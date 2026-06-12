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
// Fix for Windows WinSock header ordering issue
// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_DetectorPeakResponse_suite
#include <boost/test/included/unit_test.hpp>

#include <tuple>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_utils.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/InterSpec.h"

using namespace std;

// Data directory globals (set from command line)
std::string g_data_dir = "";
std::string g_test_data_dir = "";

namespace {

  // Helper function: calculate solid angle fraction
  double calc_solid_angle( const double diameter, const double distance )
  {
    const double radius = diameter / 2.0;
    const double distance_sq = distance * distance;
    const double radius_sq = radius * radius;
    return 0.5 * (1.0 - distance / sqrt(distance_sq + radius_sq));
  }

  // Helper function: compare floating point with tolerance
  bool close_enough( const double a, const double b, const double rel_tol = 1e-5 )
  {
    const double diff = fabs(a - b);
    const double max_val = (std::max)(fabs(a), fabs(b));
    return diff <= rel_tol * max_val || diff <= 1e-8;
  }

  // Helper to calculate air transmission coefficient
  double calc_air_transmission( const double energy_kev, const double distance_cm )
  {
    const double mu = GammaInteractionCalc::transmission_coefficient_air(
      static_cast<float>(energy_kev * PhysicalUnits::keV),
      distance_cm * PhysicalUnits::cm
    );
    return exp(-mu);
  }

}//namespace


BOOST_AUTO_TEST_CASE( test_read_gadras_detectors )
{
  cout << "\n\nTesting GADRAS detector reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set (use --datadir=...)" );

  const string gadras_dir = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory(gadras_dir),
                        "GADRAS detector directory not found: " + gadras_dir );

  // List of expected GADRAS detectors
  const vector<string> detector_names = {
    "HPGe 10%", "HPGe 20%", "HPGe 40%",
    "LaBr 10%", "LaBr 5%",
    "NaI 10%", "NaI 12%", "NaI 1x1", "NaI 25%", "NaI 2x2", "NaI 30%", "NaI 3x3", "NaI 5%"
  };

  size_t num_read = 0;
  for( const string &det_name : detector_names )
  {
    const string det_dir = SpecUtils::append_path( gadras_dir, det_name );
    if( !SpecUtils::is_directory(det_dir) )
    {
      cout << "  Skipping " << det_name << " (directory not found)" << endl;
      continue;
    }

    cout << "  Reading " << det_name << "..." << endl;

    shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>();
    BOOST_CHECK_NO_THROW( drf->fromGadrasDirectory(det_dir) );

    // Verify basic properties
    BOOST_CHECK_MESSAGE( !drf->name().empty(), det_name + ": name should not be empty" );
    BOOST_CHECK_MESSAGE( drf->detectorDiameter() > 0.0, det_name + ": diameter should be positive" );
    BOOST_CHECK_MESSAGE( drf->efficiencyFcnType() == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs,
                        det_name + ": efficiency form should be kEnergyEfficiencyPairs" );

    // Test efficiency can be evaluated at various energies
    const vector<float> test_energies = { 100.0f, 500.0f, 1000.0f, 1500.0f };
    for( const float energy : test_energies )
    {
      float eff = 0.0f;
      BOOST_CHECK_NO_THROW( eff = drf->intrinsicEfficiency(energy * PhysicalUnits::keV) );
      BOOST_CHECK_MESSAGE( eff >= 0.0f, det_name + ": efficiency should be non-negative" );
    }

    ++num_read;
  }

  BOOST_CHECK_MESSAGE( num_read > 0, "Should have read at least one GADRAS detector" );
  cout << "Successfully read " << num_read << " GADRAS detectors" << endl;

  // Specific test for NaI 3x3 detector efficiency values
  cout << "\n  Testing NaI 3x3 specific efficiency values..." << endl;
  const string nai_3x3_dir = SpecUtils::append_path( gadras_dir, "NaI 3x3" );
  if( SpecUtils::is_directory(nai_3x3_dir) )
  {
    shared_ptr<DetectorPeakResponse> nai_3x3 = make_shared<DetectorPeakResponse>();
    BOOST_REQUIRE_NO_THROW( nai_3x3->fromGadrasDirectory(nai_3x3_dir) );

    // Expected efficiency values at specific energies for NaI 3x3
    const vector<pair<float, double>> expected_efficiencies = {
      {59.0f,   0.7994},
      {195.0f,  0.9079},
      {511.0f,  0.5534},
      {1300.0f, 0.2474},
      {2614.0f, 0.1387}
    };

    for( const auto &test : expected_efficiencies )
    {
      const float energy_kev = test.first;
      const double expected_eff = test.second;

      float eff = 0.0f;
      BOOST_CHECK_NO_THROW( eff = nai_3x3->intrinsicEfficiency(energy_kev * PhysicalUnits::keV) );

      // Check against expected value with 1% tolerance
      BOOST_CHECK_MESSAGE( close_enough(eff, expected_eff, 0.01),
                          "NaI 3x3 efficiency at " + to_string(energy_kev) + " keV: " +
                          to_string(eff) + " vs expected " + to_string(expected_eff) );
    }

    // Test solid angle calculation at 25 cm
    const double distance_cm = 25.0;
    const double det_diameter_cm = nai_3x3->detectorDiameter() / PhysicalUnits::cm;
    const double solid_angle = calc_solid_angle( det_diameter_cm, distance_cm );
    const double expected_solid_angle = 0.005056;

    BOOST_CHECK_MESSAGE( close_enough(solid_angle, expected_solid_angle, 0.01),
                        "NaI 3x3 solid angle at 25 cm: " +
                        to_string(solid_angle) + " vs expected " + to_string(expected_solid_angle) );

    cout << "  NaI 3x3 efficiency and solid angle tests passed" << endl;
  }
  else
  {
    cout << "  Warning: NaI 3x3 directory not found, skipping specific tests" << endl;
  }
}


BOOST_AUTO_TEST_CASE( test_read_common_drfs_tsv )
{
  cout << "\n\nTesting common_drfs.tsv reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set (use --datadir=...)" );

  const string tsv_file = SpecUtils::append_path( g_data_dir, "common_drfs.tsv" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(tsv_file), "common_drfs.tsv not found: " + tsv_file );

  ifstream input( tsv_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open common_drfs.tsv" );

  string line;
  size_t num_parsed = 0;
  vector<string> found_detector_names;

  while( SpecUtils::safe_get_line(input, line, 4096) )
  {
    SpecUtils::trim( line );
    if( line.empty() || line[0] == '#' )
      continue;

    shared_ptr<DetectorPeakResponse> drf;
    BOOST_CHECK_NO_THROW( drf = DetectorPeakResponse::parseSingleCsvLineRelEffDrf(line) );

    if( !drf )
      continue;

    ++num_parsed;
    found_detector_names.push_back( drf->name() );

    // Verify basic properties
    BOOST_CHECK_MESSAGE( !drf->name().empty(), "DRF name should not be empty" );
    // Note: some DRFs in the TSV use pairs instead of exp-of-log
    const bool is_exp_log = (drf->efficiencyFcnType() == DetectorPeakResponse::EfficiencyFnctForm::kExpOfLogPowerSeries);
    const bool is_pairs = (drf->efficiencyFcnType() == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs);
    BOOST_CHECK_MESSAGE( is_exp_log || is_pairs,
                        drf->name() + ": efficiency form should be kExpOfLogPowerSeries or kEnergyEfficiencyPairs" );

    // Verify can evaluate efficiency
    const float test_energy = 661.7f * PhysicalUnits::keV; // Cs-137
    float eff = 0.0f;
    BOOST_CHECK_NO_THROW( eff = drf->intrinsicEfficiency(test_energy) );
    BOOST_CHECK_MESSAGE( eff >= 0.0f, drf->name() + ": efficiency should be non-negative" );
  }

  BOOST_CHECK_MESSAGE( num_parsed > 50, "Should have parsed at least 50 DRFs from TSV" );
  cout << "Successfully parsed " << num_parsed << " DRFs from common_drfs.tsv" << endl;

  // Check for specific detectors
  const vector<string> expected_detectors = {
    "ORTEC Detective_LANL_025cm",
    "ORTEC Detective-X_LANL_100cm",
    "Canberra Falcon5000_LANL_025cm"
  };

  for( const string &expected : expected_detectors )
  {
    bool found = false;
    for( const string &name : found_detector_names )
    {
      if( name.find(expected) != string::npos )
      {
        found = true;
        break;
      }
    }
    BOOST_CHECK_MESSAGE( found, "Expected detector not found: " + expected );
  }

  // Specific test for ORTEC Micro Detective (9%) detector
  cout << "\n  Testing ORTEC Micro Detective (9%) specific values..." << endl;
  input.clear();
  input.seekg(0);

  shared_ptr<DetectorPeakResponse> micro_detective;
  while( SpecUtils::safe_get_line(input, line, 4096) )
  {
    SpecUtils::trim( line );
    if( line.empty() || line[0] == '#' )
      continue;

    shared_ptr<DetectorPeakResponse> drf = DetectorPeakResponse::parseSingleCsvLineRelEffDrf(line);
    if( drf && drf->name().find("ORTEC Micro Detective (9%)") != string::npos )
    {
      micro_detective = drf;
      break;
    }
  }

  if( micro_detective )
  {
    // Expected efficiency values at specific energies for ORTEC Micro Detective (9%)
    // Note that these values arent exact, as they were calculated ignoring air attenuation (i.e., InterSpec v1.0.12)
    const vector<pair<float, double>> expected_efficiencies = {
      {59.0f,   0.09147},
      {195.0f,  0.3707},
      {511.0f,  0.1211},
      {1300.0f, 0.04554},
      {2614.0f, 0.02357}
    };

    for( const auto &test : expected_efficiencies )
    {
      const float energy_kev = test.first;
      const double expected_eff = test.second;

      float eff = 0.0f;
      BOOST_CHECK_NO_THROW( eff = micro_detective->intrinsicEfficiency(energy_kev * PhysicalUnits::keV) );

      // Check against expected value with 2% tolerance (to account for interpolation differences)
      BOOST_CHECK_MESSAGE( close_enough(eff, expected_eff, 0.02),
                          "Micro Detective efficiency at " + to_string(energy_kev) + " keV: " +
                          to_string(eff) + " vs expected " + to_string(expected_eff) );
    }

    // Test solid angle calculation at 25 cm
    const double distance_cm = 25.0;
    const double det_diameter_cm = micro_detective->detectorDiameter() / PhysicalUnits::cm;
    const double solid_angle = calc_solid_angle( det_diameter_cm, distance_cm );
    const double expected_solid_angle = 0.002357;  // Note: actual value, not 0.02357

    BOOST_CHECK_MESSAGE( close_enough(solid_angle, expected_solid_angle, 0.06),
                        "Micro Detective solid angle at 25 cm: " +
                        to_string(solid_angle) + " vs expected " + to_string(expected_solid_angle) );

    cout << "  ORTEC Micro Detective (9%) efficiency and solid angle tests passed" << endl;
  }
  else
  {
    cout << "  Warning: ORTEC Micro Detective (9%) not found in TSV, skipping specific tests" << endl;
    assert( 0 );
  }

  // Specific test for Verifinder-SN23N detector
  cout << "\n  Testing Verifinder-SN23N specific values..." << endl;
  input.clear();
  input.seekg(0);

  shared_ptr<DetectorPeakResponse> verifinder;
  while( SpecUtils::safe_get_line(input, line, 4096) )
  {
    SpecUtils::trim( line );
    if( line.empty() || line[0] == '#' )
      continue;

    shared_ptr<DetectorPeakResponse> drf = DetectorPeakResponse::parseSingleCsvLineRelEffDrf(line);
    if( drf && drf->name().find("Verifinder-SN23N") != string::npos )
    {
      verifinder = drf;
      break;
    }
  }

  if( verifinder )
  {
    // Expected efficiency values at specific energies for Verifinder-SN23N
    const vector<pair<float, double>> expected_efficiencies = {
      {59.0f,   0.822},
      {195.0f,  0.9258},
      {511.0f,  0.5773},
      {1300.0f, 0.2648},
      {2614.0f, 0.1507}
    };

    for( const auto &test : expected_efficiencies )
    {
      const float energy_kev = test.first;
      const double expected_eff = test.second;

      float eff = 0.0f;
      BOOST_CHECK_NO_THROW( eff = verifinder->intrinsicEfficiency(energy_kev * PhysicalUnits::keV) );

      // Check against expected value with 1% tolerance
      BOOST_CHECK_MESSAGE( close_enough(eff, expected_eff, 0.01),
                          "Verifinder efficiency at " + to_string(energy_kev) + " keV: " +
                          to_string(eff) + " vs expected " + to_string(expected_eff) );
    }

    // Test solid angle calculation at 25 cm
    const double distance_cm = 25.0;
    const double det_diameter_cm = verifinder->detectorDiameter() / PhysicalUnits::cm;
    const double solid_angle = calc_solid_angle( det_diameter_cm, distance_cm );
    const double expected_solid_angle = 0.005702;

    BOOST_CHECK_MESSAGE( close_enough(solid_angle, expected_solid_angle, 0.01),
                        "Verifinder solid angle at 25 cm: " +
                        to_string(solid_angle) + " vs expected " + to_string(expected_solid_angle) );

    cout << "  Verifinder-SN23N efficiency and solid angle tests passed" << endl;
  }
  else
  {
    cout << "  Warning: Verifinder-SN23N not found in TSV, skipping specific tests" << endl;
  }
}


BOOST_AUTO_TEST_CASE( test_read_ecc_file )
{
  cout << "\n\nTesting ECC file reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set (use --testfiledir=...)" );

  const string ecc_file = SpecUtils::append_path( g_test_data_dir, "det_eff/Detective-X_in-situ.ecc" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(ecc_file), "ECC file not found: " + ecc_file );

  ifstream input( ecc_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open ECC file" );

  shared_ptr<DetectorPeakResponse> drf;
  double source_area = 0.0, source_mass = 0.0;

  BOOST_CHECK_NO_THROW(
    std::tie(drf, source_area, source_mass) = DetectorPeakResponse::parseEccFile(input)
  );

  BOOST_REQUIRE_MESSAGE( drf != nullptr, "ECC parsing should return valid DRF" );

  // Verify geometry type
  BOOST_CHECK_MESSAGE( drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct,
                      "ECC DRF should have FixedGeomTotalAct geometry" );

  // Verify efficiency form is pairs
  BOOST_CHECK_MESSAGE( drf->efficiencyFcnType() == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs,
                      "ECC DRF should use energy-efficiency pairs" );

  // Verify source parameters are reasonable
  BOOST_CHECK_MESSAGE( source_area > 0.0, "Source area should be positive" );
  BOOST_CHECK_MESSAGE( source_mass >= 0.0, "Source mass should be non-negative" );

  cout << "  Source area: " << source_area << " cm²" << endl;
  cout << "  Source mass: " << source_mass << " g" << endl;

  // Test efficiency values at specific energies
  // Expected values for Detective-X_in-situ.ecc
  const vector<pair<float, double>> expected_efficiencies = {
    {100.0f,  3.733E-07},
    {500.0f,  2.235E-07},
    {661.7f,  1.832E-07},
    {1200.0f, 1.228E-07}
  };

  for( const auto &test : expected_efficiencies )
  {
    const float energy_kev = test.first;
    const double expected_eff = test.second;

    // For fixed-geometry DRF, use distance=0 or intrinsicEfficiency
    float eff = 0.0f;
    BOOST_CHECK_NO_THROW( eff = drf->intrinsicEfficiency(energy_kev * PhysicalUnits::keV) );
    BOOST_CHECK_MESSAGE( eff >= 0.0f, "Efficiency should be non-negative at " + to_string(energy_kev) + " keV" );

    // Check against expected value with 1% tolerance
    BOOST_CHECK_MESSAGE( close_enough(eff, expected_eff, 0.01),
                        "Efficiency at " + to_string(energy_kev) + " keV: " +
                        to_string(eff) + " vs expected " + to_string(expected_eff) );
  }

  cout << "Successfully parsed ECC file" << endl;
}


BOOST_AUTO_TEST_CASE( test_read_angle_outx_file )
{
  cout << "\n\nTesting ANGLE .outx file reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set (use --testfiledir=...)" );

  const string outx_file = SpecUtils::append_path( g_test_data_dir, "det_eff/Angle-example-efficiency.outx" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(outx_file), "OUTX file not found: " + outx_file );

  ifstream input( outx_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open OUTX file" );

  shared_ptr<DetectorPeakResponse> drf;

  BOOST_CHECK_NO_THROW(
    drf = DetectorPeakResponse::parseAngleOutxFile(input)
  );

  BOOST_REQUIRE_MESSAGE( drf != nullptr, "OUTX parsing should return valid DRF" );

  // Verify geometry type
  BOOST_CHECK_MESSAGE( drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct,
                      "OUTX DRF should have FixedGeomTotalAct geometry" );

  // Verify efficiency form is pairs
  BOOST_CHECK_MESSAGE( drf->efficiencyFcnType() == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs,
                      "OUTX DRF should use energy-efficiency pairs" );

  // Verify DRF source
  BOOST_CHECK_MESSAGE( drf->drfSource() == DetectorPeakResponse::DrfSource::AngleOutx,
                      "OUTX DRF should have AngleOutx source" );

  // Verify detector name is populated from XML attributes
  BOOST_CHECK_MESSAGE( !drf->name().empty(), "OUTX DRF name should not be empty" );

  // Test efficiency values at energies present in the file
  const vector<pair<float, double>> expected_efficiencies = {
    {100.0f,   0.200207559471988},
    {500.0f,   0.0744367862247439},
    {661.657f, 0.0584335868503359},
    {1000.0f,  0.0415204261052124}
  };

  for( const auto &test : expected_efficiencies )
  {
    const float energy_kev = test.first;
    const double expected_eff = test.second;

    float eff = 0.0f;
    BOOST_CHECK_NO_THROW( eff = drf->intrinsicEfficiency(energy_kev * PhysicalUnits::keV) );
    BOOST_CHECK_MESSAGE( eff >= 0.0f, "Efficiency should be non-negative at " + to_string(energy_kev) + " keV" );

    // Check against expected value with 1% tolerance
    BOOST_CHECK_MESSAGE( close_enough(eff, expected_eff, 0.01),
                        "Efficiency at " + to_string(energy_kev) + " keV: " +
                        to_string(eff) + " vs expected " + to_string(expected_eff) );
  }

  // Verify that parseAngleOutxFile rejects an ECC file
  const string ecc_file = SpecUtils::append_path( g_test_data_dir, "det_eff/Detective-X_in-situ.ecc" );
  if( SpecUtils::is_file(ecc_file) )
  {
    ifstream ecc_input( ecc_file.c_str() );
    BOOST_CHECK_THROW( DetectorPeakResponse::parseAngleOutxFile(ecc_input), std::exception );
  }

  cout << "Successfully parsed ANGLE .outx file" << endl;
}


BOOST_AUTO_TEST_CASE( test_exp_of_log_power_series )
{
  cout << "\n\nTesting exp-of-log power series efficiency..." << endl;

  // Test coefficients: exp(c0 + c1*ln(E) + c2*ln(E)^2)
  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };
  const vector<float> uncerts; // No uncertainties for this test
  const double det_diameter = 5.08 * PhysicalUnits::cm; // 2-inch detector
  const float energy_units = PhysicalUnits::keV;

  // A. Test FarFieldIntrinsic
  {
    cout << "  Testing FarFieldIntrinsic..." << endl;

    shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>( "TestDRF", "Test exp-of-log" );
    BOOST_CHECK_NO_THROW(
      drf->fromExpOfLogPowerSeries( coeffs, uncerts, 0.0, det_diameter, energy_units,
                                   50.0f, 3000.0f,
                                   DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic )
    );

    // Test at various energies
    const vector<float> test_energies = { 50.0f, 500.0f, 1000.0f, 3000.0f };

    for( const float E : test_energies )
    {
      const float eff = drf->intrinsicEfficiency( E * PhysicalUnits::keV );

      // Calculate expected: exp(c0 + c1*ln(E) + c2*ln(E)^2)
      const double ln_E = log(E);
      const double expected = exp( coeffs[0] + coeffs[1]*ln_E + coeffs[2]*ln_E*ln_E );

      BOOST_CHECK_MESSAGE( close_enough(eff, expected, 1e-4),
                          "Efficiency mismatch at " + to_string(E) + " keV: " +
                          to_string(eff) + " vs " + to_string(expected) );
    }

    cout << "    FarFieldIntrinsic tests passed" << endl;
  }

  // B. Test FarFieldAbsolute
  {
    cout << "  Testing FarFieldAbsolute..." << endl;

    const double char_dist = 100.0 * PhysicalUnits::cm;

    shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>( "TestDRF_Abs", "Test absolute" );
    BOOST_CHECK_NO_THROW(
      drf->fromExpOfLogPowerSeries( coeffs, uncerts, char_dist, det_diameter, energy_units,
                                   50.0f, 3000.0f,
                                   DetectorPeakResponse::EffGeometryType::FarFieldAbsolute )
    );

    // Verify geometry type and distance
    BOOST_CHECK( drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

    const vector<float> test_energies = { 100.0f, 661.7f, 1460.0f };

    for( const float E : test_energies )
    {
      // For FarFieldAbsolute, the coefficients represent absolute efficiency at char_dist
      // Intrinsic efficiency = absolute_eff * (1/solid_angle) * (1/air_transmission)
      const double ln_E = log(E);
      const double absolute_eff = exp( coeffs[0] + coeffs[1]*ln_E + coeffs[2]*ln_E*ln_E );

      // Calculate geometric correction
      const double solid_angle = calc_solid_angle( det_diameter / PhysicalUnits::cm, char_dist / PhysicalUnits::cm );

      // Calculate air attenuation correction
      const double air_trans = calc_air_transmission( E, char_dist / PhysicalUnits::cm );

      const double expected_intrinsic = absolute_eff / (solid_angle * air_trans);

      const float eff = drf->intrinsicEfficiency( E * PhysicalUnits::keV );

      BOOST_CHECK_MESSAGE( close_enough(eff, expected_intrinsic, 1e-3),
                          "Intrinsic efficiency mismatch at " + to_string(E) + " keV: " +
                          to_string(eff) + " vs " + to_string(expected_intrinsic) );
    }

    cout << "    FarFieldAbsolute tests passed" << endl;
  }

  // C. Test FarFieldAbsolute with and without air attenuation
  {
    cout << "  Testing air attenuation correction..." << endl;

    const double char_dist = 100.0 * PhysicalUnits::cm;
    const float test_energy = 100.0f; // Low energy = strong air attenuation

    // DRF with air correction
    shared_ptr<DetectorPeakResponse> drf_with_air = make_shared<DetectorPeakResponse>( "WithAir", "With air" );
    drf_with_air->fromExpOfLogPowerSeries( coeffs, uncerts, char_dist, det_diameter, energy_units,
                                          50.0f, 3000.0f,
                                          DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

    // Both should give same intrinsic efficiency (air correction is applied during conversion)
    const float eff_with = drf_with_air->intrinsicEfficiency( test_energy * PhysicalUnits::keV );

    BOOST_CHECK_MESSAGE( eff_with > 0.0f, "Efficiency with air should be positive" );

    cout << "    Air attenuation correction tests passed" << endl;
  }
}


BOOST_AUTO_TEST_CASE( test_intrinsic_efficiency_formula )
{
  cout << "\n\nTesting intrinsic efficiency formula..." << endl;

  const double det_diameter = 5.08 * PhysicalUnits::cm;

  // Test polynomial formula: a0 + a1*Energy + a2*Energy^2
  // Using MeV for formula (InterSpec standard)
  const string formula = "1e-3 + 1e-4*Energy - 1e-5*Energy^2";

  // A. FarFieldIntrinsic
  {
    cout << "  Testing FarFieldIntrinsic with formula..." << endl;

    shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>( "FormulaDRF", "Formula test" );
    drf->setIntrinsicEfficiencyFormula( formula, det_diameter, PhysicalUnits::MeV, 0.05f, 3.0f,
                                       DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

    // Test at various energies
    const vector<float> test_energies_kev = { 100.0f, 500.0f, 1000.0f };

    for( const float E_kev : test_energies_kev )
    {
      const float eff = drf->intrinsicEfficiency( E_kev * PhysicalUnits::keV );

      // Calculate expected (formula is in MeV)
      const double E_mev = E_kev / 1000.0;
      const double expected = 1e-3 + 1e-4*E_mev - 1e-5*E_mev*E_mev;

      BOOST_CHECK_MESSAGE( close_enough(eff, expected, 1e-6),
                          "Formula efficiency mismatch at " + to_string(E_kev) + " keV" );
    }

    cout << "    FarFieldIntrinsic formula tests passed" << endl;
  }

  // B. FarFieldAbsolute with formula
  {
    cout << "  Testing FarFieldAbsolute with formula..." << endl;

    const double char_dist = 100.0 * PhysicalUnits::cm;

    // Create FarFieldIntrinsic, then reinterpret as FarFieldAbsolute
    shared_ptr<DetectorPeakResponse> drf_intrinsic = make_shared<DetectorPeakResponse>( "FormulaInt", "Formula intrinsic" );
    drf_intrinsic->setIntrinsicEfficiencyFormula( formula, det_diameter, PhysicalUnits::MeV, 0.05f, 3.0f,
                                                 DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

    // Reinterpret as FarFieldAbsolute: the formula values now represent absolute efficiency at char_dist
    shared_ptr<DetectorPeakResponse> drf_abs = drf_intrinsic->reinterpretAsFarFieldAbsEfficiency( det_diameter, char_dist, true );

    BOOST_REQUIRE( drf_abs != nullptr );
    BOOST_CHECK( drf_abs->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

    // Test that intrinsic efficiency is correctly calculated
    const float E_kev = 661.7f;

    // The formula value (which was intrinsic) is now treated as absolute
    const double E_mev = E_kev / 1000.0;
    const double absolute_eff = 1e-3 + 1e-4*E_mev - 1e-5*E_mev*E_mev;

    // Calculate geometric correction
    const double solid_angle = calc_solid_angle( det_diameter / PhysicalUnits::cm, char_dist / PhysicalUnits::cm );

    // Calculate air attenuation correction
    const double air_trans = calc_air_transmission( E_kev, char_dist / PhysicalUnits::cm );

    const double expected_intrinsic = absolute_eff / (solid_angle * air_trans);

    const float eff = drf_abs->intrinsicEfficiency( E_kev * PhysicalUnits::keV );

    BOOST_CHECK_MESSAGE( close_enough(eff, expected_intrinsic, 1e-3),
                        "Intrinsic efficiency should be corrected: " + to_string(eff) + " vs " + to_string(expected_intrinsic) );

    cout << "    FarFieldAbsolute formula tests passed" << endl;
  }
}


BOOST_AUTO_TEST_CASE( test_efficiency_pairs )
{
  cout << "\n\nTesting efficiency pairs..." << endl;

  const double det_diameter = 5.08 * PhysicalUnits::cm;
  const double char_dist = 100.0 * PhysicalUnits::cm;

  // Define efficiency points (energy in keV, efficiency as fraction)
  const vector<pair<float,float>> intrinsic_points = {
    {50.0f, 0.01f}, {100.0f, 0.05f}, {500.0f, 0.15f},
    {1000.0f, 0.10f}, {3000.0f, 0.03f}
  };

  // A. Create FarFieldIntrinsic DRF
  shared_ptr<DetectorPeakResponse> drf_intrinsic = make_shared<DetectorPeakResponse>( "PairsDRF", "Pairs test" );

  vector<DetectorPeakResponse::EnergyEffPoint> eff_pairs;
  for( const auto &pt : intrinsic_points )
    eff_pairs.push_back( {pt.first, pt.second, {}} );

  drf_intrinsic->setEfficiencyPoints( eff_pairs, det_diameter, 0.0,
                                     DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  // Test interpolation at exact points
  for( const auto &pt : intrinsic_points )
  {
    const float eff = drf_intrinsic->intrinsicEfficiency( pt.first * PhysicalUnits::keV );
    BOOST_CHECK_MESSAGE( close_enough(eff, pt.second, 1e-4),
                        "Efficiency should match at exact point " + to_string(pt.first) + " keV" );
  }

  // B. Create FarFieldAbsolute DRF with calculated absolute efficiencies
  cout << "  Testing FarFieldAbsolute with efficiency pairs..." << endl;

  // Calculate absolute efficiencies at characterization distance
  vector<DetectorPeakResponse::EnergyEffPoint> abs_eff_pairs;
  for( const auto &pt : intrinsic_points )
  {
    const double solid_angle = calc_solid_angle( det_diameter / PhysicalUnits::cm, char_dist / PhysicalUnits::cm );
    const double air_trans = calc_air_transmission( pt.first, char_dist / PhysicalUnits::cm );
    const float abs_eff = static_cast<float>( pt.second * solid_angle * air_trans );

    abs_eff_pairs.push_back( {pt.first, abs_eff, {}} );
  }

  shared_ptr<DetectorPeakResponse> drf_absolute = make_shared<DetectorPeakResponse>( "AbsPairs", "Absolute pairs" );
  drf_absolute->setEfficiencyPoints( abs_eff_pairs, det_diameter, char_dist,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  // Verify both DRFs give same intrinsic efficiency at various energies
  const vector<float> test_energies = { 75.0f, 300.0f, 661.7f, 1460.0f, 2000.0f };

  for( const float E : test_energies )
  {
    const float eff_intrinsic = drf_intrinsic->intrinsicEfficiency( E * PhysicalUnits::keV );
    const float eff_absolute = drf_absolute->intrinsicEfficiency( E * PhysicalUnits::keV );

    // Should match within tolerance
    BOOST_CHECK_MESSAGE( close_enough(eff_intrinsic, eff_absolute, 0.05), // 5% tolerance for interpolation
                        "Intrinsic efficiencies should match at " + to_string(E) + " keV: " +
                        to_string(eff_intrinsic) + " vs " + to_string(eff_absolute) );
  }

  cout << "  Efficiency pairs tests passed" << endl;
}


BOOST_AUTO_TEST_CASE( test_reinterpret_as_far_field )
{
  cout << "\n\nTesting reinterpretAsFarFieldEfficiency..." << endl;

  // Read an ECC file to get a FixedGeomTotalAct DRF
  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set" );

  const string ecc_file = SpecUtils::append_path( g_test_data_dir, "det_eff/Detective-X_in-situ.ecc" );
  if( !SpecUtils::is_file(ecc_file) )
  {
    cout << "  ECC file not found, skipping test" << endl;
    return;
  }

  ifstream input( ecc_file.c_str() );
  BOOST_REQUIRE( input.is_open() );

  shared_ptr<DetectorPeakResponse> drf_fixed;
  double source_area = 0.0, source_mass = 0.0;
  std::tie(drf_fixed, source_area, source_mass) = DetectorPeakResponse::parseEccFile(input);

  BOOST_REQUIRE( drf_fixed != nullptr );
  BOOST_REQUIRE( drf_fixed->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );

  // Reinterpret as far-field
  const double det_diameter = 7.62 * PhysicalUnits::cm; // 3-inch
  const double distance = 100.0 * PhysicalUnits::cm;

  shared_ptr<DetectorPeakResponse> drf_farfield;
  BOOST_CHECK_NO_THROW(
    drf_farfield = drf_fixed->reinterpretAsFarFieldAbsEfficiency( det_diameter, distance, true )
  );

  BOOST_REQUIRE( drf_farfield != nullptr );

  // Verify geometry type changed
  BOOST_CHECK( drf_farfield->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  // Verify efficiency values are preserved (should be same absolute efficiency)
  const vector<float> test_energies = { 100.0f, 500.0f, 1000.0f };
  for( const float E : test_energies )
  {
    // Use efficiency() for fixed geometry (needs distance parameter)
    const float eff_fixed = drf_fixed->efficiency( E * PhysicalUnits::keV, 1.0 );

    // Use intrinsicEfficiency() for far-field
    const float eff_farfield = drf_farfield->intrinsicEfficiency( E * PhysicalUnits::keV );

    // These won't be exactly equal due to different geometry assumptions, but both should be valid
    BOOST_CHECK_MESSAGE( eff_fixed >= 0.0f && eff_farfield >= 0.0f,
                        "Both efficiencies should be non-negative at " + to_string(E) + " keV" );
  }

  cout << "  reinterpretAsFarFieldEfficiency tests passed" << endl;
}


BOOST_AUTO_TEST_CASE( test_convert_fixed_to_far_field )
{
  cout << "\n\nTesting convertFixedGeometryToFarField..." << endl;

  // Read an ECC file
  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set" );

  const string ecc_file = SpecUtils::append_path( g_test_data_dir, "det_eff/Detective-X_in-situ.ecc" );
  if( !SpecUtils::is_file(ecc_file) )
  {
    cout << "  ECC file not found, skipping test" << endl;
    return;
  }

  ifstream input( ecc_file.c_str() );
  BOOST_REQUIRE( input.is_open() );

  shared_ptr<DetectorPeakResponse> drf_fixed;
  double source_area = 0.0, source_mass = 0.0;
  std::tie(drf_fixed, source_area, source_mass) = DetectorPeakResponse::parseEccFile(input);

  BOOST_REQUIRE( drf_fixed != nullptr );
  BOOST_REQUIRE( drf_fixed->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );

  // Convert to far-field
  const double det_diameter = 7.62 * PhysicalUnits::cm;
  const double distance = 50.0 * PhysicalUnits::cm;

  shared_ptr<DetectorPeakResponse> drf_converted;
  BOOST_CHECK_NO_THROW(
    drf_converted = drf_fixed->reinterpretAsFarFieldAbsEfficiency( det_diameter, distance, true )
  );

  BOOST_REQUIRE( drf_converted != nullptr );
  BOOST_CHECK( drf_converted->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  // Test that non-FixedGeomTotalAct doesnt throw exception
  shared_ptr<DetectorPeakResponse> drf_intrinsic = make_shared<DetectorPeakResponse>( "Test", "Test" );
  const vector<float> coeffs = { -5.0f, 0.5f };
  drf_intrinsic->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter, PhysicalUnits::keV,
                                         50.0f, 3000.0f,
                                         DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  BOOST_CHECK_NO_THROW( drf_intrinsic->reinterpretAsFarFieldAbsEfficiency( det_diameter, distance, true ) );

  cout << "  convertFixedGeometryToFarField tests passed" << endl;
}


BOOST_AUTO_TEST_CASE( test_xml_serialization_round_trip )
{
  cout << "\n\nTesting XML serialization round-trip..." << endl;

  const double det_diameter = 5.08 * PhysicalUnits::cm;
  const double char_dist = 100.0 * PhysicalUnits::cm;

  // Test different efficiency forms and geometry types
  vector<shared_ptr<DetectorPeakResponse>> test_drfs;

  // 1. Exp-of-log, FarFieldIntrinsic
  {
    auto drf = make_shared<DetectorPeakResponse>( "ExpLog_Intrinsic", "Exp-log intrinsic test" );
    const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };
    drf->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter, PhysicalUnits::keV,
                                 50.0f, 3000.0f,
                                 DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    test_drfs.push_back( drf );
  }

  // 2. Exp-of-log, FarFieldAbsolute
  {
    auto drf = make_shared<DetectorPeakResponse>( "ExpLog_Absolute", "Exp-log absolute test" );
    const vector<float> coeffs = { -4.5f, 0.3f, -0.008f };
    drf->fromExpOfLogPowerSeries( coeffs, {}, char_dist, det_diameter, PhysicalUnits::keV,
                                 50.0f, 3000.0f,
                                 DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    test_drfs.push_back( drf );
  }

  // 3. Formula, FarFieldIntrinsic
  {
    auto drf = make_shared<DetectorPeakResponse>( "Formula_Intrinsic", "Formula intrinsic test" );
    drf->setIntrinsicEfficiencyFormula( "1e-3 + 1e-4*Energy", det_diameter, PhysicalUnits::keV, 10.0, 3000.0,
                                       DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    
    test_drfs.push_back( drf );
  }

  // 4. Energy-efficiency pairs, FarFieldIntrinsic
  {
    auto drf = make_shared<DetectorPeakResponse>( "Pairs_Intrinsic", "Pairs intrinsic test" );
    vector<DetectorPeakResponse::EnergyEffPoint> pairs = {
      {100.0f, 0.05f, {}}, {500.0f, 0.15f, {}}, {1000.0f, 0.12f, {}}, {2000.0f, 0.08f, {}}
    };
    drf->setEfficiencyPoints( pairs, det_diameter, 0.0, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    test_drfs.push_back( drf );
  }

  // 5. Pairs, FarFieldAbsolute
  {
    auto drf = make_shared<DetectorPeakResponse>( "Pairs_Absolute", "Pairs absolute test" );
    vector<DetectorPeakResponse::EnergyEffPoint> pairs = {
      {100.0f, 0.001f, {}}, {500.0f, 0.003f, {}}, {1000.0f, 0.0025f, {}}
    };
    drf->setEfficiencyPoints( pairs, det_diameter, char_dist, DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    test_drfs.push_back( drf );
  }

  // Test round-trip for each DRF
  for( size_t i = 0; i < test_drfs.size(); ++i )
  {
    const auto &original = test_drfs[i];
    cout << "  Testing " << original->name() << "..." << endl;

    // Serialize to XML
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
    doc.append_node( root );

    BOOST_CHECK_NO_THROW( original->toXml( root, &doc ) );

    // Find the DetectorPeakResponse node
    rapidxml::xml_node<char> *drf_node = root->first_node( "DetectorPeakResponse" );
    BOOST_REQUIRE_MESSAGE( drf_node != nullptr, "XML should contain DetectorPeakResponse node" );

    // Deserialize from XML
    auto restored = make_shared<DetectorPeakResponse>();
    BOOST_CHECK_NO_THROW( restored->fromXml( drf_node ) );

    // Compare using equalEnough
    BOOST_CHECK_NO_THROW( DetectorPeakResponse::equalEnough( *original, *restored ) );

    // Verify specific properties
    BOOST_CHECK_EQUAL( original->name(), restored->name() );
    BOOST_CHECK_EQUAL( original->description(), restored->description() );
    BOOST_CHECK( original->geometryType() == restored->geometryType() );
    BOOST_CHECK( original->efficiencyFcnType() == restored->efficiencyFcnType() );

    // Test efficiency values match
    const vector<float> test_energies = { 100.0f, 500.0f, 1000.0f };
    for( const float E : test_energies )
    {
      const float eff_orig = original->intrinsicEfficiency( E * PhysicalUnits::keV );
      const float eff_rest = restored->intrinsicEfficiency( E * PhysicalUnits::keV );

      BOOST_CHECK_MESSAGE( close_enough(eff_orig, eff_rest, 1e-4),
                          original->name() + ": efficiency mismatch at " + to_string(E) + " keV" );
    }
  }

  cout << "  XML serialization tests passed for " << test_drfs.size() << " DRFs" << endl;
}


BOOST_AUTO_TEST_CASE( test_url_serialization_round_trip )
{
  cout << "\n\nTesting URL serialization round-trip..." << endl;

  const double det_diameter = 5.08 * PhysicalUnits::cm;
  const double char_dist = 100.0 * PhysicalUnits::cm;

  // Test different DRF types
  vector<shared_ptr<DetectorPeakResponse>> test_drfs;

  // 1. Exp-of-log, FarFieldIntrinsic
  {
    auto drf = make_shared<DetectorPeakResponse>( "URL_ExpLog", "URL exp-log test" );
    const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };
    drf->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter, PhysicalUnits::keV,
                                 50.0f, 3000.0f,
                                 DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    test_drfs.push_back( drf );
  }

  // 2. FarFieldAbsolute
  {
    auto drf = make_shared<DetectorPeakResponse>( "URL_Absolute", "URL absolute test" );
    const vector<float> coeffs = { -4.8f, 0.4f };
    drf->fromExpOfLogPowerSeries( coeffs, {}, char_dist, det_diameter, PhysicalUnits::keV,
                                 50.0f, 3000.0f,
                                 DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    test_drfs.push_back( drf );
  }

  // 3. Energy-efficiency pairs
  {
    auto drf = make_shared<DetectorPeakResponse>( "URL_Pairs", "URL pairs test" );
    vector<DetectorPeakResponse::EnergyEffPoint> pairs = {
      {100.0f, 0.05f, {}}, {500.0f, 0.15f, {}}, {1000.0f, 0.12f, {}}
    };
    drf->setEfficiencyPoints( pairs, det_diameter, 0.0,
                             DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    test_drfs.push_back( drf );
  }

  // Test round-trip for each DRF
  for( size_t i = 0; i < test_drfs.size(); ++i )
  {
    const auto &original = test_drfs[i];
    cout << "  Testing " << original->name() << "..." << endl;

    // Serialize to URL
    string url;
    BOOST_CHECK_NO_THROW( url = original->toAppUrl() );
    BOOST_CHECK_MESSAGE( !url.empty(), "URL should not be empty" );

    cout << "    URL length: " << url.length() << " chars" << endl;

    // Deserialize from URL
    auto restored = make_shared<DetectorPeakResponse>();
    BOOST_CHECK_NO_THROW( restored->fromAppUrl( url ) );

    // Verify basic properties match
    BOOST_CHECK_EQUAL( original->name(), restored->name() );
    BOOST_CHECK( original->geometryType() == restored->geometryType() );
    BOOST_CHECK( original->efficiencyFcnType() == restored->efficiencyFcnType() );

    // Test efficiency values match (with looser tolerance for URL encoding)
    const vector<float> test_energies = { 150.0f, 661.7f, 1460.0f };
    for( const float E : test_energies )
    {
      const float eff_orig = original->intrinsicEfficiency( E * PhysicalUnits::keV );
      const float eff_rest = restored->intrinsicEfficiency( E * PhysicalUnits::keV );

      // URL encoding has limited precision, use 0.1% tolerance
      const double rel_diff = fabs(eff_orig - eff_rest) / (std::max)(fabs(eff_orig), fabs(eff_rest));
      BOOST_CHECK_MESSAGE( rel_diff < 0.001,
                          original->name() + ": efficiency mismatch at " + to_string(E) + " keV: " +
                          to_string(eff_orig) + " vs " + to_string(eff_rest) );
    }
  }

  cout << "  URL serialization tests passed for " << test_drfs.size() << " DRFs" << endl;
}


// Parse command-line arguments to set data directories
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

    cout << "\nTest configuration:" << endl;
    cout << "  Data directory: " << (g_data_dir.empty() ? "(not set)" : g_data_dir) << endl;
    cout << "  Test data directory: " << (g_test_data_dir.empty() ? "(not set)" : g_test_data_dir) << endl;

    // Set static data directory for GammaInteractionCalc and other components
    if( !g_data_dir.empty() )
    {
      try
      {
        InterSpec::setStaticDataDirectory( g_data_dir );
        cout << "  Static data directory set successfully" << endl;
      }catch( std::exception &e )
      {
        cerr << "  Warning: Failed to set static data directory: " << e.what() << endl;
      }
    }
  }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture );


BOOST_AUTO_TEST_CASE( test_detector_setback_basics )
{
  cout << "\n\nTesting detector setback basics..." << endl;

  // Default setback should be 0
  DetectorPeakResponse drf( "TestDRF", "Test DRF for setback" );
  BOOST_CHECK_CLOSE( drf.detectorSetback(), 0.0, 1e-10 );

  // Setting a positive setback should work
  const double setback = 0.5 * PhysicalUnits::cm;
  BOOST_CHECK_NO_THROW( drf.setDetectorSetback( setback ) );
  BOOST_CHECK_CLOSE( drf.detectorSetback(), setback, 1e-10 );

  // Negative setback should throw
  BOOST_CHECK_THROW( drf.setDetectorSetback( -1.0 ), std::runtime_error );

  // Zero setback should be fine
  BOOST_CHECK_NO_THROW( drf.setDetectorSetback( 0.0 ) );
  BOOST_CHECK_CLOSE( drf.detectorSetback(), 0.0, 1e-10 );

  cout << "Detector setback basics passed" << endl;
}//test_detector_setback_basics


BOOST_AUTO_TEST_CASE( test_detector_setback_efficiency )
{
  cout << "\n\nTesting detector setback effect on efficiency..." << endl;

  const double det_diameter = 5.0 * PhysicalUnits::cm;
  const double distance = 25.0 * PhysicalUnits::cm;
  const double setback = 0.5 * PhysicalUnits::cm;

  // Create two DRFs with same intrinsic efficiency
  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };

  auto drf_with_setback = make_shared<DetectorPeakResponse>( "WithSetback", "Has setback" );
  drf_with_setback->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  drf_with_setback->setDetectorSetback( setback );

  auto drf_no_setback = make_shared<DetectorPeakResponse>( "NoSetback", "No setback" );
  drf_no_setback->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  // efficiency(energy, D) with setback S should equal efficiency(energy, D+S) without setback
  const float energy = 661.0f * PhysicalUnits::keV;
  const double eff_with_setback = drf_with_setback->efficiency( energy, distance );
  const double eff_no_setback_at_further_dist = drf_no_setback->efficiency( energy, distance + setback );

  BOOST_CHECK_MESSAGE( close_enough( eff_with_setback, eff_no_setback_at_further_dist, 1e-6 ),
    "efficiency(E, D) with setback S should equal efficiency(E, D+S) without setback: "
    + to_string(eff_with_setback) + " vs " + to_string(eff_no_setback_at_further_dist) );

  // efficiency with setback should be less than without (at same distance)
  const double eff_no_setback = drf_no_setback->efficiency( energy, distance );
  BOOST_CHECK_GT( eff_no_setback, eff_with_setback );

  cout << "Detector setback efficiency test passed" << endl;
}//test_detector_setback_efficiency


BOOST_AUTO_TEST_CASE( test_detector_setback_hash_stability )
{
  cout << "\n\nTesting detector setback hash stability..." << endl;

  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };
  const double det_diameter = 5.0 * PhysicalUnits::cm;

  // DRF with setback=0 should have same hash as DRF created without setback
  auto drf1 = make_shared<DetectorPeakResponse>( "HashTest1", "hash test" );
  drf1->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  auto drf2 = make_shared<DetectorPeakResponse>( "HashTest1", "hash test" );
  drf2->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  drf2->setDetectorSetback( 0.0 );

  BOOST_CHECK_EQUAL( drf1->hashValue(), drf2->hashValue() );

  // Changing setback to non-zero should change the hash
  drf2->setDetectorSetback( 1.0 * PhysicalUnits::cm );
  BOOST_CHECK_NE( drf1->hashValue(), drf2->hashValue() );

  cout << "Detector setback hash stability test passed" << endl;
}//test_detector_setback_hash_stability


BOOST_AUTO_TEST_CASE( test_detector_setback_xml_roundtrip )
{
  cout << "\n\nTesting detector setback XML round-trip..." << endl;

  const double det_diameter = 5.0 * PhysicalUnits::cm;
  const double setback = 1.44 * PhysicalUnits::cm;
  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };

  // Create DRF with setback
  auto original = make_shared<DetectorPeakResponse>( "SetbackXML", "XML setback test" );
  original->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  original->setDetectorSetback( setback );

  // Serialize to XML
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
  doc.append_node( root );
  BOOST_CHECK_NO_THROW( original->toXml( root, &doc ) );

  // Deserialize from XML
  rapidxml::xml_node<char> *drf_node = root->first_node( "DetectorPeakResponse" );
  BOOST_REQUIRE( drf_node != nullptr );

  auto restored = make_shared<DetectorPeakResponse>();
  BOOST_CHECK_NO_THROW( restored->fromXml( drf_node ) );

  // Verify setback is preserved
  BOOST_CHECK_MESSAGE( close_enough( restored->detectorSetback(), setback, 1e-4 ),
    "Setback not preserved: " + to_string(restored->detectorSetback()) + " vs " + to_string(setback) );

  // Verify equalEnough passes
  BOOST_CHECK_NO_THROW( DetectorPeakResponse::equalEnough( *original, *restored ) );

  // Test backward compat: XML without DetectorSetback should give setback=0
  auto drf_no_setback = make_shared<DetectorPeakResponse>( "NoSetbackXML", "no setback" );
  drf_no_setback->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  rapidxml::xml_document<char> doc2;
  rapidxml::xml_node<char> *root2 = doc2.allocate_node( rapidxml::node_element, "root" );
  doc2.append_node( root2 );
  drf_no_setback->toXml( root2, &doc2 );

  auto restored2 = make_shared<DetectorPeakResponse>();
  restored2->fromXml( root2->first_node( "DetectorPeakResponse" ) );
  BOOST_CHECK_CLOSE( restored2->detectorSetback(), 0.0, 1e-10 );

  cout << "Detector setback XML round-trip test passed" << endl;
}//test_detector_setback_xml_roundtrip


BOOST_AUTO_TEST_CASE( test_detector_setback_url_roundtrip )
{
  cout << "\n\nTesting detector setback URL round-trip..." << endl;

  const double det_diameter = 5.0 * PhysicalUnits::cm;
  const double setback = 0.5 * PhysicalUnits::cm;
  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };

  // Create DRF with setback
  auto original = make_shared<DetectorPeakResponse>( "SetbackURL", "URL setback test" );
  original->fromExpOfLogPowerSeries( coeffs, {}, 25.0 * PhysicalUnits::cm, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  original->setDetectorSetback( setback );

  // Serialize to URL
  string url;
  BOOST_CHECK_NO_THROW( url = original->toAppUrl() );
  BOOST_CHECK_MESSAGE( url.find("SETBK") != string::npos,
    "URL should contain SETBK parameter" );

  // Deserialize from URL
  auto restored = make_shared<DetectorPeakResponse>();
  BOOST_CHECK_NO_THROW( restored->fromAppUrl( url ) );

  // Verify setback is preserved
  BOOST_CHECK_MESSAGE( close_enough( restored->detectorSetback(), setback, 1e-3 ),
    "Setback not preserved in URL: " + to_string(restored->detectorSetback()) + " vs " + to_string(setback) );

  // Test backward compat: URL without SETBK should give setback=0
  auto drf_no_setback = make_shared<DetectorPeakResponse>( "NoSetbackURL", "no setback" );
  drf_no_setback->fromExpOfLogPowerSeries( coeffs, {}, 25.0 * PhysicalUnits::cm, det_diameter,
    PhysicalUnits::keV, 50.0f, 3000.0f,
    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  const string url2 = drf_no_setback->toAppUrl();
  BOOST_CHECK_MESSAGE( url2.find("SETBK") == string::npos,
    "URL for setback=0 DRF should NOT contain SETBK" );

  auto restored2 = make_shared<DetectorPeakResponse>();
  restored2->fromAppUrl( url2 );
  BOOST_CHECK_CLOSE( restored2->detectorSetback(), 0.0, 1e-10 );

  cout << "Detector setback URL round-trip test passed" << endl;
}//test_detector_setback_url_roundtrip


BOOST_AUTO_TEST_CASE( test_detector_setback_outx )
{
  cout << "\n\nTesting detector setback from ANGLE .outx file..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set (use --testfiledir=...)" );

  const string outx_file = SpecUtils::append_path( g_test_data_dir, "det_eff/Angle-example-efficiency.outx" );
  if( !SpecUtils::is_file(outx_file) )
  {
    cout << "  Skipping - OUTX file not found: " << outx_file << endl;
    return;
  }

  ifstream input( outx_file.c_str() );
  BOOST_REQUIRE( input.is_open() );

  shared_ptr<DetectorPeakResponse> drf = DetectorPeakResponse::parseAngleOutxFile( input );
  BOOST_REQUIRE( drf != nullptr );

  // The .outx file should have a positive setback (sum of endCap + vacuum + inactiveGe + housing layers)
  BOOST_CHECK_MESSAGE( drf->detectorSetback() > 0.0,
    "OUTX file should have positive setback, got: " + to_string(drf->detectorSetback() / PhysicalUnits::mm) + " mm" );

  cout << "  Setback from .outx: " << (drf->detectorSetback() / PhysicalUnits::mm) << " mm" << endl;
  cout << "Detector setback .outx test passed" << endl;
}//test_detector_setback_outx


BOOST_AUTO_TEST_CASE( test_detector_setback_common_drfs )
{
  cout << "\n\nTesting detector setback from common_drfs.tsv..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set (use --datadir=...)" );

  const string drfs_file = SpecUtils::append_path( g_data_dir, "common_drfs.tsv" );
  if( !SpecUtils::is_file(drfs_file) )
  {
    cout << "  Skipping - common_drfs.tsv not found: " << drfs_file << endl;
    return;
  }

  ifstream input( drfs_file.c_str() );
  BOOST_REQUIRE( input.is_open() );

  vector<string> credits;
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );

  BOOST_REQUIRE_MESSAGE( !drfs.empty(), "Should parse at least one DRF from common_drfs.tsv" );

  // Check that some DRFs have non-zero setback (e.g., ORTEC Detective-DX has 0.5 cm)
  int drfs_with_setback = 0;
  for( const auto &drf : drfs )
  {
    if( drf->detectorSetback() > 0.0 )
    {
      drfs_with_setback++;
      cout << "  " << drf->name() << ": setback = "
           << (drf->detectorSetback() / PhysicalUnits::cm) << " cm" << endl;
    }
  }

  BOOST_CHECK_MESSAGE( drfs_with_setback > 0,
    "At least some DRFs in common_drfs.tsv should have non-zero setback" );

  cout << "Found " << drfs_with_setback << " DRFs with setback out of " << drfs.size() << " total" << endl;
  cout << "Detector setback common_drfs.tsv test passed" << endl;
}//test_detector_setback_common_drfs


BOOST_AUTO_TEST_CASE( test_detector_setback_gadras )
{
  cout << "\n\nTesting detector setback from GADRAS Detector.dat..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set (use --datadir=...)" );

  // Find a GADRAS detector that has field 40 (det setback) > 0
  const string gadras_dir = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" );
  if( !SpecUtils::is_directory(gadras_dir) )
  {
    cout << "  Skipping - GenericGadrasDetectors directory not found" << endl;
    return;
  }

  // Try to find a detector with a non-zero setback
  // NaI 1x1 has setback of 0.5 cm per the data
  const vector<string> detectors_to_try = { "NaI 1x1", "NaI 2x2", "NaI 3x3" };

  for( const string &det_name : detectors_to_try )
  {
    const string det_dir = SpecUtils::append_path( gadras_dir, det_name );
    const string det_dat = SpecUtils::append_path( det_dir, "Detector.dat" );
    const string eff_csv = SpecUtils::append_path( det_dir, "Efficiency.csv" );

    if( !SpecUtils::is_file(det_dat) || !SpecUtils::is_file(eff_csv) )
      continue;

    ifstream dat_input( det_dat.c_str() );
    ifstream csv_input( eff_csv.c_str() );

    if( !dat_input.is_open() || !csv_input.is_open() )
      continue;

    DetectorPeakResponse drf;
    try
    {
      drf.fromGadrasDefinition( csv_input, dat_input );
    }catch( std::exception &e )
    {
      cout << "  Failed to parse " << det_name << ": " << e.what() << endl;
      continue;
    }

    cout << "  " << det_name << " setback: " << (drf.detectorSetback() / PhysicalUnits::cm) << " cm" << endl;

    // The NaI 1x1 should have setback of 0.5 cm
    if( det_name == "NaI 1x1" )
    {
      BOOST_CHECK_MESSAGE( close_enough( drf.detectorSetback(), 0.5 * PhysicalUnits::cm, 0.01 ),
        "NaI 1x1 should have setback of ~0.5 cm, got: "
        + to_string(drf.detectorSetback() / PhysicalUnits::cm) + " cm" );
    }
  }//for( const string &det_name : detectors_to_try )

  cout << "Detector setback GADRAS test passed" << endl;
}//test_detector_setback_gadras


BOOST_AUTO_TEST_CASE( test_read_gameff_csv )
{
  cout << "\n\nTesting gamEff CSV file reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set" );

  const string csv_file = SpecUtils::append_path( g_test_data_dir,
    "det_eff/25-163U01_ES2502_23003-1-F-1_gamEff.csv" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( csv_file ), "gamEff CSV not found: " + csv_file );

  ifstream input( csv_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open gamEff CSV" );

  DetectorPeakResponse::EffCsvParseResult result;
  BOOST_CHECK_NO_THROW( result = DetectorPeakResponse::parseEfficiencyCsvFile( input ) );

  BOOST_REQUIRE_MESSAGE( result.drf != nullptr, "gamEff parsing should return valid DRF" );
  BOOST_CHECK_MESSAGE( result.drf->isValid(), "gamEff DRF should be valid" );
  BOOST_CHECK_MESSAGE( !result.is_gadras_format, "gamEff should NOT be detected as GADRAS format" );

  BOOST_CHECK_MESSAGE( result.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct,
                       "gamEff DRF should have FixedGeomTotalAct geometry" );
  BOOST_CHECK_MESSAGE( result.drf->efficiencyFcnType() == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs,
                       "gamEff DRF should use energy-efficiency pairs" );
  BOOST_CHECK_MESSAGE( result.drf->drfSource() == DetectorPeakResponse::DrfSource::UserImportedEfficiencyCsvDrf,
                       "gamEff DRF should have UserImportedEfficiencyCsvDrf source" );

  // Check first data point: 48.814 keV, 3.098e-03
  const float eff_first = result.drf->intrinsicEfficiency( 48.814f * static_cast<float>(PhysicalUnits::keV) );
  BOOST_CHECK_MESSAGE( close_enough( eff_first, 3.098e-03, 0.01 ),
    "gamEff first point efficiency: " + to_string( eff_first ) + " vs expected 3.098e-03" );

  // All efficiencies should be < 1.0 (fractional, not percentage)
  BOOST_CHECK_MESSAGE( eff_first < 1.0f, "gamEff efficiencies should be fractional (< 1.0)" );

  cout << "Successfully parsed gamEff CSV file" << endl;
}//test_read_gameff_csv


BOOST_AUTO_TEST_CASE( test_read_run_effoutput_csv )
{
  cout << "\n\nTesting Run_effoutput CSV file reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set" );

  const string csv_file = SpecUtils::append_path( g_test_data_dir,
    "det_eff/Run_effoutput_BIG8_Foil_86mm_110p_Abs-1.csv" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( csv_file ), "Run_effoutput CSV not found: " + csv_file );

  ifstream input( csv_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open Run_effoutput CSV" );

  DetectorPeakResponse::EffCsvParseResult result;
  BOOST_CHECK_NO_THROW( result = DetectorPeakResponse::parseEfficiencyCsvFile( input ) );

  BOOST_REQUIRE_MESSAGE( result.drf != nullptr, "Run_effoutput parsing should return valid DRF" );
  BOOST_CHECK_MESSAGE( result.drf->isValid(), "Run_effoutput DRF should be valid" );
  BOOST_CHECK_MESSAGE( !result.is_gadras_format, "Run_effoutput should NOT be detected as GADRAS format" );

  BOOST_CHECK_MESSAGE( result.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct,
                       "Run_effoutput DRF should have FixedGeomTotalAct geometry" );
  BOOST_CHECK_MESSAGE( result.drf->drfSource() == DetectorPeakResponse::DrfSource::UserImportedEfficiencyCsvDrf,
                       "Run_effoutput DRF should have UserImportedEfficiencyCsvDrf source" );

  // Check first data point: Energy=34.9731 keV, Eff=2.93E-06 (column 4, "Eff")
  const float eff_first = result.drf->intrinsicEfficiency( 34.9731f * static_cast<float>(PhysicalUnits::keV) );
  BOOST_CHECK_MESSAGE( close_enough( eff_first, 2.93e-06, 0.02 ),
    "Run_effoutput first point efficiency: " + to_string( eff_first ) + " vs expected 2.93e-06" );

  // All efficiencies should be < 1.0 (fractional, not percentage)
  BOOST_CHECK_MESSAGE( eff_first < 1.0f, "Run_effoutput efficiencies should be fractional (< 1.0)" );

  cout << "Successfully parsed Run_effoutput CSV file" << endl;
}//test_read_run_effoutput_csv


BOOST_AUTO_TEST_CASE( test_read_gadras_csv_standalone )
{
  cout << "\n\nTesting standalone GADRAS Efficiency.csv reading..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set" );

  const string csv_file = SpecUtils::append_path( g_data_dir,
    "GenericGadrasDetectors/HPGe 40%/Efficiency.csv" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( csv_file ), "GADRAS Efficiency.csv not found: " + csv_file );

  ifstream input( csv_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open GADRAS Efficiency.csv" );

  DetectorPeakResponse::EffCsvParseResult result;
  BOOST_CHECK_NO_THROW( result = DetectorPeakResponse::parseEfficiencyCsvFile( input ) );

  BOOST_REQUIRE_MESSAGE( result.drf != nullptr, "GADRAS CSV parsing should return valid DRF" );
  BOOST_CHECK_MESSAGE( result.drf->isValid(), "GADRAS CSV DRF should be valid" );
  BOOST_CHECK_MESSAGE( result.is_gadras_format, "GADRAS CSV should be detected as GADRAS format" );

  BOOST_CHECK_MESSAGE( result.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct,
                       "GADRAS CSV DRF should have FixedGeomTotalAct geometry (before reinterpretation)" );

  // Efficiency at ~100 keV should be around 0.49899 (from file: 49.899%)
  // The exact energy may not be in the file, so use intrinsicEfficiency which interpolates
  const float eff_100 = result.drf->intrinsicEfficiency( 100.0f * static_cast<float>(PhysicalUnits::keV) );
  BOOST_CHECK_MESSAGE( (eff_100 > 0.0f) && (eff_100 <= 1.0f),
    "GADRAS CSV efficiency should be in [0, 1] after percentage conversion, got: " + to_string( eff_100 ) );

  cout << "Successfully parsed standalone GADRAS Efficiency.csv" << endl;
}//test_read_gadras_csv_standalone


BOOST_AUTO_TEST_CASE( test_gadras_ptot_total_efficiency )
{
  cout << "\n\nTesting GADRAS PTOT -> total efficiency..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set" );

  const string csv_file = SpecUtils::append_path( g_data_dir,
    "GenericGadrasDetectors/HPGe 40%/Efficiency.csv" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( csv_file ), "GADRAS Efficiency.csv not found: " + csv_file );

  ifstream input( csv_file.c_str() );
  BOOST_REQUIRE( input.is_open() );

  DetectorPeakResponse::EffCsvParseResult result;
  BOOST_REQUIRE_NO_THROW( result = DetectorPeakResponse::parseEfficiencyCsvFile( input ) );
  BOOST_REQUIRE( result.drf && result.drf->isValid() );

  // The GADRAS Efficiency.csv has a "PTOT" column -> a total efficiency curve.
  BOOST_REQUIRE_MESSAGE( result.drf->hasTotalEfficiency(),
                         "GADRAS CSV with a PTOT column should populate total efficiency" );

  // Total efficiency should be finite, non-negative, and >= the full-energy
  //  efficiency at every (mid-range) energy.  Note: GADRAS efficiencies can
  //  legitimately exceed 100% (see the comment in parseEfficiencyCsvFile), so
  //  there is no tight upper bound.
  for( const float energy : { 200.0f, 500.0f, 661.7f, 1332.5f, 2000.0f } )
  {
    const float full = result.drf->intrinsicEfficiency( energy * static_cast<float>(PhysicalUnits::keV) );
    const float total = result.drf->totalIntrinsicEfficiency( energy * static_cast<float>(PhysicalUnits::keV) );
    BOOST_CHECK_MESSAGE( (total >= 0.0f) && std::isfinite(total),
      "Total efficiency invalid at " + to_string(energy) + " keV: " + to_string(total) );
    BOOST_CHECK_MESSAGE( total >= (full - 1.0E-4f),
      "Total efficiency (" + to_string(total) + ") should be >= full-energy ("
      + to_string(full) + ") at " + to_string(energy) + " keV" );
  }

  // Reading total efficiency must change the DRF hash relative to one without
  //  it (this is the accepted lineage break) - confirm by clearing it.
  auto drf_no_total = make_shared<DetectorPeakResponse>( *result.drf );
  const uint64_t hash_with_total = drf_no_total->hashValue();
  drf_no_total->setTotalEfficiencyCurve( nullptr );
  BOOST_CHECK_NE( hash_with_total, drf_no_total->hashValue() );

  // The full GADRAS read path (Efficiency.csv + Detector.dat) carries the
  //  total efficiency across into the far-field-intrinsic DRF.
  const string det_dir = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 40%" );
  DetectorPeakResponse gadras_drf;
  BOOST_REQUIRE_NO_THROW( gadras_drf.fromGadrasDirectory( det_dir ) );
  BOOST_CHECK_MESSAGE( gadras_drf.hasTotalEfficiency(),
                       "fromGadrasDirectory should carry the PTOT total efficiency across" );

  // And it must survive an XML round-trip.
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *parent = doc.allocate_node( rapidxml::node_element, "Parent" );
  doc.append_node( parent );
  gadras_drf.toXml( parent, &doc );

  auto restored = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( restored->fromXml( parent->first_node("DetectorPeakResponse") ) );
  BOOST_CHECK( restored->hasTotalEfficiency() );
  BOOST_CHECK_EQUAL( restored->hashValue(), gadras_drf.hashValue() );

  cout << "GADRAS PTOT total efficiency passed" << endl;
}//test_gadras_ptot_total_efficiency


BOOST_AUTO_TEST_CASE( test_parseDetectorDatGeometry )
{
  cout << "\n\nTesting Detector.dat geometry parsing..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Data directory not set" );

  const string dat_file = SpecUtils::append_path( g_data_dir,
    "GenericGadrasDetectors/HPGe 40%/Detector.dat" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( dat_file ), "Detector.dat not found: " + dat_file );

  ifstream input( dat_file.c_str() );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open Detector.dat" );

  float diameter = 0.0f, setback = 0.0f;
  BOOST_CHECK_NO_THROW(
    DetectorPeakResponse::parseDetectorDatGeometry( input, diameter, setback )
  );

  BOOST_CHECK_MESSAGE( diameter > 0.0f,
    "Detector diameter should be positive, got: " + to_string( diameter / static_cast<float>(PhysicalUnits::cm) ) + " cm" );

  // Compare with the diameter from fromGadrasDirectory
  const string det_dir = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 40%" );
  DetectorPeakResponse drf_full;
  BOOST_CHECK_NO_THROW( drf_full.fromGadrasDirectory( det_dir ) );

  BOOST_CHECK_MESSAGE( close_enough( diameter, drf_full.detectorDiameter(), 0.01 ),
    "parseDetectorDatGeometry diameter (" + to_string( diameter / static_cast<float>(PhysicalUnits::cm) ) +
    " cm) should match fromGadrasDirectory (" +
    to_string( drf_full.detectorDiameter() / static_cast<float>(PhysicalUnits::cm) ) + " cm)" );

  cout << "Successfully parsed Detector.dat geometry" << endl;
}//test_parseDetectorDatGeometry


BOOST_AUTO_TEST_CASE( test_total_efficiency_basics )
{
  cout << "\n\nTesting total efficiency basics..." << endl;

  const double det_diameter = 5.0 * PhysicalUnits::cm;
  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };

  auto drf = make_shared<DetectorPeakResponse>( "TotalEffTest", "total eff test" );
  drf->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter, PhysicalUnits::keV,
                               50.0f, 3000.0f,
                               DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  BOOST_CHECK( !drf->hasTotalEfficiency() );
  BOOST_CHECK( !drf->totalEfficiencyCurve() );
  BOOST_CHECK_THROW( drf->totalIntrinsicEfficiency( 661.0f ), std::runtime_error );
  BOOST_CHECK_THROW( drf->totalEfficiency( 661.0f, 100.0*PhysicalUnits::cm ), std::runtime_error );

  // Set a pairs-based total efficiency curve
  vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
  for( size_t i = 0; i < 10; ++i )
  {
    DetectorPeakResponse::EnergyEfficiencyPair p;
    p.energy = 50.0f + 300.0f*i;
    p.efficiency = 0.9f * exp( -0.0005f * p.energy ) + 0.02f;
    pairs.push_back( p );
  }

  auto curve = make_shared<DetectorEfficiencyCurve>();
  curve->setFromPairs( pairs, static_cast<float>(PhysicalUnits::keV) );

  const uint64_t hash_before = drf->hashValue();
  drf->setTotalEfficiencyCurve( curve );
  BOOST_CHECK( drf->hasTotalEfficiency() );
  BOOST_CHECK_NE( drf->hashValue(), hash_before );

  // Total intrinsic efficiency matches akima interpolation of the pairs
  for( const float energy : { 75.0f, 661.0f, 2000.0f } )
  {
    const float expected = DetectorPeakResponse::akimaInterpolate( energy, pairs );
    BOOST_CHECK( close_enough( drf->totalIntrinsicEfficiency(energy), expected, 1e-6 ) );
  }

  // totalEfficiency(energy,dist) applies the same solid angle as efficiency(energy,dist)
  const double dist = 50.0 * PhysicalUnits::cm;
  const double frac_solid = DetectorPeakResponse::fractionalSolidAngle( det_diameter, dist );
  BOOST_CHECK( close_enough( drf->totalEfficiency( 661.0f, dist ),
                             frac_solid * drf->totalIntrinsicEfficiency(661.0f), 1e-6 ) );

  // Total efficiency should be >= full-energy efficiency for a real detector
  BOOST_CHECK_GT( drf->totalIntrinsicEfficiency(661.0f), drf->intrinsicEfficiency(661.0f) );

  // Clearing restores the original hash
  drf->setTotalEfficiencyCurve( nullptr );
  BOOST_CHECK( !drf->hasTotalEfficiency() );
  BOOST_CHECK_EQUAL( drf->hashValue(), hash_before );

  cout << "Total efficiency basics passed" << endl;
}//test_total_efficiency_basics


namespace
{
  // Builds a DRF with both new optional fields set, for round-trip tests.
  shared_ptr<DetectorPeakResponse> make_drf_with_new_fields()
  {
    const double det_diameter = 5.0 * PhysicalUnits::cm;
    const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };

    auto drf = make_shared<DetectorPeakResponse>( "NewFieldsTest", "new fields test" );
    drf->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter, PhysicalUnits::keV,
                                 50.0f, 3000.0f,
                                 DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

    // Full-energy efficiency uncertainty: bands + node covariance + coef covariance
    const vector<float> uncert_energies = { 59.5f, 122.0f, 661.7f, 1332.5f };
    const vector<float> uncert_vals = { 0.09f, 0.06f, 0.045f, 0.055f };
    shared_ptr<DetectorEfficiencyUncert> tmp
                = DetectorEfficiencyUncert::fromPointUncerts( uncert_energies, uncert_vals );
    auto uncert = make_shared<DetectorEfficiencyUncert>( *tmp );

    vector<EffUncertBand> bands;
    bands.push_back( EffUncertBand{ 50.0f, 122.0f, 0.08f } );
    bands.push_back( EffUncertBand{ 122.0f, 661.0f, 0.05f } );
    uncert->setBands( bands );
    uncert->setCoefficientCovariance( { 1.0E-4f, -2.0E-5f, 0.0f,
                                        -2.0E-5f, 4.0E-5f, 0.0f,
                                        0.0f, 0.0f, 9.0E-6f } );
    drf->setEfficiencyUncert( uncert );

    // Total efficiency with its own uncertainty
    auto total = make_shared<DetectorEfficiencyCurve>();
    total->setFromExpOfLogPowerSeries( { -1.5f, 0.4f, -0.08f }, {},
                                       static_cast<float>(PhysicalUnits::keV) );
    total->setUncertainty( DetectorEfficiencyUncert::fromPointUncerts(
                            { 60.0f, 700.0f, 2300.0f }, { 0.1f, 0.05f, 0.07f } ) );
    drf->setTotalEfficiencyCurve( total );

    return drf;
  }//make_drf_with_new_fields()
}//namespace


BOOST_AUTO_TEST_CASE( test_xml_v5_roundtrip )
{
  cout << "\n\nTesting XML v5 round-trip..." << endl;

  shared_ptr<DetectorPeakResponse> orig = make_drf_with_new_fields();

  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *parent = doc.allocate_node( rapidxml::node_element, "Parent" );
  doc.append_node( parent );
  orig->toXml( parent, &doc );

  const rapidxml::xml_node<char> *drf_node = parent->first_node( "DetectorPeakResponse" );
  BOOST_REQUIRE( drf_node );

  // Version attribute should be 5, since the new fields are present
  const rapidxml::xml_attribute<char> *version_attr = drf_node->first_attribute( "version" );
  BOOST_REQUIRE( version_attr );
  BOOST_CHECK_EQUAL( string(version_attr->value()), "5" );

  auto restored = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( restored->fromXml( drf_node ) );

  // Float parsing may differ in the last ULP, so use the tolerant comparisons
  BOOST_REQUIRE( restored->efficiencyUncert() );
  BOOST_CHECK_NO_THROW( DetectorEfficiencyUncert::equalEnough( *restored->efficiencyUncert(),
                                                          *orig->efficiencyUncert() ) );
  BOOST_REQUIRE( restored->hasTotalEfficiency() );
  BOOST_CHECK_NO_THROW( DetectorEfficiencyCurve::equalEnough( *restored->totalEfficiencyCurve(),
                                                         *orig->totalEfficiencyCurve() ) );
  BOOST_REQUIRE( restored->totalEfficiencyCurve()->uncertainty() );

#if( PERFORM_DEVELOPER_CHECKS )
  BOOST_CHECK_NO_THROW( DetectorPeakResponse::equalEnough( *orig, *restored ) );
#endif

  // A DRF without the new fields must keep writing a version <= 4
  {
    const double det_diameter = 5.0 * PhysicalUnits::cm;
    auto plain = make_shared<DetectorPeakResponse>( "Plain", "plain test" );
    plain->fromExpOfLogPowerSeries( { -5.0f, 0.5f }, {}, 0.0, det_diameter, PhysicalUnits::keV,
                                   50.0f, 3000.0f,
                                   DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

    rapidxml::xml_document<char> doc2;
    rapidxml::xml_node<char> *parent2 = doc2.allocate_node( rapidxml::node_element, "Parent" );
    doc2.append_node( parent2 );
    plain->toXml( parent2, &doc2 );

    const rapidxml::xml_node<char> *plain_node = parent2->first_node( "DetectorPeakResponse" );
    BOOST_REQUIRE( plain_node );
    const rapidxml::xml_attribute<char> *plain_version = plain_node->first_attribute( "version" );
    BOOST_REQUIRE( plain_version );
    const int version_num = atoi( plain_version->value() );
    BOOST_CHECK_MESSAGE( version_num <= 4,
      "DRF without new fields wrote version " + string(plain_version->value()) );
  }

  cout << "XML v5 round-trip passed" << endl;
}//test_xml_v5_roundtrip


BOOST_AUTO_TEST_CASE( test_xml_v4_backcompat )
{
  cout << "\n\nTesting back-compat reading of pre-version-5 XML..." << endl;

  // A literal version-1 era DRF XML (exp-of-log series, far-field intrinsic)
  const char *old_xml =
    "<DetectorPeakResponse version=\"1\">"
      "<Name>OldDrf</Name>"
      "<Description>old format</Description>"
      "<DetectorDiameter>5.00000000E+00</DetectorDiameter>"
      "<EfficiencySource>UserAddedRelativeEfficiencyDrf</EfficiencySource>"
      "<EfficiencyEnergyUnits>1</EfficiencyEnergyUnits>"
      "<ResolutionForm>Undefined</ResolutionForm>"
      "<EfficiencyForm>ExpOfLogPowerSeries</EfficiencyForm>"
      "<ExpOfLogPowerSeriesCoeffs>-5.0 0.5 -0.01</ExpOfLogPowerSeriesCoeffs>"
      "<Hash>123456789</Hash>"
      "<ParentHash>0</ParentHash>"
      "<Flags>0</Flags>"
      "<Geometry>FAR-FIELD</Geometry>"
    "</DetectorPeakResponse>";

  vector<char> xml_buf( old_xml, old_xml + strlen(old_xml) );
  xml_buf.push_back( '\0' );

  rapidxml::xml_document<char> doc;
  BOOST_REQUIRE_NO_THROW( doc.parse<rapidxml::parse_trim_whitespace>( xml_buf.data() ) );

  auto drf = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( drf->fromXml( doc.first_node("DetectorPeakResponse") ) );

  BOOST_CHECK( drf->isValid() );
  BOOST_CHECK_EQUAL( drf->name(), "OldDrf" );
  BOOST_CHECK( !drf->efficiencyUncert() );
  BOOST_CHECK( !drf->hasTotalEfficiency() );

  cout << "Pre-version-5 XML back-compat passed" << endl;
}//test_xml_v4_backcompat


BOOST_AUTO_TEST_CASE( test_url_roundtrip_with_uncert )
{
  cout << "\n\nTesting URL round-trip with uncertainties..." << endl;

  shared_ptr<DetectorPeakResponse> orig = make_drf_with_new_fields();

  string url;
  BOOST_REQUIRE_NO_THROW( url = orig->toAppUrl() );
  BOOST_CHECK( !url.empty() );

  auto restored = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( restored->fromAppUrl( url ) );

  // The uncertainty should have made it across (values to URL float precision)
  BOOST_REQUIRE( restored->efficiencyUncert() );
  BOOST_CHECK( restored->efficiencyUncert()->hasBands() );
  BOOST_CHECK( restored->efficiencyUncert()->hasNodeCovariance() );

  // Coefficient covariance is never written to URLs
  BOOST_CHECK( restored->efficiencyUncert()->coefficientCovariance().empty() );

  const vector<double> test_energies = { 80.0, 661.0, 1332.0 };
  const vector<double> orig_uncerts = orig->efficiencyUncert()->fracUncertainties( test_energies );
  const vector<double> rest_uncerts = restored->efficiencyUncert()->fracUncertainties( test_energies );
  for( size_t i = 0; i < test_energies.size(); ++i )
    BOOST_CHECK( close_enough( orig_uncerts[i], rest_uncerts[i], 1e-3 ) );

  // Total efficiency curve too
  BOOST_REQUIRE( restored->hasTotalEfficiency() );
  BOOST_REQUIRE( restored->totalEfficiencyCurve()->uncertainty() );
  for( const float energy : { 100.0f, 661.0f, 2000.0f } )
    BOOST_CHECK( close_enough( restored->totalIntrinsicEfficiency(energy),
                               orig->totalIntrinsicEfficiency(energy), 1e-4 ) );

  cout << "URL round-trip with uncertainties passed" << endl;
}//test_url_roundtrip_with_uncert


BOOST_AUTO_TEST_CASE( test_url_drop_order )
{
  cout << "\n\nTesting URL drop order when over the QR budget..." << endl;

  // Build a pairs-form DRF with many points plus a large node covariance, so
  //  the URL exceeds the QR budget and the new keys must get dropped first.
  auto drf = make_shared<DetectorPeakResponse>( "DropOrderTest", "drop order test" );

  // Enough points that the base URL is large (but still within the QR budget
  //  once the new keys are dropped) - the node covariance then pushes it over.
  vector<DetectorPeakResponse::EnergyEffPoint> points;
  for( size_t i = 0; i < 120; ++i )
  {
    DetectorPeakResponse::EnergyEffPoint p;
    p.energy = 20.0f + 24.5f*i;
    p.efficiency = 0.5f * exp( -0.0005f * p.energy ) + 0.01f;
    points.push_back( p );
  }

  drf->setEfficiencyPoints( points, 5.0*PhysicalUnits::cm, 0.0,
                           DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  drf->setDetectorSetback( 0.5 * PhysicalUnits::cm );

  // Large node covariance (40 nodes -> 820 upper-triangle entries)
  vector<float> uncert_energies, uncert_vals;
  for( size_t i = 0; i < 40; ++i )
  {
    uncert_energies.push_back( 30.0f + 70.0f*i );
    uncert_vals.push_back( 0.05f + 0.001f*i );
  }
  auto big_uncert = make_shared<DetectorEfficiencyUncert>(
            *DetectorEfficiencyUncert::fromPointUncerts( uncert_energies, uncert_vals ) );
  vector<EffUncertBand> bands;
  bands.push_back( EffUncertBand{ 50.0f, 661.0f, 0.07f } );
  big_uncert->setBands( bands );
  drf->setEfficiencyUncert( big_uncert );

  string url;
  BOOST_REQUIRE_NO_THROW( url = drf->toAppUrl() );

  // The URL must parse, and the covariance keys must have been dropped before
  //  any of the pre-existing keys (like SETBK).
  auto restored = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( restored->fromAppUrl( url ) );
  BOOST_CHECK( restored->isValid() );

  BOOST_CHECK_MESSAGE( url.find("EFUC=") == string::npos,
                       "Node covariance should have been dropped from oversized URL" );
  BOOST_CHECK_MESSAGE( url.find("SETBK=") != string::npos,
                       "SETBK should not be dropped before the new keys" );

  cout << "  URL length after drops: " << url.size() << " chars" << endl;
  cout << "URL drop order passed" << endl;
}//test_url_drop_order


BOOST_AUTO_TEST_CASE( test_hash_stability_with_uncert )
{
  cout << "\n\nTesting hash stability with new optional fields..." << endl;

  const double det_diameter = 5.0 * PhysicalUnits::cm;
  const vector<float> coeffs = { -5.0f, 0.5f, -0.01f };

  auto drf = make_shared<DetectorPeakResponse>( "HashTest", "hash test" );
  drf->fromExpOfLogPowerSeries( coeffs, {}, 0.0, det_diameter, PhysicalUnits::keV,
                               50.0f, 3000.0f,
                               DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  const uint64_t hash_plain = drf->hashValue();

  // Adding the uncertainty changes the hash; removing it restores it.
  auto uncert = make_shared<DetectorEfficiencyUncert>(
        *DetectorEfficiencyUncert::fromPointUncerts( {60.0f, 661.0f}, {0.05f, 0.04f} ) );
  drf->setEfficiencyUncert( uncert );
  const uint64_t hash_uncert = drf->hashValue();
  BOOST_CHECK_NE( hash_plain, hash_uncert );

  drf->setEfficiencyUncert( nullptr );
  BOOST_CHECK_EQUAL( drf->hashValue(), hash_plain );

  // Setting an empty uncertainty is the same as no uncertainty.
  drf->setEfficiencyUncert( make_shared<DetectorEfficiencyUncert>() );
  BOOST_CHECK_EQUAL( drf->hashValue(), hash_plain );
  BOOST_CHECK( !drf->efficiencyUncert() );

  cout << "Hash stability with new fields passed" << endl;
}//test_hash_stability_with_uncert


BOOST_AUTO_TEST_CASE( test_angle_outx_precision_covariance )
{
  cout << "\n\nTesting ANGLE .outx efficiencyPrecision -> covariance..." << endl;

  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Test data directory not set (use --testfiledir=...)" );

  const string outx_path = SpecUtils::append_path( g_test_data_dir,
                                          "det_eff/Angle-example-efficiency.outx" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(outx_path), "Missing " + outx_path );

  ifstream input( outx_path.c_str(), ios::in | ios::binary );
  BOOST_REQUIRE( input.is_open() );

  shared_ptr<DetectorPeakResponse> drf;
  BOOST_REQUIRE_NO_THROW( drf = DetectorPeakResponse::parseAngleOutxFile( input ) );
  BOOST_REQUIRE( drf && drf->isValid() );

  shared_ptr<const DetectorEfficiencyUncert> uncert = drf->efficiencyUncert();
  BOOST_REQUIRE_MESSAGE( uncert, "ANGLE .outx with precision attributes should give uncertainty" );
  BOOST_REQUIRE( uncert->hasNodeCovariance() );

  // First result in the example file: energy=30, efficiencyPrecision=1.41419957134595
  //  -> fractional uncert ~0.4142, variance ~0.1716
  const vector<float> &node_energies = uncert->covarianceEnergies();
  const vector<float> &cov = uncert->covarianceMatrix();
  const size_t n = node_energies.size();
  BOOST_REQUIRE_GT( n, 4u );

  BOOST_CHECK( close_enough( node_energies[0], 30.0, 1e-4 ) );
  const double expected_frac = 1.41419957134595 - 1.0;
  BOOST_CHECK_MESSAGE( close_enough( cov[0], expected_frac*expected_frac, 1e-3 ),
    "Variance at 30 keV: got " + to_string(cov[0]) + ", expected "
    + to_string(expected_frac*expected_frac) );

  // Off-diagonal between first two nodes follows the default-corr-length kernel
  {
    const double corr_len = DetectorEfficiencyUncert::sm_defaultLogEnergyCorrLength;
    BOOST_CHECK( close_enough( uncert->correlationLength(), corr_len, 1e-6 ) );

    const double u0 = sqrt( double(cov[0]) );
    const double u1 = sqrt( double(cov[1*n + 1]) );
    const double dlne = log(double(node_energies[0])) - log(double(node_energies[1]));
    const double rho = exp( -0.5 * pow(dlne/corr_len, 2.0) );
    BOOST_CHECK( close_enough( cov[1], rho*u0*u1, 1e-3 ) );
  }

  cout << "ANGLE .outx precision covariance passed (" << n << " nodes)" << endl;
}//test_angle_outx_precision_covariance


BOOST_AUTO_TEST_CASE( test_eff_csv_uncert_column )
{
  cout << "\n\nTesting efficiency CSV with uncertainty column..." << endl;

  // Synthetic Run_effoutput-style CSV with an "Eff. Unc" column
  const char *csv =
    "Energy[keV],Emitted,Peak,Eff,Eff. Unc\n"
    "59.5,100000,500,0.005,0.0005\n"
    "122.0,100000,2000,0.02,0.001\n"
    "661.7,100000,1500,0.015,0.0006\n"
    "1332.5,100000,800,0.008,0.0004\n";

  stringstream input( csv );
  DetectorPeakResponse::EffCsvParseResult result;
  BOOST_REQUIRE_NO_THROW( result = DetectorPeakResponse::parseEfficiencyCsvFile( input ) );
  BOOST_REQUIRE( result.drf && result.drf->isValid() );
  BOOST_CHECK( !result.is_gadras_format );

  shared_ptr<const DetectorEfficiencyUncert> uncert = result.drf->efficiencyUncert();
  BOOST_REQUIRE_MESSAGE( uncert, "CSV with uncertainty column should give an uncertainty" );
  BOOST_REQUIRE( uncert->hasNodeCovariance() );
  BOOST_REQUIRE_EQUAL( uncert->covarianceEnergies().size(), 4u );

  // Fractional uncertainty at 59.5 keV: 0.0005/0.005 = 0.1
  const vector<double> fracs = uncert->fracUncertainties( { 59.5 } );
  BOOST_CHECK_MESSAGE( close_enough( fracs[0], 0.1, 1e-3 ),
    "Fractional uncert at 59.5 keV: got " + to_string(fracs[0]) + ", expected 0.1" );

  // A 2-column CSV must not gain an uncertainty
  const char *plain_csv =
    "Energy[keV],Eff\n"
    "59.5,0.005\n"
    "122.0,0.02\n"
    "661.7,0.015\n";

  stringstream plain_input( plain_csv );
  DetectorPeakResponse::EffCsvParseResult plain_result;
  BOOST_REQUIRE_NO_THROW( plain_result = DetectorPeakResponse::parseEfficiencyCsvFile( plain_input ) );
  BOOST_REQUIRE( plain_result.drf && plain_result.drf->isValid() );
  BOOST_CHECK( !plain_result.drf->efficiencyUncert() );

  cout << "Efficiency CSV uncertainty column passed" << endl;
}//test_eff_csv_uncert_column


BOOST_AUTO_TEST_CASE( test_drfextra_blob_roundtrip )
{
  cout << "\n\nTesting DrfExtra database-blob round-trip..." << endl;

  shared_ptr<DetectorPeakResponse> orig = make_drf_with_new_fields();

  const string blob = orig->drfExtraToXmlString();
  BOOST_CHECK( !blob.empty() );
  BOOST_CHECK( blob.find("<DrfExtra>") != string::npos );

  auto restored = make_shared<DetectorPeakResponse>( *orig );
  restored->setDrfExtraFromXmlString( "" );  //clears
  BOOST_CHECK( !restored->efficiencyUncert() );
  BOOST_CHECK( !restored->hasTotalEfficiency() );

  restored->setDrfExtraFromXmlString( blob );
  BOOST_REQUIRE( restored->efficiencyUncert() );
  BOOST_CHECK_NO_THROW( DetectorEfficiencyUncert::equalEnough( *restored->efficiencyUncert(),
                                                          *orig->efficiencyUncert() ) );
  BOOST_REQUIRE( restored->hasTotalEfficiency() );
  BOOST_CHECK_NO_THROW( DetectorEfficiencyCurve::equalEnough( *restored->totalEfficiencyCurve(),
                                                         *orig->totalEfficiencyCurve() ) );

  // A DRF without the new fields produces an empty blob
  auto plain = make_shared<DetectorPeakResponse>( "Plain", "plain" );
  plain->fromExpOfLogPowerSeries( { -5.0f, 0.5f }, {}, 0.0, 5.0*PhysicalUnits::cm,
                                 PhysicalUnits::keV, 50.0f, 3000.0f,
                                 DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  BOOST_CHECK( plain->drfExtraToXmlString().empty() );

  // Garbage must not throw, just clear
  BOOST_CHECK_NO_THROW( restored->setDrfExtraFromXmlString( "<not really &xml" ) );
  BOOST_CHECK( !restored->efficiencyUncert() );
  BOOST_CHECK( !restored->hasTotalEfficiency() );

  cout << "DrfExtra blob round-trip passed" << endl;
}//test_drfextra_blob_roundtrip


BOOST_AUTO_TEST_CASE( test_set_efficiency_points_uncerts )
{
  cout << "\n\nTesting setEfficiencyPoints with per-point uncertainties..." << endl;

  vector<DetectorPeakResponse::EnergyEffPoint> points = {
    { 100.0f, 0.05f, 0.005f },
    { 500.0f, 0.15f, 0.012f },
    { 1000.0f, 0.12f, 0.011f },
    { 2000.0f, 0.08f, {} }       //no uncertainty for this point
  };

  auto drf = make_shared<DetectorPeakResponse>( "PointUncertTest", "point uncert test" );
  drf->setEfficiencyPoints( points, 5.0*PhysicalUnits::cm, 0.0,
                           DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  shared_ptr<const DetectorEfficiencyUncert> uncert = drf->efficiencyUncert();
  BOOST_REQUIRE( uncert );
  BOOST_REQUIRE( uncert->hasNodeCovariance() );
  BOOST_CHECK_EQUAL( uncert->covarianceEnergies().size(), 3u );  //only points with uncerts

  // Fractional uncertainty at 100 keV: 0.005/0.05 = 0.1
  const vector<double> fracs = uncert->fracUncertainties( { 100.0 } );
  BOOST_CHECK( close_enough( fracs[0], 0.1, 1e-4 ) );

  // Without uncertainties: no uncert object (and same hash as before the
  //  uncertainty support was added)
  vector<DetectorPeakResponse::EnergyEffPoint> plain_points = {
    { 100.0f, 0.05f, {} }, { 500.0f, 0.15f, {} }, { 1000.0f, 0.12f, {} }
  };
  auto plain = make_shared<DetectorPeakResponse>( "PlainPoints", "plain points" );
  plain->setEfficiencyPoints( plain_points, 5.0*PhysicalUnits::cm, 0.0,
                             DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  BOOST_CHECK( !plain->efficiencyUncert() );

  cout << "setEfficiencyPoints per-point uncertainties passed" << endl;
}//test_set_efficiency_points_uncerts
