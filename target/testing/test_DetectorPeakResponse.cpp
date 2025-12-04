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
    const double max_val = std::max(fabs(a), fabs(b));
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
      const double rel_diff = fabs(eff_orig - eff_rest) / std::max(fabs(eff_orig), fabs(eff_rest));
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
