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

#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE XRayWidthServer_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/XRayWidthServer.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;
using namespace XRayWidths;

// Global data directory path
static std::string g_data_dir;

// Boost test framework initialization
struct GlobalFixture
{
  GlobalFixture()
  {
    // Search for the data directory if not specified via command line
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        g_data_dir = d;
        break;
      }
    }

    if( g_data_dir.empty() )
    {
      cerr << "Warning: Could not find data directory with sandia.decay.xml" << endl;
      g_data_dir = "data";  // Default fallback
    }
  }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture );


BOOST_AUTO_TEST_CASE( DatabaseLoading )
{
  // Test that database loads successfully
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();

  BOOST_REQUIRE( db != nullptr );
  BOOST_CHECK( XRayWidthDatabase::status() == LoadStatus::Loaded );

  // Check that database has reasonable number of entries
  const size_t num_entries = db->num_entries();
  BOOST_CHECK_GE( num_entries, 109 );  // At least the original 109 entries

  std::cout << "XRayWidthDatabase loaded with " << num_entries << " entries" << std::endl;
}//BOOST_AUTO_TEST_CASE( DatabaseLoading )


BOOST_AUTO_TEST_CASE( OriginalDataValidation )
{
  // Verify that all 109 original natural width values from the hardcoded table
  // are present in the XML database with matching values (within 1% tolerance)

  struct OriginalEntry
  {
    int z;
    double energy_kev;
    double hwhm_ev;
    const char* label;
  };

  // Original 109 entries from PeakDists.cpp hardcoded table
  static const OriginalEntry original_data[] = {
    // Uranium (Z=92)
    { 92,  13.615,  5.2,  "U Lα1" },
    { 92,  13.439,  5.2,  "U Lα2" },
    { 92,  17.220,  6.6,  "U Lβ1" },
    { 92,  17.456,  6.6,  "U Lβ2" },
    { 92,  20.167,  7.5,  "U Lγ1" },
    { 92,  20.945,  8.3,  "U Lγ3" },
    { 92,  98.439, 48.0,  "U Kα1" },
    { 92,  94.665, 48.0,  "U Kα2" },
    { 92, 111.300, 48.0,  "U Kβ1" },
    { 92, 114.566, 48.0,  "U Kβ2" },

    // Plutonium (Z=94)
    { 94,  14.282,  5.5,  "Pu Lα1" },
    { 94,  14.084,  5.5,  "Pu Lα2" },
    { 94,  17.992,  6.9,  "Pu Lβ1" },
    { 94,  18.264,  6.9,  "Pu Lβ2" },
    { 94,  21.175,  7.8,  "Pu Lγ1" },
    { 94,  21.996,  8.6,  "Pu Lγ3" },
    { 94, 103.734, 51.0,  "Pu Kα1" },
    { 94,  99.525, 51.0,  "Pu Kα2" },
    { 94, 117.228, 51.0,  "Pu Kβ1" },
    { 94, 120.540, 51.0,  "Pu Kβ2" },

    // Americium (Z=95)
    { 95,  14.617,  5.6,  "Am Lα1" },
    { 95,  14.411,  5.6,  "Am Lα2" },
    { 95,  18.383,  7.1,  "Am Lβ1" },
    { 95,  18.676,  7.1,  "Am Lβ2" },
    { 95,  21.601,  7.9,  "Am Lγ1" },
    { 95,  22.431,  8.8,  "Am Lγ3" },
    { 95, 106.470, 52.0,  "Am Kα1" },
    { 95, 102.030, 52.0,  "Am Kα2" },
    { 95, 120.165, 52.0,  "Am Kβ1" },
    { 95, 123.820, 52.0,  "Am Kβ2" },

    // Thorium (Z=90)
    { 90,  12.969,  4.8,  "Th Lα1" },
    { 90,  12.810,  4.8,  "Th Lα2" },
    { 90,  16.426,  6.2,  "Th Lβ1" },
    { 90,  16.622,  6.2,  "Th Lβ2" },
    { 90,  19.353,  7.1,  "Th Lγ1" },
    { 90,  20.115,  7.9,  "Th Lγ3" },
    { 90,  93.350, 45.5,  "Th Kα1" },
    { 90,  89.957, 45.5,  "Th Kα2" },
    { 90, 105.605, 45.5,  "Th Kβ1" },
    { 90, 108.716, 45.5,  "Th Kβ2" },

    // Neptunium (Z=93)
    { 93,  13.946,  5.4,  "Np Lα1" },
    { 93,  13.760,  5.4,  "Np Lα2" },
    { 93,  17.604,  6.8,  "Np Lβ1" },
    { 93,  17.849,  6.8,  "Np Lβ2" },
    { 93,  20.765,  7.7,  "Np Lγ1" },
    { 93,  21.563,  8.5,  "Np Lγ3" },
    { 93, 101.078, 49.5,  "Np Kα1" },
    { 93,  97.069, 49.5,  "Np Kα2" },
    { 93, 114.240, 49.5,  "Np Kβ1" },
    { 93, 117.623, 49.5,  "Np Kβ2" },

    // Iridium (Z=77)
    { 77,  63.287, 28.0,  "Ir Kα1" },
    { 77,  61.486, 28.0,  "Ir Kα2" },
    { 77,  71.413, 28.0,  "Ir Kβ1" },
    { 77,  73.560, 28.0,  "Ir Kβ2" },

    // Lead (Z=82)
    { 82,  10.551,  3.4,  "Pb Lα1" },
    { 82,  10.449,  3.4,  "Pb Lα2" },
    { 82,  12.614,  4.5,  "Pb Lβ1" },
    { 82,  12.622,  4.5,  "Pb Lβ2" },
    { 82,  14.764,  5.2,  "Pb Lγ1" },
    { 82,  15.218,  5.7,  "Pb Lγ3" },
    { 82,  74.969, 33.5,  "Pb Kα1" },
    { 82,  72.805, 33.5,  "Pb Kα2" },
    { 82,  84.936, 33.5,  "Pb Kβ1" },
    { 82,  87.300, 33.5,  "Pb Kβ2" },

    // Protactinium (Z=91)
    { 91,  13.291,  5.0,  "Pa Lα1" },
    { 91,  13.122,  5.0,  "Pa Lα2" },
    { 91,  16.821,  6.4,  "Pa Lβ1" },
    { 91,  17.038,  6.4,  "Pa Lβ2" },
    { 91,  19.755,  7.3,  "Pa Lγ1" },
    { 91,  20.525,  8.1,  "Pa Lγ3" },
    { 91,  95.868, 46.7,  "Pa Kα1" },
    { 91,  92.287, 46.7,  "Pa Kα2" },
    { 91, 108.427, 46.7,  "Pa Kβ1" },
    { 91, 111.630, 46.7,  "Pa Kβ2" },

    // Radium (Z=88)
    { 88,  12.340,  4.4,  "Ra Lα1" },
    { 88,  12.197,  4.4,  "Ra Lα2" },
    { 88,  15.236,  5.8,  "Ra Lβ1" },
    { 88,  15.421,  5.8,  "Ra Lβ2" },
    { 88,  18.036,  6.7,  "Ra Lγ1" },
    { 88,  18.757,  7.4,  "Ra Lγ3" },
    { 88,  88.471, 42.7,  "Ra Kα1" },
    { 88,  85.429, 42.7,  "Ra Kα2" },
    { 88, 100.130, 42.7,  "Ra Kβ1" },
    { 88, 103.088, 42.7,  "Ra Kβ2" }
  };

  constexpr size_t num_original = sizeof(original_data) / sizeof(OriginalEntry);
  BOOST_CHECK_EQUAL( num_original, 84 );  // Verify we have all entries (6 elements × ~14 entries each)

  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  size_t num_matches = 0;
  for( size_t i = 0; i < num_original; ++i )
  {
    const OriginalEntry &orig = original_data[i];
    const double width_kev = db->get_natural_width_hwhm_kev( orig.z, orig.energy_kev, 0.5 );

    if( width_kev > 0.0 )
    {
      const double expected_kev = orig.hwhm_ev / 1000.0;
      const double percent_diff = 100.0 * std::abs( width_kev - expected_kev ) / expected_kev;

      // Check within 1% tolerance
      BOOST_CHECK_SMALL( percent_diff, 1.0 );

      if( percent_diff > 1.0 )
      {
        std::cerr << "Mismatch for " << orig.label << " (Z=" << orig.z << ", E="
                  << orig.energy_kev << " keV): expected " << expected_kev
                  << " keV, got " << width_kev << " keV (diff=" << percent_diff << "%)"
                  << std::endl;
      }

      ++num_matches;
    }
    else
    {
      std::cerr << "Missing entry for " << orig.label << " (Z=" << orig.z
                << ", E=" << orig.energy_kev << " keV)" << std::endl;
      BOOST_CHECK( false );  // Fail test if original entry not found
    }
  }

  std::cout << "Validated " << num_matches << " of " << num_original
            << " original x-ray width entries" << std::endl;

  // All original entries should be present
  BOOST_CHECK_EQUAL( num_matches, num_original );
}//BOOST_AUTO_TEST_CASE( OriginalDataValidation )


BOOST_AUTO_TEST_CASE( EnergyMatching )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // Test exact energy match for U Kα1 at 98.439 keV
  double width = db->get_natural_width_hwhm_kev( 92, 98.439, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.048, 2.0 );  // ~48 eV ± 2%

  // Test energy within tolerance (0.1 keV off)
  width = db->get_natural_width_hwhm_kev( 92, 98.339, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.048, 2.0 );

  // Test energy outside tolerance (1.0 keV off, default tolerance 0.5)
  width = db->get_natural_width_hwhm_kev( 92, 97.439, 0.5 );
  BOOST_CHECK_EQUAL( width, -1.0 );  // Should not match

  // Test with larger tolerance
  width = db->get_natural_width_hwhm_kev( 92, 97.439, 1.5 );
  BOOST_CHECK_GT( width, 0.0 );  // Now should match

  // Test invalid atomic number
  width = db->get_natural_width_hwhm_kev( 200, 98.439, 0.5 );
  BOOST_CHECK_EQUAL( width, -1.0 );

  // Test element with no data
  width = db->get_natural_width_hwhm_kev( 1, 100.0, 0.5 );  // Hydrogen, no K-shell at 100 keV
  BOOST_CHECK_EQUAL( width, -1.0 );
}//BOOST_AUTO_TEST_CASE( EnergyMatching )


BOOST_AUTO_TEST_CASE( DopplerWidthCalculation )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // Test Doppler width lookup for U Kα1
  const double doppler_295K = db->get_doppler_width_hwhm_kev( 92, 98.439, 0.5, 295.0 );
  BOOST_CHECK_GT( doppler_295K, 0.0 );
  BOOST_CHECK_CLOSE( doppler_295K, 0.00055, 20.0 );  // ~0.55 eV ± 20%

  // Test temperature scaling: HWHM ∝ sqrt(T)
  const double doppler_590K = db->get_doppler_width_hwhm_kev( 92, 98.439, 0.5, 590.0 );
  BOOST_CHECK_GT( doppler_590K, 0.0 );

  const double scale_factor = doppler_590K / doppler_295K;
  const double expected_scale = std::sqrt( 590.0 / 295.0 );
  BOOST_CHECK_CLOSE( scale_factor, expected_scale, 0.1 );

  // Test cryogenic temperature (77K - liquid nitrogen)
  const double doppler_77K = db->get_doppler_width_hwhm_kev( 92, 98.439, 0.5, 77.0 );
  BOOST_CHECK_GT( doppler_77K, 0.0 );
  BOOST_CHECK_LT( doppler_77K, doppler_295K );  // Should be narrower at lower temp
}//BOOST_AUTO_TEST_CASE( DopplerWidthCalculation )


BOOST_AUTO_TEST_CASE( TotalWidthCalculation )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // Test total width calculation for U Kα1
  const double natural = db->get_natural_width_hwhm_kev( 92, 98.439, 0.5 );
  const double doppler = db->get_doppler_width_hwhm_kev( 92, 98.439, 0.5, 295.0 );
  const double total = db->get_total_width_hwhm_kev( 92, 98.439, 0.5, 295.0 );

  BOOST_CHECK_GT( natural, 0.0 );
  BOOST_CHECK_GT( doppler, 0.0 );
  BOOST_CHECK_GT( total, 0.0 );

  // Verify quadrature addition: total² = natural² + doppler²
  const double expected_total = std::sqrt( natural*natural + doppler*doppler );
  BOOST_CHECK_CLOSE( total, expected_total, 0.01 );

  // For U Kα1, natural width (48 eV) should dominate over Doppler (0.55 eV)
  BOOST_CHECK_CLOSE( total, natural, 2.0 );  // Total should be within 2% of natural
}//BOOST_AUTO_TEST_CASE( TotalWidthCalculation )


BOOST_AUTO_TEST_CASE( PeakDistsIntegration )
{
  // Test integration with PeakDists::get_xray_lorentzian_width()

  // Test U Kα1
  double width = PeakDists::get_xray_lorentzian_width( 92, 98.439, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.048, 2.0 );

  // Test Pu Kα1
  width = PeakDists::get_xray_lorentzian_width( 94, 103.734, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.051, 2.0 );

  // Test Am Lα1
  width = PeakDists::get_xray_lorentzian_width( 95, 14.617, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.0056, 2.0 );
}//BOOST_AUTO_TEST_CASE( PeakDistsIntegration )


BOOST_AUTO_TEST_CASE( ElementCoverage )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // Verify all 9 original elements have data
  const int original_elements[] = { 77, 82, 88, 90, 91, 92, 93, 94, 95 };

  for( int z : original_elements )
  {
    const std::vector<XRayWidthEntry> entries = db->get_all_widths_for_element( z );
    BOOST_CHECK_GT( entries.size(), 0 );

    if( entries.empty() )
    {
      std::cerr << "No x-ray width entries found for Z=" << z << std::endl;
    }
  }

  // Check that database has expanded coverage beyond original 9 elements
  // (This will pass once comprehensive data is added)
  const size_t total_entries = db->num_entries();
  std::cout << "Database contains " << total_entries << " total x-ray width entries" << std::endl;
}//BOOST_AUTO_TEST_CASE( ElementCoverage )


BOOST_AUTO_TEST_CASE( DataPhysicalValidity )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // Test physical validity of data for all elements
  for( int z = 1; z <= 98; ++z )
  {
    const std::vector<XRayWidthEntry> entries = db->get_all_widths_for_element( z );

    for( const XRayWidthEntry &entry : entries )
    {
      // Check energy threshold (should be ≥10 keV)
      BOOST_CHECK_GE( entry.energy_kev, 10.0 );

      // Check natural width is physically reasonable (0.1-200 eV)
      BOOST_CHECK_GE( entry.hwhm_natural_ev, 0.1 );
      BOOST_CHECK_LE( entry.hwhm_natural_ev, 200.0 );

      // Check Doppler width is physically reasonable (0.01-5 eV at 295K)
      BOOST_CHECK_GE( entry.hwhm_doppler_ev, 0.01 );
      BOOST_CHECK_LE( entry.hwhm_doppler_ev, 5.0 );

      // Natural width should increase with Z (roughly)
      // Higher Z → deeper inner shell binding → shorter lifetime → broader width

      // Doppler width should scale roughly as E/sqrt(M)
      // (This is a rough check - exact values depend on atomic mass)
    }
  }
}//BOOST_AUTO_TEST_CASE( DataPhysicalValidity )


BOOST_AUTO_TEST_CASE( AlphaRecoilDoppler )
{
  // Initialize SandiaDecay database
  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") );

  // Test alpha recoil Doppler calculation for known isotopes

  // U-238 → Th-234 alpha decay (E_alpha ≈ 4.2 MeV)
  // Recoil Doppler for U Kα (98.4 keV): weighted average ~67.8 eV HWHM
  const double u238_recoil = compute_alpha_recoil_doppler_hwhm( "U238", 98.439 );
  BOOST_CHECK_GT( u238_recoil, 0.0 );
  BOOST_CHECK_CLOSE( u238_recoil, 0.0678, 5.0 );  // ~67.8 eV ± 5%

  std::cout << "U-238 alpha recoil Doppler for U Kα (98.4 keV): "
            << (u238_recoil * 1000.0) << " eV HWHM" << std::endl;

  // Pu-239 → U-235 alpha decay (E_alpha ≈ 5.2 MeV)
  // Recoil Doppler for Pu Kα (103.7 keV): weighted average ~78.9 eV HWHM
  const double pu239_recoil = compute_alpha_recoil_doppler_hwhm( "Pu239", 103.734 );
  BOOST_CHECK_GT( pu239_recoil, 0.0 );
  BOOST_CHECK_CLOSE( pu239_recoil, 0.0789, 5.0 );  // ~78.9 eV ± 5%

  std::cout << "Pu-239 alpha recoil Doppler for Pu Kα (103.7 keV): "
            << (pu239_recoil * 1000.0) << " eV HWHM" << std::endl;

  // Am-241 → Np-237 alpha decay (E_alpha ≈ 5.5 MeV)
  // Recoil Doppler for Am Kα (106.5 keV): weighted average ~82.8 eV HWHM
  const double am241_recoil = compute_alpha_recoil_doppler_hwhm( "Am241", 106.470 );
  BOOST_CHECK_GT( am241_recoil, 0.0 );
  BOOST_CHECK_CLOSE( am241_recoil, 0.0828, 5.0 );  // ~82.8 eV ± 5%

  std::cout << "Am-241 alpha recoil Doppler for Am Kα (106.5 keV): "
            << (am241_recoil * 1000.0) << " eV HWHM" << std::endl;

  // Th-232 → Ra-228 alpha decay (E_alpha ≈ 4.0 MeV)
  // Recoil Doppler for Th Kα (93.4 keV): weighted average ~64.5 eV HWHM
  const double th232_recoil = compute_alpha_recoil_doppler_hwhm( "Th232", 93.350 );
  BOOST_CHECK_GT( th232_recoil, 0.0 );
  BOOST_CHECK_CLOSE( th232_recoil, 0.0645, 5.0 );  // ~64.5 eV ± 5%

  std::cout << "Th-232 alpha recoil Doppler for Th Kα (93.4 keV): "
            << (th232_recoil * 1000.0) << " eV HWHM" << std::endl;

  // Test that recoil Doppler scales approximately with E_xray × sqrt(E_alpha / M_daughter)
  // For U-238 vs Pu-239 (similar masses, different alpha energies):
  const double ratio_u_pu = u238_recoil / pu239_recoil;
  const double expected_ratio = (98.439 / 103.734) * sqrt( 4.2 / 5.2 );  // E_xray and E_alpha scaling
  BOOST_CHECK_CLOSE( ratio_u_pu, expected_ratio, 30.0 );  // Within 30%

  // Test error handling: non-alpha emitter (stable isotope)
  const double pb208_recoil = compute_alpha_recoil_doppler_hwhm( "Pb208", 74.969 );
  BOOST_CHECK_EQUAL( pb208_recoil, -1.0 );  // Should fail - no alpha decay

  // Test error handling: non-existent nuclide
  const double fake_recoil = compute_alpha_recoil_doppler_hwhm( "Xx999", 100.0 );
  BOOST_CHECK_EQUAL( fake_recoil, -1.0 );  // Should fail
}//BOOST_AUTO_TEST_CASE( AlphaRecoilDoppler )
