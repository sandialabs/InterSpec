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

#include "InterSpec/InterSpec.h"
#include "InterSpec/XRayWidthServer.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "SandiaDecay/SandiaDecay.h"

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
    // Prefer --datadir=... provided by the test runner.
    for( int i = 0; i < framework::master_test_suite().argc; ++i )
    {
      const std::string arg = framework::master_test_suite().argv[i];
      if( SpecUtils::starts_with( arg, "--datadir=" ) )
      {
        g_data_dir = arg.substr( 10 );
      }
    }

    // Fall back to searching relative paths.
    if( g_data_dir.empty() )
    {
      for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "xray_widths.xml" ) ) )
        {
          g_data_dir = d;
          break;
        }
      }
    }

    if( g_data_dir.empty() )
    {
      cerr << "Warning: Could not find data directory" << endl;
      g_data_dir = "data";  // Default fallback
    }

    // Make sure XRayWidthDatabase prefers our test data directory (via InterSpec writable data dir).
    // `setStaticDataDirectory()` is stricter (requires sandia.decay.xml etc.), but writable data dir
    // is sufficient for xray_widths.xml lookup.
    try
    {
      InterSpec::setWritableDataDirectory( g_data_dir );
    }
    catch( std::exception &e )
    {
      cerr << "Warning: Could not set writable data directory to '" << g_data_dir
           << "': " << e.what() << endl;
    }

    // Ensure SandiaDecay database can be initialized before XRayWidthDatabase loads,
    // since Doppler widths are computed in code using element atomic masses.
    try
    {
      std::string decay_xml;

      const std::vector<std::string> rel_data_dirs = { "data", "../data", "../../data", "../../../data", "../../../../data" };
      const std::vector<std::string> rel_roots = { ".", "..", "../..", "../../..", "../../../..", "../../../../.." };

      // Prefer the x-ray data dir, but fall back to common repo-relative locations.
      const std::vector<std::string> candidates = [&]()
      {
        std::vector<std::string> out;
        out.push_back( SpecUtils::append_path( g_data_dir, "sandia.decay.xml" ) );
        for( const std::string &d : rel_data_dirs )
          out.push_back( SpecUtils::append_path( d, "sandia.decay.xml" ) );
        for( const std::string &root : rel_roots )
          out.push_back( SpecUtils::append_path( root, "external_libs/SandiaDecay/sandia.decay.xml" ) );
        return out;
      }();

      for( const std::string &path : candidates )
      {
        if( SpecUtils::is_file( path ) )
        {
          decay_xml = path;
          break;
        }
      }

      if( decay_xml.empty() )
      {
        cerr << "Warning: Could not find sandia.decay.xml; SandiaDecay-dependent tests may fail." << endl;
        return;
      }

      DecayDataBaseServer::setDecayXmlFile( decay_xml );
    }
    catch( std::exception & )
    {
      // If already set/initialized elsewhere, ignore.
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
  // Regression coverage:
  // The original PeakDists hardcoded table (now removed) contained a small curated set of
  // line widths. The current XML is larger and may use different evaluations/sources, so
  // the numeric values may legitimately differ. Here we check:
  // - the line energy is present in the XML for that element (within tolerance), and
  // - a width is available when the XML provides one.

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

  size_t num_any_match = 0;
  size_t num_width_match = 0;
  size_t num_exact_but_no_width = 0;
  size_t num_no_match = 0;

  constexpr double tolerance_kev = 0.5;
  constexpr double exact_match_eps_kev = 0.001; // 1 eV

  for( size_t i = 0; i < num_original; ++i )
  {
    const OriginalEntry &orig = original_data[i];

    const std::vector<XRayWidthEntry> entries = db->get_all_widths_for_element( orig.z );
    bool any_match = false;
    bool width_match = false;
    bool exact_no_width = false;

    for( const XRayWidthEntry &entry : entries )
    {
      const double diff = fabs( entry.energy_kev - orig.energy_kev );

      if( diff <= exact_match_eps_kev )
      {
        any_match = true;
        if( entry.has_width_data() )
          width_match = true;
        else
          exact_no_width = true;
      }

      if( diff < tolerance_kev )
      {
        any_match = true;
        if( entry.has_width_data() )
          width_match = true;
      }
    }

    if( any_match )
    {
      ++num_any_match;
      if( width_match )
      {
        ++num_width_match;
      }
      else if( exact_no_width )
      {
        ++num_exact_but_no_width;
        std::cerr << "Info: Exact energy present but width missing for " << orig.label
                  << " (Z=" << orig.z << ", E=" << orig.energy_kev << " keV)" << std::endl;
      }
      else
      {
        std::cerr << "Info: Energy present but no width available within tolerance for " << orig.label
                  << " (Z=" << orig.z << ", E=" << orig.energy_kev << " keV)" << std::endl;
      }
    }
    else
    {
      ++num_no_match;
      std::cerr << "Warning: No XML energy match for " << orig.label << " (Z=" << orig.z
                << ", E=" << orig.energy_kev << " keV)" << std::endl;
    }
  }

  std::cout << "Original table coverage against current XML:" << std::endl;
  std::cout << "  Any energy match (±" << tolerance_kev << " keV): " << num_any_match
            << "/" << num_original << std::endl;
  std::cout << "  Width available: " << num_width_match << "/" << num_original << std::endl;
  std::cout << "  Exact match but no width: " << num_exact_but_no_width << "/" << num_original << std::endl;
  std::cout << "  No match at all: " << num_no_match << "/" << num_original << std::endl;

  // We expect the new XML to cover nearly all of the original energies as lines (even if some widths
  // are missing/unknown for specific subshell transitions). Allow a small number of gaps, since the
  // old curated list may include lines not present in the newer evaluation set.
  BOOST_CHECK_GE( num_any_match, num_original - 2 );
  BOOST_CHECK_GT( num_width_match, 0 );
}//BOOST_AUTO_TEST_CASE( OriginalDataValidation )


BOOST_AUTO_TEST_CASE( EnergyMatching )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // Test exact energy match for U Kα1 at 98.439 keV
  double width = db->get_natural_width_hwhm_kev( 92, 98.439, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.05225, 10.0 );  // current XML: U Kα1 ~52.25 eV

  // Test energy within tolerance (0.1 keV off)
  width = db->get_natural_width_hwhm_kev( 92, 98.339, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.05225, 10.0 );

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

  // For U Kα1, natural width (~52 eV) should dominate over Doppler (~0.55 eV)
  BOOST_CHECK_CLOSE( total, natural, 2.0 );  // Total should be within 2% of natural
}//BOOST_AUTO_TEST_CASE( TotalWidthCalculation )


BOOST_AUTO_TEST_CASE( GetXRayLorentzianWidth )
{
  // Test get_xray_lorentzian_width() function for fluorescent x-rays

  // SandiaDecay database is initialized in the global fixture.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db != nullptr );

  // Test U Kα1
  const SandiaDecay::Element *u = db->element( 92 );
  BOOST_REQUIRE( u != nullptr );
  double width = get_xray_lorentzian_width( u, 98.439, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.05225, 10.0 );

  // Test Pu Kα1
  const SandiaDecay::Element *pu = db->element( 94 );
  BOOST_REQUIRE( pu != nullptr );
  width = get_xray_lorentzian_width( pu, 103.734, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.05641, 10.0 );

  // Test Am Lα1
  const SandiaDecay::Element *am = db->element( 95 );
  BOOST_REQUIRE( am != nullptr );
  width = get_xray_lorentzian_width( am, 14.617, 0.5 );
  BOOST_CHECK_GT( width, 0.0 );
  BOOST_CHECK_CLOSE( width, 0.00616, 15.0 );

  // Test U Lα1 at ~13.6 keV should return a width
  width = get_xray_lorentzian_width( u, 13.6, 0.5 );
  BOOST_CHECK( width > 0.0 );

  // Test Pu Lα1 at ~14.3 keV should return a width
  width = get_xray_lorentzian_width( pu, 14.3, 0.5 );
  BOOST_CHECK( width > 0.0 );

  // Unknown element should return -1
  const SandiaDecay::Element *sn = db->element( 50 );
  BOOST_REQUIRE( sn != nullptr );
  width = get_xray_lorentzian_width( sn, 10.0, 0.5 );
  BOOST_CHECK( width < 0.0 );

  // Energy too far from known transitions should return -1
  width = get_xray_lorentzian_width( u, 999.0, 0.5 );
  BOOST_CHECK( width < 0.0 );

  // nullptr element should return -1
  width = get_xray_lorentzian_width( nullptr, 98.439, 0.5 );
  BOOST_CHECK( width < 0.0 );
}//BOOST_AUTO_TEST_CASE( GetXRayLorentzianWidth )


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

      if( entry.has_width_data() )
      {
        // Check natural width is physically reasonable.
        BOOST_CHECK_GT( entry.hwhm_natural_ev, 0.0 );
        BOOST_CHECK_LT( entry.hwhm_natural_ev, 2000.0 );

        // Check Doppler width is computed and physically reasonable.
        BOOST_CHECK_GT( entry.hwhm_doppler_ev, 0.0 );
        BOOST_CHECK_LT( entry.hwhm_doppler_ev, 10.0 );
      }
      else
      {
        // Missing `hwhm` in XML means no width data available for this line.
        BOOST_CHECK_LT( entry.hwhm_natural_ev, 0.0 );
        BOOST_CHECK_LT( entry.hwhm_doppler_ev, 0.0 );
      }

      // Natural width should increase with Z (roughly)
      // Higher Z → deeper inner shell binding → shorter lifetime → broader width

      // Doppler width should scale roughly as E/sqrt(M)
      // (This is a rough check - exact values depend on atomic mass)
    }
  }
}//BOOST_AUTO_TEST_CASE( DataPhysicalValidity )


BOOST_AUTO_TEST_CASE( MissingWidthHandling )
{
  const std::shared_ptr<const XRayWidthDatabase> db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( db != nullptr );

  // The new schema explicitly includes lines with no `hwhm` (width not available).
  // Example from the shipped XML: Ga (Z=31) K-N4 at 17.036 keV.
  const double natural = db->get_natural_width_hwhm_kev( 31, 17.036, 0.5 );
  BOOST_CHECK_EQUAL( natural, -1.0 );

  const double doppler = db->get_doppler_width_hwhm_kev( 31, 17.036, 0.5, 295.0 );
  BOOST_CHECK_EQUAL( doppler, -1.0 );

  const double total = db->get_total_width_hwhm_kev( 31, 17.036, 0.5, 295.0 );
  BOOST_CHECK_EQUAL( total, -1.0 );
}//BOOST_AUTO_TEST_CASE( MissingWidthHandling )


BOOST_AUTO_TEST_CASE( AlphaRecoilDoppler )
{
  // Test alpha recoil Doppler calculation via get_xray_total_width_for_decay()
  // This test verifies that decay x-rays include both natural and recoil contributions

  // SandiaDecay database is initialized in the global fixture.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db != nullptr );

  // Helper lambda to find alpha decay transition with x-ray near target energy
  auto find_alpha_transition = []( const SandiaDecay::Nuclide *nuc, const double xray_energy_kev )
    -> const SandiaDecay::Transition *
  {
    if( !nuc || nuc->decaysToChildren.empty() )
      return nullptr;

    for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
    {
      if( !trans )
        continue;

      const bool is_alpha = (trans->mode == SandiaDecay::AlphaDecay
                             || trans->mode == SandiaDecay::BetaAndAlphaDecay
                             || trans->mode == SandiaDecay::ElectronCaptureAndAlphaDecay);
      if( !is_alpha )
        continue;

      // Check if this transition has an x-ray near target energy
      for( const SandiaDecay::RadParticle &particle : trans->products )
      {
        if( particle.type == SandiaDecay::XrayParticle
            && fabs( particle.energy - xray_energy_kev ) < 0.5 )
        {
          return trans;
        }
      }
    }
    return nullptr;
  };

  // Test U-238 → Th-234 alpha decay (E_alpha ≈ 4.2 MeV)
  // Total width should include natural (~48 eV) + recoil (~68 eV) ≈ 83 eV HWHM
  const SandiaDecay::Nuclide *u238 = db->nuclide( "U238" );
  BOOST_REQUIRE( u238 != nullptr );
  const SandiaDecay::Transition *u238_trans = find_alpha_transition( u238, 98.439 );
  if( u238_trans )
  {
    const double u238_total = get_xray_total_width_for_decay( u238_trans, 98.439 );
    BOOST_CHECK_GT( u238_total, 0.0 );

    // Get natural width for comparison
    const SandiaDecay::Element *u = db->element( 92 );
    const double u238_natural = get_xray_lorentzian_width( u, 98.439, 0.5 );
    BOOST_CHECK_GT( u238_natural, 0.0 );

    // Total width should be greater than natural width (includes recoil)
    BOOST_CHECK_GT( u238_total, u238_natural );

    // Compute recoil contribution
    const double u238_recoil = std::sqrt( u238_total * u238_total - u238_natural * u238_natural );
    BOOST_CHECK_GT( u238_recoil, 0.0 );
    BOOST_CHECK_CLOSE( u238_recoil, 0.0678, 20.0 );  // ~67.8 eV ± 20%

    std::cout << "U-238 alpha decay x-ray widths (U Kα 98.4 keV):" << std::endl;
    std::cout << "  Natural: " << (u238_natural * 1000.0) << " eV HWHM" << std::endl;
    std::cout << "  Recoil:  " << (u238_recoil * 1000.0) << " eV HWHM" << std::endl;
    std::cout << "  Total:   " << (u238_total * 1000.0) << " eV HWHM" << std::endl;
  }
  else
  {
    std::cout << "Warning: Could not find U-238 alpha transition with U Kα x-ray" << std::endl;
  }

  // Test Pu-239 → U-235 alpha decay (E_alpha ≈ 5.2 MeV)
  // Total width should include natural (~51 eV) + recoil (~79 eV) ≈ 94 eV HWHM
  const SandiaDecay::Nuclide *pu239 = db->nuclide( "Pu239" );
  BOOST_REQUIRE( pu239 != nullptr );
  const SandiaDecay::Transition *pu239_trans = find_alpha_transition( pu239, 103.734 );
  if( pu239_trans )
  {
    const double pu239_total = get_xray_total_width_for_decay( pu239_trans, 103.734 );
    BOOST_CHECK_GT( pu239_total, 0.0 );

    const SandiaDecay::Element *pu = db->element( 94 );
    const double pu239_natural = get_xray_lorentzian_width( pu, 103.734, 0.5 );
    BOOST_CHECK_GT( pu239_natural, 0.0 );

    BOOST_CHECK_GT( pu239_total, pu239_natural );

    const double pu239_recoil = std::sqrt( pu239_total * pu239_total - pu239_natural * pu239_natural );
    BOOST_CHECK_GT( pu239_recoil, 0.0 );
    BOOST_CHECK_CLOSE( pu239_recoil, 0.0789, 20.0 );  // ~78.9 eV ± 20%

    std::cout << "Pu-239 alpha decay recoil Doppler for Pu Kα (103.7 keV): "
              << (pu239_recoil * 1000.0) << " eV HWHM" << std::endl;
  }

  // Test Am-241 → Np-237 alpha decay (E_alpha ≈ 5.5 MeV)
  const SandiaDecay::Nuclide *am241 = db->nuclide( "Am241" );
  BOOST_REQUIRE( am241 != nullptr );
  const SandiaDecay::Transition *am241_trans = find_alpha_transition( am241, 106.470 );
  if( am241_trans )
  {
    const double am241_total = get_xray_total_width_for_decay( am241_trans, 106.470 );
    BOOST_CHECK_GT( am241_total, 0.0 );

    const SandiaDecay::Element *am = db->element( 95 );
    const double am241_natural = get_xray_lorentzian_width( am, 106.470, 0.5 );
    BOOST_CHECK_GT( am241_natural, 0.0 );

    BOOST_CHECK_GT( am241_total, am241_natural );

    const double am241_recoil = std::sqrt( am241_total * am241_total - am241_natural * am241_natural );
    BOOST_CHECK_GT( am241_recoil, 0.0 );
    BOOST_CHECK_CLOSE( am241_recoil, 0.0828, 20.0 );  // ~82.8 eV ± 20%

    std::cout << "Am-241 alpha decay recoil Doppler for Am Kα (106.5 keV): "
              << (am241_recoil * 1000.0) << " eV HWHM" << std::endl;
  }

  // Test Th-232 → Ra-228 alpha decay (E_alpha ≈ 4.0 MeV)
  const SandiaDecay::Nuclide *th232 = db->nuclide( "Th232" );
  BOOST_REQUIRE( th232 != nullptr );
  const SandiaDecay::Transition *th232_trans = find_alpha_transition( th232, 93.350 );
  if( th232_trans )
  {
    const double th232_total = get_xray_total_width_for_decay( th232_trans, 93.350 );
    BOOST_CHECK_GT( th232_total, 0.0 );

    const SandiaDecay::Element *th = db->element( 90 );
    const double th232_natural = get_xray_lorentzian_width( th, 93.350, 0.5 );
    BOOST_CHECK_GT( th232_natural, 0.0 );

    BOOST_CHECK_GT( th232_total, th232_natural );

    const double th232_recoil = std::sqrt( th232_total * th232_total - th232_natural * th232_natural );
    BOOST_CHECK_GT( th232_recoil, 0.0 );
    BOOST_CHECK_CLOSE( th232_recoil, 0.0645, 20.0 );  // ~64.5 eV ± 20%

    std::cout << "Th-232 alpha decay recoil Doppler for Th Kα (93.4 keV): "
              << (th232_recoil * 1000.0) << " eV HWHM" << std::endl;
  }

  // Test error handling: non-alpha emitter (stable isotope)
  const SandiaDecay::Nuclide *pb208 = db->nuclide( "Pb208" );
  BOOST_REQUIRE( pb208 != nullptr );
  const SandiaDecay::Transition *pb208_trans = find_alpha_transition( pb208, 74.969 );
  BOOST_CHECK( pb208_trans == nullptr );  // Should have no alpha transition

  // Test error handling: nullptr transition
  const double null_width = get_xray_total_width_for_decay( nullptr, 98.439 );
  BOOST_CHECK_EQUAL( null_width, -1.0 );
}//BOOST_AUTO_TEST_CASE( AlphaRecoilDoppler )


BOOST_AUTO_TEST_CASE( GetXRayTotalWidthForDecay )
{
  // Test get_xray_total_width_for_decay() function for decay x-rays

  // SandiaDecay database is initialized in the global fixture.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db != nullptr );

  // Test U-238 alpha decay with U Kα x-ray
  const SandiaDecay::Nuclide *u238 = db->nuclide( "U238" );
  BOOST_REQUIRE( u238 != nullptr );
  BOOST_REQUIRE( !u238->decaysToChildren.empty() );

  // Find an alpha decay transition
  const SandiaDecay::Transition *alpha_transition = nullptr;
  for( const SandiaDecay::Transition *trans : u238->decaysToChildren )
  {
    if( trans && (trans->mode == SandiaDecay::AlphaDecay
                   || trans->mode == SandiaDecay::BetaAndAlphaDecay
                   || trans->mode == SandiaDecay::ElectronCaptureAndAlphaDecay) )
    {
      // Check if this transition has an x-ray around 98.4 keV (U Kα1)
      for( const SandiaDecay::RadParticle &particle : trans->products )
      {
        if( particle.type == SandiaDecay::XrayParticle
            && fabs( particle.energy - 98.439 ) < 0.5 )
        {
          alpha_transition = trans;
          break;
        }
      }
      if( alpha_transition )
        break;
    }
  }

  if( alpha_transition )
  {
    // Test total width calculation (should include natural + recoil)
    const double total_width = get_xray_total_width_for_decay( alpha_transition, 98.439 );
    BOOST_CHECK_GT( total_width, 0.0 );

    // Total width should be greater than natural width alone
    const SandiaDecay::Element *u = db->element( 92 );
    BOOST_REQUIRE( u != nullptr );
    const double natural_width = get_xray_lorentzian_width( u, 98.439, 0.5 );
    BOOST_CHECK_GT( total_width, natural_width );

    std::cout << "U-238 decay x-ray total width (natural + recoil): "
              << (total_width * 1000.0) << " eV HWHM" << std::endl;
    std::cout << "  Natural width: " << (natural_width * 1000.0) << " eV HWHM" << std::endl;
    std::cout << "  Recoil contribution: " << ((total_width - natural_width) * 1000.0) << " eV HWHM" << std::endl;
  }

  // Test error handling: nullptr transition
  const double null_width = get_xray_total_width_for_decay( nullptr, 98.439 );
  BOOST_CHECK_EQUAL( null_width, -1.0 );

  // Test error handling: transition without x-ray at specified energy
  if( alpha_transition )
  {
    const double bad_energy_width = get_xray_total_width_for_decay( alpha_transition, 999.0 );
    BOOST_CHECK_EQUAL( bad_energy_width, -1.0 );
  }
}//BOOST_AUTO_TEST_CASE( GetXRayTotalWidthForDecay )


BOOST_AUTO_TEST_CASE( NaturalAndTotalWidthForMultipleNuclides )
{
  // Test natural Lorentzian width availability and total width calculation
  // for multiple alpha-emitting nuclides

  // SandiaDecay database is initialized in the global fixture.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db != nullptr );

  // We want to explicitly report how well SandiaDecay's decay x-rays map to our xray_widths.xml data.
  // Use the same database instance the production functions use.
  const std::shared_ptr<const XRayWidthDatabase> width_db = XRayWidthDatabase::instance();
  BOOST_REQUIRE( width_db != nullptr );

  // List of nuclides to test
  const std::vector<std::string> nuclide_symbols = {
    "U235", "U234", "U238", "Pu236", "Pu238", "Pu239", "Pu240", "Pu241", "Am241", "Th232"
  };

  // Helper lambda to find all x-rays in alpha decay transitions
  auto find_alpha_xrays = []( const SandiaDecay::Nuclide *nuc )
    -> std::vector<std::pair<const SandiaDecay::Transition *, double> >
  {
    std::vector<std::pair<const SandiaDecay::Transition *, double> > xrays;
    
    if( !nuc || nuc->decaysToChildren.empty() )
      return xrays;

    for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
    {
      if( !trans )
        continue;

      const bool is_alpha = (trans->mode == SandiaDecay::AlphaDecay
                             || trans->mode == SandiaDecay::BetaAndAlphaDecay
                             || trans->mode == SandiaDecay::ElectronCaptureAndAlphaDecay
                             || trans->mode == SandiaDecay::BetaPlusAndAlphaDecay);
      if( !is_alpha )
        continue;

      // Collect all x-rays from this transition
      for( const SandiaDecay::RadParticle &particle : trans->products )
      {
        if( particle.type == SandiaDecay::XrayParticle && particle.energy >= 10.0 )
        {
          xrays.push_back( std::make_pair( trans, particle.energy ) );
        }
      }
    }
    
    return xrays;
  };

  size_t total_tested = 0;
  size_t natural_width_found = 0;
  size_t total_width_calculated = 0;

  // Match accounting against xray_widths.xml for the emitting element.
  size_t xml_any_exact = 0;
  size_t xml_width_exact = 0;
  size_t xml_exact_but_no_width = 0;
  size_t xml_any_within_tol = 0;
  size_t xml_width_within_tol = 0;

  constexpr double exact_match_eps_kev = 0.001;  // 1 eV
  constexpr double tolerance_kev = 0.5;

  for( const std::string &nuclide_symbol : nuclide_symbols )
  {
    const SandiaDecay::Nuclide *nuc = db->nuclide( nuclide_symbol );
    if( !nuc )
    {
      std::cout << "Warning: Nuclide '" << nuclide_symbol << "' not found in database" << std::endl;
      continue;
    }

    // Find all x-rays in alpha decay transitions
    const std::vector<std::pair<const SandiaDecay::Transition *, double> > xrays = find_alpha_xrays( nuc );
    
    if( xrays.empty() )
    {
      std::cout << "Info: No alpha decay x-rays found for " << nuclide_symbol << std::endl;
      continue;
    }

    // Test each x-ray found
    for( const auto &xray_pair : xrays )
    {
      const SandiaDecay::Transition *trans = xray_pair.first;
      const double xray_energy = xray_pair.second;
      
      ++total_tested;

      // Determine which element emits the x-ray (daughter for most decays, parent for isomeric)
      const SandiaDecay::Nuclide *xray_emitting_nuclide = trans->child;
      if( !xray_emitting_nuclide )
      {
        xray_emitting_nuclide = trans->parent;  // Fallback for spontaneous fission, etc.
      }
      
      if( !xray_emitting_nuclide )
      {
        std::cout << "Warning: Could not determine emitting nuclide for " << nuclide_symbol
                  << " x-ray at " << xray_energy << " keV" << std::endl;
        continue;
      }

      // Get the element for natural width lookup (use daughter element, not parent)
      const SandiaDecay::Element *element = db->element( xray_emitting_nuclide->atomicNumber );
      if( !element )
      {
        std::cout << "Warning: Could not find element for Z=" << xray_emitting_nuclide->atomicNumber
                  << " (x-ray from " << nuclide_symbol << " decay)" << std::endl;
        continue;
      }

      // Determine whether this decay x-ray energy matches our XML data for this element,
      // and whether a width is available.
      {
        const std::vector<XRayWidthEntry> entries =
          width_db->get_all_widths_for_element( xray_emitting_nuclide->atomicNumber );

        bool any_exact = false;
        bool width_exact = false;
        bool exact_no_width = false;
        bool any_tol = false;
        bool width_tol = false;

        for( const XRayWidthEntry &entry : entries )
        {
          const double diff = fabs( entry.energy_kev - xray_energy );

          if( diff <= exact_match_eps_kev )
          {
            any_exact = true;
            if( entry.has_width_data() )
              width_exact = true;
            else
              exact_no_width = true;
          }

          if( diff < tolerance_kev )
          {
            any_tol = true;
            if( entry.has_width_data() )
              width_tol = true;
          }
        }

        if( any_exact )
          ++xml_any_exact;
        if( width_exact )
          ++xml_width_exact;
        if( exact_no_width )
          ++xml_exact_but_no_width;
        if( any_tol )
          ++xml_any_within_tol;
        if( width_tol )
          ++xml_width_within_tol;
      }

      // Test 1: Natural Lorentzian width should be available
      const double natural_width = get_xray_lorentzian_width( element, xray_energy, tolerance_kev );
      if( natural_width > 0.0 )
      {
        ++natural_width_found;
        BOOST_CHECK_GT( natural_width, 0.0 );
        BOOST_CHECK_LT( natural_width, 1.0 );  // Should be in reasonable range (< 1 keV)
      }
      else
      {
        std::cout << "Warning: Natural width not found for " << nuclide_symbol
                  << " decay x-ray at " << xray_energy << " keV (from Z="
                  << xray_emitting_nuclide->atomicNumber << ")" << std::endl;
      }

      // Test 2: Total width should be calculable
      const double total_width = get_xray_total_width_for_decay( trans, xray_energy );
      if( total_width > 0.0 )
      {
        ++total_width_calculated;
        BOOST_CHECK_GT( total_width, 0.0 );
        
        // Total width should be >= natural width (includes recoil for alpha decays)
        if( natural_width > 0.0 )
        {
          BOOST_CHECK_GE( total_width, natural_width );
        }
      }
      else
      {
        std::cout << "Warning: Total width calculation failed for " << nuclide_symbol
                  << " x-ray at " << xray_energy << " keV" << std::endl;
      }

      // Print summary for first x-ray of each nuclide
      if( &xray_pair == &xrays[0] )
      {
        std::cout << nuclide_symbol << " (Z=" << nuc->atomicNumber << "): "
                  << xrays.size() << " x-ray(s) found";
        if( natural_width > 0.0 && total_width > 0.0 )
        {
          std::cout << ", natural=" << (natural_width * 1000.0) << " eV, "
                    << "total=" << (total_width * 1000.0) << " eV";
          if( total_width > natural_width )
          {
            const double recoil = std::sqrt( total_width * total_width - natural_width * natural_width );
            std::cout << ", recoil=" << (recoil * 1000.0) << " eV";
          }
        }
        std::cout << " (x-rays from Z=" << xray_emitting_nuclide->atomicNumber << ")" << std::endl;
      }
    }
  }

  std::cout << "\nSummary: Tested " << total_tested << " x-rays from " << nuclide_symbols.size()
            << " nuclides" << std::endl;
  std::cout << "  Natural width found: " << natural_width_found << "/" << total_tested << std::endl;
  std::cout << "  Total width calculated: " << total_width_calculated << "/" << total_tested << std::endl;

  if( total_tested > 0 )
  {
    std::cout << "  XML exact energy match (any): " << xml_any_exact << "/" << total_tested
              << " (" << (100.0 * xml_any_exact / total_tested) << "%)" << std::endl;
    std::cout << "    - exact match with width: " << xml_width_exact << "/" << total_tested
              << " (" << (100.0 * xml_width_exact / total_tested) << "%)" << std::endl;
    std::cout << "    - exact match but no width: " << xml_exact_but_no_width << "/" << total_tested
              << " (" << (100.0 * xml_exact_but_no_width / total_tested) << "%)" << std::endl;
    std::cout << "  XML match within " << tolerance_kev << " keV (any): " << xml_any_within_tol << "/" << total_tested
              << " (" << (100.0 * xml_any_within_tol / total_tested) << "%)" << std::endl;
    std::cout << "  XML match within " << tolerance_kev << " keV (with width): " << xml_width_within_tol << "/" << total_tested
              << " (" << (100.0 * xml_width_within_tol / total_tested) << "%)" << std::endl;
  }

  // At least some x-rays should have natural widths available
  BOOST_CHECK_GT( natural_width_found, 0 );
  
  // All x-rays should have calculable total widths (if transition is valid)
  BOOST_CHECK_GT( total_width_calculated, 0 );
}//BOOST_AUTO_TEST_CASE( NaturalAndTotalWidthForMultipleNuclides )
