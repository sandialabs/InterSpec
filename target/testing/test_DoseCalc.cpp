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
#include <memory>
#include <vector>
#include <iostream>


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DoseCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/DoseCalc.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DoseCalcWidget.h"
#include "InterSpec/GadrasShieldScatter.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;
using namespace boost::unit_test;


static string find_data_dir( const string file )
{
  // The "InterSpec/data" directory that contains all the cross-sections, reactions, etc
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  std::string datadir;
  
  for( int i = 1; i < argc; ++i )
  {
    const std::string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    const string trial_paths[] = {
      "data",
      "../data",
      "../../data",
      "../../../data",
      "external_libs/SandiaDecay/",
      "../external_libs/SandiaDecay/",
      "../../external_libs/SandiaDecay/",
      "../../../external_libs/SandiaDecay/",
      "/Users/wcjohns/rad_ana/InterSpec/data/external_libs/SandiaDecay/",
      "/Users/wcjohns/rad_ana/InterSpec/data"
    };
    
    for( const auto &d : trial_paths )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, file) ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  return datadir;
}//void find_data_dir()


// One-time setup of the static data dir, decay DB, and scatter table.
// All test cases share a single GadrasShieldScatter instance.
static GadrasShieldScatter *get_scatter()
{
  static unique_ptr<GadrasShieldScatter> s_scatter;
  static bool s_init = false;

  if( !s_init )
  {
    s_init = true;

    const string data_dir = find_data_dir( "sandia.shieldscatter.db" );
    const string scatter_file = SpecUtils::append_path( data_dir, "sandia.shieldscatter.db" );

    const string decay_xml_dir = find_data_dir( "sandia.decay.xml" );
    const string sandia_decay_file = SpecUtils::append_path( decay_xml_dir, "sandia.decay.xml" );

    BOOST_REQUIRE_MESSAGE( !sandia_decay_file.empty()
                          && SpecUtils::is_file( sandia_decay_file ),
                          "Error finding sandia.decay.xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( scatter_file ),
                          "Error finding sandia.shieldscatter.db at " + scatter_file );

    BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( sandia_decay_file ) );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( data_dir ) );

    BOOST_REQUIRE_NO_THROW( s_scatter.reset( new GadrasShieldScatter( scatter_file ) ) );
  }

  return s_scatter.get();
}


// Compute dose for `nuc` aged `age`, viewed at `distance` through a slab of
// areal density `areal_density` and effective Z `atomic_number`, and return
// the relative difference vs. `expected_dose` (all in PhysicalUnits internals).
// 100 µCi source activity matches the convention used in
// DoseCalcWidget::runtime_sanity_checks.
static void check_dose( const string &nuclabel,
                        const double age,
                        const float distance,
                        const float areal_density,
                        const float atomic_number,
                        const double expected_dose )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Nuclide database not initiated" );

  const SandiaDecay::Nuclide *nuc = db->nuclide( nuclabel );
  BOOST_REQUIRE_MESSAGE( nuc, "Could not retrieve '" + nuclabel + "' from nuclide database" );

  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, 100.0E-6 * SandiaDecay::curie, age );

  vector<float> energies, intensities;
  for( const auto &p : mix.photons( 0 ) )
  {
    energies.push_back( p.energy );
    intensities.push_back( p.numPerSecond );
  }

  const GadrasShieldScatter * const scatter = get_scatter();
  BOOST_REQUIRE( scatter );

  const double computed = DoseCalc::gamma_dose_with_shielding(
                              energies, intensities,
                              areal_density, atomic_number,
                              distance, *scatter );

  // 2% tolerance, computed relative to the larger of the two values to match
  // the in-app runtime_sanity_checks convention.
  const double denom = std::max( std::abs( computed ), std::abs( expected_dose ) );
  const double rel_diff = ( denom > 0.0 ) ? std::abs( computed - expected_dose ) / denom : 0.0;

  BOOST_CHECK_MESSAGE( rel_diff <= 0.02,
    nuclabel << " AN=" << atomic_number
      << " AD=" << ( areal_density * PhysicalUnits::cm2 / PhysicalUnits::gram ) << " g/cm2"
      << " age=" << PhysicalUnits::printToBestTimeUnits( age )
      << " dist=" << PhysicalUnits::printToBestLengthUnits( distance )
      << ": computed=" << PhysicalUnits::printToBestEquivalentDoseRateUnits( computed, 4, false )
      << " expected=" << PhysicalUnits::printToBestEquivalentDoseRateUnits( expected_dose, 4, false )
      << " rel_diff=" << rel_diff );
}


// End-to-end smoke check: run the same routine the widget triggers at load
// time, so a single failure here surfaces the broader "data files / build
// look broken" failure mode.
BOOST_AUTO_TEST_CASE( RuntimeSanityChecks )
{
  GadrasShieldScatter * const scatter = get_scatter();
  BOOST_REQUIRE( scatter );

  try
  {
    DoseCalcWidget::runtime_sanity_checks( scatter );
  }catch( std::exception &e )
  {
    BOOST_CHECK_MESSAGE( false, e.what() );
  }
}//BOOST_AUTO_TEST_CASE( RuntimeSanityChecks )


// ------------------------------------------------------------------------
// Per-scenario dose checks. Reference values are mirrored from
// DoseCalcWidget::runtime_sanity_checks (commit c3f40d0e). Splitting them
// out lets boost.test report exactly which scenario regresses.
// ------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Co60_NoShield_100cm )
{
  check_dose( "Co60",
              0.5 * PhysicalUnits::year,
              100.0f * PhysicalUnits::cm,
              0.0f,
              26.0f,
              115.4210E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Cs137_NoShield_100cm )
{
  check_dose( "Cs137",
              0.5 * PhysicalUnits::year,
              100.0f * PhysicalUnits::cm,
              0.0f,
              26.0f,
              29.6987E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Cs137_5gcm2_Fe_100cm )
{
  check_dose( "Cs137",
              0.5 * PhysicalUnits::year,
              100.0f * PhysicalUnits::cm,
              5.0f * static_cast<float>( PhysicalUnits::gram / PhysicalUnits::cm2 ),
              26.0f,
              25.4660E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Cs137_50gcm2_Fe_100cm )
{
  check_dose( "Cs137",
              0.5 * PhysicalUnits::year,
              100.0f * PhysicalUnits::cm,
              50.0f * static_cast<float>( PhysicalUnits::gram / PhysicalUnits::cm2 ),
              26.0f,
              2.9512E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( U238_NoShield_100cm )
{
  check_dose( "U238",
              20.0 * PhysicalUnits::year,
              100.0f * PhysicalUnits::cm,
              0.0f,
              60.0f,
              1.8638E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( U238_2gcm2_100cm )
{
  check_dose( "U238",
              20.0 * PhysicalUnits::year,
              100.0f * PhysicalUnits::cm,
              2.0f * static_cast<float>( PhysicalUnits::gram / PhysicalUnits::cm2 ),
              60.0f,
              0.8020508E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Na22_13gcm2_10cm )
{
  check_dose( "Na22",
              0.0,
              10.0f * PhysicalUnits::cm,
              13.0f * static_cast<float>( PhysicalUnits::gram / PhysicalUnits::cm2 ),
              5.0f,
              7.5218E-3 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Na22_NoShield_10cm )
{
  check_dose( "Na22",
              0.0,
              10.0f * PhysicalUnits::cm,
              0.0f,
              5.0f,
              10.81E-3 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( F18_NoShield_200cm )
{
  check_dose( "F18",
              0.0,
              200.0f * PhysicalUnits::cm,
              0.0f,
              5.0f,
              13.2870E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Ba133_NoShield_10cm )
{
  check_dose( "Ba133",
              7.0 * 24.0 * 3600.0 * PhysicalUnits::second,
              10.0f * PhysicalUnits::cm,
              0.0f,
              82.0f,
              2.4430E-3 * PhysicalUnits::rem / PhysicalUnits::hour );
}

BOOST_AUTO_TEST_CASE( Ba133_10gcm2_10cm )
{
  check_dose( "Ba133",
              7.0 * 24.0 * 3600.0 * PhysicalUnits::second,
              10.0f * PhysicalUnits::cm,
              10.0f * static_cast<float>( PhysicalUnits::gram / PhysicalUnits::cm2 ),
              82.0f,
              134.0646E-6 * PhysicalUnits::rem / PhysicalUnits::hour );
}


// ------------------------------------------------------------------------
// Dose-formatting checks. printToBestEquivalentDoseRateUnits uses snprintf
// with %.2f, whose trailing-zero behaviour can vary between platforms; the
// values below were chosen so the formatted output is platform-stable.
// ------------------------------------------------------------------------

namespace
{
  void check_print_dose( const double dose_rate_value, const string &expected, const bool use_sv )
  {
    const double scale = use_sv ? PhysicalUnits::sievert : PhysicalUnits::rem;
    const double dose_in_internal = dose_rate_value * scale / PhysicalUnits::hour;
    const string actual = PhysicalUnits::printToBestEquivalentDoseRateUnits( dose_in_internal, 2, use_sv );
    BOOST_CHECK_EQUAL( actual, expected );
  }
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Whole )
{
  check_print_dose( 1.234, "1.23 sv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Rem_Whole )
{
  check_print_dose( 1.234, "1.23 rem/hr", false );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_RoundUp )
{
  check_print_dose( 8.236, "8.24 sv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Nano )
{
  check_print_dose( 8.23688E-7, "823.69 nsv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Micro )
{
  check_print_dose( 8.23688E-5, "82.37 usv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Milli )
{
  check_print_dose( 8.23688E-2, "82.37 msv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Hundreds )
{
  check_print_dose( 8.23688E2, "823.69 sv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Kilo )
{
  check_print_dose( 8.23688E5, "823.69 ksv/hr", true );
}

BOOST_AUTO_TEST_CASE( PrintDose_Sv_Mega )
{
  check_print_dose( 8.23688E6, "8.24 Msv/hr", true );
}
