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

#define BOOST_TEST_MODULE test_DecayActivityCalc_suite
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

#include <Wt/WString.h>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayActivityCalc.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace DecayActivityCalc;

// Data directory global (set from the command line).
std::string g_data_dir = "";

namespace
{
  /** Reads the data directory from the command line and points the decay database at it.
   Runs once, on first use, and must happen before the first DecayDataBaseServer::database() call.
   */
  const SandiaDecay::SandiaDecayDataBase *database()
  {
    static bool s_have_init = false;
    if( !s_have_init )
    {
      s_have_init = true;

      const int argc = boost::unit_test::framework::master_test_suite().argc;
      char **argv = boost::unit_test::framework::master_test_suite().argv;

      for( int i = 1; i < argc; ++i )
      {
        const string arg = argv[i];
        if( arg.find("--datadir=") == 0 )
          g_data_dir = arg.substr( 10 );
      }

      SpecUtils::ireplace_all( g_data_dir, "%20", " " );

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

      BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Could not find the data directory (--datadir=...)" );

      const string decay_xml = SpecUtils::append_path( g_data_dir, "sandia.decay.xml" );
      BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(decay_xml), "No sandia.decay.xml at " + decay_xml );

      BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( decay_xml ) );
      BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
    }//if( !s_have_init )

    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE( db );
    return db;
  }//database()


  const SandiaDecay::Nuclide *nuclide( const string &name )
  {
    const SandiaDecay::Nuclide * const nuc = database()->nuclide( name );
    BOOST_REQUIRE_MESSAGE( nuc, "No nuclide " + name );
    return nuc;
  }


  /** Splits "path?query" as the GUI does, and parses it. */
  ToolState parseUrl( const string &url, const bool defaultUseBq, vector<Wt::WString> &errors )
  {
    const string::size_type qpos = url.find( '?' );
    const string path = url.substr( 0, qpos );
    const string query = (qpos == string::npos) ? string() : url.substr( qpos + 1 );
    return fromAppUrl( path, query, defaultUseBq, errors );
  }


  bool close( const double a, const double b, const double rel = 1.0E-6 )
  {
    return fabs(a - b) <= rel * std::max( fabs(a), fabs(b) );
  }
}//namespace


BOOST_AUTO_TEST_CASE( url_round_trip )
{
  const SandiaDecay::Nuclide * const ba133 = nuclide( "Ba133" );
  const SandiaDecay::Nuclide * const cs137 = nuclide( "Cs137" );

  ToolState state;
  state.tab = Tab::Chain;
  state.gridLines = false;
  state.logY = true;
  state.chartType = ChartType::Gamma;
  state.displayUnit = "kBq";
  state.timeSpan = "5 hl";

  SourceNuclide ba;
  ba.nuclide = ba133;
  ba.activity = 3.0 * PhysicalUnits::microCi;
  ba.activityStr = "3 uCi";
  ba.useCurie = true;
  state.nuclides.push_back( ba );

  SourceNuclide cs;
  cs.nuclide = cs137;
  cs.activity = 37.0 * PhysicalUnits::kBq;
  cs.activityStr = "37 kBq";
  cs.useCurie = false;
  cs.age = 2.0 * PhysicalUnits::year;
  state.nuclides.push_back( cs );

  const string url = toAppUrl( state );
  BOOST_TEST_MESSAGE( "Round-trip URL: " << url );

  // The default units should not matter, as the URL says them.
  for( const bool defaultUseBq : { false, true } )
  {
    vector<Wt::WString> errors;
    const ToolState parsed = parseUrl( url, defaultUseBq, errors );

    BOOST_CHECK( errors.empty() );
    BOOST_CHECK( parsed.tab == Tab::Chain );
    BOOST_CHECK( parsed.gridLines == false );
    BOOST_CHECK( parsed.logY == true );
    BOOST_CHECK( parsed.chartType == ChartType::Gamma );
    BOOST_CHECK_EQUAL( parsed.displayUnit, "kBq" );
    BOOST_CHECK_EQUAL( parsed.timeSpan, "5 hl" );
    BOOST_REQUIRE_EQUAL( parsed.nuclides.size(), 2 );

    for( size_t i = 0; i < 2; ++i )
    {
      const SourceNuclide &expected = state.nuclides[i];
      const SourceNuclide &got = parsed.nuclides[i];
      BOOST_CHECK( got.nuclide == expected.nuclide );
      BOOST_CHECK( close( got.activity, expected.activity ) );
      BOOST_CHECK_EQUAL( got.activityStr, expected.activityStr );
      BOOST_CHECK_EQUAL( got.useCurie, expected.useCurie );
      BOOST_CHECK( (expected.age == 0.0) ? (got.age == 0.0) : close( got.age, expected.age ) );
    }
  }//for( const bool defaultUseBq : { false, true } )

  // Display units differing only by case must survive exactly.
  for( const string unit : { "mBq", "MBq", "microCi", DecayActivityCalc::sm_arbitraryDisplayUnit } )
  {
    state.displayUnit = unit;
    vector<Wt::WString> errors;
    const ToolState parsed = parseUrl( toAppUrl( state ), false, errors );
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK_EQUAL( parsed.displayUnit, unit );
  }

  // Each tab and chart type round-trips
  for( const Tab tab : { Tab::Chart, Tab::Chain, Tab::Calc } )
  {
    for( const ChartType type : { ChartType::Activity, ChartType::Gamma, ChartType::Beta, ChartType::Alpha } )
    {
      state.tab = tab;
      state.chartType = type;
      vector<Wt::WString> errors;
      const ToolState parsed = parseUrl( toAppUrl( state ), false, errors );
      BOOST_CHECK( errors.empty() );
      BOOST_CHECK( parsed.tab == tab );
      BOOST_CHECK( parsed.chartType == type );
    }
  }
}//BOOST_AUTO_TEST_CASE( url_round_trip )


BOOST_AUTO_TEST_CASE( documented_examples )
{
  {// "interspec://decay/chain?nuclide=U238&activity=3uCi&initialage=20y&timespan=22y&actunits=ci"
    vector<Wt::WString> errors;
    const ToolState state = fromAppUrl( "chain", "nuclide=U238&activity=3uCi&initialage=20y&timespan=22y&actunits=ci",
                                        true, errors );
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK( state.tab == Tab::Chain );
    BOOST_CHECK_EQUAL( state.timeSpan, "22y" );
    BOOST_CHECK_EQUAL( state.displayUnit, "mCi" ); // From "actunits", as there is no "dispunits"
    BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
    BOOST_CHECK( state.nuclides[0].nuclide == nuclide("U238") );
    BOOST_CHECK( close( state.nuclides[0].activity, 3.0*PhysicalUnits::microCi ) );
    BOOST_CHECK_EQUAL( state.nuclides[0].activityStr, "3uCi" );
    BOOST_CHECK( state.nuclides[0].useCurie );
    BOOST_CHECK( close( state.nuclides[0].age, 20.0*PhysicalUnits::year ) );
  }

  {// "interspec://decay/chart?nuc=Ba133&act=3uCi&nuc=Cs137&act=2uci&actunits=ci&timespan=20y"
    vector<Wt::WString> errors;
    const ToolState state = fromAppUrl( "chart", "nuc=Ba133&act=3uCi&nuc=Cs137&act=2uci&actunits=ci&timespan=20y",
                                        true, errors );
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK( state.tab == Tab::Chart );
    BOOST_CHECK_EQUAL( state.timeSpan, "20y" );
    BOOST_REQUIRE_EQUAL( state.nuclides.size(), 2 );
    BOOST_CHECK( state.nuclides[0].nuclide == nuclide("Ba133") );
    BOOST_CHECK( close( state.nuclides[0].activity, 3.0*PhysicalUnits::microCi ) );
    BOOST_CHECK( state.nuclides[1].nuclide == nuclide("Cs137") );
    BOOST_CHECK( close( state.nuclides[1].activity, 2.0*PhysicalUnits::microCi ) );
    BOOST_CHECK( state.nuclides[0].age == 0.0 );
  }

  {// "interspec://decay/chain?nuc=232-th"; activity defaults to 1 mCi, in the default units.
    for( const bool useBq : { false, true } )
    {
      vector<Wt::WString> errors;
      const ToolState state = fromAppUrl( "chain", "nuc=232-th", useBq, errors );
      BOOST_CHECK( errors.empty() );
      BOOST_CHECK( state.tab == Tab::Chain );
      BOOST_CHECK( state.displayUnit.empty() );
      BOOST_CHECK( state.timeSpan.empty() );
      BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
      BOOST_CHECK( state.nuclides[0].nuclide == nuclide("Th232") );
      BOOST_CHECK( close( state.nuclides[0].activity, PhysicalUnits::mCi ) );
      BOOST_CHECK_EQUAL( state.nuclides[0].useCurie, !useBq );
      BOOST_CHECK_EQUAL( state.nuclides[0].activityStr, (useBq ? "37 MBq" : "1 mCi") );
    }
  }
}//BOOST_AUTO_TEST_CASE( documented_examples )


BOOST_AUTO_TEST_CASE( lenient_parsing_and_errors )
{
  {// Keys and values are case-insensitive; an unknown key is one error, and the rest still applies
    vector<Wt::WString> errors;
    const ToolState state = fromAppUrl( "CHAIN", "NUC=Co60&ChartType=Gamma&logy=TRUE&Grid=Yes&foo=bar",
                                        false, errors );
    BOOST_CHECK_EQUAL( errors.size(), 1 );
    BOOST_CHECK( state.tab == Tab::Chain );
    BOOST_CHECK( state.chartType == ChartType::Gamma );
    BOOST_CHECK( state.logY );
    BOOST_CHECK( state.gridLines );
    BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
    BOOST_CHECK( state.nuclides[0].nuclide == nuclide("Co60") );
  }

  {// An empty path is the chart; a bad path is an error
    vector<Wt::WString> errors;
    BOOST_CHECK( fromAppUrl( "", "nuc=Co60", false, errors ).tab == Tab::Chart );
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK( fromAppUrl( "nonsense", "?nuc=Co60", false, errors ).tab == Tab::Chart );
    BOOST_CHECK_EQUAL( errors.size(), 1 );
  }

  {// An unparsable activity falls back to 1 mCi, with an activity string that matches it
    vector<Wt::WString> errors;
    const ToolState state = fromAppUrl( "chart", "nuc=Cs137&act=garbage", false, errors );
    BOOST_CHECK_EQUAL( errors.size(), 1 );
    BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
    BOOST_CHECK( close( state.nuclides[0].activity, PhysicalUnits::mCi ) );
    BOOST_CHECK( close( PhysicalUnits::stringToActivity(state.nuclides[0].activityStr), PhysicalUnits::mCi, 1.0E-3 ) );
  }

  {// Units written in the activity win over the default (and, with "dispunits", "actunits"); the
   //  order of a single nuclides arguments does not matter.
    for( const string query : { "nuc=Cs137&act=10kBq", "act=10kBq&nuc=Cs137",
                                "dispunits=MBq&actunits=ci&act=10kBq&nuc=Cs137" } )
    {
      vector<Wt::WString> errors;
      const ToolState state = fromAppUrl( "chart", query, false, errors );
      BOOST_CHECK( errors.empty() );
      BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
      BOOST_CHECK( state.nuclides[0].nuclide == nuclide("Cs137") );
      BOOST_CHECK( close( state.nuclides[0].activity, 10.0*PhysicalUnits::kBq ) );
      BOOST_CHECK_EQUAL( state.nuclides[0].activityStr, "10kBq" );
      BOOST_CHECK( !state.nuclides[0].useCurie );
    }
  }

  {// URLs from before "dispunits" wrote "act" in curies, with the users units only in "actunits"
    vector<Wt::WString> errors;
    const ToolState state = fromAppUrl( "chart", "grid=0&logy=0&actunits=bq&charttype=activity"
                                        "&timespan=301.000000 y&nuc=Cs137&act=27.0270 uCi", false, errors );
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK_EQUAL( state.displayUnit, "MBq" );
    BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
    BOOST_CHECK( !state.nuclides[0].useCurie );
    BOOST_CHECK( close( state.nuclides[0].activity, PhysicalUnits::MBq, 1.0E-4 ) );
    // Reformatted into the users units (the old URL rounded it to "999.999 kBq")
    BOOST_CHECK( activityStringIsCurie( state.nuclides[0].activityStr ) == false );
    BOOST_CHECK( close( PhysicalUnits::stringToActivity(state.nuclides[0].activityStr), state.nuclides[0].activity, 1.0E-3 ) );
  }

  {// Display units: exact match, else an unambiguous case-insensitive one
    vector<Wt::WString> errors;
    BOOST_CHECK_EQUAL( fromAppUrl( "chart", "dispunits=kbq", false, errors ).displayUnit, "kBq" );
    BOOST_CHECK_EQUAL( fromAppUrl( "chart", "dispunits=ARBITRARY", false, errors ).displayUnit,
                       DecayActivityCalc::sm_arbitraryDisplayUnit );
    BOOST_CHECK_EQUAL( fromAppUrl( "chart", "dispunits=uCi", false, errors ).displayUnit, "microCi" ); // by value
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK( fromAppUrl( "chart", "dispunits=mbq", false, errors ).displayUnit.empty() ); // mBq or MBq?
    BOOST_CHECK_EQUAL( errors.size(), 1 );
  }

  {// Bad and stable nuclides, ages, units, chart types, and time spans are each an error
    vector<Wt::WString> errors;
    const ToolState state = fromAppUrl( "chart",
                     "nuc=Xx999&nuc=Fe56&nuc=Am241&age=-3y&actunits=furlongs&charttype=neutron&timespan=forever",
                     false, errors );
    BOOST_CHECK_EQUAL( errors.size(), 6 );
    BOOST_REQUIRE_EQUAL( state.nuclides.size(), 1 );
    BOOST_CHECK( state.nuclides[0].nuclide == nuclide("Am241") );
    BOOST_CHECK( state.nuclides[0].age == 0.0 );
    BOOST_CHECK( state.chartType == ChartType::Activity );
    BOOST_CHECK( state.timeSpan.empty() );
  }

  {// Half-life time spans need a nuclide
    vector<Wt::WString> errors;
    BOOST_CHECK_EQUAL( fromAppUrl( "chart", "nuc=Co60&timespan=5hl", false, errors ).timeSpan, "5hl" );
    BOOST_CHECK( errors.empty() );
    BOOST_CHECK( fromAppUrl( "chart", "timespan=5hl", false, errors ).timeSpan.empty() );
    BOOST_CHECK_EQUAL( errors.size(), 1 );
  }
}//BOOST_AUTO_TEST_CASE( lenient_parsing_and_errors )


BOOST_AUTO_TEST_CASE( activity_strings )
{
  BOOST_CHECK( activityStringIsCurie( "3 uCi" ) == true );
  BOOST_CHECK( activityStringIsCurie( "2 microCi" ) == true );
  BOOST_CHECK( activityStringIsCurie( "5 curie" ) == true );
  BOOST_CHECK( activityStringIsCurie( "10kBq" ) == false );
  BOOST_CHECK( activityStringIsCurie( "1.5 microBq" ) == false ); // the 'c' in "micro" is not curies
  BOOST_CHECK( activityStringIsCurie( "1e3 becquerel" ) == false );
  BOOST_CHECK( !activityStringIsCurie( "7" ).has_value() );

  SourceNuclide src;
  src.nuclide = nuclide( "Co60" );
  src.activity = PhysicalUnits::mCi;
  src.useCurie = true;
  normalizeActivityString( src );
  BOOST_CHECK_EQUAL( src.activityStr, "1 mCi" );

  src.activityStr = "5 uCi"; // Does not match the activity
  normalizeActivityString( src );
  BOOST_CHECK_EQUAL( src.activityStr, "1 mCi" );

  src.activityStr = "1.0 millicurie"; // Matches, so is kept as the user wrote it
  normalizeActivityString( src );
  BOOST_CHECK_EQUAL( src.activityStr, "1.0 millicurie" );

  src.activityStr.clear();
  src.useCurie = false;
  src.activity = 37.5 * PhysicalUnits::kBq;
  normalizeActivityString( src );
  BOOST_CHECK_EQUAL( src.activityStr, "37.5 kBq" );

  BOOST_CHECK_EQUAL( printTimeForUrl( 40.42 * PhysicalUnits::hour ), "40.42 h" );
  BOOST_CHECK_EQUAL( printTimeForUrl( 2.0 * PhysicalUnits::year ), "2 y" );
  const double longTime = 2.5E10 * PhysicalUnits::year;  // printed with an exponent
  BOOST_CHECK( close( PhysicalUnits::stringToTimeDurationPossibleHalfLife( printTimeForUrl(longTime), -1.0 ), longTime ) );

  // Formatting keeps the activity to within the 0.1% the GUI checks.
  src.activityStr.clear();
  src.activity = 1.23456789 * PhysicalUnits::MBq;
  normalizeActivityString( src );
  BOOST_CHECK( close( PhysicalUnits::stringToActivity(src.activityStr), src.activity, 1.0E-3 ) );
}//BOOST_AUTO_TEST_CASE( activity_strings )


BOOST_AUTO_TEST_CASE( time_for_activity_fraction )
{
  {// A single nuclide with no progeny: log2(10) half-lives
    const SandiaDecay::Nuclide * const co60 = nuclide( "Co60" );
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( co60, PhysicalUnits::curie );
    const double t = findTimeForActivityFrac( mix, 0.1 );
    BOOST_CHECK_CLOSE( t, std::log2(10.0) * co60->halfLife, 1.0 ); // percent
  }

  {// With progeny building in, the total activity is what falls to the fraction
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide("Cs137"), PhysicalUnits::curie );
    const double t = findTimeForActivityFrac( mix, 0.1 );
    BOOST_CHECK_GT( t, 0.0 );
    BOOST_CHECK_SMALL( mix.totalActivity(t) / mix.totalActivity(0.0) - 0.1, 0.0011 );
  }

  {// No activity
    SandiaDecay::NuclideMixture mix;
    BOOST_CHECK_EQUAL( findTimeForActivityFrac( mix, 0.1 ), 0.0 );
  }

  {// Invalid inputs
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide("Co60"), PhysicalUnits::curie );
    BOOST_CHECK_THROW( findTimeForActivityFrac( mix, 0.0 ), std::exception );
    BOOST_CHECK_THROW( findTimeForActivityFrac( mix, 1.5 ), std::exception );
    BOOST_CHECK_THROW( findTimeForActivityFrac( mix, 0.1, 10.0, 1.0 ), std::exception );
  }
}//BOOST_AUTO_TEST_CASE( time_for_activity_fraction )


BOOST_AUTO_TEST_CASE( particle_rates_match_mixture )
{
  // The chart must agree with the summary and CSV, which use NuclideMixture; e.g., both skip the
  //  gamma records of pure E0 transitions.
  const SandiaDecay::SandiaDecayDataBase * const db = database();
  const double activity = 1000.0 * PhysicalUnits::becquerel;

  auto sum = []( const vector<SandiaDecay::EnergyRatePair> &lines ) -> double {
    double total = 0.0;
    for( const SandiaDecay::EnergyRatePair &line : lines )
      total += line.numPerSecond;
    return total;
  };

  size_t num_checked = 0, num_mismatch = 0, num_with_e0 = 0;
  for( const SandiaDecay::Nuclide * const nuc : db->nuclides() )
  {
    if( !nuc || nuc->isStable() || nuc->decaysToChildren.empty() )
      continue;

    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuc, activity );

    // The mixtures decay solution gives the activity at t=0 to ~1E-8, so compare at its value.
    const double act0 = mix.activity( 0.0, nuc );

    const double gammas = sum( mix.gammas( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy, false ) );
    const double betas = sum( mix.betas( 0.0 ) );
    const double alphas = sum( mix.alphas( 0.0 ) );

    // NuclideMixture keeps rates as floats, so allow float precision.
    auto matches = [&]( const double expected, const SandiaDecay::ProductType type ) -> bool {
      const double rate = particleRate( *nuc, act0, type );
      return fabs(expected - rate) <= (1.0E-5 * fabs(expected) + 1.0E-9 * activity);
    };

    const bool ok = matches( gammas, SandiaDecay::GammaParticle )
                    && matches( betas, SandiaDecay::BetaParticle )
                    && matches( alphas, SandiaDecay::AlphaParticle );
    if( !ok && (num_mismatch++ < 10) )
      BOOST_ERROR( "Particle rate for " << nuc->symbol << " does not match NuclideMixture" );

    ++num_checked;

    bool has_e0 = false;
    for( const SandiaDecay::Transition * const trans : nuc->decaysToChildren )
      for( const SandiaDecay::RadParticle &particle : trans->products )
        has_e0 |= ((particle.type == SandiaDecay::GammaParticle) && particle.e0_verified);
    num_with_e0 += has_e0;
  }//for( loop over nuclides )

  BOOST_TEST_MESSAGE( "Checked particle rates of " << num_checked << " nuclides, " << num_with_e0
                      << " of which have pure E0 gamma records" );
  BOOST_WARN_GT( num_with_e0, 0 ); // Else the E0 filtering is not exercised
  BOOST_CHECK_GT( num_checked, 1000 );
  BOOST_CHECK_EQUAL( num_mismatch, 0 );
}//BOOST_AUTO_TEST_CASE( particle_rates_match_mixture )
