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

/** Tests of the Activity/Shielding fit "Back decay" option: quoting the fit activities at a time
 other than the measurement (ShieldingSourceFitCalc::compute_back_decay), including fitted nuclides
 that feed other fitted nuclides (fission products), nuclides too many half-lives back, its XML
 serialization, and its use through a full fit and through the batch activity fit.
 */

#include "InterSpec_config.h"

#include <map>
#include <set>
#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <iostream>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ShieldSourceBackDecay_suite
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayBatchCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/SourceFitNuclideDisplay.h"


using namespace std;
using ShieldingSourceFitCalc::BackDecaySpec;
using ShieldingSourceFitCalc::SourceFitDef;
using ShieldingSourceFitCalc::BackDecayStatus;
using ShieldingSourceFitCalc::BackDecayResults;
using ShieldingSourceFitCalc::BackDecaySourceResult;


namespace
{
  std::string g_test_file_dir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char ** const argv = boost::unit_test::framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        g_test_file_dir = arg.substr( 14 );
    }//for( loop over arguments )

    SpecUtils::ireplace_all( datadir, "%20", " " );
    SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }//for( loop over candidate dirs )
    }//if( datadir.empty() )

    if( g_test_file_dir.empty() )
    {
      for( const char *d : { "test_data", "../test_data", "../../test_data" } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path(d, "ActShieldBackDecay") ) )
        {
          g_test_file_dir = d;
          break;
        }
      }//for( loop over candidate dirs )
    }//if( g_test_file_dir.empty() )

    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(datadir, "sandia.decay.xml") ),
                           "sandia.decay.xml not found in '" << datadir << "'" );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
    BOOST_REQUIRE( DecayDataBaseServer::database() );
  }//void set_data_dir()


  const SandiaDecay::Nuclide *nuc( const char *symbol )
  {
    set_data_dir();
    const SandiaDecay::Nuclide * const n = DecayDataBaseServer::database()->nuclide( symbol );
    BOOST_REQUIRE_MESSAGE( n, "No nuclide " << symbol );
    return n;
  }


  SourceFitDef make_src( const char *symbol, const double activity, const double age,
                         const double act_uncert = -1.0 )
  {
    SourceFitDef src;
    src.nuclide = nuc( symbol );
    src.activity = activity;
    src.fitActivity = true;
    src.age = age;
    src.fitAge = false;
    src.ageDefiningNuc = nullptr;
    src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
    if( act_uncert > 0.0 )
      src.activityUncertainty = act_uncert;
    return src;
  }//make_src(...)


  BackDecaySpec duration_spec( const double dt )
  {
    BackDecaySpec spec;
    spec.enabled = true;
    spec.mode = BackDecaySpec::Mode::DurationBeforeMeas;
    spec.duration = dt;
    return spec;
  }


  SpecUtils::time_point_t parse_time( const string &str )
  {
    const SpecUtils::time_point_t t = SpecUtils::time_from_string( str );
    BOOST_REQUIRE_MESSAGE( !SpecUtils::is_special(t), "Could not parse '" << str << "'" );
    return t;
  }


  const BackDecaySourceResult &result_for( const BackDecayResults &res, const char *symbol )
  {
    const BackDecaySourceResult * const src = res.source( nuc(symbol) );
    BOOST_REQUIRE_MESSAGE( src, "No back-decay result for " << symbol );
    return *src;
  }


  /** Activity of `progeny` at time `t` in a mixture that is pure `parent` of activity
   `parent_act` at t = 0. */
  double ingrowth( const char *parent, const double parent_act, const char *progeny, const double t )
  {
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuc(parent), parent_act );
    return static_cast<double>( mix.activity( t, nuc(progeny) ) );
  }


  /** Activity of `progeny` at the measurement, in a source that was pure `parent`, with activity
   `parent_act_now` at the measurement, an `age` ago. */
  double ingrowth_at_age( const char *parent, const double parent_act_now, const char *progeny,
                          const double age )
  {
    if( age <= 0.0 )
      return 0.0;
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc(parent), parent_act_now, age );
    return static_cast<double>( mix.activity( 0.0, nuc(progeny) ) );
  }


  /** A fission "release": at the reference time there is `parent_ref` of a parent and, from its
   independent fission yield, `progeny_ref` of its progeny; the measurement is `dt` later.

   Gives the parent and progeny source activities the activity fit would have (with interference
   correction, i.e., the progeny source being what the parent's source does not grow in), for a
   parent source age `fit_parent_age`, then checks back-decaying gives back the release activities.
   */
  void check_release( const char *parent, const double parent_ref, const char *progeny,
                      const double progeny_ref, const double dt, const double fit_parent_age )
  {
    const SandiaDecay::Nuclide * const p = nuc( parent ), * const d = nuc( progeny );

    SandiaDecay::NuclideMixture truth;
    truth.addNuclideByActivity( p, parent_ref );
    truth.addNuclideByActivity( d, progeny_ref );

    const double parent_now = static_cast<double>( truth.activity( dt, p ) );
    const double progeny_total_now = static_cast<double>( truth.activity( dt, d ) );

    // What the fit attributes to the progeny's own source: all of the progeny, less what the
    //  parent source, of the fit's age, grows in.
    const double progeny_src_now = progeny_total_now
                                    - ingrowth_at_age( parent, parent_now, progeny, fit_parent_age );

    const vector<SourceFitDef> srcs{ make_src( parent, parent_now, fit_parent_age ),
                                     make_src( progeny, progeny_src_now, 0.0 ) };
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    BOOST_REQUIRE_MESSAGE( res.error.empty(), res.error );

    const BackDecaySourceResult &pres = result_for( res, parent );
    const BackDecaySourceResult &dres = result_for( res, progeny );
    BOOST_REQUIRE( pres.status == BackDecayStatus::Ok );
    BOOST_REQUIRE( dres.status == BackDecayStatus::Ok );
    BOOST_CHECK( dres.has_fitted_ancestor );
    BOOST_CHECK( !pres.has_fitted_ancestor );

    BOOST_CHECK_MESSAGE( fabs(pres.activity - parent_ref) < 1.0E-8*parent_ref,
                         parent << " back-decayed to " << pres.activity << ", expected " << parent_ref );
    // If the fit says the parent is older than the time since release, the fit's own split of the
    //  progeny is kept (the parent then already had in-grown progeny at the reference time), and
    //  the progeny's own share just decays back; otherwise the release activity is recovered.
    //  (A fit activity can not be negative, so neither is the progeny's own share then.)
    const double expected_progeny = (fit_parent_age > dt)
                                     ? std::max(0.0, progeny_src_now) * std::pow( 2.0, dt/d->halfLife )
                                     : progeny_ref;
    BOOST_CHECK_MESSAGE( fabs(dres.activity - expected_progeny) < 1.0E-6*std::max(1.0, fabs(expected_progeny)),
                         progeny << " (parent age " << fit_parent_age/PhysicalUnits::hour << " h) back-decayed to "
                         << dres.activity << ", expected " << expected_progeny );

    // When the parent's age is shorter than the back-decay time, it is taken as pure at the
    //  reference time, and the user is told.
    const bool expect_note = (fit_parent_age < dt);
    BOOST_CHECK_MESSAGE( (dres.note.find("taken as pure") != string::npos) == expect_note,
                         "Unexpected note for " << progeny << ": '" << dres.note << "'" );
  }//check_release(...)
}//namespace


BOOST_AUTO_TEST_CASE( SingleNuclideScaling )
{
  set_data_dir();

  const double act = 1.0*PhysicalUnits::MBq, unc = 0.05*PhysicalUnits::MBq;
  const SourceFitDef cs137 = make_src( "Cs137", act, 5.0*PhysicalUnits::year, unc );
  const double hl = cs137.nuclide->halfLife;

  // Ten years before the measurement.
  const double dt = 10.0*PhysicalUnits::year;
  const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( {cs137}, duration_spec(dt),
                                                                SpecUtils::time_point_t{}, true );
  BOOST_REQUIRE( res.error.empty() );
  BOOST_REQUIRE_EQUAL( res.sources.size(), 1 );
  const BackDecaySourceResult &r = res.sources[0];
  const double factor = std::pow( 2.0, dt/hl );
  BOOST_CHECK( r.status == BackDecayStatus::Ok );
  BOOST_CHECK_CLOSE( r.activity, act*factor, 1.0E-8 );
  BOOST_CHECK_CLOSE( r.ratio, factor, 1.0E-8 );
  BOOST_CHECK_CLOSE( r.half_lives, dt/hl, 1.0E-10 );
  BOOST_REQUIRE( r.activity_uncertainty.has_value() );
  BOOST_CHECK_CLOSE( *r.activity_uncertainty, unc*factor, 1.0E-6 );
  BOOST_CHECK( r.note.empty() );
  BOOST_CHECK( res.warnings.empty() );
  BOOST_CHECK( SpecUtils::is_special(res.ref_time) );  //no measurement start time given

  // Later than the measurement: decays forward.
  const BackDecayResults fwd = ShieldingSourceFitCalc::compute_back_decay( {cs137}, duration_spec(-dt),
                                                                SpecUtils::time_point_t{}, true );
  BOOST_REQUIRE( fwd.error.empty() );
  BOOST_CHECK_CLOSE( fwd.sources[0].activity, act/factor, 1.0E-8 );
  BOOST_CHECK_CLOSE( *fwd.sources[0].activity_uncertainty, unc/factor, 1.0E-8 );
  BOOST_CHECK( fwd.sources[0].half_lives < 0.0 );

  // At the measurement time.
  const BackDecayResults same = ShieldingSourceFitCalc::compute_back_decay( {cs137}, duration_spec(0.0),
                                                                SpecUtils::time_point_t{}, true );
  BOOST_CHECK_CLOSE( same.sources[0].activity, act, 1.0E-10 );

  // A date/time: two days before the measurement start.
  const SpecUtils::time_point_t meas_start = parse_time( "2020-01-01T00:00:00" );
  BackDecaySpec date_spec;
  date_spec.enabled = true;
  date_spec.mode = BackDecaySpec::Mode::DateTime;
  date_spec.date = parse_time( "2019-12-30T00:00:00" );

  const SourceFitDef br82 = make_src( "Br82", act, 0.0 );
  const BackDecayResults dated = ShieldingSourceFitCalc::compute_back_decay( {br82}, date_spec, meas_start, true );
  BOOST_REQUIRE( dated.error.empty() );
  BOOST_CHECK_CLOSE( dated.time_before_meas, 2.0*PhysicalUnits::day, 1.0E-10 );
  BOOST_CHECK( dated.ref_time == date_spec.date );
  BOOST_CHECK_CLOSE( dated.sources[0].activity, act*std::pow(2.0, 2.0*PhysicalUnits::day/br82.nuclide->halfLife), 1.0E-8 );
  BOOST_CHECK( !dated.sources[0].activity_uncertainty.has_value() );  //no input uncertainty

  // A duration gives the corresponding reference time.
  const BackDecayResults dur = ShieldingSourceFitCalc::compute_back_decay( {br82}, duration_spec(2.0*PhysicalUnits::day),
                                                                meas_start, true );
  BOOST_CHECK( dur.ref_time == date_spec.date );
  BOOST_CHECK_EQUAL( ShieldingSourceFitCalc::back_decay_time_str(dur.ref_time), "2019-12-30 00:00:00" );

  // A date/time can not be used without the measurement start time; that is an error for the
  //  back-decay as a whole, but not an exception.
  const BackDecayResults no_start = ShieldingSourceFitCalc::compute_back_decay( {br82}, date_spec,
                                                                     SpecUtils::time_point_t{}, true );
  BOOST_CHECK( !no_start.error.empty() );
  BOOST_CHECK( no_start.sources.empty() );
}//BOOST_AUTO_TEST_CASE( SingleNuclideScaling )


// What can be computed is; a nuclide too many half-lives back gets an error of its own.
BOOST_AUTO_TEST_CASE( PerNuclideLimits )
{
  set_data_dir();

  // The activation products in the in-situ example data.  30 days is ~20 half-lives of Br82, but
  //  ~48 of Na24 and ~58 of K42.
  const vector<SourceFitDef> srcs{ make_src("Br82", 1000.0, 0.0, 10.0), make_src("Na24", 50.0, 0.0, 5.0),
                                   make_src("K42", 20.0, 0.0, 4.0) };
  const double dt = 30.0*PhysicalUnits::day;

  BackDecayResults res;
  BOOST_REQUIRE_NO_THROW( res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                   SpecUtils::time_point_t{}, true ) );
  BOOST_REQUIRE( res.error.empty() );

  const BackDecaySourceResult &br82 = result_for( res, "Br82" );
  BOOST_CHECK( br82.status == BackDecayStatus::Ok );
  BOOST_CHECK_CLOSE( br82.activity, 1000.0*std::pow(2.0, dt/nuc("Br82")->halfLife), 1.0E-6 );
  BOOST_CHECK( br82.activity_uncertainty.has_value() );

  for( const char *symbol : { "Na24", "K42" } )
  {
    const BackDecaySourceResult &r = result_for( res, symbol );
    BOOST_CHECK_MESSAGE( r.status == BackDecayStatus::DecayedAway, symbol << " should not be determinable" );
    BOOST_CHECK_EQUAL( r.activity, 0.0 );
    BOOST_CHECK( !r.activity_uncertainty.has_value() );
    BOOST_CHECK_MESSAGE( r.note.find("half-lives") != string::npos, "Note: " << r.note );
  }

  // ~100 half-lives of Na24, but only ~42 of Br82.
  const double dt100 = 100.0*nuc("Na24")->halfLife;
  const BackDecayResults res100 = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt100),
                                                                     SpecUtils::time_point_t{}, true );
  BOOST_CHECK( result_for(res100, "Br82").status == BackDecayStatus::Ok );
  BOOST_CHECK( result_for(res100, "Na24").status == BackDecayStatus::DecayedAway );
  BOOST_CHECK_CLOSE( result_for(res100, "Na24").half_lives, 100.0, 1.0E-8 );
}//BOOST_AUTO_TEST_CASE( PerNuclideLimits )


// With fitted ancestors at least as old as the back-decay time, each source just decays back on its own.
BOOST_AUTO_TEST_CASE( CoupledAgesLongerThanBackDecay )
{
  set_data_dir();

  const double dt = 2.0*PhysicalUnits::day;
  const double te_now = 1.0E4, i_now = 300.0;
  const vector<SourceFitDef> srcs{ make_src("Te132", te_now, 5.0*PhysicalUnits::day),
                                   make_src("I132", i_now, 0.0) };

  const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                SpecUtils::time_point_t{}, true );
  BOOST_REQUIRE( res.error.empty() );
  const BackDecaySourceResult &te = result_for( res, "Te132" ), &i132 = result_for( res, "I132" );
  BOOST_REQUIRE( te.status == BackDecayStatus::Ok );
  BOOST_REQUIRE( i132.status == BackDecayStatus::Ok );
  BOOST_CHECK_CLOSE( te.activity, te_now*std::pow(2.0, dt/nuc("Te132")->halfLife), 1.0E-8 );
  BOOST_CHECK_CLOSE( i132.activity, i_now*std::pow(2.0, dt/nuc("I132")->halfLife), 1.0E-5 );
  BOOST_CHECK( i132.note.empty() );  //Te132's age covers the back-decay time

  // Decaying that reference-time state forward reproduces the measurement: Te132 aged 3 days at
  //  the reference time, plus pure I132.
  SandiaDecay::NuclideMixture ref_state;
  ref_state.addAgedNuclideByActivity( nuc("Te132"), te.activity, 3.0*PhysicalUnits::day );
  ref_state.addNuclideByActivity( nuc("I132"), i132.activity );
  const double te_fwd = static_cast<double>( ref_state.activity( dt, nuc("Te132") ) );
  const double i_fwd = static_cast<double>( ref_state.activity( dt, nuc("I132") ) );
  BOOST_CHECK_CLOSE( te_fwd, te_now, 1.0E-6 );
  BOOST_CHECK_CLOSE( i_fwd, i_now + ingrowth_at_age("Te132", te_now, "I132", 5.0*PhysicalUnits::day), 1.0E-5 );
}//BOOST_AUTO_TEST_CASE( CoupledAgesLongerThanBackDecay )


// The case prone to issues: a fission product that is both produced directly, and fed by another
//  fitted fission product.
BOOST_AUTO_TEST_CASE( FissionReleaseCoupledNuclides )
{
  set_data_dir();

  struct Chain{ const char *parent; double parent_ref; const char *progeny; double progeny_ref; double dt; };
  const vector<Chain> chains{
    { "Te132", 2.0E4, "I132",  5.0E3, 6.0*PhysicalUnits::hour },
    { "Te132", 2.0E4, "I132",  5.0E3, 20.0*PhysicalUnits::hour },
    { "Mo99",  1.0E4, "Tc99m", 3.0E3, 12.0*PhysicalUnits::hour },
    { "Ba140", 1.0E4, "La140", 2.0E3, 3.0*PhysicalUnits::day },
    { "Zr95",  1.0E4, "Nb95",  1.0E3, 20.0*PhysicalUnits::day },
    { "Ce144", 1.0E4, "Pr144", 4.0E3, 3.0*PhysicalUnits::hour }
  };

  for( const Chain &c : chains )
  {
    // Fit parent age: the full time since release (the fit's split matches the truth), no age (all
    //  of the progeny attributed to its own source), and in-between.
    for( const double age_frac : { 1.0, 0.0, 0.5, 3.0 } )
      check_release( c.parent, c.parent_ref, c.progeny, c.progeny_ref, c.dt, age_frac*c.dt );
  }//for( const Chain &c : chains )

  // Over 1.5 days Pr144 (17 minutes) is ~125 half-lives: it can not be determined, but Ce144 can.
  {
    const double dt = 1.5*PhysicalUnits::day;
    const vector<SourceFitDef> srcs{ make_src("Ce144", 1.0E4, 0.0), make_src("Pr144", 1.0E4, 0.0) };
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    BOOST_CHECK( result_for(res, "Ce144").status == BackDecayStatus::Ok );
    BOOST_CHECK( result_for(res, "Pr144").status == BackDecayStatus::DecayedAway );
    BOOST_CHECK_CLOSE( result_for(res, "Ce144").activity, 1.0E4*std::pow(2.0, dt/nuc("Ce144")->halfLife), 1.0E-8 );
  }

  // Zr95 -> Nb95m -> Nb95: over 200 days Nb95m (3.6 d) is ~55 half-lives, and its unknown amount
  //  then would have grown into Nb95, so Nb95 can not be determined either; Zr95 can.
  {
    const double dt = 200.0*PhysicalUnits::day;
    const vector<SourceFitDef> srcs{ make_src("Zr95", 1.0E4, 0.0), make_src("Nb95m", 50.0, 0.0),
                                     make_src("Nb95", 5.0E3, 0.0) };
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    BOOST_REQUIRE( res.error.empty() );
    BOOST_CHECK( result_for(res, "Zr95").status == BackDecayStatus::Ok );
    BOOST_CHECK( result_for(res, "Nb95m").status == BackDecayStatus::DecayedAway );
    BOOST_CHECK( result_for(res, "Nb95").status == BackDecayStatus::DependsOnUndetermined );
    BOOST_CHECK_MESSAGE( result_for(res, "Nb95").note.find("Nb95m") != string::npos,
                         "Note: " << result_for(res, "Nb95").note );
    BOOST_CHECK_CLOSE( result_for(res, "Zr95").activity, 1.0E4*std::pow(2.0, dt/nuc("Zr95")->halfLife), 1.0E-8 );
  }

  // The three of them over 10 days is fine, and decaying back then forward is consistent.
  {
    const double dt = 10.0*PhysicalUnits::day;
    SandiaDecay::NuclideMixture truth;
    truth.addNuclideByActivity( nuc("Zr95"), 1.0E4 );
    truth.addNuclideByActivity( nuc("Nb95m"), 200.0 );
    truth.addNuclideByActivity( nuc("Nb95"), 3.0E3 );

    // All three sources at age zero: each fitted activity is the nuclide's total.
    const vector<SourceFitDef> srcs{ make_src("Zr95", static_cast<double>(truth.activity(dt, nuc("Zr95"))), 0.0),
                                     make_src("Nb95m", static_cast<double>(truth.activity(dt, nuc("Nb95m"))), 0.0),
                                     make_src("Nb95", static_cast<double>(truth.activity(dt, nuc("Nb95"))), 0.0) };
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    BOOST_CHECK_CLOSE( result_for(res, "Zr95").activity, 1.0E4, 1.0E-6 );
    BOOST_CHECK_CLOSE( result_for(res, "Nb95m").activity, 200.0, 1.0E-4 );
    BOOST_CHECK_CLOSE( result_for(res, "Nb95").activity, 3.0E3, 1.0E-4 );
  }
}//BOOST_AUTO_TEST_CASE( FissionReleaseCoupledNuclides )


// Without interference correction a fitted activity is the nuclide's whole activity, so in-growth
//  from a fitted parent must be subtracted even when the parent's age covers the back-decay time.
BOOST_AUTO_TEST_CASE( InterferenceCorrectionOff )
{
  set_data_dir();

  const double dt = 6.0*PhysicalUnits::hour;
  const double te_ref = 2.0E4, i_ref = 5.0E3;
  SandiaDecay::NuclideMixture truth;
  truth.addNuclideByActivity( nuc("Te132"), te_ref );
  truth.addNuclideByActivity( nuc("I132"), i_ref );
  const double te_now = static_cast<double>( truth.activity( dt, nuc("Te132") ) );
  const double i_total_now = static_cast<double>( truth.activity( dt, nuc("I132") ) );

  const vector<SourceFitDef> srcs{ make_src("Te132", te_now, dt), make_src("I132", i_total_now, 0.0) };

  const BackDecayResults off = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                SpecUtils::time_point_t{}, false );
  BOOST_CHECK_CLOSE( result_for(off, "I132").activity, i_ref, 1.0E-5 );
  bool warned = false;
  for( const string &w : off.warnings )
    warned |= (w.find("Interference correction is off") != string::npos);
  BOOST_CHECK( warned );

  // The same input with interference correction on would mean Te132's in-growth is in addition to
  //  the 'I132' activity, so the result is different (and not the truth).
  const BackDecayResults on = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                               SpecUtils::time_point_t{}, true );
  BOOST_CHECK( fabs(result_for(on, "I132").activity - i_ref) > 0.5*i_ref );
}//BOOST_AUTO_TEST_CASE( InterferenceCorrectionOff )


// Uncertainty propagation, and a negative (unclamped) result.
BOOST_AUTO_TEST_CASE( CoupledUncertaintyAndNegative )
{
  set_data_dir();

  const double dt = 6.0*PhysicalUnits::hour;
  const double te_now = 1.0E4, i_now = 1.6E4, te_unc = 200.0, i_unc = 150.0;
  const double s_i = std::pow( 2.0, -dt/nuc("I132")->halfLife );  //surviving fraction of I132

  // Te132 with an age covering the back-decay time: the I132 result does not depend on Te132's activity.
  {
    const vector<SourceFitDef> srcs{ make_src("Te132", te_now, 2.0*PhysicalUnits::day, te_unc),
                                     make_src("I132", 500.0, 0.0, i_unc) };
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    BOOST_REQUIRE( result_for(res, "I132").activity_uncertainty.has_value() );
    BOOST_CHECK_CLOSE( *result_for(res, "I132").activity_uncertainty, i_unc/s_i, 1.0E-4 );
  }

  // Te132 pure at the reference time (age zero): I132 then = (A_I - g*A_Te)/s, where g is the I132
  //  that unit Te132 at the measurement has grown in since the reference time.
  {
    const vector<SourceFitDef> srcs{ make_src("Te132", te_now, 0.0, te_unc),
                                     make_src("I132", i_now, 0.0, i_unc) };
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    const double te_ref = te_now * std::pow( 2.0, dt/nuc("Te132")->halfLife );
    const double g = ingrowth( "Te132", te_ref, "I132", dt ) / te_now;
    const double expected = (i_now - g*te_now) / s_i;
    const double expected_unc = std::sqrt( i_unc*i_unc + g*g*te_unc*te_unc ) / s_i;

    const BackDecaySourceResult &r = result_for( res, "I132" );
    BOOST_REQUIRE( r.status == BackDecayStatus::Ok );
    BOOST_CHECK_CLOSE( r.activity, expected, 1.0E-4 );
    BOOST_REQUIRE( r.activity_uncertainty.has_value() );
    BOOST_CHECK_CLOSE( *r.activity_uncertainty, expected_unc, 1.0E-3 );

    bool warned_corr = false;
    for( const string &w : res.warnings )
      warned_corr |= (w.find("neglect correlations") != string::npos);
    BOOST_CHECK( warned_corr );

    // Less I132 than Te132 grows in: a negative value, with a note, rather than being clamped.
    const double too_little = 0.5*g*te_now;
    const vector<SourceFitDef> neg_srcs{ make_src("Te132", te_now, 0.0, te_unc),
                                         make_src("I132", too_little, 0.0, i_unc) };
    const BackDecayResults neg = ShieldingSourceFitCalc::compute_back_decay( neg_srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    const BackDecaySourceResult &nr = result_for( neg, "I132" );
    BOOST_REQUIRE( nr.status == BackDecayStatus::Ok );
    BOOST_CHECK( nr.activity < 0.0 );
    BOOST_CHECK_CLOSE( nr.activity, (too_little - g*te_now)/s_i, 1.0E-4 );
    BOOST_CHECK_MESSAGE( nr.note.find("In-growth from Te132") != string::npos, "Note: " << nr.note );
    BOOST_CHECK_MESSAGE( nr.note.find("inconsistent") != string::npos, "Note: " << nr.note );
  }
}//BOOST_AUTO_TEST_CASE( CoupledUncertaintyAndNegative )


BOOST_AUTO_TEST_CASE( Serialization )
{
  set_data_dir();

  // Not used: nothing is written, so the XML is exactly as before this option existed.
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *base = doc.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
    doc.append_node( base );
    ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
    options.serialize( base );
    BOOST_CHECK( !base->first_node("BackDecay") );

    ShieldingSourceFitCalc::ShieldingSourceFitOptions read_back;
    read_back.deSerialize( base );
    BOOST_CHECK( read_back.back_decay.is_default() );
  }

  // Round trip both modes - both values are kept, whichever mode is active.
  for( const bool use_date : { false, true } )
  {
    for( const bool enabled : { false, true } )
    {
      ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
      options.back_decay.enabled = enabled;
      options.back_decay.mode = use_date ? BackDecaySpec::Mode::DateTime : BackDecaySpec::Mode::DurationBeforeMeas;
      options.back_decay.duration = 3.25*PhysicalUnits::day + 0.123456*PhysicalUnits::second;
      options.back_decay.date = parse_time( "2019-12-28T12:34:56.789" );

      rapidxml::xml_document<char> doc;
      rapidxml::xml_node<char> *base = doc.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
      doc.append_node( base );
      options.serialize( base );
      BOOST_REQUIRE( base->first_node("BackDecay") );

      // Through text, as it would be in a file.
      string xml;
      rapidxml::print( std::back_inserter(xml), doc, 0 );
      rapidxml::xml_document<char> read_doc;
      read_doc.parse<rapidxml::parse_non_destructive>( &xml[0] );

      ShieldingSourceFitCalc::ShieldingSourceFitOptions read_back;
      BOOST_REQUIRE_NO_THROW( read_back.deSerialize( read_doc.first_node() ) );
      BOOST_CHECK( read_back.back_decay == options.back_decay );
      BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( options, read_back ) );
    }//for( enabled )
  }//for( use_date )

  // The <ShieldingSourceFit> version is only bumped when the XML actually has back-decay in it.
  for( const bool use_back_decay : { false, true } )
  {
    auto config = make_shared<GammaInteractionCalc::ShieldSourceConfig>();
    config->distance = 1.0*PhysicalUnits::meter;
    config->geometry = GammaInteractionCalc::GeometryType::Spherical;
    if( use_back_decay )
      config->options.back_decay = duration_spec( 1.0*PhysicalUnits::day );

    ShieldingSourceDisplay::ShieldingSourceDisplayState state;
    state.config = config;

    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *node = nullptr;
    BOOST_REQUIRE_NO_THROW( node = state.serialize( &doc ) );
    BOOST_REQUIRE( node );
    const string version = SpecUtils::xml_value_str( node->first_attribute("version") );
    BOOST_CHECK_EQUAL( version, (use_back_decay ? "0.2" : "0.1") );
  }//for( use_back_decay )
}//BOOST_AUTO_TEST_CASE( Serialization )


namespace
{
  /** The decay transition, and index of the gamma in it, of `emitter` closest to `energy` keV. */
  pair<const SandiaDecay::Transition *, int> find_gamma( const SandiaDecay::Nuclide *emitter, const double energy )
  {
    const SandiaDecay::Transition *best = nullptr;
    int best_index = -1;
    double best_diff = 1.0;
    for( const SandiaDecay::Transition *trans : emitter->decaysToChildren )
    {
      for( size_t i = 0; i < trans->products.size(); ++i )
      {
        const SandiaDecay::RadParticle &part = trans->products[i];
        if( (part.type == SandiaDecay::GammaParticle) && (fabs(part.energy - energy) < best_diff) )
        {
          best = trans;
          best_index = static_cast<int>( i );
          best_diff = fabs( part.energy - energy );
        }
      }//for( loop over products )
    }//for( loop over transitions )

    BOOST_REQUIRE_MESSAGE( best, "No " << energy << " keV gamma from " << emitter->symbol );
    return { best, best_index };
  }//find_gamma(...)


  /** Counts in `live_time` from every photon of `truth` at `time` within `window` keV of `mean`,
   detected with `detector` at `distance` - i.e., what the fit predicts for a peak at `mean` whose
   sigma times the cluster width is `window`. */
  double truth_counts( const SandiaDecay::NuclideMixture &truth, const double time, const double mean,
                       const double window, const double live_time,
                       const DetectorPeakResponse &detector, const double distance )
  {
    double rate = 0.0;
    for( const SandiaDecay::EnergyRatePair &erp : truth.photons( time ) )
    {
      if( fabs(erp.energy - mean) < window )
        rate += erp.numPerSecond * detector.efficiency( erp.energy, distance );
    }
    return rate * live_time / PhysicalUnits::second;
  }//truth_counts(...)


  struct SynthFit
  {
    std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results;
  };


  /** Fits Te132, I132, Ba140 and La140 to peaks synthesized from a "release" `dt` before the
   measurement, with the parents' ages being `parent_age`, and back-decay to the release enabled.

   If `background_counts` is positive, a background peak with that many counts (in the background
   live time) is put under the 667.7 keV I132 peak, and also added to the foreground peak.
   */
  SynthFit synthetic_fit( const double dt, const double parent_age, const double background_counts,
                          const bool background_peak_subtract )
  {
    const double live_time = 1000.0*PhysicalUnits::second;
    const double back_live_time = 2000.0*PhysicalUnits::second;
    const double distance = 100.0*PhysicalUnits::cm;
    const double sigma = 0.3;
    const double cluster_sigma = 1.25;

    // The release: parents, and each progeny's independent fission yield.
    SandiaDecay::NuclideMixture truth;
    truth.addNuclideByActivity( nuc("Te132"), 2.0E4*PhysicalUnits::bq );
    truth.addNuclideByActivity( nuc("I132"), 5.0E3*PhysicalUnits::bq );
    truth.addNuclideByActivity( nuc("Ba140"), 1.0E4*PhysicalUnits::bq );
    truth.addNuclideByActivity( nuc("La140"), 2.0E3*PhysicalUnits::bq );

    // A (nearly) 100% efficient detector.
    auto detector = make_shared<DetectorPeakResponse>();
    detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance, 5*PhysicalUnits::cm, PhysicalUnits::keV,
                                       0, 3000*PhysicalUnits::keV,
                                       DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

    auto foreground = make_shared<SpecUtils::Measurement>();
    auto spec = make_shared<vector<float>>( vector<float>{0.0f, 1.0f, 5.0f, 2.0f, 3.5f} );
    foreground->set_gamma_counts( spec, static_cast<float>(live_time), static_cast<float>(live_time) );
    foreground->set_start_time( parse_time("2020-01-01T00:00:00") );

    struct PeakSpec{ double energy; const char *parent; const char *emitter; };
    const vector<PeakSpec> peak_specs{
      { 228.16, "Te132", "Te132" },
      { 667.71, "I132", "I132" }, { 772.60, "I132", "I132" }, { 954.55, "I132", "I132" },
      { 537.26, "Ba140", "Ba140" }, { 162.66, "Ba140", "Ba140" },
      { 1596.21, "La140", "La140" }, { 487.02, "La140", "La140" }, { 815.77, "La140", "La140" }
    };

    std::deque<std::shared_ptr<const PeakDef>> fore_peaks;
    auto back_peaks = make_shared<std::deque<std::shared_ptr<const PeakDef>>>();
    for( const PeakSpec &ps : peak_specs )
    {
      double area = truth_counts( truth, dt, ps.energy, cluster_sigma*sigma, live_time, *detector, distance );
      if( (background_counts > 0.0) && (fabs(ps.energy - 667.71) < 0.01) )
      {
        area += background_counts * live_time / back_live_time;

        auto back_peak = make_shared<PeakDef>( ps.energy, sigma, background_counts );
        back_peak->setPeakAreaUncert( sqrt(background_counts) );
        back_peaks->push_back( back_peak );
      }//if( this peak has background under it )

      auto peak = make_shared<PeakDef>( ps.energy, sigma, area );
      peak->setPeakAreaUncert( sqrt(area) );
      const pair<const SandiaDecay::Transition *, int> gamma = find_gamma( nuc(ps.emitter), ps.energy );
      peak->setNuclearTransition( nuc(ps.parent), gamma.first, gamma.second, PeakDef::SourceGammaType::NormalGamma );
      peak->useForShieldingSourceFit( true );
      fore_peaks.push_back( peak );
    }//for( const PeakSpec &ps : peak_specs )

    shared_ptr<SpecUtils::Measurement> background;
    if( background_counts > 0.0 )
    {
      background = make_shared<SpecUtils::Measurement>();
      auto back_spec = make_shared<vector<float>>( vector<float>{0.0f, 1.0f, 5.0f, 2.0f, 3.5f} );
      background->set_gamma_counts( back_spec, static_cast<float>(back_live_time), static_cast<float>(back_live_time) );
    }

    // The progeny get their default age, as in the GUI - which for La140 brings its microsecond
    //  Ce140m progeny into equilibrium, as it is in the truth spectrum.
    vector<SourceFitDef> srcs;
    for( const char *symbol : { "Te132", "I132", "Ba140", "La140" } )
    {
      const bool is_parent = (string(symbol) == "Te132") || (string(symbol) == "Ba140");
      const double age = is_parent ? parent_age : PeakDef::defaultDecayTime( nuc(symbol) );
      srcs.push_back( make_src( symbol, 1.0E3*PhysicalUnits::bq, age ) );
    }

    ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
    options.attenuate_for_air = false;
    options.multiple_nucs_contribute_to_peaks = true;
    options.photopeak_cluster_sigma = cluster_sigma;
    options.background_peak_subtract = background_peak_subtract;
    options.back_decay = duration_spec( dt );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    chi_input.config.distance = distance;
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
    chi_input.config.sources = srcs;
    chi_input.config.options = options;
    chi_input.detector = detector;
    chi_input.foreground = foreground;
    chi_input.background = background;
    chi_input.foreground_peaks = fore_peaks;
    chi_input.background_peaks = background ? back_peaks : nullptr;
    chi_input.supplemental_options.compute = false;

    const pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ShieldingSourceFitCalc::FitParameters> fcn_pars
                  = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
    const auto input_pars = make_shared<ShieldingSourceFitCalc::FitParameters>( fcn_pars.second );
    auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();

    SynthFit answer;
    answer.results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
    ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, input_pars, progress, [](){}, answer.results, [](){} );

    return answer;
  }//synthetic_fit(...)
}//namespace


// Through an actual fit: activities at the release time come back to the truth, whether the
//  parent ages in the fit were set to the time since release or left at zero.
BOOST_AUTO_TEST_CASE( SyntheticFitEndToEnd )
{
  set_data_dir();

  const double dt = 6.0*PhysicalUnits::hour;
  const map<string,double> release{ {"Te132", 2.0E4}, {"I132", 5.0E3}, {"Ba140", 1.0E4}, {"La140", 2.0E3} };

  struct Case{ double parent_age; double back_counts; bool back_sub; };
  for( const Case &c : { Case{dt, 0.0, false}, Case{0.0, 0.0, false}, Case{dt, 5.0E5, true} } )
  {
    const SynthFit fit = synthetic_fit( dt, c.parent_age, c.back_counts, c.back_sub );
    const ShieldingSourceFitCalc::ModelFitResults &res = *fit.results;
    BOOST_REQUIRE( res.successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );
    BOOST_REQUIRE_MESSAGE( res.back_decay.has_value(), "Fit did not compute back-decay" );
    BOOST_REQUIRE( res.back_decay->error.empty() );
    BOOST_CHECK( res.foreground_start_time == parse_time("2020-01-01T00:00:00") );
    BOOST_CHECK( res.back_decay->ref_time == parse_time("2019-12-31T18:00:00") );

    // The synthetic data is exactly what the model predicts.
    BOOST_CHECK_MESSAGE( res.chi2 < 1.0E-3, "Synthetic fit chi2=" << res.chi2 );


    for( const pair<const string,double> &truth : release )
    {
      const BackDecaySourceResult &r = result_for( *res.back_decay, truth.first.c_str() );
      BOOST_REQUIRE( r.status == BackDecayStatus::Ok );
      BOOST_REQUIRE( r.activity_uncertainty.has_value() );
      BOOST_CHECK_MESSAGE( fabs(r.activity - truth.second) < 1.0E-3*truth.second,
                           truth.first << " (parent age " << c.parent_age/PhysicalUnits::hour << " h, back counts "
                           << c.back_counts << ") at release: " << r.activity << " +- "
                           << *r.activity_uncertainty << ", expected " << truth.second );
    }//for( loop over truth )
  }//for( loop over cases )

  // Background under an I132 peak, but not subtracted: the fit (and so the back-decay) is off, but
  //  each source is still back-decayed consistently - with the parents' ages covering the back-decay
  //  time, each is just its own activity decayed back.
  {
    const SynthFit fit = synthetic_fit( dt, dt, 5.0E5, false );
    const ShieldingSourceFitCalc::ModelFitResults &res = *fit.results;
    BOOST_REQUIRE( res.successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );
    BOOST_REQUIRE( res.back_decay.has_value() );
    for( const SourceFitDef &src : res.fit_src_info )
    {
      const BackDecaySourceResult * const r = res.back_decay->source( src.nuclide );
      BOOST_REQUIRE( r );
      BOOST_CHECK_CLOSE( r->activity, src.activity*std::pow(2.0, dt/src.nuclide->halfLife), 1.0E-4 );
    }
    BOOST_CHECK( result_for(*res.back_decay, "I132").activity > 1.05*5.0E3 );
  }
}//BOOST_AUTO_TEST_CASE( SyntheticFitEndToEnd )


namespace
{
  string in_situ_dir()
  {
    set_data_dir();
    const string dir = SpecUtils::append_path( g_test_file_dir, "ActShieldBackDecay" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory(dir), "No test data directory '" << dir << "'" );
    return dir;
  }


  /** The in-situ exemplar (its Activity/Shielding model fits Br82, Na24 and K42, with a fixed
   geometry per-m2 detector response), with its model modified. */
  shared_ptr<SpecMeas> load_exemplar( const std::optional<BackDecaySpec> &back_decay,
                                      const bool background_peak_subtract )
  {
    const string path = SpecUtils::append_path( in_situ_dir(), "exemplar.n42" );
    auto meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( meas->load_file( path, SpecUtils::ParserType::Auto ), "Could not load " << path );

    rapidxml::xml_document<char> * const doc = meas->shieldingSourceModel();
    BOOST_REQUIRE( doc );
    rapidxml::xml_node<char> * const base = doc->first_node( "ShieldingSourceFit" );
    BOOST_REQUIRE( base );

    rapidxml::xml_node<char> * const back_sub = base->first_node( "BackgroundPeakSubtraction" );
    BOOST_REQUIRE( back_sub );
    back_sub->value( background_peak_subtract ? "1" : "0" );

    if( back_decay.has_value() )
      back_decay->serialize( base );

    return meas;
  }//load_exemplar(...)


  vector<string> location_files()
  {
    vector<string> files;
    for( int i = 1; i <= 5; ++i )
      files.push_back( SpecUtils::append_path( in_situ_dir(), "location" + std::to_string(i) + ".n42" ) );
    return files;
  }


  /** The survey files, freshly loaded (a batch fit may change their energy calibration), with the
   start times a survey of the sites over a few days could have: location<i> starts `24 + 7*(i-1)`
   hours after `base` - except location2, whose start time ("22020-01-01T00:00:00Z") is not a valid
   date, and is left that way.  `start_offsets` gets each valid start time, relative to `base`.
   */
  vector<shared_ptr<SpecMeas>> survey_files( const SpecUtils::time_point_t base,
                                             map<string,double> &start_offsets )
  {
    start_offsets.clear();
    vector<shared_ptr<SpecMeas>> files;
    const vector<string> paths = location_files();
    for( size_t i = 0; i < paths.size(); ++i )
    {
      auto meas = make_shared<SpecMeas>();
      BOOST_REQUIRE( meas->load_file( paths[i], SpecUtils::ParserType::Auto ) );
      const string fname = SpecUtils::filename( paths[i] );
      if( fname != "location2.n42" )
      {
        const int hours = 24 + 7*static_cast<int>(i);
        for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
          meas->set_start_time( base + std::chrono::hours(hours), m );
        start_offsets[fname] = hours * PhysicalUnits::hour;
      }//if( a valid start time )
      files.push_back( meas );
    }//for( loop over files )
    return files;
  }//survey_files(...)


  BatchActivity::BatchActivityFitOptions batch_options()
  {
    BatchActivity::BatchActivityFitOptions options;
    options.fit_all_peaks = false;
    options.to_stdout = false;
    // The survey files' energy calibrations differ from the exemplar's by a few keV at 1.5 MeV;
    //  without a refit, the Na24 and K42 lines, and most Br82 lines above 1 MeV, are not fit.
    options.refit_energy_cal = true;
    options.use_exemplar_energy_cal = false;
    options.write_n42_with_results = false;
    options.show_nonfit_peaks = false;
    options.overwrite_output_files = false;
    options.create_csv_output = false;
    options.create_json_output = false;
    options.concatenate_to_n42 = false;
    options.use_existing_background_peaks = false;
    options.use_exemplar_energy_cal_for_background = false;
    options.peak_stat_threshold = 0.0;
    options.peak_hypothesis_threshold = 0.0;
    options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::None;
    options.template_include_dir = "default";
    options.report_templates = { "txt", "html" };
    options.summary_report_templates = { "csv", "html" };
    options.use_bq = true;
    options.hard_background_sub = false;
    return options;
  }//batch_options()


  BatchActivity::BatchActivityFitSummary run_batch( const shared_ptr<SpecMeas> &exemplar,
                                                    const BatchActivity::BatchActivityFitOptions &options,
                                                    vector<shared_ptr<SpecMeas>> cached_files = {} )
  {
    BatchActivity::BatchActivityFitSummary summary;
    BOOST_REQUIRE_NO_THROW( BatchActivity::fit_activities_in_files( "exemplar.n42", exemplar, {1},
                                                    location_files(), cached_files, options, &summary ) );
    BOOST_REQUIRE_EQUAL( summary.file_results.size(), 5 );
    BOOST_REQUIRE_EQUAL( summary.file_json.size(), 5 );
    BOOST_REQUIRE_EQUAL( summary.file_reports.size(), 5 );
    for( const BatchActivity::BatchActivityFitResult &fr : summary.file_results )
      BOOST_REQUIRE_MESSAGE( (fr.m_result_code == BatchActivity::BatchActivityFitResult::ResultCode::Success)
                             && fr.m_fit_results, fr.m_filename << " fit failed: " << fr.m_error_msg );
    return summary;
  }//run_batch(...)


  bool contains( const vector<string> &strs, const string &substr )
  {
    for( const string &s : strs )
    {
      if( s.find(substr) != string::npos )
        return true;
    }
    return false;
  }


  /** The rows of the summary CSV's "Back-decayed activities" table, by {filename, nuclide} (an
   error row has an empty nuclide); each row is its 14 comma-separated fields, then the quoted note.
   */
  map<pair<string,string>,vector<string>> summary_csv_back_decay_rows( const string &csv )
  {
    map<pair<string,string>,vector<string>> rows;
    const string heading = "Back-decayed activities\n";
    size_t pos = csv.find( heading );
    BOOST_REQUIRE_MESSAGE( pos != string::npos, "The summary CSV has no back-decay table" );
    pos = csv.find( '\n', pos + heading.size() );  //Skip the column headings
    BOOST_REQUIRE( pos != string::npos );
    ++pos;

    while( (pos < csv.size()) && (csv[pos] != '\n') && (csv[pos] != '\r') )
    {
      size_t end = csv.find( '\n', pos );
      if( end == string::npos )
        end = csv.size();
      string line = csv.substr( pos, end - pos );
      pos = end + 1;
      SpecUtils::trim( line );

      vector<string> fields;
      size_t start = 0;
      for( int i = 0; i < 14; ++i )
      {
        const size_t comma = line.find( ',', start );
        BOOST_REQUIRE_MESSAGE( comma != string::npos, "Short back-decay CSV row: " << line );
        string field = line.substr( start, comma - start );
        SpecUtils::trim( field );
        fields.push_back( field );
        start = comma + 1;
      }//for( loop over fields before the note )
      fields.push_back( line.substr(start) );

      BOOST_CHECK_MESSAGE( !rows.count( {fields[0], fields[1]} ), "Duplicate CSV row: " << line );
      rows[{fields[0], fields[1]}] = fields;
    }//while( loop over table rows )

    return rows;
  }//summary_csv_back_decay_rows(...)


  typedef map<string,map<string,BackDecayStatus>> SurveyStatuses;  //{filename: {nuclide: status}}

  /** Checks a batch fit of the survey files.  Each file in `time_before_meas` must have that
   back-decay time, and (as none of the nuclides are related) each nuclide's value must be its
   fitted activity decayed by it, determinable only if enough of it survives that time; every
   other file must have a back-decay error.  The per-file JSON and the summary CSV must have the
   same values.  Returns the status of each fit nuclide.
   */
  SurveyStatuses check_survey_batch( const BatchActivity::BatchActivityFitSummary &summary,
                                     const map<string,double> &time_before_meas )
  {
    SurveyStatuses statuses;
    BOOST_REQUIRE( !summary.summary_reports.empty() );  //CSV, then HTML
    const map<pair<string,string>,vector<string>> csv_rows
                                        = summary_csv_back_decay_rows( summary.summary_reports[0] );
    size_t num_expected_rows = 0;

    for( size_t i = 0; i < summary.file_results.size(); ++i )
    {
      const string fname = SpecUtils::filename( summary.file_results[i].m_filename );
      const ShieldingSourceFitCalc::ModelFitResults &res = *summary.file_results[i].m_fit_results;
      BOOST_REQUIRE_MESSAGE( res.back_decay.has_value(), fname << " has no back-decay results" );
      const BackDecayResults &bd = *res.back_decay;

      const nlohmann::json json = nlohmann::json::parse( summary.file_json[i] );
      BOOST_REQUIRE_MESSAGE( json.contains("BackDecay") && json.contains("Sources"), fname << " JSON" );
      const nlohmann::json &json_bd = json["BackDecay"];
      const nlohmann::json &json_srcs = json["Sources"];

      ++num_expected_rows;  //An error row, or (at least) one row per source

      const auto dt_pos = time_before_meas.find( fname );
      if( dt_pos == end(time_before_meas) )
      {
        BOOST_CHECK_MESSAGE( !bd.error.empty(), fname << " should have a back-decay error" );
        BOOST_CHECK( json_bd["HasError"].get<bool>() );
        const auto row = csv_rows.find( {fname, ""} );
        BOOST_CHECK_MESSAGE( (row != end(csv_rows)) && (row->second[13] == "Error"),
                             "No summary CSV error row for " << fname );
        BOOST_CHECK( contains( summary.file_reports[i], "Back-decay could not be done" ) );
        continue;
      }//if( back-decay can not be done for this file )

      const double dt = dt_pos->second;
      BOOST_REQUIRE_MESSAGE( bd.error.empty(), fname << ": " << bd.error );
      BOOST_CHECK_CLOSE( bd.time_before_meas, dt, 1.0E-6 );
      BOOST_CHECK_CLOSE( json_bd["TimeBeforeMeas_s"].get<double>(), dt/PhysicalUnits::second, 1.0E-6 );
      BOOST_CHECK_EQUAL( json_bd["IsAfterMeas"].get<bool>(), (dt < 0.0) );
      BOOST_CHECK( contains( summary.file_reports[i], (dt < 0.0) ? " after the measurement start"
                                                                 : " before the measurement start" ) );

      BOOST_REQUIRE_EQUAL( bd.sources.size(), res.fit_src_info.size() );
      BOOST_REQUIRE_EQUAL( json_srcs.size(), bd.sources.size() );
      num_expected_rows += bd.sources.size() - 1;

      for( size_t j = 0; j < bd.sources.size(); ++j )
      {
        const BackDecaySourceResult &r = bd.sources[j];
        const SourceFitDef &src = res.fit_src_info[j];
        BOOST_REQUIRE( r.nuclide && (r.nuclide == src.nuclide) );
        const string &symbol = r.nuclide->symbol;
        const string status = ShieldingSourceFitCalc::backDecayStatusName( r.status );
        statuses[fname][symbol] = r.status;

        const double half_lives = dt / r.nuclide->halfLife;
        const double factor = std::pow( 2.0, half_lives );
        const bool expect_ok = (1.0/factor) >= DecayBatchCalc::sm_min_back_decay_surviving_fraction;
        BOOST_CHECK_MESSAGE( (r.status == BackDecayStatus::Ok) == expect_ok,
                             fname << " " << symbol << " is " << status << " at " << half_lives << " half-lives" );
        BOOST_CHECK_CLOSE( r.half_lives, half_lives, 1.0E-6 );
        BOOST_CHECK_CLOSE( r.activity_now, src.activity, 1.0E-9 );

        const nlohmann::json &json_src = json_srcs[j];
        BOOST_CHECK_EQUAL( json_src["Nuclide"].get<string>(), symbol );
        BOOST_CHECK_EQUAL( json_src["BackDecayStatus"].get<string>(), status );
        BOOST_CHECK_EQUAL( json_src["BackDecayIsValid"].get<bool>(), (r.status == BackDecayStatus::Ok) );

        const auto row = csv_rows.find( {fname, symbol} );
        BOOST_REQUIRE_MESSAGE( row != end(csv_rows), "No summary CSV row for " << fname << " " << symbol );
        const vector<string> &fields = row->second;
        BOOST_CHECK_EQUAL( fields[13], status );

        if( r.status == BackDecayStatus::Ok )
        {
          BOOST_CHECK_CLOSE( r.activity, src.activity*factor, 1.0E-6 );
          BOOST_REQUIRE( src.activityUncertainty.has_value() && r.activity_uncertainty.has_value() );
          BOOST_CHECK_CLOSE( *r.activity_uncertainty, (*src.activityUncertainty)*factor, 1.0E-4 );

          BOOST_CHECK_CLOSE( json_src["BackDecayActivity_bq"].get<double>(), r.activity/PhysicalUnits::bq, 1.0E-9 );
          BOOST_CHECK_CLOSE( json_src["BackDecayActivityUncert_bq"].get<double>(),
                             (*r.activity_uncertainty)/PhysicalUnits::bq, 1.0E-9 );
          BOOST_CHECK_CLOSE( std::stod(fields[8]), r.activity/PhysicalUnits::bq, 1.0E-9 );
          BOOST_CHECK_CLOSE( std::stod(fields[9]), (*r.activity_uncertainty)/PhysicalUnits::bq, 1.0E-9 );
        }else
        {
          BOOST_CHECK( fields[8].empty() && fields[9].empty() );
          BOOST_CHECK( contains( summary.file_reports[i], "    " + symbol + ": " + r.note ) );
        }
      }//for( loop over sources )
    }//for( loop over file results )

    BOOST_CHECK_EQUAL( csv_rows.size(), num_expected_rows );
    BOOST_CHECK( contains( summary.summary_reports, "Back-decayed activities" ) );
    BOOST_CHECK( summary.summary_json.find("\"AnyBackDecay\":true") != string::npos );

    return statuses;
  }//check_survey_batch(...)


  /** The time before measurement of each survey file, for a reference time `ref_offset` after the
   survey's base time. */
  map<string,double> survey_time_before_meas( const map<string,double> &start_offsets, const double ref_offset )
  {
    map<string,double> answer;
    for( const pair<const string,double> &name_offset : start_offsets )
      answer[name_offset.first] = name_offset.second - ref_offset;
    return answer;
  }


  /** Writes the reports to `$BACK_DECAY_TEST_DUMP_DIR`, if that is an existing directory. */
  void maybe_dump_reports( const BatchActivity::BatchActivityFitSummary &summary, const string &prefix )
  {
    const char * const dump_dir = std::getenv( "BACK_DECAY_TEST_DUMP_DIR" );
    if( !dump_dir || !SpecUtils::is_directory(dump_dir) )
      return;

    const auto write_file = [dump_dir, &prefix]( const string &name, const string &contents ){
      std::ofstream out( SpecUtils::append_path( dump_dir, prefix + name ).c_str(), std::ios::out | std::ios::binary );
      out << contents;
    };
    for( size_t i = 0; i < summary.summary_reports.size(); ++i )
      write_file( "summary" + std::to_string(i), summary.summary_reports[i] );  // (CSV, then HTML)
    for( size_t i = 0; i < summary.file_reports.size(); ++i )
    {
      for( size_t j = 0; j < summary.file_reports[i].size(); ++j )
        write_file( "file" + std::to_string(i) + "_" + std::to_string(j), summary.file_reports[i][j] );
    }
  }//maybe_dump_reports(...)
}//namespace


// The batch activity fit of the in-situ survey (each site fits Br82, and some of Na24 and K42),
//  quoting the sites at a common time.  The files get start times over a few days; location2's is
//  invalid, so it can only be back-decayed by a duration.
BOOST_AUTO_TEST_CASE( InSituBatchBackDecay )
{
  set_data_dir();

  const SpecUtils::time_point_t base = parse_time( "2019-12-30T00:00:00" );
  const auto date_spec = [base]( const double hours_after_base ) -> BackDecaySpec {
    BackDecaySpec spec;
    spec.enabled = true;
    spec.mode = BackDecaySpec::Mode::DateTime;
    spec.date = base + std::chrono::seconds( static_cast<long long>(3600.0*hours_after_base) );
    return spec;
  };

  // A date before the survey, with and without background peak subtraction, and with hard
  //  background subtraction.
  struct Variant{ bool back_peak_sub; bool use_background; bool hard_sub; };
  for( const Variant &v : { Variant{false, false, false}, Variant{true, true, false}, Variant{false, true, true} } )
  {
    BatchActivity::BatchActivityFitOptions options = batch_options();
    if( v.use_background )
    {
      options.background_subtract_file = SpecUtils::append_path( in_situ_dir(), "Background.n42" );
      options.hard_background_sub = v.hard_sub;
    }

    map<string,double> start_offsets;
    vector<shared_ptr<SpecMeas>> files = survey_files( base, start_offsets );
    const BatchActivity::BatchActivityFitSummary summary
                                  = run_batch( load_exemplar(date_spec(0.0), v.back_peak_sub), options, files );
    maybe_dump_reports( summary, "variant" + std::to_string( v.back_peak_sub ) + std::to_string( v.use_background )
                                 + std::to_string( v.hard_sub ) + "_" );

    const SurveyStatuses statuses = check_survey_batch( summary, survey_time_before_meas(start_offsets, 0.0) );

    // Location3 (the exemplar's site) has all three nuclides, and location1 does not have Na24, so
    //  the summary table has differing nuclides by site.
    BOOST_REQUIRE( statuses.count("location3.n42") && statuses.count("location1.n42") );
    BOOST_CHECK_EQUAL( statuses.at("location3.n42").size(), 3 );
    BOOST_CHECK( !statuses.at("location1.n42").count("Na24") );
  }//for( loop over variants )


  // A reference time in the middle of the survey: the earlier sites are decayed forward to it, the
  //  later ones back.
  {
    const double ref_hours = 40.0;
    map<string,double> start_offsets;
    vector<shared_ptr<SpecMeas>> files = survey_files( base, start_offsets );
    const BatchActivity::BatchActivityFitSummary summary
                              = run_batch( load_exemplar(date_spec(ref_hours), false), batch_options(), files );
    maybe_dump_reports( summary, "midsurvey_" );

    const map<string,double> dts = survey_time_before_meas( start_offsets, ref_hours*PhysicalUnits::hour );
    BOOST_REQUIRE( (dts.at("location1.n42") < 0.0) && (dts.at("location5.n42") > 0.0) );
    check_survey_batch( summary, dts );
  }


  // Weeks before the survey, K42 (12.36 h half-life) can still be determined for the first site
  //  (46.1 half-lives back), but not the later ones (47.3 half-lives back and more), while Na24
  //  (15.0 h) and Br82 (35.3 h) can be for all the sites.
  {
    const double ref_hours = 24.0 - 570.0;
    map<string,double> start_offsets;
    vector<shared_ptr<SpecMeas>> files = survey_files( base, start_offsets );
    const BatchActivity::BatchActivityFitSummary summary
                              = run_batch( load_exemplar(date_spec(ref_hours), false), batch_options(), files );
    maybe_dump_reports( summary, "weeks_" );

    const SurveyStatuses statuses
                  = check_survey_batch( summary, survey_time_before_meas(start_offsets, ref_hours*PhysicalUnits::hour) );

    size_t num_k42_undetermined = 0;
    for( const pair<const string,map<string,BackDecayStatus>> &file : statuses )
    {
      for( const pair<const string,BackDecayStatus> &nuc : file.second )
      {
        if( nuc.first != "K42" )
          BOOST_CHECK_MESSAGE( nuc.second == BackDecayStatus::Ok, file.first << " " << nuc.first );
        else if( file.first == "location1.n42" )
          BOOST_CHECK( nuc.second == BackDecayStatus::Ok );
        else
          num_k42_undetermined += (nuc.second == BackDecayStatus::DecayedAway);
      }
    }//for( loop over files )
    BOOST_CHECK( statuses.count("location1.n42") && statuses.at("location1.n42").count("K42") );
    BOOST_CHECK_GE( num_k42_undetermined, 2 );
  }


  // A duration needs no start time, so location2 gets results too.  Thirty days is ~20 half-lives
  //  of Br82, but Na24 (48) and K42 (58) can not be determined.
  {
    const double dt = 30.0*PhysicalUnits::day;
    const BatchActivity::BatchActivityFitSummary summary
                     = run_batch( load_exemplar(duration_spec(dt), false), batch_options() );
    maybe_dump_reports( summary, "duration_" );

    map<string,double> dts;
    for( const string &path : location_files() )
      dts[SpecUtils::filename(path)] = dt;
    const SurveyStatuses statuses = check_survey_batch( summary, dts );

    set<string> undetermined;
    for( const pair<const string,map<string,BackDecayStatus>> &file : statuses )
    {
      for( const pair<const string,BackDecayStatus> &nuc : file.second )
      {
        BOOST_CHECK_EQUAL( (nuc.second == BackDecayStatus::Ok), (nuc.first == "Br82") );
        if( nuc.second != BackDecayStatus::Ok )
          undetermined.insert( nuc.first );
      }
    }//for( loop over files )
    BOOST_CHECK( undetermined.count("Na24") && undetermined.count("K42") );
  }


  // Not used: the reports and JSON are as they were before back-decay existed.
  {
    const BatchActivity::BatchActivityFitSummary summary = run_batch( load_exemplar({}, false), batch_options() );
    BOOST_CHECK( !contains( summary.summary_reports, "ack-decay" ) );
    // (The test-data directory name has "BackDecay" in it, so look for the JSON keys.)
    for( const char *key : { "\"BackDecay\"", "\"AnyBackDecay\"", "\"HasBackDecay\"", "\"BackDecayEnabled\"" } )
    {
      BOOST_CHECK_MESSAGE( summary.summary_json.find(key) == string::npos, "Summary JSON has " << key );
      for( const string &json : summary.file_json )
        BOOST_CHECK_MESSAGE( json.find(key) == string::npos, "File JSON has " << key );
    }
    for( const vector<string> &rpts : summary.file_reports )
      BOOST_CHECK( !contains( rpts, "ack-decay" ) );
    for( const BatchActivity::BatchActivityFitResult &fr : summary.file_results )
      BOOST_CHECK( !fr.m_fit_results->back_decay.has_value() );
  }
}//BOOST_AUTO_TEST_CASE( InSituBatchBackDecay )


// The bundled per-fit reports render the back-decayed activities the same both ways InterSpec renders
//  them: the GUI calculation log parses the template files with trim_blocks on, batch uses the
//  pre-registered copies parsed with it off (see "Whitespace after {% ... %} tags" in the README).
BOOST_AUTO_TEST_CASE( ReportsRenderBothWays )
{
  set_data_dir();

  const double dt = 6.0*PhysicalUnits::hour;
  const SynthFit fit = synthetic_fit( dt, 0.0, 0.0, false );  //parents age zero -> notes on progeny
  BOOST_REQUIRE( fit.results->back_decay.has_value() );

  nlohmann::json data;
  data["Filename"] = "test.n42";
  data["HasFitResults"] = true;
  data["HasErrorMessage"] = false;
  for( const char *asset : { "D3_JS", "SpectrumChart_JS", "SpectrumChart_CSS",
                             "ShieldingSourceFitPlot_JS", "ShieldingSourceFitPlot_CSS" } )
    data[asset] = "";
  BatchInfoLog::add_exe_info_to_json( data );
  BatchInfoLog::shield_src_fit_results_to_json( *fit.results, nullptr, true, data );

  BOOST_REQUIRE( data.contains("BackDecay") );
  BOOST_CHECK_EQUAL( data["BackDecay"]["RefTime"].get<string>(), "2019-12-31 18:00:00" );
  BOOST_CHECK( data["BackDecay"]["AnyCoupled"].get<bool>() );
  BOOST_CHECK( data["ActShieldFitSetup"]["FitOptions"]["BackDecayEnabled"].get<bool>() );
  for( const nlohmann::json &src : data["Sources"] )
  {
    BOOST_REQUIRE( src.contains("HasBackDecay") );
    BOOST_CHECK( src["BackDecayIsValid"].get<bool>() );
    BOOST_CHECK( src.contains("BackDecayActivity_bq") && src.contains("BackDecayActivityUncert_bq") );
  }

  const string tmplt_dir = BatchInfoLog::default_template_dir();
  BatchPeak::BatchPeakFitOptions options;
  const BatchInfoLog::TemplateRenderType type = BatchInfoLog::TemplateRenderType::ActShieldIndividual;

  inja::Environment batch_env = BatchInfoLog::get_default_inja_env( tmplt_dir );
  inja::Environment gui_env = BatchInfoLog::get_default_inja_env( tmplt_dir );

  const string heading = "18:00:00, 6.0000 h before the measurement start (2020-01-01 00:00:00):";
  for( const bool is_gui : { false, true } )
  {
    string txt, html;
    BOOST_REQUIRE_NO_THROW( txt = is_gui ? gui_env.render_file( "std_fit_log.tmplt.txt", data )
                                         : BatchInfoLog::render_template( "txt", batch_env, type, options, data ) );
    BOOST_REQUIRE_NO_THROW( html = is_gui ? gui_env.render_file( "act_fit.tmplt.html", data )
                                          : BatchInfoLog::render_template( "html", batch_env, type, options, data ) );

    BOOST_CHECK_MESSAGE( txt.find(heading) != string::npos, (is_gui ? "GUI" : "batch") << " txt heading missing:\n" << txt );
    BOOST_CHECK_MESSAGE( html.find(heading) != string::npos, (is_gui ? "GUI" : "batch") << " html heading missing" );

    // Each nuclide on its own line, its note (if any) on the next.
    for( const char *symbol : { "Te132", "I132", "Ba140", "La140" } )
      BOOST_CHECK_MESSAGE( txt.find( string("\n    ") + symbol + ": " ) != string::npos,
                           (is_gui ? "GUI" : "batch") << " txt has no line for " << symbol );
    BOOST_CHECK_MESSAGE( txt.find( "\n      The age of Te132 is shorter than the back-decay time" ) != string::npos,
                         (is_gui ? "GUI" : "batch") << " txt has no note line for I132" );
    BOOST_CHECK( html.find( "BackDecaySummary" ) != string::npos );
  }//for( is_gui )
}//BOOST_AUTO_TEST_CASE( ReportsRenderBothWays )


// Without interference correction a fit activity is all of the nuclide's activity: the other fitted
//  nuclides are taken as pure at the reference time (earlier) or the measurement (later), so the
//  result is continuous through the measurement time, whatever the parent's age.
BOOST_AUTO_TEST_CASE( InterferenceOffContinuityAndForward )
{
  set_data_dir();

  const double te_now = 1.0E4, i_total_now = 1.2E4;
  const double te_age = 2.0*PhysicalUnits::day;  //longer than the back-decay times used below
  const vector<SourceFitDef> srcs{ make_src("Te132", te_now, te_age), make_src("I132", i_total_now, 0.0) };
  const double lambda_i = nuc("I132")->decayConstant(), lambda_te = nuc("Te132")->decayConstant();

  const auto i132_at = [&]( const double dt, const bool interference ) -> double {
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, interference );
    BOOST_REQUIRE( res.error.empty() );
    const BackDecaySourceResult &r = result_for( res, "I132" );
    BOOST_REQUIRE( r.status == BackDecayStatus::Ok );
    return r.activity;
  };

  for( const bool interference : { false, true } )
  {
    const double at_zero = i132_at( 0.0, interference );
    BOOST_CHECK_CLOSE( at_zero, i_total_now, 1.0E-10 );
    BOOST_CHECK_MESSAGE( fabs(i132_at( 1.0*PhysicalUnits::second, interference ) - at_zero) < 1.0E-3*at_zero,
                         "Not continuous just before the measurement (interference=" << interference << ")" );
    BOOST_CHECK_MESSAGE( fabs(i132_at( -1.0*PhysicalUnits::second, interference ) - at_zero) < 1.0E-3*at_zero,
                         "Not continuous just after the measurement (interference=" << interference << ")" );
  }//for( interference )

  // An hour earlier, without interference correction: Te132 taken as pure then, so only the I132 it
  //  grew in during that hour is subtracted.
  const double hour = 1.0*PhysicalUnits::hour;
  const double te_ref = te_now * std::exp( lambda_te * hour );
  const double expected_back = (i_total_now - ingrowth( "Te132", te_ref, "I132", hour )) * std::exp( lambda_i * hour );
  BOOST_CHECK_CLOSE( i132_at( hour, false ), expected_back, 1.0E-5 );

  // An hour later, without interference correction: the I132 Te132 grows in by then is added.
  const double expected_fwd = i_total_now * std::exp( -lambda_i * hour ) + ingrowth( "Te132", te_now, "I132", hour );
  BOOST_CHECK_CLOSE( i132_at( -hour, false ), expected_fwd, 1.0E-6 );

  // ... while with interference correction the later value is just I132's own share decaying.
  BOOST_CHECK_CLOSE( i132_at( -hour, true ), i_total_now * std::exp( -lambda_i * hour ), 1.0E-8 );

  // Back then forward again gets back to the measurement.
  {
    const BackDecayResults back = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(hour),
                                                                 SpecUtils::time_point_t{}, false );
    vector<SourceFitDef> ref_srcs = srcs;
    for( SourceFitDef &src : ref_srcs )
      src.activity = back.source( src.nuclide )->activity;
    const BackDecayResults fwd = ShieldingSourceFitCalc::compute_back_decay( ref_srcs, duration_spec(-hour),
                                                                SpecUtils::time_point_t{}, false );
    BOOST_CHECK_CLOSE( result_for(fwd, "Te132").activity, te_now, 1.0E-6 );
    BOOST_CHECK_CLOSE( result_for(fwd, "I132").activity, i_total_now, 1.0E-5 );
  }

  // The warning says which way the in-growth went.
  const BackDecayResults later = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(-hour),
                                                                SpecUtils::time_point_t{}, false );
  bool warned = false;
  for( const string &w : later.warnings )
    warned |= (w.find("is added") != string::npos);
  BOOST_CHECK( warned );
}//BOOST_AUTO_TEST_CASE( InterferenceOffContinuityAndForward )


// The limit is at a surviving fraction of 2^-47; also very long ages, and far-off reference times.
BOOST_AUTO_TEST_CASE( LimitsAndRanges )
{
  set_data_dir();

  const SourceFitDef na24 = make_src( "Na24", 100.0, 0.0, 1.0 );
  const double hl = na24.nuclide->halfLife;
  const BackDecayResults just_ok = ShieldingSourceFitCalc::compute_back_decay( {na24}, duration_spec(46.9*hl),
                                                                    SpecUtils::time_point_t{}, true );
  BOOST_CHECK( just_ok.sources[0].status == BackDecayStatus::Ok );
  BOOST_CHECK_CLOSE( just_ok.sources[0].activity, 100.0*std::pow(2.0, 46.9), 1.0E-6 );
  const BackDecayResults just_not = ShieldingSourceFitCalc::compute_back_decay( {na24}, duration_spec(47.1*hl),
                                                                     SpecUtils::time_point_t{}, true );
  BOOST_CHECK( just_not.sources[0].status == BackDecayStatus::DecayedAway );

  // A parent age of thousands of its half-lives (only possible by typing it in) is limited, rather
  //  than underflowing and making the parent look undeterminable.
  {
    const vector<SourceFitDef> srcs{ make_src("Te132", 1.0E4, 10.0*PhysicalUnits::year), make_src("I132", 500.0, 0.0) };
    const double dt = 1.0*PhysicalUnits::day;
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                  SpecUtils::time_point_t{}, true );
    const BackDecaySourceResult &te = result_for( res, "Te132" );
    BOOST_CHECK( te.status == BackDecayStatus::Ok );
    BOOST_CHECK_CLOSE( te.activity, 1.0E4*std::pow(2.0, dt/nuc("Te132")->halfLife), 1.0E-8 );
    BOOST_CHECK_MESSAGE( te.note.find("longer than can be evaluated") != string::npos, "Note: " << te.note );
    BOOST_CHECK( result_for(res, "I132").status == BackDecayStatus::Ok );
  }

  // A reference time thousands of years from the measurement can not be written as a date, so is
  //  not given as one - but the activities are still computed from the duration.
  const SpecUtils::time_point_t meas_start = parse_time( "2020-01-01T00:00:00" );
  const SourceFitDef u235 = make_src( "U235", 1.0E3, 0.0 );
  for( const double years : { 5000.0, 50000.0, -8000.0 } )
  {
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( {u235}, duration_spec(years*PhysicalUnits::year),
                                                                  meas_start, true );
    BOOST_CHECK_MESSAGE( SpecUtils::is_special(res.ref_time), years << " years gave a reference time" );
    BOOST_CHECK( ShieldingSourceFitCalc::back_decay_time_str(res.ref_time).empty() );
    BOOST_CHECK( res.sources[0].status == BackDecayStatus::Ok );
  }
  const BackDecayResults recent = ShieldingSourceFitCalc::compute_back_decay( {u235}, duration_spec(100.0*PhysicalUnits::year),
                                                                   meas_start, true );
  BOOST_CHECK( !SpecUtils::is_special(recent.ref_time) );

  // A <BackDecay> node this version can not read is skipped, rather than failing the model.
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *base = doc.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
    doc.append_node( base );
    ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
    options.serialize( base );
    rapidxml::xml_node<char> *bd = doc.allocate_node( rapidxml::node_element, "BackDecay" );
    bd->append_attribute( doc.allocate_attribute( "Enabled", "1" ) );
    bd->append_attribute( doc.allocate_attribute( "Mode", "SomeFutureMode" ) );
    base->append_node( bd );

    ShieldingSourceFitCalc::ShieldingSourceFitOptions read_back;
    BOOST_CHECK_NO_THROW( read_back.deSerialize( base ) );
    BOOST_CHECK( !read_back.back_decay.enabled );
  }
}//BOOST_AUTO_TEST_CASE( LimitsAndRanges )


// The "Show Decay" link starts the decay tool from each source's activity and age at the earlier of
//  the reference time and the measurement; decaying that must give back what was fit, or the
//  back-decayed values: at the measurement each fitted nuclide's total is what the fit sources hold
//  (with interference correction a source of age `a` holds what it grew in that time; without, each
//  fit activity is all of the nuclide), and at a reference time after the measurement each value is
//  the back-decayed one.
BOOST_AUTO_TEST_CASE( DecayToolStartReproducesFit )
{
  set_data_dir();
  const double day = PhysicalUnits::day, bq = PhysicalUnits::bq;
  const SpecUtils::time_point_t meas = parse_time( "2020-01-10T00:00:00" );

  // A release-like Ba140 + La140 (Ba140 of age zero, two days, and older than the back-decay
  //  time), a fission pair with an unrelated old source, and a three-member chain.
  const vector<vector<SourceFitDef>> cases = {
    { make_src("Ba140", 762.0*bq, 0.0),     make_src("La140", 769.5*bq, 0.0) },
    { make_src("Ba140", 762.0*bq, 2.0*day), make_src("La140", 320.3*bq, 0.0) },
    { make_src("Ba140", 762.0*bq, 20.0*day), make_src("La140", 950.0*bq, 0.0) },
    { make_src("Te132", 72.3*bq, 0.0), make_src("I132", 74.5*bq, 0.0), make_src("Cs137", 10.0*bq, 3000.0*day) },
    { make_src("Zr95", 1000.0*bq, 10.0*day), make_src("Nb95m", 5.0*bq, 0.0), make_src("Nb95", 200.0*bq, 0.0) }
  };

  // The decay tool's mixture of one source.
  const auto source_mix = []( const SourceFitNuclideDisplay::DecayToolSource &src ) {
    SandiaDecay::NuclideMixture mix;
    if( src.age > 0.0 )
      mix.addAgedNuclideByActivity( src.nuclide, src.activity, src.age );
    else
      mix.addNuclideByActivity( src.nuclide, src.activity );
    return mix;
  };

  // (NuclideMixture::activity(t, nuc) throws for a nuclide not in the mixture's chains.)
  const auto activity_of = []( const SandiaDecay::NuclideMixture &mix, const double t,
                               const SandiaDecay::Nuclide *nuc ) -> double {
    for( const SandiaDecay::NuclideActivityPair &p : mix.activity( t ) )
    {
      if( p.nuclide == nuc )
        return p.activity;
    }
    return 0.0;
  };

  size_t num_checked = 0;
  for( const vector<SourceFitDef> &srcs : cases )
  {
    for( const bool interference : { true, false } )
    {
      for( const double dt : { 5.0*day, 1.5*day, -3.0*day } )
      {
        const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( srcs, duration_spec(dt),
                                                                               meas, interference );
        BOOST_REQUIRE( res.error.empty() && (res.sources.size() == srcs.size()) );

        const vector<SourceFitNuclideDisplay::DecayToolSource> start
                                                 = SourceFitNuclideDisplay::decayToolStart( res );

        // Undetermined or negative values are left out; then the rest can not reproduce the fit.
        bool all_usable = true;
        for( const BackDecaySourceResult &r : res.sources )
          all_usable = all_usable && (r.status == BackDecayStatus::Ok) && (r.activity > 0.0);
        BOOST_REQUIRE_EQUAL( all_usable, (start.size() == srcs.size()) );
        if( !all_usable )
          continue;
        ++num_checked;

        const string desc = srcs[0].nuclide->symbol + "+... dt=" + std::to_string(dt/day) + " d, interference "
                            + (interference ? "on" : "off");
        vector<SandiaDecay::NuclideMixture> mixes;
        for( const SourceFitNuclideDisplay::DecayToolSource &s : start )
          mixes.push_back( source_mix(s) );

        if( dt >= 0.0 )
        {
          // Started at the reference time; at the measurement, each fitted nuclide's total.
          for( const SourceFitDef &x : srcs )
          {
            double total = 0.0, expected = 0.0;
            for( const SandiaDecay::NuclideMixture &mix : mixes )
              total += activity_of( mix, dt, x.nuclide );
            for( const SourceFitDef &src : srcs )
            {
              if( interference )
              {
                SandiaDecay::NuclideMixture fit_mix;
                if( src.age > 0.0 )
                  fit_mix.addAgedNuclideByActivity( src.nuclide, src.activity, src.age );
                else
                  fit_mix.addNuclideByActivity( src.nuclide, src.activity );
                expected += activity_of( fit_mix, 0.0, x.nuclide );
              }else if( src.nuclide == x.nuclide )
              {
                expected += src.activity;
              }
            }//for( loop over sources )

            BOOST_CHECK_MESSAGE( fabs(total - expected) < 1.0E-6*expected,
                                 desc << ": " << x.nuclide->symbol << " at measurement " << total
                                 << ", but fit " << expected );
          }//for( loop over fitted nuclides )
        }else
        {
          // Started at the measurement; at the (later) reference time, the back-decayed values -
          //  each source's own share with interference correction, or each nuclide's total without.
          for( size_t i = 0; i < srcs.size(); ++i )
          {
            const BackDecaySourceResult &r = res.sources[i];
            double value = 0.0;
            for( size_t j = 0; j < mixes.size(); ++j )
            {
              if( !interference || (j == i) )
                value += activity_of( mixes[j], -dt, r.nuclide );
            }
            BOOST_CHECK_MESSAGE( fabs(value - r.activity) < 1.0E-6*fabs(r.activity),
                                 desc << ": " << r.nuclide->symbol << " at reference time " << value
                                 << ", but back-decayed " << r.activity );
          }//for( loop over sources )
        }//if( dt >= 0.0 ) / else
      }//for( loop over dt )
    }//for( interference )
  }//for( loop over cases )

  BOOST_CHECK_GE( num_checked, 20 );

  // A parent older than the back-decay time keeps its age less that time.
  {
    const BackDecayResults res = ShieldingSourceFitCalc::compute_back_decay( cases[2], duration_spec(5.0*day),
                                                                           meas, true );
    BOOST_CHECK_CLOSE( result_for( res, "Ba140" ).age_at_ref, 15.0*day, 1.0E-9 );
    BOOST_CHECK_EQUAL( result_for( res, "La140" ).age_at_ref, 0.0 );
  }
}//BOOST_AUTO_TEST_CASE( DecayToolStartReproducesFit )
