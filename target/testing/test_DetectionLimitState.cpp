/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
 Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain
 rights in this software.
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

// Must be defined before Windows.h (or any header that includes it) is included; see CLAUDE.md.
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#include <string>
#include <memory>
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE TestDetectionLimitState
#include <boost/test/included/unit_test.hpp>

#include <Wt/WApplication>
#include <Wt/Test/WTestEnvironment>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DetectionLimitTool.h"
#include "InterSpec/DetectionLimitSimple.h"

using namespace std;
using namespace boost::unit_test;


/** The state a detection-limit URI encodes has to keep decoding forever: both tools persist it in
 the database and in QR codes.  These cases pin the grammar and, above all, the migration of the
 retired "assume no signal is present" continuum option - which must never be reinterpreted
 silently, because it measured about 40% coverage where 95% was claimed.

 Deliberately tests the pure helpers rather than the widgets: the migration decision is the part
 that can silently change what a saved number means, and it needs no Wt session to check.
 */

namespace
{
  /** Convenience wrapper: decodes a token and returns the whole outcome. */
  struct Decoded
  {
    bool recognized = false;
    DetectionLimitCalc::DeconContinuumNorm norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
    DetectionLimitCalc::DeconMeasurementModel model
          = DetectionLimitCalc::DeconMeasurementModel::CurrentSpectrum;
    bool migrated = false;
  };

  /** Decodes twice from OPPOSITE seeds and requires the two runs to agree.

   Seeding once cannot detect an output the function never writes: whichever seed you pick happens
   to equal the expected value for some token, and the assertion on it silently goes vacuous.  (It
   did: seeding `BackgroundReference`/`migrated=true` made both of the assertions that
   `DeprecatedTokensMigrateVisibly` exists to make unfalsifiable.)  Two seedings have no such blind
   spot - an unwritten output shows up as a disagreement whatever the token.
   */
  Decoded decode( const string &token )
  {
    Decoded a, b;

    a.norm = DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange;
    a.model = DetectionLimitCalc::DeconMeasurementModel::BackgroundReference;
    a.migrated = true;
    a.recognized = DetectionLimitCalc::decode_continuum_norm_token( token, a.norm, a.model,
                                                                    a.migrated );

    b.norm = DetectionLimitCalc::DeconContinuumNorm::FixedByEdges;
    b.model = DetectionLimitCalc::DeconMeasurementModel::CurrentSpectrum;
    b.migrated = false;
    b.recognized = DetectionLimitCalc::decode_continuum_norm_token( token, b.norm, b.model,
                                                                    b.migrated );

    BOOST_REQUIRE_MESSAGE( a.recognized == b.recognized,
                          "decode of '" << token << "' depended on the seed values" );

    // Only meaningful when recognized: the contract for an unrecognized token is to leave the
    //  outputs untouched, so the two seedings are *supposed* to differ there.
    if( a.recognized )
    {
      BOOST_REQUIRE_MESSAGE( a.norm == b.norm,
                            "decode left `norm` unwritten for '" << token << "'" );
      BOOST_REQUIRE_MESSAGE( a.model == b.model,
                            "decode left `model` unwritten for '" << token << "'" );
      BOOST_REQUIRE_MESSAGE( a.migrated == b.migrated,
                            "decode left `migrated` unwritten for '" << token << "'" );
    }

    return a;
  }
}//namespace


BOOST_AUTO_TEST_CASE( Ver1ContinuumTokensStillDecode )
{
  using namespace DetectionLimitCalc;

  // Simple MDA, VER=1 vocabulary.  "UNKNOWN" and "FIXED" describe continuum treatments that still
  // exist, so they decode unchanged and raise no notice.
  const Decoded unknown = decode( "UNKNOWN" );
  BOOST_REQUIRE( unknown.recognized );
  BOOST_CHECK( unknown.norm == DeconContinuumNorm::Floating );
  BOOST_CHECK( unknown.model == DeconMeasurementModel::CurrentSpectrum );
  BOOST_CHECK( !unknown.migrated );

  const Decoded fixed = decode( "FIXED" );
  BOOST_REQUIRE( fixed.recognized );
  BOOST_CHECK( fixed.norm == DeconContinuumNorm::FixedByEdges );
  BOOST_CHECK( fixed.model == DeconMeasurementModel::CurrentSpectrum );
  BOOST_CHECK( !fixed.migrated );

  // The old decoders matched by prefix, and abbreviated tokens are in the wild.
  BOOST_CHECK( decode("UNK").norm == DeconContinuumNorm::Floating );
  BOOST_CHECK( decode("FIX").norm == DeconContinuumNorm::FixedByEdges );
  BOOST_CHECK( decode("NOS").migrated );

  // Case-insensitive, as every other token in these URIs is.
  BOOST_CHECK( decode("unknown").norm == DeconContinuumNorm::Floating );
  BOOST_CHECK( decode("Fixed").norm == DeconContinuumNorm::FixedByEdges );
}//BOOST_AUTO_TEST_CASE( Ver1ContinuumTokensStillDecode )


BOOST_AUTO_TEST_CASE( DeprecatedTokensMigrateVisibly )
{
  using namespace DetectionLimitCalc;

  // The whole point of the migration: NOSIG (Simple, v1) and FULL (Tool, v1) were the retired
  // "assume no signal is present" option.  It was never a continuum treatment - it asserted the
  // spectrum is signal-free and predicted a future measurement from it - so it becomes a floating
  // continuum on a background-reference measurement, and the caller is told.
  for( const char * const token : { "NOSIG", "FULL", "nosig", "full", "NOS" } )
  {
    const Decoded d = decode( token );
    BOOST_REQUIRE_MESSAGE( d.recognized, "token '" << token << "' was not recognized" );
    BOOST_CHECK_MESSAGE( d.norm == DeconContinuumNorm::Floating,
                        "token '" << token << "' did not migrate to a floating continuum" );
    BOOST_CHECK_MESSAGE( d.model == DeconMeasurementModel::BackgroundReference,
                        "token '" << token << "' did not migrate to a background reference" );
    BOOST_CHECK_MESSAGE( d.migrated,
                        "token '" << token << "' migrated SILENTLY - the notice would not fire" );
  }

  // ... and the treatments that were not retired must NOT raise a notice, or the warning becomes
  // noise the user learns to ignore.
  for( const char * const token : { "UNKNOWN", "FIXED", "FLOAT", "EDGES" } )
    BOOST_CHECK_MESSAGE( !decode(token).migrated, "token '" << token << "' raised a spurious notice" );
}//BOOST_AUTO_TEST_CASE( DeprecatedTokensMigrateVisibly )


BOOST_AUTO_TEST_CASE( Ver2TokensRoundTrip )
{
  using namespace DetectionLimitCalc;

  for( const DeconContinuumNorm norm : { DeconContinuumNorm::Floating,
                                         DeconContinuumNorm::FixedByEdges } )
  {
    const string token = continuum_norm_token( norm );

    const Decoded d = decode( token );
    BOOST_REQUIRE_MESSAGE( d.recognized, "emitted token '" << token << "' does not decode" );
    BOOST_CHECK( d.norm == norm );
    BOOST_CHECK( d.model == DeconMeasurementModel::CurrentSpectrum );
    BOOST_CHECK( !d.migrated );
  }

  BOOST_CHECK_EQUAL( continuum_norm_token( DeconContinuumNorm::Floating ), "FLOAT" );
  BOOST_CHECK_EQUAL( continuum_norm_token( DeconContinuumNorm::FixedByEdges ), "EDGES" );

  // The deprecated value must never be written back out - a state carrying it is migrated on
  // decode, so re-emitting it would resurrect the retired option in a freshly saved URI.
  BOOST_CHECK_THROW( continuum_norm_token( DeconContinuumNorm::FixedByFullRange ), std::exception );
}//BOOST_AUTO_TEST_CASE( Ver2TokensRoundTrip )


BOOST_AUTO_TEST_CASE( UnknownTokenIsRejectedNotGuessed )
{
  using namespace DetectionLimitCalc;

  // An unrecognized token must leave the outputs alone so the caller can log and keep its default,
  // rather than have a typo silently select a different statistical treatment.
  for( const char * const token : { "BOGUS", "", "F", "NO", "FI" } )
  {
    DeconContinuumNorm norm = DeconContinuumNorm::FixedByEdges;
    DeconMeasurementModel model = DeconMeasurementModel::BackgroundReference;
    bool migrated = true;

    const bool ok = decode_continuum_norm_token( token, norm, model, migrated );
    BOOST_CHECK_MESSAGE( !ok, "token '" << token << "' was accepted but should not be" );
    if( !ok )
    {
      BOOST_CHECK( norm == DeconContinuumNorm::FixedByEdges );
      BOOST_CHECK( model == DeconMeasurementModel::BackgroundReference );
      BOOST_CHECK( migrated );
    }
  }
}//BOOST_AUTO_TEST_CASE( UnknownTokenIsRejectedNotGuessed )


BOOST_AUTO_TEST_CASE( ComboIndexMappingIsGuarded )
{
  using namespace DetectionLimitCalc;

  // Simple MDA's continuum combo used to run index 1 -> FixedByFullRange, index 2 -> FixedByEdges,
  // the reverse of the enum for two of three values, with no static_assert guarding it.  Both
  // widgets and both URL decoders now share this mapping.
  BOOST_CHECK_EQUAL( num_selectable_continuum_norms(), 2 );

  BOOST_CHECK( continuum_norm_from_index(0) == DeconContinuumNorm::Floating );
  BOOST_CHECK( continuum_norm_from_index(1) == DeconContinuumNorm::FixedByEdges );
  BOOST_CHECK_EQUAL( index_from_continuum_norm( DeconContinuumNorm::Floating ), 0 );
  BOOST_CHECK_EQUAL( index_from_continuum_norm( DeconContinuumNorm::FixedByEdges ), 1 );

  // The deprecated value has no combo slot; asking for one means a stored state reached a widget
  // without being migrated, which must fail loudly rather than display as something else.
  BOOST_CHECK_THROW( index_from_continuum_norm( DeconContinuumNorm::FixedByFullRange ),
                     std::exception );

  // Out-of-range indices throw rather than silently selecting the first entry.
  BOOST_CHECK_THROW( continuum_norm_from_index(-1), std::exception );
  BOOST_CHECK_THROW( continuum_norm_from_index(2), std::exception );

  BOOST_CHECK( measurement_model_from_index(0) == DeconMeasurementModel::CurrentSpectrum );
  BOOST_CHECK( measurement_model_from_index(1) == DeconMeasurementModel::BackgroundReference );
  BOOST_CHECK_EQUAL( index_from_measurement_model( DeconMeasurementModel::CurrentSpectrum ), 0 );
  BOOST_CHECK_EQUAL( index_from_measurement_model( DeconMeasurementModel::BackgroundReference ), 1 );
  BOOST_CHECK_THROW( measurement_model_from_index(2), std::exception );

  BOOST_CHECK( limit_type_from_index(0) == DeconLimitType::OneSidedUpperLimit );
  BOOST_CHECK( limit_type_from_index(1) == DeconLimitType::CentralInterval );
  BOOST_CHECK_EQUAL( index_from_limit_type( DeconLimitType::OneSidedUpperLimit ), 0 );
  BOOST_CHECK_EQUAL( index_from_limit_type( DeconLimitType::CentralInterval ), 1 );
  BOOST_CHECK_THROW( limit_type_from_index(2), std::exception );
}//BOOST_AUTO_TEST_CASE( ComboIndexMappingIsGuarded )


namespace
{
/** Points InterSpec at its data directory; the round-trip cases below need a real session, and a
 session needs the decay database.  Same `--datadir=` convention the other test targets use.
 */
void set_state_test_data_dir()
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
  }

  SpecUtils::ireplace_all( datadir, "%20", " " );

  if( datadir.empty() )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }
  }//if( datadir.empty() )

  const string decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(decay_file),
                        "sandia.decay.xml not at '" << decay_file << "'" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
}//set_state_test_data_dir()
}//namespace


/** The tools' own URI round trip, which is what undo/redo actually rides on.

 The cases above test the pure token helpers.  They cannot see the thing that breaks undo/redo in
 practice: a control whose value `encodeStateToUrl()` does not write, or writes but
 `handleAppUrl()` does not read.  Either way the two states encode to the same string, so
 `render()`'s `sameAsPrev` check records no undo step and the change simply cannot be undone - and
 nothing fails loudly, which is why it needs a test rather than a review.

 The invariant checked is `encode -> decode -> encode` stability, plus, for each control, that
 changing it changes the URI.  Both need a real widget, hence a Wt test session.
 */
namespace
{
/** A minimal InterSpec session with a spectrum loaded, enough to build the two tools. */
class ToolSessionFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app = nullptr;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_lock;
  InterSpec *m_interspec = nullptr;

  ToolSessionFixture()
  {
    set_state_test_data_dir();

    string app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), ".." );
    app_root = SpecUtils::lexically_normalize_path( app_root );

    m_env.reset( new Wt::Test::WTestEnvironment( "", "", Wt::Application ) );
    m_env->setAppRoot( app_root );
    m_app = new InterSpecApp( *m_env );
    m_lock.reset( new Wt::WApplication::UpdateLock(m_app) );
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );

    // Both tools need a foreground before they will compute anything.
    const string spec = SpecUtils::append_path( InterSpec::staticDataDirectory(),
                    "reference_spectra/Common_Field_Nuclides/Detective X/Br82_Unshielded.txt" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(spec), "test spectrum not at '" << spec << "'" );

    auto meas = std::make_shared<SpecMeas>();
    BOOST_REQUIRE( meas->load_file( spec, SpecUtils::ParserType::Auto, spec ) );
    BOOST_REQUIRE( meas->num_measurements() > 0 );
    const std::shared_ptr<const SpecUtils::Measurement> m = meas->measurement_at_index( 0 );
    BOOST_REQUIRE( !!m );

    // The Detection Confidence Tool builds one row per gamma line, and it needs a detector response
    //  with resolution information to give each row a region of interest.  Without one there are no
    //  rows at all, and the per-row half of the state round trip would go untested.
    auto drf = std::make_shared<DetectorPeakResponse>();
    drf->setIntrinsicEfficiencyFormula( "exp(-343.63 + 269.10*log(x) - 83.80*log(x)^2"
                                        " + 12.44*log(x)^3 - 0.708*log(x)^4)",
                                        2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                        0.0f, 0.0f,
                                        DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    drf->setFwhmCoefficients( std::vector<float>{ 1.5f, 0.035f },
                             DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
    BOOST_REQUIRE( drf->isValid() && drf->hasResolutionInfo() );
    meas->setDetector( drf );

    m_interspec->setSpectrum( meas, {m->sample_number()}, SpecUtils::SpectrumType::Foreground, 0 );
  }

  ~ToolSessionFixture()
  {
    // `m_app` is owned by the environment - deleting it here as well crashes during teardown,
    //  after the test body has already passed, which reads confusingly as a product failure.
    m_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
    m_app = nullptr;
  }
};//class ToolSessionFixture
}//namespace


BOOST_AUTO_TEST_CASE( SimpleMdaStateRoundTrips )
{
  ToolSessionFixture fixture;

  DetectionLimitSimpleWindow * const window = fixture.m_interspec->showSimpleMdaWindow();
  BOOST_REQUIRE( window && window->tool() );
  DetectionLimitSimple * const tool = window->tool();

  // A state exercising every control the tool serializes, including the ones this increment
  //  touched: the measurement model, the continuum type, and the planned measurement time.
  const string uri = "DECON?VER=2.1&NUC=Cs137&ENERGY=661.657&DIST=100 cm&LROI=640&UROI=680"
                     "&CL=95&NSIDE=4&MODEL=BACKREF&CONTNORM=FLOAT&CONTTYPE=QUAD&SCALE=1800s"
                     "&ADV=1&ALPHA=0.02&BETA=0.1&DISTUNCERT=2 cm&EFFUNCERT=5";

  BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( uri ) );

  const string once = tool->encodeStateToUrl();
  BOOST_REQUIRE( !once.empty() );

  // Decoding what we just encoded, then encoding again, must give the same string.  Anything the
  //  encoder writes but the decoder drops (or vice versa) shows up right here.
  BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( once ) );
  const string twice = tool->encodeStateToUrl();

  BOOST_CHECK_MESSAGE( once == twice,
                      "Simple MDA state did not survive a round trip:\n  first:  " << once
                      << "\n  second: " << twice );

  // And the settings this increment cares about have to actually be in there - a round trip is
  //  stable if a field is dropped by *both* halves, so stability alone is not enough.
  BOOST_CHECK_MESSAGE( once.find("MODEL=BACKREF") != string::npos,
                      "measurement model missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("CONTTYPE=QUAD") != string::npos,
                      "continuum type missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("SCALE=") != string::npos,
                      "planned measurement time missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( SpecUtils::istarts_with(once, "DECON"),
                      "calculation method missing from encoded state: " << once );

  // Changing a setting must change the URI, or `render()` sees `sameAsPrev` and records no undo
  //  step.  Flipping the measurement model is the one this increment added behaviour to.
  const string backref_uri = tool->encodeStateToUrl();
  string current_uri = backref_uri;
  SpecUtils::ireplace_all( current_uri, "MODEL=BACKREF", "MODEL=CUR" );
  BOOST_REQUIRE( current_uri != backref_uri );

  BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( current_uri ) );
  const string after_current = tool->encodeStateToUrl();
  BOOST_CHECK_MESSAGE( after_current != backref_uri,
                      "switching the measurement model left the encoded state unchanged, so the"
                      " change cannot be undone" );
  BOOST_CHECK_MESSAGE( after_current.find("MODEL=CUR") != string::npos,
                      "measurement model did not switch: " << after_current );

  // ---- the "Advanced" section ------------------------------------------------------------------
  BOOST_CHECK_MESSAGE( once.find("ADV=1") != string::npos,
                      "advanced checkbox missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("ALPHA=") != string::npos,
                      "alpha missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("BETA=") != string::npos,
                      "beta missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("DISTUNCERT=") != string::npos,
                      "distance uncertainty missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("EFFUNCERT=") != string::npos,
                      "efficiency uncertainty missing from encoded state: " << once );

  // Turning "Advanced" off has to change the URI, or the change cannot be undone.
  {
    string adv_off = once;
    SpecUtils::ireplace_all( adv_off, "&ADV=1", "" );
    BOOST_REQUIRE( adv_off != once );

    BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( adv_off ) );
    const string after_adv_off = tool->encodeStateToUrl();
    BOOST_CHECK_MESSAGE( after_adv_off.find("ADV=1") == string::npos,
                        "advanced checkbox stayed on after being removed from the URI: "
                        << after_adv_off );
    BOOST_CHECK_MESSAGE( after_adv_off != once,
                        "turning the advanced section off left the encoded state unchanged, so the"
                        " change cannot be undone" );
  }

  // Alpha and beta track the confidence level until edited, and the ALPHA/BETA tokens are how that
  //  latch is encoded.  A state without ALPHA must come back still tracking - otherwise loading a
  //  state, then changing the confidence level, silently stops moving the error rates.  A plain
  //  round trip cannot see this: both halves agree on a value either way.
  {
    const string tracking = "CURRIE?VER=2.1&NUC=Cs137&ENERGY=661.657&DIST=100 cm&LROI=640&UROI=680"
                            "&CL=95&NSIDE=4&MODEL=CUR&ADV=1";
    BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( tracking ) );

    const string no_alpha = tool->encodeStateToUrl();
    BOOST_CHECK_MESSAGE( no_alpha.find("ADV=1") != string::npos, no_alpha );
    BOOST_CHECK_MESSAGE( no_alpha.find("ALPHA=") == string::npos,
                        "alpha latched by a state that never set it, so it has stopped following"
                        " the confidence level: " << no_alpha );
    BOOST_CHECK_MESSAGE( no_alpha.find("BETA=") == string::npos,
                        "beta latched by a state that never set it: " << no_alpha );
    BOOST_CHECK_MESSAGE( no_alpha.find("DISTUNCERT=") == string::npos,
                        "distance uncertainty survived a state that did not carry it: " << no_alpha );
    BOOST_CHECK_MESSAGE( no_alpha.find("EFFUNCERT=") == string::npos,
                        "efficiency uncertainty survived a state that did not carry it: " << no_alpha );

    // ...and it must still round trip in that un-latched condition.
    BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( no_alpha ) );
    BOOST_CHECK_MESSAGE( tool->encodeStateToUrl() == no_alpha,
                        "un-latched advanced state did not survive a round trip" );
  }

  // A VER=2 URI predates the advanced section entirely, and must decode to Advanced-off rather than
  //  inheriting whatever the tool happened to be showing.
  {
    const string legacy = "CURRIE?VER=2&NUC=Cs137&ENERGY=661.657&DIST=100 cm&LROI=640&UROI=680"
                          "&CL=95&NSIDE=4&MODEL=CUR";
    BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( legacy ) );

    const string from_legacy = tool->encodeStateToUrl();
    BOOST_CHECK_MESSAGE( from_legacy.find("ADV=") == string::npos,
                        "a pre-advanced URI left the advanced section on: " << from_legacy );
    BOOST_CHECK_MESSAGE( from_legacy.find("DISTUNCERT=") == string::npos,
                        "a pre-advanced URI left a stale distance uncertainty: " << from_legacy );
  }
}// BOOST_AUTO_TEST_CASE( SimpleMdaStateRoundTrips )


BOOST_AUTO_TEST_CASE( DetectionLimitToolStateRoundTrips )
{
  ToolSessionFixture fixture;

  DetectionLimitWindow * const window = fixture.m_interspec->createDetectionLimitTool();
  BOOST_REQUIRE( window && window->tool() );
  DetectionLimitTool * const tool = window->tool();

  // Includes a per-row entry that deviates from the row defaults on every axis the encoder tests
  //  (used-for-likelihood, side channels, continuum type), because rows at their defaults are
  //  deliberately *not* emitted - that is a URL-size optimization, not a gap.
  const string uri = "VER=1&NUC=Cs137&AGE=20y&LIM=ACT&DIST=100 cm&CL=95&AIRATTN=1"
                     "&MODEL=BACKREF&LIMTYPE=UPPER&SCALE=1800s"
                     "&ROW0=E:661.657,U:1,L:650,H:670,N:8,CN:FLOAT,CT:QUAD";

  BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( uri ) );

  const string once = tool->encodeStateToUrl();
  BOOST_REQUIRE( !once.empty() );

  BOOST_REQUIRE_NO_THROW( tool->handleAppUrl( once ) );
  const string twice = tool->encodeStateToUrl();

  BOOST_CHECK_MESSAGE( once == twice,
                      "Detection Confidence Tool state did not survive a round trip:\n  first:  "
                      << once << "\n  second: " << twice );

  BOOST_CHECK_MESSAGE( once.find("MODEL=") != string::npos,
                      "measurement model missing from encoded state: " << once );
  BOOST_CHECK_MESSAGE( once.find("LIMTYPE=") != string::npos,
                      "limit type missing from encoded state: " << once );

  // Per-row settings have no equivalent in the Simple MDA tool and are where a dropped field would
  //  be least obvious.  A row is only emitted when it deviates from its defaults, so the input
  //  above deviates deliberately; if the session produced no rows at all (no detector response,
  //  so no gamma lines to make rows from) there is nothing to assert and saying so beats a
  //  vacuous pass.
  //  The fixture installs a detector response with resolution info precisely so rows exist, so
  //  "no rows" is the regression worth failing on rather than skipping past.
  BOOST_REQUIRE_MESSAGE( once.find("ROW") != string::npos,
                        "no per-row state encoded, so the row round trip is untested: " << once );
  {
    BOOST_CHECK_MESSAGE( once.find("U:1") != string::npos,
                        "row's use-for-likelihood flag lost: " << once );
    BOOST_CHECK_MESSAGE( once.find("N:8") != string::npos,
                        "row's side-channel count lost: " << once );
    BOOST_CHECK_MESSAGE( once.find("CT:QUAD") != string::npos,
                        "row's continuum type lost: " << once );
  }
}// BOOST_AUTO_TEST_CASE( DetectionLimitToolStateRoundTrips )
