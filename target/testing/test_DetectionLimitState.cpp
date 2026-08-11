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
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE TestDetectionLimitState
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/DetectionLimitCalc.h"

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
