/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2.1 of the License, or (at your option)
 any later version.
 */

#include "InterSpec_config.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE RelActAutoTruthManifest_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;
using namespace boost::unit_test;

namespace
{
namespace fs = std::filesystem;
using json = nlohmann::json;

const array<const char *,5> sm_pu_isotopes{{ "Pu238", "Pu239", "Pu240", "Pu241", "Pu242" }};

string test_file_dir()
{
  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;

  string answer;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      answer = arg.substr( 14 );
  }

  SpecUtils::ireplace_all( answer, "%20", " " );
  if( answer.empty() )
  {
    for( const char *candidate : { "test_data", "../test_data", "../../test_data",
                                   "../../../target/testing/test_data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(candidate,
                       "rel_act_auto/pu_610_775/manifest.json") ) )
      {
        answer = candidate;
        break;
      }
    }
  }

  return answer;
}

json read_json( const fs::path &path )
{
  ifstream input( path );
  if( !input )
    throw runtime_error( "Could not open " + path.string() );
  return json::parse( input );
}

string read_text( const fs::path &path )
{
  ifstream input( path, ios::binary );
  if( !input )
    throw runtime_error( "Could not open " + path.string() );
  ostringstream contents;
  contents << input.rdbuf();
  return contents.str();
}

bool is_lower_hex_sha256( const string &value )
{
  if( value.size() != 64 )
    return false;
  for( const unsigned char c : value )
  {
    if( !std::isdigit(c) && !(c >= 'a' && c <= 'f') )
      return false;
  }
  return true;
}

std::uint32_t rotate_right( const std::uint32_t value, const unsigned int count )
{
  return (value >> count) | (value << (32u-count));
}

string sha256( const string &input )
{
  static const array<std::uint32_t,64> round_constants{{
    0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
    0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
    0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
    0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
    0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
    0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
    0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
    0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
  }};
  array<std::uint32_t,8> state{{
    0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,
    0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u
  }};

  vector<std::uint8_t> message;
  message.reserve(input.size()+72u);
  for( const unsigned char byte : input )
    message.push_back(byte);
  const std::uint64_t bit_count = static_cast<std::uint64_t>(message.size()) * 8u;
  message.push_back(0x80u);
  while( (message.size() % 64u) != 56u )
    message.push_back(0u);
  for( int shift = 56; shift >= 0; shift -= 8 )
    message.push_back(static_cast<std::uint8_t>(bit_count >> shift));

  for( size_t offset = 0; offset < message.size(); offset += 64u )
  {
    array<std::uint32_t,64> words{};
    for( size_t i = 0; i < 16u; ++i )
    {
      const size_t j = offset + 4u*i;
      words[i] = (static_cast<std::uint32_t>(message[j]) << 24u)
               | (static_cast<std::uint32_t>(message[j+1u]) << 16u)
               | (static_cast<std::uint32_t>(message[j+2u]) << 8u)
               | static_cast<std::uint32_t>(message[j+3u]);
    }
    for( size_t i = 16u; i < words.size(); ++i )
    {
      const std::uint32_t s0 = rotate_right(words[i-15u],7u)
                             ^ rotate_right(words[i-15u],18u) ^ (words[i-15u] >> 3u);
      const std::uint32_t s1 = rotate_right(words[i-2u],17u)
                             ^ rotate_right(words[i-2u],19u) ^ (words[i-2u] >> 10u);
      words[i] = words[i-16u] + s0 + words[i-7u] + s1;
    }

    std::uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
    std::uint32_t e = state[4], f = state[5], g = state[6], h = state[7];
    for( size_t i = 0; i < words.size(); ++i )
    {
      const std::uint32_t sum1 = rotate_right(e,6u) ^ rotate_right(e,11u)
                               ^ rotate_right(e,25u);
      const std::uint32_t choose = (e & f) ^ ((~e) & g);
      const std::uint32_t temp1 = h + sum1 + choose + round_constants[i] + words[i];
      const std::uint32_t sum0 = rotate_right(a,2u) ^ rotate_right(a,13u)
                               ^ rotate_right(a,22u);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temp2 = sum0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    state[0] += a; state[1] += b; state[2] += c; state[3] += d;
    state[4] += e; state[5] += f; state[6] += g; state[7] += h;
  }

  ostringstream digest;
  digest << hex << nouppercase << setfill('0');
  for( const std::uint32_t word : state )
    digest << setw(8) << word;
  return digest.str();
}

size_t maximum_gamma_channels( const SpecUtils::SpecFile &file )
{
  size_t answer = 0;
  for( const shared_ptr<const SpecUtils::Measurement> &measurement : file.measurements() )
    if( measurement )
      answer = std::max( answer, measurement->num_gamma_channels() );
  return answer;
}

shared_ptr<const SpecUtils::Measurement> largest_gamma_measurement(
                                                const SpecUtils::SpecFile &file )
{
  shared_ptr<const SpecUtils::Measurement> answer;
  for( const shared_ptr<const SpecUtils::Measurement> &measurement : file.measurements() )
    if( measurement && (!answer
        || (measurement->num_gamma_channels() > answer->num_gamma_channels())) )
      answer = measurement;
  return answer;
}
}//namespace


BOOST_AUTO_TEST_CASE( Sha256ImplementationMatchesStandardVectors )
{
  BOOST_CHECK_EQUAL( sha256(""),
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" );
  BOOST_CHECK_EQUAL( sha256("abc"),
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad" );
}


BOOST_AUTO_TEST_CASE( ManifestFilesAndProvenanceAreSelfConsistent )
{
  const string data_dir = test_file_dir();
  BOOST_REQUIRE_MESSAGE( !data_dir.empty(), "Could not locate target/testing/test_data" );
  const fs::path fixture_dir = fs::path(data_dir) / "rel_act_auto" / "pu_610_775";

  const json manifest = read_json( fixture_dir / "manifest.json" );
  BOOST_REQUIRE_EQUAL( manifest.at("schema_version").get<int>(), 1 );
  BOOST_REQUIRE_EQUAL( manifest.at("fixture_set").get<string>(), "rel_act_auto_pu_610_775" );

  const json &sources = manifest.at("sources");
  BOOST_REQUIRE_EQUAL( sources.size(), 2 );
  BOOST_CHECK_EQUAL( sources.at("jrc_upu_gamma").at("license").get<string>(), "CC-BY-4.0" );
  BOOST_CHECK_EQUAL( sources.at("idb_v2024_01").at("license").get<string>(), "CC-BY-4.0" );
  BOOST_CHECK( sources.at("idb_v2024_01").at("dataset_doi").get<string>().find("10.61092/") != string::npos );

  const fs::path copyright_path = fixture_dir
      / sources.at("jrc_upu_gamma").at("copyright_notice_file").get<string>();
  BOOST_REQUIRE( fs::is_regular_file(copyright_path) );
  BOOST_CHECK_EQUAL( fs::file_size(copyright_path), 540 );
  const string copyright_text = read_text( copyright_path );
  BOOST_CHECK( copyright_text.find("(c) European Union, 1995-2026") != string::npos );
  BOOST_CHECK( copyright_text.find("Creative Commons Attribution 4.0 International") != string::npos );
  BOOST_CHECK( copyright_text.find("changes are indicated") != string::npos );
  const string copyright_sha256
      = sources.at("jrc_upu_gamma").at("copyright_notice_sha256").get<string>();
  BOOST_REQUIRE( is_lower_hex_sha256(copyright_sha256) );
  BOOST_CHECK_EQUAL( sha256(copyright_text),copyright_sha256 );

  const json &artifacts = manifest.at("artifacts");
  BOOST_REQUIRE_EQUAL( artifacts.size(), 6 );
  for( const auto &[filename, artifact] : artifacts.items() )
  {
    BOOST_CHECK( filename.find('/') == string::npos );
    BOOST_CHECK( filename.find('\\') == string::npos );
    BOOST_CHECK( filename != "." && filename != ".." );
    BOOST_REQUIRE( sources.contains( artifact.at("source").get<string>() ) );
    BOOST_CHECK( !artifact.at("original_relative_path").get<string>().empty() );
    BOOST_CHECK( !artifact.at("changes").get<string>().empty() );
    BOOST_CHECK( is_lower_hex_sha256( artifact.at("sha256").get<string>() ) );

    const fs::path path = fixture_dir / filename;
    BOOST_REQUIRE_MESSAGE( fs::is_regular_file(path), "Missing fixture " << path );
    BOOST_CHECK_EQUAL( fs::file_size(path), artifact.at("size_bytes").get<uintmax_t>() );
    const string actual_sha256 = sha256(read_text(path));
    BOOST_CHECK_MESSAGE( actual_sha256 == artifact.at("sha256").get<string>(),
                         filename << " SHA-256 mismatch: got " << actual_sha256 );

    SpecUtils::SpecFile spectrum;
    BOOST_REQUIRE_MESSAGE( spectrum.load_file(path.string(), SpecUtils::ParserType::Auto),
                           "Could not parse fixture " << path );
    BOOST_CHECK_MESSAGE( maximum_gamma_channels(spectrum) > 0,
                         "No gamma channels found in " << path );
  }
}


BOOST_AUTO_TEST_CASE( CertificateTruthIsDerivedAndNormalized )
{
  const string data_dir = test_file_dir();
  BOOST_REQUIRE_MESSAGE( !data_dir.empty(), "Could not locate target/testing/test_data" );
  const fs::path fixture_dir = fs::path(data_dir) / "rel_act_auto" / "pu_610_775";
  const json manifest = read_json( fixture_dir / "manifest.json" );
  const json &convention = manifest.at("truth_convention");

  BOOST_CHECK_EQUAL( convention.at("certificate_uncertainty_sigma").get<int>(), 2 );
  BOOST_CHECK( !convention.at("jrc_generated_decay_corrected_field_is_truth").get<bool>() );
  const json &effective_age_convention = convention.at("effective_age_diagnostic");
  BOOST_CHECK( !effective_age_convention.at("certified_separation_truth").get<bool>() );
  BOOST_CHECK( effective_age_convention.at("derivation").get<string>().find("acquisition_truth")
               != string::npos );
  BOOST_CHECK_EQUAL( effective_age_convention.at("time_unit").get<string>(),
                     "365.25-day years" );

  const json &half_life = convention.at("half_life_seconds");
  const json &atomic_mass = convention.at("atomic_mass_u");
  const double branch = convention.at("pu241_to_am241_branching_ratio").get<double>();
  const double seconds_per_year = 365.25 * 24.0 * 60.0 * 60.0;
  const double ln2 = std::log(2.0);

  const json &cases = manifest.at("cases");
  BOOST_REQUIRE_EQUAL( cases.size(), 4 );
  size_t jrc_cases = 0;
  size_t idb_cases = 0;
  set<string> backgrounds;

  for( const json &test_case : cases )
  {
    const bool is_jrc = (test_case.at("source").get<string>() == "jrc_upu_gamma");
    jrc_cases += is_jrc ? 1 : 0;
    idb_cases += is_jrc ? 0 : 1;

    const json &certificate = test_case.at("certificate");
    const json &cert_pu = certificate.at("pu_mass_fractions_wt_percent");
    const json &cert_unc = certificate.at("pu_mass_fraction_uncertainty_wt_percent");
    const json &truth = test_case.at("acquisition_truth");
    const json &truth_pu = truth.at("pu_mass_fractions");

    BOOST_CHECK_EQUAL( certificate.at("uncertainty_sigma").get<int>(), 2 );
    BOOST_CHECK( !cert_pu.contains("Am241") );
    BOOST_CHECK( !truth_pu.contains("Am241") );

    double cert_sum = 0.0;
    double truth_sum = 0.0;
    for( const char *nuclide : sm_pu_isotopes )
    {
      BOOST_REQUIRE( cert_pu.contains(nuclide) );
      BOOST_REQUIRE( cert_unc.contains(nuclide) );
      BOOST_REQUIRE( truth_pu.contains(nuclide) );
      cert_sum += cert_pu.at(nuclide).get<double>();
      truth_sum += truth_pu.at(nuclide).get<double>();
      BOOST_CHECK( truth_pu.at(nuclide).get<double>() >= 0.0 );
      BOOST_CHECK( truth_pu.at(nuclide).get<double>() <= 1.0 );
    }
    BOOST_CHECK_SMALL( cert_sum - 100.0, 1.0e-3 );
    BOOST_CHECK_SMALL( truth_sum - 1.0, 5.0e-11 );

    const auto ref = SpecUtils::time_from_string(
                         certificate.at("reference_date_utc").get<string>() );
    const auto acq = SpecUtils::time_from_string(
                         test_case.at("acquisition_start_utc").get<string>() );
    BOOST_REQUIRE( !SpecUtils::is_special(ref) );
    BOOST_REQUIRE( !SpecUtils::is_special(acq) );
    const double elapsed_seconds = chrono::duration<double>( acq - ref ).count();
    BOOST_CHECK_CLOSE_FRACTION( elapsed_seconds / seconds_per_year,
                               truth.at("elapsed_years_365_25d").get<double>(), 5.0e-12 );

    double remaining_pu_mass = 0.0;
    array<double,5> remaining{};
    for( size_t i = 0; i < sm_pu_isotopes.size(); ++i )
    {
      const string nuclide = sm_pu_isotopes[i];
      const double initial_mass = cert_pu.at(nuclide).get<double>();
      remaining[i] = initial_mass * std::exp( -ln2 * elapsed_seconds
                                              / half_life.at(nuclide).get<double>() );
      remaining_pu_mass += remaining[i];
    }
    for( size_t i = 0; i < sm_pu_isotopes.size(); ++i )
    {
      const string nuclide = sm_pu_isotopes[i];
      BOOST_CHECK_CLOSE_FRACTION( remaining[i] / remaining_pu_mass,
                                 truth_pu.at(nuclide).get<double>(), 5.0e-9 );
    }

    const double lambda_pu241 = ln2 / half_life.at("Pu241").get<double>();
    const double lambda_am241 = ln2 / half_life.at("Am241").get<double>();
    const double initial_pu241_atoms = cert_pu.at("Pu241").get<double>()
                                      / atomic_mass.at("Pu241").get<double>();
    const double initial_am241_atoms = certificate.at("am241_to_total_pu_mass_wt_percent").get<double>()
                                      / atomic_mass.at("Am241").get<double>();
    const double final_am241_atoms = initial_am241_atoms * std::exp(-lambda_am241*elapsed_seconds)
        + branch * initial_pu241_atoms * lambda_pu241 / (lambda_am241 - lambda_pu241)
          * (std::exp(-lambda_pu241*elapsed_seconds) - std::exp(-lambda_am241*elapsed_seconds));
    const double am241_to_pu = final_am241_atoms * atomic_mass.at("Am241").get<double>()
                               / remaining_pu_mass;
    BOOST_CHECK_CLOSE_FRACTION( am241_to_pu,
                               truth.at("am241_to_total_pu_mass_ratio").get<double>(), 5.0e-10 );

    const string spectrum_name = test_case.at("spectrum_file").get<string>();
    BOOST_REQUIRE( manifest.at("artifacts").contains(spectrum_name) );
    SpecUtils::SpecFile spectrum;
    BOOST_REQUIRE( spectrum.load_file((fixture_dir / spectrum_name).string(),
                                      SpecUtils::ParserType::Auto) );
    const shared_ptr<const SpecUtils::Measurement> parsed = largest_gamma_measurement(spectrum);
    BOOST_REQUIRE( parsed );
    BOOST_CHECK_EQUAL( parsed->num_gamma_channels(),
                       test_case.at("detector").at("channels").get<size_t>() );
    const json &embedded_header = test_case.at("embedded_spe_header");
    BOOST_CHECK_LE( std::fabs(parsed->live_time()
                       - test_case.at("live_time_seconds").get<double>()),
                    embedded_header.at("live_time_absolute_tolerance_seconds").get<double>() );
    BOOST_CHECK_LE( std::fabs(parsed->real_time()
                       - test_case.at("real_time_seconds").get<double>()),
                    embedded_header.at("real_time_absolute_tolerance_seconds").get<double>() );
    const auto embedded_start = SpecUtils::time_from_string(
                       embedded_header.at("start_time_no_zone").get<string>() );
    BOOST_REQUIRE( !SpecUtils::is_special(embedded_start) );
    const double embedded_start_delta = std::fabs(
      chrono::duration<double>(parsed->start_time()-embedded_start).count() );
    BOOST_CHECK_LE( embedded_start_delta,
      0.5*embedded_header.at("start_time_resolution_seconds").get<double>() );
    double parsed_counts = 0.0;
    for( const shared_ptr<const SpecUtils::Measurement> &measurement : spectrum.measurements() )
      if( measurement && measurement->gamma_counts() )
        for( const float count : *measurement->gamma_counts() )
          parsed_counts += count;
    BOOST_CHECK_CLOSE_FRACTION( parsed_counts,
                               test_case.at("total_gamma_counts").get<double>(), 1.0e-12 );

    if( test_case.at("background_file").is_null() )
    {
      BOOST_CHECK( !is_jrc );
      BOOST_CHECK( certificate.at("separation_date_utc").is_null() );
      BOOST_CHECK( certificate.at("separation_date_note").get<string>().find("No distinct")
                   != string::npos );
      BOOST_CHECK_EQUAL( test_case.at("fit_acceptance").at("mode").get<string>(),
                         "provenance_parsing_and_robustness_only" );

      // IDB independently publishes an acquisition-date composition in the SPE metadata.  It is
      // not the propagation oracle, but it must agree with the certificate/date derivation within
      // the documented 0.02 absolute wt%-point tolerance.
      const json &published = test_case.at("source_metadata_acquisition_composition");
      BOOST_CHECK_EQUAL( published.at("units").get<string>(), "fraction" );
      const double wt_percent_tolerance
          = published.at("comparison_tolerance_absolute_wt_percent").get<double>();
      for( const char *nuclide : sm_pu_isotopes )
      {
        const double delta_wt_percent = 100.0*std::fabs(
            published.at("pu_mass_fractions").at(nuclide).get<double>()
            - truth_pu.at(nuclide).get<double>() );
        BOOST_CHECK_LE( delta_wt_percent,wt_percent_tolerance );
      }
      const double am_delta_wt_percent = 100.0*std::fabs(
          published.at("am241_to_total_pu_mass_ratio").get<double>()
          - truth.at("am241_to_total_pu_mass_ratio").get<double>() );
      BOOST_CHECK_LE( am_delta_wt_percent,wt_percent_tolerance );
    }else
    {
      BOOST_CHECK( is_jrc );
      const string background = test_case.at("background_file").get<string>();
      backgrounds.insert( background );
      BOOST_REQUIRE( manifest.at("artifacts").contains(background) );
      const json &acceptance = test_case.at("fit_acceptance");
      const json &known_age = acceptance.at("known_age_diagnostic");
      BOOST_CHECK_CLOSE_FRACTION(
        known_age.at("total_variation_max").get<double>(),
        0.05, 1.0e-12 );
      const double decay_constant_delta = lambda_pu241-lambda_am241;
      const double acquisition_am_to_pu241_mass
          = truth.at("am241_to_total_pu_mass_ratio").get<double>()
            / truth_pu.at("Pu241").get<double>();
      const double effective_age_seconds = std::log1p(
          acquisition_am_to_pu241_mass
            * atomic_mass.at("Pu241").get<double>() / atomic_mass.at("Am241").get<double>()
            * decay_constant_delta / (branch*lambda_pu241) ) / decay_constant_delta;
      BOOST_CHECK_CLOSE_FRACTION(
        effective_age_seconds/seconds_per_year,
        known_age.at("effective_age_years_365_25d").get<double>(),5.0e-12 );
      BOOST_CHECK_EQUAL( acceptance.at("one_sigma_fraction_max").size(), 5 );
      BOOST_CHECK_CLOSE_FRACTION(
        acceptance.at("relative_one_sigma_max_for_truth_fraction_at_least_0_01").get<double>(),
        1.0, 1.0e-12 );
      // The JRC certificate-date value must not masquerade as acquisition truth.
      BOOST_CHECK( std::abs( 0.01*cert_pu.at("Pu241").get<double>()
                             - truth_pu.at("Pu241").get<double>() ) > 1.0e-4 );
    }
  }

  BOOST_CHECK_EQUAL( jrc_cases, 3 );
  BOOST_CHECK_EQUAL( idb_cases, 1 );
  BOOST_CHECK_EQUAL( backgrounds.size(), 2 );
}


BOOST_AUTO_TEST_CASE( ExistingManualFixtureRecordsItsTransformations )
{
  const string data_dir = test_file_dir();
  BOOST_REQUIRE_MESSAGE( !data_dir.empty(), "Could not locate target/testing/test_data" );
  const string notice = read_text( fs::path(data_dir) / "manual_rel_eff" / "source.txt" );
  BOOST_CHECK( notice.find("IDB-v2024-01") != string::npos );
  BOOST_CHECK( notice.find("record 184") != string::npos );
  BOOST_CHECK( notice.find("CC BY 4.0") != string::npos );
  BOOST_CHECK( notice.find("converted to N42-2012") != string::npos );
  BOOST_CHECK( notice.find("0.0692044422") != string::npos );
  BOOST_CHECK( notice.find("peak fits/nuclide assignments") != string::npos );
  BOOST_CHECK( notice.find("f831779ff25eeae31af532613b4c1d76da5b55c6fcdfaa6e32f108297efe8c07") != string::npos );
  BOOST_CHECK( notice.find("014bc03cebd62c0583f76fa921695a146ca1fd7d666b9bf3720021d230edad8c") != string::npos );
}
