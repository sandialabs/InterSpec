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

#include <string>
#include <memory>
#include <sstream>

#define BOOST_TEST_MODULE PeakFitDetPrefs_suite
#include <boost/test/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"

using namespace std;


BOOST_AUTO_TEST_CASE( test_default_construction )
{
  PeakFitDetPrefs prefs;
  BOOST_CHECK( prefs.m_det_type == PeakFitUtils::CoarseResolutionType::Unknown );
  BOOST_CHECK( prefs.m_peak_skew_type == PeakDef::NoSkew );
  BOOST_CHECK( prefs.m_source == PeakFitDetPrefs::LoadingSource::Default );

  for( int i = 0; i < 4; ++i )
  {
    BOOST_CHECK( !prefs.m_lower_energy_skew[i].has_value() );
    BOOST_CHECK( !prefs.m_upper_energy_skew[i].has_value() );
  }
}//test_default_construction


BOOST_AUTO_TEST_CASE( test_equality_operators )
{
  PeakFitDetPrefs a, b;
  BOOST_CHECK( a == b );
  BOOST_CHECK( !(a != b) );

  a.m_det_type = PeakFitUtils::CoarseResolutionType::High;
  BOOST_CHECK( a != b );

  b.m_det_type = PeakFitUtils::CoarseResolutionType::High;
  BOOST_CHECK( a == b );

  a.m_lower_energy_skew[0] = 1.5;
  BOOST_CHECK( a != b );

  b.m_lower_energy_skew[0] = 1.5;
  BOOST_CHECK( a == b );

  a.m_source = PeakFitDetPrefs::LoadingSource::UserInputInGui;
  BOOST_CHECK( a != b );
}//test_equality_operators


BOOST_AUTO_TEST_CASE( test_coarse_res_str_roundtrip )
{
  using CRT = PeakFitUtils::CoarseResolutionType;

  const CRT types[] = { CRT::Low, CRT::LaBr, CRT::CZT, CRT::High, CRT::Unknown };
  for( const CRT t : types )
  {
    const char *s = PeakFitDetPrefs::to_str( t );
    BOOST_CHECK( s != nullptr );
    const CRT parsed = PeakFitDetPrefs::coarse_res_from_str( s );
    BOOST_CHECK( parsed == t );
  }

  // Test aliases
  BOOST_CHECK( PeakFitDetPrefs::coarse_res_from_str( "NaI" ) == CRT::Low );
  BOOST_CHECK( PeakFitDetPrefs::coarse_res_from_str( "LaBr3" ) == CRT::LaBr );
  BOOST_CHECK( PeakFitDetPrefs::coarse_res_from_str( "HPGe" ) == CRT::High );

  // Test invalid
  BOOST_CHECK_THROW( PeakFitDetPrefs::coarse_res_from_str( "InvalidType" ), std::exception );
}//test_coarse_res_str_roundtrip


BOOST_AUTO_TEST_CASE( test_loading_source_str_roundtrip )
{
  using LS = PeakFitDetPrefs::LoadingSource;

  const LS sources[] = { LS::Default, LS::UserInputInGui, LS::FromDetectorPeakResponse,
                          LS::DefaultForDetectorType, LS::FromSpectralData };
  for( const LS s : sources )
  {
    const char *str = PeakFitDetPrefs::to_str( s );
    BOOST_CHECK( str != nullptr );
    const LS parsed = PeakFitDetPrefs::loading_source_from_str( str );
    BOOST_CHECK( parsed == s );
  }

  BOOST_CHECK_THROW( PeakFitDetPrefs::loading_source_from_str( "BadSource" ), std::exception );
}//test_loading_source_str_roundtrip


BOOST_AUTO_TEST_CASE( test_xml_roundtrip_full )
{
  PeakFitDetPrefs orig;
  orig.m_det_type = PeakFitUtils::CoarseResolutionType::High;
  orig.m_peak_skew_type = PeakDef::DoubleSidedCrystalBall;
  orig.m_lower_energy_skew[0] = 1.5;
  orig.m_upper_energy_skew[0] = 2.5;
  orig.m_lower_energy_skew[1] = 3.0;
  // SkewPar1 for DSCB is not energy dependent, so upper not needed
  orig.m_lower_energy_skew[2] = 0.8;
  orig.m_upper_energy_skew[2] = 1.2;
  orig.m_lower_energy_skew[3] = 5.0;
  // SkewPar3 for DSCB is not energy dependent
  orig.m_source = PeakFitDetPrefs::LoadingSource::UserInputInGui;

  // Serialize to XML
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "Root" );
  doc.append_node( root );

  orig.toXml( root, &doc );

  // Print to string and re-parse
  string xml_str;
  rapidxml::print( back_inserter( xml_str ), doc, 0 );

  rapidxml::xml_document<char> doc2;
  // Need to make a non-const copy for rapidxml parsing
  vector<char> xml_copy( xml_str.begin(), xml_str.end() );
  xml_copy.push_back( '\0' );
  doc2.parse<rapidxml::parse_default>( xml_copy.data() );

  const rapidxml::xml_node<char> *root2 = doc2.first_node( "Root" );
  BOOST_REQUIRE( root2 );
  const rapidxml::xml_node<char> *pfp_node = root2->first_node( "PeakFitDetPrefs" );
  BOOST_REQUIRE( pfp_node );

  PeakFitDetPrefs parsed;
  parsed.fromXml( pfp_node );

  BOOST_CHECK( parsed == orig );
}//test_xml_roundtrip_full


BOOST_AUTO_TEST_CASE( test_xml_roundtrip_minimal )
{
  PeakFitDetPrefs orig;
  orig.m_det_type = PeakFitUtils::CoarseResolutionType::Low;
  orig.m_peak_skew_type = PeakDef::NoSkew;
  orig.m_source = PeakFitDetPrefs::LoadingSource::Default;
  // No skew params set

  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "Root" );
  doc.append_node( root );

  orig.toXml( root, &doc );

  string xml_str;
  rapidxml::print( back_inserter( xml_str ), doc, 0 );

  rapidxml::xml_document<char> doc2;
  vector<char> xml_copy( xml_str.begin(), xml_str.end() );
  xml_copy.push_back( '\0' );
  doc2.parse<rapidxml::parse_default>( xml_copy.data() );

  const rapidxml::xml_node<char> *root2 = doc2.first_node( "Root" );
  BOOST_REQUIRE( root2 );
  const rapidxml::xml_node<char> *pfp_node = root2->first_node( "PeakFitDetPrefs" );
  BOOST_REQUIRE( pfp_node );

  PeakFitDetPrefs parsed;
  parsed.fromXml( pfp_node );

  BOOST_CHECK( parsed == orig );
}//test_xml_roundtrip_minimal


BOOST_AUTO_TEST_CASE( test_url_roundtrip_full )
{
  PeakFitDetPrefs orig;
  orig.m_det_type = PeakFitUtils::CoarseResolutionType::LaBr;
  orig.m_peak_skew_type = PeakDef::Bortel;
  orig.m_lower_energy_skew[0] = 1.23;
  orig.m_upper_energy_skew[0] = 4.56;
  orig.m_source = PeakFitDetPrefs::LoadingSource::FromDetectorPeakResponse;

  const string url_parts = orig.toUrlQueryParts();
  BOOST_CHECK( !url_parts.empty() );

  // Parse back the URL parts
  // Simple parser: split by '&', then by '='
  map<string, string> parts;
  {
    istringstream iss( url_parts );
    string token;
    while( getline( iss, token, '&' ) )
    {
      const size_t eq = token.find( '=' );
      if( eq != string::npos )
        parts[token.substr( 0, eq )] = token.substr( eq + 1 );
    }
  }

  PeakFitDetPrefs parsed;
  parsed.fromUrlQueryParts( parts );

  // Source is not preserved in URL format (defaults to Default)
  BOOST_CHECK( parsed.m_det_type == orig.m_det_type );
  BOOST_CHECK( parsed.m_peak_skew_type == orig.m_peak_skew_type );

  // Check skew params
  BOOST_REQUIRE( parsed.m_lower_energy_skew[0].has_value() );
  BOOST_CHECK_CLOSE( parsed.m_lower_energy_skew[0].value(), 1.23, 1e-5 );
  BOOST_REQUIRE( parsed.m_upper_energy_skew[0].has_value() );
  BOOST_CHECK_CLOSE( parsed.m_upper_energy_skew[0].value(), 4.56, 1e-5 );
}//test_url_roundtrip_full


BOOST_AUTO_TEST_CASE( test_url_roundtrip_minimal )
{
  PeakFitDetPrefs orig;
  orig.m_det_type = PeakFitUtils::CoarseResolutionType::Unknown;
  orig.m_peak_skew_type = PeakDef::NoSkew;

  const string url_parts = orig.toUrlQueryParts();
  BOOST_CHECK( !url_parts.empty() );

  map<string, string> parts;
  {
    istringstream iss( url_parts );
    string token;
    while( getline( iss, token, '&' ) )
    {
      const size_t eq = token.find( '=' );
      if( eq != string::npos )
        parts[token.substr( 0, eq )] = token.substr( eq + 1 );
    }
  }

  PeakFitDetPrefs parsed;
  parsed.fromUrlQueryParts( parts );

  BOOST_CHECK( parsed.m_det_type == orig.m_det_type );
  BOOST_CHECK( parsed.m_peak_skew_type == orig.m_peak_skew_type );

  for( int i = 0; i < 4; ++i )
  {
    BOOST_CHECK( !parsed.m_lower_energy_skew[i].has_value() );
    BOOST_CHECK( !parsed.m_upper_energy_skew[i].has_value() );
  }
}//test_url_roundtrip_minimal


BOOST_AUTO_TEST_CASE( test_default_for_detector_type )
{
  using DT = SpecUtils::DetectorType;

  // HPGe detectors should return High
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::DetectiveEx ), "", "" );
    BOOST_REQUIRE( prefs );
    BOOST_CHECK( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::High );
    BOOST_CHECK( prefs->m_source == PeakFitDetPrefs::LoadingSource::DefaultForDetectorType );
  }

  // NaI detectors should return Low
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::IdentiFinder ), "", "" );
    BOOST_REQUIRE( prefs );
    BOOST_CHECK( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::Low );
  }

  // LaBr3 detectors
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::IdentiFinderLaBr3 ), "", "" );
    BOOST_REQUIRE( prefs );
    BOOST_CHECK( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::LaBr );
  }

  // CZT detectors
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::KromekGR1 ), "", "" );
    BOOST_REQUIRE( prefs );
    BOOST_CHECK( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::CZT );
  }

  // Unknown with no keywords should return nullptr
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::Unknown ), "", "" );
    BOOST_CHECK( !prefs );
  }

  // Unknown with HPGe keyword should return High
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::Unknown ), "Canberra", "HPGe coaxial" );
    BOOST_REQUIRE( prefs );
    BOOST_CHECK( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::High );
  }

  // Unknown with LaBr keyword
  {
    auto prefs = PeakFitDetPrefs::defaultForDetectorType( static_cast<int>( DT::Unknown ), "", "LaBr detector" );
    BOOST_REQUIRE( prefs );
    BOOST_CHECK( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::LaBr );
  }
}//test_default_for_detector_type


BOOST_AUTO_TEST_CASE( test_guess_from_spectral_data_stub )
{
  auto prefs = PeakFitDetPrefs::guessFromSpectralData( nullptr );
  BOOST_CHECK( !prefs );
}//test_guess_from_spectral_data_stub
