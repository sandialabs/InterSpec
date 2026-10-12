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

#include <regex>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#define BOOST_TEST_MODULE LlmConversationLog_suite
#include <boost/test/included/unit_test.hpp>

#include <rapidxml/rapidxml.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace std;


namespace
{
  string write_temp_file( const string &contents )
  {
    const string path = SpecUtils::temp_file_name( "llm_log_test", SpecUtils::temp_dir() );
    ofstream out( path.c_str(), ios::binary | ios::out );
    out << contents;
    return path;
  }

  string read_file( const string &path )
  {
    ifstream in( path.c_str(), ios::binary | ios::in );
    return string( (istreambuf_iterator<char>(in)), istreambuf_iterator<char>() );
  }

  const char * const sm_minimal_config =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<LlmConfig version=\"0\">\n"
    "  <LlmApi version=\"1\">\n"
    "    <Enabled>true</Enabled>\n"
    "    <ApiProvider active=\"true\">\n"
    "      <ApiEndpoint>https://api.example.com/v1/chat/completions</ApiEndpoint>\n"
    "      <BearerToken>x</BearerToken>\n"
    "      <Model active=\"true\"><Name>test-model</Name></Model>\n"
    "    </ApiProvider>\n"
    "  </LlmApi>\n"
    "  <McpServer version=\"1\"><Enabled>false</Enabled><BearerToken>x</BearerToken></McpServer>\n"
    "</LlmConfig>\n";
}//namespace


BOOST_AUTO_TEST_CASE( LogConversationsConfigRoundTrip )
{
  // A config without <LogConversations> must load with logging off.
  const string orig_path = write_temp_file( sm_minimal_config );
  LlmConfig config;
  std::tie( config.llmApi, config.mcpServer ) = LlmConfig::loadApiAndMcpConfigs( orig_path );
  SpecUtils::remove_file( orig_path );
  BOOST_CHECK( !config.llmApi.logConversations );

  // When off, the element is not written at all, so older configs round-trip unchanged.
  const string off_xml = LlmConfig::toXmlString( config );
  BOOST_CHECK( off_xml.find( "<LogConversations>" ) == string::npos );

  // When on, it is written and read back.
  config.llmApi.logConversations = true;
  const string on_xml = LlmConfig::toXmlString( config );
  BOOST_CHECK( on_xml.find( "<LogConversations>true</LogConversations>" ) != string::npos );

  const string on_path = write_temp_file( on_xml );
  const pair<LlmConfig::LlmApi, LlmConfig::McpServer> reloaded = LlmConfig::loadApiAndMcpConfigs( on_path );
  SpecUtils::remove_file( on_path );
  BOOST_CHECK( reloaded.first.logConversations );
}//LogConversationsConfigRoundTrip


BOOST_AUTO_TEST_CASE( LogFileNameFormat )
{
  const string name = LlmConversationHistory::logFileName( chrono::system_clock::now(), "aB3dE6" );
  BOOST_CHECK_MESSAGE( std::regex_match( name, std::regex( R"(\d{8}T\d{4}_aB3dE6_llm_log\.xml)" ) ),
                       "Unexpected log file name: '" << name << "'" );
}//LogFileNameFormat


BOOST_AUTO_TEST_CASE( WriteLogFileRoundTrip )
{
  LlmConversationHistory history;
  const shared_ptr<LlmInteraction> first = history.addUserMessageToMainConversation( "What is in my spectrum?" );
  history.addAssistantMessageWithThinking( "Looks like Cs-137 & <Co-60>.", "thinking...", "", first );
  LlmConversationHistory::addTokenUsage( first, 1200, 34, 1234, 1000, 5 );

  const shared_ptr<LlmInteraction> second = history.addUserMessageToMainConversation( "Thanks" );
  history.addAssistantMessageWithThinking( "You're welcome.", "", "", second );

  const string path = SpecUtils::temp_file_name( "llm_log_test", SpecUtils::temp_dir() ) + ".xml";
  history.writeLogFile( path, { {"model", "test-model"}, {"started", "2026-10-08T12:32:00-0700"} } );
  // Writing again must replace the file (not fail because it exists), and leave no temp file behind.
  history.writeLogFile( path, { {"model", "test-model"} } );
  BOOST_REQUIRE( SpecUtils::is_file( path ) );
  BOOST_CHECK( !SpecUtils::is_file( path + ".tmp" ) );

  string contents = read_file( path );
  SpecUtils::remove_file( path );

  rapidxml::xml_document<char> doc;
  BOOST_REQUIRE_NO_THROW( doc.parse<rapidxml::parse_trim_whitespace>( &contents[0] ) );

  const rapidxml::xml_node<char> * const root = XML_FIRST_NODE( (&doc), "LlmConversationLog" );
  BOOST_REQUIRE( root );
  BOOST_CHECK_EQUAL( SpecUtils::xml_value_str( XML_FIRST_ATTRIB( root, "version" ) ), "0" );
  BOOST_CHECK_EQUAL( SpecUtils::xml_value_str( XML_FIRST_ATTRIB( root, "model" ) ), "test-model" );

  vector<shared_ptr<LlmInteraction>> loaded;
  LlmConversationHistory::fromXml( root, loaded );
  BOOST_REQUIRE_EQUAL( loaded.size(), 2 );
  BOOST_CHECK_EQUAL( loaded[0]->conversationId, first->conversationId );
  BOOST_CHECK_EQUAL( loaded[1]->conversationId, second->conversationId );
  BOOST_REQUIRE_EQUAL( loaded[0]->responses.size(), 2 );

  const shared_ptr<const LlmInteractionFinalResponse> reply
                = dynamic_pointer_cast<const LlmInteractionFinalResponse>( loaded[0]->responses[1] );
  BOOST_REQUIRE( reply );
  BOOST_CHECK_EQUAL( reply->content(), "Looks like Cs-137 & <Co-60>." );

  // Token counts are now serialized too.
  BOOST_CHECK( loaded[0]->promptTokens == first->promptTokens );
  BOOST_CHECK( loaded[0]->completionTokens == first->completionTokens );
  BOOST_CHECK( loaded[0]->totalTokens == first->totalTokens );
  BOOST_CHECK( loaded[0]->cachedTokens == first->cachedTokens );
  BOOST_CHECK( loaded[0]->cacheCreationTokens == first->cacheCreationTokens );
  BOOST_CHECK( first->totalTokens.has_value() );
  BOOST_CHECK( !loaded[1]->totalTokens.has_value() );
}//WriteLogFileRoundTrip
