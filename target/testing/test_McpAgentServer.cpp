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

/* Tests for the pieces of InterSpec's MCP agent server that do not need a running server:
   env(...) references and the <AgentServer> settings in llm_config.xml, file-name and path
   confinement, and that an exported N42 carries the analysis (peaks, models, the assistant
   conversation with its images) through a write and re-load.
 */

#include "InterSpec_config.h"

#include <set>
#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

#define BOOST_TEST_MODULE McpAgentServer_suite
#include <boost/test/included/unit_test.hpp>

#include <rapidxml/rapidxml.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecServer.h"
#include "target/mcp_http/McpSessionManager.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace std;
namespace fs = std::filesystem;


namespace
{
  string write_temp_file( const string &contents, const string &ext = ".xml" )
  {
    const string path = SpecUtils::temp_file_name( "mcp_agent_test", SpecUtils::temp_dir() ) + ext;
    ofstream out( path.c_str(), ios::binary | ios::out );
    out << contents;
    return path;
  }

  /** Sets (or, with value nullptr, unsets) an environment variable. */
  void set_env( const char *name, const char *value )
  {
#ifdef _WIN32
    _putenv_s( name, value ? value : "" );
#else
    if( value )
      setenv( name, value, 1 );
    else
      unsetenv( name );
#endif
  }

  /** Has conversation histories save their users' images for the scope, as the MCP agent server does. */
  struct SaveUserImagesScope
  {
    const bool m_previous = LlmConversationHistory::saveUserImages();
    SaveUserImagesScope(){ LlmConversationHistory::setSaveUserImages( true ); }
    ~SaveUserImagesScope(){ LlmConversationHistory::setSaveUserImages( m_previous ); }
  };

  const char * const sm_env_config =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<LlmConfig version=\"0\">\n"
    "  <LlmApi version=\"1\">\n"
    "    <Enabled>true</Enabled>\n"
    "    <ApiProvider active=\"true\" apiFormat=\"env(INTERSPEC_TEST_API_FORMAT,anthropic)\">\n"
    "      <ApiEndpoint>env(INTERSPEC_TEST_ENDPOINT)</ApiEndpoint>\n"
    "      <BearerToken> env( INTERSPEC_TEST_TOKEN ) </BearerToken>\n"
    "      <Model active=\"true\" supportsImages=\"env(INTERSPEC_TEST_IMAGES,true)\">\n"
    "        <Name>env(INTERSPEC_TEST_MODEL,default-model)</Name>\n"
    "        <MaxTokens>env(INTERSPEC_TEST_MAX_TOKENS,1234)</MaxTokens>\n"
    "      </Model>\n"
    "    </ApiProvider>\n"
    "    <CompactionSystemPrompt>Say env(INTERSPEC_TEST_TOKEN) here</CompactionSystemPrompt>\n"
    "  </LlmApi>\n"
    "  <McpServer version=\"1\">\n"
    "    <Enabled>true</Enabled>\n"
    "    <BearerToken>env(INTERSPEC_TEST_MCP_TOKEN)</BearerToken>\n"
    "    <AllowAssistantPrompt>true</AllowAssistantPrompt>\n"
    "    <AgentServer>\n"
    "      <MaxSessions>env(INTERSPEC_TEST_MAX_SESSIONS,3)</MaxSessions>\n"
    "      <ExchangeDirectory>/tmp/exchange-test</ExchangeDirectory>\n"
    "      <KeepSessionFiles>true</KeepSessionFiles>\n"
    "      <MaxImageBytes>1000</MaxImageBytes>\n"
    "      <PublicBaseUrl>env(INTERSPEC_TEST_BASE_URL,)</PublicBaseUrl>\n"
    "    </AgentServer>\n"
    "  </McpServer>\n"
    "</LlmConfig>\n";
}//namespace


BOOST_AUTO_TEST_CASE( EnvReferenceResolution )
{
  set_env( "INTERSPEC_TEST_VALUE", "the value" );
  set_env( "INTERSPEC_TEST_UNSET", nullptr );
  set_env( "INTERSPEC_TEST_EMPTY", "" );

  bool was_ref = true;
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "plain text", was_ref ), "plain text" );
  BOOST_CHECK( !was_ref );

  // Only a value that is entirely a reference is resolved.
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "x env(INTERSPEC_TEST_VALUE)", was_ref ), "x env(INTERSPEC_TEST_VALUE)" );
  BOOST_CHECK( !was_ref );
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_VALUE) x", was_ref ), "env(INTERSPEC_TEST_VALUE) x" );
  BOOST_CHECK( !was_ref );

  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_VALUE)", was_ref ), "the value" );
  BOOST_CHECK( was_ref );
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "  env( INTERSPEC_TEST_VALUE , dflt )  ", was_ref ), "the value" );
  BOOST_CHECK( was_ref );

  // Defaults, including an empty one; a set-but-empty variable is used as is.
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_UNSET,some default)", was_ref ), "some default" );
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_UNSET,)", was_ref ), "" );
  BOOST_CHECK( was_ref );
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_EMPTY,dflt)", was_ref ), "" );

  // A default may itself contain commas.
  BOOST_CHECK_EQUAL( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_UNSET,a,b)", was_ref ), "a,b" );

  // Errors: unset with no default, a variable without the INTERSPEC_ prefix, or a malformed name.
  BOOST_CHECK_THROW( LlmConfig::resolveEnvReference( "env(INTERSPEC_TEST_UNSET)", was_ref ), std::exception );
  BOOST_CHECK_THROW( LlmConfig::resolveEnvReference( "env(HOME)", was_ref ), std::exception );
  BOOST_CHECK_THROW( LlmConfig::resolveEnvReference( "env(PATH,x)", was_ref ), std::exception );
  BOOST_CHECK_THROW( LlmConfig::resolveEnvReference( "env(INTERSPEC TEST)", was_ref ), std::exception );
  BOOST_CHECK_THROW( LlmConfig::resolveEnvReference( "env()", was_ref ), std::exception );
}//EnvReferenceResolution


BOOST_AUTO_TEST_CASE( EnvReferencesInConfigFile )
{
  set_env( "INTERSPEC_TEST_ENDPOINT", "https://llm.example.com/v1/messages" );
  set_env( "INTERSPEC_TEST_TOKEN", "secret-llm-token" );
  set_env( "INTERSPEC_TEST_MCP_TOKEN", "secret-mcp-token-0123456789" );
  set_env( "INTERSPEC_TEST_MODEL", nullptr );
  set_env( "INTERSPEC_TEST_API_FORMAT", nullptr );
  set_env( "INTERSPEC_TEST_IMAGES", nullptr );
  set_env( "INTERSPEC_TEST_MAX_TOKENS", nullptr );
  set_env( "INTERSPEC_TEST_MAX_SESSIONS", "5" );
  set_env( "INTERSPEC_TEST_BASE_URL", nullptr );

  const string path = write_temp_file( sm_env_config );

  bool uses_env = false;
  pair<LlmConfig::LlmApi, LlmConfig::McpServer> parsed;
  BOOST_REQUIRE_NO_THROW( parsed = LlmConfig::loadApiAndMcpConfigs( path, &uses_env ) );
  BOOST_CHECK( uses_env );

  const LlmConfig::LlmApi &api = parsed.first;
  BOOST_REQUIRE_EQUAL( api.providers.size(), 1u );
  const LlmConfig::LlmApi::ApiProvider &provider = api.providers[0];
  BOOST_CHECK_EQUAL( provider.apiEndpoint, "https://llm.example.com/v1/messages" );
  BOOST_CHECK_EQUAL( provider.bearerToken, "secret-llm-token" );
  BOOST_REQUIRE( provider.apiFormat.has_value() );
  BOOST_CHECK( provider.apiFormat.value() == LlmConfig::LlmApi::ApiFormat::Anthropic );
  BOOST_REQUIRE_EQUAL( provider.models.size(), 1u );
  BOOST_CHECK_EQUAL( provider.models[0].name, "default-model" );
  BOOST_CHECK( provider.models[0].supportsImages );
  BOOST_CHECK_EQUAL( provider.models[0].maxTokens, 1234 );

  // Text merely containing a reference is left alone.
  BOOST_CHECK_EQUAL( api.compactionSystemPrompt, "Say env(INTERSPEC_TEST_TOKEN) here" );

  const LlmConfig::McpServer &mcp = parsed.second;
  BOOST_CHECK( mcp.enabled );
  BOOST_CHECK_EQUAL( mcp.bearerToken, "secret-mcp-token-0123456789" );
  BOOST_REQUIRE( mcp.agentServer.has_value() );
  BOOST_CHECK_EQUAL( mcp.agentServer->maxSessions, 5u );
  BOOST_CHECK_EQUAL( mcp.agentServer->exchangeDirectory, "/tmp/exchange-test" );
  BOOST_CHECK( mcp.agentServer->keepSessionFiles );
  BOOST_CHECK_EQUAL( mcp.agentServer->maxImageBytes, 1000u );
  BOOST_CHECK( mcp.agentServer->publicBaseUrl.empty() );
  // Unspecified elements keep their defaults.
  BOOST_CHECK_EQUAL( mcp.agentServer->sessionIdleTimeoutSeconds, LlmConfig::McpServer::AgentServer{}.sessionIdleTimeoutSeconds );

  // The <AgentServer> element survives serialization.
  LlmConfig config;
  config.llmApi = api;
  config.mcpServer = mcp;
  const string xml = LlmConfig::toXmlString( config );
  BOOST_CHECK( xml.find( "<AgentServer>" ) != string::npos );
  const string resaved = write_temp_file( xml );
  bool resaved_uses_env = true;
  const pair<LlmConfig::LlmApi, LlmConfig::McpServer> reparsed = LlmConfig::loadApiAndMcpConfigs( resaved, &resaved_uses_env );
  BOOST_CHECK( !resaved_uses_env );
  BOOST_REQUIRE( reparsed.second.agentServer.has_value() );
  BOOST_CHECK_EQUAL( reparsed.second.agentServer->maxSessions, 5u );

  // A desktop config without <AgentServer> does not gain one.
  config.mcpServer.agentServer.reset();
  BOOST_CHECK( LlmConfig::toXmlString( config ).find( "AgentServer" ) == string::npos );

  // An unset variable without a default is an error.
  set_env( "INTERSPEC_TEST_TOKEN", nullptr );
  BOOST_CHECK_THROW( LlmConfig::loadApiAndMcpConfigs( path, &uses_env ), std::exception );

  SpecUtils::remove_file( path );
  SpecUtils::remove_file( resaved );
}//EnvReferencesInConfigFile


BOOST_AUTO_TEST_CASE( FilenameSanitizing )
{
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "spectrum.n42" ), "spectrum.n42" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "../../etc/passwd" ), "passwd" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "C:\\\\Users\\\\x\\\\my file (1).spe" ), "my_file__1_.spe" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( ".hidden" ), "hidden" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "..." ), "file" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "" ), "file" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "dir/" ), "file" );
  BOOST_CHECK_EQUAL( McpSessionManager::sanitizeFilename( "a\nb;rm -rf.n42" ), "a_b_rm_-rf.n42" );

  const string long_name = string( 300, 'x' ) + ".n42";
  const string sanitized = McpSessionManager::sanitizeFilename( long_name );
  BOOST_CHECK_EQUAL( sanitized.size(), 128u );
  BOOST_CHECK( SpecUtils::iends_with( sanitized, ".n42" ) );

  const fs::path dir = fs::path( SpecUtils::temp_dir() ) / ("mcp_unique_" + std::to_string( std::rand() ));
  fs::create_directories( dir );
  BOOST_CHECK_EQUAL( McpSessionManager::uniqueFilename( dir.string(), "a.n42" ), "a.n42" );
  ofstream( (dir / "a.n42").string() ) << "x";
  BOOST_CHECK_EQUAL( McpSessionManager::uniqueFilename( dir.string(), "a.n42" ), "a-1.n42" );
  ofstream( (dir / "a-1.n42").string() ) << "x";
  BOOST_CHECK_EQUAL( McpSessionManager::uniqueFilename( dir.string(), "a.n42" ), "a-2.n42" );
  fs::remove_all( dir );
}//FilenameSanitizing


BOOST_AUTO_TEST_CASE( SessionPathConfinement )
{
  const fs::path root = fs::canonical( SpecUtils::temp_dir() ) / ("mcp_paths_" + std::to_string( std::rand() ));
  const fs::path session = root / "2026-10-17_SessionA";
  const fs::path other = root / "2026-10-17_SessionB";
  fs::create_directories( session / "input" );
  fs::create_directories( other / "input" );
  ofstream( (session / "input" / "spec.n42").string() ) << "x";
  ofstream( (other / "input" / "secret.n42").string() ) << "x";

  const string base = session.string();

  // Allowed: anything inside the session directory, existing or not.
  BOOST_CHECK_EQUAL( McpSessionManager::resolvePathWithin( base, "input/spec.n42" ), (session / "input" / "spec.n42").string() );
  BOOST_CHECK_NO_THROW( McpSessionManager::resolvePathWithin( base, "output/new.n42" ) );
  BOOST_CHECK_NO_THROW( McpSessionManager::resolvePathWithin( base, "input/../input/spec.n42" ) );

  // Refused: empty, absolute, or escaping paths - including into another session.
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "" ), std::exception );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "/etc/passwd" ), std::exception );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, (session / "input" / "spec.n42").string() ), std::exception );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "../2026-10-17_SessionB/input/secret.n42" ), std::exception );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "input/../../2026-10-17_SessionB/input/secret.n42" ), std::exception );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "../../../../../../etc/hosts" ), std::exception );
  // A sibling whose name starts with the session's name is not inside it.
  fs::create_directories( root / "2026-10-17_SessionAX" );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "../2026-10-17_SessionAX" ), std::exception );

#if( !defined(_WIN32) )
  // A symbolic link inside the session that points outside of it is refused.
  fs::create_symlink( other / "input" / "secret.n42", session / "input" / "link.n42" );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "input/link.n42" ), std::exception );
  fs::create_directory_symlink( other, session / "input" / "linkdir" );
  BOOST_CHECK_THROW( McpSessionManager::resolvePathWithin( base, "input/linkdir/input/secret.n42" ), std::exception );
#endif

  fs::remove_all( root );
}//SessionPathConfinement


BOOST_AUTO_TEST_CASE( AnalysisN42RoundTrip )
{
  // A spectrum with user and automated-search peaks, an activity/shielding model, and an assistant
  //  conversation whose user message carries an image.
  shared_ptr<SpecUtils::Measurement> measurement = make_shared<SpecUtils::Measurement>();
  measurement->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 300.0f, 310.0f );
  shared_ptr<SpecUtils::EnergyCalibration> cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, { 0.0f, 3.0f }, {} );
  measurement->set_energy_calibration( cal );
  measurement->set_sample_number( 1 );

  SpecMeas meas;
  meas.add_measurement( measurement, true );
  const set<int> samples = meas.sample_numbers();
  BOOST_REQUIRE_EQUAL( samples.size(), 1u );

  deque<shared_ptr<const PeakDef>> user_peaks;
  user_peaks.push_back( make_shared<PeakDef>( 661.7, 1.5, 1000.0 ) );
  user_peaks.push_back( make_shared<PeakDef>( 1460.8, 2.0, 500.0 ) );
  meas.setPeaks( user_peaks, samples );

  shared_ptr<SpecMeas::PeakDeque> auto_peaks = make_shared<SpecMeas::PeakDeque>();
  auto_peaks->push_back( make_shared<PeakDef>( 609.3, 1.4, 300.0 ) );
  meas.setAutomatedSearchPeaks( samples, auto_peaks );

  {
    unique_ptr<rapidxml::xml_document<char>> model = make_unique<rapidxml::xml_document<char>>();
    rapidxml::xml_node<char> *node = model->allocate_node( rapidxml::node_element, "ShieldingSourceFit", "test-model" );
    model->append_node( node );
    meas.setShieldingSourceModel( std::move(model) );
  }

  LlmToolCall::ImageContent image;
  image.base64Data = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==";
  image.mimeType = "image/png";
  image.caption = "Photo A";
  image.widthPx = image.heightPx = 1;

  LlmConversationHistory history;
  history.addUserMessageToMainConversation( "What is in this spectrum?", { image } );
  meas.setLlmConversationHistory( samples, make_shared<vector<shared_ptr<LlmInteraction>>>( history.getConversations() ) );

  // A copy (as made by some export paths) keeps the automated-search peaks too.
  {
    SpecMeas copy;
    copy.uniqueCopyContents( meas );
    const shared_ptr<const SpecMeas::PeakDeque> copied = copy.automatedSearchPeaks( samples );
    BOOST_REQUIRE( copied );
    BOOST_CHECK_EQUAL( copied->size(), 1u );
  }

  const auto write_and_reload = [&meas]( SpecMeas &reloaded ){
    std::ostringstream n42;
    BOOST_REQUIRE( meas.write_2012_N42( n42 ) );
    const string path = write_temp_file( n42.str(), ".n42" );
    BOOST_REQUIRE( reloaded.load_file( path, SpecUtils::ParserType::N42_2012, "" ) );
    SpecUtils::remove_file( path );
  };

  // As written by an agent server: the conversation's images are kept.
  {
    const SaveUserImagesScope save_images;
    SpecMeas reloaded;
    write_and_reload( reloaded );

    const set<int> reloaded_samples = reloaded.sample_numbers();
    const shared_ptr<const SpecMeas::PeakDeque> peaks = reloaded.peaks( reloaded_samples );
    BOOST_REQUIRE( peaks );
    BOOST_CHECK_EQUAL( peaks->size(), 2u );

    const shared_ptr<const SpecMeas::PeakDeque> reloaded_auto = reloaded.automatedSearchPeaks( reloaded_samples );
    BOOST_REQUIRE( reloaded_auto );
    BOOST_CHECK_EQUAL( reloaded_auto->size(), 1u );

    BOOST_CHECK( reloaded.shieldingSourceModel() != nullptr );

    const shared_ptr<vector<shared_ptr<LlmInteraction>>> convo = reloaded.llmConversationHistory( reloaded_samples );
    BOOST_REQUIRE( convo && (convo->size() == 1u) );
    BOOST_REQUIRE( !convo->front()->responses.empty() );
    const LlmInteractionInitialRequest *request
                  = dynamic_cast<const LlmInteractionInitialRequest *>( convo->front()->responses.front().get() );
    BOOST_REQUIRE( request );
    BOOST_CHECK_EQUAL( request->content(), "What is in this spectrum?" );
    BOOST_REQUIRE_EQUAL( request->imageContent().size(), 1u );
    BOOST_CHECK_EQUAL( request->imageContent()[0].base64Data, image.base64Data );
    BOOST_CHECK_EQUAL( request->imageContent()[0].mimeType, "image/png" );
    BOOST_CHECK_EQUAL( request->imageContent()[0].caption, "Photo A" );
  }

  // As written by desktop InterSpec: unchanged from before - the images are not written.
  {
    SpecMeas reloaded;
    write_and_reload( reloaded );
    const shared_ptr<vector<shared_ptr<LlmInteraction>>> convo = reloaded.llmConversationHistory( reloaded.sample_numbers() );
    BOOST_REQUIRE( convo && (convo->size() == 1u) && !convo->front()->responses.empty() );
    const LlmInteractionInitialRequest *request
                  = dynamic_cast<const LlmInteractionInitialRequest *>( convo->front()->responses.front().get() );
    BOOST_REQUIRE( request );
    BOOST_CHECK( request->imageContent().empty() );
  }
}//AnalysisN42RoundTrip
