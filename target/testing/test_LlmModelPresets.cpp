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

#include <set>
#include <limits>
#include <string>
#include <vector>
#include <variant>
#include <filesystem>
#include <system_error>

#define BOOST_TEST_MODULE LlmModelPresets_suite
#include <boost/test/included/unit_test.hpp>

#include <nlohmann/json.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmModelPresets.h"

using namespace std;

typedef LlmConfig::LlmApi::ReasoningEffort ReasoningEffort;
// (A typedef, as the comma in the template arguments would break the BOOST_CHECK macros.)
typedef std::variant<bool,ReasoningEffort> Reasoning_t;

namespace
{
// The static data directory, from the `--datadir=` argument, or else searched for relative to the CWD.
string find_data_dir()
{
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
    for( const char *d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
      {
        datadir = d;
        break;
      }
    }
  }//if( datadir.empty() )

  return datadir;
}//find_data_dir()


const LlmModelPresets &shipped_presets()
{
  static const LlmModelPresets presets = [](){
    const string datadir = find_data_dir();
    BOOST_REQUIRE_MESSAGE( !datadir.empty(), "Could not find the data directory; pass --datadir=..." );
    return LlmModelPresets::fromXmlFile( SpecUtils::append_path( datadir, "llm_model_presets.xml" ) );
  }();
  return presets;
}//shipped_presets()


Reasoning_t reasoning_for( const string &model )
{
  const std::optional<LlmModelPresets::ModelDefaults> d = shipped_presets().defaultsFor( model );
  BOOST_REQUIRE_MESSAGE( d.has_value(), "No defaults for '" << model << "'" );
  return d->reasoning.value_or( Reasoning_t(false) );
}
}//namespace


BOOST_AUTO_TEST_CASE( GlobMatch )
{
  BOOST_CHECK( LlmModelPresets::globMatch( "gpt-5*", "gpt-5-mini" ) );
  BOOST_CHECK( LlmModelPresets::globMatch( "GPT-5*", "gpt-5" ) );
  BOOST_CHECK( LlmModelPresets::globMatch( "*flash*", "gemini-2.5-flash-lite" ) );
  BOOST_CHECK( LlmModelPresets::globMatch( "*", "" ) );
  BOOST_CHECK( LlmModelPresets::globMatch( "a*b*c", "aXXbYYc" ) );
  BOOST_CHECK( !LlmModelPresets::globMatch( "a*b*c", "aXXbYY" ) );
  BOOST_CHECK( !LlmModelPresets::globMatch( "gpt-5", "gpt-5-mini" ) );
  BOOST_CHECK( !LlmModelPresets::globMatch( "claude-*", "anthropic/claude-x" ) );
}//GlobMatch


BOOST_AUTO_TEST_CASE( ShippedProviderPresets )
{
  const LlmModelPresets &presets = shipped_presets();
  BOOST_REQUIRE( presets.providerPresets.size() >= 8 );

  set<string> labels;
  for( const LlmModelPresets::ProviderPreset &p : presets.providerPresets )
  {
    BOOST_CHECK_MESSAGE( labels.insert( p.label ).second, "Duplicate preset label '" << p.label << "'" );
    BOOST_CHECK_MESSAGE( SpecUtils::istarts_with( p.endpoint, "https://" )
                         || SpecUtils::istarts_with( p.endpoint, "http://localhost:" )
                         || SpecUtils::istarts_with( p.endpoint, "http://127.0.0.1:" ),
                         "Unexpected endpoint '" << p.endpoint << "'" );
    BOOST_CHECK_MESSAGE( p.endpoint.find( ' ' ) == string::npos, "Space in endpoint '" << p.endpoint << "'" );

    // A suggested model should get filled-in defaults.
    if( !p.suggestedModel.empty() )
      BOOST_CHECK_MESSAGE( presets.defaultsFor( p.suggestedModel, p.endpoint ).has_value(),
                           "No defaults for suggested model '" << p.suggestedModel << "' of '" << p.label << "'" );
  }
}//ShippedProviderPresets


BOOST_AUTO_TEST_CASE( ShippedModelDefaults )
{
  const LlmModelPresets &presets = shipped_presets();
  BOOST_REQUIRE( !presets.modelDefaults.empty() );

  // Unknown models get nothing, so the user's settings are left alone.
  BOOST_CHECK( !presets.defaultsFor( "totally-unknown-model" ).has_value() );
  BOOST_CHECK( !presets.defaultsFor( "" ).has_value() );

  // Provider prefixes ("openai/", "models/", "anthropic/") are looked through.
  BOOST_CHECK( presets.defaultsFor( "openai/gpt-oss-120b" ).has_value() );
  BOOST_CHECK( presets.defaultsFor( "models/gemini-2.5-flash" ).has_value() );
  BOOST_CHECK( presets.defaultsFor( "anthropic/claude-sonnet-4.5" ).has_value() );
  BOOST_CHECK( presets.defaultsFor( "meta-llama/Llama-3.3-70B-Instruct" ).has_value() );

  // Frontier models get medium effort, smaller/flash models high.
  BOOST_CHECK( reasoning_for( "gpt-5" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "gpt-5-mini" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "gpt-5.4-mini" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "gpt-6-luna" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "gpt-6-sol" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "gpt-5-pro" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "claude-opus-5-5" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "claude-haiku-5-5" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "gemini-4-argon" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "gemini-3.8-flash" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "grok-4.7" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "grok-4.3" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "grok-4.20-0309-reasoning" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "claude-sonnet-4-5" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "claude-opus-4-6" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "claude-haiku-4-5" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "gemini-2.5-pro" ) == Reasoning_t( ReasoningEffort::medium ) );
  BOOST_CHECK( reasoning_for( "gemini-2.5-flash" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "gemini-2.5-flash-lite" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "openai/gpt-oss-120b" ) == Reasoning_t( ReasoningEffort::high ) );

  // Non-reasoning models must not get reasoning turned on, or the provider may reject requests.
  BOOST_CHECK( reasoning_for( "gpt-4.1" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "gpt-4o-mini" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "gpt-5-chat-latest" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "gpt-5.2-chat-latest" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "gemini-2.0-flash" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "gemini-2.0-flash-lite" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "claude-3-5-haiku-20241022" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "anthropic/claude-3.5-sonnet" ) == Reasoning_t( false ) );
  BOOST_CHECK( reasoning_for( "claude-3-7-sonnet-latest" ) == Reasoning_t( ReasoningEffort::medium ) );

  // Grok models that reject reasoning_effort (most of them reason anyway).
  for( const char *model : { "grok-4", "grok-4-0709", "x-ai/grok-4", "grok-4-fast-reasoning",
                             "grok-4-1-fast-non-reasoning", "grok-code-fast-1", "grok-3" } )
    BOOST_CHECK_MESSAGE( reasoning_for( model ) == Reasoning_t( false ), "Reasoning set for '" << model << "'" );
  BOOST_CHECK( reasoning_for( "grok-3-mini" ) == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK( reasoning_for( "grok-4.3-fast" ) == Reasoning_t( ReasoningEffort::high ) );

  // Other vendors' flash-tier models just get reasoning turned on, with no other settings.
  BOOST_CHECK( reasoning_for( "somevendor/somevendor-v4-flash" ) == Reasoning_t( true ) );
  BOOST_CHECK( !presets.defaultsFor( "somevendor/somevendor-v4-flash" )->maxTokens.has_value() );

  // At OpenAI itself, some models need the Responses API; other providers may translate.
  typedef LlmConfig::LlmApi::ApiFormat ApiFormat;
  const string openaiChat = "https://api.openai.com/v1/chat/completions";
  for( const char *model : { "gpt-6-astra", "gpt-6.1-sol", "gpt-5-pro" } )
  {
    const std::optional<LlmModelPresets::ModelDefaults> atOpenAi = presets.defaultsFor( model, openaiChat );
    BOOST_REQUIRE( atOpenAi.has_value() );
    BOOST_CHECK_MESSAGE( atOpenAi->requiredFormat == std::optional<ApiFormat>( ApiFormat::OpenAiResponses ),
                         "No Responses requirement for '" << model << "'" );
    BOOST_CHECK( atOpenAi->reasoning.has_value() );  // the general entry's fields are still there
    BOOST_CHECK( !presets.defaultsFor( model, "https://openrouter.ai/api/v1/chat/completions" )->requiredFormat.has_value() );
  }
  BOOST_CHECK( !presets.defaultsFor( "gpt-5-mini", openaiChat )->requiredFormat.has_value() );

  // Old Claude models have small output limits.
  BOOST_CHECK_EQUAL( presets.defaultsFor( "claude-3-haiku-20240307" )->maxTokens.value_or(0), 4096 );
  BOOST_CHECK_EQUAL( presets.defaultsFor( "claude-3-5-haiku-latest" )->maxTokens.value_or(0), 8192 );

  // Image support.
  BOOST_CHECK( presets.defaultsFor( "gpt-5-mini" )->supportsImages.value_or( false ) );
  BOOST_CHECK( presets.defaultsFor( "claude-sonnet-4-5" )->supportsImages.value_or( false ) );
  BOOST_CHECK( !presets.defaultsFor( "openai/gpt-oss-120b" )->supportsImages.value_or( true ) );

  // Anthropic requires max_tokens, so every Claude entry must give one.
  for( const char *model : { "claude-opus-4-6", "claude-sonnet-4-5", "claude-haiku-4-5", "claude-fable-5-1" } )
  {
    const std::optional<LlmModelPresets::ModelDefaults> d = presets.defaultsFor( model );
    BOOST_REQUIRE_MESSAGE( d.has_value(), "No defaults for '" << model << "'" );
    BOOST_CHECK_MESSAGE( d->maxTokens.value_or( 0 ) > 0, "No max tokens for '" << model << "'" );
  }

  // Every entry should be reachable: no earlier, more general pattern may completely shadow a later one.
  //  We check this by matching each pattern's own text (with '*' removed) against the earlier patterns.
  const vector<LlmModelPresets::ModelDefaultsEntry> &entries = presets.modelDefaults;
  for( size_t i = 0; i < entries.size(); ++i )
  {
    string sample = entries[i].match;
    SpecUtils::ireplace_all( sample, "*", "" );
    for( size_t j = 0; j < i; ++j )
      BOOST_CHECK_MESSAGE( (entries[j].endpoint != entries[i].endpoint)
                           || !LlmModelPresets::globMatch( entries[j].match, sample ),
                           "Pattern '" << entries[i].match << "' is shadowed by earlier '" << entries[j].match << "'" );
  }

  // Endpoint-specific entries override just their fields of the general entry.
  const string asksage = "https://api.nnsa.asksage.ai/server/openai/v1/chat/completions";
  const std::optional<LlmModelPresets::ModelDefaults> gemGeneral = presets.defaultsFor( "google-gemini-2.5-pro" );
  const std::optional<LlmModelPresets::ModelDefaults> gemAskSage = presets.defaultsFor( "google-gemini-2.5-pro", asksage );
  BOOST_REQUIRE( gemGeneral.has_value() && gemAskSage.has_value() );
  BOOST_CHECK_EQUAL( gemAskSage->maxTokens.value_or(0), 8192 );
  BOOST_CHECK( gemGeneral->maxTokens.value_or(0) > 8192 );
  BOOST_CHECK( gemAskSage->reasoning == gemGeneral->reasoning );
  BOOST_CHECK( gemAskSage->supportsImages == gemGeneral->supportsImages );
  BOOST_CHECK_EQUAL( presets.defaultsFor( "google-gemini-2.5-pro", "https://openrouter.ai/api/v1/chat/completions" )->maxTokens.value_or(0),
                     gemGeneral->maxTokens.value_or(0) );
}//ShippedModelDefaults


BOOST_AUTO_TEST_CASE( ApplyAndCombine )
{
  LlmConfig::LlmApi::ModelInfo model;
  model.name = "some-model";
  model.temperature = 0.3;
  model.instructionVerbosity = LlmConfig::LlmApi::InstructionVerbosity::Verbose;

  LlmModelPresets::ModelDefaults table;
  table.supportsImages = true;
  table.reasoning = ReasoningEffort::high;
  table.maxTokens = 32000;
  table.contextLengthLimit = 200000;

  BOOST_CHECK( table.applyTo( model ) );
  BOOST_CHECK( model.supportsImages );
  BOOST_CHECK( model.reasoning == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK_EQUAL( model.maxTokens, 32000 );
  BOOST_CHECK_EQUAL( model.contextLengthLimit, 200000 );
  BOOST_CHECK( !model.requireToolInStateMachine );  // not set in the defaults, so untouched
  BOOST_CHECK_EQUAL( model.name, "some-model" );
  BOOST_REQUIRE( model.temperature.has_value() );
  BOOST_CHECK_CLOSE( model.temperature.value(), 0.3, 1.0E-9 );
  BOOST_CHECK( model.instructionVerbosity == LlmConfig::LlmApi::InstructionVerbosity::Verbose );
  BOOST_CHECK( !table.applyTo( model ) );  // nothing changes the second time

  // The listing is authoritative for capabilities, but the smaller token counts win.
  LlmModelPresets::ModelDefaults fetched;
  fetched.supportsImages = false;
  fetched.contextLengthLimit = 1048576;
  fetched.maxTokens = 16000;
  const LlmModelPresets::ModelDefaults both = LlmModelPresets::ModelDefaults::combine( &table, &fetched );
  BOOST_CHECK( !both.supportsImages.value() );
  BOOST_CHECK( both.reasoning == Reasoning_t( ReasoningEffort::high ) );
  BOOST_CHECK_EQUAL( both.contextLengthLimit.value(), 200000 );
  BOOST_CHECK_EQUAL( both.maxTokens.value(), 16000 );

  // A model the table does not know: the listing's technical ceilings are capped to sane request
  //  values, but smaller listed values are kept.
  LlmModelPresets::ModelDefaults bigFetched;
  bigFetched.contextLengthLimit = 1048576;
  bigFetched.maxTokens = 393216;
  const LlmModelPresets::ModelDefaults onlyFetched = LlmModelPresets::ModelDefaults::combine( nullptr, &bigFetched );
  BOOST_CHECK_EQUAL( onlyFetched.contextLengthLimit.value(), LlmModelPresets::ModelDefaults::sm_maxFetchedContextLength );
  BOOST_CHECK_EQUAL( onlyFetched.maxTokens.value(), LlmModelPresets::ModelDefaults::sm_maxFetchedMaxTokens );
  BOOST_CHECK( !onlyFetched.reasoning.has_value() );

  LlmModelPresets::ModelDefaults smallFetched;
  smallFetched.contextLengthLimit = 32768;
  smallFetched.maxTokens = 8192;
  const LlmModelPresets::ModelDefaults smallOnly = LlmModelPresets::ModelDefaults::combine( nullptr, &smallFetched );
  BOOST_CHECK_EQUAL( smallOnly.contextLengthLimit.value(), 32768 );
  BOOST_CHECK_EQUAL( smallOnly.maxTokens.value(), 8192 );
  BOOST_CHECK( LlmModelPresets::ModelDefaults::combine( nullptr, nullptr ).empty() );
}//ApplyAndCombine


BOOST_AUTO_TEST_CASE( ModelsListingCapabilities )
{
  // OpenRouter-style entry
  const nlohmann::json openrouter = nlohmann::json::parse( R"({
    "id": "google/gemini-2.5-flash", "context_length": 1048576,
    "architecture": { "input_modalities": ["text", "image"] },
    "top_provider": { "max_completion_tokens": 65535 },
    "supported_parameters": ["tools", "reasoning", "max_tokens"] })" );
  const LlmModelPresets::ModelDefaults orCaps = LlmModelPresets::defaultsFromModelsListEntry( openrouter );
  BOOST_CHECK_EQUAL( orCaps.contextLengthLimit.value_or(0), 1048576 );
  BOOST_CHECK_EQUAL( orCaps.maxTokens.value_or(0), 65535 );
  BOOST_CHECK( orCaps.supportsImages.value_or( false ) );
  BOOST_CHECK( !orCaps.reasoning.has_value() );  // supported - effort level comes from the table

  const nlohmann::json noReasoning = nlohmann::json::parse( R"({
    "id": "meta-llama/llama-3.3-70b-instruct", "context_length": 131072,
    "architecture": { "input_modalities": ["text"] }, "supported_parameters": ["tools"] })" );
  const LlmModelPresets::ModelDefaults nrCaps = LlmModelPresets::defaultsFromModelsListEntry( noReasoning );
  BOOST_CHECK( !nrCaps.supportsImages.value_or( true ) );
  BOOST_CHECK( nrCaps.reasoning == std::optional<Reasoning_t>( false ) );

  // Google native-style entry
  const nlohmann::json google = nlohmann::json::parse( R"({
    "name": "models/gemini-2.5-pro", "inputTokenLimit": 1048576, "outputTokenLimit": 65536, "thinking": true })" );
  const LlmModelPresets::ModelDefaults gCaps = LlmModelPresets::defaultsFromModelsListEntry( google );
  BOOST_CHECK_EQUAL( gCaps.contextLengthLimit.value_or(0), 1048576 );
  BOOST_CHECK_EQUAL( gCaps.maxTokens.value_or(0), 65536 );

  // vLLM / LM Studio give the context length under their own names.
  const nlohmann::json vllm = nlohmann::json::parse( R"({ "id": "Qwen/Qwen3-32B", "max_model_len": 40960 })" );
  BOOST_CHECK_EQUAL( LlmModelPresets::defaultsFromModelsListEntry( vllm ).contextLengthLimit.value_or(0), 40960 );
  const nlohmann::json lmstudio = nlohmann::json::parse( R"({ "id": "qwen3-8b", "max_context_length": 32768 })" );
  BOOST_CHECK_EQUAL( LlmModelPresets::defaultsFromModelsListEntry( lmstudio ).contextLengthLimit.value_or(0), 32768 );

  // Out-of-range or oddly typed values must neither overflow nor throw.
  const nlohmann::json huge = nlohmann::json::parse( R"({ "id": "x", "context_length": 1e30,
    "max_tokens": 18446744073709551615, "architecture": "text", "supported_parameters": "reasoning" })" );
  LlmModelPresets::ModelDefaults hugeCaps;
  BOOST_REQUIRE_NO_THROW( hugeCaps = LlmModelPresets::defaultsFromModelsListEntry( huge ) );
  BOOST_CHECK_EQUAL( hugeCaps.contextLengthLimit.value_or(0), std::numeric_limits<int>::max() );
  BOOST_CHECK_EQUAL( hugeCaps.maxTokens.value_or(0), std::numeric_limits<int>::max() );
  const nlohmann::json negative = nlohmann::json::parse( R"({ "id": "x", "context_length": -5, "max_tokens": -1.5 })" );
  BOOST_CHECK( LlmModelPresets::defaultsFromModelsListEntry( negative ).empty() );

  // OpenAI-style entry: only an id, so nothing is known.
  const nlohmann::json openai = nlohmann::json::parse( R"({ "id": "gpt-5-mini", "object": "model", "owned_by": "openai" })" );
  BOOST_CHECK( LlmModelPresets::defaultsFromModelsListEntry( openai ).empty() );
}//ModelsListingCapabilities


BOOST_AUTO_TEST_CASE( LoadWithoutUserConfig )
{
  // InterSpec no longer ships a llm_config.xml; loading without one must not throw.
  const string datadir = find_data_dir();
  BOOST_REQUIRE( !datadir.empty() );
  BOOST_REQUIRE( !SpecUtils::is_file( SpecUtils::append_path( datadir, "llm_config.xml" ) ) );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  // An empty, temporary, writable directory, so nothing can be written into the repository.
  const string writabledir = SpecUtils::temp_file_name( "test_LlmModelPresets", SpecUtils::temp_dir() );
  BOOST_REQUIRE( SpecUtils::create_directory( writabledir ) == 1 );
  BOOST_REQUIRE_NO_THROW( InterSpec::setWritableDataDirectory( writabledir ) );

  std::shared_ptr<LlmConfig> config;
  BOOST_REQUIRE_NO_THROW( config = LlmConfig::load() );
  BOOST_REQUIRE( config );

  // Unit-test builds stub in a provider so the tool tests have a config to work with.
  BOOST_CHECK( config->llmApi.isConfigured() );
  BOOST_CHECK( !config->mcpServer.enabled );

  // A default-constructed config (what non-test builds get) is not configured, and does not throw.
  const LlmConfig::LlmApi empty_api;
  BOOST_CHECK( !empty_api.isConfigured() );

  // Removing fails if anything was written into the directory.
  std::error_code ec;
  BOOST_CHECK_MESSAGE( std::filesystem::remove( writabledir, ec ), "Something was written to " << writabledir );
}//LoadWithoutUserConfig
