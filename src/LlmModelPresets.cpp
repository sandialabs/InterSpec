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

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <rapidxml/rapidxml.hpp>

#include <nlohmann/json.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/LlmModelPresets.h"

using namespace std;

namespace
{
  /** Parses an optional boolean attribute; an unrecognized value is ignored (with a warning). */
  std::optional<bool> bool_attrib( const rapidxml::xml_node<char> *node, const char *name )
  {
    const rapidxml::xml_attribute<char> * const attrib = node->first_attribute( name );
    if( !attrib )
      return std::nullopt;

    string value = SpecUtils::xml_value_str( attrib );
    SpecUtils::trim( value );
    SpecUtils::to_lower_ascii( value );
    if( (value == "true") || (value == "1") || (value == "yes") )
      return true;
    if( (value == "false") || (value == "0") || (value == "no") )
      return false;

    cerr << "llm_model_presets.xml: invalid boolean '" << value << "' for " << name << "; ignoring." << endl;
    return std::nullopt;
  }//bool_attrib(...)


  /** Parses an optional, positive, integer attribute; an invalid value is ignored (with a warning). */
  std::optional<int> int_attrib( const rapidxml::xml_node<char> *node, const char *name )
  {
    const rapidxml::xml_attribute<char> * const attrib = node->first_attribute( name );
    if( !attrib )
      return std::nullopt;

    string value = SpecUtils::xml_value_str( attrib );
    SpecUtils::trim( value );
    int result = 0;
    if( SpecUtils::parse_int( value.c_str(), value.size(), result ) && (result > 0) )
      return result;

    cerr << "llm_model_presets.xml: invalid integer '" << value << "' for " << name << "; ignoring." << endl;
    return std::nullopt;
  }//int_attrib(...)


  /** Returns a positive integer from `j[key]`, if it is present as a number (or numeric string). */
  std::optional<int> json_positive_int( const nlohmann::json &j, const char *key )
  {
    if( !j.is_object() || !j.contains(key) )
      return std::nullopt;

    // Read wide, and clamp, so out-of-range values cannot overflow the int.
    const nlohmann::json &v = j[key];
    int value = 0;
    if( v.is_number_unsigned() )
      value = static_cast<int>( std::min<uint64_t>( v.get<uint64_t>(), std::numeric_limits<int>::max() ) );
    else if( v.is_number_integer() )
      value = static_cast<int>( std::clamp<int64_t>( v.get<int64_t>(), 0, std::numeric_limits<int>::max() ) );
    else if( v.is_number_float() )
    {
      const double dval = v.get<double>();
      if( std::isfinite(dval) )
        value = static_cast<int>( std::clamp<double>( dval, 0.0, std::numeric_limits<int>::max() ) );
    }else if( v.is_string() )
    {
      const string str = v.get<string>();
      if( !SpecUtils::parse_int( str.c_str(), str.size(), value ) )
        return std::nullopt;
    }

    if( value > 0 )
      return value;
    return std::nullopt;
  }//json_positive_int(...)


  /** Returns true if `arr` is a JSON array containing (case-insensitively) the string `value`. */
  bool json_array_contains( const nlohmann::json &arr, const string &value )
  {
    if( !arr.is_array() )
      return false;

    for( const nlohmann::json &el : arr )
    {
      if( el.is_string() && SpecUtils::iequals_ascii( el.get<string>(), value ) )
        return true;
    }
    return false;
  }//json_array_contains(...)
}//namespace


bool LlmModelPresets::ModelDefaults::empty() const
{
  return !supportsImages.has_value() && !reasoning.has_value() && !maxTokens.has_value()
         && !contextLengthLimit.has_value() && !requireToolInStateMachine.has_value()
         && !requiredFormat.has_value();
}//ModelDefaults::empty()


void LlmModelPresets::ModelDefaults::mergeFrom( const ModelDefaults &other )
{
  if( other.supportsImages.has_value() )
    supportsImages = other.supportsImages;
  if( other.reasoning.has_value() )
    reasoning = other.reasoning;
  if( other.maxTokens.has_value() )
    maxTokens = other.maxTokens;
  if( other.contextLengthLimit.has_value() )
    contextLengthLimit = other.contextLengthLimit;
  if( other.requireToolInStateMachine.has_value() )
    requireToolInStateMachine = other.requireToolInStateMachine;
  if( other.requiredFormat.has_value() )
    requiredFormat = other.requiredFormat;
}//ModelDefaults::mergeFrom(...)


LlmModelPresets::ModelDefaults LlmModelPresets::ModelDefaults::combine( const ModelDefaults *table,
                                                                         const ModelDefaults *fetched )
{
  ModelDefaults result = table ? *table : ModelDefaults();
  if( !fetched )
    return result;

  const auto smaller = []( const std::optional<int> &a, const std::optional<int> &b ) -> std::optional<int> {
    if( a.has_value() && b.has_value() )
      return std::min( a.value(), b.value() );
    return a.has_value() ? a : b;
  };

  if( fetched->supportsImages.has_value() )
    result.supportsImages = fetched->supportsImages;
  // The listing can only tell us reasoning is unsupported; the effort level comes from the table.
  if( fetched->reasoning.has_value() )
    result.reasoning = fetched->reasoning;
  if( fetched->requireToolInStateMachine.has_value() )
    result.requireToolInStateMachine = fetched->requireToolInStateMachine;
  // The listing gives technical ceilings, so cap them; they can lower, but never raise, the
  //  table's suggested values.
  const auto capped = []( const std::optional<int> &value, const int cap ) -> std::optional<int> {
    return value.has_value() ? std::optional<int>( std::min( value.value(), cap ) ) : value;
  };
  result.maxTokens = smaller( result.maxTokens, capped( fetched->maxTokens, sm_maxFetchedMaxTokens ) );
  result.contextLengthLimit = smaller( result.contextLengthLimit,
                                       capped( fetched->contextLengthLimit, sm_maxFetchedContextLength ) );

  return result;
}//ModelDefaults::combine(...)


bool LlmModelPresets::ModelDefaults::applyTo( LlmConfig::LlmApi::ModelInfo &model ) const
{
  const LlmConfig::LlmApi::ModelInfo orig = model;

  if( supportsImages.has_value() )
    model.supportsImages = supportsImages.value();
  if( reasoning.has_value() )
    model.reasoning = reasoning.value();
  if( maxTokens.has_value() )
    model.maxTokens = maxTokens.value();
  if( contextLengthLimit.has_value() )
    model.contextLengthLimit = contextLengthLimit.value();
  if( requireToolInStateMachine.has_value() )
    model.requireToolInStateMachine = requireToolInStateMachine.value();

  return (orig.supportsImages != model.supportsImages)
         || (orig.reasoning != model.reasoning)
         || (orig.maxTokens != model.maxTokens)
         || (orig.contextLengthLimit != model.contextLengthLimit)
         || (orig.requireToolInStateMachine != model.requireToolInStateMachine);
}//ModelDefaults::applyTo(...)


std::shared_ptr<const LlmModelPresets> LlmModelPresets::load()
{
  string path;
  try
  {
    path = LlmConfig::configFilePath( "llm_model_presets.xml" );
    if( SpecUtils::is_file( path ) )
      return make_shared<const LlmModelPresets>( fromXmlFile( path ) );

    cerr << "LLM model presets file '" << path << "' not found; no provider presets or model defaults." << endl;
  }catch( std::exception &e )
  {
    cerr << "Failed to load LLM model presets from '" << path << "': " << e.what() << endl;
  }

  return make_shared<const LlmModelPresets>();
}//LlmModelPresets::load()


LlmModelPresets LlmModelPresets::fromXmlFile( const std::string &filename )
{
  std::vector<char> data;
  SpecUtils::load_file_data( filename.c_str(), data );

  rapidxml::xml_document<char> doc;
  doc.parse<rapidxml::parse_trim_whitespace>( &data[0] );

  const rapidxml::xml_node<char> * const root = XML_FIRST_NODE( &doc, "LlmModelPresets" );
  if( !root )
    throw runtime_error( "Missing LlmModelPresets root node" );

  XmlUtils::check_xml_version( root, sm_xmlSerializationVersion );

  LlmModelPresets presets;

  const rapidxml::xml_node<char> * const providers_node = XML_FIRST_NODE( root, "ProviderPresets" );
  XML_FOREACH_CHILD( preset_node, providers_node, "Preset" )
  {
    ProviderPreset preset;
    const rapidxml::xml_attribute<char> * const label_attrib = XML_FIRST_ATTRIB( preset_node, "label" );
    preset.label = SpecUtils::xml_value_str( label_attrib );
    preset.endpoint = XmlUtils::get_string_node_value( preset_node, "ApiEndpoint" );
    preset.suggestedModel = SpecUtils::xml_value_str( XML_FIRST_NODE(preset_node, "SuggestedModel") );
    SpecUtils::trim( preset.label );
    SpecUtils::trim( preset.endpoint );
    SpecUtils::trim( preset.suggestedModel );

    const rapidxml::xml_attribute<char> * const fmt_attrib = XML_FIRST_ATTRIB( preset_node, "apiFormat" );
    if( fmt_attrib )
      preset.apiFormat = LlmConfig::LlmApi::parseApiFormat( SpecUtils::xml_value_str(fmt_attrib) );

    const rapidxml::xml_attribute<char> * const http_attrib = XML_FIRST_ATTRIB( preset_node, "httpBackend" );
    if( http_attrib )
      preset.httpBackend = LlmConfig::LlmApi::parseHttpBackend( SpecUtils::xml_value_str(http_attrib) );

    if( preset.label.empty() || preset.endpoint.empty() )
      throw runtime_error( "A <Preset> is missing its label or <ApiEndpoint>." );

    presets.providerPresets.push_back( std::move(preset) );
  }//XML_FOREACH_CHILD( preset_node, providers_node, "Preset" )

  const rapidxml::xml_node<char> * const models_node = XML_FIRST_NODE( root, "ModelDefaults" );
  XML_FOREACH_CHILD( model_node, models_node, "Model" )
  {
    ModelDefaultsEntry entry;
    const rapidxml::xml_attribute<char> * const match_attrib = XML_FIRST_ATTRIB( model_node, "match" );
    entry.match = SpecUtils::xml_value_str( match_attrib );
    SpecUtils::trim( entry.match );
    if( entry.match.empty() )
      throw runtime_error( "A <Model> under <ModelDefaults> is missing its 'match' attribute." );

    const rapidxml::xml_attribute<char> * const endpoint_attrib = XML_FIRST_ATTRIB( model_node, "endpoint" );
    entry.endpoint = SpecUtils::xml_value_str( endpoint_attrib );
    SpecUtils::trim( entry.endpoint );

    ModelDefaults &defaults = entry.defaults;
    defaults.supportsImages = bool_attrib( model_node, "supportsImages" );
    defaults.maxTokens = int_attrib( model_node, "maxTokens" );
    defaults.contextLengthLimit = int_attrib( model_node, "contextLength" );
    defaults.requireToolInStateMachine = bool_attrib( model_node, "requireToolInStateMachine" );

    const rapidxml::xml_attribute<char> * const format_attrib = XML_FIRST_ATTRIB( model_node, "requiresFormat" );
    if( format_attrib )
    {
      try
      {
        defaults.requiredFormat = LlmConfig::LlmApi::parseApiFormat( SpecUtils::xml_value_str(format_attrib) );
      }catch( std::exception &e )
      {
        cerr << "llm_model_presets.xml: invalid requiresFormat for '" << entry.match << "': " << e.what() << endl;
      }
    }

    const rapidxml::xml_attribute<char> * const reasoning_attrib = XML_FIRST_ATTRIB( model_node, "reasoning" );
    if( reasoning_attrib )
    {
      defaults.reasoning = LlmConfig::LlmApi::parseReasoning( SpecUtils::xml_value_str(reasoning_attrib) );
      if( !defaults.reasoning.has_value() )
        cerr << "llm_model_presets.xml: invalid reasoning value for '" << entry.match << "'; ignoring." << endl;
    }

    presets.modelDefaults.push_back( std::move(entry) );
  }//XML_FOREACH_CHILD( model_node, models_node, "Model" )

  return presets;
}//LlmModelPresets::fromXmlFile(...)


std::optional<LlmModelPresets::ModelDefaults> LlmModelPresets::defaultsFor( const std::string &modelId,
                                                                             const std::string &endpoint ) const
{
  string id = modelId;
  SpecUtils::trim( id );
  if( id.empty() )
    return std::nullopt;

  const size_t slash_pos = id.find_last_of( '/' );
  const string base_id = (slash_pos == string::npos) ? string() : id.substr( slash_pos + 1 );
  const auto matches = [&id,&base_id]( const string &pattern ) -> bool {
    return globMatch( pattern, id ) || (!base_id.empty() && globMatch( pattern, base_id ));
  };

  const ModelDefaults *general = nullptr, *specific = nullptr;
  for( const ModelDefaultsEntry &entry : modelDefaults )
  {
    if( entry.endpoint.empty() )
    {
      if( !general && matches( entry.match ) )
        general = &entry.defaults;
    }else if( !specific && !endpoint.empty() && globMatch( entry.endpoint, endpoint ) && matches( entry.match ) )
    {
      specific = &entry.defaults;
    }
  }//for( const ModelDefaultsEntry &entry : modelDefaults )

  if( !general && !specific )
    return std::nullopt;

  ModelDefaults result = general ? *general : ModelDefaults();
  if( specific )
    result.mergeFrom( *specific );
  return result;
}//LlmModelPresets::defaultsFor(...)


bool LlmModelPresets::globMatch( const std::string &pattern, const std::string &str )
{
  // Iterative wildcard match with single-star backtracking; O(n*m) worst case, but these are short.
  size_t p = 0, s = 0, star_p = string::npos, star_s = 0;

  const auto lower = []( const char c ) -> char {
    return ((c >= 'A') && (c <= 'Z')) ? static_cast<char>(c - 'A' + 'a') : c;
  };

  while( s < str.size() )
  {
    if( (p < pattern.size()) && (pattern[p] == '*') )
    {
      star_p = p++;
      star_s = s;
    }else if( (p < pattern.size()) && (lower(pattern[p]) == lower(str[s])) )
    {
      ++p;
      ++s;
    }else if( star_p != string::npos )
    {
      p = star_p + 1;
      s = ++star_s;
    }else
    {
      return false;
    }
  }//while( s < str.size() )

  while( (p < pattern.size()) && (pattern[p] == '*') )
    ++p;

  return (p == pattern.size());
}//LlmModelPresets::globMatch(...)


LlmModelPresets::ModelDefaults LlmModelPresets::defaultsFromModelsListEntry( const nlohmann::json &entry )
{
  ModelDefaults defaults;
  if( !entry.is_object() )
    return defaults;

  // Context window: OpenRouter "context_length", Google native "inputTokenLimit", Anthropic
  //  "max_input_tokens", vLLM "max_model_len", LM Studio "max_context_length", others.
  for( const char *key : { "context_length", "inputTokenLimit", "max_input_tokens", "max_model_len",
                           "max_context_length", "context_window" } )
  {
    defaults.contextLengthLimit = json_positive_int( entry, key );
    if( defaults.contextLengthLimit.has_value() )
      break;
  }

  // Max output tokens: OpenRouter "top_provider.max_completion_tokens", Google "outputTokenLimit",
  //  Anthropic "max_tokens", others.
  if( entry.contains("top_provider") )
    defaults.maxTokens = json_positive_int( entry["top_provider"], "max_completion_tokens" );
  for( const char *key : { "outputTokenLimit", "max_output_tokens", "max_completion_tokens", "max_tokens" } )
  {
    if( defaults.maxTokens.has_value() )
      break;
    defaults.maxTokens = json_positive_int( entry, key );
  }

  // Image input: OpenRouter lists input modalities.
  if( entry.contains("architecture") && entry["architecture"].is_object()
      && entry["architecture"].contains("input_modalities") )
  {
    defaults.supportsImages = json_array_contains( entry["architecture"]["input_modalities"], "image" );
  }

  // Reasoning: OpenRouter lists the request parameters it supports; Google gives a "thinking" bool.
  //  We only learn that reasoning is *possible*, so leave the effort level to the preset table.
  if( entry.contains("supported_parameters") && entry["supported_parameters"].is_array()
      && !json_array_contains( entry["supported_parameters"], "reasoning" ) )
  {
    defaults.reasoning = false;
  }
  if( entry.contains("thinking") && entry["thinking"].is_boolean() && !entry["thinking"].get<bool>() )
    defaults.reasoning = false;

  return defaults;
}//LlmModelPresets::defaultsFromModelsListEntry(...)
