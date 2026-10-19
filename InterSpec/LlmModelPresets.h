#ifndef LLM_MODEL_PRESETS_H
#define LLM_MODEL_PRESETS_H
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

#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <optional>

#include <nlohmann/json_fwd.hpp>

#include "InterSpec/LlmConfig.h"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );


/** Known LLM providers and per-model default settings, used by the "LLM Provider Settings" window
 (LlmConfigWindow) to pre-fill an endpoint from a list of well-known providers, and to fill in a
 model's capabilities (image support, max tokens, context limit, reasoning, ...) when the user picks
 or types a model name.

 Read from "llm_model_presets.xml" (see `LlmConfig::configFilePath` - a copy in the writable data
 directory overrides the distributed one), so the list can be updated without a rebuild.  These are
 only conveniences, so a missing or invalid file just means no presets/defaults.
 */
struct LlmModelPresets
{
  static const int sm_xmlSerializationVersion = 1;

  typedef LlmConfig::LlmApi::ApiFormat ApiFormat;
  typedef LlmConfig::LlmApi::HttpBackend HttpBackend;
  typedef LlmConfig::LlmApi::ReasoningEffort ReasoningEffort;

  /** A well-known provider endpoint the user can select from a drop-down. */
  struct ProviderPreset
  {
    std::string label;          ///< Displayed to the user, e.g. "OpenAI (Chat Completions)"
    std::string endpoint;       ///< Full endpoint URL, e.g. "https://api.openai.com/v1/chat/completions"
    std::optional<ApiFormat> apiFormat;     ///< Unset means auto-detect from the URL
    std::optional<HttpBackend> httpBackend; ///< e.g. Native for providers that block browser CORS
    std::string suggestedModel; ///< Optional model to fill in if the provider has no model yet
  };//struct ProviderPreset

  /** Default settings for a model; only fields that are set get applied. */
  struct ModelDefaults
  {
    std::optional<bool> supportsImages;
    std::optional<std::variant<bool,ReasoningEffort>> reasoning;
    std::optional<int> maxTokens;
    std::optional<int> contextLengthLimit;
    std::optional<bool> requireToolInStateMachine;

    /** The only wire format the model works with (at a given provider); not a model setting, so
     `applyTo(...)` ignores it - it is only used to warn the user. */
    std::optional<LlmConfig::LlmApi::ApiFormat> requiredFormat;

    bool empty() const;

    /** Overwrites fields of this object with the fields that are set in `other`. */
    void mergeFrom( const ModelDefaults &other );

    /** The most a provider's "/models" listing can set max-tokens or the context limit to.

     Listings report a model's technical ceilings (e.g., 384K output tokens), which make poor
     request values: some providers check credit or rate limits against the requested max tokens,
     or require prompt + max tokens to fit in the context window.  These match the ceilings the
     llm_model_presets.xml table uses for the models it knows.
     */
    static constexpr int sm_maxFetchedMaxTokens = 64000;
    static constexpr int sm_maxFetchedContextLength = 200000;

    /** Combines the defaults from llm_model_presets.xml with what a provider's "/models" listing
     reported (either may be null).  The listing is authoritative for image/reasoning support, but
     for the context and max-output token counts the smaller value is used, since the table values
     may intentionally be below a model's true limit (e.g., to keep long conversations affordable);
     the listing's values are also capped at `sm_maxFetchedMaxTokens`/`sm_maxFetchedContextLength`.
     */
    static ModelDefaults combine( const ModelDefaults *table, const ModelDefaults *fetched );

    /** Sets the fields of `model` that are set in this object; never touches name, temperature,
     or instruction verbosity.  Returns true if any field of `model` changed.
     */
    bool applyTo( LlmConfig::LlmApi::ModelInfo &model ) const;
  };//struct ModelDefaults

  /** A model-name pattern, and the defaults for models that match it. */
  struct ModelDefaultsEntry
  {
    std::string match;     ///< Case-insensitive glob; '*' matches any sequence of characters
    /** If non-empty, a glob the provider's endpoint URL must match; such entries override some
     fields of the general (no `endpoint`) entry - e.g., for a provider with lower token limits. */
    std::string endpoint;
    ModelDefaults defaults;
  };//struct ModelDefaultsEntry

  std::vector<ProviderPreset> providerPresets;
  std::vector<ModelDefaultsEntry> modelDefaults;  ///< In file order; first match wins


  /** Loads "llm_model_presets.xml" from the writable or static data directory.
   Never throws; returns an empty object (and prints a warning) if the file is missing or invalid.
   */
  static std::shared_ptr<const LlmModelPresets> load();

  /** Parses the given file; throws std::exception on error. */
  static LlmModelPresets fromXmlFile( const std::string &filename );

  /** Returns the defaults for a model, or an empty optional if nothing is known about it.

   The defaults of the first general entry (one without an `endpoint`) whose pattern matches
   `modelId` are used, with the fields of the first matching endpoint-specific entry (whose
   `endpoint` glob also matches `endpoint`) overriding them.

   Each pattern is tried against the full id, and against the part after the last '/', so that
   "openai/gpt-oss-120b", "models/gemini-2.5-flash", or OpenRouter's "anthropic/claude-sonnet-4.5"
   can all be matched by patterns written for the bare model names.
   */
  std::optional<ModelDefaults> defaultsFor( const std::string &modelId,
                                            const std::string &endpoint = std::string() ) const;

  /** Case-insensitive glob match, where '*' matches any (possibly empty) sequence of characters. */
  static bool globMatch( const std::string &pattern, const std::string &str );

  /** Extracts what capability information an entry of a provider's "/models" listing gives.

   Providers vary in what they report: OpenRouter gives `context_length`,
   `top_provider.max_completion_tokens`, `architecture.input_modalities`, and `supported_parameters`;
   Google's native API gives `inputTokenLimit`, `outputTokenLimit`, and `thinking`; others give
   `max_input_tokens`/`max_tokens`/`context_window`, or nothing beyond the id.
   */
  static ModelDefaults defaultsFromModelsListEntry( const nlohmann::json &entry );
};//struct LlmModelPresets

#endif //LLM_MODEL_PRESETS_H
