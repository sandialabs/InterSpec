#ifndef LLM_API_PROTOCOL_H
#define LLM_API_PROTOCOL_H
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

#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <optional>

#include <nlohmann/json.hpp>

#include "InterSpec/LlmConfig.h"  // for LlmConfig::LlmApi and LlmConfig::LlmApi::ApiFormat

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

// Forward declarations
struct LlmInteraction;

/** How the model should be allowed/forced to call tools for a request. */
enum class ToolChoice
{
  Auto,     // Model decides whether to call a tool (OpenAI "auto" / Anthropic {type:"auto"}).
  Required  // Model must call a tool (OpenAI "required" / Anthropic {type:"any"}).
};//enum class ToolChoice

/** A provider-agnostic tool definition, as produced from the InterSpec tool registry.

 Each protocol translates this into its own on-the-wire shape (OpenAI nests under "function";
 OpenAI Responses flattens; Anthropic uses "input_schema").
 */
struct NormalizedTool
{
  std::string name;
  std::string description;
  nlohmann::json parameters;  // JSON-schema object describing the tool's parameters
};//struct NormalizedTool

/** A normalized view of an LLM response, independent of the provider wire format.

 `LlmApiProtocol::parseResponse()` fills this in; `LlmInterface::handleApiResponse()` then drives
 the (format-agnostic) orchestration from it.  Note that transport- and API-level *errors* are not
 represented here: they are detected generically (top-level "error" field) before parseResponse()
 is ever called, and malformed payloads cause parseResponse() to throw (also handled upstream).
 */
struct ParsedLlmResponse
{
  struct ToolCall
  {
    std::string id;        // The provider's tool-call id (OpenAI id / Anthropic toolu_ / Responses call_)
    std::string name;      // Tool name
    std::string argsRaw;   // Tool arguments as a JSON *string* (passed through to the tool executor verbatim)
    nlohmann::json extra;  // Optional provider-specific passthrough (e.g. OpenAI "extra_content"); null if absent
  };//struct ToolCall

  std::string content;            // Assistant visible text (concatenated across text blocks)
  std::string thinkingContent;    // Human-readable reasoning/thinking text (for display)
  std::string thinkingSignature;  // Anthropic thinking-block signature (must be echoed back unchanged); empty otherwise
  std::string reasoningContent;   // OpenAI-style `reasoning_content` passthrough
  std::string reasoningDetails;   // OpenRouter `reasoning_details` array, as a JSON string (echoed back unchanged)
  std::vector<ToolCall> toolCalls;

  std::optional<size_t> promptTokens, completionTokens, totalTokens;
  std::optional<size_t> cachedTokens;          // Prompt tokens served from cache (cache reads)
  std::optional<size_t> cacheCreationTokens;   // Prompt tokens written to cache (Anthropic cache_creation_input_tokens)
  std::string stopReason;               // Normalized stop reason ("" if unknown/unavailable)
};//struct ParsedLlmResponse


/** Abstracts the on-the-wire translation for a single LLM API "family" (wire format).

 Implementations isolate everything that differs between the OpenAI Chat Completions API, the
 newer OpenAI Responses API, and the native Anthropic Messages API: how a conversation's turns are
 serialized, how the request body is assembled (system-prompt placement, token caps, tool schema,
 tool_choice, reasoning/thinking), how a raw response is parsed, and which HTTP headers carry auth.

 Instances are stateless and cheap; create one per `LlmInterface` from the active provider's
 resolved `apiFormat()`.
 */
class LlmApiProtocol
{
public:
  virtual ~LlmApiProtocol() = default;

  /** Serialize the turns of the given conversations into this format's message/input array.

   The returned array does NOT include the system prompt (that is placed by buildRequestBody()).
   For MainAgent requests `convos` is the prior history through the current conversation; for
   sub-agents it is just the single conversation.
   */
  virtual nlohmann::json serializeConversations(
      const std::vector<std::shared_ptr<LlmInteraction>> &convos ) const = 0;

  /** Append an ephemeral (not-persisted-in-history) user-authored note to the messages array.

   Used for the state-machine reminder.  Placement differs per format (e.g. OpenAI appends to a
   trailing role:"tool" message; Anthropic appends a text block to the trailing tool_result user
   message; otherwise a new user message is pushed).
   */
  virtual void appendEphemeralUserNote( nlohmann::json &messages, const std::string &note ) const = 0;

  /** Assemble the full request body: model + token caps + reasoning + system placement + the
   already-serialized messages + tools + tool_choice.
   */
  virtual nlohmann::json buildRequestBody(
      const LlmConfig::LlmApi &api,
      const std::string &systemPrompt,
      const nlohmann::json &messages,
      const std::vector<NormalizedTool> &tools,
      ToolChoice toolChoice ) const = 0;

  /** Parse a raw (already-JSON-parsed) provider response into the normalized struct.
   May throw on a structurally unrecognizable payload (handled by the caller).
   */
  virtual ParsedLlmResponse parseResponse( const nlohmann::json &responseJson ) const = 0;

  /** HTTP headers (other than the JSON content-type, which the bridge always sets) for the request,
   as ordered key/value pairs.  `endpoint` is provided so format-specific endpoint quirks can be
   honored (e.g. the Ask Sage `x-access-tokens` header for OpenAI-compatible endpoints).
   */
  virtual std::vector<std::pair<std::string,std::string>> headers(
      const std::string &endpoint, const std::string &bearerToken ) const = 0;

  /** Field names, in cache-stable priority order, to emit *before* the bulk message/input array
   when serializing the request for prompt-cache reuse (see serializeRequestForCaching).
   */
  virtual std::vector<std::string> cachePriorityFields() const = 0;

  /** Factory: create the protocol implementation for a given wire format. */
  static std::unique_ptr<const LlmApiProtocol> create( LlmConfig::LlmApi::ApiFormat fmt );
};//class LlmApiProtocol

#endif // LLM_API_PROTOCOL_H
