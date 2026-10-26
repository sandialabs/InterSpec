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
#include <vector>
#include <variant>
#include <algorithm>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/LlmApiProtocol.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace std;
using json = nlohmann::json;

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

namespace
{
/** Anthropic requires `max_tokens`; used when the active model leaves MaxTokens unset (0). */
const int sm_anthropic_default_max_tokens = 4096;

/** The family ("opus", "sonnet", "haiku", "fable", ...) and version of a Claude model id. */
struct ClaudeModelVersion
{
  std::string family;
  int version = 0;  // 100*major + minor (e.g., 405 for 4.5); 0 if not recognized
};


/** Parses ids like "claude-sonnet-4-5", "claude-opus-4-20250514" (= 4.0), OpenRouter's
 "anthropic/claude-sonnet-4.5", Bedrock's "anthropic.claude-haiku-4-5-20251001-v1:0", Vertex's
 "claude-opus-4-5@20251101", and the Claude 3 style "claude-3-5-haiku" (version before family, which
 Ask Sage also uses, e.g. "claude-35-sonnet").
 */
ClaudeModelVersion claude_model_version( std::string model )
{
  ClaudeModelVersion result;
  SpecUtils::to_lower_ascii( model );
  const size_t pos = model.find( "claude-" );
  if( pos == std::string::npos )
    return result;

  std::string rest = model.substr( pos + 7 );
  std::replace( begin(rest), end(rest), '.', '-' );
  std::replace( begin(rest), end(rest), '@', '-' );
  std::vector<std::string> tokens;
  SpecUtils::split( tokens, rest, "-" );

  const auto is_number = []( const std::string &s ) -> bool {
    return !s.empty() && (s.size() <= 2)
           && std::all_of( begin(s), end(s), []( char c ){ return (c >= '0') && (c <= '9'); } );
  };

  // Version is either first ("3-5-haiku"), or follows the family ("sonnet-4-5").
  size_t ver_index = 0;
  if( !tokens.empty() && is_number(tokens[0]) )
  {
    if( tokens.size() > 1 )
      result.family = (is_number(tokens[1]) && (tokens.size() > 2)) ? tokens[2] : tokens[1];
  }else if( tokens.size() > 1 )
  {
    result.family = tokens[0];
    ver_index = 1;
  }

  if( (ver_index == 0) && !tokens.empty() && (tokens[0].size() == 2) && is_number(tokens[0]) )
  {
    // Ask Sage style "claude-35-sonnet"
    result.version = 100*(tokens[0][0] - '0') + (tokens[0][1] - '0');
  }else if( (ver_index < tokens.size()) && is_number(tokens[ver_index]) )
  {
    const int major = std::stoi( tokens[ver_index] );
    // A minor version is one or two digits; a date ("4-20250514") or nothing means ".0".
    const bool has_minor = ((ver_index + 1) < tokens.size()) && is_number(tokens[ver_index + 1]);
    const int minor = has_minor ? std::stoi( tokens[ver_index + 1] ) : 0;
    result.version = 100*major + minor;
  }

  return result;
}//claude_model_version(...)


/** Whether a Claude model accepts thinking:{type:"adaptive"} (the model sizes its own thinking).
 Adaptive thinking arrived with Opus/Sonnet 4.6; Haiku 4.5 and all earlier models reject it with a
 400, and need the explicit {type:"enabled",budget_tokens:N} form instead - which is in turn
 rejected by Opus 4.7+, Sonnet 5+, Haiku 5.5, and Fable.  Unrecognized ids get the modern form. */
bool claude_supports_adaptive_thinking( const ClaudeModelVersion &claude )
{
  if( claude.version == 0 )
    return true;
  if( claude.family == "haiku" )
    return (claude.version >= 500);
  return (claude.version >= 406);
}


/** Whether a Claude model rejects a `temperature` (and other sampling parameters): Opus 4.7+,
 Sonnet 5+, Haiku 5+, and Fable / Mythos. */
bool claude_rejects_temperature( const ClaudeModelVersion &claude )
{
  if( (claude.family == "fable") || (claude.family == "mythos") )
    return true;
  if( claude.family == "opus" )
    return (claude.version >= 407);
  return (claude.version >= 500);
}


/** Whether a Claude model rejects a forced tool_choice ("any" or "tool"): Opus/Sonnet 5.5+, and
 Fable / Mythos 5.1+. */
bool claude_rejects_forced_tool_choice( const ClaudeModelVersion &claude )
{
  if( (claude.family == "fable") || (claude.family == "mythos") )
    return (claude.version >= 501);
  if( (claude.family == "opus") || (claude.family == "sonnet") )
    return (claude.version >= 505);
  return false;
}


/** Extract a non-negative integer token count from a usage field, if present and numeric.
 Read as 64-bit (or via double for float-encoded counts) so large counts do not overflow a 32-bit int. */
std::optional<size_t> get_usage_count( const json &usage, const char *key )
{
  if( !usage.contains(key) || !usage[key].is_number() )
    return std::nullopt;

  const json &v = usage[key];
  const long long count = v.is_number_float() ? static_cast<long long>( v.get<double>() )
                                              : v.get<long long>();
  if( count < 0 )
    return std::nullopt;
  return static_cast<size_t>( count );
}//get_usage_count(...)


/** Attach a 5-minute ephemeral cache breakpoint to the last content block of `msg`, promoting
 string content to a single text block first so it can carry the marker.  Only user turns are
 marked: that keeps the breakpoint off thinking blocks (which the API disallows cache_control on)
 and matches where a request normally ends before an assistant reply.
 Returns whether a breakpoint was set.
 */
bool set_anthropic_cache_breakpoint( json &msg )
{
  if( !msg.is_object() || (msg.value("role","") != "user") || !msg.contains("content") )
    return false;

  json &content = msg["content"];
  if( content.is_string() )
  {
    // Promote string content to a single text block so it can carry the cache breakpoint.
    json textBlock;
    textBlock["type"] = "text";
    textBlock["text"] = content.get<string>();
    textBlock["cache_control"]["type"] = "ephemeral";
    content = json::array({ textBlock });
    return true;
  }

  if( content.is_array() && !content.empty() && content.back().is_object() )
  {
    content.back()["cache_control"]["type"] = "ephemeral";
    return true;
  }

  return false;
}//set_anthropic_cache_breakpoint(...)


/** Place the rolling conversation cache breakpoints, so the growing prefix is cached incrementally
 across tool-call rounds and follow-up turns.  Two markers are placed:

   - the "head", on the end of the newest user turn.  This writes the entry that the *next* request
     reads back.
   - the "anchor", on the end of the previous user turn.  This is exactly where the previous
     request's head marker wrote its entry, so reading it is an exact-position hit.  The anchor
     matters because a breakpoint only searches back a bounded number of content blocks for a
     usable entry, and a single round of parallel tool calls can add enough blocks to overshoot.

 IMPORTANT: call this *before* appending any text that is not persisted into history (i.e. the
 state-machine reminder from appendEphemeralUserNote()).  A breakpoint written onto a block the
 next request will not contain can never be read back - it only ever pays the cache-write premium.
 */
void add_anthropic_cache_breakpoints( json &messages )
{
  if( !messages.is_array() )
    return;

  vector<size_t> userMsgIdx;
  for( size_t i = 0; i < messages.size(); ++i )
  {
    if( !messages[i].is_object() || (messages[i].value("role","") != "user")
        || !messages[i].contains("content") )
      continue;  // note: contains() guard matters - operator[] below would insert a null "content"

    // Normalize *every* user turn's string content to a one-block array, not just the two we are
    // about to mark.  Marking is what promotes string content to a block array, so marking only
    // some turns would make a given message's wire shape oscillate between rounds (block array
    // while it is the head/anchor, bare string once it is neither).  Whether the provider treats
    // the two shapes as the same cached content is an assumption we do not need to make, and if it
    // did not the conversation would silently lose its cache with no error - only cost.
    json &content = messages[i]["content"];
    if( content.is_string() )
    {
      json textBlock;
      textBlock["type"] = "text";
      textBlock["text"] = content.get<string>();
      content = json::array({ textBlock });
    }

    userMsgIdx.push_back( i );
  }//for( i : messages )

  if( userMsgIdx.empty() )
    return;

  set_anthropic_cache_breakpoint( messages[userMsgIdx.back()] );  // head

  if( userMsgIdx.size() >= 2 )
    set_anthropic_cache_breakpoint( messages[userMsgIdx[userMsgIdx.size() - 2]] );  // anchor
}//add_anthropic_cache_breakpoints(...)


/** Normalize string `content` into an array of text blocks (block `type` given by `textType`), so two
 messages with the same role can be concatenated.  Array content passes through; null/empty becomes [].
 */
json content_to_block_array( const json &content, const char *textType )
{
  if( content.is_array() )
    return content;
  json arr = json::array();
  if( content.is_string() && !content.get<std::string>().empty() )
  {
    json b;
    b["type"] = textType;
    b["text"] = content.get<std::string>();
    arr.push_back( b );
  }
  return arr;
}//content_to_block_array(...)


/** Merge consecutive messages that share the same role.  The Anthropic Messages API rejects
 consecutive same-role turns (and the Responses API mishandles them); this can otherwise happen at
 the summarization seam (summary `user` immediately followed by the next conversation's `user`
 InitialRequest).  Items without a "role" (e.g. Responses function_call / function_call_output) are
 never merged.  `textType` is the per-format text-block type ("text" for Anthropic, "input_text" for
 the Responses input).
 */
void merge_adjacent_same_role( json &messages, const char *textType )
{
  if( !messages.is_array() )
    return;

  json out = json::array();
  for( const json &msg : messages )
  {
    const std::string role = msg.is_object() ? msg.value("role","") : std::string();
    if( !role.empty() && !out.empty() && out.back().is_object()
        && (out.back().value("role","") == role) )
    {
      json &prev = out.back()["content"];
      const json cur = msg.contains("content") ? msg["content"] : json();
      if( prev.is_string() && cur.is_string() )
      {
        prev = prev.get<std::string>() + "\n\n" + cur.get<std::string>();
      }else
      {
        json merged = content_to_block_array( prev, textType );
        for( json &b : content_to_block_array( cur, textType ) )
          merged.push_back( b );
        prev = merged;
      }
    }else
    {
      out.push_back( msg );
    }
  }//for( msg : messages )

  messages = out;
}//merge_adjacent_same_role(...)


//==============================================================================
//  OpenAI Chat Completions  (the legacy format, and all OpenAI-compatible providers)
//==============================================================================
class OpenAiChatProtocol : public LlmApiProtocol
{
public:
  nlohmann::json serializeConversations(
      const std::vector<std::shared_ptr<LlmInteraction>> &convos ) const override
  {
    // Delegate to the existing, well-exercised OpenAI-Chat serializer so this format stays
    // byte-for-byte identical to the pre-refactor behavior.
    json messages = json::array();
    for( const std::shared_ptr<LlmInteraction> &c : convos )
    {
      if( c )
        LlmConversationHistory::addConversationToLlmApiHistory( *c, messages );
    }
    return messages;
  }


  void appendEphemeralUserNote( nlohmann::json &messages, const std::string &note ) const override
  {
    // Some providers (e.g. Mistral) require an assistant turn immediately after tool-result
    // messages (role "tool"), so we append the note to the last tool-result message instead of
    // inserting a new user message there.
    if( !messages.empty() && (messages.back().value("role","") == "tool") )
    {
      json &last = messages.back();
      if( last["content"].is_string() )
      {
        last["content"] = last["content"].get<string>() + "\n\n" + note;
      }else if( last["content"].is_array() )
      {
        json textBlock;
        textBlock["type"] = "text";
        textBlock["text"] = note;
        last["content"].push_back( textBlock );
      }
    }else
    {
      json msg;
      msg["role"] = "user";
      msg["content"] = note;
      messages.push_back( msg );
    }
  }//appendEphemeralUserNote(...)


  /** No-op: the OpenAI Chat Completions API caches automatically by longest-common-prefix, with no
   client-side markers, so a trailing ephemeral note simply falls outside the matched prefix.
   */
  void placeConversationCacheBreakpoints( nlohmann::json &/*messages*/ ) const override
  {
  }


  nlohmann::json buildRequestBody( const LlmConfig::LlmApi &api,
                                   const std::string &systemPrompt,
                                   const nlohmann::json &messages,
                                   const std::vector<NormalizedTool> &tools,
                                   ToolChoice toolChoice ) const override
  {
    json request;
    request["model"] = api.model();

    // Use max_completion_tokens for newer OpenAI models, max_tokens for others.
    if( api.maxTokens() > 0 )
    {
      const string modelName = api.model();
      if( (modelName.find("gpt-") != string::npos)
         || (modelName.find("o1") != string::npos)
         || (modelName.find("o3") != string::npos)
         || (modelName.find("o4") != string::npos)
         || (modelName.find("gemini") != string::npos) )
      {
        request["max_completion_tokens"] = api.maxTokens();
      }else
      {
        request["max_tokens"] = api.maxTokens();
      }
    }//if( api.maxTokens() > 0 )

    // Reasoning configuration based on the variant type.
    const std::variant<bool,LlmConfig::LlmApi::ReasoningEffort> &reasoning = api.reasoning();
    if( std::holds_alternative<bool>(reasoning) )
    {
      if( std::get<bool>(reasoning) )
        request["reasoning"]["enabled"] = true;  // OpenRouter style
    }else if( std::holds_alternative<LlmConfig::LlmApi::ReasoningEffort>(reasoning) )
    {
      switch( std::get<LlmConfig::LlmApi::ReasoningEffort>(reasoning) )
      {
        case LlmConfig::LlmApi::ReasoningEffort::low:    request["reasoning_effort"] = "low";    break;
        case LlmConfig::LlmApi::ReasoningEffort::medium: request["reasoning_effort"] = "medium"; break;
        case LlmConfig::LlmApi::ReasoningEffort::high:   request["reasoning_effort"] = "high";   break;
      }
    }

    if( api.temperature().has_value() )
      request["temperature"] = api.temperature().value();

    // The OpenAI Chat format carries the system prompt as the first message.
    json msgs = messages;
    if( !systemPrompt.empty() )
    {
      json systemMsg;
      systemMsg["role"] = "system";
      systemMsg["content"] = systemPrompt;
      msgs.insert( msgs.begin(), systemMsg );
    }
    request["messages"] = msgs;

    if( !tools.empty() )
    {
      json toolsArr = json::array();
      for( const NormalizedTool &t : tools )
      {
        json toolDef;
        toolDef["type"] = "function";
        toolDef["function"]["name"] = t.name;
        toolDef["function"]["description"] = t.description;
        toolDef["function"]["parameters"] = t.parameters;
        toolsArr.push_back( toolDef );
      }
      request["tools"] = toolsArr;
      request["tool_choice"] = (toolChoice == ToolChoice::Required) ? "required" : "auto";
    }

    return request;
  }//buildRequestBody(...)


  ParsedLlmResponse parseResponse( const nlohmann::json &responseJson ) const override
  {
    ParsedLlmResponse out;

    if( responseJson.contains("usage") && responseJson["usage"].is_object() )
    {
      const json &usage = responseJson["usage"];
      out.promptTokens = get_usage_count( usage, "prompt_tokens" );
      out.completionTokens = get_usage_count( usage, "completion_tokens" );
      out.totalTokens = get_usage_count( usage, "total_tokens" );

      if( usage.contains("prompt_tokens_details") && usage["prompt_tokens_details"].is_object() )
        out.cachedTokens = get_usage_count( usage["prompt_tokens_details"], "cached_tokens" );
    }

    auto parseOpenAiToolCalls = [&out]( const json &toolCalls )
    {
      for( const json &tc : toolCalls )
      {
        ParsedLlmResponse::ToolCall call;
        call.id = tc.value("id","");
        if( tc.contains("function") && tc["function"].is_object() )
        {
          const json &fn = tc["function"];
          call.name = fn.value("name","");
          if( fn.contains("arguments") && fn["arguments"].is_string() )
            call.argsRaw = fn["arguments"].get<string>();
        }
        if( tc.contains("extra_content") )
          call.extra = tc["extra_content"];
        out.toolCalls.push_back( std::move(call) );
      }
    };//parseOpenAiToolCalls

    if( responseJson.contains("choices") && responseJson["choices"].is_array()
       && !responseJson["choices"].empty() )
    {
      const json &choice = responseJson["choices"][0];
      if( choice.contains("message") && choice["message"].is_object() )
      {
        const json &message = choice["message"];

        // content: string, or a Claude-style content-block array (thinking + text)
        if( message.contains("content") )
        {
          if( message["content"].is_string() )
          {
            out.content = message["content"].get<string>();
          }else if( message["content"].is_array() )
          {
            for( const json &block : message["content"] )
            {
              const string blockType = block.value("type","");
              if( blockType == "thinking" )
              {
                out.thinkingContent = block.value("thinking","");
                out.thinkingSignature = block.value("signature","");
              }else if( blockType == "text" )
              {
                out.content += block.value("text","");
              }
            }
          }
        }//if( message.contains("content") )

        if( message.contains("reasoning_content") && message["reasoning_content"].is_string() )
        {
          out.reasoningContent = message["reasoning_content"].get<string>();
          if( out.thinkingContent.empty() )
            out.thinkingContent = out.reasoningContent;
        }

        if( message.contains("reasoning") && message["reasoning"].is_string() )
        {
          const string reasoning = message["reasoning"].get<string>();
          if( out.thinkingContent.empty() )
            out.thinkingContent = reasoning;
        }

        // OpenRouter reasoning_details: store the array verbatim (must be passed back unmodified),
        // and use its text only for display when we have nothing else.
        if( message.contains("reasoning_details") && message["reasoning_details"].is_array() )
        {
          out.reasoningDetails = message["reasoning_details"].dump();
          if( out.thinkingContent.empty() )
          {
            for( const json &detail : message["reasoning_details"] )
            {
              if( !detail.is_object() )
                continue;
              const string detailType = detail.value("type","");
              if( (detailType == "reasoning.text") && detail.contains("text") && detail["text"].is_string() )
              {
                const string text = detail["text"].get<string>();
                if( !text.empty() )
                  out.thinkingContent += (out.thinkingContent.empty() ? "" : "\n") + text;
              }else if( (detailType == "reasoning.summary") && detail.contains("summary") && detail["summary"].is_string() )
              {
                const string summary = detail["summary"].get<string>();
                if( !summary.empty() )
                  out.thinkingContent += (out.thinkingContent.empty() ? "" : "\n") + summary;
              }
            }
          }
        }//if( reasoning_details )

        if( message.contains("tool_calls") && message["tool_calls"].is_array() )
          parseOpenAiToolCalls( message["tool_calls"] );
      }//if( choice.contains("message") )
    }//if( has choices )

    // Ask Sage `/server/query` style: top-level tool_calls and message fields.
    if( out.toolCalls.empty() && responseJson.contains("tool_calls")
       && responseJson["tool_calls"].is_array() && !responseJson["tool_calls"].empty() )
    {
      parseOpenAiToolCalls( responseJson["tool_calls"] );
    }

    if( out.content.empty() && responseJson.contains("message")
       && responseJson["message"].is_string() )
    {
      out.content = responseJson["message"].get<string>();
    }

    return out;
  }//parseResponse(...)


  std::vector<std::pair<std::string,std::string>> headers(
      const std::string &endpoint, const std::string &bearerToken ) const override
  {
    std::vector<std::pair<std::string,std::string>> hdrs;
    if( !bearerToken.empty() )
    {
      hdrs.emplace_back( "Authorization", "Bearer " + bearerToken );

      // Ask Sage's `/server/query` endpoint additionally expects the token in x-access-tokens.
      // Other providers reject the extra header, so only set it for the endpoint that needs it.
      if( SpecUtils::icontains( endpoint, ".ai/server/query" ) )
        hdrs.emplace_back( "x-access-tokens", bearerToken );
    }
    return hdrs;
  }//headers(...)


  std::vector<std::string> cachePriorityFields() const override
  {
    return { "model", "temperature", "max_completion_tokens", "max_tokens",
             "reasoning", "reasoning_effort", "tools", "tool_choice" };
  }
};//class OpenAiChatProtocol


//==============================================================================
//  Anthropic Messages API
//==============================================================================
class AnthropicProtocol : public LlmApiProtocol
{
  static void addConversation( const LlmInteraction &conv, json &messages )
  {
    for( const std::shared_ptr<LlmInteractionTurn> &response : conv.responses )
    {
      if( response->excludeFromHistory() )
        continue;

      switch( response->type() )
      {
        case LlmInteractionTurn::Type::InitialRequest:
        {
          const LlmInteractionInitialRequest *req
            = dynamic_cast<const LlmInteractionInitialRequest*>( response.get() );
          if( !req )
            break;
          // Anthropic has no system role in `messages` (the agent system prompt is top-level), so
          // both user and (app-generated) system queries become user turns here.
          json msg;
          msg["role"] = "user";
          if( !req->imageContent().empty() )
          {
            json arr = json::array();
            json textBlock;
            textBlock["type"] = "text";
            textBlock["text"] = req->content();
            arr.push_back( textBlock );
            for( const LlmToolCall::ImageContent &img : req->imageContent() )
            {
              // Only benchmark questions set a caption, so pasted/uploaded images produce the
              //  exact same payload as before.
              if( !img.caption.empty() )
              {
                json capBlock;
                capBlock["type"] = "text";
                capBlock["text"] = img.caption;
                arr.push_back( capBlock );
              }

              json imgBlock;
              imgBlock["type"] = "image";
              imgBlock["source"] = { {"type","base64"}, {"media_type", img.mimeType}, {"data", img.base64Data} };
              arr.push_back( imgBlock );
            }
            msg["content"] = arr;
          }else
          {
            msg["content"] = req->content();
          }
          messages.push_back( msg );
          break;
        }//InitialRequest

        case LlmInteractionTurn::Type::FinalLlmResponse:
        {
          const LlmInteractionFinalResponse *resp
            = dynamic_cast<const LlmInteractionFinalResponse*>( response.get() );
          if( !resp )
            break;
          json msg;
          msg["role"] = "assistant";
          const string &sig = resp->thinkingSignature();
          const string &think = resp->thinkingContent();
          if( !sig.empty() && !think.empty() )
          {
            json arr = json::array();
            json thinkingBlock;
            thinkingBlock["type"] = "thinking";
            thinkingBlock["thinking"] = think;
            thinkingBlock["signature"] = sig;
            arr.push_back( thinkingBlock );
            if( !resp->content().empty() )
            {
              json textBlock;
              textBlock["type"] = "text";
              textBlock["text"] = resp->content();
              arr.push_back( textBlock );
            }
            msg["content"] = arr;
          }else
          {
            // Skip a FinalLlmResponse with no text and no (signed) thinking: Anthropic rejects an
            // assistant message with empty content (can happen for an empty model turn, or a
            // signatureless/redacted thinking response).
            if( resp->content().empty() )
              break;
            msg["content"] = resp->content();
          }
          messages.push_back( msg );
          break;
        }//FinalLlmResponse

        case LlmInteractionTurn::Type::ToolCall:
        {
          const LlmToolRequest *req = dynamic_cast<const LlmToolRequest*>( response.get() );
          if( !req )
            break;
          json msg;
          msg["role"] = "assistant";
          json arr = json::array();
          const string &sig = req->thinkingSignature();
          const string &think = req->thinkingContent();
          if( !sig.empty() && !think.empty() )
          {
            json thinkingBlock;
            thinkingBlock["type"] = "thinking";
            thinkingBlock["thinking"] = think;
            thinkingBlock["signature"] = sig;
            arr.push_back( thinkingBlock );
          }
          for( const LlmToolCall &tc : req->toolCalls() )
          {
            json toolUse;
            toolUse["type"] = "tool_use";
            toolUse["id"] = tc.invocationId;
            toolUse["name"] = tc.toolName;
            toolUse["input"] = tc.toolParameters.is_null() ? json::object() : tc.toolParameters;
            arr.push_back( toolUse );
          }
          msg["content"] = arr;
          messages.push_back( msg );
          break;
        }//ToolCall

        case LlmInteractionTurn::Type::ToolResult:
        {
          const LlmToolResults *res = dynamic_cast<const LlmToolResults*>( response.get() );
          if( !res )
            break;
          // All tool results for the preceding assistant tool_use blocks go in ONE user message.
          json arr = json::array();
          for( const LlmToolCall &tc : res->toolCalls() )
          {
            json toolResult;
            toolResult["type"] = "tool_result";
            toolResult["tool_use_id"] = tc.invocationId;
            if( tc.imageContent.has_value() )
            {
              json contentArr = json::array();
              if( !tc.content.empty() )
              {
                json textBlock;
                textBlock["type"] = "text";
                textBlock["text"] = tc.content;
                contentArr.push_back( textBlock );
              }
              json imgBlock;
              imgBlock["type"] = "image";
              imgBlock["source"] = { {"type","base64"}, {"media_type", tc.imageContent->mimeType},
                                     {"data", tc.imageContent->base64Data} };
              contentArr.push_back( imgBlock );
              toolResult["content"] = contentArr;
            }else
            {
              // Anthropic rejects an empty tool_result content - substitute a placeholder.
              toolResult["content"] = tc.content.empty() ? string("(no output)") : tc.content;
            }
            if( tc.status == LlmToolCall::CallStatus::Error )
              toolResult["is_error"] = true;
            arr.push_back( toolResult );
          }
          json msg;
          msg["role"] = "user";
          msg["content"] = arr;
          messages.push_back( msg );
          break;
        }//ToolResult

        case LlmInteractionTurn::Type::Error:
        {
          const LlmInteractionError *err = dynamic_cast<const LlmInteractionError*>( response.get() );
          if( !err )
            break;
          json msg;
          msg["role"] = "assistant";
          msg["content"] = "Error: " + err->errorMessage();
          messages.push_back( msg );
          break;
        }//Error

        case LlmInteractionTurn::Type::AutoReply:
        {
          const LlmInteractionAutoReply *ar = dynamic_cast<const LlmInteractionAutoReply*>( response.get() );
          if( !ar )
            break;
          json msg;
          msg["role"] = "user";
          msg["content"] = ar->content();
          messages.push_back( msg );
          break;
        }//AutoReply

        case LlmInteractionTurn::Type::ConversationSummary:
        {
          const LlmConversationSummary *sum = dynamic_cast<const LlmConversationSummary*>( response.get() );
          if( !sum )
            break;
          json msg;
          msg["role"] = "user";
          msg["content"] = "[System: Summary of previous conversation history]\n\n" + sum->summary();
          messages.push_back( msg );
          break;
        }//ConversationSummary
      }//switch( response->type() )
    }//for( responses )
  }//addConversation(...)

public:
  nlohmann::json serializeConversations(
      const std::vector<std::shared_ptr<LlmInteraction>> &convos ) const override
  {
    json messages = json::array();
    for( const std::shared_ptr<LlmInteraction> &c : convos )
    {
      if( c )
        addConversation( *c, messages );
    }
    // Anthropic rejects consecutive same-role turns (e.g. summary user -> InitialRequest user).
    merge_adjacent_same_role( messages, "text" );
    return messages;
  }


  void appendEphemeralUserNote( nlohmann::json &messages, const std::string &note ) const override
  {
    // NOTE: placeConversationCacheBreakpoints() must already have run (LlmInterface calls it right
    // after serializeConversations()).  This note is deliberately not persisted into history, so it
    // must land *after* the last breakpoint: a cached prefix ending in text the next request cannot
    // reproduce can never be read back, and only ever pays the cache-write premium.

    // Anthropic has no role:"tool"; tool results, auto-replies, and user queries are all user
    // messages.  Merge the note into a trailing user turn (rather than pushing a second user
    // message) so the note stays attached to that turn and we never emit consecutive user messages.
    if( !messages.empty() && (messages.back().value("role","") == "user")
       && messages.back().contains("content") )
    {
      json &content = messages.back()["content"];
      json textBlock;
      textBlock["type"] = "text";
      textBlock["text"] = note;

      if( content.is_array() )
      {
        content.push_back( textBlock );
      }else if( content.is_string() )
      {
        // Promote the string content to a block array so we can carry both pieces.
        json existing;
        existing["type"] = "text";
        existing["text"] = content.get<string>();
        content = json::array({ existing, textBlock });
      }else
      {
        content = json::array({ textBlock });
      }
    }else
    {
      json msg;
      msg["role"] = "user";
      msg["content"] = note;
      messages.push_back( msg );
    }
  }//appendEphemeralUserNote(...)


  void placeConversationCacheBreakpoints( nlohmann::json &messages ) const override
  {
    add_anthropic_cache_breakpoints( messages );
  }


  nlohmann::json buildRequestBody( const LlmConfig::LlmApi &api,
                                   const std::string &systemPrompt,
                                   const nlohmann::json &messages,
                                   const std::vector<NormalizedTool> &tools,
                                   ToolChoice toolChoice ) const override
  {
    json request;
    request["model"] = api.model();

    // max_tokens is REQUIRED by Anthropic.
    request["max_tokens"] = (api.maxTokens() > 0) ? api.maxTokens() : sm_anthropic_default_max_tokens;

    if( !systemPrompt.empty() )
    {
      // Use the block-array form of `system` so it can carry a 5-minute ephemeral cache breakpoint.
      // Tools render before system, so this single breakpoint caches the large, stable tools+system
      // prefix together (re-read cheaply on every subsequent request in the conversation).
      json sysBlock;
      sysBlock["type"] = "text";
      sysBlock["text"] = systemPrompt;
      sysBlock["cache_control"]["type"] = "ephemeral";
      request["system"] = json::array({ sysBlock });
    }

    // Conversation breakpoints are NOT placed here: only the caller knows whether this request is
    // part of a continuing conversation, and where the not-persisted ephemeral note ends.  See
    // placeConversationCacheBreakpoints().
    request["messages"] = messages;

    // Map any configured reasoning to Claude extended thinking.  Claude 4.6+ models take the
    // "adaptive" type (the model sizes its own thinking), with the effort level in
    // output_config.effort; older ones (and Haiku 4.5) need the explicit
    // {type:"enabled",budget_tokens:N} form - gate on model capability so reasoning works across
    // models.  When thinking is on, `temperature` must be omitted, and the newest models reject it
    // outright, so only set temperature otherwise.
    const ClaudeModelVersion claude = claude_model_version( api.model() );
    const bool allowTemperature = api.temperature().has_value() && !claude_rejects_temperature( claude );
    const std::variant<bool,LlmConfig::LlmApi::ReasoningEffort> &reasoning = api.reasoning();
    bool thinkingEnabled = false;
    LlmConfig::LlmApi::ReasoningEffort effort = LlmConfig::LlmApi::ReasoningEffort::medium;
    if( std::holds_alternative<bool>(reasoning) )
    {
      thinkingEnabled = std::get<bool>(reasoning);
    }else
    {
      thinkingEnabled = true;  // any explicit effort level enables thinking
      effort = std::get<LlmConfig::LlmApi::ReasoningEffort>(reasoning);
    }

    if( thinkingEnabled && claude_supports_adaptive_thinking( claude ) )
    {
      request["thinking"]["type"] = "adaptive";

      // Plain "true" leaves the model's default effort.
      if( std::holds_alternative<LlmConfig::LlmApi::ReasoningEffort>(reasoning) )
      {
        switch( effort )
        {
          case LlmConfig::LlmApi::ReasoningEffort::low:    request["output_config"]["effort"] = "low";    break;
          case LlmConfig::LlmApi::ReasoningEffort::medium: request["output_config"]["effort"] = "medium"; break;
          case LlmConfig::LlmApi::ReasoningEffort::high:   request["output_config"]["effort"] = "high";   break;
        }
      }
    }else if( thinkingEnabled )
    {
      // Explicit extended-thinking budget (supported across thinking-capable models incl. Haiku).
      // Anthropic requires 1024 <= budget_tokens < max_tokens; scale by effort and clamp to fit.
      const int maxTok = request["max_tokens"].get<int>();
      int budget = 4096;  // medium / unspecified default
      if( effort == LlmConfig::LlmApi::ReasoningEffort::low )
        budget = 2048;
      else if( effort == LlmConfig::LlmApi::ReasoningEffort::high )
        budget = 8192;
      if( budget >= maxTok )
        budget = maxTok - 1;
      if( budget < 1024 )
        budget = 1024;

      if( budget < maxTok )  // only enable if a valid budget fits under max_tokens
      {
        request["thinking"]["type"] = "enabled";
        request["thinking"]["budget_tokens"] = budget;
      }else if( allowTemperature )
      {
        request["temperature"] = api.temperature().value();  // max_tokens too small for thinking
      }
    }else if( allowTemperature )
    {
      request["temperature"] = api.temperature().value();
    }

    if( !tools.empty() )
    {
      json toolsArr = json::array();
      for( const NormalizedTool &t : tools )
      {
        json toolDef;
        toolDef["name"] = t.name;
        toolDef["description"] = t.description;
        // Anthropic rejects oneOf/anyOf/allOf at the top level of input_schema.
        toolDef["input_schema"] = strip_top_level_schema_combinators( t.parameters );
        toolsArr.push_back( toolDef );
      }
      request["tools"] = toolsArr;
      // Anthropic requires tool_choice to be "auto" (or "none") when extended thinking is enabled,
      // and the newest models reject a forced choice outright, so don't force "any" in those cases
      // (it would be rejected with a 400).
      const bool forceTool = (toolChoice == ToolChoice::Required) && !thinkingEnabled
                             && !claude_rejects_forced_tool_choice( claude );
      request["tool_choice"]["type"] = forceTool ? "any" : "auto";
    }

    return request;
  }//buildRequestBody(...)


  ParsedLlmResponse parseResponse( const nlohmann::json &responseJson ) const override
  {
    ParsedLlmResponse out;

    if( responseJson.contains("usage") && responseJson["usage"].is_object() )
    {
      const json &usage = responseJson["usage"];
      out.completionTokens = get_usage_count( usage, "output_tokens" );
      out.cachedTokens = get_usage_count( usage, "cache_read_input_tokens" );
      out.cacheCreationTokens = get_usage_count( usage, "cache_creation_input_tokens" );

      // Anthropic's `input_tokens` counts only the *uncached* input; cache reads and cache writes
      // are reported separately and are additive.  Fold them in so promptTokens means the TOTAL
      // input processed - matching OpenAI's `prompt_tokens` (which already includes cached tokens).
      // This keeps the cached/input ratio sensible (<=100%) and the input count comparable across
      // providers, and makes cachedTokens a true subset of promptTokens.
      out.promptTokens = get_usage_count( usage, "input_tokens" );
      if( out.promptTokens.has_value() || out.cachedTokens.has_value() || out.cacheCreationTokens.has_value() )
        out.promptTokens = out.promptTokens.value_or(0) + out.cachedTokens.value_or(0) + out.cacheCreationTokens.value_or(0);

      if( out.promptTokens.has_value() && out.completionTokens.has_value() )
        out.totalTokens = out.promptTokens.value() + out.completionTokens.value();
    }

    out.stopReason = responseJson.value("stop_reason","");

    if( responseJson.contains("content") && responseJson["content"].is_array() )
    {
      for( const json &block : responseJson["content"] )
      {
        const string blockType = block.value("type","");
        if( blockType == "text" )
        {
          out.content += block.value("text","");
        }else if( blockType == "thinking" )
        {
          out.thinkingContent = block.value("thinking","");
          out.thinkingSignature = block.value("signature","");
        }else if( blockType == "tool_use" )
        {
          ParsedLlmResponse::ToolCall call;
          call.id = block.value("id","");
          call.name = block.value("name","");
          call.argsRaw = block.contains("input") ? block["input"].dump() : string("{}");
          out.toolCalls.push_back( std::move(call) );
        }
      }
    }//if( content array )

    return out;
  }//parseResponse(...)


  std::vector<std::pair<std::string,std::string>> headers(
      const std::string &endpoint, const std::string &bearerToken ) const override
  {
    // The native Anthropic API (api.anthropic.com) authenticates with x-api-key, requires a version
    // header, and needs the "dangerous-direct-browser-access" opt-in so the browser's fetch() is
    // permitted cross-origin.  Those three headers are specific to api.anthropic.com.
    if( SpecUtils::icontains( endpoint, "api.anthropic.com" ) )
    {
      return {
        { "x-api-key", bearerToken },
        { "anthropic-version", "2023-06-01" },
        { "anthropic-dangerous-direct-browser-access", "true" }
      };
    }

    // Other providers that merely speak the Anthropic Messages *format* (e.g. OpenRouter's
    // /v1/messages) authenticate with a Bearer token, and their CORS policy rejects the
    // Anthropic-specific headers above in preflight - so send only the Bearer token here.
    std::vector<std::pair<std::string,std::string>> hdrs;
    if( !bearerToken.empty() )
      hdrs.emplace_back( "Authorization", "Bearer " + bearerToken );
    return hdrs;
  }//headers(...)


  std::vector<std::string> cachePriorityFields() const override
  {
    return { "model", "max_tokens", "system", "thinking", "output_config", "tools", "tool_choice" };
  }
};//class AnthropicProtocol


//==============================================================================
//  OpenAI Responses API
//==============================================================================
class OpenAiResponsesProtocol : public LlmApiProtocol
{
  static void addConversation( const LlmInteraction &conv, json &input )
  {
    // Echo any captured reasoning items (id + encrypted_content) ahead of the items they produced -
    // required by the Responses API for reasoning models when store=false (else the follow-up turn
    // 400s with the reasoning items missing).  reasoningDetails carries them as a JSON array string.
    const auto echo_reasoning = []( const LlmInteractionTurn &turn, json &out ){
      const std::string &rd = turn.reasoningDetails();
      if( rd.empty() )
        return;
      try {
        const json items = json::parse( rd );
        if( items.is_array() )
        {
          for( const json &it : items )
            out.push_back( it );
        }
      } catch( ... ) { /* not Responses reasoning items (e.g. cross-provider history) - skip */ }
    };

    for( const std::shared_ptr<LlmInteractionTurn> &response : conv.responses )
    {
      if( response->excludeFromHistory() )
        continue;

      switch( response->type() )
      {
        case LlmInteractionTurn::Type::InitialRequest:
        {
          const LlmInteractionInitialRequest *req
            = dynamic_cast<const LlmInteractionInitialRequest*>( response.get() );
          if( !req )
            break;
          json item;
          item["role"] = "user";
          if( !req->imageContent().empty() )
          {
            json arr = json::array();
            json textBlock;
            textBlock["type"] = "input_text";
            textBlock["text"] = req->content();
            arr.push_back( textBlock );
            for( const LlmToolCall::ImageContent &img : req->imageContent() )
            {
              // Only benchmark questions set a caption, so pasted/uploaded images produce the
              //  exact same payload as before.
              if( !img.caption.empty() )
              {
                json capBlock;
                capBlock["type"] = "input_text";
                capBlock["text"] = img.caption;
                arr.push_back( capBlock );
              }

              json imgBlock;
              imgBlock["type"] = "input_image";
              imgBlock["image_url"] = "data:" + img.mimeType + ";base64," + img.base64Data;
              arr.push_back( imgBlock );
            }
            item["content"] = arr;
          }else
          {
            item["content"] = req->content();
          }
          input.push_back( item );
          break;
        }//InitialRequest

        case LlmInteractionTurn::Type::FinalLlmResponse:
        {
          const LlmInteractionFinalResponse *resp
            = dynamic_cast<const LlmInteractionFinalResponse*>( response.get() );
          if( !resp )
            break;
          echo_reasoning( *resp, input );
          json item;
          item["role"] = "assistant";
          item["content"] = resp->content();
          input.push_back( item );
          break;
        }//FinalLlmResponse

        case LlmInteractionTurn::Type::ToolCall:
        {
          const LlmToolRequest *req = dynamic_cast<const LlmToolRequest*>( response.get() );
          if( !req )
            break;
          // Replay the reasoning items that preceded these tool calls (required when store=false),
          // then emit each tool call as its own function_call input item.
          echo_reasoning( *req, input );
          for( const LlmToolCall &tc : req->toolCalls() )
          {
            json item;
            item["type"] = "function_call";
            item["call_id"] = tc.invocationId;
            item["name"] = tc.toolName;
            item["arguments"] = tc.toolParameters.is_null() ? string("{}") : tc.toolParameters.dump();
            input.push_back( item );
          }
          break;
        }//ToolCall

        case LlmInteractionTurn::Type::ToolResult:
        {
          const LlmToolResults *res = dynamic_cast<const LlmToolResults*>( response.get() );
          if( !res )
            break;
          for( const LlmToolCall &tc : res->toolCalls() )
          {
            json item;
            item["type"] = "function_call_output";
            item["call_id"] = tc.invocationId;
            item["output"] = tc.content;
            input.push_back( item );
          }
          break;
        }//ToolResult

        case LlmInteractionTurn::Type::Error:
        {
          const LlmInteractionError *err = dynamic_cast<const LlmInteractionError*>( response.get() );
          if( !err )
            break;
          json item;
          item["role"] = "assistant";
          item["content"] = "Error: " + err->errorMessage();
          input.push_back( item );
          break;
        }//Error

        case LlmInteractionTurn::Type::AutoReply:
        {
          const LlmInteractionAutoReply *ar = dynamic_cast<const LlmInteractionAutoReply*>( response.get() );
          if( !ar )
            break;
          json item;
          item["role"] = "user";
          item["content"] = ar->content();
          input.push_back( item );
          break;
        }//AutoReply

        case LlmInteractionTurn::Type::ConversationSummary:
        {
          const LlmConversationSummary *sum = dynamic_cast<const LlmConversationSummary*>( response.get() );
          if( !sum )
            break;
          json item;
          item["role"] = "user";
          item["content"] = "[System: Summary of previous conversation history]\n\n" + sum->summary();
          input.push_back( item );
          break;
        }//ConversationSummary
      }//switch( response->type() )
    }//for( responses )
  }//addConversation(...)

public:
  nlohmann::json serializeConversations(
      const std::vector<std::shared_ptr<LlmInteraction>> &convos ) const override
  {
    json input = json::array();
    for( const std::shared_ptr<LlmInteraction> &c : convos )
    {
      if( c )
        addConversation( *c, input );
    }
    // Merge consecutive same-role message items (function_call/output items have no role and are
    // left untouched).
    merge_adjacent_same_role( input, "input_text" );
    return input;
  }


  void appendEphemeralUserNote( nlohmann::json &input, const std::string &note ) const override
  {
    json item;
    item["role"] = "user";
    item["content"] = note;
    input.push_back( item );
  }


  /** No-op: the OpenAI Responses API caches automatically by longest-common-prefix. */
  void placeConversationCacheBreakpoints( nlohmann::json &/*input*/ ) const override
  {
  }


  nlohmann::json buildRequestBody( const LlmConfig::LlmApi &api,
                                   const std::string &systemPrompt,
                                   const nlohmann::json &messages,
                                   const std::vector<NormalizedTool> &tools,
                                   ToolChoice toolChoice ) const override
  {
    json request;
    request["model"] = api.model();

    if( !systemPrompt.empty() )
      request["instructions"] = systemPrompt;

    request["input"] = messages;

    if( api.maxTokens() > 0 )
      request["max_output_tokens"] = api.maxTokens();

    const std::variant<bool,LlmConfig::LlmApi::ReasoningEffort> &reasoning = api.reasoning();
    if( std::holds_alternative<bool>(reasoning) )
    {
      if( std::get<bool>(reasoning) )
        request["reasoning"]["effort"] = "medium";
    }else if( std::holds_alternative<LlmConfig::LlmApi::ReasoningEffort>(reasoning) )
    {
      switch( std::get<LlmConfig::LlmApi::ReasoningEffort>(reasoning) )
      {
        case LlmConfig::LlmApi::ReasoningEffort::low:    request["reasoning"]["effort"] = "low";    break;
        case LlmConfig::LlmApi::ReasoningEffort::medium: request["reasoning"]["effort"] = "medium"; break;
        case LlmConfig::LlmApi::ReasoningEffort::high:   request["reasoning"]["effort"] = "high";   break;
      }
    }

    // With store=false (below) and a reasoning model, the Responses API requires the reasoning items
    // produced before each tool call to be replayed on the following turn - otherwise it 400s
    // ("reasoning item ... without its required following item").  Ask for the encrypted reasoning so
    // serializeConversations() can echo it back across tool rounds (see addConversation / parseResponse).
    if( request.contains("reasoning") )
      request["include"] = json::array({ "reasoning.encrypted_content" });

    if( api.temperature().has_value() )
      request["temperature"] = api.temperature().value();

    if( !tools.empty() )
    {
      json toolsArr = json::array();
      for( const NormalizedTool &t : tools )
      {
        // The Responses API flattens function tools (name/description/parameters at the top level).
        json toolDef;
        toolDef["type"] = "function";
        toolDef["name"] = t.name;
        toolDef["description"] = t.description;
        toolDef["parameters"] = t.parameters;
        toolsArr.push_back( toolDef );
      }
      request["tools"] = toolsArr;
      request["tool_choice"] = (toolChoice == ToolChoice::Required) ? "required" : "auto";
    }

    // Stay stateless: do not persist server-side response state (we resend full input each turn).
    request["store"] = false;

    return request;
  }//buildRequestBody(...)


  ParsedLlmResponse parseResponse( const nlohmann::json &responseJson ) const override
  {
    ParsedLlmResponse out;

    if( responseJson.contains("usage") && responseJson["usage"].is_object() )
    {
      const json &usage = responseJson["usage"];
      out.promptTokens = get_usage_count( usage, "input_tokens" );
      out.completionTokens = get_usage_count( usage, "output_tokens" );
      out.totalTokens = get_usage_count( usage, "total_tokens" );
      if( !out.totalTokens.has_value() && out.promptTokens.has_value() && out.completionTokens.has_value() )
        out.totalTokens = out.promptTokens.value() + out.completionTokens.value();

      if( usage.contains("input_tokens_details") && usage["input_tokens_details"].is_object() )
        out.cachedTokens = get_usage_count( usage["input_tokens_details"], "cached_tokens" );
    }

    out.stopReason = responseJson.value("status","");

    if( responseJson.contains("output") && responseJson["output"].is_array() )
    {
      json reasoningItems = json::array();
      for( const json &item : responseJson["output"] )
      {
        const string itemType = item.value("type","");
        if( itemType == "message" )
        {
          if( item.contains("content") && item["content"].is_array() )
          {
            for( const json &block : item["content"] )
            {
              const string blockType = block.value("type","");
              if( (blockType == "output_text") || (blockType == "text") )
                out.content += block.value("text","");
            }
          }else if( item.contains("content") && item["content"].is_string() )
          {
            out.content += item["content"].get<string>();
          }
        }else if( itemType == "function_call" )
        {
          ParsedLlmResponse::ToolCall call;
          call.id = item.contains("call_id") ? item.value("call_id","") : item.value("id","");
          call.name = item.value("name","");
          if( item.contains("arguments") )
            call.argsRaw = item["arguments"].is_string() ? item["arguments"].get<string>() : item["arguments"].dump();
          out.toolCalls.push_back( std::move(call) );
        }else if( itemType == "reasoning" )
        {
          // Capture the reasoning item verbatim (incl. id + encrypted_content, requested via
          // include=["reasoning.encrypted_content"]) so it can be echoed back across tool rounds when
          // store=false (see addConversation).  Also surface any summary text for display.
          reasoningItems.push_back( item );
          if( item.contains("summary") && item["summary"].is_array() )
          {
            for( const json &s : item["summary"] )
            {
              if( s.is_object() && s.contains("text") && s["text"].is_string() )
                out.thinkingContent += s["text"].get<string>();
            }
          }
        }
      }//for( output items )

      // Carry the reasoning items (as a JSON array string) on the turn's reasoningDetails so they are
      // replayed before the items they produced on the next request.
      if( !reasoningItems.empty() )
        out.reasoningDetails = reasoningItems.dump();
    }//if( output array )

    if( out.content.empty() && responseJson.contains("output_text")
       && responseJson["output_text"].is_string() )
    {
      out.content = responseJson["output_text"].get<string>();
    }

    return out;
  }//parseResponse(...)


  std::vector<std::pair<std::string,std::string>> headers(
      const std::string & /*endpoint*/, const std::string &bearerToken ) const override
  {
    std::vector<std::pair<std::string,std::string>> hdrs;
    if( !bearerToken.empty() )
      hdrs.emplace_back( "Authorization", "Bearer " + bearerToken );
    return hdrs;
  }//headers(...)


  std::vector<std::string> cachePriorityFields() const override
  {
    return { "model", "max_output_tokens", "reasoning", "instructions", "tools", "tool_choice", "store" };
  }
};//class OpenAiResponsesProtocol

}//anonymous namespace


std::unique_ptr<const LlmApiProtocol> LlmApiProtocol::create( LlmConfig::LlmApi::ApiFormat fmt )
{
  switch( fmt )
  {
    case LlmConfig::LlmApi::ApiFormat::OpenAiChat:      return std::make_unique<OpenAiChatProtocol>();
    case LlmConfig::LlmApi::ApiFormat::OpenAiResponses: return std::make_unique<OpenAiResponsesProtocol>();
    case LlmConfig::LlmApi::ApiFormat::Anthropic:       return std::make_unique<AnthropicProtocol>();
  }
  throw std::logic_error( "LlmApiProtocol::create: unknown ApiFormat" );
}//LlmApiProtocol::create(...)


nlohmann::json strip_top_level_schema_combinators( nlohmann::json schema )
{
  if( schema.is_object() )
  {
    for( const char *key : { "oneOf", "anyOf", "allOf", "not", "if", "then", "else" } )
      schema.erase( key );
  }
  return schema;
}//strip_top_level_schema_combinators(...)
