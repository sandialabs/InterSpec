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

#include <chrono>
#include <string>
#include <iostream>

#define BOOST_TEST_MODULE LlmInterface_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/LlmInterface.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;
using namespace boost::unit_test;
using json = nlohmann::json;


BOOST_AUTO_TEST_SUITE( JsonParsing )

// Test Case 1: Valid JSON (fast path - no repair needed)
BOOST_AUTO_TEST_CASE( ValidJsonNoRepair )
{
  // Simple object
  {
    const string jsonStr = R"({"key": "value"})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Nested object
  {
    const string jsonStr = R"({"nested": {"key": "value"}})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["nested"]["key"].get<string>(), "value" );
  }

  // Array
  {
    const string jsonStr = R"([1, 2, 3])";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed.size(), 3 );
    BOOST_CHECK_EQUAL( parsed[0].get<int>(), 1 );
    BOOST_CHECK_EQUAL( parsed[1].get<int>(), 2 );
    BOOST_CHECK_EQUAL( parsed[2].get<int>(), 3 );
  }

  // Complex structure
  {
    const string jsonStr = R"({"array": [1, 2, {"nested": "value"}], "bool": true})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["array"][2]["nested"].get<string>(), "value" );
    BOOST_CHECK_EQUAL( parsed["bool"].get<bool>(), true );
  }
}

// Test Case 2: Missing closing quote
BOOST_AUTO_TEST_CASE( MissingClosingQuote )
{
  // Simple case
  {
    const string jsonStr = R"({"key": "value)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Multiple fields, last one incomplete
  {
    const string jsonStr = R"({"a":"b","c":"d)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["a"].get<string>(), "b" );
    BOOST_CHECK_EQUAL( parsed["c"].get<string>(), "d" );
  }

  // The actual failing case from the user's issue
  {
    const string jsonStr = R"({"editAction":"SetSource","energy":45.58,"stringValue":"Eu152)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["editAction"].get<string>(), "SetSource" );
    BOOST_CHECK_CLOSE( parsed["energy"].get<double>(), 45.58, 0.001 );
    BOOST_CHECK_EQUAL( parsed["stringValue"].get<string>(), "Eu152" );
  }
}

// Test Case 3: Missing closing brace
BOOST_AUTO_TEST_CASE( MissingClosingBrace )
{
  // Simple object
  {
    const string jsonStr = R"({"key": "value")";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Nested objects
  {
    const string jsonStr = R"({"a":{"b":"c")";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["a"]["b"].get<string>(), "c" );
  }

  // Multiple levels
  {
    const string jsonStr = R"({"level1":{"level2":{"level3":"value")";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["level1"]["level2"]["level3"].get<string>(), "value" );
  }
}

// Test Case 4: Missing closing bracket
BOOST_AUTO_TEST_CASE( MissingClosingBracket )
{
  // Simple array
  {
    const string jsonStr = R"([1, 2, 3)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed.size(), 3 );
    BOOST_CHECK_EQUAL( parsed[0].get<int>(), 1 );
    BOOST_CHECK_EQUAL( parsed[1].get<int>(), 2 );
    BOOST_CHECK_EQUAL( parsed[2].get<int>(), 3 );
  }

  // Array of objects
  {
    const string jsonStr = R"([{"a":"b"},{"c":"d")";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed.size(), 2 );
    BOOST_CHECK_EQUAL( parsed[0]["a"].get<string>(), "b" );
    BOOST_CHECK_EQUAL( parsed[1]["c"].get<string>(), "d" );
  }
}

// Test Case 5: Combined issues (missing quote and brace/bracket)
BOOST_AUTO_TEST_CASE( CombinedIssues )
{
  // Missing quote and brace
  {
    const string jsonStr = R"({"key":"val)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "val" );
  }

  // Nested with multiple issues
  {
    const string jsonStr = R"({"a":[1,2,{"b":"c)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["a"].size(), 3 );
    BOOST_CHECK_EQUAL( parsed["a"][0].get<int>(), 1 );
    BOOST_CHECK_EQUAL( parsed["a"][1].get<int>(), 2 );
    BOOST_CHECK_EQUAL( parsed["a"][2]["b"].get<string>(), "c" );
  }

  // Complex real-world case
  {
    const string jsonStr = R"({"tool":"analyze","params":{"spectrum":"fg","energy":661.7,"window":[655,668)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["tool"].get<string>(), "analyze" );
    BOOST_CHECK_CLOSE( parsed["params"]["energy"].get<double>(), 661.7, 0.001 );
    BOOST_CHECK_EQUAL( parsed["params"]["window"][0].get<int>(), 655 );
    BOOST_CHECK_EQUAL( parsed["params"]["window"][1].get<int>(), 668 );
  }
}

// Test Case 6: Escaped characters
BOOST_AUTO_TEST_CASE( EscapedCharacters )
{
  // Escaped quote in incomplete string
  {
    const string jsonStr = R"({"text":"He said \"hi)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["text"].get<string>(), "He said \"hi" );
  }

  // Backslash
  {
    const string jsonStr = R"({"path":"C:\\Users\\file)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["path"].get<string>(), "C:\\Users\\file" );
  }

  // Multiple escapes
  {
    const string jsonStr = R"({"regex":"\\d+\\s*\\w+)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["regex"].get<string>(), "\\d+\\s*\\w+" );
  }
}

// Test Case 7: Edge cases
BOOST_AUTO_TEST_CASE( EdgeCases )
{
  // Empty string
  {
    const string jsonStr = "";
    BOOST_CHECK_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ), json::parse_error );
  }

  // Just opening brace
  {
    const string jsonStr = "{";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK( parsed.is_object() );
    BOOST_CHECK( parsed.empty() );
  }

  // Just opening bracket
  {
    const string jsonStr = "[";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK( parsed.is_array() );
    BOOST_CHECK( parsed.empty() );
  }

  // Deeply nested
  {
    const string jsonStr = "[[[[1,2";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK( parsed.is_array() );
    BOOST_CHECK_EQUAL( parsed[0][0][0].size(), 2 );
  }
}

// Test Case 8: Unfixable cases (should still throw proper errors)
BOOST_AUTO_TEST_CASE( UnfixableCases )
{
  // Invalid syntax (wrong bracket type matching)
  {
    const string jsonStr = R"({"bad": [})";
    BOOST_CHECK_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ), json::parse_error );
  }

  // Unquoted key
  {
    const string jsonStr = R"({key: "value"})";
    BOOST_CHECK_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ), json::parse_error );
  }

  // Invalid number
  {
    const string jsonStr = R"({"num": 123abc})";
    BOOST_CHECK_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ), json::parse_error );
  }
}

// Test Case 9: Existing sanitization (UTF-8 BOM, whitespace)
BOOST_AUTO_TEST_CASE( ExistingSanitization )
{
  // UTF-8 BOM
  {
    const string jsonStr = "\xEF\xBB\xBF{\"key\":\"value\"}";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Leading/trailing whitespace
  {
    const string jsonStr = "  \n\t{\"key\":\"value\"}\n  ";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( jsonStr ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( jsonStr );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }
}

// Test Case 10: Performance test (no regression for valid JSON)
BOOST_AUTO_TEST_CASE( PerformanceNoRegression )
{
  // Large valid JSON should not have performance penalty
  json large_obj;
  for( int i = 0; i < 1000; ++i )
  {
    large_obj["field_" + std::to_string( i )] = "value_" + std::to_string( i );
  }

  const string json_str = large_obj.dump();

  // Time the parsing (should be fast - single parse call)
  const auto start = std::chrono::steady_clock::now();
  BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( json_str ) );
  const auto end = std::chrono::steady_clock::now();

  const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( end - start );
  BOOST_TEST_MESSAGE( "Large valid JSON parsed in " << duration.count() << "ms" );

  // Should be very fast (< 100ms) since no repair attempted
  BOOST_CHECK_LT( duration.count(), 100 );
}

// Test Case 11: Repair log functionality
BOOST_AUTO_TEST_CASE( RepairLogOutput )
{
  // Test that repair log is populated correctly
  {
    const string jsonStr = R"({"key":"val)";
    string repairLog;
    const string repaired = LlmInterfaceTests::repairIncompleteJson( jsonStr, &repairLog );

    BOOST_CHECK( !repairLog.empty() );
    BOOST_CHECK( repairLog.find( "closing quote" ) != string::npos );
    BOOST_CHECK( repairLog.find( "closing brace" ) != string::npos );

    // Verify repaired string is valid JSON
    BOOST_REQUIRE_NO_THROW( json::parse( repaired ) );
  }

  // Test no repair log for valid JSON
  {
    const string jsonStr = R"({"key":"val"})";
    string repairLog;
    const string repaired = LlmInterfaceTests::repairIncompleteJson( jsonStr, &repairLog );

    BOOST_CHECK( repairLog.empty() );
    BOOST_CHECK_EQUAL( repaired, jsonStr );
  }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================================
// New test suites for the additional repair strategies added to lenientlyParseJson
// and sanitizeJsonString.
//
// Several test cases below were inspired by the test suite of the open-source
// TypeScript library "ai-json-fixer" (https://github.com/aotakeda/ai-json-fixer,
// MIT licence).  The specific cases were translated to C++ and adapted for the
// InterSpec repair pipeline.
// ============================================================================


// ---------------------------------------------------------------------------
// Trailing junk removal
// Some LLMs append sentinel tokens or plain text after the closing delimiter.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( TrailingContentRemoval )

BOOST_AUTO_TEST_CASE( TrailingTextAfterObject )
{
  // Inspired by ai-json-fixer trailing-content-removal tests
  {
    const string s = R"({"status": "complete", "count": 42})" "\n\nThat's the JSON response you requested.";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["status"].get<string>(), "complete" );
    BOOST_CHECK_EQUAL( parsed["count"].get<int>(), 42 );
  }

  {
    const string s = R"({"data": {"user": "john", "age": 30}} This represents the user data structure.)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["data"]["user"].get<string>(), "john" );
  }

  // Brackets inside string values must not confuse the depth tracker
  {
    const string s = R"({"text": "Use [brackets] and {braces} carefully"} Note the special characters.)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["text"].get<string>(), "Use [brackets] and {braces} carefully" );
  }

  // Escaped quotes inside strings must not confuse the depth tracker
  {
    const string s = "{\"message\": \"He said \\\"hello\\\" to me\"} - Contains escaped quotes";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["message"].get<string>(), "He said \"hello\" to me" );
  }
}

BOOST_AUTO_TEST_CASE( TrailingTextAfterArray )
{
  // Inspired by ai-json-fixer trailing-content-removal tests
  {
    const string s = "[1, 2, 3, 4, 5]\n\nThese are the first five numbers.";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed.size(), 5 );
  }

  {
    const string s = R"([{"id": 1}, {"id": 2}] - Two items total)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed.size(), 2 );
    BOOST_CHECK_EQUAL( parsed[0]["id"].get<int>(), 1 );
  }
}

BOOST_AUTO_TEST_CASE( LlmSentinelTokens )
{
  // The original motivating case: LLM appends "<|call|>" after the JSON object
  {
    const string s = "{\n"
      "\"context\": \"Foreground spectrum currently loaded.\",\n"
      "\"task\": \"Perform energy calibration.\"\n"
      "}    <|call|>";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["task"].get<string>(), "Perform energy calibration." );
  }

  // Another common sentinel
  {
    const string s = R"({"action": "done"})" " <|end|>";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["action"].get<string>(), "done" );
  }
}

BOOST_AUTO_TEST_SUITE_END() // TrailingContentRemoval


// ---------------------------------------------------------------------------
// Leading garbage removal
// Some LLMs prepend XML-like tags, markdown artefacts, or other non-JSON text.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( LeadingGarbageRemoval )

BOOST_AUTO_TEST_CASE( XmlTagPrefix )
{
  // Inspired by ai-json-fixer markdown-extraction tests
  {
    const string s = R"(<tool_call>{"name": "calibrate", "args": {}}</tool_call>)";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["name"].get<string>(), "calibrate" );
  }

  {
    const string s = "<json>\n{\"key\": \"value\"}\n</json>";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }
}

BOOST_AUTO_TEST_CASE( CombinedLeadingAndTrailing )
{
  // Leading tag + trailing sentinel
  {
    const string s = "<tool_call>{\"action\": \"fit_peaks\"}    <|call|>";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["action"].get<string>(), "fit_peaks" );
  }
}

BOOST_AUTO_TEST_SUITE_END() // LeadingGarbageRemoval


// ---------------------------------------------------------------------------
// Markdown code-fence extraction
// LLMs commonly wrap JSON in ```json ... ``` blocks.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( MarkdownExtraction )

BOOST_AUTO_TEST_CASE( LabelledJsonFence )
{
  // Inspired by ai-json-fixer markdown-extraction tests
  {
    const string s = "```json\n{\"key\": \"value\"}\n```";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Uppercase label
  {
    const string s = "```JSON\n{\"key\": \"value\"}\n```";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Fence with surrounding prose (as an LLM would emit)
  {
    const string s = "Here is the result:\n```json\n{\"status\": \"ok\"}\n```\nI hope that helps!";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["status"].get<string>(), "ok" );
  }

  // Nested structure inside fence
  {
    const string s = "```json\n{\"nuclides\": [\"Cs137\", \"Co60\"], \"count\": 2}\n```";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["nuclides"][0].get<string>(), "Cs137" );
    BOOST_CHECK_EQUAL( parsed["count"].get<int>(), 2 );
  }
}

BOOST_AUTO_TEST_CASE( UnlabelledFence )
{
  // Unlabelled fence whose body starts with '{' — should be extracted
  {
    const string s = "```\n{\"key\": \"value\"}\n```";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }

  // Unlabelled fence whose body starts with '[' — should be extracted
  {
    const string s = "```\n[1, 2, 3]\n```";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed.size(), 3 );
  }
}

BOOST_AUTO_TEST_CASE( FenceWithCRLF )
{
  // Windows-style line endings inside a fence
  {
    const string s = "```json\r\n{\"key\": \"value\"}\r\n```";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["key"].get<string>(), "value" );
  }
}

BOOST_AUTO_TEST_SUITE_END() // MarkdownExtraction


// ---------------------------------------------------------------------------
// Missing comma insertion
// Some LLMs omit the comma between object fields or array elements.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( MissingCommaInsertion )

BOOST_AUTO_TEST_CASE( MissingCommaBetweenObjectFields )
{
  // Inspired by ai-json-fixer missing-comma-detection tests
  {
    const string s = "{\n  \"a\": 1\n  \"b\": 2\n}";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["a"].get<int>(), 1 );
    BOOST_CHECK_EQUAL( parsed["b"].get<int>(), 2 );
  }

  // String values
  {
    const string s = "{\n  \"name\": \"Alice\"\n  \"role\": \"analyst\"\n}";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["name"].get<string>(), "Alice" );
    BOOST_CHECK_EQUAL( parsed["role"].get<string>(), "analyst" );
  }

  // Boolean and null values
  {
    const string s = "{\n  \"active\": true\n  \"data\": null\n  \"count\": 5\n}";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["active"].get<bool>(), true );
    BOOST_CHECK( parsed["data"].is_null() );
    BOOST_CHECK_EQUAL( parsed["count"].get<int>(), 5 );
  }
}

BOOST_AUTO_TEST_CASE( MissingCommaBetweenArrayElements )
{
  // Inspired by ai-json-fixer missing-comma-detection tests
  {
    const string s = "[\n  {\"id\": 1}\n  {\"id\": 2}\n]";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed.size(), 2 );
    BOOST_CHECK_EQUAL( parsed[0]["id"].get<int>(), 1 );
    BOOST_CHECK_EQUAL( parsed[1]["id"].get<int>(), 2 );
  }
}

BOOST_AUTO_TEST_CASE( AlreadyValidJsonNotModified )
{
  // Correctly-comma'd JSON must survive unchanged
  {
    const string s = "{\n  \"a\": 1,\n  \"b\": 2\n}";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["a"].get<int>(), 1 );
    BOOST_CHECK_EQUAL( parsed["b"].get<int>(), 2 );
  }

  // String values that contain spaces or JSON-like text must not be corrupted
  {
    const string s = "{\n  \"desc\": \"a b c\",\n  \"note\": \"x y z\"\n}";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["desc"].get<string>(), "a b c" );
    BOOST_CHECK_EQUAL( parsed["note"].get<string>(), "x y z" );
  }

  // Empty containers
  {
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( "{}" ) );
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( "[]" ) );
  }
}

BOOST_AUTO_TEST_SUITE_END() // MissingCommaInsertion


// ---------------------------------------------------------------------------
// Regression: valid JSON must always parse correctly and identically
// Covers all JSON value types to ensure no new repair phase corrupts them.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( ValidJsonRegression )

BOOST_AUTO_TEST_CASE( AllPrimitiveTypes )
{
  // String
  {
    const json parsed = LlmInterfaceTests::lenientlyParseJson( R"("hello")" );
    BOOST_CHECK_EQUAL( parsed.get<string>(), "hello" );
  }

  // Integer
  {
    const json parsed = LlmInterfaceTests::lenientlyParseJson( "42" );
    BOOST_CHECK_EQUAL( parsed.get<int>(), 42 );
  }

  // Float
  {
    const json parsed = LlmInterfaceTests::lenientlyParseJson( "3.14" );
    BOOST_CHECK_CLOSE( parsed.get<double>(), 3.14, 0.001 );
  }

  // Boolean
  {
    const json t = LlmInterfaceTests::lenientlyParseJson( "true" );
    const json f = LlmInterfaceTests::lenientlyParseJson( "false" );
    BOOST_CHECK_EQUAL( t.get<bool>(), true );
    BOOST_CHECK_EQUAL( f.get<bool>(), false );
  }

  // Null
  {
    const json parsed = LlmInterfaceTests::lenientlyParseJson( "null" );
    BOOST_CHECK( parsed.is_null() );
  }
}

BOOST_AUTO_TEST_CASE( SpecialStringContents )
{
  // String with backticks (must not trigger markdown extraction)
  {
    const string s = R"({"cmd": "echo `date`"})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["cmd"].get<string>(), "echo `date`" );
  }

  // String that looks like trailing garbage but is a value
  {
    const string s = R"({"token": "<|end|>"})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["token"].get<string>(), "<|end|>" );
  }

  // String with embedded braces and brackets
  {
    const string s = R"({"expr": "{a} + [b]"})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["expr"].get<string>(), "{a} + [b]" );
  }

  // Unicode
  {
    const string s = "{\"msg\": \"caf\xC3\xA9\"}";  // café in UTF-8
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK( !parsed["msg"].get<string>().empty() );
  }
}

BOOST_AUTO_TEST_CASE( RealWorldToolCallArguments )
{
  // Typical tool-call argument objects that must round-trip perfectly
  {
    const string s = R"({
  "context": "Foreground spectrum currently loaded. Need to improve energy calibration.",
  "task": "Perform energy calibration using NORM peaks (K-40 1460 keV and Th-232 2614 keV)."
})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK( parsed["task"].get<string>().find( "K-40" ) != string::npos );
  }

  {
    const string s = R"({"editAction":"SetSource","energy":45.58,"stringValue":"Eu152"})";
    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );
    BOOST_CHECK_EQUAL( parsed["editAction"].get<string>(), "SetSource" );
    BOOST_CHECK_CLOSE( parsed["energy"].get<double>(), 45.58, 0.001 );
    BOOST_CHECK_EQUAL( parsed["stringValue"].get<string>(), "Eu152" );
  }
}

BOOST_AUTO_TEST_SUITE_END() // ValidJsonRegression
