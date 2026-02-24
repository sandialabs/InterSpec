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

  {
    const string s = R"({
  "choices": [
    {
      "finish_reason": "stop",
      "index": 0,
      "message": {
        "content": "**Validation Summary – Co‑60**\n\n1. **Energy match**  \n   *Peak 1173 keV* – fitted energy 1173.12 keV, FWHM ≈ 2.35 keV → within 1 FWHM of the Co‑60 line (1173.23 keV).  \n   *Peak 1332 keV* – fitted energy 1332.26 keV, FWHM ≈ 2.49 keV → within 1 FWHM of the Co‑60 line (1332.49 keV).\n\n2. **Source information** (`source_info`)  \n   – Co‑60 is a common industrial/sterilization gamma source.  \n   – Associated nuclides (Mn‑54, Zn‑65, Ag‑110m, etc.) are *not* observed, which is fine because they are only expected if the source is a mixed activation product; a pure Co‑60 calibration/source is typical.\n\n3. **Corroborating peaks** (`get_detected_peaks`)  \n   – Both primary gammas are present with comparable amplitudes (≈ 3.8 × 10⁶ counts each).  \n   – The intensity ratio is ~1 : 1, matching Co‑60's 100 % branching for both lines.  \n   – No additional strong Co‑60 lines are expected in this energy range, so the absence of other peaks is not a concern.\n\n4. **Escape‑peak consideration**  \n   – A possible single‑escape of the 1173 keV line would appear near 662 keV, but a very strong, distinct peak at 661.66 keV is present (amplitude 1.5 × 10⁶ cps). This is the well‑known Cs‑137 662 keV line, not an escape peak, and its intensity is far higher than would be expected from a Co‑60 escape (≈ 1 % of the parent). Therefore the 661 keV peak is **not** an escape artifact.\n\n5. **Relative‑efficiency check**  \n   – Only the two Co‑60 peaks are available; the relative‑efficiency tool requires ≥ 3 peaks for a polynomial fit, so a formal fit is unnecessary. The two‑peak agreement is already decisive.\n\n6. **Ultimate‑parent rule**  \n   – Co‑60 is itself a long‑lived parent (half‑life 5.27 y). No higher‑order parent is required.\n\n**Conclusion:** All validation criteria are satisfied. The identification of the 1173 keV and 1332 keV peaks as **Co‑60** is robust.\n\n**Next step:** Commit the source to the analysis state.\n\n---  \n\n**Transitioning workflow:**  \n\n```json\n{\n  \"notes\": \"Co-60 validated and ready to be added as a source peak.\",\n  \"state\": \"ADD_SOURCE_PEAK\"\n}\n```",
        "provider_specific_fields": {
          "reasoning_content": "We need to validate Co-60 identification. Steps:\n\n- Check source_info: done (Co-60 is plausible; associated sources include Mn-54 etc.)\n\n- Check other peaks: Co-60 also has 1173 and 1332 keV pairs (they are present). Check if any other peaks of Co-60 expected? Co-60 gamma lines: only these two main lines ~1173 and 1332 keV, plus maybe 44 keV? Actually Co-60 also has low-energy 0.0? It also emits a 0.0? Not relevant. So we have both lines, good.\n\n- Check ratio of amplitudes: Expect ratio close to 1 (branching 100% each). Observed amplitude: 3.84M and 3.74M, fairly similar.\n\n- Use relative efficiency: need enough peaks; we only have two; the rule says need at least as many peaks as equation order (default eqn order 3). Can't do rel eff. Could set eqn order lower? Could override options? The tool requires at least as many peaks as order; but we could request eqn_order=1 to fit with 2 peaks. But maybe we can skip relative efficiency for Co-60 due to limited peaks. However guidelines: if >2 peaks, use relative efficiency. Here only 2, okay.\n\n- Check for other sources overlapping: any other peaks near these energies? Not likely.\n\n- Check for escape peaks: 1173 and 1332 could have single-escape at ~1022 (511 subtractions) ~662 keV? Actually single escape is energy -511 keV: 1173-511=662 keV (close to 661.66 peak, which is present). Could that 661.66 be single escape of 1173? But 662 is a known background line from Bi-214 (609) ??? Actually 662 is Compton edge? Wait 661.66 is characteristic gamma of Cs-137 (662 keV). That is strong. Could also be sum of 511+151? But likely it's Cs-137. So need to check possible escape peaks.\n\nCheck escape_peak_check for 662? But we can note that 661.66 peak likely Cs-137.\n\nThus Co-60 identification seems solid.\n\nNow we must transition to ADD_SOURCE_PEAK.\n\n"
        },
        "reasoning_content": "We need to validate Co-60 identification. Steps:\n\n- Check source_info: done (Co-60 is plausible; associated sources include Mn-54 etc.)\n\n- Check other peaks: Co-60 also has 1173 and 1332 keV pairs (they are present). Check if any other peaks of Co-60 expected? Co-60 gamma lines: only these two main lines ~1173 and 1332 keV, plus maybe 44 keV? Actually Co-60 also has low-energy 0.0? It also emits a 0.0? Not relevant. So we have both lines, good.\n\n- Check ratio of amplitudes: Expect ratio close to 1 (branching 100% each). Observed amplitude: 3.84M and 3.74M, fairly similar.\n\n- Use relative efficiency: need enough peaks; we only have two; the rule says need at least as many peaks as equation order (default eqn order 3). Can't do rel eff. Could set eqn order lower? Could override options? The tool requires at least as many peaks as order; but we could request eqn_order=1 to fit with 2 peaks. But maybe we can skip relative efficiency for Co-60 due to limited peaks. However guidelines: if >2 peaks, use relative efficiency. Here only 2, okay.\n\n- Check for other sources overlapping: any other peaks near these energies? Not likely.\n\n- Check for escape peaks: 1173 and 1332 could have single-escape at ~1022 (511 subtractions) ~662 keV? Actually single escape is energy -511 keV: 1173-511=662 keV (close to 661.66 peak, which is present). Could that 661.66 be single escape of 1173? But 662 is a known background line from Bi-214 (609) ??? Actually 662 is Compton edge? Wait 661.66 is characteristic gamma of Cs-137 (662 keV). That is strong. Could also be sum of 511+151? But likely it's Cs-137. So need to check possible escape peaks.\n\nCheck escape_peak_check for 662? But we can note that 661.66 peak likely Cs-137.\n\nThus Co-60 identification seems solid.\n\nNow we must transition to ADD_SOURCE_PEAK.\n\n",
        "role": "assistant"
      }
    }
  ],
  "created": 177197,
  "id": "chatcmpl-...",
  "model": "openai/gpt-oss-120b",
  "object": "chat.completion",
  "usage": {
    "completion_tokens": 1130,
    "prompt_tokens": 16936,
    "total_tokens": 18066
  }
})";

    BOOST_REQUIRE_NO_THROW( LlmInterfaceTests::lenientlyParseJson( s ) );
    const json parsed = LlmInterfaceTests::lenientlyParseJson( s );

    BOOST_REQUIRE( parsed.contains("choices") && parsed["choices"].is_array() && (parsed["choices"].size() == 1) );

    BOOST_REQUIRE( parsed["choices"][0].contains("index") && parsed["choices"][0]["index"].is_number() );
    BOOST_CHECK_EQUAL( parsed["choices"][0]["index"].get<int>(), 0 );

    BOOST_REQUIRE( parsed["choices"][0].contains("finish_reason") && parsed["choices"][0]["finish_reason"].is_string() );
    BOOST_CHECK_EQUAL( parsed["choices"][0]["finish_reason"].get<string>(), "stop" );
  }
}

BOOST_AUTO_TEST_SUITE_END() // ValidJsonRegression
