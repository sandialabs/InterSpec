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
#include <stdexcept>

#define BOOST_TEST_MODULE LlmBenchmarkJudge_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/LlmBenchmarkRunner.h"

using namespace std;

// Parsing of the judge LLM's verdict (LlmBenchmarkRunner::parseJudgeResponse); the judge's reply
//  format varies by model, so these are replies actually seen in benchmark runs.

typedef LlmBenchmarkRunner::JudgeVerdict JudgeVerdict;


BOOST_AUTO_TEST_CASE( PlainJson )
{
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse(
    "{\"correct\": true, \"score\": 3, \"maxScore\": 4, \"reason\": \"Most points met.\"}", 4 );
  BOOST_CHECK( v.correct );
  BOOST_CHECK_EQUAL( v.score, 3 );
  BOOST_CHECK_EQUAL( v.maxScore, 4 );
  BOOST_CHECK_EQUAL( v.reason, "Most points met." );
}


BOOST_AUTO_TEST_CASE( MissingMaxScoreUsesRubric )
{
  // Seen from Sonnet 5.5 as judge: no "maxScore"; it used to default to 1, giving "4/1"
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse(
    "{\"correct\": true, \"score\": 4, \"reason\": \"All criteria met.\"}", 4 );
  BOOST_CHECK_EQUAL( v.score, 4 );
  BOOST_CHECK_EQUAL( v.maxScore, 4 );
}


BOOST_AUTO_TEST_CASE( RubricOverridesJudgeMax )
{
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse(
    "{\"correct\": true, \"score\": 5, \"maxScore\": 6, \"reason\": \"x\"}", 4 );
  BOOST_CHECK_EQUAL( v.maxScore, 4 );
  BOOST_CHECK_EQUAL( v.score, 4 );   // clamped to the rubric total
}


BOOST_AUTO_TEST_CASE( NoRubricUsesJudgeMax )
{
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse(
    "{\"correct\": false, \"score\": 0, \"maxScore\": 1, \"reason\": \"Wrong nuclide.\"}", 0 );
  BOOST_CHECK( !v.correct );
  BOOST_CHECK_EQUAL( v.score, 0 );
  BOOST_CHECK_EQUAL( v.maxScore, 1 );

  const JudgeVerdict w = LlmBenchmarkRunner::parseJudgeResponse( "{\"correct\": true, \"reason\": \"ok\"}", 0 );
  BOOST_CHECK_EQUAL( w.maxScore, 1 );
  BOOST_CHECK_EQUAL( w.score, 1 );   // no score given: correct means full marks
}


BOOST_AUTO_TEST_CASE( ReasoningBeforeJson )
{
  // Seen from Sonnet 5.5 as judge: criterion-by-criterion reasoning, then the JSON
  const string reply = "Criterion 1 (2 points): The answer names Be. Met.\n"
    "Criterion 2 (2 points): Gives 4438 keV {from 12C*}. Met.\n\n"
    "{\"correct\": true, \"score\": 4, \"maxScore\": 4, \"reason\": \"Identifies Be and the 4438 keV line.\"}";
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse( reply, 4 );
  BOOST_CHECK_EQUAL( v.score, 4 );
  BOOST_CHECK_EQUAL( v.reason, "Identifies Be and the 4438 keV line." );
}


BOOST_AUTO_TEST_CASE( MarkdownFenceAndTrailingText )
{
  const string reply = "```json\n{\"correct\": false, \"score\": 1, \"maxScore\": 4, \"reason\": \"Partial.\"}\n```\nDone.";
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse( reply, 4 );
  BOOST_CHECK( !v.correct );
  BOOST_CHECK_EQUAL( v.score, 1 );
}


BOOST_AUTO_TEST_CASE( LenientNumbers )
{
  BOOST_CHECK_EQUAL( LlmBenchmarkRunner::parseJudgeResponse( "{\"score\": 2.6, \"reason\": \"\"}", 4 ).score, 3 );
  BOOST_CHECK_EQUAL( LlmBenchmarkRunner::parseJudgeResponse( "{\"score\": \"3\", \"reason\": \"\"}", 4 ).score, 3 );
  BOOST_CHECK_EQUAL( LlmBenchmarkRunner::parseJudgeResponse( "{\"score\": -1, \"reason\": \"\"}", 4 ).score, 0 );
}


BOOST_AUTO_TEST_CASE( NestedObjectInVerdict )
{
  const string reply = "{\"correct\": true, \"score\": 2, \"maxScore\": 2, \"detail\": {\"c1\": 1, \"c2\": 1}, \"reason\": \"ok\"}";
  const JudgeVerdict v = LlmBenchmarkRunner::parseJudgeResponse( reply, 2 );
  BOOST_CHECK_EQUAL( v.score, 2 );
}


BOOST_AUTO_TEST_CASE( NoVerdictThrows )
{
  BOOST_CHECK_THROW( LlmBenchmarkRunner::parseJudgeResponse( "I could not grade this.", 4 ), std::runtime_error );
  BOOST_CHECK_THROW( LlmBenchmarkRunner::parseJudgeResponse( "{\"note\": \"no verdict\"}", 4 ), std::runtime_error );
  BOOST_CHECK_THROW( LlmBenchmarkRunner::parseJudgeResponse( "{\"score\": 3, \"reason\": ", 4 ), std::runtime_error );
}
