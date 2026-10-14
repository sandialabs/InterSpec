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
#include <fstream>
#include <stdexcept>

#define BOOST_TEST_MODULE LlmKnowledgeBase_suite
#include <boost/test/included/unit_test.hpp>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/Filesystem.h"

#include "InterSpec/LlmKnowledgeBase.h"

using namespace std;
using json = nlohmann::json;

namespace
{
void write_file( const string &path, const string &content )
{
  ofstream out( path.c_str(), ios::binary | ios::out );
  BOOST_REQUIRE( out.is_open() );
  out << content;
}

// Creates a temporary llm_knowledge directory with one skill that has two references, and one without.
string make_knowledge_dir()
{
  const string base = SpecUtils::temp_file_name( "test_llm_knowledge", SpecUtils::temp_dir() );
  BOOST_REQUIRE( SpecUtils::create_directory( base ) == 1 );

  const string skill_dir = SpecUtils::append_path( base, "detector-models" );
  BOOST_REQUIRE( SpecUtils::create_directory( skill_dir ) == 1 );
  const string ref_dir = SpecUtils::append_path( skill_dir, "references" );
  BOOST_REQUIRE( SpecUtils::create_directory( ref_dir ) == 1 );

  write_file( SpecUtils::append_path( skill_dir, "SKILL.md" ),
              "---\nname: detector-models\ndescription: Model specific detector info.\n---\n\n# Detectors\nIndex body.\n" );
  write_file( SpecUtils::append_path( ref_dir, "identifinder-family.md" ),
              "---\ndescription: \"identiFINDER variants\"\n---\n# ignored heading\nNaI 1.4x2 in.\n" );
  write_file( SpecUtils::append_path( ref_dir, "detective-family.md" ),
              "# Detective family details\nHPGe 65x50 mm.\n" );
  write_file( SpecUtils::append_path( ref_dir, "notes.txt" ), "not a reference" );

  const string other_dir = SpecUtils::append_path( base, "Neutron_Notes" );
  BOOST_REQUIRE( SpecUtils::create_directory( other_dir ) == 1 );
  write_file( SpecUtils::append_path( other_dir, "SKILL.md" ), "Just a body, no frontmatter.\n" );

  return base;
}//make_knowledge_dir()
}//namespace


BOOST_AUTO_TEST_CASE( skill_references )
{
  const string base = make_knowledge_dir();

  const vector<LlmKnowledgeBase::SkillInfo> skills = LlmKnowledgeBase::loadSkillsFromKnowledgeDir( base );
  BOOST_REQUIRE_EQUAL( skills.size(), 2 );

  // Sorted by path, so "Neutron_Notes" comes first
  const LlmKnowledgeBase::SkillInfo &plain = skills[0];
  const LlmKnowledgeBase::SkillInfo &det = skills[1];
  BOOST_CHECK_EQUAL( plain.canonical_name, "neutron-notes" );
  BOOST_CHECK( plain.references.empty() );
  BOOST_CHECK_EQUAL( det.canonical_name, "detector-models" );
  BOOST_CHECK_EQUAL( det.description, "Model specific detector info." );

  BOOST_REQUIRE_EQUAL( det.references.size(), 2 );
  BOOST_CHECK_EQUAL( det.references[0].name, "detective-family" );
  BOOST_CHECK_EQUAL( det.references[0].summary, "Detective family details" );
  BOOST_CHECK_EQUAL( det.references[1].name, "identifinder-family" );
  BOOST_CHECK_EQUAL( det.references[1].summary, "identiFINDER variants" );

  // Schema: no `reference` parameter without references; enum of names otherwise
  const json plain_schema = LlmKnowledgeBase::skillToolParameterSchema( plain );
  BOOST_CHECK( plain_schema["properties"].empty() );

  const json det_schema = LlmKnowledgeBase::skillToolParameterSchema( det );
  BOOST_REQUIRE( det_schema["properties"].contains( "reference" ) );
  const json &names = det_schema["properties"]["reference"]["enum"];
  BOOST_REQUIRE_EQUAL( names.size(), 2 );
  BOOST_CHECK_EQUAL( names[0].get<string>(), "detective-family" );

  // Calling without a reference gives body plus index
  const json index = LlmKnowledgeBase::executeSkillTool( det, json::object() );
  BOOST_CHECK( index["body"].get<string>().find( "Index body." ) != string::npos );
  BOOST_REQUIRE( index.contains( "available_references" ) );
  BOOST_CHECK_EQUAL( index["available_references"].size(), 2 );

  const json plain_res = LlmKnowledgeBase::executeSkillTool( plain, json::object() );
  BOOST_CHECK( !plain_res.contains( "available_references" ) );

  // Loading a reference (case-insensitive) returns its full content
  const json ref = LlmKnowledgeBase::executeSkillTool( det, json{ {"reference", "IdentiFINDER-Family"} } );
  BOOST_CHECK_EQUAL( ref["reference"].get<string>(), "identifinder-family" );
  BOOST_CHECK( ref["body"].get<string>().find( "NaI 1.4x2 in." ) != string::npos );

  // Unknown names, and anything path-like, are rejected
  BOOST_CHECK_THROW( LlmKnowledgeBase::executeSkillTool( det, json{ {"reference", "nope"} } ), std::exception );
  BOOST_CHECK_THROW( LlmKnowledgeBase::executeSkillTool( det, json{ {"reference", "../SKILL"} } ), std::exception );
  BOOST_CHECK_THROW( LlmKnowledgeBase::executeSkillTool( det, json{ {"reference", "notes"} } ), std::exception );

  // Tool names and the invoke_KnowledgeBase description suffix
  BOOST_CHECK_EQUAL( LlmKnowledgeBase::skillToolNameFromCanonical( det.canonical_name ), "knowledge_detector_models" );

  const string suffix = LlmKnowledgeBase::invokeDescriptionSuffix( skills, false );
  BOOST_CHECK( suffix.find( "detector-models" ) != string::npos );
  BOOST_CHECK( suffix.find( "Neutron_Notes" ) != string::npos ); //no frontmatter name, so the directory name is used
  BOOST_CHECK( suffix.find( "Model specific detector info." ) == string::npos ); //descriptions are not included
  BOOST_CHECK( suffix.find( "JavaScript" ) != string::npos );
  BOOST_CHECK( suffix.find( "RAG" ) == string::npos );
  BOOST_CHECK( LlmKnowledgeBase::invokeDescriptionSuffix( skills, true ).find( "RAG" ) != string::npos );
  BOOST_CHECK( LlmKnowledgeBase::invokeDescriptionSuffix( {}, true ).find( "curated reference knowledge" ) == string::npos );

  SpecUtils::remove_file( SpecUtils::append_path( SpecUtils::append_path( base, "detector-models" ), "references/notes.txt" ) );
}//BOOST_AUTO_TEST_CASE( skill_references )
