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
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

#define BOOST_TEST_MODULE LlmKnowledgeBase_suite
#include <boost/test/included/unit_test.hpp>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
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


// Checks every skill shipped in `<datadir>/llm_knowledge` against the loader's limits, so a
// mis-authored skill (multi-line frontmatter, oversized reference, ...) fails here, not silently in use.
BOOST_AUTO_TEST_CASE( shipped_skills_well_formed )
{
  const string datadir = find_data_dir();
  BOOST_REQUIRE_MESSAGE( !datadir.empty(), "Could not find the static data directory; pass --datadir=..." );

  const string knowledge_dir = SpecUtils::append_path( datadir, "llm_knowledge" );
  const vector<LlmKnowledgeBase::SkillInfo> skills = LlmKnowledgeBase::loadSkillsFromKnowledgeDir( knowledge_dir );
  if( skills.empty() )
  {
    BOOST_TEST_MESSAGE( "No skills found in '" << knowledge_dir << "'." );
    return;
  }

  for( const LlmKnowledgeBase::SkillInfo &skill : skills )
  {
    BOOST_TEST_CONTEXT( "Skill '" << skill.canonical_name << "' (" << skill.file_path << ")" )
    {
      // Descriptions must be a single frontmatter line; YAML block scalars (e.g., ">-") are not parsed.
      BOOST_CHECK( !skill.description.empty() );
      BOOST_CHECK( skill.description.find( "Loaded from " ) != 0 );
      BOOST_CHECK( skill.description != ">-" && skill.description != ">" && skill.description != "|" );
      BOOST_CHECK_LE( skill.description.size(), 400 );
      BOOST_CHECK( !skill.body.empty() );

      // LLM APIs commonly limit tool names to 64 characters
      BOOST_CHECK_LE( LlmKnowledgeBase::skillToolNameFromCanonical( skill.canonical_name ).size(), 64 );

      // The overview call is made for every question, so keep it small (body + reference index)
      const json index = LlmKnowledgeBase::executeSkillTool( skill, json::object() );
      BOOST_CHECK_LE( index.dump().size(), 16*1024 );

      for( const LlmKnowledgeBase::SkillReference &ref : skill.references )
      {
        BOOST_TEST_CONTEXT( "Reference '" << ref.name << "'" )
        {
          // An empty or overly long summary usually means a malformed frontmatter block
          BOOST_CHECK( !ref.summary.empty() );
          BOOST_CHECK_LE( ref.summary.size(), 150 );

          const json content = LlmKnowledgeBase::executeSkillTool( skill, json{ {"reference", ref.name} } );
          const string body = content["body"].get<string>();
          BOOST_CHECK( !body.empty() );
          BOOST_CHECK( body.find( "[... reference truncated ...]" ) == string::npos );
        }
      }//for( const SkillReference &ref : skill.references )
    }
  }//for( const SkillInfo &skill : skills )
}//BOOST_AUTO_TEST_CASE( shipped_skills_well_formed )


BOOST_AUTO_TEST_CASE( skill_filter )
{
  const string datadir = find_data_dir();
  BOOST_REQUIRE_MESSAGE( !datadir.empty(), "Could not find the static data directory; pass --datadir=..." );
  InterSpec::setStaticDataDirectory( datadir );

  BOOST_CHECK( !LlmKnowledgeBase::skillFilter().has_value() );
  const vector<LlmKnowledgeBase::SkillInfo> all = LlmKnowledgeBase::loadSkills();
  if( all.size() < 2 )
  {
    BOOST_TEST_MESSAGE( "Fewer than two skills found; not testing the skill filter." );
    return;
  }

  // A subset: only the named skills, by canonical name.
  const string wanted = all[1].canonical_name;
  LlmKnowledgeBase::setSkillFilter( set<string>{ wanted, "no-such-skill" } );
  const vector<LlmKnowledgeBase::SkillInfo> subset = LlmKnowledgeBase::loadSkills();
  BOOST_REQUIRE_EQUAL( subset.size(), 1 );
  BOOST_CHECK_EQUAL( subset[0].canonical_name, wanted );

  // None.
  LlmKnowledgeBase::setSkillFilter( set<string>{} );
  BOOST_CHECK( LlmKnowledgeBase::loadSkills().empty() );

  // Back to all.
  LlmKnowledgeBase::setSkillFilter( std::nullopt );
  BOOST_CHECK_EQUAL( LlmKnowledgeBase::loadSkills().size(), all.size() );
}//BOOST_AUTO_TEST_CASE( skill_filter )
