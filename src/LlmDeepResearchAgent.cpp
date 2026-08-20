#include "InterSpec_config.h"
#include "InterSpec/LlmDeepResearchAgent.h"

#if( USE_LLM_INTERFACE )

#include <map>
#include <cctype>
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <boost/process.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmToolRegistry.h"

using namespace std;

namespace
{
std::string query_deep_research_endpoint( const std::string &queryText, const std::string &deep_query_url )
{
  try
  {
    // Construct the JSON payload
    nlohmann::json payload = {
      {"messages", {{{"content", queryText}}}},
      {"corpora", {"ldrd_llm_gamma_spec"}}
    };
    std::string jsonStr = payload.dump();

    // Setup the curl command arguments
    std::vector<std::string> args = {
      "-X", "POST",
      deep_query_url,
      "-H", "accept: application/json",
      "-H", "Content-Type: application/json",
      "-d", jsonStr
    };

    // Execute curl and capture stdout into a stream
    boost::process::ipstream is;
    boost::process::child c(boost::process::search_path("curl"), args, boost::process::std_out > is);

    std::string output;
    std::string line;
    while( std::getline(is, line) )
      output += line;

    c.wait(); // Ensure the process finishes

    if( (c.exit_code() == 6) || (c.exit_code() == 7) )
      throw runtime_error( "The deep-research service is down - this probably means this tool-call will not be useful in the future either." );

    if( c.exit_code() != 0 )
      throw std::runtime_error( "curl command failed with exit code " + std::to_string(c.exit_code()) );

    // Parse the result
    nlohmann::json result = nlohmann::json::parse(output);

    if( !result.contains("response") )
      throw std::runtime_error( "'response' key not found in JSON output from deep-research endpoint." );

    return result["response"].get<std::string>();
  }catch( const nlohmann::json::parse_error &e )
  {
    throw std::runtime_error( std::string("JSON Parse Error: ") + e.what() );
  }catch( const std::exception &e )
  {
    throw std::runtime_error( std::string("Error: ") + e.what() );
  }
}//std::string query_deep_research_endpoint(...)


std::string trim_copy( std::string value )
{
  SpecUtils::trim( value );
  return value;
}


bool read_text_file( const std::string &path, std::string &contents )
{
#ifdef _WIN32
  const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16( path );
  std::ifstream input( wpath.c_str(), std::ios_base::binary | std::ios_base::in );
#else
  std::ifstream input( path.c_str(), std::ios_base::binary | std::ios_base::in );
#endif

  if( !input.is_open() )
    return false;

  std::stringstream buffer;
  buffer << input.rdbuf();
  contents = buffer.str();
  return true;
}//bool read_text_file(...)


std::string canonicalize_skill_name( std::string value )
{
  SpecUtils::to_lower_ascii( value );

  std::string canonical;
  canonical.reserve( value.size() );

  bool prev_was_dash = false;
  for( char c : value )
  {
    if( std::isalnum( static_cast<unsigned char>(c) ) )
    {
      canonical.push_back( c );
      prev_was_dash = false;
    }else if( !prev_was_dash )
    {
      canonical.push_back( '-' );
      prev_was_dash = true;
    }
  }//for( char c : value )

  while( !canonical.empty() && (canonical.front() == '-') )
    canonical.erase( canonical.begin() );

  while( !canonical.empty() && (canonical.back() == '-') )
    canonical.pop_back();

  return canonical;
}//std::string canonicalize_skill_name(...)


bool parse_skill_markdown( const std::string &content,
                           const std::string &path,
                           LlmDeepResearch::SkillInfo &skill,
                           std::string &warning )
{
  warning.clear();

  std::string text = content;
  if( (text.size() >= 3)
     && (static_cast<unsigned char>(text[0]) == 0xEF)
     && (static_cast<unsigned char>(text[1]) == 0xBB)
     && (static_cast<unsigned char>(text[2]) == 0xBF) )
  {
    text.erase( 0, 3 );
  }

  std::vector<std::string> lines;
  {
    std::istringstream input( text );
    std::string line;
    while( std::getline( input, line ) )
    {
      if( !line.empty() && (line.back() == '\r') )
        line.pop_back();
      lines.push_back( line );
    }
  }

  const std::string parent_dir_name = SpecUtils::filename( SpecUtils::parent_path(path) );
  std::string skill_name = parent_dir_name.empty() ? "unnamed-skill" : parent_dir_name;
  std::string description;
  std::string body = text;

  if( !lines.empty() && (trim_copy(lines[0]) == "---") )
  {
    size_t frontmatter_end = std::numeric_limits<size_t>::max();
    for( size_t i = 1; i < lines.size(); ++i )
    {
      if( trim_copy(lines[i]) == "---" )
      {
        frontmatter_end = i;
        break;
      }
    }

    if( frontmatter_end == std::numeric_limits<size_t>::max() )
    {
      warning = "frontmatter starts with '---' but missing closing delimiter";
    }else
    {
      for( size_t i = 1; i < frontmatter_end; ++i )
      {
        std::string line = trim_copy( lines[i] );
        if( line.empty() || (line[0] == '#') )
          continue;

        const size_t colon_pos = line.find( ':' );
        if( colon_pos == std::string::npos )
          continue;

        std::string key = trim_copy( line.substr(0, colon_pos) );
        std::string value = trim_copy( line.substr(colon_pos + 1) );
        SpecUtils::to_lower_ascii( key );

        if( (value.size() >= 2)
           && (((value.front() == '"') && (value.back() == '"'))
               || ((value.front() == '\'') && (value.back() == '\''))) )
        {
          value = value.substr( 1, value.size() - 2 );
        }

        if( key == "name" )
        {
          if( !value.empty() )
            skill_name = value;
        }else if( key == "description" )
        {
          description = value;
        }
      }//for( lines in frontmatter )

      std::ostringstream body_stream;
      for( size_t i = frontmatter_end + 1; i < lines.size(); ++i )
      {
        if( i > (frontmatter_end + 1) )
          body_stream << '\n';
        body_stream << lines[i];
      }
      body = body_stream.str();
    }
  }//if( frontmatter starts with --- )

  body = trim_copy( body );
  if( body.empty() )
  {
    warning = "SKILL.md body was empty";
    return false;
  }

  if( description.empty() )
    description = "Loaded from " + path;

  std::string canonical_name = canonicalize_skill_name( skill_name );
  if( canonical_name.empty() )
    canonical_name = canonicalize_skill_name( parent_dir_name );
  if( canonical_name.empty() )
    canonical_name = "unnamed-skill";

  skill.canonical_name = canonical_name;
  skill.skill_name = skill_name;
  skill.description = description;
  skill.body = body;
  skill.file_path = path;

  return true;
}//bool parse_skill_markdown(...)


std::string find_case_insensitive_key( const std::string &desired_key, const nlohmann::json &obj )
{
  std::string desired_lower = desired_key;
  SpecUtils::to_lower_ascii( desired_lower );

  for( auto it = obj.begin(); it != obj.end(); ++it )
  {
    std::string key = it.key();
    SpecUtils::to_lower_ascii( key );
    if( key == desired_lower )
      return it.key();
  }

  return desired_key;
}//std::string find_case_insensitive_key(...)
}//namespace


namespace LlmDeepResearch
{
std::string skillToolNameFromCanonical( std::string canonical_name )
{
  for( char &c : canonical_name )
  {
    if( c == '-' )
      c = '_';
  }

  if( canonical_name.empty() )
    canonical_name = "unnamed_skill";

  return "deepresearch_skill_" + canonical_name;
}//std::string skillToolNameFromCanonical(...)


std::vector<SkillInfo> loadSkills()
{
  std::map<std::string, SkillInfo> skills_by_name;

  auto load_skills_from_data_dir = [&skills_by_name]( const std::string &data_dir, const bool overwrite_existing )
  {
    if( data_dir.empty() )
      return;

    const std::string knowledge_dir = SpecUtils::append_path( data_dir, "llm_knowledge" );
    if( !SpecUtils::is_directory(knowledge_dir) )
      return;

    std::vector<std::string> candidate_files = SpecUtils::recursive_ls( knowledge_dir, ".md" );
    std::sort( begin(candidate_files), end(candidate_files) );

    for( const std::string &path : candidate_files )
    {
      if( SpecUtils::filename(path) != "SKILL.md" )
        continue;

      std::string content;
      if( !read_text_file(path, content) )
      {
        cerr << "Warning: Failed reading DeepResearch SKILL.md '" << path << "'" << endl;
        continue;
      }

      SkillInfo skill;
      std::string warning;
      if( !parse_skill_markdown(content, path, skill, warning) )
      {
        cerr << "Warning: Failed parsing DeepResearch SKILL.md '" << path << "': " << warning << endl;
        continue;
      }

      if( !warning.empty() )
        cerr << "Warning: DeepResearch SKILL.md '" << path << "': " << warning << endl;

      const auto pos = skills_by_name.find( skill.canonical_name );
      if( pos != end(skills_by_name) )
      {
        if( !overwrite_existing )
        {
          cerr << "Warning: Duplicate DeepResearch skill '" << skill.canonical_name << "' at '" << path
               << "' ignored in favor of '" << pos->second.file_path << "'" << endl;
          continue;
        }

        cerr << "Info: DeepResearch skill '" << skill.canonical_name << "' from '" << path
             << "' overriding '" << pos->second.file_path << "'" << endl;
      }

      skills_by_name[skill.canonical_name] = skill;
    }//for( const std::string &path : candidate_files )
  };//load_skills_from_data_dir(...)

  load_skills_from_data_dir( InterSpec::staticDataDirectory(), false );

  try
  {
    load_skills_from_data_dir( InterSpec::writableDataDirectory(), true );
  }catch( std::exception &e )
  {
    cerr << "Warning: Could not scan writable data directory for DeepResearch skills: " << e.what() << endl;
  }

  std::vector<SkillInfo> loaded_skills;
  loaded_skills.reserve( skills_by_name.size() );
  for( const auto &entry : skills_by_name )
    loaded_skills.push_back( entry.second );

  return loaded_skills;
}//std::vector<SkillInfo> loadSkills()


void registerDeepResearchTools( const std::string &deep_research_url,
                                const std::function<void(const LlmTools::SharedTool&)> &register_tool,
                                bool &loaded_any_skill )
{
  if( !register_tool )
    throw std::invalid_argument( "registerDeepResearchTools requires a valid register callback." );

  loaded_any_skill = false;

  if( !deep_research_url.empty() )
  {
    LlmTools::SharedTool query_tool;
    query_tool.name = "query_deep_research_endpoint";
    query_tool.description = "Query the optional deep research HTTP endpoint for supplemental domain information.";
    query_tool.parameters_schema = nlohmann::json::parse(R"({
      "type": "object",
      "properties": {
        "question": {
          "type": "string",
          "description": "Natural-language research question to send to the remote deep research endpoint."
        }
      },
      "required": ["question"]
    })");
    query_tool.availableForAgents = {AgentType::DeepResearch};
    query_tool.executor = [deep_research_url]( const nlohmann::json& params,
                                               InterSpec* interspec,
                                               std::shared_ptr<LlmInteraction>,
                                               LlmConversationHistory* ) -> nlohmann::json {
      (void)interspec;
      return executeQueryDeepResearchEndpoint( params, deep_research_url );
    };

    register_tool( query_tool );
  }//if( !deep_research_url.empty() )

  const std::vector<SkillInfo> deep_research_skills = loadSkills();
  loaded_any_skill = !deep_research_skills.empty();
  for( const SkillInfo &skill : deep_research_skills )
  {
    LlmTools::SharedTool skill_tool;
    skill_tool.name = skillToolNameFromCanonical( skill.canonical_name );
    skill_tool.description = "Load DeepResearch skill '" + skill.skill_name + "': " + skill.description;
    skill_tool.parameters_schema = nlohmann::json::parse(R"({
      "type": "object",
      "properties": {}
    })");
    skill_tool.availableForAgents = {AgentType::DeepResearch};
    skill_tool.executor = [skill]( const nlohmann::json& params,
                                   InterSpec* interspec,
                                   std::shared_ptr<LlmInteraction>,
                                   LlmConversationHistory* ) -> nlohmann::json {
      (void)params;
      (void)interspec;

      nlohmann::json answer;
      answer["skill_name"] = skill.skill_name;
      answer["description"] = skill.description;
      answer["source_file"] = skill.file_path;
      answer["body"] = skill.body;
      return answer;
    };

    register_tool( skill_tool );
  }//for( const SkillInfo &skill : deep_research_skills )
}//void registerDeepResearchTools(...)


nlohmann::json executeQueryDeepResearchEndpoint( const nlohmann::json& params, const std::string &deep_research_url )
{
  const std::string question_key = find_case_insensitive_key( "question", params );
  if( !params.contains(question_key) || !params[question_key].is_string() )
  {
    cerr << "executeQueryDeepResearchEndpoint error: The 'question' parameter must be present, and a string; params: " << params.dump() << endl << endl;
    throw runtime_error( "The 'question' parameter must be present, and a string." );
  }

  const std::string question = params[question_key];

  cout << "In executeQueryDeepResearchEndpoint; question: " << question << endl << endl;

  nlohmann::json result;

  try
  {
    const std::string query_result = query_deep_research_endpoint( question, deep_research_url );
    cout << "executeQueryDeepResearchEndpoint answer: " << query_result << endl << endl;
    result["success"] = true;
    result["answer"] = query_result;
  }catch( std::exception &e )
  {
    result["success"] = false;
    result["error"] = e.what();
  }

  return result;
}//nlohmann::json executeQueryDeepResearchEndpoint(...)
}//namespace LlmDeepResearch

#endif //#if( USE_LLM_INTERFACE )
