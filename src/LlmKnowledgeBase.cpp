#include "InterSpec_config.h"
#include "InterSpec/LlmKnowledgeBase.h"

#if( USE_LLM_INTERFACE )

#include <map>
#include <mutex>
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
/** Defaults for the remote RAG service, used when `<KnowledgeBaseRagModel>` or
 `<KnowledgeBaseRagCorpora>` are not given in llm_config.xml.  The service has renamed its corpus more
 than once, which is why these are overridable from the config rather than only being defined here.
 */
const char * const sm_default_rag_model = "openai/gpt-oss-120b";
const char * const sm_default_rag_corpus = "ldrd_llm_gamma_spec_labeled";


std::string query_rag_endpoint( const std::string &queryText,
                                          const std::string &rag_query_url,
                                          const std::string &model,
                                          const std::vector<std::string> &corpora )
{
  try
  {
    // Construct the JSON payload; the endpoint requires all three of these fields, and validates
    //  `corpora` against a fixed list of names, rejecting the request with a HTTP 422 otherwise.
    nlohmann::json payload;
    payload["query"] = queryText;
    payload["model"] = model;
    payload["corpora"] = nlohmann::json::array();
    for( const std::string &corpus : corpora )
      payload["corpora"].push_back( corpus );

    const std::string jsonStr = payload.dump();

    // Setup the curl command arguments.  We append the HTTP status code so a failed request
    //  reports the actual server error, rather than us trying to parse the error body as an answer.
    //  The service does multiple LLM round-trips per query, so it needs a generous timeout.
    const std::vector<std::string> args = {
      "-s",
      "-X", "POST",
      rag_query_url,
      "-H", "accept: application/json",
      "-H", "Content-Type: application/json",
      "--max-time", "600",
      "-w", "\n%{http_code}",
      "-d", jsonStr
    };

    // Execute curl and capture stdout into a stream
    boost::process::ipstream is;
    boost::process::child c(boost::process::search_path("curl"), args, boost::process::std_out > is);

    // Keep the lines seperate so we can pull the status code `-w` appended off the end; the body
    //  itself is JSON, so re-joining it with newlines does not change its meaning.
    std::vector<std::string> lines;
    std::string line;
    while( std::getline(is, line) )
    {
      if( !line.empty() && (line.back() == '\r') )
        line.pop_back();
      lines.push_back( line );
    }

    c.wait(); // Ensure the process finishes

    if( (c.exit_code() == 6) || (c.exit_code() == 7) )
      throw runtime_error( "The RAG service is down - this probably means this tool-call will not be useful in the future either." );

    if( c.exit_code() != 0 )
      throw std::runtime_error( "curl command failed with exit code " + std::to_string(c.exit_code()) );

    // Last line is the HTTP status code from `-w`; everything before it is the response body.
    int http_status = 0;
    if( !lines.empty() )
    {
      if( !SpecUtils::parse_int(lines.back().c_str(), lines.back().size(), http_status) )
        throw std::runtime_error( "Could not parse HTTP status code from RAG endpoint response." );
      lines.pop_back();
    }

    std::string output;
    for( size_t i = 0; i < lines.size(); ++i )
      output += (i ? "\n" : "") + lines[i];

    // A non-2xx means the request itself was rejected (e.g. a HTTP 422 when the endpoint changes
    //  the fields it expects), so surface the server's explanation instead of hunting for an answer.
    if( (http_status < 200) || (http_status > 299) )
      throw std::runtime_error( "RAG endpoint returned HTTP " + std::to_string(http_status)
                                + " for payload " + jsonStr + "; response: " + output );

    // Parse the result
    nlohmann::json result = nlohmann::json::parse(output);

    if( !result.contains("response") )
      throw std::runtime_error( "'response' key not found in JSON output from RAG endpoint." );

    return result["response"].get<std::string>();
  }catch( const nlohmann::json::parse_error &e )
  {
    throw std::runtime_error( std::string("JSON Parse Error: ") + e.what() );
  }catch( const std::exception &e )
  {
    throw std::runtime_error( std::string("Error: ") + e.what() );
  }
}//std::string query_rag_endpoint(...)


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
                           LlmKnowledgeBase::SkillInfo &skill,
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


/** Returns a one-line summary for a skill reference file: its frontmatter `description` if present,
 otherwise the text of its first Markdown heading, otherwise empty.
 */
std::string reference_summary( const std::string &content )
{
  std::istringstream input( content );
  std::string line;
  bool in_frontmatter = false, first_line = true;
  while( std::getline( input, line ) )
  {
    line = trim_copy( line );
    if( first_line && ((line == "---") || (line == "\xEF\xBB\xBF---")) ) //allow for UTF-8 BOM
    {
      in_frontmatter = true;
      first_line = false;
      continue;
    }
    first_line = false;

    if( in_frontmatter )
    {
      if( line == "---" )
      {
        in_frontmatter = false;
        continue;
      }

      if( SpecUtils::istarts_with( line, "description:" ) )
      {
        std::string value = trim_copy( line.substr( 12 ) );
        if( (value.size() >= 2) && ((value.front() == '"') || (value.front() == '\''))
           && (value.back() == value.front()) )
          value = value.substr( 1, value.size() - 2 );
        if( !value.empty() )
          return value;
      }
      continue;
    }//if( in_frontmatter )

    if( !line.empty() && (line[0] == '#') )
      return trim_copy( line.substr( line.find_first_not_of( '#' ) ) );
  }//while( std::getline( input, line ) )

  return "";
}//std::string reference_summary(...)


/** Finds the `references/*.md` files in a skill's directory. */
std::vector<LlmKnowledgeBase::SkillReference> find_skill_references( const std::string &skill_md_path )
{
  std::vector<LlmKnowledgeBase::SkillReference> refs;

  const std::string ref_dir = SpecUtils::append_path( SpecUtils::parent_path( skill_md_path ), "references" );
  if( !SpecUtils::is_directory( ref_dir ) )
    return refs;

  std::vector<std::string> files = SpecUtils::ls_files_in_directory( ref_dir, ".md" );
  std::sort( begin(files), end(files) );

  for( const std::string &path : files )
  {
    std::string content;
    if( !read_text_file( path, content ) )
    {
      cerr << "Warning: Failed reading KnowledgeBase skill reference '" << path << "'" << endl;
      continue;
    }

    LlmKnowledgeBase::SkillReference ref;
    ref.name = SpecUtils::filename( path );
    if( SpecUtils::iends_with( ref.name, ".md" ) )
      ref.name = ref.name.substr( 0, ref.name.size() - 3 );
    ref.summary = reference_summary( content );
    ref.file_path = path;

    if( !ref.name.empty() )
      refs.push_back( ref );
  }//for( const std::string &path : files )

  return refs;
}//find_skill_references(...)


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


namespace LlmKnowledgeBase
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

  return "knowledge_" + canonical_name;
}//std::string skillToolNameFromCanonical(...)


std::vector<SkillInfo> loadSkillsFromKnowledgeDir( const std::string &knowledge_dir )
{
  std::vector<SkillInfo> answer;
  if( knowledge_dir.empty() || !SpecUtils::is_directory(knowledge_dir) )
    return answer;

  std::vector<std::string> candidate_files = SpecUtils::recursive_ls( knowledge_dir, ".md" );
  std::sort( begin(candidate_files), end(candidate_files) );

  for( const std::string &path : candidate_files )
  {
    if( SpecUtils::filename(path) != "SKILL.md" )
      continue;

    std::string content;
    if( !read_text_file(path, content) )
    {
      cerr << "Warning: Failed reading KnowledgeBase SKILL.md '" << path << "'" << endl;
      continue;
    }

    SkillInfo skill;
    std::string warning;
    if( !parse_skill_markdown(content, path, skill, warning) )
    {
      cerr << "Warning: Failed parsing KnowledgeBase SKILL.md '" << path << "': " << warning << endl;
      continue;
    }

    if( !warning.empty() )
      cerr << "Warning: KnowledgeBase SKILL.md '" << path << "': " << warning << endl;

    const auto pos = std::find_if( begin(answer), end(answer), [&skill]( const SkillInfo &s ){
      return s.canonical_name == skill.canonical_name;
    } );
    if( pos != end(answer) )
    {
      cerr << "Warning: Duplicate KnowledgeBase skill '" << skill.canonical_name << "' at '" << path
           << "' ignored in favor of '" << pos->file_path << "'" << endl;
      continue;
    }

    skill.references = find_skill_references( path );
    answer.push_back( std::move(skill) );
  }//for( const std::string &path : candidate_files )

  return answer;
}//std::vector<SkillInfo> loadSkillsFromKnowledgeDir( const std::string &knowledge_dir )


std::vector<std::string> knowledgeDirectories()
{
  // Within each data directory, look in "llm_knowledge" and in the "llm_knowledge" directory of a
  //  checked-out knowledge repository (e.g. a git clone placed directly in the data directory).
  const char * const sub_dirs[] = { "llm_knowledge", "gamma_spec_skills/llm_knowledge" };

  std::vector<std::string> data_dirs;
  data_dirs.push_back( InterSpec::staticDataDirectory() );
  try
  {
    data_dirs.push_back( InterSpec::writableDataDirectory() );
  }catch( std::exception & )
  {
    // Writable directory not set (e.g., batch mode or some unit tests)
  }

  // A deployment's mounted config directory (--llm-config-dir) comes last, so it takes precedence.
  const std::string config_dir = LlmConfig::configDirectory();
  if( !config_dir.empty() )
    data_dirs.push_back( config_dir );

  std::vector<std::string> answer;
  for( const std::string &data_dir : data_dirs )
  {
    if( data_dir.empty() )
      continue;

    for( const char * const sub_dir : sub_dirs )
    {
      const std::string dir = SpecUtils::lexically_normalize_path( SpecUtils::append_path( data_dir, sub_dir ) );
      if( SpecUtils::is_directory( dir ) && (std::find( begin(answer), end(answer), dir ) == end(answer)) )
        answer.push_back( dir );
    }
  }//for( const std::string &data_dir : data_dirs )

  return answer;
}//std::vector<std::string> knowledgeDirectories()


namespace
{
  std::mutex ns_skill_filter_mutex;
  std::optional<std::set<std::string>> ns_skill_filter;
}//namespace


void setSkillFilter( std::optional<std::set<std::string>> canonical_names )
{
  std::lock_guard<std::mutex> lock( ns_skill_filter_mutex );
  ns_skill_filter = std::move( canonical_names );
}


std::optional<std::set<std::string>> skillFilter()
{
  std::lock_guard<std::mutex> lock( ns_skill_filter_mutex );
  return ns_skill_filter;
}


std::vector<SkillInfo> loadSkills()
{
  const std::optional<std::set<std::string>> filter = skillFilter();
  if( filter && filter->empty() )
    return {};

  std::map<std::string, SkillInfo> skills_by_name;

  // Directories are in increasing precedence, so a later skill with the same name replaces an earlier one.
  for( const std::string &dir : knowledgeDirectories() )
  {
    for( SkillInfo &skill : loadSkillsFromKnowledgeDir( dir ) )
    {
      const auto pos = skills_by_name.find( skill.canonical_name );
      if( pos != end(skills_by_name) )
        cerr << "Info: KnowledgeBase skill '" << skill.canonical_name << "' from '" << skill.file_path
             << "' overriding '" << pos->second.file_path << "'" << endl;
      skills_by_name[skill.canonical_name] = std::move(skill);
    }
  }//for( const std::string &dir : knowledgeDirectories() )

  std::vector<SkillInfo> loaded_skills;
  loaded_skills.reserve( skills_by_name.size() );
  for( const auto &entry : skills_by_name )
  {
    if( !filter || filter->count( entry.first ) )
      loaded_skills.push_back( entry.second );
  }

  return loaded_skills;
}//std::vector<SkillInfo> loadSkills()


std::string invokeDescriptionSuffix( const std::vector<SkillInfo> &skills, const bool has_rag )
{
  // Only the skill names are listed, to keep this description (which every calling agent sees on every
  //  turn) short; the full skill descriptions are on the `knowledge_<name>` tools, which only the
  //  KnowledgeBase agent sees, and where the actual skill selection happens.
  std::string answer;
  if( !skills.empty() )
  {
    answer += "\n\nIt has curated reference knowledge on these topics: ";
    for( size_t i = 0; i < skills.size(); ++i )
      answer += (i ? ", " : "") + skills[i].skill_name;
    answer += ".";
  }

  answer += "\n\nIt can also run JavaScript calculations (e.g., decay, unit, or dose arithmetic) to apply this knowledge.";

  if( has_rag )
    answer += "\nIt can also search (RAG) over a corpus of reference documents for topics not covered above.";

  return answer;
}//std::string invokeDescriptionSuffix(...)


nlohmann::json skillToolParameterSchema( const SkillInfo &skill )
{
  nlohmann::json schema = { {"type", "object"}, {"properties", nlohmann::json::object()} };

  if( !skill.references.empty() )
  {
    nlohmann::json names = nlohmann::json::array();
    for( const SkillReference &ref : skill.references )
      names.push_back( ref.name );

    schema["properties"]["reference"] = {
      {"type", "string"},
      {"description", "Optional. Name of one of this skill's reference documents to load, instead of the main"
                      " skill body. Call first without this parameter to see what each reference covers."},
      {"enum", names}
    };
  }//if( !skill.references.empty() )

  return schema;
}//nlohmann::json skillToolParameterSchema( const SkillInfo &skill )


nlohmann::json executeSkillTool( const SkillInfo &skill, const nlohmann::json &params )
{
  // Generous cap, to keep a single (mis-authored) reference from blowing out a small models context
  const size_t max_reference_chars = 48*1024;

  std::string requested;
  if( params.is_object() && params.contains("reference") && params["reference"].is_string() )
    requested = trim_copy( params["reference"].get<std::string>() );

  nlohmann::json answer;
  answer["skill_name"] = skill.skill_name;

  if( requested.empty() )
  {
    answer["description"] = skill.description;
    answer["source_file"] = skill.file_path;
    answer["body"] = skill.body;

    if( !skill.references.empty() )
    {
      nlohmann::json refs = nlohmann::json::array();
      for( const SkillReference &ref : skill.references )
        refs.push_back( { {"name", ref.name}, {"summary", ref.summary} } );
      answer["available_references"] = refs;
      answer["note"] = "Call this tool again with {\"reference\": \"<name>\"} to load only the references relevant to the question.";
    }

    return answer;
  }//if( requested.empty() )

  // Only exact (case-insensitive) matches to known references are allowed - never a raw path.
  const auto pos = std::find_if( begin(skill.references), end(skill.references),
    [&requested]( const SkillReference &ref ){ return SpecUtils::iequals_ascii( ref.name, requested ); } );

  if( pos == end(skill.references) )
  {
    std::string valid;
    for( const SkillReference &ref : skill.references )
      valid += (valid.empty() ? "" : ", ") + ref.name;
    throw std::runtime_error( "Unknown reference '" + requested + "' for skill '" + skill.skill_name
                              + "'; valid references are: " + (valid.empty() ? std::string("(none)") : valid) );
  }

  std::string content;
  if( !read_text_file( pos->file_path, content ) )
    throw std::runtime_error( "Failed to read reference '" + pos->name + "' of skill '" + skill.skill_name + "'." );

  if( content.size() > max_reference_chars )
  {
    content = content.substr( 0, max_reference_chars );
    SpecUtils::utf8_limit_str_size( content, max_reference_chars );
    content += "\n\n[... reference truncated ...]";
  }

  answer["reference"] = pos->name;
  answer["summary"] = pos->summary;
  answer["source_file"] = pos->file_path;
  answer["body"] = content;

  return answer;
}//nlohmann::json executeSkillTool(...)


void registerKnowledgeBaseTools( const std::vector<SkillInfo> &skills,
                                 const std::string &rag_url,
                                 const std::string &rag_model,
                                 const std::vector<std::string> &rag_corpora,
                                 const std::function<void(const LlmTools::SharedTool&)> &register_tool )
{
  if( !register_tool )
    throw std::invalid_argument( "registerKnowledgeBaseTools requires a valid register callback." );

  if( !rag_url.empty() )
  {
    LlmTools::SharedTool query_tool;
    query_tool.name = "query_rag_endpoint";
    query_tool.description = "Query the optional retrieval-augmented-generation (RAG) HTTP endpoint (searches a corpus of reference documents) for supplemental domain information.";
    query_tool.parameters_schema = nlohmann::json::parse(R"({
      "type": "object",
      "properties": {
        "question": {
          "type": "string",
          "description": "Natural-language research question to send to the remote RAG endpoint."
        }
      },
      "required": ["question"]
    })");
    query_tool.availableForAgents = {AgentType::KnowledgeBase};
    query_tool.executor = [rag_url, rag_model, rag_corpora](
                                               const nlohmann::json& params,
                                               InterSpec* interspec,
                                               std::shared_ptr<LlmInteraction>,
                                               LlmConversationHistory* ) -> nlohmann::json {
      (void)interspec;
      return executeQueryRagEndpoint( params, rag_url,
                                               rag_model, rag_corpora );
    };

    register_tool( query_tool );
  }//if( !rag_url.empty() )

  for( const SkillInfo &skill : skills )
  {
    LlmTools::SharedTool skill_tool;
    skill_tool.name = skillToolNameFromCanonical( skill.canonical_name );
    skill_tool.description = "Load knowledge skill '" + skill.skill_name + "': " + skill.description;
    if( !skill.references.empty() )
      skill_tool.description += " (Has " + std::to_string(skill.references.size())
                                + " optional detailed references, loadable via the 'reference' parameter.)";
    skill_tool.parameters_schema = skillToolParameterSchema( skill );
    skill_tool.availableForAgents = {AgentType::KnowledgeBase};
    skill_tool.executor = [skill]( const nlohmann::json& params,
                                   InterSpec* interspec,
                                   std::shared_ptr<LlmInteraction>,
                                   LlmConversationHistory* ) -> nlohmann::json {
      (void)interspec;
      return executeSkillTool( skill, params );
    };

    register_tool( skill_tool );
  }//for( const SkillInfo &skill : skills )
}//void registerKnowledgeBaseTools(...)


nlohmann::json executeQueryRagEndpoint( const nlohmann::json& params,
                                                 const std::string &rag_url,
                                                 const std::string &rag_model,
                                                 const std::vector<std::string> &rag_corpora )
{
  const std::string question_key = find_case_insensitive_key( "question", params );
  if( !params.contains(question_key) || !params[question_key].is_string() )
  {
    cerr << "executeQueryRagEndpoint error: The 'question' parameter must be present, and a string; params: " << params.dump() << endl << endl;
    throw runtime_error( "The 'question' parameter must be present, and a string." );
  }

  const std::string question = params[question_key];

  cout << "In executeQueryRagEndpoint; question: " << question << endl << endl;

  const std::string model = rag_model.empty() ? sm_default_rag_model
                                                        : rag_model;
  const std::vector<std::string> corpora = rag_corpora.empty()
                                     ? std::vector<std::string>{ sm_default_rag_corpus }
                                     : rag_corpora;

  nlohmann::json result;

  try
  {
    const std::string query_result = query_rag_endpoint( question, rag_url, model, corpora );
    cout << "executeQueryRagEndpoint answer: " << query_result << endl << endl;
    result["success"] = true;
    result["answer"] = query_result;
  }catch( std::exception &e )
  {
    result["success"] = false;
    result["error"] = e.what();
  }

  return result;
}//nlohmann::json executeQueryRagEndpoint(...)
}//namespace LlmKnowledgeBase

#endif //#if( USE_LLM_INTERFACE )
