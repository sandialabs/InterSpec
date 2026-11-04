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

/** The stdio MCP server: an MCP client launches this executable and exchanges newline-delimited
 JSON-RPC messages with it over stdin/stdout.  It hosts one InterSpec session, with no browser and
 no web server, in a `Wt::Test::WTestEnvironment`.  See README.md in this directory.

 Threads: the main thread reads stdin and hands each message to the `LlmMcpDispatcher`, which runs
 tools in the session through `WServer::schedule()` (on Wt's IO threads, which take the session
 lock); the main thread never holds the session lock while waiting.  Replies are written by
 whichever thread has them, one whole line at a time.

 With `--run-benchmark`, it is not an MCP server: it runs LLM benchmarks through the assistant
 (`LlmToolGui`/`LlmBenchmarkRunner`) in the session, the same way, prints a summary of each, and exits.
 */

#include "InterSpec_config.h"

#include <map>
#include <set>
#include <mutex>
#include <atomic>
#include <string>
#include <vector>
#include <cstdio>
#include <thread>
#include <memory>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <variant>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <functional>
#include <condition_variable>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#else
#include <unistd.h>
#endif

#include <boost/asio/post.hpp>
#include <boost/program_options.hpp>

#include <Wt/WServer.h>
#include <Wt/WLogger.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>
#include <Wt/Test/WTestEnvironment.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SerialToDetectorModel.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/LlmKnowledgeBase.h"
#include "InterSpec/LlmMcpDispatcher.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/LlmBenchmarkRunner.h"
#include "InterSpec/SpectrumChartImage.h"
#include "InterSpec/DecayDataBaseServer.h"

#include "QuickJsSandbox.h"
#include "PlutoPaintDevice.h"

#ifdef _WIN32
// After the Boost.Asio (via Wt) headers, which need <winsock2.h> before <windows.h>.
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>  // for TerminateProcess
#endif

static_assert( BUILD_AS_MCP_STDIO_SERVER, "This file is only for the stdio MCP server build" );

using namespace std;
using json = nlohmann::json;


namespace
{
  /** Keeps the original stdout for our own output (the MCP protocol, or a benchmark report), and
   sends everything else written to stdout (by InterSpec, Wt, or libraries) to stderr instead, so it
   can not get mixed into it.
   */
  class ReservedStdout
  {
  public:
    ReservedStdout()
    {
#ifdef _WIN32
      // A stream the process was started without has _fileno() -2; passing that to _setmode, _dup
      //  or _dup2 would end the process from the CRT's invalid-parameter handler.
      if( _fileno(stdout) < 0 )
        throw std::runtime_error( "stdout must be connected." );
      _setmode( _fileno(stdout), _O_BINARY );
      fflush( stdout );
      m_fd = _dup( _fileno(stdout) );

      // With no stderr, stray output goes nowhere.
      const int nul_fd = (_fileno(stderr) < 0) ? _open( "NUL", _O_WRONLY ) : -1;
      const int sink_fd = (_fileno(stderr) < 0) ? nul_fd : _fileno(stderr);
      const bool redirected = (sink_fd >= 0) && (_dup2( sink_fd, _fileno(stdout) ) == 0);
      if( nul_fd >= 0 )
        _close( nul_fd );
      if( (m_fd < 0) || !redirected )
        throw std::runtime_error( "Could not redirect stdout." );
#else
      fflush( stdout );
      m_fd = dup( STDOUT_FILENO );
      if( (m_fd < 0) || (dup2( STDERR_FILENO, STDOUT_FILENO ) < 0) )
        throw std::runtime_error( "Could not redirect stdout." );
#endif
    }//ReservedStdout()

    /** Writes all of `text` (in one piece, with respect to other calls); returns false if the reader
     has gone away. */
    bool write( const std::string &text )
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      size_t written = 0;
      while( written < text.size() )
      {
#ifdef _WIN32
        const int n = _write( m_fd, text.data() + written, static_cast<unsigned int>(text.size() - written) );
#else
        const ssize_t n = ::write( m_fd, text.data() + written, text.size() - written );
#endif
        if( n <= 0 )
          return false;
        written += static_cast<size_t>( n );
      }
      return true;
    }//write(...)

  private:
    int m_fd = -1;
    std::mutex m_mutex;
  };//class ReservedStdout


  /** Writes protocol messages, one per line, to the original stdout (see ReservedStdout). */
  class ProtocolWriter
  {
  public:
    ProtocolWriter()
    {
#ifdef _WIN32
      if( _fileno(stdin) < 0 )
        throw std::runtime_error( "stdin and stdout must be connected to the MCP client." );
      _setmode( _fileno(stdin), _O_BINARY );
#endif
    }

    /** Writes `message` as one line; returns false if the client has gone away. */
    bool send( const json &message )
    {
      return m_out.write( LlmMcpDispatcher::dumpMessage( message ) + "\n" );
    }

  private:
    ReservedStdout m_out;
  };//class ProtocolWriter


  /** `run_javascript` calls waiting for QuickJS, by id.  Their callbacks are only created, called,
   and destroyed in the session; only the id and the (plain JSON) result cross to the worker thread.
   */
  std::mutex ns_js_calls_mutex;
  uint64_t ns_next_js_call = 0;
  std::map<uint64_t, LlmTools::SharedTool::AsyncCallback> ns_js_calls;

  /** Runs `run_javascript` with QuickJS on a Wt IO thread (it may take up to 30 s), and gives the
   result back in the session.  Called in the session.
   */
  void run_javascript( const json &params, LlmTools::SharedTool::AsyncCallback callback )
  {
    Wt::WServer * const server = Wt::WServer::instance();
    Wt::WApplication * const app = Wt::WApplication::instance();
    if( !server || !app )
    {
      callback( std::string( "No session to run JavaScript in." ) );
      return;
    }

    uint64_t id = 0;
    {
      std::lock_guard<std::mutex> lock( ns_js_calls_mutex );
      id = ++ns_next_js_call;
      ns_js_calls[id] = std::move( callback );
    }

    const std::string session_id = app->sessionId();
    boost::asio::post( server->ioService(), [id, session_id, params](){
      std::variant<json, std::string> result = QuickJsSandbox::run( params );

      Wt::WServer * const server = Wt::WServer::instance();
      if( !server )
        return;

      server->schedule( std::chrono::milliseconds(1), session_id, [id, result](){
        LlmTools::SharedTool::AsyncCallback done;
        {
          std::lock_guard<std::mutex> lock( ns_js_calls_mutex );
          const auto pos = ns_js_calls.find( id );
          if( pos == end(ns_js_calls) )
            return;
          done = std::move( pos->second );
          ns_js_calls.erase( pos );
        }
        done( result );
      } );
    } );
  }//run_javascript(...)


  /** The first of `dirs` that exists (and contains `must_contain`, if given), or empty. */
  std::string first_existing_dir( const std::vector<std::string> &dirs, const std::string &must_contain = "" )
  {
    for( const std::string &dir : dirs )
    {
      if( dir.empty() || !SpecUtils::is_directory( dir ) )
        continue;
      if( must_contain.empty() || SpecUtils::is_file( SpecUtils::append_path( dir, must_contain ) ) )
        return dir;
    }
    return std::string();
  }//first_existing_dir(...)


  /** InterSpec's static data directory, if installed next to this executable: as for a wxWidgets
   or Electron install of InterSpec, a `make install`, or this build directory.
   */
  std::string default_static_data_dir()
  {
    std::string exe_dir;
    try
    {
      exe_dir = SpecUtils::parent_path( AppUtils::current_exe_path() );
    }catch( std::exception & )
    {
      return std::string();
    }

    return first_existing_dir( {
      SpecUtils::append_path( exe_dir, "data" ),
      SpecUtils::append_path( exe_dir, "resources/app/data" ),
      SpecUtils::append_path( exe_dir, "../data" ),
      SpecUtils::append_path( exe_dir, "../share/interspec/data" ),
      SpecUtils::append_path( exe_dir, "LibInterSpec/data" )
    }, "sandia.decay.xml" );
  }//default_static_data_dir()


  /** Sets up the per-user data InterSpec reads at startup, as `InterSpecServer::start_server` does
   for the GUI builds (the user database is in memory for this build, so is not involved).
   */
  void set_user_data_dir( const std::string &user_data_dir )
  {
    InterSpec::setWritableDataDirectory( user_data_dir );

    for( const std::string &dir : { user_data_dir, InterSpec::staticDataDirectory() } )
    {
      const std::vector<std::string> serial_db = SpecUtils::ls_files_in_directory( dir, "serial_to_model.csv" );
      if( !serial_db.empty() )
      {
        SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
        break;
      }
    }

    const std::string user_decay = SpecUtils::append_path( user_data_dir, "sandia.decay.xml" );
    if( SpecUtils::is_file( user_decay ) )
      DecayDataBaseServer::setDecayXmlFile( user_decay );

    const std::string user_reaction = SpecUtils::append_path( user_data_dir, "sandia.reactiongamma.xml" );
    if( SpecUtils::is_file( user_reaction ) )
      ReactionGammaServer::set_xml_file_location( user_reaction );
  }//set_user_data_dir(...)


  /** The `instructions` sent to the MCP client with its `initialize` reply: mcp_stdio_instructions.md,
   from the LLM config, user data, or static data directory (see LlmConfig::configFilePath); empty
   if there is none.
   */
  std::string load_mcp_instructions()
  {
    const std::string path = LlmConfig::configFilePath( "mcp_stdio_instructions.md" );
    try
    {
      std::vector<char> data;
      SpecUtils::load_file_data( path.c_str(), data );  // null-terminates `data`
      return SpecUtils::trim_copy( std::string( data.data() ) );
    }catch( std::exception & )
    {
      std::cerr << "Warning: no instructions for the MCP client; could not read '" << path << "'." << std::endl;
    }
    return std::string();
  }//load_mcp_instructions()


  /** Exit codes of `--run-benchmark` / `--list-benchmarks`. */
  const int sm_exit_ok = 0;              // Every benchmark ran to the end
  const int sm_exit_setup_error = 1;     // Bad arguments, or the setup failed
  const int sm_exit_benchmark_error = 2; // A benchmark could not be run (or was interrupted)


  /** Set by the SIGINT handler; the main thread then cancels the running benchmark. */
  std::atomic<bool> ns_interrupted{ false };

  void handle_sigint( int )
  {
    ns_interrupted = true;
    std::signal( SIGINT, SIG_DFL );  // A second Ctrl-C ends the program right away
  }


  /** What the main thread learns of a benchmark run, from the session. */
  struct BenchmarkRunState
  {
    std::mutex mutex;
    std::condition_variable cv;
    bool started = false;    // Set once the run was started in the session (or failed to be)
    bool startFailed = false;
    std::string startError;
    bool finished = false;
    bool completed = false;
    std::string resultsPath;
    BenchmarkResults results;
  };//struct BenchmarkRunState


  std::string percent_str( const int num, const int denom )
  {
    if( denom <= 0 )
      return std::string();
    std::ostringstream strm;
    strm << " (" << std::fixed << std::setprecision( 1 ) << (100.0 * num / denom) << "%)";
    return strm.str();
  }//percent_str(...)


  /** Plain-text summary of a benchmark run, for the report on stdout. */
  std::string benchmark_summary( const BenchmarkResults &results, const std::string &results_path,
                                 const bool completed )
  {
    std::ostringstream out;

    std::string skills;
    for( const std::string &skill : results.skills )
      skills += (skills.empty() ? "" : ",") + skill;

    const auto opt_str = []( const std::optional<size_t> &v ) -> std::string {
      return v.has_value() ? std::to_string( v.value() ) : std::string( "-" );
    };

    out << "Benchmark:    " << results.benchmarkName << (completed ? "" : "  [INCOMPLETE]") << "\n"
        << "File:         " << results.benchmarkFilePath << "\n"
        << "Model:        " << results.modelName << "\n"
        << "Judge model:  " << (results.judgeModelName.empty() ? std::string( "-" ) : results.judgeModelName) << "\n"
        << "Skills:       " << (skills.empty() ? std::string( "none" ) : skills) << "\n"
        << "Score:        " << results.totalScore << "/" << results.maxTotalScore
                            << percent_str( results.totalScore, results.maxTotalScore ) << "\n"
        << "Correct:      " << results.correctCount << "/" << results.gradedQuestions
                            << percent_str( results.correctCount, results.gradedQuestions ) << "\n"
        << "Errors:       " << results.errorCount << "\n"
        << "Questions:    " << results.totalQuestions << "\n";

    for( const BenchmarkQuestionResult &r : results.questionResults )
    {
      const bool skipped = SpecUtils::starts_with( r.llmAnswer, "[SKIPPED]" );
      std::string score = "-";
      if( r.graded && !r.hadError )
        score = std::to_string( r.score ) + "/" + std::to_string( r.maxScore );

      std::string marker;
      if( r.hadError )
        marker = "ERROR";
      else if( skipped )
        marker = "SKIP";
      else if( !r.graded )
        marker = "UNGRADED";
      else
        marker = r.correct ? "ok" : "wrong";

      const std::string id = r.problemId + " Q" + std::to_string( r.questionPart );
      out << "  " << std::left << std::setw( 26 ) << id << std::right
          << " " << std::setw( 7 ) << score
          << "  " << std::left << std::setw( 8 ) << marker << std::right
          << "  retries=" << r.errorRetries
          << "  tokens=" << opt_str( r.totalTokens ) << "\n";
    }//for( const BenchmarkQuestionResult &r : results.questionResults )

    out << "Results JSON: " << (results_path.empty() ? std::string( "-" ) : results_path) << "\n\n";

    return out.str();
  }//benchmark_summary(...)


  /** Parses `--skills`: "all" (no filter), "none", or comma-separated canonical skill names.
   @throws std::runtime_error, listing the available skill names, if a name is not known.
   */
  std::optional<std::set<std::string>> parse_skills_arg( const std::string &arg )
  {
    const std::string value = SpecUtils::trim_copy( arg );
    if( value.empty() || SpecUtils::iequals_ascii( value, "all" ) )
      return std::nullopt;
    if( SpecUtils::iequals_ascii( value, "none" ) )
      return std::set<std::string>{};

    std::set<std::string> available;
    for( const LlmKnowledgeBase::SkillInfo &skill : LlmKnowledgeBase::loadSkills() )
      available.insert( skill.canonical_name );

    std::vector<std::string> names;
    SpecUtils::split( names, value, "," );

    std::set<std::string> answer, unknown;
    for( std::string name : names )
    {
      SpecUtils::trim( name );
      if( name.empty() )
        continue;
      if( available.count( name ) )
        answer.insert( name );
      else
        unknown.insert( name );
    }

    if( !unknown.empty() || answer.empty() )
    {
      std::string msg = "Unknown skill name(s) for --skills:";
      for( const std::string &name : unknown )
        msg += " '" + name + "'";
      msg += ".  Available skills:";
      for( const std::string &name : available )
        msg += " " + name;
      if( available.empty() )
        msg += " (none found)";
      throw std::runtime_error( msg );
    }

    return answer;
  }//parse_skills_arg(...)
}//namespace


int main( int argc, char **argv )
{
#ifdef _WIN32
  AppUtils::getUtf8Args( argc, argv );
#else
  // A client that goes away must not kill us mid-write; the write fails instead.
  std::signal( SIGPIPE, SIG_IGN );
#endif

  std::string static_data_dir, user_data_dir, llm_config_dir, log_file;
  std::vector<std::string> benchmark_files;
  std::string skills_arg = "all", model_override, judge_model, benchmark_output_dir;

  namespace po = boost::program_options;
  po::options_description cl_desc( "Allowed options" );
  cl_desc.add_options()
  ("help,h",  "produce this help message")
  ("datadir", po::value<std::string>(&static_data_dir),
   "InterSpec's static data directory (nuclear data, detector responses, ...); by default, the"
   " 'data' directory of an InterSpec installed next to this executable."
   )
  ("userdatadir", po::value<std::string>(&user_data_dir),
   "InterSpec's per-user data directory (LLM configuration, knowledge base, detector responses,"
   " ...); by default, the one the installed InterSpec app uses, if it exists."
   )
  ("llm-config-dir", po::value<std::string>(&llm_config_dir),
   "Optional: a directory with llm_config.xml (and llm_tools_config.xml, ...), used before the user"
   " data and static data directories."
   )
  ("log-file", po::value<std::string>(&log_file),
   "Optional: append log output (otherwise written to stderr) to this file."
   )
  ;

  po::options_description bench_desc( "LLM benchmark options (instead of being an MCP server)" );
  bench_desc.add_options()
  ("run-benchmark", po::value<std::vector<std::string>>(&benchmark_files)->composing(),
   "Run this *_llm_benchmark.xml through the LLM assistant, print a summary of the results to"
   " stdout, and exit; may be given more than once (run in order, in one session).  Progress goes"
   " to stderr.  Exit code: 0 if every benchmark ran to the end (wrong answers and per-question"
   " errors included), 1 for bad arguments or setup, 2 if a benchmark could not be run."
   )
  ("skills", po::value<std::string>(&skills_arg),
   "With --run-benchmark: the KnowledgeBase skills the assistant gets: 'all' (default), 'none' (no"
   " KnowledgeBase agent, even with a RAG endpoint), or a comma-separated list of skill names."
   )
  ("model", po::value<std::string>(&model_override),
   "With --run-benchmark: the model (of the active provider in llm_config.xml) to test, for this"
   " run only; need not be listed in llm_config.xml.  Default: the active model."
   )
  ("judge-model", po::value<std::string>(&judge_model),
   "With --run-benchmark: the model (of the same provider) that grades the answers, for this run"
   " only.  Default: the <Model judge=\"true\"> in llm_config.xml, or else the model being tested."
   )
  ("output-dir", po::value<std::string>(&benchmark_output_dir),
   "With --run-benchmark: where to write the results JSON and conversation logs.  Default:"
   " <userdatadir>/llm_benchmarks; required if there is no user data directory."
   )
  ("list-benchmarks",
   "List the benchmark XML files found (in test_datasets/ and the KnowledgeBase skill directories),"
   " one '<name><tab><path>' per line on stdout, and exit."
   )
  ;
  cl_desc.add( bench_desc );

  po::variables_map cl_vm;
  try
  {
    po::store( po::command_line_parser(argc,argv).options(cl_desc).run(), cl_vm );
    po::notify( cl_vm );

    if( benchmark_files.empty() )
    {
      for( const char *opt : { "skills", "model", "judge-model", "output-dir" } )
      {
        if( cl_vm.count( opt ) )
          throw std::runtime_error( std::string( "--" ) + opt + " is only used with --run-benchmark" );
      }
    }
  }catch( std::exception &e )
  {
    std::cerr << "Command line argument error: " << e.what() << std::endl << cl_desc << std::endl;
    return 1;
  }

  if( cl_vm.count("help") )
  {
    std::cerr << "InterSpec MCP server (stdio).  An MCP client launches this program, and talks to it"
              << " over stdin/stdout.\n"
              << "Or, with --run-benchmark, it runs LLM benchmarks through the InterSpec LLM assistant.\n"
              << cl_desc << std::endl;
    return 0;
  }

  const bool run_benchmarks = !benchmark_files.empty();
  const bool list_benchmarks = (cl_vm.count( "list-benchmarks" ) > 0);
  const bool batch_mode = run_benchmarks || list_benchmarks;

#ifdef _WIN32
  // The path is UTF-8 (see getUtf8Args), which narrow-path functions would read as the ANSI code page.
  if( !log_file.empty() && !_wfreopen( SpecUtils::convert_from_utf8_to_utf16( log_file ).c_str(), L"a", stderr ) )
#else
  if( !log_file.empty() && !std::freopen( log_file.c_str(), "a", stderr ) )
#endif
  {
    std::cerr << "Could not open log file '" << log_file << "'" << std::endl;
    return batch_mode ? sm_exit_setup_error : 2;
  }

  // From here on, stdout carries only the protocol (or the benchmark report).
  std::unique_ptr<ProtocolWriter> writer;
  std::unique_ptr<ReservedStdout> report;
  try
  {
    if( batch_mode )
      report = std::make_unique<ReservedStdout>();
    else
      writer = std::make_unique<ProtocolWriter>();
  }catch( std::exception &e )
  {
    std::cerr << e.what() << std::endl;
    return batch_mode ? sm_exit_setup_error : 3;
  }

  std::cerr << "InterSpec stdio MCP server, version " << InterSpec_VERSION << ", compiled "
            << __DATE__ << std::endl;

  try
  {
    if( static_data_dir.empty() )
      static_data_dir = default_static_data_dir();
    if( static_data_dir.empty()
        || !SpecUtils::is_file( SpecUtils::append_path( static_data_dir, "sandia.decay.xml" ) ) )
      throw std::runtime_error( "InterSpec's data directory was not found; specify it with --datadir." );
    InterSpec::setStaticDataDirectory( static_data_dir );

    if( user_data_dir.empty() )
      user_data_dir = first_existing_dir( AppUtils::candidate_user_data_dirs() );

    if( !user_data_dir.empty() )
    {
      if( !SpecUtils::is_directory( user_data_dir ) )
        throw std::runtime_error( "The user data directory '" + user_data_dir + "' does not exist." );
      set_user_data_dir( user_data_dir );
    }else
    {
      // Things that would be in the user data directory are then looked for in the static one.
      InterSpec::setWritableDataDirectory( static_data_dir );
    }

    if( !llm_config_dir.empty() )
    {
      if( !SpecUtils::is_directory( llm_config_dir ) )
        throw std::runtime_error( "The LLM config directory '" + llm_config_dir + "' does not exist." );
      LlmConfig::setConfigDirectory( llm_config_dir );
    }
  }catch( std::exception &e )
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return batch_mode ? sm_exit_setup_error : 4;
  }

  std::cerr << "Static data directory: " << static_data_dir << "\n"
            << "User data directory: " << (user_data_dir.empty() ? std::string("(none)") : user_data_dir)
            << std::endl;

  if( list_benchmarks )
  {
    std::string listing;
    for( const std::pair<std::string,std::string> &bench : LlmBenchmarkRunner::findAvailableBenchmarks() )
      listing += bench.first + "\t" + bench.second + "\n";
    report->write( listing );
    if( !run_benchmarks )
      return sm_exit_ok;
  }//if( list_benchmarks )

  // Benchmark-run arguments that do not need the LLM config.
  if( run_benchmarks )
  {
    try
    {
      for( std::string &file : benchmark_files )
      {
        if( !SpecUtils::is_absolute_path( file ) )
          SpecUtils::make_canonical_path( file );
        if( !SpecUtils::is_file( file ) )
          throw std::runtime_error( "The benchmark file '" + file + "' does not exist." );
      }

      // Results must not go into an installed data directory, which is what the writable data
      //  directory is set to when there is no user data directory.
      if( benchmark_output_dir.empty() )
      {
        if( user_data_dir.empty() )
          throw std::runtime_error( "There is no user data directory to write results to; specify"
                                    " --output-dir (or --userdatadir)." );
        benchmark_output_dir = SpecUtils::append_path( user_data_dir, "llm_benchmarks" );
      }
      if( !SpecUtils::is_absolute_path( benchmark_output_dir ) )
        SpecUtils::make_canonical_path( benchmark_output_dir );
      if( !SpecUtils::is_directory( benchmark_output_dir ) && !SpecUtils::create_directory( benchmark_output_dir ) )
        throw std::runtime_error( "Could not create the output directory '" + benchmark_output_dir + "'." );

      // Before the session (and so its tool registries) is created.
      LlmKnowledgeBase::setSkillFilter( parse_skills_arg( skills_arg ) );
    }catch( std::exception &e )
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return sm_exit_setup_error;
    }
  }//if( run_benchmarks )

  // Launching this program is the opt-in to MCP, so <McpServer><Enabled> is not required (see
  //  LlmConfig::load); the config is still needed for the tool definitions, and the assistant
  //  settings, if the user has a llm_config.xml.
  std::shared_ptr<const LlmConfig> llm_config;
  try
  {
    llm_config = LlmConfig::load();
    if( !llm_config )
      throw std::runtime_error( "it could not be loaded." );
  }catch( std::exception &e )
  {
    std::cerr << "Error loading the LLM configuration: " << e.what() << std::endl;
    return batch_mode ? sm_exit_setup_error : 5;
  }

  if( run_benchmarks )
  {
    // Benchmarks drive the assistant directly, so do not need <McpServer><AllowAssistantPrompt>.
    try
    {
      if( !llm_config->llmApi.isConfigured() )
        throw std::runtime_error( "no LLM provider is configured (llm_config.xml needs <LlmApi><Enabled>"
                                  " and an <ApiProvider> with an endpoint)." );
      if( !NativeHttp::available() )
        throw std::runtime_error( "this build has no HTTPS client for LLM requests (on Linux, libcurl.so.4"
                                  " must be installed)." );

      // For this run only - never written back to disk.
      if( !model_override.empty() )
      {
        const std::shared_ptr<LlmConfig> config = std::make_shared<LlmConfig>( *llm_config );
        config->llmApi.selectModel( model_override );
        llm_config = config;
      }
    }catch( std::exception &e )
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return sm_exit_setup_error;
    }

    // The session's assistant (LlmToolGui) takes its config from here, rather than reading the file.
    InterSpecServer::set_llm_config( llm_config );

    std::cerr << "LLM requests are made with " << NativeHttp::backendName() << "; model "
              << llm_config->llmApi.model() << "." << std::endl;
  }//if( run_benchmarks )

  // The assistant (when the config allows MCP clients to use it) makes its LLM requests from here.
  if( !run_benchmarks && llm_config->llmApi.isConfigured()
      && (llm_config->mcpServer.allowAssistantPrompt || llm_config->mcpServer.allowAgentInvocation) )
  {
    if( NativeHttp::available() )
      std::cerr << "The assistant will make LLM requests with " << NativeHttp::backendName() << "." << std::endl;
    else
      std::cerr << "Warning: this build has no HTTPS client for LLM requests, so the assistant will fail."
                << std::endl;
  }//if( the assistant is offered )

  // The session.  The environment owns the app (through its session), and deletes it when it is
  //  destroyed.
  // InterSpec_resources/ (its text, report templates, ...) is next to the data directory.
  const std::string app_root = SpecUtils::lexically_normalize_path( SpecUtils::append_path( static_data_dir, ".." ) );
  const std::string app_text = SpecUtils::append_path( app_root, "InterSpec_resources/app_text/InterSpec.xml" );
  if( !SpecUtils::is_file( app_text ) )
    std::cerr << "Warning: '" << app_text << "' was not found, so InterSpec's messages will show as"
              << " ??key?? placeholders; InterSpec_resources/ should be next to the data directory." << std::endl;

  std::unique_ptr<Wt::Test::WTestEnvironment> env;
  std::string session_id;
  try
  {
    env = std::make_unique<Wt::Test::WTestEnvironment>( "", "", Wt::EntryPointType::Application );
    Wt::WServer::instance()->logger().configure( "-* +error +warning +secure" );
    env->setAppRoot( app_root );
    // Files read relative to the doc root would otherwise be looked for in the working directory
    //  the client started us in.
    env->setDocRoot( app_root );

    InterSpecApp * const app = new InterSpecApp( *env );
    if( !app->viewer() )
      throw std::runtime_error( "the InterSpec session could not be created." );
    session_id = app->sessionId();

    // The environment took the session lock on this thread; release it, so work scheduled into the
    //  session can run.
    env->endRequest();
  }catch( std::exception &e )
  {
    std::cerr << "Error starting the session: " << e.what() << std::endl;
    return run_benchmarks ? sm_exit_setup_error : 6;
  }

  // Chart images (get_spectrum_image, get_time_chart_image) are drawn here, as there is no browser.
  SpectrumChartImage::setRasterRenderer( []( const int width, const int height, const std::string &format,
                                             const std::function<void( Wt::WPaintDevice & )> &paint ){
    PlutoPaintDevice device( width, height );
    paint( device );
    return device.encode( format );
  } );

  // The assistant's run_javascript tool, as there is no browser for its sandbox.
  LlmTools::setJavaScriptRunner( &run_javascript );

  // When done: let whatever is running in the session finish (taking the session lock waits for
  //  it), and tear down.  A calculation still running (an isotopics fit, an assistant turn) could
  //  hold that up for minutes; nothing needs saving (the database is in memory), so do not wait
  //  more than a few seconds for it.
  const auto exit_within_seconds = []( const int exit_code ){
    std::thread( [exit_code](){
      std::this_thread::sleep_for( std::chrono::seconds( 5 ) );
      std::cerr << "Still busy; exiting without waiting for the running calculation." << std::endl;
#ifdef _WIN32
      // _Exit() ends in ExitProcess(), whose DLL detach notifications could wait on threads it has
      //  already killed (e.g., one holding the loader lock); TerminateProcess skips them.
      ::TerminateProcess( ::GetCurrentProcess(), exit_code );
#endif
      std::_Exit( exit_code );
    } ).detach();
  };//exit_within_seconds

  if( run_benchmarks )
  {
    std::signal( SIGINT, &handle_sigint );

    BenchmarkRunOptions run_options;
    run_options.interactive = false;
    run_options.outputDir = benchmark_output_dir;
    run_options.judgeModel = judge_model;

    int exit_code = sm_exit_ok;
    for( const std::string &xml_path : benchmark_files )
    {
      if( ns_interrupted )
        break;

      const std::shared_ptr<BenchmarkRunState> state = std::make_shared<BenchmarkRunState>();
      const auto start_failed = [state]( const std::string &why ){
        std::lock_guard<std::mutex> lock( state->mutex );
        state->started = state->startFailed = true;
        state->startError = why;
        state->cv.notify_all();
      };

      std::cerr << "Running benchmark " << xml_path << std::endl;
      Wt::WServer::instance()->schedule( std::chrono::milliseconds(0), session_id,
        [state, xml_path, run_options, start_failed](){
          InterSpec * const viewer = InterSpec::instance();
          if( !viewer )
            return start_failed( "there is no InterSpec session." );
          if( !viewer->currentLlmTool() )
            viewer->createLlmTool();
          LlmToolGui * const gui = viewer->currentLlmTool();
          if( !gui || !gui->llmInterface() )
            return start_failed( "the LLM assistant could not be created." );

          const bool started = gui->startBenchmark( xml_path, run_options,
            [state]( const LlmBenchmarkRunner &runner ){
              std::lock_guard<std::mutex> lock( state->mutex );
              state->finished = true;
              state->completed = runner.completed();
              state->resultsPath = runner.resultsFilePath();
              state->results = runner.results();
              state->cv.notify_all();
          } );
          if( !started )
            return start_failed( "a benchmark is already running." );

          std::lock_guard<std::mutex> lock( state->mutex );
          state->started = true;
          state->cv.notify_all();
        },
        [start_failed](){ start_failed( "the session has gone away." ); } );

      // Wait for the run, never holding the session lock; on Ctrl-C, cancel it (the checkpoint is
      //  written after every question, so the partial results are kept).
      bool cancel_sent = false;
      std::unique_lock<std::mutex> lock( state->mutex );
      while( !state->finished && !state->startFailed )
      {
        state->cv.wait_for( lock, std::chrono::milliseconds( 200 ) );
        if( ns_interrupted && !cancel_sent && state->started && !state->finished )
        {
          cancel_sent = true;
          std::cerr << "Interrupted; cancelling the benchmark." << std::endl;
          lock.unlock();
          Wt::WServer::instance()->schedule( std::chrono::milliseconds(0), session_id, [](){
            InterSpec * const viewer = InterSpec::instance();
            LlmToolGui * const gui = viewer ? viewer->currentLlmTool() : nullptr;
            if( gui )
              gui->cancelBenchmark();
          } );
          lock.lock();
        }
      }//while( waiting for the run )

      if( state->startFailed )
      {
        std::cerr << "Error: could not run benchmark '" << xml_path << "': " << state->startError << std::endl;
        exit_code = sm_exit_benchmark_error;
        continue;
      }

      report->write( benchmark_summary( state->results, state->resultsPath, state->completed ) );
      if( !state->completed )
      {
        std::cerr << "Benchmark '" << xml_path << "' did not run to the end." << std::endl;
        exit_code = sm_exit_benchmark_error;
      }
    }//for( const std::string &xml_path : benchmark_files )

    if( ns_interrupted )
      exit_code = sm_exit_benchmark_error;

    exit_within_seconds( exit_code );
    env->startRequest();
    env.reset();

    return exit_code;
  }//if( run_benchmarks )

  LlmMcpDispatcher::Options options;
  options.session_argument = false;
  options.browser_available = false;
  options.instructions = load_mcp_instructions();
  std::unique_ptr<LlmMcpDispatcher> dispatcher;
  try
  {
    dispatcher = std::make_unique<LlmMcpDispatcher>( llm_config, options,
                                    [session_id]( const std::string & ){ return session_id; } );
  }catch( std::exception &e )
  {
    std::cerr << "Error setting up the MCP tools: " << e.what() << std::endl;
    return 7;
  }

  LlmMcpDispatcher::RequestContext context;
  context.can_stream = true;

  ProtocolWriter * const out = writer.get();
  const auto send = [out]( const json &message ){ out->send( message ); };

  std::string line;
  while( std::getline( std::cin, line ) )
  {
    if( !line.empty() && (line.back() == '\r') )
      line.pop_back();
    if( line.find_first_not_of( " \t" ) == std::string::npos )
      continue;

    json message;
    try
    {
      message = json::parse( line );
    }catch( std::exception &e )
    {
      send( LlmMcpDispatcher::rpcError( nullptr, -32700, std::string( "Parse error: " ) + e.what() ) );
      continue;
    }

    const LlmMcpDispatcher::Outcome outcome = dispatcher->handle( message, context );
    if( outcome.pending )
      LlmMcpDispatcher::deliverWhenReady( outcome.pending, send );
    else if( outcome.reply )
      send( *outcome.reply );
  }//while( reading messages )

  // The client closed stdin: stop waiting on replies no one will read, and tear down.
  std::cerr << "stdin closed; shutting down." << std::endl;
  exit_within_seconds( 0 );
  LlmMcpDispatcher::cancelDeliveries();
  dispatcher.reset();

  env->startRequest();
  env.reset();

  return 0;
}//int main( int argc, char **argv )
