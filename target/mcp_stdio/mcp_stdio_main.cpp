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
 */

#include "InterSpec_config.h"

#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <cstdio>
#include <thread>
#include <memory>
#include <chrono>
#include <csignal>
#include <variant>
#include <iostream>
#include <stdexcept>
#include <functional>

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
#include "SpecUtils/SerialToDetectorModel.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmMcpDispatcher.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/NativeHttpClient.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/SpectrumChartImage.h"
#include "InterSpec/DecayDataBaseServer.h"

#include "QuickJsSandbox.h"
#include "PlutoPaintDevice.h"

static_assert( BUILD_AS_MCP_STDIO_SERVER, "This file is only for the stdio MCP server build" );

using namespace std;
using json = nlohmann::json;


namespace
{
  /** Writes protocol messages to the original stdout; everything else written to stdout (by
   InterSpec, Wt, or libraries) goes to stderr instead, so it can not corrupt the protocol.
   */
  class ProtocolWriter
  {
  public:
    ProtocolWriter()
    {
#ifdef _WIN32
      _setmode( _fileno(stdin), _O_BINARY );
      _setmode( _fileno(stdout), _O_BINARY );
      fflush( stdout );
      m_fd = _dup( _fileno(stdout) );
      if( (m_fd < 0) || (_dup2( _fileno(stderr), _fileno(stdout) ) != 0) )
        throw std::runtime_error( "Could not redirect stdout." );
#else
      fflush( stdout );
      m_fd = dup( STDOUT_FILENO );
      if( (m_fd < 0) || (dup2( STDERR_FILENO, STDOUT_FILENO ) < 0) )
        throw std::runtime_error( "Could not redirect stdout." );
#endif
    }//ProtocolWriter()

    /** Writes `message` as one line; returns false if the client has gone away. */
    bool send( const json &message )
    {
      const std::string line = LlmMcpDispatcher::dumpMessage( message ) + "\n";

      std::lock_guard<std::mutex> lock( m_mutex );
      size_t written = 0;
      while( written < line.size() )
      {
#ifdef _WIN32
        const int n = _write( m_fd, line.data() + written, static_cast<unsigned int>(line.size() - written) );
#else
        const ssize_t n = write( m_fd, line.data() + written, line.size() - written );
#endif
        if( n <= 0 )
          return false;
        written += static_cast<size_t>( n );
      }
      return true;
    }//send(...)

  private:
    int m_fd = -1;
    std::mutex m_mutex;
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

  po::variables_map cl_vm;
  try
  {
    po::store( po::command_line_parser(argc,argv).options(cl_desc).run(), cl_vm );
    po::notify( cl_vm );
  }catch( std::exception &e )
  {
    std::cerr << "Command line argument error: " << e.what() << std::endl << cl_desc << std::endl;
    return 1;
  }

  if( cl_vm.count("help") )
  {
    std::cerr << "InterSpec MCP server (stdio).  An MCP client launches this program, and talks to it"
              << " over stdin/stdout.\n" << cl_desc << std::endl;
    return 0;
  }

  if( !log_file.empty() && !std::freopen( log_file.c_str(), "a", stderr ) )
  {
    std::cerr << "Could not open log file '" << log_file << "'" << std::endl;
    return 2;
  }

  // From here on, stdout carries only the protocol.
  std::unique_ptr<ProtocolWriter> writer;
  try
  {
    writer = std::make_unique<ProtocolWriter>();
  }catch( std::exception &e )
  {
    std::cerr << e.what() << std::endl;
    return 3;
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
    return 4;
  }

  std::cerr << "Static data directory: " << static_data_dir << "\n"
            << "User data directory: " << (user_data_dir.empty() ? std::string("(none)") : user_data_dir)
            << std::endl;

  // Launching this program is the opt-in to MCP, so <McpServer><Enabled> is not required; the
  //  config is still needed for the tool definitions, and the assistant settings.
  std::shared_ptr<const LlmConfig> llm_config;
  try
  {
    llm_config = LlmConfig::load();
    if( !llm_config )
      throw std::runtime_error( "no llm_config.xml was found." );
  }catch( std::exception &e )
  {
    std::cerr << "Error loading the LLM configuration: " << e.what() << std::endl;
    return 5;
  }

  // The assistant (when the config allows MCP clients to use it) makes its LLM requests from here.
  if( llm_config->llmApi.enabled && (llm_config->mcpServer.allowAssistantPrompt
                                     || llm_config->mcpServer.allowAgentInvocation) )
  {
    if( NativeHttp::available() )
      std::cerr << "The assistant will make LLM requests with " << NativeHttp::backendName() << "." << std::endl;
    else
      std::cerr << "Warning: this build has no HTTPS client for LLM requests, so the assistant will fail."
                << std::endl;
  }//if( the assistant is offered )

  // The session.  The environment owns the app (through its session), and deletes it when it is
  //  destroyed.
  std::unique_ptr<Wt::Test::WTestEnvironment> env;
  std::string session_id;
  try
  {
    env = std::make_unique<Wt::Test::WTestEnvironment>( "", "", Wt::EntryPointType::Application );
    Wt::WServer::instance()->logger().configure( "-* +error +warning +secure" );
    env->setAppRoot( SpecUtils::lexically_normalize_path( SpecUtils::append_path( static_data_dir, ".." ) ) );

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
    return 6;
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

  LlmMcpDispatcher::Options options;
  options.session_argument = false;
  options.browser_available = false;
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

  // The client closed stdin: stop waiting on replies no one will read, let whatever is running in
  //  the session finish (taking the session lock waits for it), and tear down.  A calculation still
  //  running (an isotopics fit, an assistant turn) could hold that up for minutes; nothing needs
  //  saving (the database is in memory), so do not wait more than a few seconds for it.
  std::cerr << "stdin closed; shutting down." << std::endl;
  std::thread( [](){
    std::this_thread::sleep_for( std::chrono::seconds( 5 ) );
    std::cerr << "Still busy; exiting without waiting for the running calculation." << std::endl;
    std::_Exit( 0 );
  } ).detach();
  LlmMcpDispatcher::cancelDeliveries();
  dispatcher.reset();

  env->startRequest();
  env.reset();

  return 0;
}//int main( int argc, char **argv )
