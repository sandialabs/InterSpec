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

#include <ctime>
#include <regex>
#include <cstring>
#include <atomic>
#include <thread>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <system_error>

#if( !defined(_WIN32) )
#include <spawn.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
#endif

#include <Wt/WServer.h>
#include <Wt/WRandom.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecServer.h"
#include "target/mcp_http/McpSessionManager.h"

#if( !defined(_WIN32) )
extern char **environ;
#endif

using namespace std;
using json = nlohmann::json;
namespace fs = std::filesystem;


namespace
{
  std::mutex ns_instance_mutex;
  std::shared_ptr<McpSessionManager> ns_instance;

  const char * const sm_marker_filename = ".interspec-session";

  /** How long a session may take to come up before it is abandoned. */
  const std::chrono::seconds sm_start_timeout( 60 );

  /** How long after a session ends its directory is removed, if its app has not ended by then. */
  const std::chrono::seconds sm_dir_removal_grace( 30 );

  /** How long a record of an ended session is kept, so clients get "session ended" rather than
   "unknown session". */
  const std::chrono::minutes sm_ended_record_lifetime( 60 );

  /** Kept short, as the server's io_service does not finish stopping until pending timers fire. */
  const std::chrono::seconds sm_reaper_interval( 10 );

  std::string utc_date_string( const std::chrono::system_clock::time_point &tp )
  {
    const std::time_t t = std::chrono::system_clock::to_time_t( tp );
    std::tm tm_utc{};
#ifdef _WIN32
    gmtime_s( &tm_utc, &t );
#else
    gmtime_r( &t, &tm_utc );
#endif
    char buffer[32] = { '\0' };
    std::strftime( buffer, sizeof(buffer), "%Y-%m-%d", &tm_utc );
    return buffer;
  }//utc_date_string(...)


  std::string utc_iso_string( const std::chrono::system_clock::time_point &tp )
  {
    const std::time_t t = std::chrono::system_clock::to_time_t( tp );
    std::tm tm_utc{};
#ifdef _WIN32
    gmtime_s( &tm_utc, &t );
#else
    gmtime_r( &t, &tm_utc );
#endif
    char buffer[32] = { '\0' };
    std::strftime( buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &tm_utc );
    return buffer;
  }//utc_iso_string(...)


  /** True if `path` is equal to, or inside of, `base`; both must be lexically normal and absolute. */
  bool path_is_within( const fs::path &base, const fs::path &path )
  {
    const auto mismatch = std::mismatch( base.begin(), base.end(), path.begin(), path.end() );
    return (mismatch.first == base.end());
  }


  /** Session directory names are "<YYYY-MM-DD>_<alphanumeric handle>". */
  bool is_session_dir_name( const std::string &name )
  {
    static const std::regex pattern( "^[0-9]{4}-[0-9]{2}-[0-9]{2}_[A-Za-z0-9]+$" );
    return std::regex_match( name, pattern );
  }
}//namespace


#if( !defined(_WIN32) )
/** A headless Chromium driven over the DevTools protocol on a pipe (`--remote-debugging-pipe`), so
 there is no DevTools network port for page JavaScript (which includes LLM-written code) to reach.

 Messages are JSON, terminated by a NUL; the browser reads commands on its fd 3, and writes replies
 and events on its fd 4.  A reader thread dispatches replies to their callbacks, and reports
 targets (tabs) that go away, and the browser exiting.
 */
class CdpBrowser
{
public:
  using TargetGoneFcn = std::function<void(const std::string &target_id)>;
  using ExitFcn = std::function<void()>;

  /** Starts the browser; @throws std::exception on failure. */
  CdpBrowser( const std::string &executable, const std::vector<std::string> &extra_args,
              TargetGoneFcn on_target_gone, ExitFcn on_exit )
    : m_pid( -1 ),
      m_to_browser( -1 ),
      m_from_browser( -1 ),
      m_running( false ),
      m_next_id( 1 ),
      m_on_target_gone( std::move(on_target_gone) ),
      m_on_exit( std::move(on_exit) )
  {
    // A throwaway profile, so nothing is shared with, or left for, a later run.
    const string tmp_dir = SpecUtils::temp_dir();
    m_profile_dir = SpecUtils::append_path( tmp_dir, "interspec-browser-" + Wt::WRandom::generateId(12) );
    std::error_code ec;
    fs::create_directories( m_profile_dir, ec );
    if( ec )
      throw runtime_error( "Could not create browser profile directory: " + ec.message() );

    int to_pipe[2] = { -1, -1 }, from_pipe[2] = { -1, -1 };
    if( (pipe( to_pipe ) != 0) || (pipe( from_pipe ) != 0) )
      throw runtime_error( "Could not create pipes for the browser." );

    // Move every end above fd 4 (with close-on-exec), so the dup2()s below always make a new
    //  descriptor - dup2() onto the same number is a no-op that would keep close-on-exec set.
    const auto move_high = []( int &fd ){
      const int high = fcntl( fd, F_DUPFD_CLOEXEC, 10 );
      close( fd );
      fd = high;
    };
    for( int *fd : { &to_pipe[0], &to_pipe[1], &from_pipe[0], &from_pipe[1] } )
      move_high( *fd );

    if( (to_pipe[0] < 0) || (to_pipe[1] < 0) || (from_pipe[0] < 0) || (from_pipe[1] < 0) )
      throw runtime_error( "Could not set up pipes for the browser." );

    vector<string> args;
    args.push_back( executable );
    const string exe_name = SpecUtils::filename( executable );
    if( exe_name.find( "headless" ) == string::npos )
      args.push_back( "--headless=new" );
    args.push_back( "--remote-debugging-pipe" );
    args.push_back( "--user-data-dir=" + m_profile_dir );
    args.push_back( "--no-first-run" );
    args.push_back( "--no-default-browser-check" );
    args.push_back( "--disable-gpu" );
    args.push_back( "--disable-extensions" );
    args.push_back( "--disable-sync" );
    args.push_back( "--disable-background-networking" );
    args.push_back( "--disable-component-update" );
    args.push_back( "--disable-default-apps" );
    args.push_back( "--disable-breakpad" );
    args.push_back( "--disable-client-side-phishing-detection" );
    args.push_back( "--disable-domain-reliability" );
    args.push_back( "--no-pings" );
    args.push_back( "--mute-audio" );
    // Sessions are in background tabs, but their timers and server-push must keep running.
    args.push_back( "--disable-background-timer-throttling" );
    args.push_back( "--disable-renderer-backgrounding" );
    args.push_back( "--disable-backgrounding-occluded-windows" );
    args.push_back( "--disable-features=IntensiveWakeUpThrottling,CalculateNativeWinOcclusion,"
                    "Translate,OptimizationHints,MediaRouter,DialMediaRouteProvider" );
    // Everything goes to a proxy that is not there - including loopback, which Chromium would
    //  otherwise never proxy ("<-loopback>") - except InterSpec itself (see `extra_args`), so page
    //  JavaScript can not reach the internet, the internal network, or other local services.
    //  WebRTC UDP would bypass the proxy, so it is limited to proxied connections (i.e., none).
    args.push_back( "--proxy-server=http://127.0.0.1:9" );
    args.push_back( "--force-webrtc-ip-handling-policy=disable_non_proxied_udp" );
    args.push_back( "--window-size=1600,1000" );
#if( defined(__linux__) )
    // Containers default to a 64 MB /dev/shm, which is not enough for Chromium.
    args.push_back( "--disable-dev-shm-usage" );
#endif
    const char * const no_sandbox = std::getenv( "INTERSPEC_BROWSER_NO_SANDBOX" );
    if( no_sandbox && (string(no_sandbox) == "1") )
      args.push_back( "--no-sandbox" );
    args.insert( end(args), begin(extra_args), end(extra_args) );
    args.push_back( "about:blank" );

    vector<char *> argv;
    for( string &arg : args )
      argv.push_back( &arg[0] );
    argv.push_back( nullptr );

    // The browser gets our environment, minus our own INTERSPEC_* variables (API keys and such).
    vector<string> env_strs;
    for( char **env = environ; env && *env; ++env )
    {
      if( !SpecUtils::starts_with( *env, "INTERSPEC_" ) )
        env_strs.push_back( *env );
    }
    vector<char *> envp;
    for( string &env : env_strs )
      envp.push_back( &env[0] );
    envp.push_back( nullptr );

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init( &actions );
    posix_spawn_file_actions_adddup2( &actions, to_pipe[0], 3 );
    posix_spawn_file_actions_adddup2( &actions, from_pipe[1], 4 );
    posix_spawn_file_actions_addopen( &actions, 0, "/dev/null", O_RDONLY, 0 );
    const char * const verbose = std::getenv( "INTERSPEC_BROWSER_VERBOSE" );
    if( !verbose || (string(verbose) != "1") )
    {
      posix_spawn_file_actions_addopen( &actions, 1, "/dev/null", O_WRONLY, 0 );
      posix_spawn_file_actions_addopen( &actions, 2, "/dev/null", O_WRONLY, 0 );
    }

    // Its own process group, so the renderer processes can be signalled along with it.
    posix_spawnattr_t attr;
    posix_spawnattr_init( &attr );
    posix_spawnattr_setflags( &attr, POSIX_SPAWN_SETPGROUP );
    posix_spawnattr_setpgroup( &attr, 0 );

    pid_t pid = -1;
    const int spawn_status = (executable.find('/') != string::npos)
                     ? posix_spawn( &pid, executable.c_str(), &actions, &attr, argv.data(), envp.data() )
                     : posix_spawnp( &pid, executable.c_str(), &actions, &attr, argv.data(), envp.data() );
    posix_spawn_file_actions_destroy( &actions );
    posix_spawnattr_destroy( &attr );

    close( to_pipe[0] );
    close( from_pipe[1] );

    if( spawn_status != 0 )
    {
      close( to_pipe[1] );
      close( from_pipe[0] );
      throw runtime_error( "Could not start the browser '" + executable + "': "
                           + std::string( strerror(spawn_status) ) );
    }

    m_pid = pid;
    m_to_browser = to_pipe[1];
    m_from_browser = from_pipe[0];
    m_running = true;

    m_reader = std::thread( [this](){ readLoop(); } );

    // Get told of tabs that are destroyed or crash.
    send( "Target.setDiscoverTargets", json{ {"discover", true} }, nullptr );
  }//CdpBrowser constructor


  ~CdpBrowser()
  {
    shutdown();

    if( m_reader.joinable() )
    {
      if( m_reader.get_id() == std::this_thread::get_id() )
        m_reader.detach();
      else
        m_reader.join();
    }

    std::error_code ec;
    fs::remove_all( m_profile_dir, ec );
  }//~CdpBrowser()


  bool running() const
  {
    return m_running.load();
  }


  /** Opens a tab at `url`; `callback` is called with the target id, or an error, from the reader
   thread.
   */
  void createTarget( const std::string &url,
                     std::function<void(const std::string &target_id, const std::string &error)> callback )
  {
    // A background tab, rather than a window each, which costs about twice the memory.  Its page is
    //  "hidden", so requestAnimationFrame does not fire, but the chart captures race that with a
    //  timer, and timers are not throttled (see the command-line flags).
    send( "Target.createTarget", json{ {"url", url}, {"background", true} },
      [callback]( const json &reply ){
        if( reply.contains("result") && reply["result"].contains("targetId") )
          callback( reply["result"]["targetId"].get<string>(), "" );
        else
          callback( "", reply.contains("error") ? reply["error"].dump() : string("no reply") );
    } );
  }//createTarget(...)


  void closeTarget( const std::string &target_id )
  {
    if( !target_id.empty() )
      send( "Target.closeTarget", json{ {"targetId", target_id} }, nullptr );
  }


  /** Stops the browser: asks it to close, then terminates its process group. */
  void shutdown()
  {
    const bool was_running = m_running.exchange( false );

    if( was_running )
      send( "Browser.close", json::object(), nullptr, true );

    {
      std::lock_guard<std::mutex> lock( m_write_mutex );
      if( m_to_browser >= 0 )
      {
        close( m_to_browser );
        m_to_browser = -1;
      }
    }

    if( m_pid > 0 )
    {
      // Give it a moment to exit on its own, then make sure.
      for( int i = 0; i < 20; ++i )
      {
        int status = 0;
        if( waitpid( m_pid, &status, WNOHANG ) == m_pid )
        {
          m_pid = -1;
          break;
        }
        std::this_thread::sleep_for( std::chrono::milliseconds(50) );
      }

      if( m_pid > 0 )
      {
        kill( -m_pid, SIGKILL );
        int status = 0;
        waitpid( m_pid, &status, 0 );
        m_pid = -1;
      }
    }//if( m_pid > 0 )
  }//shutdown()


private:
  /** Sends a command; `callback`, if given, gets the reply (on the reader thread), or an error
   (possibly before this function returns, but never with a lock held).
   */
  void send( const std::string &method, const json &params,
             std::function<void(const json &)> callback, const bool even_if_stopped = false )
  {
    std::unique_lock<std::mutex> lock( m_write_mutex );
    if( (m_to_browser < 0) || (!m_running.load() && !even_if_stopped) )
    {
      lock.unlock();
      if( callback )
        callback( json{ {"error", "browser not running"} } );
      return;
    }

    const int id = m_next_id++;
    if( callback )
    {
      std::lock_guard<std::mutex> pending_lock( m_pending_mutex );
      m_pending[id] = std::move( callback );
    }

    string message = json{ {"id", id}, {"method", method}, {"params", params} }.dump();
    message.push_back( '\0' );

    size_t written = 0;
    while( written < message.size() )
    {
      const ssize_t nwrote = ::write( m_to_browser, message.data() + written, message.size() - written );
      if( nwrote < 0 )
      {
        if( errno == EINTR )
          continue;
        break;
      }
      written += static_cast<size_t>( nwrote );
    }//while( written < message.size() )
  }//send(...)


  void readLoop()
  {
    string buffer;
    char chunk[16384];

    while( true )
    {
      const ssize_t nread = ::read( m_from_browser, chunk, sizeof(chunk) );
      if( nread < 0 )
      {
        if( errno == EINTR )
          continue;
        break;
      }
      if( nread == 0 )
        break;

      buffer.append( chunk, static_cast<size_t>(nread) );

      size_t end_pos = 0;
      while( (end_pos = buffer.find( '\0' )) != string::npos )
      {
        const string message = buffer.substr( 0, end_pos );
        buffer.erase( 0, end_pos + 1 );
        handleMessage( message );
      }
    }//while( true )

    close( m_from_browser );
    m_from_browser = -1;

    const bool was_running = m_running.exchange( false );

    // Fail anything still waiting on a reply.
    map<int,std::function<void(const json &)>> pending;
    {
      std::lock_guard<std::mutex> lock( m_pending_mutex );
      pending.swap( m_pending );
    }
    for( auto &id_fcn : pending )
      id_fcn.second( json{ {"error", "browser exited"} } );

    if( was_running && m_on_exit )
      m_on_exit();
  }//readLoop()


  void handleMessage( const std::string &message )
  {
    json msg;
    try
    {
      msg = json::parse( message );
    }catch( std::exception & )
    {
      return;
    }

    if( msg.contains( "id" ) && msg["id"].is_number_integer() )
    {
      std::function<void(const json &)> callback;
      {
        std::lock_guard<std::mutex> lock( m_pending_mutex );
        const auto pos = m_pending.find( msg["id"].get<int>() );
        if( pos != end(m_pending) )
        {
          callback = std::move( pos->second );
          m_pending.erase( pos );
        }
      }
      if( callback )
        callback( msg );
      return;
    }//if( a reply )

    const string method = msg.value( "method", "" );
    if( ((method == "Target.targetDestroyed") || (method == "Target.targetCrashed"))
        && msg.contains("params") && m_on_target_gone )
    {
      const string target_id = msg["params"].value( "targetId", "" );
      if( !target_id.empty() )
        m_on_target_gone( target_id );
    }
  }//handleMessage(...)


  pid_t m_pid;
  int m_to_browser;
  int m_from_browser;
  std::atomic<bool> m_running;
  std::string m_profile_dir;

  std::mutex m_write_mutex;
  int m_next_id;

  std::mutex m_pending_mutex;
  std::map<int,std::function<void(const json &)>> m_pending;

  std::thread m_reader;
  const TargetGoneFcn m_on_target_gone;
  const ExitFcn m_on_exit;
};//class CdpBrowser


namespace
{
  /** Finds the browser executable: the configured one, else the first of the usual names/places. */
  std::string find_browser( const std::string &configured )
  {
    const auto is_executable = []( const std::string &path ) -> bool {
      return !path.empty() && SpecUtils::is_file( path ) && (access( path.c_str(), X_OK ) == 0);
    };

    const auto search_path = [&is_executable]( const std::string &name ) -> std::string {
      if( name.find('/') != string::npos )
        return is_executable( name ) ? name : string();

      const char * const path_env = std::getenv( "PATH" );
      vector<string> dirs;
      SpecUtils::split( dirs, path_env ? path_env : "/usr/bin:/bin", ":" );
      for( const string &dir : dirs )
      {
        const string candidate = SpecUtils::append_path( dir, name );
        if( is_executable( candidate ) )
          return candidate;
      }
      return string();
    };

    if( !configured.empty() )
    {
      const string found = search_path( configured );
      if( found.empty() )
        throw runtime_error( "The configured browser '" + configured + "' was not found, or is not executable." );
      return found;
    }

    const char * const candidates[] = {
      "chromium-headless-shell",
      "/usr/lib/chromium/chromium-headless-shell",
      "/usr/lib/chromium/headless_shell",
      "chromium",
      "chromium-browser",
      "google-chrome",
      "google-chrome-stable",
      "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
      "/Applications/Chromium.app/Contents/MacOS/Chromium"
    };

    for( const char * const candidate : candidates )
    {
      const string found = search_path( candidate );
      if( !found.empty() )
        return found;
    }

    throw runtime_error( "No Chromium/Chrome browser found; set <BrowserPath> in the <AgentServer>"
                         " element of llm_config.xml." );
  }//find_browser(...)
}//namespace
#else
class CdpBrowser
{
};
#endif //#if( !defined(_WIN32) )


McpSessionManager::McpSessionManager( const AgentSettings &settings, const std::string &app_host,
                                      const int app_port )
  : m_settings( settings ),
    m_app_host( app_host ),
    m_app_port( app_port ),
    m_running( true )
{
#if( defined(_WIN32) )
  throw runtime_error( "Running as an MCP agent server is not supported on Windows." );
#else
  if( m_settings.exchangeDirectory.empty() )
    throw runtime_error( "No exchange directory specified." );

  std::error_code ec;
  fs::create_directories( m_settings.exchangeDirectory, ec );
  if( ec )
    throw runtime_error( "Could not create exchange directory '" + m_settings.exchangeDirectory
                         + "': " + ec.message() );

  m_exchange_root = fs::canonical( m_settings.exchangeDirectory, ec ).string();
  if( ec )
    throw runtime_error( "Invalid exchange directory '" + m_settings.exchangeDirectory + "': " + ec.message() );

  if( !SpecUtils::can_rw_in_directory( m_exchange_root ) )
    throw runtime_error( "Can not write to exchange directory '" + m_exchange_root + "'." );

  // Fail at startup, rather than at the first create_session, if there is no browser.
  find_browser( m_settings.browserPath );

  if( m_settings.maxSessions < 1 )
    throw runtime_error( "<MaxSessions> must be at least 1." );
#endif
}//McpSessionManager constructor


McpSessionManager::~McpSessionManager()
{
}


void McpSessionManager::start( const AgentSettings &settings, const std::string &app_host, const int app_port )
{
  std::shared_ptr<McpSessionManager> manager( new McpSessionManager( settings, app_host, app_port ) );
  manager->sweepOrphanedDirectories();

  {
    std::lock_guard<std::mutex> lock( ns_instance_mutex );
    if( ns_instance )
      throw logic_error( "McpSessionManager already started." );
    ns_instance = manager;
  }

  manager->scheduleReaper();

  cout << "MCP agent server: exchange directory '" << manager->m_exchange_root << "', at most "
       << settings.maxSessions << " sessions." << endl;
}//start(...)


void McpSessionManager::stop()
{
  std::shared_ptr<McpSessionManager> manager;
  {
    std::lock_guard<std::mutex> lock( ns_instance_mutex );
    manager.swap( ns_instance );
  }

  if( !manager )
    return;

  vector<std::function<void()>> listeners;
  std::shared_ptr<CdpBrowser> browser;
  {
    std::lock_guard<std::mutex> lock( manager->m_mutex );
    manager->m_running = false;
    for( auto &handle_session : manager->m_sessions )
    {
      vector<std::function<void()>> these = manager->endSessionLocked( handle_session.second,
                                                                       "the server is shutting down" );
      listeners.insert( end(listeners), begin(these), end(these) );
      manager->removeSessionDirectory( handle_session.second );
    }
    browser.swap( manager->m_browser );
  }

  for( const std::function<void()> &fcn : listeners )
    fcn();

#if( !defined(_WIN32) )
  if( browser )
    browser->shutdown();
#endif
}//stop()


std::shared_ptr<McpSessionManager> McpSessionManager::instance()
{
  std::lock_guard<std::mutex> lock( ns_instance_mutex );
  return ns_instance;
}


const McpSessionManager::AgentSettings &McpSessionManager::settings() const
{
  return m_settings;
}


std::string McpSessionManager::createSession()
{
#if( defined(_WIN32) )
  throw runtime_error( "Not supported on Windows." );
#else
  string handle, url;
  std::shared_ptr<CdpBrowser> browser;

  // Destroyed after our lock is released: its destructor joins its reader thread, which may be
  //  waiting on our lock (to end the sessions of the browser).
  std::shared_ptr<CdpBrowser> dead_browser;

  {
    std::lock_guard<std::mutex> lock( m_mutex );

    if( !m_running )
      throw runtime_error( "The server is shutting down." );

    size_t num_open = 0;
    for( const auto &handle_session : m_sessions )
      num_open += (handle_session.second.info.state != State::Ended);
    if( num_open >= m_settings.maxSessions )
      throw runtime_error( "The maximum number of sessions (" + std::to_string(m_settings.maxSessions)
                           + ") are open; close one with close_session, or try again later." );

    do
    {
      handle = Wt::WRandom::generateId( 24 );
    }while( m_sessions.count( handle ) );

    // The browser is started first, so a failure leaves no session directory behind.
    if( m_browser && !m_browser->running() )
      dead_browser.swap( m_browser );

    if( !m_browser )
    {
      // Callbacks come from the browser's reader thread; they hold the manager only weakly.
      const std::weak_ptr<McpSessionManager> weak_self = shared_from_this();
      // Only InterSpec's own port bypasses the browser's (dead) proxy.
      const string port = std::to_string( m_app_port );
      const uint64_t generation = ++m_browser_generation;
      const vector<string> extra_args{
        "--proxy-bypass-list=<-loopback>;*.localhost:" + port + ";localhost:" + port + ";" + m_app_host + ":" + port
      };
      m_browser = std::make_shared<CdpBrowser>( find_browser( m_settings.browserPath ), extra_args,
        [weak_self]( const std::string &target_id ){
          const std::shared_ptr<McpSessionManager> self = weak_self.lock();
          if( self )
            self->handleTargetGone( target_id );
        },
        [weak_self, generation](){
          const std::shared_ptr<McpSessionManager> self = weak_self.lock();
          if( self )
            self->handleBrowserExited( generation );
        } );
    }//if( !m_browser )

    Session session;
    session.info.handle = handle;
    session.info.created = std::chrono::system_clock::now();
    session.info.expires = session.info.created + std::chrono::seconds( m_settings.sessionMaxLifetimeSeconds );
    session.created = session.last_activity = std::chrono::steady_clock::now();
    session.info.dir_name = utc_date_string( session.info.created ) + "_" + handle;

    const fs::path dir = fs::path( m_exchange_root ) / session.info.dir_name;
    const auto fail = [&dir]( const string &message ){
      std::error_code remove_ec;
      fs::remove_all( dir, remove_ec );
      throw runtime_error( message );
    };

    std::error_code ec;
    for( const char * const sub_dir : { "input", "output", "logs" } )
    {
      fs::create_directories( dir / sub_dir, ec );
      if( ec )
        fail( "Could not create session directory: " + ec.message() );
    }

    {// The marker is what allows the directory to be deleted later.
      ofstream marker( (dir / sm_marker_filename).string() );
      marker << json{ {"handle", handle}, {"created", utc_iso_string(session.info.created)} }.dump() << "\n";
      if( !marker )
        fail( "Could not write the session directory marker file." );
    }

    session.info.dir = fs::canonical( dir, ec ).string();
    if( ec )
      fail( "Could not resolve session directory: " + ec.message() );

    session.browser_generation = m_browser_generation;

    session.info.upload_token = Wt::WRandom::generateId( 32 );
    m_file_tokens[session.info.upload_token] = FileToken{ handle, true, "" };

    InterSpecServer::add_allowed_session_token( handle.c_str(),
                                                InterSpecServer::SessionType::ExternalBrowserInstance );

    m_sessions[handle] = std::move( session );
    browser = m_browser;
    url = sessionAppUrl( handle );
  }

  // Not under our lock: on failure the callback is called right away, and ends the session.
  const std::weak_ptr<McpSessionManager> weak_self = shared_from_this();
  browser->createTarget( url, [weak_self,handle]( const std::string &target_id, const std::string &error ){
    const std::shared_ptr<McpSessionManager> self = weak_self.lock();
    if( !self )
      return;

    if( !error.empty() )
    {
      cerr << "MCP agent server: could not open a browser tab for a session: " << error << endl;
      self->endSession( handle, "the browser could not open a tab for the session" );
      return;
    }

    std::lock_guard<std::mutex> lock( self->m_mutex );
    const auto pos = self->m_sessions.find( handle );
    if( pos == end(self->m_sessions) )
      return;
    pos->second.cdp_target_id = target_id;
    // If the session was ended before its tab came up, close the tab now.
    if( (pos->second.info.state == State::Ended) && self->m_browser )
      self->m_browser->closeTarget( target_id );
  } );

  return handle;
#endif
}//createSession()


std::string McpSessionManager::sessionAppUrl( const std::string &handle )
{
  // Each session's long-polling (server-push) request stays open, and a browser allows only 6
  //  connections per host - so a 7th session would stall until another's poll times out.  Giving
  //  each session its own host name avoids this; "*.localhost" is loopback, and never proxied.
  string host = m_app_host;
  if( (host == "127.0.0.1") || (host == "localhost") )
    host = "s" + std::to_string( ++m_num_sessions_created ) + ".localhost";

  return "http://" + host + ":" + std::to_string( m_app_port ) + "/?externalid=" + handle + "&restore=no";
}//sessionAppUrl(...)


std::optional<McpSessionManager::SessionInfo> McpSessionManager::session( const std::string &handle ) const
{
  std::lock_guard<std::mutex> lock( m_mutex );
  const auto pos = m_sessions.find( handle );
  if( pos == end(m_sessions) )
    return std::nullopt;
  return pos->second.info;
}


void McpSessionManager::endSession( const std::string &handle, const std::string &reason )
{
  vector<std::function<void()>> listeners;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_sessions.find( handle );
    if( pos == end(m_sessions) )
      return;
    listeners = endSessionLocked( pos->second, reason );
  }

  for( const std::function<void()> &fcn : listeners )
    fcn();
}//endSession(...)


std::vector<std::function<void()>> McpSessionManager::endSessionLocked( Session &session,
                                                                         const std::string &reason )
{
  if( session.info.state == State::Ended )
    return {};

  session.info.state = State::Ended;
  session.info.end_reason = reason;
  session.ended = std::chrono::steady_clock::now();

  m_file_tokens.erase( session.info.upload_token );
  for( const string &token : session.download_tokens )
    m_file_tokens.erase( token );
  session.download_tokens.clear();

  InterSpecServer::remove_allowed_session_token( session.info.handle.c_str() );

  if( !session.info.wt_session_id.empty() && !session.app_ended )
  {
    // A delayed schedule() rather than post(): the app's teardown (which can wait seconds for a
    //  fit to stop) would otherwise run on the strand that every session's posted events share.
    Wt::WServer * const server = Wt::WServer::instance();
    if( server )
      server->schedule( std::chrono::milliseconds(1), session.info.wt_session_id, [](){
        if( wApp )
          wApp->quit();
      } );
  }//if( app may still be running )

#if( !defined(_WIN32) )
  if( m_browser && !session.cdp_target_id.empty() )
    m_browser->closeTarget( session.cdp_target_id );
#endif

  // If there was never an app, there is nothing that could still be writing to the directory.
  if( session.info.wt_session_id.empty() || session.app_ended )
    removeSessionDirectory( session );

  cout << "MCP agent server: a session ended (" << reason << ")." << endl;

  vector<std::function<void()>> listeners;
  listeners.swap( session.listeners );
  return listeners;
}//endSessionLocked(...)


void McpSessionManager::touch( const std::string &handle )
{
  std::lock_guard<std::mutex> lock( m_mutex );
  const auto pos = m_sessions.find( handle );
  if( pos != end(m_sessions) )
    pos->second.last_activity = std::chrono::steady_clock::now();
}


void McpSessionManager::callStarted( const std::string &handle )
{
  std::lock_guard<std::mutex> lock( m_mutex );
  const auto pos = m_sessions.find( handle );
  if( pos != end(m_sessions) )
  {
    pos->second.active_calls += 1;
    pos->second.last_activity = std::chrono::steady_clock::now();
  }
}//callStarted(...)


void McpSessionManager::callFinished( const std::string &handle )
{
  std::lock_guard<std::mutex> lock( m_mutex );
  const auto pos = m_sessions.find( handle );
  if( pos != end(m_sessions) )
  {
    pos->second.active_calls = std::max( 0, pos->second.active_calls - 1 );
    pos->second.last_activity = std::chrono::steady_clock::now();
  }
}//callFinished(...)


void McpSessionManager::addSessionListener( const std::string &handle, std::function<void()> fn )
{
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_sessions.find( handle );
    if( (pos != end(m_sessions)) && (pos->second.info.state == State::Starting) )
    {
      pos->second.listeners.push_back( std::move(fn) );
      return;
    }
  }

  fn();
}//addSessionListener(...)


void McpSessionManager::appStarted( const std::string &handle, const std::string &wt_session_id )
{
  vector<std::function<void()>> listeners;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_sessions.find( handle );
    if( pos == end(m_sessions) )
      return;

    Session &session = pos->second;
    if( !session.info.wt_session_id.empty() && (session.info.wt_session_id != wt_session_id) )
      return;

    session.info.wt_session_id = wt_session_id;
    session.last_activity = std::chrono::steady_clock::now();

    if( session.info.state == State::Starting )
    {
      session.info.state = State::Live;
      listeners.swap( session.listeners );
    }else if( session.info.state == State::Ended )
    {
      // Ended while the tab was loading; the app should go too.
      Wt::WServer * const server = Wt::WServer::instance();
      if( server )  //(not post(); see endSessionLocked())
        server->schedule( std::chrono::milliseconds(1), wt_session_id, [](){ if( wApp ) wApp->quit(); } );
    }
  }

  for( const std::function<void()> &fcn : listeners )
    fcn();
}//appStarted(...)


void McpSessionManager::appEnded( const std::string &handle, const std::string &wt_session_id )
{
  vector<std::function<void()>> listeners;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_sessions.find( handle );
    if( pos == end(m_sessions) )
      return;

    Session &session = pos->second;
    if( !session.info.wt_session_id.empty() && (session.info.wt_session_id != wt_session_id) )
      return;
    session.app_ended = true;
    listeners = endSessionLocked( session, "the session's application ended" );
  }

  for( const std::function<void()> &fcn : listeners )
    fcn();

  // The app is still being destroyed (we are called from its destructor), and may write a final
  //  log file; so remove the directory a little later.
  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
  {
    const std::weak_ptr<McpSessionManager> weak_self = shared_from_this();
    server->ioService().schedule( std::chrono::seconds(2), [weak_self,handle](){
      const std::shared_ptr<McpSessionManager> self = weak_self.lock();
      if( !self )
        return;
      std::lock_guard<std::mutex> lock( self->m_mutex );
      const auto pos = self->m_sessions.find( handle );
      if( pos != end(self->m_sessions) )
        self->removeSessionDirectory( pos->second );
    } );
  }//if( server )
}//appEnded(...)


bool McpSessionManager::isManagedHandle( const std::string &handle ) const
{
  std::lock_guard<std::mutex> lock( m_mutex );
  return !handle.empty() && m_sessions.count( handle );
}


std::string McpSessionManager::resolveSessionPath( const std::string &handle, const std::string &relpath ) const
{
  string session_dir;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_sessions.find( handle );
    if( pos == end(m_sessions) )
      throw runtime_error( "Unknown session." );
    if( pos->second.info.state == State::Ended )
      throw runtime_error( "The session has ended (" + pos->second.info.end_reason + ")." );
    session_dir = pos->second.info.dir;
  }

  return resolvePathWithin( session_dir, relpath );
}//resolveSessionPath(...)


std::string McpSessionManager::resolvePathWithin( const std::string &base_dir, const std::string &relpath )
{
  // Deliberately generic: a client probing the file system learns nothing from the message.
  static const char * const invalid_path_msg = "Invalid path: give a path relative to the session"
                                               " directory, like \"input/<filename>\".";

  if( base_dir.empty() || relpath.empty() || (relpath.find('\0') != string::npos) )
    throw runtime_error( invalid_path_msg );

  const fs::path rel( relpath );
  if( rel.is_absolute() || rel.has_root_name() || rel.has_root_directory() )
    throw runtime_error( invalid_path_msg );

  std::error_code ec;
  const fs::path base( base_dir );
  const fs::path full = fs::weakly_canonical( base / rel, ec );
  if( ec || !path_is_within( base, full ) )
    throw runtime_error( invalid_path_msg );

  return full.string();
}//resolvePathWithin(...)


std::string McpSessionManager::sanitizeFilename( const std::string &name )
{
  string filename = name;
  const size_t last_sep = filename.find_last_of( "/\\" );
  if( last_sep != string::npos )
    filename = filename.substr( last_sep + 1 );

  for( char &c : filename )
  {
    const bool ok = ((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z'))
                    || ((c >= '0') && (c <= '9')) || (c == '.') || (c == '_') || (c == '-');
    if( !ok )
      c = '_';
  }

  while( !filename.empty() && (filename.front() == '.') )
    filename.erase( begin(filename) );

  if( filename.size() > 128 )
  {
    const size_t dot = filename.find_last_of( '.' );
    const string ext = ((dot != string::npos) && ((filename.size() - dot) <= 16)) ? filename.substr(dot) : string();
    filename = filename.substr( 0, 128 - ext.size() ) + ext;
  }

  if( filename.empty() || (filename.find_first_not_of( "._-" ) == string::npos) )
    filename = "file";

  return filename;
}//sanitizeFilename(...)


std::string McpSessionManager::uniqueFilename( const std::string &directory, const std::string &filename )
{
  const string name = sanitizeFilename( filename );

  std::error_code ec;
  if( !fs::exists( fs::symlink_status( fs::path(directory) / name, ec ) ) )
    return name;

  const size_t dot = name.find_last_of( '.' );
  const string stem = (dot == string::npos) ? name : name.substr( 0, dot );
  const string ext = (dot == string::npos) ? string() : name.substr( dot );

  for( int i = 1; i < 10000; ++i )
  {
    const string candidate = stem + "-" + std::to_string(i) + ext;
    if( !fs::exists( fs::symlink_status( fs::path(directory) / candidate, ec ) ) )
      return candidate;
  }

  return stem + "-" + Wt::WRandom::generateId(8) + ext;
}//uniqueFilename(...)


uint64_t McpSessionManager::sessionBytesUsed( const std::string &handle ) const
{
  string dir;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_sessions.find( handle );
    if( pos == end(m_sessions) )
      return 0;
    dir = pos->second.info.dir;
  }

  uint64_t total = 0;
  std::error_code ec;
  for( fs::recursive_directory_iterator iter( dir, ec ), iter_end; !ec && (iter != iter_end); iter.increment(ec) )
  {
    if( iter->is_regular_file( ec ) && !iter->is_symlink( ec ) )
      total += static_cast<uint64_t>( iter->file_size( ec ) );
  }

  return total;
}//sessionBytesUsed(...)


std::string McpSessionManager::addDownload( const std::string &handle, const std::string &relpath )
{
  std::lock_guard<std::mutex> lock( m_mutex );
  const auto pos = m_sessions.find( handle );
  if( (pos == end(m_sessions)) || (pos->second.info.state == State::Ended) )
    throw runtime_error( "The session has ended." );

  for( const string &token : pos->second.download_tokens )
  {
    const auto token_pos = m_file_tokens.find( token );
    if( (token_pos != end(m_file_tokens)) && (token_pos->second.relpath == relpath) )
      return token;
  }

  string token;
  do
  {
    token = Wt::WRandom::generateId( 32 );
  }while( m_file_tokens.count( token ) );

  m_file_tokens[token] = FileToken{ handle, false, relpath };
  pos->second.download_tokens.insert( token );

  return token;
}//addDownload(...)


std::optional<McpSessionManager::FileToken> McpSessionManager::lookupFileToken( const std::string &token ) const
{
  std::lock_guard<std::mutex> lock( m_mutex );
  const auto pos = m_file_tokens.find( token );
  if( pos == end(m_file_tokens) )
    return std::nullopt;
  return pos->second;
}


McpSessionManager::HealthInfo McpSessionManager::health() const
{
  HealthInfo info;
  std::lock_guard<std::mutex> lock( m_mutex );
  info.max_sessions = m_settings.maxSessions;
  for( const auto &handle_session : m_sessions )
  {
    info.live_sessions += (handle_session.second.info.state == State::Live);
    info.starting_sessions += (handle_session.second.info.state == State::Starting);
  }
#if( !defined(_WIN32) )
  info.browser_running = (m_browser && m_browser->running());
#endif
  return info;
}//health()


std::string McpSessionManager::currentHandle()
{
  InterSpecApp * const app = dynamic_cast<InterSpecApp *>( Wt::WApplication::instance() );
  if( !app )
    return string();

  const std::shared_ptr<McpSessionManager> manager = instance();
  const string &token = app->externalToken();
  return (manager && manager->isManagedHandle( token )) ? token : string();
}//currentHandle()


void McpSessionManager::sweepOrphanedDirectories()
{
  if( m_settings.keepSessionFiles )
    return;

  size_t num_removed = 0;
  std::error_code ec;
  for( fs::directory_iterator iter( m_exchange_root, ec ), iter_end; !ec && (iter != iter_end); iter.increment(ec) )
  {
    const fs::path path = iter->path();
    std::error_code status_ec;
    if( !fs::is_directory( fs::symlink_status( path, status_ec ) ) )
      continue;
    if( !is_session_dir_name( path.filename().string() ) )
      continue;
    if( !fs::is_regular_file( fs::symlink_status( path / sm_marker_filename, status_ec ) ) )
      continue;

    std::error_code remove_ec;
    fs::remove_all( path, remove_ec );
    num_removed += !remove_ec;
  }//for( loop over exchange directory )

  if( num_removed )
    cout << "MCP agent server: removed " << num_removed << " session directories left by an earlier run." << endl;
}//sweepOrphanedDirectories()


void McpSessionManager::removeSessionDirectory( Session &session )
{
  if( session.dir_removed )
    return;
  session.dir_removed = true;

  if( m_settings.keepSessionFiles || session.info.dir.empty() )
    return;

  // Only ever delete something we created: a directly-contained, session-named directory that
  //  holds our marker file.
  const fs::path dir( session.info.dir );
  std::error_code ec;
  if( (dir.parent_path() != fs::path(m_exchange_root))
      || (dir.filename().string() != session.info.dir_name)
      || !is_session_dir_name( session.info.dir_name )
      || !fs::is_directory( fs::symlink_status( dir, ec ) )
      || !fs::is_regular_file( fs::symlink_status( dir / sm_marker_filename, ec ) ) )
  {
    cerr << "MCP agent server: not removing session directory '" << session.info.dir
         << "', as it does not look like one we created." << endl;
    return;
  }

  fs::remove_all( dir, ec );
  if( ec )
    cerr << "MCP agent server: error removing session directory: " << ec.message() << endl;
}//removeSessionDirectory(...)


void McpSessionManager::scheduleReaper()
{
  Wt::WServer * const server = Wt::WServer::instance();
  if( !server )
    return;

  const std::weak_ptr<McpSessionManager> weak_self = shared_from_this();
  server->ioService().schedule( sm_reaper_interval, [weak_self](){
    const std::shared_ptr<McpSessionManager> self = weak_self.lock();
    if( !self )
      return;
    {
      std::lock_guard<std::mutex> lock( self->m_mutex );
      if( !self->m_running )
        return;
    }
    self->reap();
    self->scheduleReaper();
  } );
}//scheduleReaper()


void McpSessionManager::reap()
{
  const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  const std::chrono::seconds idle_timeout( m_settings.sessionIdleTimeoutSeconds );
  const std::chrono::seconds max_lifetime( m_settings.sessionMaxLifetimeSeconds );

  vector<std::function<void()>> listeners;
  std::shared_ptr<CdpBrowser> dead_browser;
  {
    std::lock_guard<std::mutex> lock( m_mutex );

    for( auto iter = begin(m_sessions); iter != end(m_sessions); )
    {
      Session &session = iter->second;
      string reason;

      switch( session.info.state )
      {
        case State::Starting:
          if( (now - session.created) > sm_start_timeout )
            reason = "the session did not start in time";
          break;

        case State::Live:
          if( (m_settings.sessionIdleTimeoutSeconds > 0) && (session.active_calls == 0)
              && ((now - session.last_activity) > idle_timeout) )
            reason = "the session was idle for too long";
          break;

        case State::Ended:
          if( !session.dir_removed && ((now - session.ended) > sm_dir_removal_grace) )
            removeSessionDirectory( session );
          break;
      }//switch( session.info.state )

      if( reason.empty() && (session.info.state != State::Ended)
          && (m_settings.sessionMaxLifetimeSeconds > 0) && ((now - session.created) > max_lifetime) )
        reason = "the session reached its maximum lifetime";

      if( !reason.empty() )
      {
        vector<std::function<void()>> these = endSessionLocked( session, reason );
        listeners.insert( end(listeners), begin(these), end(these) );
      }

      if( (session.info.state == State::Ended) && session.dir_removed
          && ((now - session.ended) > sm_ended_record_lifetime) )
        iter = m_sessions.erase( iter );
      else
        ++iter;
    }//for( loop over sessions )

#if( !defined(_WIN32) )
    // Dropped here, not on the browser's own reader thread (whose exit callback flagged it).
    if( m_browser && !m_browser->running() )
      dead_browser.swap( m_browser );
#endif
  }

  for( const std::function<void()> &fcn : listeners )
    fcn();

  dead_browser.reset();
}//reap()


void McpSessionManager::handleTargetGone( const std::string &target_id )
{
  string handle;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    for( const auto &handle_session : m_sessions )
    {
      if( handle_session.second.cdp_target_id == target_id )
        handle = handle_session.first;
    }
  }

  if( !handle.empty() )
    endSession( handle, "the session's browser tab closed" );
}//handleTargetGone(...)


void McpSessionManager::handleBrowserExited( const uint64_t browser_generation )
{
  cerr << "MCP agent server: the browser exited; ending its sessions." << endl;

  vector<std::function<void()>> listeners;
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    for( auto &handle_session : m_sessions )
    {
      // (A replacement browser may already have sessions.)
      if( handle_session.second.browser_generation != browser_generation )
        continue;

      vector<std::function<void()>> these = endSessionLocked( handle_session.second, "the browser exited" );
      listeners.insert( end(listeners), begin(these), end(these) );
    }
  }

  for( const std::function<void()> &fcn : listeners )
    fcn();
}//handleBrowserExited()
