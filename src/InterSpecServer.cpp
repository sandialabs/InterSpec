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

#include <map>
#include <mutex>
#include <string>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>


#if __APPLE__
#include "TargetConditionals.h"
#endif
#include <Wt/WApplication>
#include <Wt/WServer>
#include <Wt/WResource>
#include <Wt/WIOService>
#include <Wt/Http/Client>
#include <Wt/Utils>
#include <Wt/Http/Response>
#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/WtSqlTraits>
#include <Wt/Dbo/backend/Sqlite3>
#include <Wt/Json/Array>
#include <Wt/Json/Object>
#include <Wt/Json/Parser>


#include <boost/filesystem.hpp>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "SpecUtils/SerialToDetectorModel.h"
#include "InterSpec/DataBaseVersionUpgrade.h"


using namespace std;

// Namespace to track allowed sessions and not
namespace
{
struct SessionState
{
  enum class State
  {
    Invalid             = 0,
    AuthorizedNotLoaded = 1,
    AuthorizedLoaded    = 2,
    NoLongerAuthorized  = 3,
    Dead                = 4
  };//enum class State
  
  State current_state;
  
  InterSpecServer::SessionType session_type;
  
  /** Either a file-system path, or app-url to open when the session is first constructed. */
  std::string initial_file_to_open;
  
  std::chrono::system_clock::time_point auth_time;
  std::chrono::system_clock::time_point load_time;
  std::chrono::system_clock::time_point deauth_time;
  std::chrono::system_clock::time_point destruct_time;
  
  SessionState()
  : current_state( State::Invalid ),
    auth_time(),
    session_type(InterSpecServer::SessionType::ExternalBrowserInstance),
    load_time(),
    deauth_time(),
    destruct_time()
  {
  }
};//struct SessionState


std::mutex ns_sessions_mutex;
std::map<string,SessionState> ns_sessions;
bool ns_allow_untokened_sessions = true;
  
  Wt::WApplication *create_application( const Wt::WEnvironment &env )
  {
    return new InterSpecApp( env );
  }// Wt::WApplication *createApplication(const Wt::WEnvironment& env)
}

namespace InterSpecServer
{
  int sm_portServedOn = -1;
  std::string sm_urlServedOn = "";
  std::mutex sm_servedOnMutex;
  
  Wt::WServer *ns_server = nullptr;
  std::mutex ns_servermutex;
  
  
  int portBeingServedOn()
  {
    std::lock_guard<std::mutex> lock( sm_servedOnMutex );
    return sm_portServedOn;
  }
  
  std::string urlBeingServedOn()
  {
    std::lock_guard<std::mutex> lock( sm_servedOnMutex );
    return sm_urlServedOn;
  }
  
  
  std::string getWtConfigXml( int argc, char *argv[] )
  {
    for( int i = 0; i < (argc-1); ++i )
    {
      if( strcmp(argv[i],"-c") == 0 )
        return argv[i+1];
      if( strcmp(argv[i+1],"--config") == 0 )
        return argv[i];
    }//for( int i = 0; i < argc; ++i )
    
#if( BUILD_OSX_APP )
    if( boost::filesystem::exists( "data/config/wt_config_osx.xml" ) )
      return (boost::filesystem::current_path() / boost::filesystem::path("data/config/wt_config_osx.xml")).string<std::string>();
#endif
    
#if( ANDROID )
    if( boost::filesystem::exists( "data/config/wt_config_android.xml" ) )
      return (boost::filesystem::current_path() / boost::filesystem::path("data/config/wt_config_android.xml")).string<std::string>();
#endif
    
#if( IOS )
    if( boost::filesystem::exists( "data/config/wt_config_ios.xml" ) )
      return (boost::filesystem::current_path() / boost::filesystem::path("data/config/wt_config_ios.xml")).string<std::string>();
#endif
    
    if( boost::filesystem::exists( "wt_config.xml" ) )
      return (boost::filesystem::current_path() / "wt_config.xml").string<std::string>();
    
    //should search the "PATH" here
    
    return WT_CONFIG_XML;
  }//std::string getWtConfigXml(int argc, char *argv[])
  
  
  
  bool changeToBaseDir( int argc, char *argv[] )
  {
    for( int i = 1; i < argc; ++i )
    {
      string arg = argv[i];
      
      if( SpecUtils::starts_with(arg, "--basedir=") )
      {
        std::string workdir = arg.substr( 10 );
        if( workdir.size() && (workdir[0]=='\"' || workdir[0]=='\'') )
          workdir = workdir.substr(1);
        if( workdir.size()
           && (workdir[workdir.size()-1]=='\"' || workdir[workdir.size()-1]=='\'') )
          workdir = workdir.substr(0, workdir.size()-1 );
        
        try
        {
          boost::filesystem::current_path( workdir );
          cout << "Changed to cwd: " << workdir << endl;
          return true;
        }catch(...)
        {
          cerr << "Unable to change to directory: '" << workdir << "'" << endl;
          return false;
        }
      }
    }//for( int i = 0; i < argc; ++i )
    
    return false;
  }//bool changeToBaseDir( argc, argv )
  
  
  Wt::WApplication *createAppForServer( const Wt::WEnvironment &env,
                                       Wt::WApplication::ApplicationCreator appCreator )
  {
    Wt::WApplication *app = appCreator( env );
    return app;
  }//WApplication *createAppForServer(...)
  
  
  void startServer( int argc, char *argv[],
                                Wt::WApplication::ApplicationCreator createApplication )
  {
    changeToBaseDir( argc, argv );
    const string xml_config_path = getWtConfigXml( argc, argv );
    
    std::lock_guard<std::mutex> serverlock( ns_servermutex );
    if( ns_server )
    {
      std::cerr << "experimental_startServer: already running" << std::endl;
      return;
    }

    Wt::WString::setDefaultEncoding( Wt::UTF8 );   

    //If we are in an Apple Sandbox, we cant write to /tmp (especially done when
    //  spooling files), so we will check if we passed in an argument of a tmp
    //  directory we can write to.
    for( int i = 0; i < (argc-1); ++i )
    {
      if( strcmp(argv[i],"--tempdir") == 0 )
      {
#ifdef WIN32
		  _putenv_s("TMPDIR", argv[i + 1]);
		  _putenv_s("WT_TMP_DIR", argv[i + 1]);
#else
		  setenv("TMPDIR", argv[i + 1], 1);
		  setenv("WT_TMP_DIR", argv[i + 1], 1);
#endif // WIN32
        cerr << "Setting tmpdir=" << argv[i+1] << endl;
        break;
      }
    }//for( int i = 0; i < argc; ++i )
    
    
    ns_server = new Wt::WServer( argv[0], xml_config_path );
    char httpaddr_param_name[]  = "--http-addr";
    char httpaddr_param_value[] = "127.0.0.1";
    char httpport_param_name[]  = "--http-port";
    char httpport_param_value[] = "0";           //Assign port automatically
    char docroot_param_name[]   = "--docroot";
    char docroot_param_value[]  = ".";
    char accesslog_param_value[]  = "--accesslog=-";  //quite down printing all the GET and POST and such
    
    char *argv_wthttp[] = { argv[0],
      httpaddr_param_name, httpaddr_param_value,
      httpport_param_name, httpport_param_value,
      docroot_param_name, docroot_param_value,
      accesslog_param_value
    };
    const int argc_wthttp = sizeof(argv_wthttp)/sizeof(argv_wthttp[0]);
    
    ns_server->setServerConfiguration( argc_wthttp, argv_wthttp, WTHTTP_CONFIGURATION );
    
    ns_server->addEntryPoint( Wt::Application, boost::bind( &createAppForServer, _1, createApplication ) );
    
    if( ns_server->start() )
    {
      // Start initializing DecayDataBaseServer.
      //  Previous to 20230203, we did this in the beginning of InterSpecApp constructor.
      //  On a 2019 macBook pro, release build of the native app, it took about 485 ms from the
      //  start of starting the server, to when the decay database was done initializing.
      //  The first request for the database had to wait about 120 ms for it to be ready; I
      //  think this was another background thread (the InterSpec::fillMaterialDb() function),
      //  but the second request for the database was a GUI render thread, which had to wait
      //  40 to 80 ms.
      //  If we do it here, its 260 ms of the equivalent timespan, and the GUI thread never has
      //  to wait to access the database.
      //  Using the minimized coincidence version of sandia.decay.xml increases parse time
      //  by about 170 ms.
      ns_server->ioService().boost::asio::io_service::post( &DecayDataBaseServer::initialize );
      
      const int port = ns_server->httpPort();
      std::string this_url = "http://127.0.0.1:" + boost::lexical_cast<string>(port);
      
      {
        std::lock_guard<std::mutex> lock( sm_servedOnMutex );
        sm_portServedOn = port;
        sm_urlServedOn = this_url;
      }
    }else
    {
      throw std::runtime_error( "Failed to start Wt server" );
    }//if( server.start() )
  }//void startServer()
  

  
  void startServerNodeAddon( string name,
                             std::string basedir,
                             const std::string xml_config_path,
                             unsigned short int server_port_num )
  {
    std::lock_guard<std::mutex> serverlock( ns_servermutex );
    if( ns_server )
    {
      std::cerr << "experimental_startServer: already running" << std::endl;
      return;
    }
    
    if( name.empty() )
      name = "\0";
    if( basedir.empty() )
      basedir = ".";
    
    Wt::WString::setDefaultEncoding( Wt::UTF8 );
    
    
    ns_server = new Wt::WServer( name, xml_config_path );
    char *exe_param_name  = &(name[0]);
    char httpaddr_param_name[]  = "--http-addr";
    char httpaddr_param_value[] = "127.0.0.1";
    char httpport_param_name[]  = "--http-port";
    string port_str = std::to_string( static_cast<int>(server_port_num) );
    assert( !port_str.empty() );
    if( port_str.empty() )
      port_str = "0";
    char *httpport_param_value = &(port_str[0]);
    //char httpport_param_value[] = "0";           //Assign port automatically
    char docroot_param_name[]   = "--docroot";
    char *docroot_param_value  = &(basedir[0]);
    //char approot_param_name[]   = "--approot";
    //char *approot_param_value  = &(basedir[0]);
    char accesslog_param_value[]  = "--accesslog=-";  //quite down printing all the GET and POST and such

    
    char *argv_wthttp[] = { exe_param_name,
      httpaddr_param_name, httpaddr_param_value,
      httpport_param_name, httpport_param_value,
      docroot_param_name, docroot_param_value,
      accesslog_param_value
    };
    const int argc_wthttp = sizeof(argv_wthttp)/sizeof(argv_wthttp[0]);
    
    ns_server->setServerConfiguration( argc_wthttp, argv_wthttp, WTHTTP_CONFIGURATION );
    
    ns_server->addEntryPoint( Wt::Application, boost::bind( &createAppForServer, _1, create_application ) );
    
    if( ns_server->start() )
    {
      // See remarks in startServer() on performance and reason for this next call
      ns_server->ioService().boost::asio::io_service::post( &DecayDataBaseServer::initialize );
      
      const int port = ns_server->httpPort();
      assert( !server_port_num || (server_port_num == port) );
      std::string this_url = "http://127.0.0.1:" + boost::lexical_cast<string>(port);
      
      {
        std::lock_guard<std::mutex> lock( sm_servedOnMutex );
        sm_portServedOn = port;
        sm_urlServedOn = this_url;
      }
    }else
    {
      throw std::runtime_error( "Failed to start Wt server" );
    }//if( server.start() )
  }//startServerNodeAddon(...)
  

int start_server( const char *process_name, const char *userdatadir,
                  const char *basedir, const char *xml_config_path, 
                  unsigned short int server_port )
{
  //Using a relative path should get us in less trouble than an absolute path
  //  on Windows.  Although havent yet tested (20190902) with network drives and such on Windows.
  //Should check if basedir is a relative or absolut path before this next step
//#warning "Need to check out how this setting basedir for Electron target works on Windows network drives and such"
  string cwd, relbasedir;

  try
  {
    cwd = SpecUtils::get_working_path();
    relbasedir = basedir; //SpecUtils::fs_relative( cwd, basedir );
    cerr << "cwd='" << cwd << "'" << endl;
    cerr << "relbasedir='" << relbasedir << "'" << endl;
    cerr << "userdatadir='" << userdatadir << "'" << endl;
    //if( relbasedir.size() >= strlen(basedir) )
    //  relbasedir = basedir;
    if( relbasedir.empty() )
      relbasedir = ".";
  }catch( std::exception &e )
  {
    cerr << "When setting directories to serve, caught: " << e.what() << endl;
    return -1;
  }
  
  try
  {
    //ToDo: refactor this userdatadir stuff into function inside InterSpecServer
    //      or InterSpecApp or something.
    if( !SpecUtils::create_directory(userdatadir) )
      throw std::runtime_error( "Failed to create directory '" + string(userdatadir) + "' for user data." );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -2;
  }
  
  try
  {
    const string preffile = SpecUtils::append_path( userdatadir, "InterSpecUserData.db" );
    
    cout << "Will set user preferences file to: '" << preffile << "'" << endl;
    DataBaseUtils::setPreferenceDatabaseFile( preffile );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -3;
  }
  
  try
  {
    const string directories_to_try[] = { userdatadir, 
      SpecUtils::append_path(relbasedir,"data"), 
      SpecUtils::append_path(relbasedir,"data_ouo") 
    };

    for( const string &trial_dir : directories_to_try )
    {
      const auto serial_db = SpecUtils::ls_files_in_directory( trial_dir, "serial_to_model.csv" );
      if( !serial_db.empty() )
      {
        SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
        cout << "Will use serial_to_model.csv from directory '" << trial_dir << "'" << endl;
        break;
      }
    }//
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -4;
  }
  
  
  try
  {
    InterSpec::setWritableDataDirectory( userdatadir );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -5;
  }
  
  try
  {
    cout << "Will make sure preferences database is up to date" << endl;
    DataBaseVersionUpgrade::checkAndUpgradeVersion();
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -6;
  }
  
  
  try
  {
    InterSpec::setStaticDataDirectory( SpecUtils::append_path(relbasedir,"data") );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -8;
  }

  try
  {
    const string user_decay = SpecUtils::append_path(userdatadir, "sandia.decay.xml");
    if( SpecUtils::is_file(user_decay) )
      DecayDataBaseServer::setDecayXmlFile(user_decay);
  }catch (std::exception& e)
  {
    cerr << e.what() << endl;
    return -9;
  }


  try
  {
    const string user_reaction = SpecUtils::append_path(userdatadir, "sandia.reactiongamma.xml");
    if( SpecUtils::is_file(user_reaction) )
      ReactionGammaServer::set_xml_file_location(user_reaction);
  }catch (std::exception& e)
  {
    cerr << e.what() << endl;
    return -10;
  }
  

    //ToDo: should look into using '--approot' Wt Argument.
    
    //try
    //{
    //  boost::filesystem::current_path( basedir );
    //  cout << "Changed to cwd: " << basedir << endl;
    //}catch ( std::exception &e )
    //{
    //  cerr << "Unable to change to directory: '" << basedir << "' :" << e.what() << endl;
    //}

  try
  {
    InterSpecServer::startServerNodeAddon( process_name, relbasedir, xml_config_path, server_port );
  }catch( std::exception &e )
  {
    std::cerr << "\n\nCaught exception trying to start InterSpec server:\n\t"
    << e.what() << std::endl << std::endl;
    return -11;
  }
  
  return InterSpecServer::portBeingServedOn();
}//int interspec_start_server( int argc, char *argv[] )
  
  void killServer()
  {
    std::lock_guard<std::mutex> serverlock( ns_servermutex );
    
    if( ns_server )
    {
      std::cerr << "About to stop server" << std::endl;
      ns_server->stop();
      delete ns_server;
      ns_server = nullptr;
      std::cerr << "Stopped and killed server" << std::endl;
    }
  }//void killServer()
  

  int wait_for_shutdown()
  {
    {
      std::lock_guard<std::mutex> serverlock( ns_servermutex );
      
      if( !ns_server )
        return -1;
    }
    
    return Wt::WServer::waitForShutdown();
  }//int wait_for_shutdown()

  
  int add_allowed_session_token( const char *session_id, const SessionType session_type )
  {
    //Returns zero if hadnt been seen before, 1 if authorized but not seen yet, 2 if authorized and loaded, 3 if deauthorized session, 4 if dead session
    SessionState newsession;
    newsession.session_type = session_type;
    newsession.auth_time = std::chrono::system_clock::now();
    newsession.current_state = SessionState::State::AuthorizedNotLoaded;
    
    lock_guard<mutex> lock( ns_sessions_mutex );
    
    const auto pos = ns_sessions.find( session_id );
    if( pos != end(ns_sessions) )
    {
      switch( pos->second.current_state )
      {
        case SessionState::State::Invalid:             return -1;  //shouldnt happen
        case SessionState::State::AuthorizedNotLoaded: return 1;
        case SessionState::State::AuthorizedLoaded:    return 2;
        case SessionState::State::NoLongerAuthorized:  return 3;
        case SessionState::State::Dead:                return 4;
      }//
      assert( 0 );
    }//if( we already have this token )
    
    ns_sessions[session_id] = newsession;
    
    // Do some very niave cleanup for super long running sessions, will probably never actually hit
    if( ns_sessions.size() > 500 )
    {
      vector<string> to_remove;
      for( auto &keyval : ns_sessions )
      {
        if( keyval.second.current_state == SessionState::State::Dead )
          to_remove.push_back( keyval.first );
      }
      
      for( const auto &key : to_remove )
        ns_sessions.erase( key );
    }//if( ns_sessions.size() > 500 )
    
    return 0;
  }//void add_allowed_session_token( const char *session_id )
  
  
  int remove_allowed_session_token( const char *session_token )
  {
    lock_guard<mutex> lock( ns_sessions_mutex );
    auto pos = ns_sessions.find( session_token );
    if( pos == end(ns_sessions) )
      return -1;
    
    switch( pos->second.current_state )
    {
      case SessionState::State::Invalid:
        return -2;
        
      case SessionState::State::AuthorizedNotLoaded:
        pos->second.current_state = SessionState::State::NoLongerAuthorized;
        pos->second.deauth_time = std::chrono::system_clock::now();
        return 1;
        
      case SessionState::State::AuthorizedLoaded:
      case SessionState::State::Dead:
        pos->second.current_state = SessionState::State::NoLongerAuthorized;
        pos->second.deauth_time = std::chrono::system_clock::now();
        return 0;
        
      case SessionState::State::NoLongerAuthorized:
        pos->second.current_state = SessionState::State::NoLongerAuthorized;
        pos->second.deauth_time = std::chrono::system_clock::now();
        return 2;
    }//switch( pos->current_state )
    
    assert( 0 );
    
    return -1;
  }//int interspec_remove_allowed_session_token( const char *session_token )
  

  bool allow_untokened_sessions()
  {
    lock_guard<mutex> lock( ns_sessions_mutex );
    return ns_allow_untokened_sessions;
  }//bool allow_untokened_sessions()


  void set_require_tokened_sessions( const bool require )
  {
    lock_guard<mutex> lock( ns_sessions_mutex );
    ns_allow_untokened_sessions = !require;
  }


  int set_session_loaded( const char *session_token )
  {
    //Returns 0 if had been authorized, 1 if had been seen, and 2 if not authorized, dead
    
    lock_guard<mutex> lock( ns_sessions_mutex );
    auto pos = ns_sessions.find( session_token );
    if( pos == end(ns_sessions) )
      return 2;
    
    switch( pos->second.current_state )
    {
      case SessionState::State::Invalid:
        return 3;
        
      case SessionState::State::AuthorizedNotLoaded:
        pos->second.current_state = SessionState::State::AuthorizedLoaded;
        pos->second.load_time = std::chrono::system_clock::now();
        return 0;
        
      case SessionState::State::AuthorizedLoaded:
        return 1;
        
      case SessionState::State::Dead:
      case SessionState::State::NoLongerAuthorized:
        return 2;
    }//switch( pos->current_state )
    
    assert( 0 );
    
    return -1;
  }//int set_session_loaded( const char *session_token )

  
  std::pair<bool,SessionType> session_type( const char *session_token )
  {
    lock_guard<mutex> lock( ns_sessions_mutex );
    auto pos = ns_sessions.find( session_token );
    if( pos == end(ns_sessions) )
      return std::pair<bool,SessionType>( false, SessionType::ExternalBrowserInstance );
    
    return std::pair<bool,SessionType>( true, pos->second.session_type );
  }

  
  void set_session_destructing( const char *session_token )
  {
    lock_guard<mutex> lock( ns_sessions_mutex );
    auto pos = ns_sessions.find( session_token );
    if( pos != end(ns_sessions) )
    {
      pos->second.current_state = SessionState::State::Dead;
      pos->second.destruct_time = std::chrono::system_clock::now();
    }
  }//void set_session_destructing( const char *session_token )


  int session_status( const char *session_token )
  {
    lock_guard<mutex> lock( ns_sessions_mutex );
    auto pos = ns_sessions.find( session_token );
    if( pos == end(ns_sessions) )
      return 0;
    
    static_assert( static_cast<int>(SessionState::State::AuthorizedNotLoaded) == 1, "" );
    static_assert( static_cast<int>(SessionState::State::AuthorizedLoaded) == 2, "" );
    static_assert( static_cast<int>(SessionState::State::NoLongerAuthorized) == 3, "" );
    static_assert( static_cast<int>(SessionState::State::Dead) == 4, "" );
    
    return static_cast<int>( pos->second.current_state );
  }//int session_status( const char *session_token )


  int open_file_in_session( const char *sessionToken, const char *files_json )
  {
#ifndef _WIN32
  #warning "Need to actually test interspec_open_file"
#endif    
    //Are we guaranteed to receeve the entire message at once?
    cerr << "Opening files not tested!" << endl;
    
    vector<string> files, appurls;
    
    try
    {
      Wt::Json::Value result;
      Wt::Json::parse( files_json, result, true );
      
      if( result.type() != Wt::Json::ArrayType )
        throw runtime_error( "Json passed in was not an array" );
      
      const Wt::Json::Array &jsonfiles = result;
      for( const auto &val : jsonfiles )
      {
        const string valstr = val.orIfNull("");
        if( SpecUtils::is_file( valstr ) )
        {
          files.push_back( valstr );
        }else if( SpecUtils::istarts_with( valstr, "interspec://" ) )
        {
          appurls.push_back( valstr );
        }else if( SpecUtils::istarts_with( valstr, "raddata://g" ) )
        {
          // TODO: should try to detect if potentually url-encoded to help decide things; I'm really not sure if how we get URLs is consistent on the various platforms, or maybe I amd the one causing this problem on macOS
#if __APPLE__
          appurls.push_back( valstr );
#else
          appurls.push_back( Wt::Utils::urlDecode( valstr ) );
#endif
        }else
        {
          cerr << "File '" << valstr << "' is not a file" << endl;
        }
      }//for( const auto &val : jsonfiles )
    }catch( std::exception &e )
    {
      cerr << "Failed to parse '" << files_json
      << "' as valid JSON. Issue: " << e.what() << endl;
      return -1;
    }//try / catch
    
    
    //Look for session with 'externalid' and open file...
    string externalid = sessionToken;
    
    InterSpecApp *app = InterSpecApp::instanceFromExtenalToken( externalid );
    if( app )
    {
      int numopened = 0;
      Wt::WApplication::UpdateLock applock( app );
      
      for( const string &filename : files )
      {
        if( app->userOpenFromFileSystem( filename ) )
          numopened += 1;
        else
          cerr << "InterSpec failed to open file: '" << filename << "'" << endl;
      }
      
      for( const string &url : appurls )
      {
        if( app->handleAppUrl( url ) )
          numopened += 1;
        else
          cerr << "InterSpec failed to open app-url: '" << url << "'"  << endl;
      }

      app->triggerUpdate();
      
      return numopened;
    }else
    {
      //I dont know why we would get here... but lets deal with it JIC
      cerr << "There is no app with externalid=" << externalid << endl;
      
#if( IOS )
      return -2;
#endif
      
      Wt::WServer *server = Wt::WServer::instance();
      if( server )
      {
        cerr << "Will ask ALL current sesssions to open file." << endl;
        
        server->postAll( std::bind( [files](){
          Wt::WApplication *wtap = wApp;
          if( !wtap )
          {
            //Not sure why this happens some times.
            cerr << "No WApplication::instance() in postAll(...)" << endl;
            return;
          }
          
          InterSpecApp *app = dynamic_cast<InterSpecApp *>( wtap );
          assert( app );
          
          InterSpec *interspec = app->viewer();
          assert( interspec );
          
          for( auto filename : files )
          {
            if( !app->userOpenFromFileSystem( filename ) )
              cerr << "InterSpec failed to open file filename" << endl;
          }
          
          app->triggerUpdate();
        }) );
      }else
      {
        cerr << "There is no server running, not opening file." << endl;
      }
    }//if( app ) / else
    
    return -2;
  }

  bool pass_app_url_to_session( const char *session_token, const std::string &url )
  {
    bool used = false;
    InterSpecApp *app = InterSpecApp::instanceFromExtenalToken( session_token );
    if( app )
    {
      // TODO: need to figure out which platforms go through here, and if they are encoded.
#if( ANDROID )
      const std::string unencoded = Wt::Utils::urlDecode(url);
      const std::string &unecodedUrl = unencoded;
#else
      const std::string &unecodedUrl = url;
#endif

      Wt::WApplication::UpdateLock applock( app );
      used = app->handleAppUrl( unecodedUrl );
      app->triggerUpdate();
    }else
    {
      cerr << "Failed to find session with token '" << session_token << "'" << endl;
    }

    return used;
  }//pass_app_url_to_session(....)


void set_file_to_open_on_load( const char *session_token, const std::string file_path )
{
  lock_guard<mutex> lock( ns_sessions_mutex );
  auto pos = ns_sessions.find( session_token );
  if( pos == end(ns_sessions) )
    throw runtime_error( "set_file_to_open_on_load: specified token is not present." );
  
  if( pos->second.current_state != SessionState::State::AuthorizedNotLoaded )
  {
    assert( 0 );
    throw runtime_error( "set_file_to_open_on_load: specified session has already been loaded" );
  }
  
  pos->second.initial_file_to_open = file_path;
}//void set_file_to_open_on_load(...)


std::string file_to_open_on_load( const std::string &session_token )
{
  lock_guard<mutex> lock( ns_sessions_mutex );
  auto pos = ns_sessions.find( session_token );
  if( pos == end(ns_sessions) )
  {
    assert( 0 );
    return "";
  }
  
  return pos->second.initial_file_to_open;
}//file_to_open_on_load(...)

}//namespace InterSpecServer

