//
//  wtServer.cpp
//  InterSpec
//
//  Created by Johnson, William C on 6/25/17.
//  Copyright © 2017 Sandia National Laboratories. All rights reserved.
//

#include "InterSpec_config.h"

#include "InterSpec/InterSpecServer.h"

#include <map>
#include <mutex>
#include <string>
#include <cstdlib>
#include <stdlib.h>

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


#include <stdlib.h>
#include <iostream>
#include <fstream>


#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>


#include "InterSpec/InterSpecApp.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DbToFilesystemLink.h"


using namespace std;

namespace
{
  //The externalid of the primary session (the main electron window)
  string ns_externalid;
  std::mutex ns_externalid_mutex;

  
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
  
  Wt::WServer *ns_server = 0;
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
      
      if( UtilityFunctions::starts_with(arg, "--basedir=") )
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
//#warning "Need to add a (optional) number that the requestor must provide in the URL in order to be served anything - this should probably be the externalid argument"
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
    
    char *argv_wthttp[] = { argv[0],
      httpaddr_param_name, httpaddr_param_value,
      httpport_param_name, httpport_param_value,
      docroot_param_name, docroot_param_value
    };
    const int argc_wthttp = sizeof(argv_wthttp)/sizeof(argv_wthttp[0]);
    
    ns_server->setServerConfiguration( argc_wthttp, argv_wthttp, WTHTTP_CONFIGURATION );
    
    ns_server->addEntryPoint( Wt::Application, boost::bind( &createAppForServer, _1, createApplication ) );
    
    if( ns_server->start() )
    {
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
                             const std::string xml_config_path )
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
    char httpport_param_value[] = "0";           //Assign port automatically
    char docroot_param_name[]   = "--docroot";
    char *docroot_param_value  = &(basedir[0]);
    //char approot_param_name[]   = "--approot";
    //char *approot_param_value  = &(basedir[0]);
    
    
    char *argv_wthttp[] = { exe_param_name,
      httpaddr_param_name, httpaddr_param_value,
      httpport_param_name, httpport_param_value,
      docroot_param_name, docroot_param_value
    };
    const int argc_wthttp = sizeof(argv_wthttp)/sizeof(argv_wthttp[0]);
    
    ns_server->setServerConfiguration( argc_wthttp, argv_wthttp, WTHTTP_CONFIGURATION );
    
    ns_server->addEntryPoint( Wt::Application, boost::bind( &createAppForServer, _1, create_application ) );
    
    if( ns_server->start() )
    {
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
  }
  
  
  void killServer()
  {
    std::lock_guard<std::mutex> serverlock( ns_servermutex );
    
    if( ns_server )
    {
      std::cerr << "About to stop server" << std::endl;
      ns_server->stop();
      delete ns_server;
      ns_server = 0;
      std::cerr << "Stopped and killed server" << std::endl;
    }
  }//void experimental_killServer()
  
  
  void add_allowed_session_token( const char *session_id )
  {
#warning "Need to make add_session_id() handle more than just single session ID"
    lock_guard<mutex> lock(ns_externalid_mutex);
    ns_externalid = session_id;
  }//void interspec_add_allowed_session_token( const char *session_id )
  
  
  int remove_allowed_session_token( const char *session_token )
  {
#warning "Need to make interspec_remove_allowed_session_token() take into account if sesion with token had been seen"
    lock_guard<mutex> lock(ns_externalid_mutex);
    if( session_token == ns_externalid )
    {
      ns_externalid = "";
      return 0;
    }
    
    return -1;
  }//int interspec_remove_allowed_session_token( const char *session_token )
  
  std::string external_id()
  {
    lock_guard<mutex> lock(ns_externalid_mutex);
    return ns_externalid;
  }

  int open_file_in_session( const char *sessionToken, const char *files_json )
  {
#warning "Need to actually test interspec_open_file"
    
    //Are we garunteed to recieve the entire message at once?
    cerr << "Opening files not tested!" << endl;
    
    vector<string> files;
    
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
        if( UtilityFunctions::is_file(valstr) )
          files.push_back(valstr);
        else
          cerr << "File '" << valstr << "' is not a file" << endl;
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
      
      for( auto filename : files )
      {
        if( app->userOpenFromFileSystem( filename ) )
          numopened += 1;
        else
          cerr << "InterSpec failed to open file filename" << endl;
      }
      
      app->triggerUpdate();
      
      return numopened;
    }else
    {
      //I dont know why we would get here... but lets deal with it JIC
      cerr << "There is no app with externalid=" << externalid << endl;
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
}//namespace InterSpecServer

