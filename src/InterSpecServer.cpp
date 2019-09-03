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
#include <string>
#include <cstdlib>
#include <stdlib.h>

#if __APPLE__
#include "TargetConditionals.h"
#endif
#include <boost/thread/thread.hpp>
#include <boost/static_assert.hpp>
#include <boost/algorithm/string.hpp>
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
#include <stdlib.h>
#include <iostream>
#include <fstream>


#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>


#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DbToFilesystemLink.h"


using namespace std;

namespace
{
  Wt::WApplication *create_application( const Wt::WEnvironment &env )
  {
    return new InterSpecApp( env );
  }// Wt::WApplication *createApplication(const Wt::WEnvironment& env)
}

namespace InterSpecServer
{
  int sm_portServedOn = -1;
  std::string sm_urlServedOn = "";
  boost::mutex sm_servedOnMutex;
  
  Wt::WServer *ns_server = 0;
  boost::mutex ns_servermutex;
  
  
  int portBeingServedOn()
  {
    boost::lock_guard<boost::mutex> lock( sm_servedOnMutex );
    return sm_portServedOn;
  }
  
  std::string urlBeingServedOn()
  {
    boost::lock_guard<boost::mutex> lock( sm_servedOnMutex );
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
      if( boost::algorithm::starts_with( arg, "--basedir=" ) )
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
    
    boost::lock_guard<boost::mutex> serverlock( ns_servermutex );
    if( ns_server )
    {
      std::cerr << "experimental_startServer: already running" << std::endl;
      return;
    }
    

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
        boost::lock_guard<boost::mutex> lock( sm_servedOnMutex );
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
    boost::lock_guard<boost::mutex> serverlock( ns_servermutex );
    if( ns_server )
    {
      std::cerr << "experimental_startServer: already running" << std::endl;
      return;
    }
    
    if( name.empty() )
      name = "\0";
    if( basedir.empty() )
      basedir = ".";
    
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
        boost::lock_guard<boost::mutex> lock( sm_servedOnMutex );
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
    boost::lock_guard<boost::mutex> serverlock( ns_servermutex );
    
    if( ns_server )
    {
      std::cerr << "About to stop server" << std::endl;
      ns_server->stop();
      delete ns_server;
      ns_server = 0;
      std::cerr << "Stopped and killed server" << std::endl;
    }
  }//void experimental_killServer()
  
  
  
}//namespace InterSpecServer

