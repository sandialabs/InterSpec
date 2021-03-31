/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#if __APPLE__
#include "TargetConditionals.h"
#endif
#if( defined(_WIN64) || defined(_WIN32) )
#define WIN32_LEAN_AND_MEAN 1
#include <Windows.h>
#include <Shellapi.h>
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



#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "DesktopWithNativeBrowser/DesktopWithNativeBrowser.h"


using namespace std;

namespace
{
  int PortServedOn = -1;
  std::string UrlServedOn = "";
  boost::mutex ServedOnMutex;
    
  //And now for our experimental style server - e.g. we dont put the WServer
  //  in its own thread.  I dont really know if this is thread or otherwise safe
  
  
  Wt::WServer *ns_server = 0;
  boost::mutex ns_servermutex;
  
//  std::vector<const InterSpecApp *> AppInstances;
//  boost::mutex AppInstancesMutex;
  
//  void appDestroyed( const InterSpecApp *app )
//  {
//    boost::lock_guard<boost::mutex> lock( AppInstancesMutex );
//    std::vector<const InterSpecApp *>::iterator pos 
//                       = find( AppInstances.begin(), AppInstances.end(), app );
//    if( pos != AppInstances.end() )
//     AppInstances.erase( pos );
//  }//void appDestroyed(...)
  
  Wt::WApplication *createAppForServer( const Wt::WEnvironment &env, 
  	                            Wt::WApplication::ApplicationCreator appCreator )
  {
    Wt::WApplication *app = appCreator( env );
//    InterSpecApp *viewerApp = dynamic_cast<InterSpecApp *>( app );
    
//    if( viewerApp )
 //   {
//      {
//        boost::lock_guard<boost::mutex> lock( AppInstancesMutex );
//        AppInstances.push_back( viewerApp );
//      }
      
//      {
//      	Wt::WApplication::UpdateLock lock( app );
//        viewerApp->destructing().connect( boost::bind( &appDestroyed, _1 ) );
//      }
//    }//if( viewerApp )
    
    return app;
  }//WApplication *createAppForServer(...)
}//namespace


//As an experiment to try to get spectra files from
#define CREATE_REACHBACK_SERVER 0

#if( CREATE_REACHBACK_SERVER )
#include <Wt/WResource>
#include <Wt/Http/Client>
#include <Wt/Http/Response>

#include <boost/foreach.hpp>
#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH


namespace
{
  class DetectorUploadResource : public Wt::WResource
  {
  public:
    DetectorUploadResource( Wt::WObject *parent = 0 )
    : WResource( parent )
    {
    }
    
    ~DetectorUploadResource()
    {
      beingDeleted();
    }
    
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response)
    {
      cerr << "Recieved request from " << request.clientAddress()
      << " with length " << request.contentLength() << endl;
      
      namespace pt = boost::posix_time;
      const pt::ptime server_local_time = pt::second_clock::local_time();
      const pt::ptime utc_time = pt::second_clock::universal_time();
      const string dt_str = pt::to_iso_string( server_local_time );
      
      stringstream error_msg;
      bool valid_request = true;
      const Wt::Http::ParameterMap &parameters = request.getParameterMap();
      
      //      if( parameters.find( "serialnumber" ) == parameters.end() )
      //      {
      //        valid_request = false;
      //        error_msg << "Request did not have a serial number parameter." << endl;
      //      }//if( parameters.find("serialnumber") == parameters.end() )
      
      const Wt::Http::UploadedFileMap &uploadedFiled = request.uploadedFiles();
      const size_t nfiles = uploadedFiled.size();
      if( nfiles != 1 )
      {
        valid_request = false;
        error_msg << "Request had " << nfiles << " files; not the expected 1." << endl;
      }//if( did not have an uploaded file )
      
      if( request.tooLarge() )
      {
        valid_request = false;
        error_msg << "Request was was rejected due to request being to large "
        << request.tooLarge() << endl;
      }//if( request.tooLarge() )
      
      if( !boost::algorithm::iequals( "POST", request.method() ) )
      {
        valid_request = false;
        error_msg << "Request was of type " << request.method() << "; not POST." << endl;
      }//if( not a POST command )
      
      if( !valid_request )
      {
        cerr << "Request Rejected" << endl
        << "\tServer Local Time: " << dt_str << endl
        << "\tIP Address: " << request.clientAddress() << endl
        << "\tContent Length: " << request.contentLength() << endl
        << "\tMethod: " << request.method() << endl
        << "\treason: " << error_msg.str()
        << endl << endl;
        
        response.addHeader( "X-ReachbackError", "Delivery Failed" );
        //    response.out() << error_msg.str();
        response.setStatus( 500 ); //generic internal server error
        return;
      }//if( !valid_request )
      
      //If we made it here, we have a valid request
      if( parameters.find( "serialnumber" ) == parameters.end() )
      {
        Wt::Http::ParameterValues serials = parameters.find("serialnumber")->second;
        if( serials.size() )
          cout << "Detector serial number claimed to be " << serials[0] << endl;
      }
      
      foreach( const Wt::Http::UploadedFileMap::value_type &a, uploadedFiled )
      {
        const string clientName = a.second.clientFileName();
        const string spoolName = a.second.spoolFileName();
        cout << "Recieved file '" << clientName
             << "' at filesystem location '" << spoolName << "'" << endl;
      }
      
      response.addHeader( "X-ReachbackSuccess", "Delivery Successful" );
      
      response.out() << "200 OK\r\n";
      response.setStatus( 200 );
    }//handleRequest(...)
  };//class DetectorUploadResource
  
  //Protected by ns_servermutex
  DetectorUploadResource *ns_reachback_server = 0;
}//namespace
#endif  //CREATE_REACHBACK_SERVER


boost::mutex g_instance_mutex;
DesktopWithNativeBrowser::InstanceState g_other_instance_state 
                              = DesktopWithNativeBrowser::InstanceStateNotKnown;
std::string g_other_instance_url = "";




  
namespace DesktopWithNativeBrowser
{  
  int portBeingServedOn()
  {
    boost::lock_guard<boost::mutex> lock( ServedOnMutex );
    return PortServedOn;
  }
  
  std::string urlBeingServedOn()
  {
    boost::lock_guard<boost::mutex> lock( ServedOnMutex );
    return UrlServedOn;
  }

  
//  std::vector<const InterSpecApp *> runningInstances()
//  {
//    std::vector<const InterSpecApp *> instances;
    
//    {
//      boost::lock_guard<boost::mutex> lock( AppInstancesMutex );
//      instances = AppInstances;
//    }
    
//    return instances;
//  }//std::vector<const InterSpecApp *> runningInstances()
  
	
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
    
#if( BUILD_AS_QT4_APP || BUILD_AS_QT5_APP )
    if( boost::filesystem::exists( "data/config/wt_config_qt.xml" ) )
      return (boost::filesystem::current_path() / boost::filesystem::path("data/config/wt_config_qt.xml")).string<std::string>();
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
  
  
  
  bool shouldForceServe( int argc, char *argv[] )
  {
    for( int i = 1; i < argc; ++i )
    {
      if( strcmp(argv[i],"--forceserve") == 0 )
        return true;
    }
    
    return false;
  }//bool shouldForceServe( int argc, char *argv[] )
  
  bool shouldOpenBlankTab( int argc, char *argv[] )
  {
    for( int i = 1; i < argc; ++i )
    {
      if( strcmp(argv[i],"--nobrowsertab") == 0 )
        return false;
    }
    
    return true;
  }//bool shouldOpenBlankTab( int argc, char *argv[] )
  
  boost::filesystem::path fileToOpen( int argc, char *argv[] )
  {
    try
    {
      if( argc < 2 )
        return "";
      
      std::string filename = argv[1];
      if( boost::algorithm::starts_with(filename, "'" )
          || boost::algorithm::starts_with(filename, "\"" ) )
        filename = filename.substr(1);
      if( boost::algorithm::ends_with(filename, "'" )
         || boost::algorithm::ends_with(filename, "\"" ) )
        filename = filename.substr(0,filename.size()-1);
                                           
      if( boost::filesystem::exists( filename ) )
        return filename;
      
      std::cerr << "The file '" << filename << "' doesnt exist" << std::endl;
      return "";
    }catch( std::exception &e )
    {
      std::cerr << "fileToOpen(...) caught: " << e.what() << std::endl;
    }//try / catch
    
    return "";
  }//path fileToOpen(...)
  
  
  boost::filesystem::path instanceInfoFile( const char *exename )
  {
    return boost::filesystem::temp_directory_path()
          / (boost::filesystem::path(exename).filename().string() + "info");
  }
  
  
#if( !defined(IOS) && !defined(ANDROID) )
  void openUrl( std::string url, boost::filesystem::path file )
  {
    if( !file.empty() )
    {  
      DbToFilesystemLink::FileIdToLocation filerequest;
      filerequest.m_foregroundFilePath = file.string<std::string>();
      const int index = DbToFilesystemLink::addFileToOpenToDatabase( filerequest );
      if( index >= 0 )
        url += "?specfile=" + boost::lexical_cast<std::string>( index );
    }//if( !file.empty() )
    
#if( defined(_WIN64) || defined(_WIN32) )
    //MS Windows product here
    string command;
    
    if( boost::filesystem::exists( "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe" ) )
      command = "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe";
   else if( boost::filesystem::exists( "C:/Program Files/Google/Chrome/Application/chrome.exe" ) )
      command = "C:/Program Files/Google/Chrome/Application/chrome.exe";
   else if( boost::filesystem::exists( "C:/Program Files (x86)/Mozilla Firefox/firefox.exe" ) )
      command = "C:/Program Files (x86)/Mozilla Firefox/firefox.exe";
    else if( boost::filesystem::exists( "C:/Program Files/Mozilla Firefox/firefox.exe" ) )
      command = "C:/Program Files/Mozilla Firefox/firefox.exe";

    if( command.empty() )
      ShellExecute( NULL, "open", url.c_str(), NULL, NULL, SW_SHOWNORMAL );
    else
      ShellExecute( NULL, "open", command.c_str(), url.c_str(), NULL, SW_SHOWNORMAL );
      
#elif __APPLE__ && __MACH__
    //Apple product here
//#if TARGET_OS_IPHONE || TARGET_IPHONE_SIMULATOR
//    BOOST_STATIC_ASSERT_MSG(false, "DesktopRunInDefaultBrowser not supported on iOS");
//#elif TARGET_OS_MAC
    system( ("open "+url).c_str() );
//#else
//    BOOST_STATIC_ASSERT_MSG(false, "Unsupported apple operating system");
//#endif
#elif __ANDROID__
    //Android
    BOOST_STATIC_ASSERT_MSG(false, "DesktopRunInDefaultBrowser not supported on Android");
#elif __linux__ || __linux || __unix
    //Linux
    system("x-www-browser " + url);
    //        system("firefox -new-tab " + url);
#endif
  }//void openUrl(...)
#endif //#if( !defined(IOS) && !defined(ANDROID) )  



  void experimental_startServer( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator createApplication )
  {
    changeToBaseDir( argc, argv );
    const string xml_config_path = getWtConfigXml( argc, argv );
    
    boost::lock_guard<boost::mutex> serverlock( ns_servermutex );
    if( ns_server )
    {
      std::cerr << "experimental_startServer: already running" << std::endl;
      return;
    }

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
    
    
#if( CREATE_REACHBACK_SERVER )
    ns_reachback_server = new DetectorUploadResource();
    ns_server->addResource( ns_reachback_server, "/upload" );
    cout << "Added a DetectorUploadResource to /upload" << endl;
#endif
  
    
      if( ns_server->start() )
      {
        const int port = ns_server->httpPort();
        std::string this_url = "http://127.0.0.1:" + boost::lexical_cast<string>(port);

        {
          boost::lock_guard<boost::mutex> lock( ServedOnMutex );
          PortServedOn = port;
          UrlServedOn = this_url;
          
#if( CREATE_REACHBACK_SERVER )
          cout << "UrlServedOn=" << UrlServedOn << ", port is " << PortServedOn << endl;
#endif
        }
      }//if( server.start() )
  }//void experimental_startServer()
  
  
  void experimental_killServer()
  {
    boost::lock_guard<boost::mutex> serverlock( ns_servermutex );
    
    if( ns_server )
    {
      std::cerr << "About to stop server" << std::endl;
      ns_server->stop();
      delete ns_server;
      ns_server = 0;
      std::cerr << "Stopped and killed server" << std::endl;
      
#if( CREATE_REACHBACK_SERVER )
      if( ns_reachback_server )
        delete ns_reachback_server;
      ns_reachback_server = 0;
      std::cerr << "Deleted rechback server" << std::endl;
#endif
    }
  }//void experimental_killServer()
  
  
  int runInDefaultBrowser( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator createApplication )
  {
    using namespace std;
    using namespace Wt;
    int return_code = EXIT_SUCCESS;
    
    string this_url = "";
    
    try
    {      
      const string xml_config_path = getWtConfigXml( argc, argv );
      
      cout << "Using wt_config file: " << xml_config_path << endl;
      
      WServer server( argv[0], xml_config_path );
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
      
      server.setServerConfiguration( argc_wthttp, argv_wthttp, WTHTTP_CONFIGURATION );
      
      try
      {
        server.addEntryPoint( Wt::Application, boost::bind( &createAppForServer, _1, createApplication ) );
        
        ServerAliveResource serverAliveRes;
        server.addResource( &serverAliveRes, "/checkalive" );
        
        if( server.start() )
        {
          const int port = server.httpPort();
          this_url = "http://127.0.0.1:" + boost::lexical_cast<string>(port);

          {
            boost::lock_guard<boost::mutex> lock( ServedOnMutex );
            PortServedOn = port;
            UrlServedOn = this_url;
  	  }
          
          boost::filesystem::path exe_info_file = instanceInfoFile( argv[0] );
          
          {//Begin code block to add this_url to info file
            //Should get a mutex on this file here!
            std::ofstream info_file( exe_info_file.c_str() );
            if( info_file.is_open() )
              info_file << this_url;
            else
              std::cerr << "Failed to open instance info file: "
                        << exe_info_file << std::endl;
          }//End code block to add this_url to info file
          
#if( !defined(IOS) && !defined(ANDROID) )
          boost::filesystem::path specfile = fileToOpen( argc, argv );
          
          if( !specfile.empty() || shouldOpenBlankTab(argc, argv) ) 
            openUrl( this_url, specfile );
#endif //#if( !defined(IOS) && !defined(ANDROID) )

//Should add a custom signal handler here, as well as a cusomt loop to
//  waite untill some global value changes to false
//Or could try to lock a mutex, and then lock it again so that someone has
//  to unlock it inbetween the two instance of locking, or maybye

          int sig = WServer::waitForShutdown(argv[0]);
          
          server.log("notice") << "Shutdown (signal = " << sig << ")";
          server.stop();
        }//if( server.start() )
      }catch (std::exception& e)
      {
        server.log("fatal") << e.what();
        return_code = EXIT_FAILURE;
      }
    }catch( WServer::Exception &e )
    {
      cerr << "WServer exception: " << e.what() << endl;
      return_code = EXIT_FAILURE;
    }catch( std::exception &e )
    {
      cerr << "Fatal exception: " << e.what() << endl;
      return_code = EXIT_FAILURE;
    }
    
    {
      boost::lock_guard<boost::mutex> lock( ServedOnMutex );
      PortServedOn = -1;
      UrlServedOn = "";
    }

    
    try
    {
      //XXX - in principle should be carefull and only remove the file if this
      //      url is the only one in it
      if( !this_url.empty() )
        boost::filesystem::remove( instanceInfoFile(argv[0]) );
    }catch( std::exception &e )
    {
      std::cerr << "Caught error removing instance info file for url "
                << this_url << std::endl;
    }//try / catch to delete the info file
    
    return return_code;
  }//int runInDefaultBrowser(...)
  
  
  void handleHttpCallback( boost::system::error_code err,
                           const Wt::Http::Message &response,
                           Wt::Http::Client *client,
                           int argc, char *argv[],
                           std::string url,
                           Wt::WApplication::ApplicationCreator appCreator )
  {
    boost::lock_guard<boost::mutex> lock( g_instance_mutex );
    
    if( client )
      delete client;
    
    if( !err && (response.status() == 200) )
    {
      std::cerr << "handleHttpCallback(...) successfuly got instance!"
                << "Response body: " << response.body() << std::endl;
      g_other_instance_state = OtherInstanceRunning;
      g_other_instance_url = url;
      return;
    }//if( we were successful )
    
    g_other_instance_state = NoOtherInstanceRunning;
  }//void handleHttpCallback(...)
                          
  
  bool checkForExistingInstanceAndCallIt( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator appCreator  )
  {
    //See also: http://rosettacode.org/wiki/Determine_if_only_one_instance_is_running
    using namespace Wt;
    using namespace std;
    
    const boost::filesystem::path exe_info_file = instanceInfoFile( argv[0] );
    
    if( boost::filesystem::exists( exe_info_file ) )
    {
      //look for the last non-empty line in this file, then we'll check if this
      //  is a working app
      string url, line;
      ifstream address( exe_info_file.c_str() );
      
      while( getline( address, line ) )
        if( boost::algorithm::starts_with(line, "http://" ) )
          url = line;
      
      if( url.empty() )
        return false;
      
      if( !url.empty() )
      {
        Http::Client *client = new Http::Client();
        client->setTimeout( 1 );
        client->setMaximumResponseSize( 1024 );
        
        Http::Message msg;
//        msg.addHeader( "InterSpecHeader", "InterSpecMessage" );
        client->done().connect( boost::bind( &handleHttpCallback, _1, _2,
                                      client, argc, argv, url, appCreator ) );
        const bool posted = client->post( url + "/checkalive", msg );
        if( posted )
          return true;
        
        cerr << "Failed to post check to " << url << endl;
        delete client;
        return false;
      }//if( !url.empty() )
    }//if( boost::filesystem::exists( exe_info_file ) )
    
    return false;
  }//void checkForExistingInstanceAndCallIt(...)
  
  int DesktopRunInDefaultBrowser( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator createApplication )
  {
    using namespace std;
    using namespace Wt;
    bool startNewSession = true;
    
    try
    {
      //Our current working directory may not be where we want, lets change this
      //  --currently only tested on Windows (probably the same on Linux, I 
      //    still have yet to decide on OSX
#if( defined(_WIN64) || defined(_WIN32) )
      boost::filesystem::path exe = argv[0];
      boost::filesystem::path exe_dir = exe.parent_path();
      if( !exe_dir.empty() )
        boost::filesystem::current_path( exe_dir );
#endif //#if( defined(_WIN64) || defined(_WIN32) )

      //We have to create a temporary server, so checkForExistingInstanceAndCallIt(...)
      //  can communicate and make sure the other possible servers are alive.
      const string xml_config_path = getWtConfigXml( argc, argv );
      WServer tempServer( argv[0], xml_config_path );
      
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
      
      changeToBaseDir( argc, argv );
      const bool forceRun = shouldForceServe( argc, argv );
      
      if( !forceRun )
      {
        tempServer.setServerConfiguration( argc_wthttp, argv_wthttp, WTHTTP_CONFIGURATION );
      
        if( tempServer.start() )
        {
          const bool isPotentiallyRunning
               = checkForExistingInstanceAndCallIt( argc, argv, createApplication );
    
          if( isPotentiallyRunning )
          {
            int running = -1;
            int numtries = 0;
            while( running<0 && numtries < 110 )  //110 sections of at least 10 ms, is 1.1 seconds
            {
              boost::lock_guard<boost::mutex> lock( g_instance_mutex );
        
              switch( g_other_instance_state )
              {
                case InstanceStateNotKnown:
                  break;
                case OtherInstanceRunning:
                  running = 1;
                  break;
                case NoOtherInstanceRunning:
                  running = 0;
                  break;
              }//switch( g_other_instance_state )
        
              ++numtries;
              boost::this_thread::sleep( boost::posix_time::milliseconds(10) );
            }//while( numtries < 11 )
      
            if( running == 1 )
            {
#if( !defined(IOS) && !defined(ANDROID) )
              boost::filesystem::path specfile = fileToOpen( argc, argv );
              std::cout << " Opening " << g_other_instance_url << " for file " 
                        << specfile << std::endl;
              if( !specfile.empty() || shouldOpenBlankTab(argc, argv) )
                openUrl( g_other_instance_url, specfile );
#endif //#if( !defined(IOS) && !defined(ANDROID) )
              startNewSession = false;
            }//if( running == 1 )
          }//if( isPotentiallyRunning )
        }//if( tempServer.start() )
      }//if( forceRun ) / else
    }catch( std::exception &e )
    {
      cerr << "Caught exception makeing temporary server in"
           << " DesktopRunInDefaultBrowser: " << e.what() << endl;
    }//try / catch
    
    if( startNewSession )
      return runInDefaultBrowser( argc, argv, createApplication );
    
    return EXIT_SUCCESS;
  }//int DesktopRunInDefaultBrowser(...)
}//namespace
