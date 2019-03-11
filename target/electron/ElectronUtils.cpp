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

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <csignal>
//#include <unistd.h>
#include <sys/stat.h>

#include <mutex>
#include <atomic>
#include <string>
#include <thread>
#include <chrono>
#include <memory>
#include <iostream>
#include <condition_variable>

#include <boost/asio.hpp>
#include <boost/asio/signal_set.hpp>

#include <Wt/WServer>
#include <Wt/Json/Array>
#include <Wt/Json/Object>
#include <Wt/Json/Parser>

#include "websocketpp/client.hpp"
#include "websocketpp/config/asio_no_tls_client.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/ResourceUpdate.h"
#include "InterSpec/InterSpecServer.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "InterSpec/MassAttenuationTool.h"
#include "target/electron/ElectronUtils.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

using namespace std;
typedef websocketpp::client<websocketpp::config::asio_client> ws_client_t;


namespace
{
  //The externalid of the primary session (the main electron window)
  string ns_externalid;
  std::mutex ns_externalid_mutex;
  
  std::mutex ns_run_mutex;
  std::condition_variable ns_run_cv;
  bool ns_keep_running = true;
  
  
  Wt::WApplication *createApplication( const Wt::WEnvironment &env )
  {
    return new InterSpecApp( env );
  }// Wt::WApplication *createApplication(const Wt::WEnvironment& env)
  
  
  //Returns true if should keep running
  bool handle_parent_message( const std::string &recieved )
  {
    const bool shutdown = (recieved.find( "command=exit" ) != string::npos);
    if( shutdown )
      return false;
    
    if( UtilityFunctions::istarts_with(recieved, "openfile=") )
    {
      //Are we garunteed to recieve the entire message at once?
      cerr << "Opening files not tested!" << endl;
      
      vector<string> files;
      
      try
      {
        Wt::Json::Value result;
        Wt::Json::parse( recieved.substr(9), result, true );
        
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
        cerr << "Failed to parse '" << recieved.substr(9)
        << "' as valid JSON. Issue: " << e.what() << endl;
      }//try / catch
      
      
      //Look for session with 'externalid' and open file...
      string externalid;
      {
        lock_guard<mutex> lock(ns_externalid_mutex);
        externalid = ns_externalid;
      }
      
      InterSpecApp *app = InterSpecApp::instanceFromExtenalIdString( externalid );
      if( app )
      {
        Wt::WApplication::UpdateLock applock( app );
        
        for( auto filename : files )
        {
          if( !app->userOpenFromFileSystem( filename ) )
            cerr << "InterSpec failed to open file filename" << endl;
        }
        
        app->triggerUpdate();
      }else
      {
        cerr << "There is no app with externalid=" << externalid
        << "; not opening files" << endl;
      }//if( app ) / else
    }else
    {
      cerr << "Unrecognized command: " << recieved << endl;
    }

    return true;
  }//void handle_parent_message( const std::string &msg )
  
  void do_ipc_shutdown()
  {
    cout << "Shutting dow IPC mechanism" << endl;
    std::unique_lock<std::mutex> run_lock( ns_run_mutex );
    ns_keep_running = false;
    run_lock.unlock();
    ns_run_cv.notify_all();
    cout << "Done shutting down IPC" << endl;
  }
  
}//namespace

namespace ElectronUtils
{

int run_app( int argc, char *argv[] )
{
  //Port used to communicate with Electron App.
  int ipc_port = 0;
  string userdatadir;
  
  for( int i = 1; i < argc-1; ++i )
  {
    if( argv[i] == string("--ipc") )
    {
      try
      {
        ipc_port = std::stoi( argv[i+1] );
      }catch( std::exception & )
      {
        std::cerr << "IntializeError: --ipc value of '" << argv[i+1]
                  << "' is not a valid integer." << endl;
        return EXIT_FAILURE;
      }
    }else if( argv[i] == string("--userdatadir") )
    {
      userdatadir = argv[i+1];
    }else if( argv[i] == string("--externalid") )
    {
      lock_guard<mutex> lock(ns_externalid_mutex);
      ns_externalid = argv[i+1];
	}else if (argv[i] == string("--basedir"))
	{
    InterSpec::setStaticDataDirectory( argv[i + 1] );

		try
		{
		  boost::filesystem::current_path(argv[i + 1]);
		  cout << "Changed to cwd: " << argv[i + 1] << endl;
		}catch ( std::exception &e )
		{
		  cerr << "Unable to change to directory: '" << argv[i + 1] << "' :" << e.what() << endl;
		}
	}
  }//for( int i = 1; i < argc-1; ++i )
  
  if( !ipc_port )
  {
    std::cerr << "IntializeError: no \"--ipc\" specified." << endl;
    return EXIT_FAILURE;
  }
  
  if( userdatadir.empty() )
  {
    std::cerr << "IntializeError: no \"--userdatadir\" specified." << endl;
    return EXIT_FAILURE;
  }
  
  if( !UtilityFunctions::create_directory(userdatadir) )
  {
    cerr << "IntializeError: Failed to create directory '"
         << userdatadir << "' for user data." << endl;
    return EXIT_FAILURE;
  }
  
  {//begin lock on ns_externalid_mutex
    lock_guard<mutex> lock(ns_externalid_mutex);
    if( ns_externalid.empty() )
    {
      std::cerr << "IntializeError: no \"--externalid\" specified." << endl;
      return EXIT_FAILURE;
    }
    cout << "Using externalid: " << ns_externalid << endl;
  }//end lock on ns_externalid_mutex
  
  
  const string preffile = UtilityFunctions::append_path( userdatadir, "InterSpecUserData.db" );
  
  cout << "Will set user preferences file to: '" << preffile << "'" << endl;
  DataBaseUtils::setPreferenceDatabaseFile( preffile );
  DbToFilesystemLink::setFileNumToFilePathDBNameBasePath( userdatadir );
  ResourceUpdate::setUserDataDirectory( userdatadir );
  InterSpec::setWritableDataDirectory( userdatadir );


  cout << "Will make sure preferences database is up to date" << endl;
  DataBaseVersionUpgrade::checkAndUpgradeVersion();
  
  ResourceUpdate::setupGlobalPrefsFromDb();
  
  cout << "Using port " << ipc_port << " for interprocess communication" << endl;
  
  try
  {
    InterSpecServer::startServer( argc, argv, &createApplication );
  }catch( std::runtime_error &e )
  {
    std::cerr << "IntializeError: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  
  //const int interspec_port = InterSpecServer::portBeingServedOn();
  std::string const interspec_url = InterSpecServer::urlBeingServedOn();
  
  
  ws_client_t ws_client;
  websocketpp::connection_hdl ws_handle;
  websocketpp::lib::asio::io_service ws_io_service;
  
  ws_client.clear_access_channels(websocketpp::log::alevel::all);
  ws_client.set_access_channels(websocketpp::log::alevel::connect);
  ws_client.set_access_channels(websocketpp::log::alevel::disconnect);
  ws_client.set_access_channels(websocketpp::log::alevel::app);
  
  websocketpp::log::level elog_level =( 0x0
                                       | websocketpp::log::elevel::devel
                                       | websocketpp::log::elevel::library
                                       | websocketpp::log::elevel::info
                                       | websocketpp::log::elevel::warn
                                       | websocketpp::log::elevel::rerror
                                       | websocketpp::log::elevel::fatal
                                       );
  ws_client.get_elog().set_channels( elog_level );
  
  ws_client.init_asio( &ws_io_service );

#if( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  boost::asio::signal_set signals(ws_io_service, SIGINT, SIGTERM, SIGABRT);
#else
  boost::asio::signal_set signals( ws_io_service, SIGINT, SIGTERM, SIGHUP );
#endif

  auto handle_interupt = []( const boost::system::error_code& error, int signal_number ) {
    if( error )  // m_signals is destructing or canceled or something maybe?
      return;
    do_ipc_shutdown();
  };//handle_interupt(...)
  
  signals.async_wait( [&handle_interupt](const boost::system::error_code& error, int signal_number){
    cerr << "boost::asio::signal_set, recived error: " << error.message() << " (" << signal_number << ")" << endl;
    handle_interupt(error,signal_number);
  });
  
  ws_client.set_interrupt_handler( [&handle_interupt]( websocketpp::connection_hdl hdl ){
    cerr << "WebSocket recieved interupt" << endl;
    handle_interupt( boost::system::error_code{}, 0 );
  });
  
  //The open handler is called after the WebSocket handshake is complete and
  //  the connection is considered OPEN.
  ws_client.set_open_handler( [&ws_client,&ws_handle,&interspec_url](websocketpp::connection_hdl hdl){
    cout << "WebSocket connection opened in c++" << endl;
    ws_handle = hdl;
    websocketpp::lib::error_code ec;
    ws_client.send( hdl, "InterSpecUrl=" + interspec_url, websocketpp::frame::opcode::TEXT, ec );
    if( ec )
    {
      cerr << "IntializeError: Failed to send InterSpecUrl after open" << endl;
      do_ipc_shutdown();
    }
  } );
  
  //The close handler is called immediately after the connection is closed.
  ws_client.set_close_handler( [&ws_handle,&ws_client](websocketpp::connection_hdl hdl){
    cout << "WebSocket close connection handler" << endl;
    ws_handle.reset();
    do_ipc_shutdown();
    ws_client.stop_perpetual();
    ws_client.stop();
  } );
  
  //The fail handler is called whenever the connection fails while the
  //  handshake is bring processed.
  ws_client.set_fail_handler( [&ws_handle,&ws_client](websocketpp::connection_hdl hdl){
    cerr << "IntializeError: WebSocket handshake failed" << endl;
    ws_handle.reset();
    do_ipc_shutdown();
    ws_client.stop_perpetual();
    ws_client.stop();
  } );

  ws_client.set_message_handler( []( websocketpp::connection_hdl hdl, ws_client_t::message_ptr msg ){
    if( !handle_parent_message( msg->get_payload() ) )
    {
      cerr << "Recieved shutdown message, will close connection." << endl;
      do_ipc_shutdown();
    }
  } );
  
  
  std::thread runner( [&ws_io_service](){ ws_io_service.run(); } );
  
  
  websocketpp::lib::error_code ec;
  
  {//begin scope of con
#ifdef __linux__
    ws_client_t::connection_ptr con = ws_client.get_connection("http://127.0.0.1:" + to_string(ipc_port), ec);
#else
    ws_client_t::connection_ptr con = ws_client.get_connection("http://[::1]:" + to_string(ipc_port), ec);
#endif
  
    cout << "Got con" << endl;
    if( ec )
    {
      cerr << "Get Connection Error: " << ec.message() << endl;
    
      return EXIT_FAILURE;
    }
  
    ws_handle = con->get_handle();
    cout << "Got handle" << endl;
  
    ws_client.connect( con );
    cout << "Queded connection" << endl;
  }//end scope of con
  
  
  {
    std::unique_lock<std::mutex> run_lock( ns_run_mutex );
    while( ns_keep_running )
      ns_run_cv.wait( run_lock, []{return !ns_keep_running;} );
  }
  
  cout << "Will kill InterSpec server" << endl;
  InterSpecServer::killServer();
  
  ws_client.send( ws_handle, "ServerKilled", websocketpp::frame::opcode::TEXT, ec );
  if( ec )
  {
    cerr << "ShutdownError: Failed to send 'ServerKilled' message: " << ec.message() << endl;
    do_ipc_shutdown();
  }
  
  cout << "Will stop ws_client" << endl;
  if( ws_client.is_listening() )
    ws_client.stop_listening( ec );
  
  if( ws_client.get_con_from_hdl(ws_handle,ec) )
    ws_client.close( ws_handle, websocketpp::close::status::going_away, "", ec );
  
  
  cout << "Will join IPC main thread" << endl;
  runner.join();
  
  cout << "Exiting native sub-process code" << endl;
  return EXIT_SUCCESS;
}//int run_app( int argc, char *argv[] )

  std::string external_id()
  {
    lock_guard<mutex> lock(ns_externalid_mutex);
    return ns_externalid;
  }
}//namespace ElectronUtils
