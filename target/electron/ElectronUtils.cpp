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
#include <boost/filesystem.hpp>
#include <boost/asio/signal_set.hpp>

#include <Wt/WRandom>
#include <Wt/WServer>
#include <Wt/Json/Array>
#include <Wt/Json/Object>
#include <Wt/Json/Parser>

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/ResourceUpdate.h"
#include "InterSpec/InterSpecServer.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "InterSpec/MassAttenuationTool.h"
#include "target/electron/ElectronUtils.h"
#include "SpecUtils/SerialToDetectorModel.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

using namespace std;


namespace
{
  //The externalid of the primary session (the main electron window)
  string ns_externalid;
  std::mutex ns_externalid_mutex;
  

}//namespace

namespace ElectronUtils
{

std::string external_id()
{
  lock_guard<mutex> lock(ns_externalid_mutex);
  return ns_externalid;
}
  
bool requestNewCleanSession()
{
  const string newexternalid = Wt::WRandom::generateId( 16 );
  
#warning "Need to implement calling back into JS for requestNewCleanSession"
  /*
  std::cout << "requestNewCleanSession()"<< endl;
  
  {
    lock_guard<mutex> lock(ns_externalid_mutex);
    ns_externalid = newexternalid;
  }
  
  std::unique_lock<std::mutex> run_lock( ns_send_message_mutex );
  
  if( !ns_send_message )
  {
    std::cerr << "requestNewCleanSession(): No active connection" << std::endl;
    return false;
  }
  
  try
  {
    ns_send_message( "NewCleanSession:" + newexternalid );
  }catch( std::exception &e )
  {
    cerr << "requestNewCleanSession(): caught exception sendign WS message: " << e.what() << endl;
    return false;
  }
  
  std::cout << "requestNewCleanSession() successfully sent" << endl;
*/
  return true;
}//void requestNewCleanSession()
  
  
bool notifyNodeJsOfNewSessionLoad( const std::string sessionid )
{
#warning "Need to implement calling back into JS for notifyNodeJsOfNewSessionLoad"
/*
 std::unique_lock<std::mutex> run_lock( ns_send_message_mutex );
 
  if( !ns_send_message )
  {
    std::cerr << "notifyNodeJsOfNewSessionLoad(): No active connection" << std::endl;
    return false;
  }
  
  try
  {
    ns_send_message( "SessionFinishedLoading:" + sessionid );
  }catch( std::exception &e )
  {
    cerr << "notifyNodeJsOfNewSessionLoad(): caught exception sendign WS message: " << e.what() << endl;
    return false;
  }
  
  std::cout << "notifyNodeJsOfNewSessionLoad() successfully sent" << endl;
*/
  
  return true;
}//bool notifyNodeJsOfNewSessionLoad( const std::string sessionid )
}//namespace ElectronUtils


#if( BUILD_AS_ELECTRON_APP )

int interspec_start_server( const char *process_name, const char *userdatadir,
                            const char *basedir, const char *xml_config_path )
{
  try
  {
    //ToDo: refactor this userdatadir stuff into function inside InterSpecServer
    //      or InterSpecApp or something.
    if( !UtilityFunctions::create_directory(userdatadir) )
      throw std::runtime_error( "Failed to create directory '" + string(userdatadir) + "' for user data." );
  
    const string preffile = UtilityFunctions::append_path( userdatadir, "InterSpecUserData.db" );
    
    cout << "Will set user preferences file to: '" << preffile << "'" << endl;
    DataBaseUtils::setPreferenceDatabaseFile( preffile );
    DbToFilesystemLink::setFileNumToFilePathDBNameBasePath( userdatadir );
    const auto serial_db = UtilityFunctions::ls_files_in_directory( userdatadir, "serial_to_model.csv" );
    if( !serial_db.empty() )
      SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
#if( ENABLE_RESOURCE_UPDATES )
    ResourceUpdate::setUserDataDirectory( userdatadir );
#endif
    InterSpec::setWritableDataDirectory( userdatadir );
    
    cout << "Will make sure preferences database is up to date" << endl;
    DataBaseVersionUpgrade::checkAndUpgradeVersion();
    
#if( ENABLE_RESOURCE_UPDATES )
    ResourceUpdate::setupGlobalPrefsFromDb();
#endif
    
    InterSpec::setStaticDataDirectory( UtilityFunctions::append_path(basedir,"data") );
    
    //try
    //{
    //  boost::filesystem::current_path( basedir );
    //  cout << "Changed to cwd: " << basedir << endl;
    //}catch ( std::exception &e )
    //{
    //  cerr << "Unable to change to directory: '" << basedir << "' :" << e.what() << endl;
    //}
    
    InterSpecServer::startServerNodeAddon( process_name, basedir, xml_config_path );
  }catch( std::exception &e )
  {
    std::cerr << "\n\nCaught exception trying to start InterSpec server:\n\t"
    << e.what() << std::endl << std::endl;
    return -1;
  }
  
  return InterSpecServer::portBeingServedOn();
}//int interspec_start_server( int argc, char *argv[] )

void interspec_add_session_id( const char *session_id )
{
#warning "Need to make add_session_id() handle more than just single session ID"
  lock_guard<mutex> lock(ns_externalid_mutex);
  ns_externalid = session_id;
}//void interspec_add_session_id( const char *session_id )

int interspec_open_file( const char *sessionToken, const char *files_json )
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
  }//try / catch
    
    
  //Look for session with 'externalid' and open file...
  string externalid = sessionToken;
    
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
  
  return -1;
}

void interspec_kill_server()
{
  InterSpecServer::killServer();
}//void interspec_kill_server()


#endif //#if( BUILD_AS_ELECTRON_APP )
