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

#include <string>
#include <memory>

#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "InterSpec/MassAttenuationTool.h"
#include "target/electron/ElectronUtils.h"
#include "SpecUtils/SerialToDetectorModel.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

using namespace std;

namespace ElectronUtils
{

  
#if( USING_ELECTRON_NATIVE_MENU )
bool requestNewCleanSession()
{
  auto app = dynamic_cast<InterSpecApp *>(wApp);
  
  const string oldexternalid = app ? app->externalToken() : string();
  if( !oldexternalid.empty() )
  {
    //should check ns_externalid==oldexternalid
    string js;
    //Speed up loading by defering calls to Menu.setApplicationMenu() until app
    //  is fully reloaded.
    js += "$(window).data('HaveTriggeredMenuUpdate',null);";
    
    //Have electron reload the page.
    js += "ipcRenderer.send('NewCleanSession','" + app->externalToken() + "');";
    
    //Just in case the page reload doesnt go through, make sure menus will get updated eventually
    //  (this shouldnt be necassary, right?)
    js += "setTimeout(function(){$(window).data('HaveTriggeredMenuUpdate',true);},5000);";
    
    wApp->doJavaScript(js);
    return true;
  }else
  {
    cerr << "requestNewCleanSession(): failed; couldnt get external token." << endl;
  }

  return false;
}//void requestNewCleanSession()
#endif //USING_ELECTRON_NATIVE_MENU
  
bool notifyNodeJsOfNewSessionLoad()
{
  auto app = dynamic_cast<InterSpecApp *>(wApp);
  if( !app )
  {
    cerr << "Error: notifyNodeJsOfNewSessionLoad: wApp is null!!!" << endl;
    return false;
  }

  app->doJavaScript( "ipcRenderer.send('SessionFinishedLoading','" + app->externalToken() + "');" );
  app->triggerUpdate();
  
  return true;
}//bool notifyNodeJsOfNewSessionLoad( const std::string sessionid )
}//namespace ElectronUtils


#if( BUILD_AS_ELECTRON_APP )

int interspec_start_server( const char *process_name, const char *userdatadir,
                            const char *basedir, const char *xml_config_path )
{
  //Using a relative path should get us in less trouble than an absolute path
  //  on Windows.  Although havent yet tested (20190902) with network drives and such on Windows.
  //Should check if basedir is a relative or absolut path before this next step
//#warning "Need to check out how this setting basedir for Electron tartget works on WIndows network drives and such"
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
    DbToFilesystemLink::setFileNumToFilePathDBNameBasePath( userdatadir );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -3;
  }
  
  try
  {
    const auto serial_db = SpecUtils::ls_files_in_directory( userdatadir, "serial_to_model.csv" );
    if( !serial_db.empty() )
      SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -4;
  }
  
  
  try
  {
#if( ENABLE_RESOURCE_UPDATES )
    ResourceUpdate::setUserDataDirectory( userdatadir );
#endif
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
#if( ENABLE_RESOURCE_UPDATES )
    ResourceUpdate::setupGlobalPrefsFromDb();
#endif
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -7;
  }
  
  try
  {
    InterSpec::setStaticDataDirectory( SpecUtils::append_path(relbasedir,"data") );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    return -8;
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
    InterSpecServer::startServerNodeAddon( process_name, relbasedir, xml_config_path );
  }catch( std::exception &e )
  {
    std::cerr << "\n\nCaught exception trying to start InterSpec server:\n\t"
    << e.what() << std::endl << std::endl;
    return -9;
  }
  
  return InterSpecServer::portBeingServedOn();
}//int interspec_start_server( int argc, char *argv[] )


void interspec_add_allowed_session_token( const char *session_id )
{
  InterSpecServer::add_allowed_session_token( session_id );
}//void interspec_add_allowed_session_token( const char *session_id )


int interspec_remove_allowed_session_token( const char *session_token )
{
  return InterSpecServer::remove_allowed_session_token( session_token );
}//int interspec_remove_allowed_session_token( const char *session_token )


int interspec_open_file( const char *session_token, const char *files_json )
{
  return InterSpecServer::open_file_in_session( session_token, files_json );
}


bool interspec_using_electron_menus()
{
#if( USING_ELECTRON_HTML_MENU )
  return false;
#elif( USING_ELECTRON_NATIVE_MENU )
  return true;
#else
  return false;
#endif
}

void interspec_kill_server()
{
  InterSpecServer::killServer();
}//void interspec_kill_server()


#endif //#if( BUILD_AS_ELECTRON_APP )
