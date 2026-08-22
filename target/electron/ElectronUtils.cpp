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

#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SerialToDetectorModel.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

#include "target/electron/NativeMenu.h"
#include "target/electron/ElectronUtils.h"
#include "target/electron/InterSpecAddOn.h"


using namespace std;

namespace ElectronUtils
{

bool requestNewCleanSession()
{
  auto app = dynamic_cast<InterSpecApp *>(wApp);
  
  const string oldexternalid = app ? app->externalToken() : string();
  if( !oldexternalid.empty() )
  {
    //Have electron reload the page.
    ElectronUtils::send_nodejs_message("NewCleanSession", "");
    
    return true;
  }else
  {
    cerr << "requestNewCleanSession(): failed; couldnt get external token." << endl;
  }

  return false;
}//void requestNewCleanSession()

  
bool notifyNodeJsOfNewSessionLoad()
{
  auto app = dynamic_cast<InterSpecApp *>(wApp);
  if( !app )
  {
    cerr << "Error: notifyNodeJsOfNewSessionLoad: wApp is null!!!" << endl;
    return false;
  }

  // send_nodejs_message(...) is a no-op if there is no session token to address the message to.
  ElectronUtils::send_nodejs_message("SessionFinishedLoading", "");
  app->triggerUpdate();
  
  return true;
}//bool notifyNodeJsOfNewSessionLoad( const std::string sessionid )


void send_nodejs_message( const std::string msg_name, const std::string msg_data )
{
  auto app = dynamic_cast<InterSpecApp *>(wApp);
  if( !app )
  {
    cerr << "Error: send_nodejs_message: wApp is null!!!" << endl;
    return;
  }
  
  // The token is the address main.js routes on - it matches it against `window.appSessionToken`
  //  over its open windows - so an empty one could only produce a message matching no window.
  const string session_token = app->externalToken();
  if( session_token.empty() )
  {
    cerr << "Error: send_nodejs_message: no session token for '" << msg_name << "'." << endl;
    return;
  }
  
  Wt::WServer *server = Wt::WServer::instance();
  assert( server );
  if( !server )
  {
    cerr << "Error: send_nodejs_message: WServer::instance() is null!!!" << endl;
    return;
  }
  
  server->ioService().post( [=](){
    InterSpecAddOn::send_nodejs_message( session_token, msg_name, msg_data );
  } );
}//void send_nodejs_message(...)


bool handle_message_from_nodejs( const std::string &session_token,
                                const std::string &msg_name, const std::string &msg_data )
{
  InterSpecApp *app = InterSpecApp::instanceFromExtenalToken( session_token );
  
  if( !app )
  {
    // We will get here if the app-instance with this token hasnt yet loaded
    cerr << "Failed to find app instance for token='" << session_token << "'" << endl;
    
    //assert( 0 );
    
    return false;
  }//if( !app )
  
  Wt::WApplication::UpdateLock lock( app );
  if( !lock )
  {
    cerr << "Failed to get WApplication::UpdateLock lock token='" << session_token << "'" << endl;
    return false;
  }
  
// TODO: maybe we should make a own function for each of the cases below; maybe make things both
//       clearer here, as well as in main.js
  // The window-state messages only drive the HTML titlebar, which does not exist when the window
  //  has a native frame (macOS, or any build using the native menus), so the JS may not be loaded.
  if( msg_name == "OnMaximize" )
  {
    app->doJavaScript( "if(Wt.WT.TitleBarChangeMaximized)Wt.WT.TitleBarChangeMaximized(true);" );
  }else if( msg_name == "OnUnMaximize" )
  {
    app->doJavaScript( "if(Wt.WT.TitleBarChangeMaximized)Wt.WT.TitleBarChangeMaximized(false);" );
  }else if( msg_name == "OnBlur" )
  {
    app->doJavaScript( "var _tb=document.querySelector('.app-titlebar');if(_tb)_tb.classList.add('inactive');" );
  }else if( msg_name == "OnFocus" )
  {
    app->doJavaScript( "var _tb=document.querySelector('.app-titlebar');if(_tb)_tb.classList.remove('inactive');" );
#if( USING_ELECTRON_NATIVE_MENU )
  }else if( msg_name == "MenuItemClicked" )
  {
    handleElectronMenuActivation( app, msg_data );
#endif
  }else if( (msg_name == "OnEnterFullScreen") || (msg_name == "OnLeaveFullScreen") )
  {
    // main.js sends these, but the HTML titlebar does not currently do anything differently in
    //  full-screen; accepted so they do not log as unrecognized.
  }else
  {
    cerr << "Unrecognized msg_name from nodejs: '" << msg_name << "'" << endl;
    return false;
  }
  
  app->triggerUpdate();
  return true;
}//handle_message_from_nodejs(...)


bool browse_for_directory( const std::string &window_title,
                           const std::string &window_message,
                           std::function<void(std::string)> callback )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( Wt::WApplication::instance() );
  
  if( !app )
    throw runtime_error( "ElectronUtils::browse_for_directory(): must be called from within Wt event-loop." );
  
  if( !InterSpecApp::isPrimaryWindowInstance() )
  {
    cerr << "Browse for directory should only be called from a primary instance\n";
    assert( 0 );
    return false;
  }
  //session_token
  assert( callback );
  if( !callback )
    return false;
  
  const string session_id = app->sessionId();
  
  std::function<void(std::string)> wrapped_callback = [session_id,callback](string result_path){
    Wt::WServer *server = Wt::WServer::instance();
    if( !server ){
      cerr << "browse_for_directory callback wrapper: WServer no longer available." << endl;
      assert( 0 );
      return;
    }
    
    server->post(session_id, [callback,result_path](){
      Wt::WApplication *app = Wt::WApplication::instance();
      if( !app )
        return;
      
      callback( result_path );
      app->triggerUpdate();
    });
  };//wrapped_callback(...)
  
  const string token = app->externalToken();
  
  auto worker = [=](){
    InterSpecAddOn::browse_for_directory( token, window_title, window_message, wrapped_callback );
  };
  
  Wt::WServer *server = Wt::WServer::instance();
  assert( server );

  //Deliberately the raw asio post rather than WIOService::post: the latter dispatches through
  //  WIOService's strand, and `worker` blocks for as long as the user has the directory dialog
  //  open - which would stall every other posted job (including send_nodejs_message) until they
  //  dismiss it.
  Wt::WIOService &io = server->ioService();
  io.boost::asio::io_service::post( worker );
  
  return true;
}//bool browse_for_directory(...)

}//namespace ElectronUtils



int interspec_start_server( const char *process_name, const char *userdatadir,
                            const char *basedir, const char *xml_config_path,
                            const unsigned short int server_port_num )
{
  return InterSpecServer::start_server( process_name, userdatadir, basedir, xml_config_path, server_port_num );
}//int interspec_start_server( int argc, char *argv[] )


void interspec_set_require_session_token( const bool require_token )
{
  InterSpecServer::set_require_tokened_sessions( require_token );
}


void interspec_set_max_undo_steps( const int max_items )
{
  UndoRedoManager::setMaxUndoRedoSteps( max_items );
}


void interspec_add_allowed_primary_session_token( const char *session_token )
{
  InterSpecServer::add_allowed_session_token( session_token, InterSpecServer::SessionType::PrimaryAppInstance );
}//void interspec_add_allowed_primary_session_token( const char *session_id )


void interspec_add_allowed_external_session_token( const char *session_token )
{
  InterSpecServer::add_allowed_session_token( session_token, InterSpecServer::SessionType::ExternalBrowserInstance );
}


int interspec_remove_allowed_session_token( const char *session_token )
{
  return InterSpecServer::remove_allowed_session_token( session_token );
}//int interspec_remove_allowed_session_token( const char *session_token )


int interspec_session_is_alive( const char *session_token )
{
  const int status = InterSpecServer::session_status( session_token );
  
  return (status == 2);
}

int interspec_open_file( const char *session_token, const char *files_json )
{
  return InterSpecServer::open_file_in_session( session_token, files_json );
}

bool interspec_open_app_url( const char *session_token, const char *files_json )
{
  return InterSpecServer::pass_app_url_to_session( session_token, files_json );
}

bool interspec_set_initial_file_to_open( const char *session_token, const char *file_path )
{
  try
  {
    InterSpecServer::set_file_to_open_on_load( session_token, file_path );
  }catch( std::exception &e )
  {
    cerr << "interspec_set_initial_file_to_open: " << e.what() << endl;
    return false;
  }
  
  return true;
}


bool interspec_using_electron_menus()
{
#if( USING_ELECTRON_NATIVE_MENU )
  return true;
#else
  return false;
#endif
}


void interspec_kill_server()
{
  InterSpecServer::killServer();
}//void interspec_kill_server()
