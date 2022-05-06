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

#include <Wt/WServer>
#include <Wt/WIOService>
#include <Wt/WApplication>

#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "target/electron/ElectronUtils.h"
#include "target/electron/InterSpecAddOn.h"
#include "SpecUtils/SerialToDetectorModel.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

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
    
#if( USING_ELECTRON_NATIVE_MENU )
    //should check ns_externalid==oldexternalid
    string js;
    //Speed up loading by deferring calls to Menu.setApplicationMenu() until app
    //  is fully reloaded.
    js += "$(window).data('HaveTriggeredMenuUpdate',null);";
    
    //Just in case the page reload doesnt go through, make sure menus will get updated eventually
    //  (this shouldnt be necassary, right?)
    js += "setTimeout(function(){$(window).data('HaveTriggeredMenuUpdate',true);},5000);";
    
    wApp->doJavaScript(js);
#endif
    
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

  const string oldexternalid = app->externalToken();
  if( !oldexternalid.empty() )
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
  
  const string session_token = app->externalToken();
  
  Wt::WServer *server = Wt::WServer::instance();
  assert( server );
  
  Wt::WIOService &io = server->ioService();
  io.boost::asio::io_service::post( [=](){
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
  if( msg_name == "..." )
  {
    
  }
#if( !USE_ELECTRON_NATIVE_MENU )
  else if( msg_name == "OnMaximize" )
  {
    cout << "\n\nOnMaximize\n\n";
    app->doJavaScript( "Wt.WT.TitleBarChangeMaximized(true);" );
  }else if( msg_name == "OnUnMaximize" )
  {
    cout << "\n\nOnUnMaximize\n\n";
    app->doJavaScript( "Wt.WT.TitleBarChangeMaximized(false);" );
  }else if( msg_name == "OnBlur" )
  {
    app->doJavaScript( "$('.app-titlebar').addClass('inactive');" );
  }else if( msg_name == "OnFocus" )
  {
    app->doJavaScript( "$('.app-titlebar').removeClass('inactive');" );
  }else if( msg_name == "OnLeaveFullScreen" )
  {
    cout << "Left to fullscreen." << endl;
  }else if( msg_name == "OnEnterFullScreen" )
  {
    cout << "Went to fullscreen." << endl;
  }else
#endif
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

  Wt::WIOService &io = server->ioService();
  io.boost::asio::io_service::post( worker );
  
  return true;
}//bool browse_for_directory(...)

}//namespace ElectronUtils



int interspec_start_server( const char *process_name, const char *userdatadir,
                            const char *basedir, const char *xml_config_path )
{
  return InterSpecServer::start_server( process_name, userdatadir, basedir, xml_config_path );
}//int interspec_start_server( int argc, char *argv[] )


void interspec_set_require_session_token( const bool require_token )
{
  InterSpecServer::set_require_tokened_sessions( require_token );
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
