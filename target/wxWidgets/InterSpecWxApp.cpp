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


#ifdef WX_PRECOMP
#include "wx/wxprec.h"
#else
#include "wx/wx.h"
#endif


#include "wx/ipc.h"
#include "wx/timer.h"
#include "wx/utils.h"
#include "wx/config.h"
#include "wx/filefn.h"
#include "wx/msgdlg.h"
#include "wx/cmdline.h"
#include "wx/filename.h"
#ifndef __APPLE__
#include "wx/snglinst.h"
#endif
#include "wx/stdpaths.h"
#include "wx/hyperlink.h"
#include "wx/webview.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpecWxApp.h"
#include "InterSpecWebFrame.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/UndoRedoManager.h"

#include "InterSpecWxUtils.h"


/*
// Running batch commands not finished being implemented 
//  to work reasonably within wxWidgets.
#if( USE_BATCH_CLI_TOOLS )
#include "wx/textctrl.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/BatchCommandLine.h"
#endif
*/

namespace 
{
  bool sm_command_line_parsed = false;  //just for debugging to make sure I understand program flow

  // Some variable to hold command line switch options
#ifndef __APPLE__
  bool sm_single_instance = true;
#endif
  bool sm_try_restore = true;
  bool sm_test_load_only = false;
  bool sm_require_session_token = true;
  bool sm_open_dev_console = false;
  std::string sm_proxy_config;
  long sm_server_port = 0;
  long sm_max_runtime_seconds = 0;
/*
#if( USE_BATCH_CLI_TOOLS )
  bool sm_run_batch_command = false;
#endif
*/

  std::string sm_base_dir = ".";
  std::string sm_user_data_dir;


  bool sm_overide_rc = false;
  int sm_rc_override_value = 0;


  void set_data_dirs()
  {
    const wxStandardPaths &paths = wxStandardPaths::Get();

    sm_base_dir = ".";
    sm_user_data_dir = paths.GetUserDataDir().utf8_string();

    const std::string exe_path = paths.GetExecutablePath().utf8_string();
    const std::string exe_parent = SpecUtils::parent_path( exe_path );

    // TODO: use wxFileName to resolve links, etc

    const std::string test_file = SpecUtils::append_path( "InterSpec_resources", "InterSpec.css" );
    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = exe_parent;

    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = paths.GetResourcesDir().utf8_string();

    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = SpecUtils::parent_path( exe_parent );

    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = SpecUtils::parent_path( sm_base_dir );

    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = SpecUtils::parent_path( sm_base_dir );

    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = SpecUtils::parent_path( sm_base_dir );

    if( !SpecUtils::is_file( SpecUtils::append_path( sm_base_dir, test_file ) ) )
      sm_base_dir = ".";
  }//void set_data_dirs()


  class MaxRuntimeTimer : public wxTimer 
  {
    InterSpecWxApp *m_app;

  public:
    MaxRuntimeTimer( InterSpecWxApp *app, int num_seconds )
      : wxTimer(),
      m_app( app )
    {
      StartOnce( 1000 * num_seconds );
    }

    virtual void Notify() override
    {
      if( m_app )
        m_app->close_all_windows_and_exit();
    }
  };//class MaxRuntimeTimer



  class CheckLoadTimer : public wxTimer
  {
    InterSpecWxApp *m_app;

  public:
    CheckLoadTimer( InterSpecWxApp *app, int num_seconds )
      : wxTimer(),
      m_app( app )
    {
      StartOnce( 1000 * num_seconds );
    }

    virtual void Notify() override
    {
      sm_overide_rc = true;
      sm_rc_override_value = -11;

      if( m_app )
        m_app->close_all_windows_and_exit();
    }
  };//class CheckLoadTimer

  // TODO: these timers should probably be member variables of InterSpecWxApp
  std::unique_ptr<MaxRuntimeTimer> sm_max_runtime_timer;

  std::unique_ptr<CheckLoadTimer> sm_check_load_timer;
}//namespace


#ifndef __APPLE__
class IpcServer;

/** Returns a per-user, per-build IPC service name that both client and server
 must use the same value of, so the IPC actually connects.

 On Windows this is the DDE service name.  On Linux/Unix it is a filesystem
 path used as a Unix-domain socket - wxWidgets' IPC framework accepts either a
 numeric port or a filename for the service argument, and a filename avoids the
 multi-user / multi-build collisions of a fixed TCP port.
 */
static wxString compute_ipc_service_name()
{
  const wxString suffix = wxGetUserId()
                          + "-" + std::to_string( AppUtils::compile_date_as_int() );
#ifdef _WIN32
  return "InterSpecIPC-" + suffix;
#else
  // Prefer XDG_RUNTIME_DIR (per-user, tmpfs, auto-cleaned at logout); fall back
  // to /tmp.  The uid in the filename guards against shared /tmp collisions.
  wxString dir;
  if( !wxGetEnv( "XDG_RUNTIME_DIR", &dir ) || dir.IsEmpty() )
    dir = "/tmp";
  return dir + "/InterSpecIPC-" + suffix + ".sock";
#endif
}//wxString compute_ipc_service_name()


class IpcConnection : public wxConnection
{
public:
  IpcConnection(void) : wxConnection() { }
  ~IpcConnection(void) { }

  virtual bool OnExec(const wxString& topic, const wxString& data)
  {
    wxLogMessage("Received IPC message on server topic='%s', data='%s'", topic, data);

    InterSpecWxApp * const app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance());
    assert( app );
    if( !app )
      return false;

    if( topic == "FileToOpen" )
    {
      app->handle_open_file_message(data.utf8_string());
    }
    else if( topic == "OpenNewWindow" )
    {
      app->new_app_window();
    }
    else
    {
      // OnAcceptConnection should already have rejected unknown topics.
      assert( 0 );
      return false;
    }

    return true;
  }
};//IpcConnection


class IpcServer : public wxServer
{
public:
  IpcServer() : wxServer()
  {
  }

  virtual wxConnectionBase* OnAcceptConnection(const wxString& topic) wxOVERRIDE
  {
    if( (topic != "FileToOpen") && (topic != "OpenNewWindow") )
    {
      wxLogWarning( "InterSpec: rejecting IPC connection for unknown topic '%s'", topic );
      return nullptr;
    }
    return new IpcConnection();
  }
};


class IpcClient : public wxClient
{
public:
  IpcClient(void) { }
  wxConnectionBase* OnMakeConnection(void)
  {
    return new IpcConnection();
  }
};
#endif //#ifndef __APPLE__

InterSpecWxApp::InterSpecWxApp() :
    wxApp(),
    m_url(""),
    m_command_line_args{},
    m_frames{},
    m_active_frame(nullptr),
    m_checker(nullptr),
    m_ipc_server(nullptr)
  {
  }


  void InterSpecWxApp::OnInitCmdLine(wxCmdLineParser& parser)
  {
    wxApp::OnInitCmdLine(parser);

    parser.AddSwitch( "n", "no-restore", "Do not restore previous working state", wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddSwitch( "m", "mult-instance", 
      "Normally only one instance of the InterSpec executable runs, and when you open a new file or window,"
      " the original instance of InterSpec will actually open the file or new Window.\n"
      "This option removes that relationship so this instance of the application wont affect other instances, and vice-versa.", 
      wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddSwitch( "s", "no-token", 
      "Allow external web-browser sessions without a token.\n"
      "Normally, when the 'Use in external browser' option is chosen, a token is generated and used to control access.\n"
      "With this option enabled any user or program on your computer with access to `localhost` can then start an InterSpec"
      " session - note that the operating system does not allow external connections to `localhost`.", 
      wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddOption( "p", "port", 
      "Port to serve application on.\n"
      "Useful along with 'no-token' option to allow creating arbitrary number of sessions in your web-broswer.\n"
      "Value given should be above 1024, or if zero is specified (default), a random port will be chosen.", 
      wxCMD_LINE_VAL_NUMBER, wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddOption( "t", "max-run-time", 
      "Maximum number of seconds to run before exiting.\n"
      "Recomend to use in combination with the 'multi-instance' option.", 
      wxCMD_LINE_VAL_NUMBER, wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddSwitch( "l", "test-load", 
      "Test that the applicaiton successfully loads and then exits.", 
      wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddSwitch( "d", "dev", 
      "Opens html/js developer console.", 
      wxCMD_LINE_PARAM_OPTIONAL );

    parser.AddLongOption( "proxy",
      "Proxy configuration to use (only applicable to maps tool); valid values are:"
      " empty (defualt), 'direct', 'auto_detect', 'system', "
      "or any other string any other string will be interpreted as the 'proxyRules'.",
      wxCMD_LINE_VAL_STRING, wxCMD_LINE_PARAM_OPTIONAL );
    
    parser.AddParam("File or URI to open", wxCMD_LINE_VAL_STRING, wxCMD_LINE_PARAM_OPTIONAL | wxCMD_LINE_PARAM_MULTIPLE);

/*
#if( USE_BATCH_CLI_TOOLS )
    parser.AddLongSwitch( "batch-peak-fit", "Batch-fit peaks.", wxCMD_LINE_PARAM_OPTIONAL );
    parser.AddLongSwitch( "batch-act-fit", "Batch shielding/source fit.", wxCMD_LINE_PARAM_OPTIONAL );
#endif
*/
  }


  bool InterSpecWxApp::OnCmdLineParsed(wxCmdLineParser& parser)
  {
    if (!wxApp::OnCmdLineParsed(parser))
      return false;

    for (size_t i = 0; i < parser.GetParamCount(); ++i)
      m_command_line_args.push_back(parser.GetParam(i));

    
    if( parser.FoundSwitch( "test-load" ) == wxCMD_SWITCH_ON )
      sm_test_load_only = true;

#ifndef __APPLE__
    if( parser.FoundSwitch("mult-instance") == wxCMD_SWITCH_ON )
      sm_single_instance = false;
#endif

    if( parser.FoundSwitch( "no-restore" ) == wxCMD_SWITCH_ON )
      sm_try_restore = false;
  
    if( parser.FoundSwitch( "no-token" ) == wxCMD_SWITCH_ON )
      sm_require_session_token = false;

    if( parser.Found( "port", &sm_server_port ) )
    {
      if( (sm_server_port < 0) 
        || (sm_server_port > std::numeric_limits<unsigned short int>::max()) )
      {
        wxLogMessage( "Invalid negative server port number specified (%i); setting to zero.", sm_server_port );
        sm_server_port = 0;
      }
    }//if( parser.Found( "port", &port_val ) )

    wxString proxy;
    if( parser.Found( "proxy", &proxy ) )
      sm_proxy_config = proxy.utf8_string();

    if( parser.Found( "max-run-time", &sm_max_runtime_seconds ) && (sm_max_runtime_seconds > 0) )
      sm_max_runtime_timer.reset( new MaxRuntimeTimer( this, sm_max_runtime_seconds ) );

    if( parser.FoundSwitch( "dev" ) == wxCMD_SWITCH_ON )
      sm_open_dev_console = true;

/*
#if( USE_BATCH_CLI_TOOLS )
    if( (parser.FoundSwitch( "batch-peak-fit" ) == wxCMD_SWITCH_ON )
      || (parser.FoundSwitch( "batch-act-fit" ) == wxCMD_SWITCH_ON) )
    {
      sm_run_batch_command = true;
    }//if( a command-line batch run )
#endif
*/

    sm_command_line_parsed = true;

    return true;
  }//bool OnCmdLineParsed(wxCmdLineParser& parser)



  void InterSpecWxApp::handle_open_file_message(const std::string& message)
  {
    wxLogMessage("Will try to open file '%s'",message.c_str());

    // wxGetActiveWindow() doesnt usually work
    InterSpecWebFrame* active_frame = dynamic_cast<InterSpecWebFrame *>( wxGetActiveWindow() );
    
    if (!active_frame)
      active_frame = m_active_frame;

    if (!active_frame && !m_frames.empty())
      active_frame = m_frames.back();

    if (!active_frame)
    {
      wxMessageBox("Sorry, there was an issue finding an issue to open the spectrum file - sorry!", "Error Opening File");
      return;
    }

    const wxString& token = active_frame->app_token();

    // Expected payload: ["/path/to/some/file","/some/other/file","interspec://drf/define?..."]
    // Validate at this layer so a malformed message doesn't trigger a user-facing
    // "issue opening" dialog further down — it's a developer/IPC bug, not a file issue.
    try
    {
      const nlohmann::json parsed = nlohmann::json::parse( message );
      if( !parsed.is_array() )
        throw std::runtime_error( "IPC payload was not a JSON array" );
    }catch( std::exception &e )
    {
      wxLogError( "InterSpec: ignoring malformed IPC open-file payload: %s", e.what() );
      return;
    }

    const int status = InterSpecServer::open_file_in_session(token.utf8_string().c_str(), message.c_str());

    active_frame->Show();
    active_frame->Raise();

    if (status < 0)
    {
      wxLogMessage("Error opening spectrum file in session, code: %i", status);
      wxMessageBox("There may have been an issue opening the requested resource - sorry!", "Error Opening File/URI");
    }
  }//void handle_open_file_message(const std::string& message)


  void InterSpecWxApp::handle_frame_closing(InterSpecWebFrame* frame)
  {
    const auto iter = std::find(std::begin(m_frames), std::end(m_frames), frame);
    assert(iter != std::end(m_frames));
    if (iter != std::end(m_frames))
      m_frames.erase(iter);
    if (m_active_frame == frame)
      m_active_frame = nullptr;
  }


  void InterSpecWxApp::handle_frame_focus(InterSpecWebFrame* frame)
  {
    assert(frame);
    assert(std::find(std::begin(m_frames), std::end(m_frames), frame) != std::end(m_frames));
    m_active_frame = frame;
  }//void handle_frame_focus(InterSpecWebFrame* frame);


  void InterSpecWxApp::session_loaded( InterSpecWebFrame *frame )
  {
    if( sm_test_load_only )
    {
      sm_overide_rc = true;
      sm_rc_override_value = 0;
      if( sm_check_load_timer )
      {
        sm_check_load_timer->Stop();
        sm_check_load_timer.reset();
      }

      close_all_windows_and_exit();
    }//if( sm_test_load_only )
  }//void session_loaded( InterSpecWebFrame *frame );


  void InterSpecWxApp::close_all_windows_and_exit()
  {
    std::vector<InterSpecWebFrame *> frames_copy = m_frames;
    for( auto frame : frames_copy )
      frame->Close();

    if( frames_copy.empty() && wxAppConsoleBase::IsMainLoopRunning() )
      ExitMainLoop();
  }


  
  void InterSpecWxApp::handle_javascript_error( const std::string error_msg, const std::string app_token )
  {
    // static functions
    auto app = dynamic_cast<InterSpecWxApp *>(wxApp::GetInstance());
    assert( app );
    if( !app )
    {
      wxLogError( "InterSpecWxApp::handle_javascript_error: failed to get wxApp" );
      return;
    }


    // Pass this call/info to main UI thread
    wxWindow *topWindow = app->GetTopWindow();
    assert( topWindow );
    if( !topWindow )
    {
      wxLogError( "InterSpecWxApp::handle_javascript_error: failed to get top window" );
      return;
    }


    topWindow->GetEventHandler()->CallAfter( [=](){ 
      app->handle_javascript_error_internal( error_msg, app_token ); 
    } ); 
  }//void handle_javascript_error( const std::string &error_msg, const std::string app_token );

  const std::string &InterSpecWxApp::proxy_config()
  {
    return sm_proxy_config;
  }

  bool InterSpecWxApp::try_restore_session()
  {
    return sm_try_restore;
  }

  bool InterSpecWxApp::require_session_token()
  {
    return sm_require_session_token;
  }

  bool InterSpecWxApp::open_dev_console()
  {
    return sm_open_dev_console;
  }

  unsigned short int InterSpecWxApp::server_port()
  {
    return static_cast<unsigned short int>(sm_server_port);
  }


  void InterSpecWxApp::handle_javascript_error_internal( const std::string &error_msg, const std::string &app_token )
  {
    wxLogMessage( "Have JS Error; msg='%s', session='%s'", error_msg, app_token );

    std::cout << "JS Error: " << error_msg << std::endl;
    if( sm_test_load_only )
    {
      sm_overide_rc = true;
      sm_rc_override_value = -12;
      close_all_windows_and_exit();

      return;
    }//if( sm_test_load_only )


    InterSpecWebFrame *frame = nullptr;
    for( InterSpecWebFrame *f : m_frames )
    {
      if( f->app_token() == app_token )
      {
        frame = f;
        break;
      }
    }//for( InterSpecWebFrame *f : m_frames )


    if( frame )
      frame->Raise();
    
    wxString caption = "Javascript application error";
    wxString message = "There was a Javascript error - the application will now close.";
    
    wxMessageDialog dialog( GetTopWindow(), message, caption, wxOK | wxICON_ERROR | wxSTAY_ON_TOP );
    dialog.SetExtendedMessage( error_msg + 
      "\n\nPlease report to InterSpec@sandia.gov, along with what you were doing when this error occured." );
    dialog.ShowModal();

    // TODO: Make a custom dialog with a URL the user can click on to send an
    //       email bug report (mailto: with percent-encoded body).

    if( frame )
      frame->Close();
    else
      close_all_windows_and_exit();
  }//void handle_javascript_error_internal( const std::string &error_msg, const std::string &app_token );


  void InterSpecWxApp::new_app_window()
  {
    wxLogMessage("Creating new app window" );
    InterSpecWebFrame* frame = new InterSpecWebFrame(m_url, true, "");
    m_frames.push_back(frame);
    m_active_frame = frame;
    frame->Show();
    frame->Raise();
  }//void new_app_window()



//#include <fstream>
//  std::unique_ptr<std::ofstream> g_stdbuf, g_errbuf;


#ifndef __APPLE__
  bool InterSpecWxApp::check_single_instance()
  {
    const wxString service = compute_ipc_service_name();

    // Make sure were the only instance running.
    m_checker = new wxSingleInstanceChecker();

    // By default wxWidgets uses the name `GetAppName() + '-' + wxGetUserId()` - however,
    // the Electron version of the app uses the same thing, so we'll modify this one a
    //  little by appending build date, so this way if you want, you can have two different
    //  builds of InterSpec running.  The IPC service name above is keyed off the same
    //  per-user/per-build suffix so client and server always agree.
    const wxString app_name = GetAppName() + '-' + wxGetUserId()
                              + "-" + std::to_string(AppUtils::compile_date_as_int());
    const bool did_create = m_checker->Create( app_name );

    // `did_create` will be true, even if another instance of the program is running; it only indicated
    //  a failure to allocate a Windows named mutex or whatever (which I assume is excedingly rare to
    //  happen, but I didnt actually check)

    if( did_create && m_checker->IsAnotherRunning() )
    {
      size_t n_valid_args = 0;
      nlohmann::json msg_json = nlohmann::json::array();
      for( size_t i = 0; i < m_command_line_args.size(); ++i )
      {
        std::string arg = m_command_line_args[i].utf8_string();

        if( SpecUtils::istarts_with( arg, "interspec:" )
          || SpecUtils::istarts_with( arg, "raddata://g" ) )
        {
          n_valid_args += 1;
        }else
        {
          wxFileName fname( m_command_line_args[i] );
          if( fname.IsOk() && fname.IsFileReadable() )
            n_valid_args += 1;

          if( fname.IsOk() && fname.IsFileReadable() && !fname.IsAbsolute() )
            arg = fname.GetAbsolutePath().utf8_string();
        }//if( not a interspec:// URI )

        msg_json.push_back( arg );
      }//for (size_t i = 0; i < m_command_line_args.size(); ++i)

      // Open a new window if there are no file/URI arguments, otherwise hand
      // the files/URIs off to the running instance.
      const char *topic = n_valid_args ? "FileToOpen" : "OpenNewWindow";
      const std::string message = msg_json.dump();

      IpcClient client;
      IpcConnection *connection = static_cast<IpcConnection *>(
        client.MakeConnection( "localhost", service, topic ) );

      if( connection )
      {
        connection->Execute( message.c_str(), message.size() + 1 );
        connection->Disconnect();
        delete connection;

        delete m_checker; // OnExit() won't be called when we return false
        m_checker = nullptr;
        return false;
      }

      // Connection to the running instance failed - the other process is alive
      // (the single-instance checker says so) but its IPC channel is unreachable
      // (stale, crashed listener, or similar).  Fall through and start as an
      // independent instance so the user isn't locked out.  Server creation
      // below will likely also fail (the other instance still owns the name),
      // and that path is now also handled gracefully.
      wxLogError( "InterSpec: another instance is running but its IPC channel "
                  "is unreachable (service '%s'). Starting an independent instance.",
                  service );
    }//if( did_create && m_checker->IsAnotherRunning() )


    m_ipc_server = new IpcServer();

#ifndef _WIN32
    // On Unix the service name is a Unix-domain socket path.  A leftover socket
    // from a prior crashed run will make Create() fail.  The single-instance
    // checker has already confirmed no live peer (or we wouldn't be here), so
    // it's safe to clear a stale file.
    if( wxFileName::FileExists( service ) )
      wxRemoveFile( service );
#endif

    if( !m_ipc_server->Create( service ) )
    {
      wxLogError( "InterSpec: failed to create IPC server on '%s'; "
                  "open-file forwarding from secondary instances will not work.",
                  service );
      delete m_ipc_server;
      m_ipc_server = nullptr;
    }

    return true;
  }//bool check_single_instance();
#endif //#ifndef __APPLE__

  

  bool InterSpecWxApp::OnInit()
  {
   //   g_stdbuf.reset( new std::ofstream( "from_cout.txt") );
   //   g_errbuf.reset(new std::ofstream("from_cerr.txt") );
   //   std::cout.rdbuf(g_stdbuf->rdbuf());
   //   std::cerr.rdbuf(g_errbuf->rdbuf());
    
    //
    if( sm_command_line_parsed )
      throw std::runtime_error( "InterSpecWxApp::OnInit(): command line was already parsed before this function????" );

    set_data_dirs();

    
    InterSpecApp::setJavascriptErrorHandler( []( std::string errormsg, std::string app_token ){
      InterSpecWxUtils::handle_javascript_error( errormsg, app_token );
    } );

    InterSpecApp::setNativeFileSaveHandler( []( std::string data, std::string suggested_name ){
      InterSpecWxUtils::save_file_data( std::move( data ), std::move( suggested_name ) );
    } );

    // Register the wx native directory picker into LibInterSpec's
    // DirectorySelector.  Must happen on the executable side because
    // wxWidgets is statically linked separately into the dylib and the
    // executable, so wx static state (e.g. wxApp::ms_appInstance) does not
    // survive crossing that boundary.
    InterSpecWxUtils::register_native_directory_picker();


    try
    {
      const std::string data_dir = SpecUtils::append_path( sm_base_dir, "data" );
      const auto app_file_config = InterSpecServer::DesktopAppConfig::init( data_dir, sm_user_data_dir );
      sm_server_port = app_file_config.m_http_port;
      sm_proxy_config = app_file_config.m_proxy;
      sm_try_restore = app_file_config.m_allow_restore;
      sm_open_dev_console = app_file_config.m_open_dev_tools;
      sm_require_session_token = app_file_config.m_require_token;
      
      UndoRedoManager::setMaxUndoRedoSteps( app_file_config.m_max_undo_steps );
    }catch( std::exception &e )
    {
      wxLogMessage( "Error parsing app configuration file: %s", e.what() );
      wxMessageBox( e.what(), "Error parsing app configuration file" );
    }// try / catch


    if( !wxApp::OnInit() )
      return false; //OnExit wont be called.

    if( !sm_command_line_parsed )
      throw std::runtime_error( "InterSpecWxApp::OnInit(): command line not parsed after wxApp::OnInit()????" );

#ifdef NDEBUG
    const bool enableLogging = InterSpecWxApp::open_dev_console();
#else
    const bool enableLogging = true;
#endif

    if( enableLogging )
    {
      // Create a log window
      new wxLogWindow( nullptr, _( "Logging" ), true, false );
    }else
    {
      wxLog::EnableLogging( false );
    }

  /*
  // None of the methods I niavely tried to get the stdout and/or stderr to show up
  //  to the user works; should be able to use either wxLogStream or wxStreamToTextRedirector
  //  but neither works.  It should maybe also be possible to get the original console the
  //  app was launched from, and output the text there, but this looks more involved (modifying
  //  linking and such), and I'm not sure if its even doable.
  //  Perhaps there is a better way to do this, like start using libInterSpec as a shared library
  //  and just creating a script to run batch stuff.
#if( USE_BATCH_CLI_TOOLS )
    if( sm_run_batch_command )
    {
      int utf8_argc = 0;
      char **utf8_argv = nullptr;
      AppUtils::getUtf8Args( utf8_argc, utf8_argv );

      //wxLogWindow *batchLogger = new wxLogWindow( nullptr, _( "Batch" ), true, true );
      //wxLog *cout_logger = new wxLogStream( &std::cout );
      //batchLogger->SetActiveTarget( cout_logger );
      
      //wxTextCtrl *text = new wxTextCtrl( nullptr, -1, "Running batch command \n", wxDefaultPosition, wxSize(640,480), wxTE_MULTILINE | wxTE_READONLY | wxTE_NOHIDESEL );
      //wxStreamToTextRedirector redirect_cout( text, &std::cout );
      //new wxLogTextCtrl( text );

      const int rcode = BatchCommandLine::run_batch_command( utf8_argc, utf8_argv );

      AppUtils::cleanupUtf8Args( utf8_argc, utf8_argv );

      return true;
    }//if( a command-line batch run )
#endif
*/


#ifdef _WIN32
    // The Edge backend is compiled in (the #error directives in wxMain.cpp
    // and InterSpecWebFrame.cpp guarantee that), but the Microsoft Edge
    // WebView2 *Runtime* is a separate component - bundled with Win11, but
    // optional on Win10.  Detect a missing runtime up front and direct the
    // user to the installer rather than letting them face an opaque error
    // when the web view fails to construct.
    if( !wxWebView::IsBackendAvailable( wxWebViewBackendEdge ) )
    {
      const wxString download_url = "https://developer.microsoft.com/en-us/microsoft-edge/webview2/";
      wxMessageDialog dlg( nullptr,
        "InterSpec requires the Microsoft Edge WebView2 Runtime to display its user interface, "
        "but it doesn't appear to be installed on this computer.\n\n"
        "The runtime ships with Windows 11, but on older Windows 10 installs it must be "
        "installed separately.\n\n"
        "Would you like to open the download page now?",
        "Microsoft Edge WebView2 Runtime is required",
        wxYES_NO | wxICON_ERROR | wxSTAY_ON_TOP );
      dlg.SetExtendedMessage( "Download URL: " + download_url );
      if( dlg.ShowModal() == wxID_YES )
        wxLaunchDefaultBrowser( download_url );
      return false; //OnExit wont be called.
    }
#endif

#ifndef __APPLE__
    // Check if any other instance is running.
    //  If so, vmessage that instance and return from here.
    if( sm_single_instance && !check_single_instance() )
      return false; //OnExit wont be called.
#endif

    // Start Wt Server, etc
    //bool InterSpecServer::changeToBaseDir(int argc, char* argv[]);

    std::string xml_config_path = SpecUtils::append_path(sm_base_dir, "data/config/wt_config_wx.xml" );

    const int host_port = InterSpecServer::start_server( "InterSpec", 
                                                         sm_user_data_dir.c_str(), 
                                                         sm_base_dir.c_str(), 
                                                         xml_config_path.c_str(), 
                                                         sm_server_port );


    if( host_port <= 0 )
    {
      const wxStandardPaths &paths = wxStandardPaths::Get();

      wxMessageBox("Error starting InterSpec (" + std::to_string(host_port) + ")\n"
      + "\tUserDataDir: " + paths.GetUserDataDir() + "\n"
        + "\tExePath: " + paths.GetExecutablePath() + "\n"
        + "\tbase_dir: " + sm_user_data_dir,
        "Error"
      );
      return false;
    }

    assert( !sm_server_port || (sm_server_port == host_port) );
    sm_server_port = host_port;

    m_url = InterSpecServer::urlBeingServedOn();

    InterSpecServer::set_require_tokened_sessions( sm_require_session_token );
    //InterSpecServer::set_require_tokened_sessions(false);
    
    wxConfigBase* config = wxConfigBase::Get(true);
    const long num_load_attempts = config->ReadLong("/NumLoadAttempts", 0);
    config->Write("/NumLoadAttempts", num_load_attempts + 1);
    //wxLogMessage("Have attempted to load %i times.", static_cast<int>(num_load_attempts));


    const bool no_restore = ((!sm_try_restore) || (num_load_attempts >= 2));

    wxString file_to_open;
    for( size_t i = 0; i < m_command_line_args.size(); ++i )
    {
      wxFileName fname = wxFileName::FileName(m_command_line_args[i]);
      if (fname.IsOk() && fname.FileExists())
      {
        file_to_open = fname.GetAbsolutePath();
        wxLogMessage("Will open file: '%s'", file_to_open.utf8_string().c_str() );
        break;
      }
      
      std::string arg = m_command_line_args[i].utf8_string();
      if( SpecUtils::istarts_with( arg, "interspec://" )
        || SpecUtils::istarts_with( arg, "raddata://g" ) )
      {
        file_to_open = m_command_line_args[i];
          break;
      }
    }//for( size_t i = 0; i < m_command_line_args.size(); ++i )


    if( sm_test_load_only )
    {
     // Right now we'll just wait 90 seconds, and if things havent loaded, declare failure.
     //  Normally shouldnt take very long, but on some test runners it can take a while.
     // When things load, we'll stop this timer and clear it.
     // 
     //  However, there are some hooks we could probably use to fail faster - see:
     //   void InterSpecApp::unload()
     //   void InterSpecApp::prepareForEndOfSession()
     //   void InterSpecApp::finalize()
     //   for possible hooks to check failure.
      sm_check_load_timer.reset( new CheckLoadTimer( this, 90 ) );
    }//if( sm_test_load_only )

    InterSpecWebFrame* frame = new InterSpecWebFrame(m_url, no_restore, file_to_open);
    m_frames.push_back(frame);
    m_active_frame = frame;
    frame->Show();

    wxLogMessage("Serving InterSpec at: %s", m_url.c_str());

    if( sm_max_runtime_seconds > 0 )
      wxLogMessage( "Will limit execution to max time of %i seconds", sm_max_runtime_seconds );

    return true;
  }


  int InterSpecWxApp::OnRun()
  {
    const int rc = wxApp::OnRun();

    return sm_overide_rc ? sm_rc_override_value : rc;
  }



  int InterSpecWxApp::OnExit()
  {
    // We need to explicitly stop/get-rid-of these timers, or else 
    //  an exception will be thrown during shutdown, and our return
    //  value wont actually be returned.
    if( sm_max_runtime_timer )
    {
      sm_max_runtime_timer->Stop();
      sm_max_runtime_timer.reset();
    }

    if( sm_check_load_timer )
    {
      sm_check_load_timer->Stop();
      sm_check_load_timer.reset();
    }

    InterSpecServer::killServer();

    if( m_checker )
      delete m_checker;
    m_checker = nullptr;

    if( m_ipc_server )
      delete m_ipc_server;
    m_ipc_server = nullptr;

    return 0;
  }
