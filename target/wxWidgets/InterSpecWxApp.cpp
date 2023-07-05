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
#include "wx/config.h"
#include "wx/msgdlg.h"
#include "wx/cmdline.h"
#include "wx/filename.h"
#ifndef __APPLE__
#include "wx/snglinst.h"
#endif
#include "wx/stdpaths.h"
#include "wx/hyperlink.h"

#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Serializer>

#include "InterSpecWxApp.h"
#include "InterSpecWebFrame.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpecServer.h"
#include "InterSpec/UndoRedoManager.h"

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

class IpcConnection : public wxConnection
{
public:
  IpcConnection(void) : wxConnection() { }
  ~IpcConnection(void) { }

  virtual bool OnExec(const wxString& topic, const wxString& data)
  {
    // Server gets messages here when client uses the Execute
    //  TODO: send this information on to wxApp to actually use the information
    wxLogMessage("Have recieved message on server topic='%s', data='%s'", topic, data);


    assert( (topic == "FileToOpen") || (topic == "OpenNewWindow"));
    auto app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance());

    if (topic == "FileToOpen")
    {
      app->handle_open_file_message(data.utf8_string());
    }
    else if (topic == "OpenNewWindow")
    {
      app->new_app_window();
    }
    else
    {
      assert(0);
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

    parser.AddLongOption( "proxy",
      "Proxy configuration to use (only applicable to maps tool); valid values are:"
      " empty (defualt), 'direct', 'auto_detect', 'system', "
      "or any other string any other string will be interpreted as the 'proxyRules'.",
      wxCMD_LINE_VAL_STRING, wxCMD_LINE_PARAM_OPTIONAL );

    parser.AddParam("File or URI to open", wxCMD_LINE_VAL_STRING, wxCMD_LINE_PARAM_OPTIONAL | wxCMD_LINE_PARAM_MULTIPLE);
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

    //["/path/to/some/file","/some/other/file","interspec://drf/define?..."]
    // TODO: use Wt::JSON to make sure this message is proper JSON, etc
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
    auto frames_copy = m_frames;
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

    /*
    // TODO: Make a custom dialog with a URL they can click on to send an email bug report.
     wxHyperlinkCtrl *hyperlink = new wxHyperlinkCtrl( &dialog, -1,
      "Click here to email",
      "mailto:interspec@sandia.gov?subject=InterSpec%20bug:%20Fatal%20JS%20Error&body=" + Wt::Utils::urlEncode(error_msg)
      );
    wxDialog dialog( frame, -1, caption, wxDefaultPosition, wxDefaultSize, wxICON_ERROR | wxSTAY_ON_TOP, "dialogBox" );
    wxSizer *btnSizer = dialog.CreateButtonSizer( wxOK );
    wxSizer *txtSizer = dialog.CreateTextSizer( message )
    ...
    dialog.ShowModal();
    */

    if( frame )
      frame->Close();
    else
      close_all_windows_and_exit();


    wxLogMessage( "Have JS Error; msg='%s', session='%s'", error_msg, app_token );
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
    // Under Unix, the service name may be either an integer port identifier 
    //  in which case an Internet domain socket will be used for the communications, 
    //  or a valid file name (which shouldn't exist and will be deleted afterwards)
    //  in which case a Unix domain socket is created.
    wxString serverName = "InterSpecIPC";
#ifndef _WIN32
    serverName = "7072";
#endif

    // Make sure were the only instance running
    //  (is this per user, or per computer? Double check)
    m_checker = new wxSingleInstanceChecker();

    // By default wxWidgets uses the name `GetAppName() + '-' + wxGetUserId()` - however, 
    // the Electron version of the app uses the same thing, so we'll modify this one a 
    //  little by appending "-webview"
    const bool did_create = m_checker->Create( GetAppName() + '-' + wxGetUserId() + "-webview" );
    // `did_create` will be true, even if another instance of the program is running; it only indicated
    //  a failure to allocate a Windows named mutex or whatever (which I assume is excedingly rare to
    //  happen, but I didnt actually check)

    if( did_create && m_checker->IsAnotherRunning() )
    {
      // Here is where we would request the existing instance to open a file, or open a new window.
      //wxLogError(_("Another program instance is already running, aborting."));

      // Use https://docs.wxwidgets.org/3.0/overview_ipc.html to message other session

      delete m_checker; // OnExit() won't be called if we return false
      m_checker = nullptr;

      size_t n_valid_args = 0;
      Wt::Json::Array msg_json;
      for( size_t i = 0; i < m_command_line_args.size(); ++i )
      {
        std::string arg = m_command_line_args[i].utf8_string();

        if( SpecUtils::istarts_with( arg, "interspec:" ) 
          || SpecUtils::istarts_with( arg, "raddata://g" ) )
        {
          n_valid_args += 1;
        }else
        {
          // We could use SpecUtils to make canonical, but we'll use the wx stuff, just to try it out
          // arg = SpecUtils::make_canonical_path(arg)

          wxFileName fname( m_command_line_args[i] );
          if( fname.IsOk() && fname.IsFileReadable() )
            n_valid_args += 1;

          if( fname.IsOk() && fname.IsFileReadable() && !fname.IsAbsolute() )
            arg = fname.GetAbsolutePath().utf8_string();
        }//if( not a interspec:// URI )

        msg_json.push_back( Wt::Json::Value( Wt::WString::fromUTF8( arg ) ) );
      }//for (size_t i = 0; i < m_command_line_args.size(); ++i)

      // Right now we will open a new Window if there are no file or URI arguments,
      //  or we will open the files and/or URI.

      const char *topic = n_valid_args ? "FileToOpen" : "OpenNewWindow";
      std::string message = Wt::Json::serialize( msg_json, 0 );

      IpcClient client;
      IpcConnection *connection = (IpcConnection *)client.MakeConnection( "localhost", serverName, topic );

      connection->Execute( message.c_str(), message.size() + 1 );
      connection->Disconnect();
      delete connection;

      return false;  //OnExit wont be called.
    }//if( did_create && m_checker->IsAnotherRunning() )


    m_ipc_server = new IpcServer();


    if( !m_ipc_server->Create( "InterSpecIPC" ) )
      wxMessageBox( "Error", "Error creating IPC server" );

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
     // Right now we'll just wait 180 seconds, and if things havent loaded, declare failure.
     //  Normally shouldnt take very long, but on some test runners it can take a while.
     // When things load, we'll stop this timer and clear it.
     // 
     //  However, there are some hooks we could probably use to fail faster - see:
     //   void InterSpecApp::unload()
     //   void InterSpecApp::prepareForEndOfSession()
     //   void InterSpecApp::finalize()
     //   for possible hooks to check failure.
      sm_check_load_timer.reset( new CheckLoadTimer( this, 180 ) );
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
