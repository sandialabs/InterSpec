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
#include "wx/config.h"
#include "wx/msgdlg.h"
#include "wx/cmdline.h"
#include "wx/filename.h"
#include "wx/snglinst.h"
#include "wx/stdpaths.h"

#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Serializer>

#include "InterSpecWxApp.h"
#include "InterSpecWebFrame.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpecServer.h"

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






InterSpecWxApp::InterSpecWxApp() :
    wxApp(),
    m_url(""),
    m_checker(nullptr),
    m_ipc_server(nullptr)
  {
  }

  void InterSpecWxApp::OnInitCmdLine(wxCmdLineParser& parser)
  {
    wxApp::OnInitCmdLine(parser);

    parser.AddParam("File or URI to open",
      wxCMD_LINE_VAL_STRING,
      wxCMD_LINE_PARAM_OPTIONAL | wxCMD_LINE_PARAM_MULTIPLE);

    // Could put 
    //  - port number
    //  - wehter not to retore or not
    //  - allow multiple EXE instances
    //  - 
  }

  bool InterSpecWxApp::OnCmdLineParsed(wxCmdLineParser& parser)
  {
    if (!wxApp::OnCmdLineParsed(parser))
      return false;

    for (size_t i = 0; i < parser.GetParamCount(); ++i)
      m_command_line_args.push_back(parser.GetParam(i));

    return true;
  }



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


  void InterSpecWxApp::new_app_window()
  {
    wxLogMessage("Creating new app window" );
    InterSpecWebFrame* frame = new InterSpecWebFrame(m_url, true, "");
    m_frames.push_back(frame);
    m_active_frame = frame;
    frame->Show();
    frame->Raise();
  }//void new_app_window()


  namespace
  {
    /*
    class androidbuf : public std::streambuf
    {
      //A utility to redirect cout/cerr to the Android logging system so it can be
      // seen using 'adb logcat'
    public:
      enum Source { FromCout, FromCerr };
      enum { bufsize = 128 }; // ... or some other suitable buffer size
      androidbuf(Source src)
        : m_type(src), m_origSrc(nullptr)
      {
        this->setp(buffer, buffer + bufsize - 1);

        switch (m_type)
        {
        case FromCout:
          m_source = "cout";
          m_origSrc = std::cout.rdbuf(this);
          break;

        case FromCerr:
          m_source = "cerr";
          m_origSrc = std::cerr.rdbuf(this);
          break;
        }//switch( src )
      }//androidbuf( Source src )

      ~androidbuf()
      {
        //cout/cerr must be given back there original rdbufs or else there can be a
        //  problems with freeing resources
        if (!m_origSrc)
          return;
        switch (m_type)
        {
        case FromCout: std::cout.rdbuf(m_origSrc); break;
        case FromCerr: std::cerr.rdbuf(m_origSrc); break;
        }
      }//~androidbuf()

    private:
      int overflow(int c) override
      {
        if (c == traits_type::eof()) {
          *this->pptr() = traits_type::to_char_type(c);
          this->sbumpc();
        }
        return this->sync() ? traits_type::eof() : traits_type::not_eof(c);
      }//int overflow(int c)

      int sync() override
      {
        int rc = 0;
        if (this->pbase() != this->pptr())
        {
          rc = __android_log_write(ANDROID_LOG_INFO, m_source,
            std::string(this->pbase(), this->pptr()).c_str());
          this->setp(buffer, buffer + bufsize - 1);
        }
        return rc;
      }//int sync()
      char buffer[bufsize];
      const char* m_source;
      const Source m_type;
      std::streambuf* m_origSrc;
    };//class androidbuf

    */
    
  }
//#include <fstream>
//  std::unique_ptr<std::ofstream> g_stdbuf, g_errbuf;


  bool InterSpecWxApp::OnInit()
  {
   //   g_stdbuf.reset( new std::ofstream( "from_cout.txt") );
   //   g_errbuf.reset(new std::ofstream("from_cerr.txt") );
   //   std::cout.rdbuf(g_stdbuf->rdbuf());
   //   std::cerr.rdbuf(g_errbuf->rdbuf());

    if (!wxApp::OnInit())
      return false;

    // Under Unix, the service name may be either an integer port identifier in which case an Internet domain socket will be used for the communications, or a valid file name (which shouldn't exist and will be deleted afterwards) in which case a Unix domain socket is created.
    wxString serverName = "InterSpecIPC";
#ifndef _WIN32
    serverName = "7072";
#endif

    // Make sure were the only instance running
    //  (is this per user, or per computer? Double check)
    m_checker = new wxSingleInstanceChecker();
    if (m_checker->IsAnotherRunning())
    {
      // Here is where we would request the existing instance to open a file, or open a new window.
      //wxLogError(_("Another program instance is already running, aborting."));

      // Use https://docs.wxwidgets.org/3.0/overview_ipc.html to message other session

      delete m_checker; // OnExit() won't be called if we return false
      m_checker = NULL;

      size_t n_valid_args = 0;
      Wt::Json::Array msg_json;
      for (size_t i = 0; i < m_command_line_args.size(); ++i)
      {
        std::string arg = m_command_line_args[i].utf8_string();
 
        if (SpecUtils::istarts_with(arg, "interspec:"))
        {
          n_valid_args += 1;
        }else
        {
          // We could use SpecUtils to make canonical, but we'll use the wx stuff, just to try it out
          // arg = SpecUtils::make_canonical_path(arg)
          
          wxFileName fname(m_command_line_args[i]);
          if( fname.IsOk() && fname.IsFileReadable() )
            n_valid_args += 1;

          if (fname.IsOk() && fname.IsFileReadable() && !fname.IsAbsolute())
            arg = fname.GetAbsolutePath().utf8_string();
        }//if( not a interspec:// URI )
          
        msg_json.push_back(Wt::Json::Value(Wt::WString::fromUTF8(arg)));
      }//for (size_t i = 0; i < m_command_line_args.size(); ++i)

      // Right now we will open a new Window if there are no file or URI arguments,
      //  or we will open the files and/or URI.
      //  TODO: we should probably allow some command line arguments to avoid this single application stuff, or maybe some other options.

      const char* topic = n_valid_args ? "FileToOpen" : "OpenNewWindow";
      std::string message = Wt::Json::serialize(msg_json, 0);

      IpcClient client;
      IpcConnection* connection = (IpcConnection*)client.MakeConnection("localhost", serverName, topic);

      connection->Execute(message.c_str(), message.size() + 1);
      connection->Disconnect();
      delete connection;

      return false;
    }


    m_ipc_server = new IpcServer();


    if (!m_ipc_server->Create("InterSpecIPC"))
      wxMessageBox("Error", "Error creating IPC server");


    // Start Wt Server, etc
    const wxStandardPaths& paths = wxStandardPaths::Get();

    std::string base_dir = ".";
    std::string user_data_dir = paths.GetUserDataDir().utf8_string();

    const std::string exe_path = paths.GetExecutablePath().utf8_string();
    const std::string exe_parent = SpecUtils::parent_path(exe_path);

    // TODO: use wxFileName to resolve links, etc

    const std::string test_file = SpecUtils::append_path("InterSpec_resources", "InterSpec.css");
    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)) )
      base_dir = exe_parent;

    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)))
      base_dir = paths.GetResourcesDir().utf8_string();
    
    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)))
      base_dir = SpecUtils::parent_path(exe_parent);

    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)))
      base_dir = SpecUtils::parent_path(base_dir);

    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)))
      base_dir = SpecUtils::parent_path(base_dir);

    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)))
      base_dir = SpecUtils::parent_path(base_dir);

    if (!SpecUtils::is_file(SpecUtils::append_path(base_dir, test_file)))
      base_dir = ".";

    //bool InterSpecServer::changeToBaseDir(int argc, char* argv[]);

    std::string xml_config_path = SpecUtils::append_path(base_dir, "data/config/wt_config_wx.xml" );

    const int host_port = InterSpecServer::start_server("InterSpec", user_data_dir.c_str(), base_dir.c_str(), xml_config_path.c_str());


    if (host_port <= 0)
    {
      wxMessageBox("Error starting InterSpec (" + std::to_string(host_port) + ")\n"
      + "\tUserDataDir: " + paths.GetUserDataDir() + "\n"
        + "\tExePath: " + paths.GetExecutablePath() + "\n"
        + "\tbase_dir: " + base_dir,
        "Error"
      );
      return false;
    }

    m_url = InterSpecServer::urlBeingServedOn();

    InterSpecServer::set_require_tokened_sessions(true);
    //InterSpecServer::set_require_tokened_sessions(false);
    
    wxConfigBase* config = wxConfigBase::Get(true);
    const long num_load_attempts = config->ReadLong("/NumLoadAttempts", 0);
    config->Write("/NumLoadAttempts", num_load_attempts + 1);
    //wxLogMessage("Have attempted to load %i times.", static_cast<int>(num_load_attempts));

    const bool no_restore = (num_load_attempts >= 2);

    wxString file_to_open;
    for (size_t i = 0; i < m_command_line_args.size(); ++i)
    {
      wxFileName fname = wxFileName::FileName(m_command_line_args[i]);
      if (fname.IsOk() && fname.FileExists())
      {
        file_to_open = fname.GetAbsolutePath();
        wxLogMessage("Will open file: '%s'", file_to_open.utf8_string().c_str() );
        break;
      }
    }

    InterSpecWebFrame* frame = new InterSpecWebFrame(m_url, no_restore, file_to_open);
    m_frames.push_back(frame);
    m_active_frame = frame;
    frame->Show();

    wxLogMessage("Serving InterSpec at: %s", m_url.c_str());

    return true;
  }


  int InterSpecWxApp::OnExit()
  {
    InterSpecServer::killServer();


    if (m_checker)
      delete m_checker;
    m_checker = nullptr;

    if (m_ipc_server)
      delete m_ipc_server;
    m_ipc_server = nullptr;

    return 0;
  }