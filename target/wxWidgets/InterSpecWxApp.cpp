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
    assert(topic == "FileToOpen");
    auto app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance());
    app->handle_open_file_message(data.utf8_string());

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
    wxMessageBox(message, "File to Open");

    // m_active_frame
    // or
    // wxWindow* wxGetActiveWindow()

    //int InterSpecServer::open_file_in_session(const char* session_token, const char* files_json);
  }


  void InterSpecWxApp::handle_frame_closing(InterSpecWebFrame* frame)
  {
    const auto iter = std::find(std::begin(m_frames), std::end(m_frames), frame);
    assert(iter != std::end(m_frames));
    if (iter != std::end(m_frames))
      m_frames.erase(iter);
    if (m_active_frame == frame)
      m_active_frame = nullptr;
  }


  bool InterSpecWxApp::OnInit()
  {
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
      wxLogError(_("Another program instance is already running, aborting."));

      // Use https://docs.wxwidgets.org/3.0/overview_ipc.html to message other session

      delete m_checker; // OnExit() won't be called if we return false
      m_checker = NULL;


      IpcClient client;
      IpcConnection* connection = (IpcConnection*)client.MakeConnection("localhost", serverName, "FileToOpen");

      // TODO: replace this JSON encoding with using Wt::JSON or something.
      std::string message = "[";
      for (size_t i = 0; i < m_command_line_args.size(); ++i)
      {
        if (i)
          message += ",";

        std::string arg = m_command_line_args[i].utf8_string();
        if (!SpecUtils::make_canonical_path(arg))
        {
          arg = m_command_line_args[i].utf8_string();

          //std::string exe_path = ...
          //if (!make_canonical_path(arg, ))
          //  arg = m_command_line_args[i];
        }

        message += "\"" + arg + "\"";
      }
      message += "]";

      connection->Execute(message.c_str(), message.size());
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

    std::string xml_config_path = SpecUtils::append_path(base_dir, "data/config/wt_config_electron.xml" );


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