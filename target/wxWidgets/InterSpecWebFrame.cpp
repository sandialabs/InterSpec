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

#include <string>
#include <random>
#include <numeric>

#include <boost/system/error_code.hpp>

#include <Wt/WServer>
#include <Wt/WIOService>
#include <Wt/Http/Client>
#include <Wt/Http/Message>

#include "wx/fileconf.h"
#include "wx/artprov.h"
#include "wx/cmdline.h"
#include "wx/snglinst.h"
#include "wx/notifmsg.h"
#include "wx/settings.h"
#include <wx/sizer.h>
#include "wx/webview.h"
#include <wx/utils.h>

#include "wx/numdlg.h"
#include "wx/infobar.h"
#include "wx/fs_arc.h"
#include "wx/fs_mem.h"
#include "wx/stdpaths.h"
#include "wx/uri.h"
#include "wx/url.h"
#include "wx/filedlg.h"
#include "wx/wfstream.h"

#ifdef _WIN32
#include <wx/mimetype.h>
#include "wx/msw/registry.h"

#if wxUSE_WEBVIEW_EDGE
#include "wx/msw/webview_edge.h"
#else
#error "Must use Edge web-view for compiling on windows"
#endif
#endif

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpecServer.h"

#include "InterSpecWxApp.h"
#include "InterSpecWebFrame.h"

 
namespace
{
#ifdef _WIN32
  void check_url_association()
  {
    // We want to associate "interspec://" URLs with the application; we'll
    //  do this through editing the registry.
    
    // For files, we would use the wxFileType, and wxMimeTypesManager classes, similar to:
    //wxFileTypeInfo fti("", app_path, "", "Gamma Radiation Spectrum File", "n42", "pcf");
    //fti.SetShortDesc("Spectrum file"));
    //wxFileType* ft = wxTheMimeTypesManager->Associate(fti);
    //ft->SetDefaultIcon(app_path, 1);
    //delete ft;

    // Only do this function the first time an application loads
    static bool s_already_called = false;
    if (s_already_called)
      return;
    s_already_called = true;

    //wxFileName exe_filename(wxStandardPaths::Get().GetExecutablePath());
    //wxString app_path(exe_filename.GetFullPath());
    const wxString app_path = wxStandardPaths::Get().GetExecutablePath();
    
    const wxString open_command = wxString::Format("\"%s\" \"%%1\"", app_path);

    /*
    {// begin play around with file extensions
    // Right now this all kinda-sorta-seems to do something, about right, but
    //  I'm not sure of the correctness, or if it messes up other 
    //  programs/user-preferences, or whatever, so not bringing out to the
    //  user at the moment, pending further looking into other ways to do
    //  this, better (maybe the solution is to just use an installer, and
    //  have a seperate portable distribution that includes a script or 
    //  something)
      wxMimeTypesManager manager;
      wxFileType* ft = manager.GetFileTypeFromExtension("n42");
      if (!ft)
      {
        wxLogMessage("Failed to get file type information from extension.");
      }else
      {
        wxString mimeType;
        ft->GetMimeType(&mimeType);
        const wxString open_cmd = ft->GetOpenCommand("MyFile.n42");
        wxLogMessage("Command for openeing files was: '%s', mimetype -> '%s'", open_cmd, mimeType);
        if( manager.Unassociate(ft) )
          wxLogMessage("Unassociated previous file type association");
        else
          wxLogMessage("Failed to unassociate previous file type association");

        delete ft;
        ft = nullptr;
      }
      
      wxFileTypeInfo type_info("application/x.gamma-spectrum");
      type_info.SetOpenCommand(open_command);
      type_info.SetShortDesc("N42");
      type_info.SetDescription("Gamma energy spectrum files");
      type_info.AddExtension("n42");
      type_info.AddExtension("cnf"); //"pcf", ...
      
      ft = manager.Associate(type_info);
      if (ft)
        wxLogMessage("Added file type association; command for opening files: '%s'", ft->GetOpenCommand("myfile.n42"));
      else
        wxLogMessage("Failed to add ");

      delete ft;
      ft = nullptr;
    }// end play around with file extensions
    */

    wxRegKey interspec_key(wxRegKey::HKCU, "SOFTWARE\\Classes\\interspec");
    if (!interspec_key.Exists() && !interspec_key.Create(true))
    {
      wxLogError("Failed to create InterSpec key - aboring");
      return;
    }

    wxString def_value, proto_value, command_value;
    interspec_key.QueryValue("", def_value); //Its okay if this fails

    if( (def_value != "URL:interspec") && !interspec_key.SetValue("", "URL:interspec") )
    {
      wxLogError("Failed to set Reg key default value - arborting");
      return;
    }

    if (!interspec_key.QueryValue("URL Protocol", proto_value) 
        && !interspec_key.SetValue("URL Protocol", "") )
    {
      wxLogError("Failed to set Reg key 'URL Protocol' value");
      return;
    }      
    
    wxRegKey command_key(wxRegKey::HKCU, "SOFTWARE\\Classes\\interspec\\shell\\open\\command");
    if (!command_key.Exists() && !command_key.Create(true))
    {
      wxLogMessage("Failed to create 'command' Reg key");
      return;
    }

    command_key.QueryValue("", command_value); //Its okay if this fails
    
    if (command_value == open_command)
    {
      wxLogMessage("No need to update registry URL command value");
      return;
    }
    

    if( !command_key.SetValue("", open_command) )
      wxLogMessage("Failed to update URL command Reg key def value");
    else
      wxLogMessage("Updated URL 'command' Reg key def value to '%s'", open_command);
  }//void check_url_association()
#endif
}//namespace


InterSpecWebFrame::InterSpecWebFrame(const wxString& url, const bool no_restore, const wxString& file_to_open) :
  wxFrame(NULL, wxID_ANY, "InterSpec", wxDefaultPosition, wxDefaultSize, wxRESIZE_BORDER | wxSYSTEM_MENU),
  m_url( url ),
  m_token( "" ),
  m_dragging_window( false ),
  m_mouse_down_pos( 0, 0 ),
  m_currently_maximized( false )
{
  // set the frame icon
  //SetIcon(wxICON(sample));
  
  m_token = generate_token();
    
  const int add_token_status = InterSpecServer::add_allowed_session_token(m_token.c_str(), InterSpecServer::SessionType::PrimaryAppInstance);
  assert(add_token_status == 0);

  if (!file_to_open.empty())
    InterSpecServer::set_file_to_open_on_load(m_token.c_str(), file_to_open.utf8_string());

  
  //SetTitle("InterSpec");
  EnableFullScreenView(); // Enable native fullscreen API on macOS

  wxBoxSizer* topsizer = new wxBoxSizer(wxVERTICAL);

#if wxUSE_WEBVIEW_EDGE
  // Check if a fixed version of edge is present in
  // $executable_path/edge_fixed and use it
  wxFileName edgeFixedDir(wxStandardPaths::Get().GetExecutablePath());
  edgeFixedDir.SetFullName("");
  edgeFixedDir.AppendDir("edge_fixed");
  if (edgeFixedDir.DirExists())
  {
    wxWebViewEdge::MSWSetBrowserExecutableDir(edgeFixedDir.GetFullPath());
    wxLogMessage("Using fixed edge version");
  }

  //void* wxWebViewEdge::GetNativeBackend() const
#endif
  // Create the webview
  m_browser = wxWebView::New();
  
  wxString app_url = m_url + "?apptoken=" + m_token;

  // There is a bug if we dont use WebSockets, when we DONT restore state, the UI freezes about half the time; doesnt seem to ever happen when we restore a state... I have no idea
  //wxLogMessage("Not restoring state for debug purposes" );
  //app_url += "&restore=no";

  if (no_restore)
    app_url += "&restore=no";

  // The user agent isnt terrible important for us, so we'll just use a default Chrome one on windows
  //  TODO: customize for macOS
  wxString user_agent = "Mozilla / 5.0 (Windows NT 10.0; Win64; x64) AppleWebKit / 537.36 (KHTML, like Gecko) Chrome / 106.0.0.0 Safari / 537.36";

#ifdef __VISUALC__
  // If the users computer is configured to use a proxy, or auto-detect a proxy, the initial load may
  //  take quite a while, especially if they arent on a proxy (I'm guessing for things to time out).  
  // However, all content is local, (the exception to this is the google maps feature, which currently 
  // isnt supported) so we never need a proxy, so we'll just not use one.
  // To achieve this, we'll hijack the user agent string, becuase in the wxWidgets source file webview_edge.cpp,
  // the user agent is set via:
  //   options->put_AdditionalBrowserArguments( wxString::Format("--user-agent=\"%s\"", m_customUserAgent).wc_str());
  // So we'll just hijack the quotes, since I dont see another way to set browser arguments.
  user_agent += "\" --no-proxy-server";
  
  // I briefly tried using a proxy bypass list, but it didnt seem to work.
  //user_agent += "\" --proxy-bypass-list=\"127.0.0.1;localhost\"";
  //user_agent += "\" --proxy-auto-detect --proxy-bypass-list=\"127.0.0.1;localhost\"";
  //user_agent += "\" --proxy-server=direct://"; //not tested
  //user_agent += "\" --proxy-pac-url=..."; //not tested
#endif


  m_browser->SetUserAgent(user_agent);


  m_browser->Create(this, wxID_ANY, app_url, wxDefaultPosition, wxDefaultSize);
  topsizer->Add(m_browser, wxSizerFlags().Expand().Proportion(1));

  
  // Log backend information
  wxLogMessage("Backend: %s Version: %s", m_browser->GetClassInfo()->GetClassName(),
    wxWebView::GetBackendVersionInfo().ToString());
  wxLogMessage("User Agent: %s", m_browser->GetUserAgent());


  if (!m_browser->AddScriptMessageHandler("wx"))
    wxLogError("Could not add script message handler");

  SetSizer(topsizer);

  //We'll set the window size based on last-closed window, or some sensible defaults.
  //  TODO: check if any other window is open, and if so size and position relative to it.
  const long screen_x = wxSystemSettings::GetMetric(wxSYS_SCREEN_X);
  const long screen_y = wxSystemSettings::GetMetric(wxSYS_SCREEN_Y);
  wxLogMessage("Screen is: %i x %i", screen_x, screen_y);

  //wxFileName pref_file = wxFileConfig::GetLocalFile("InterSpec_app_state");
//wxFileConfig config("InterSpec", "Sandia", pref_file.GetFullPath(), "", wxCONFIG_USE_LOCAL_FILE);

  const wxSize min_win_size = FromDIP(wxSize(800, 600));
  const wxSize max_win_size = FromDIP(wxSize(1600, 1024));
  const long def_width = std::min(static_cast<long>(max_win_size.x), static_cast<long>(0.85 * screen_x));
  const long def_height = std::min(static_cast<long>(max_win_size.y), static_cast<long>(0.85 * screen_y));

  wxConfigBase* config = wxConfigBase::Get(true);
  long ww = config->ReadLong("/WindowWidth", def_width);
  long wh = config->ReadLong("/WindowHeight", def_height);
  long wx = config->ReadLong("/WindowX", static_cast<long>(0.025 * screen_x));
  long wy = config->ReadLong("/WindowY", static_cast<long>(0.025 * screen_y));
  
  wxLogMessage("Prev window: %i x %i at %i,%i", ww, wh, wx, wy);
  
  ww = std::min( std::max(ww, static_cast<long>(min_win_size.x)), static_cast<long>(screen_x) );
  wh = std::min( std::max(wh, static_cast<long>(min_win_size.y)), static_cast<long>(screen_y) );
  wx = std::max( 0l, std::min(wx, screen_x / 2) );
  wy = std::max( 0l, std::min(wy, screen_y / 2) );

  wxLogMessage("Will set window tp: %i x %i at %i,%i", ww, wh, wx, wy);

  //I'm not sure how/if screen_x and screen_y take into account GetDPIScaleFactor()
  //SetSize(FromDIP(wxSize(0.85*screen_x, 0.85*screen_y)));
  SetSize( wxSize(ww, wh) );
  SetPosition( wxPoint(wx, wy) );
    

  // Allow the right-click context menu on web page contents
  m_browser->EnableContextMenu(true);

  // Enable the "Inspect" option in the right-click context menu
  m_browser->EnableAccessToDevTools(true);

  //wxString customUserAgent = "Mozilla/5.0 (iPhone; CPU iPhone OS 13_1_3 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/13.0.1 Mobile/15E148 Safari/604.1";
  //if (!m_browser->SetUserAgent(customUserAgent))
  //  wxLogError("Could not set custom user agent");

  // Not sure what the zoom type actually means
  //m_browser->SetZoomType(wxWEBVIEW_ZOOM_TYPE_LAYOUT); //Scales entire page, including images
  //m_browser->SetZoomType(wxWEBVIEW_ZOOM_TYPE_TEXT);  // I think only scales the text

  // Coult set the zoom:
  //m_browser->SetZoom(wxWEBVIEW_ZOOM_TINY);
  //m_browser->SetZoom(wxWEBVIEW_ZOOM_SMALL);
  //m_browser->SetZoom(wxWEBVIEW_ZOOM_MEDIUM);
  //m_browser->SetZoom(wxWEBVIEW_ZOOM_LARGE);
  //m_browser->SetZoom(wxWEBVIEW_ZOOM_LARGEST);
  //m_browser->SetZoomFactor( (float)120/100);


  //m_browser->SetPage
  //(
  //  "<html><title>New Page</title>"
  //  "<body>Created using <tt>SetPage()</tt> method.</body></html>",
  //  wxString()
  //);


  //m_edit_cut->Enable();
  //m_edit_copy->Enable();
  //m_edit_paste->Enable();

  //m_edit_undo->Enable();
  //m_edit_redo->Enable();

  //m_browser->CanCut();
  //m_browser->Cut();
  //m_browser->CanCopy();
  //m_browser->Copy();
  //m_browser->CanPaste();
  //m_browser->Paste();
  //m_browser->CanUndo();
  //m_browser->Undo();
  //m_browser->CanRedo();
  //m_browser->Redo();

  //if( m_browser->HasSelection() )
  //{
  //  m_browser->ClearSelection();
  //  m_browser->DeleteSelection();
  //  m_browser->SelectAll();
  //}


  // To Enter the JavaScript code to run as the initialization script that runs before any script in the HTML document.
  //if (!m_browser->AddUserScript( "window.my_var = 5;" ) )
  //  wxLogError("Could not add user script");


  // To run JS, and get back a stringified answer asynchronously, we can do like:
  //RunScript("function f(a){return a;}f(2.34);");
  //RunScript("function f(a){return a;}f(false);");
  //RunScript("function f(){var person = new Object();person.name = 'Foo'; person.lastName = 'Bar';return person;}f();");
  // For dates, we neet to take into account tize zone; below will return "2017-10-08T21:30:40.000Z"
  //RunScript("function f(){var d = new Date('10/08/2017 21:30:40'); \
        var tzoffset = d.getTimezoneOffset() * 60000; \
        return new Date(d.getTime() - tzoffset);}f();");
    //
    // To run script asynchonously
    //   `m_browser->RunScriptAsync("function f(a){return a;}f('Hello World!');");`
    // Will then cause `InterSpecWebFrame::OnScriptResult` to recieve the result, when its done
    //
    // Note:
    // To call C++ from JS, `window.wx.postMessage('This is a message from JS to C++ land');` will then call &InterSpecWebFrame::OnScriptMessage(..)
    // 

    // Connect the webview events
  Bind(wxEVT_WEBVIEW_NAVIGATING, &InterSpecWebFrame::OnNavigationRequest, this, m_browser->GetId());
  Bind(wxEVT_WEBVIEW_NAVIGATED, &InterSpecWebFrame::OnNavigationComplete, this, m_browser->GetId());
  Bind(wxEVT_WEBVIEW_LOADED, &InterSpecWebFrame::OnDocumentLoaded, this, m_browser->GetId());
  Bind(wxEVT_WEBVIEW_ERROR, &InterSpecWebFrame::OnError, this, m_browser->GetId());
  Bind(wxEVT_WEBVIEW_NEWWINDOW, &InterSpecWebFrame::OnNewWindow, this, m_browser->GetId());
  Bind(wxEVT_WEBVIEW_TITLE_CHANGED, &InterSpecWebFrame::OnTitleChanged, this, m_browser->GetId());
  //Bind(wxEVT_WEBVIEW_FULLSCREEN_CHANGED, &InterSpecWebFrame::OnFullScreenChanged, this, m_browser->GetId()); // Doesnt seem to ever get called, instead wxEVT_MAXIMIZE cused
  // OnScriptMessage will recieve messages from JS from code like:
  //  window.wx.postMessage('This is a message from JS to C++ land');
  Bind(wxEVT_WEBVIEW_SCRIPT_MESSAGE_RECEIVED, &InterSpecWebFrame::OnScriptMessage, this, m_browser->GetId());
  Bind(wxEVT_WEBVIEW_SCRIPT_RESULT, &InterSpecWebFrame::OnScriptResult, this, m_browser->GetId());

  //Connect the idle events
  Bind(wxEVT_IDLE, &InterSpecWebFrame::OnIdle, this);

  Bind(wxEVT_CLOSE_WINDOW, &InterSpecWebFrame::handleOnClose, this);
  Bind(wxEVT_MAXIMIZE, &InterSpecWebFrame::handleWindowMaximizeChange, this); // Only sent when maximized, not when restored
  Bind(wxEVT_SIZE, &InterSpecWebFrame::handleWinowSizeChange, this); //Listen to this for restore


  // TODO: The set/kill focus events dont seem to trigger - we can probably remove these callbacks anyway, as the JS handles this anyway
  //Bind(wxEVT_SET_FOCUS, &InterSpecWebFrame::handleOnFocus, this);
  //Bind(wxEVT_KILL_FOCUS, &InterSpecWebFrame::handleFocusLost, this);t
  // 
  // The EVT_CHILD_FOCUS does seem to reliably trigger
  Bind(wxEVT_CHILD_FOCUS, &InterSpecWebFrame::handleChildFocus, this);
  Bind(wxEVT_THREAD, &InterSpecWebFrame::handle_self_event, this);

  //
  Bind(wxEVT_MOTION, &InterSpecWebFrame::handle_mouse_move, this);
  Bind(wxEVT_LEFT_UP, &InterSpecWebFrame::handle_mouse_up, this);
 }


InterSpecWebFrame::~InterSpecWebFrame()
{
}


#ifdef _WIN32
WXLRESULT InterSpecWebFrame::MSWWindowProc(WXUINT nMsg, WXWPARAM wParam, WXLPARAM lParam)
{
  //Adapted from https://stackoverflow.com/questions/41179437/wxwidgets-wxborder-none-and-wxresize-border-makes-white-area

  switch (nMsg)
  {
  case WM_NCACTIVATE:
  {
    lParam = -1;
    break;
  }
  case WM_NCCALCSIZE:
    if (wParam)
    {
      HWND hWnd = (HWND)this->GetHandle();
      WINDOWPLACEMENT wPos;
      wPos.length = sizeof(wPos);
      GetWindowPlacement(hWnd, &wPos);
      if (wPos.showCmd != SW_SHOWMAXIMIZED)
      {
        RECT borderThickness;
        SetRectEmpty(&borderThickness);
        AdjustWindowRectEx(&borderThickness,
          GetWindowLongPtr(hWnd, GWL_STYLE) & ~WS_CAPTION, FALSE, NULL);
        borderThickness.left *= -1;
        borderThickness.top *= -1;
        NCCALCSIZE_PARAMS* sz = reinterpret_cast<NCCALCSIZE_PARAMS*>(lParam);
        sz->rgrc[0].top += 1;
        sz->rgrc[0].left += borderThickness.left;
        sz->rgrc[0].right -= borderThickness.right;
        sz->rgrc[0].bottom -= borderThickness.bottom;
        return 0;
      }
    }
    break;
  }
  return wxFrame::MSWWindowProc(nMsg, wParam, lParam);
}
#endif //#ifdef _WIN32


wxString InterSpecWebFrame::generate_token()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<unsigned int> dist;
  const unsigned int token = dist(rng);
  return std::to_string(token);
}//wxString InterSpecWebFrame::generate_token()


void InterSpecWebFrame::handleOnClose(wxCloseEvent& event)
{
  // We'll save the window position so the next time we open a window 
  //  it will be in same pos/size.
  const wxPoint win_pos = GetPosition();
  const wxSize win_size = GetSize();

  wxLogMessage("On window close: pos {%i,%i}, and size: {%i,%i}", 
    win_pos.x, win_pos.y, win_size.x, win_size.y );

  wxConfigBase* config = wxConfigBase::Get(true);
  config->Write("/WindowWidth", static_cast<long>(win_size.x));
  config->Write("/WindowHeight", static_cast<long>(win_size.y));
  config->Write("/WindowX", static_cast<long>(win_pos.x));
  config->Write("/WindowY", static_cast<long>(win_pos.y));
  config->Flush();

  auto app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance());
  if (app)
    app->handle_frame_closing(this);

  Destroy();
}//OnClose


/*
void InterSpecWebFrame::handleOnFocus(wxFocusEvent& event)
{
  wxLogMessage("Got Focus");
  event.Skip();

  RunScript("$('.app-titlebar').removeClass('inactive');");
}//void handleOnFocus(wxFocusEvent& event);

void InterSpecWebFrame::handleFocusLost(wxFocusEvent& event)
{
  // We dont seem to ever really get here, except when we programitacly call minimize, 
  // or whatever - so detecting window losing/getting focus in the JS must work okay enough
  wxLogMessage("Lost Focus");
  event.Skip();

  RunScript("$('.app-titlebar').addClass('inactive');");
}//void handleFocusLost(wxFocusEvent& event);
*/


void InterSpecWebFrame::handleChildFocus(wxChildFocusEvent& evt)
{
  wxLogMessage("handleChildFocus");

  // This next line doesnt appear to be necassary
  //RunScript("$('.app-titlebar').removeClass('inactive');");

  const auto app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance());
  assert(app);
  if (app)
    app->handle_frame_focus(this);
}//void handleChildFocus(wxChildFocusEvent& evt);


void InterSpecWebFrame::handleWindowMaximizeChange(wxMaximizeEvent& evt)
{
  evt.Skip();

  const bool isMax = IsMaximized();
  assert(isMax == true);

  m_currently_maximized = isMax;
  wxLogMessage("Screen maximized changed; max=%s", isMax ? "true" : "false");
  m_browser->RunScriptAsync("Wt.WT.TitleBarChangeMaximized(" + wxString(isMax ? "true" : "false") + ");");
}//handleWindowMaximizeChange(...)


void InterSpecWebFrame::handleWinowSizeChange(wxSizeEvent& evt)
{
  evt.Skip();

  const bool isMax = IsMaximized();
  if (isMax != m_currently_maximized)
  {
    m_currently_maximized = isMax;
    wxLogMessage("Screen sized changed; max=%s", isMax ? "true" : "false");
    m_browser->RunScriptAsync("Wt.WT.TitleBarChangeMaximized(" + wxString(isMax ? "true" : "false") + ");");
  }
}//handleWinowSizeChange(...)


/**
  * Method that retrieves the current state from the web control and updates the GUI
  * the reflect this current state.
  */
void InterSpecWebFrame::UpdateState()
{
  if (m_browser->IsBusy())
  {
    // We're loading
  }


}

void InterSpecWebFrame::OnIdle(wxIdleEvent& WXUNUSED(evt))
{
  /*x
  if (m_browser->IsBusy())
  {
    wxSetCursor(wxCURSOR_ARROWWAIT);
  }
  else
  {
    wxSetCursor(wxNullCursor);
  }
  */
}





/**
  * Callback invoked when there is a request to load a new page (for instance
  * when the user clicks a link)
  */
void InterSpecWebFrame::OnNavigationRequest(wxWebViewEvent& evt)
{
  // Here we would 
  wxLogMessage("%s", "Navigation request to '" + evt.GetURL() + "' (target='" +
    evt.GetTarget() + "')");

  //If we don't want to handle navigation then veto the event and navigation
  //will not take place, we also need to stop the loading animation
  bool handle_request = true; //... 
  if (handle_request)
  {
    // blah blah bl;ah
  }
  else
  {
    evt.Veto();
    // blah blah blah
  }
}

/**
  * Callback invoked when a navigation request was accepted
  */
void InterSpecWebFrame::OnNavigationComplete(wxWebViewEvent& evt)
{
  wxLogMessage("%s", "Navigation complete; url='" + evt.GetURL() + "'");
  UpdateState();
}

/**
  * Callback invoked when a page is finished loading
  */
void InterSpecWebFrame::OnDocumentLoaded(wxWebViewEvent& evt)
{
  //Only notify if the document is the main frame, not a subframe
  if (evt.GetURL() == m_browser->GetCurrentURL())
  {
    wxLogMessage("%s", "Document loaded; url='" + evt.GetURL() + "'");
  }
  UpdateState();
}


void InterSpecWebFrame::handle_self_event(wxThreadEvent& evt)
{
  wxLogMessage("In handle_self_event" );
  const wxString str_val = evt.GetString();

  if (str_val != "Wt::Http::Message")
  {
    wxLogError("InterSpecWebFrame::handle_self_event: Unrecognized event type: %s", str_val );
    return;
  }

  Wt::Http::Message* response_ptr = evt.GetPayload< Wt::Http::Message*>();
  assert(response_ptr);
  if (!response_ptr)
  {
    wxLogError("InterSpecWebFrame::handle_self_event: unexpected nullptr msg" );
    return;
  }

  std::unique_ptr<Wt::Http::Message> response_owner(response_ptr);
  const Wt::Http::Message& response = *response_ptr;

  std::string filename = "filename";
  const std::string* disposition_value = response.getHeader("Content-Disposition");
  if (disposition_value)
  {
    //attachment;filename="Ba133_gammas.csv";filename*=UTF-8''Ba133_gammas.csv

    std::vector<std::string> fields;
    SpecUtils::split(fields, *disposition_value, ";");
    for (std::string val : fields)
    {
      SpecUtils::trim(val);
      if (SpecUtils::istarts_with(val, "filename="))
      {
        filename = val.substr(9);
        SpecUtils::erase_any_character(filename, "'\\/:?\"<>|");
      }
    }//for (std::string val : fields)

  }//if(disposition_value )

  const std::string* type_value = response.getHeader("Content-Type");
  if (type_value)
  {
    // "text/csv"
  }

  wxString ext = SpecUtils::file_extension(filename);
  if (!ext.empty())
    ext = "*" + ext;

  wxConfigBase* config = wxConfigBase::Get(true);
  wxString defaultDir = config->Read("/LastSaveDir", wxString(""));
  if (!defaultDir.empty())
  {
    wxFileName defDirFile(defaultDir);
    if (!defDirFile.IsOk() || !defDirFile.DirExists() || !defDirFile.IsDirWritable())
      defaultDir = "";
  }

  wxFileDialog saveFileDialog(this, _("Save file"), defaultDir, filename,
    ext, wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
  if (saveFileDialog.ShowModal() == wxID_CANCEL)
  {
    wxLogMessage("User canceled saving file");
    return;     // the user changed idea...
  }

  // save the current contents in the file;
  // this can be done with e.g. wxWidgets output streams:
  const wxString savePath = saveFileDialog.GetPath();
  wxFileOutputStream output_stream(savePath);
  if (!output_stream.IsOk())
  {
    wxLogError("Cannot save current contents in file '%s'.", saveFileDialog.GetPath());
    return;
  }

  std::string body = response.body();
  output_stream.Write(body.c_str(), body.size());

  if (!output_stream.IsOk())
  {
    wxLogError("Error writing output stream.");
  }
  else
  {
    wxFileName lastSaveName(savePath);
    wxString lastSavePath = lastSaveName.GetPath();
    config->Write("/LastSaveDir", lastSavePath);
    wxLogMessage("Saved file: %s", saveFileDialog.GetPath().c_str());
  }
  output_stream.Close();
}//void handle_self_event(wxThreadEvent& evt)



void InterSpecWebFrame::handle_download_response(Wt::Http::Client* client, const boost::system::error_code &err, const Wt::Http::Message& response) {
  // TODO: should we be in a the MAIN thread here?  Cause I think we're in a Wt owned Asio thread, I would guess 

  //
  std::unique_ptr<Wt::Http::Client> client_ptr(client);

  if (!err) 
  {
    if (response.status() == 200) 
    {
      // Doing the Save As dialog here in this thread seems to work, but it isnt the proper thing 
      //  to do (we should be in wxWidgets main thread), so we'll post an event to be handled in
      //  the main thread.
      wxThreadEvent* evt = new wxThreadEvent();
      Wt::Http::Message* msg_copy = new Wt::Http::Message(response);
      evt->SetString("Wt::Http::Message");
      evt->SetPayload(msg_copy);
      this->QueueEvent(evt);
    }
    else
    {
      wxLogMessage("Http::Client invalis response code: %i", response.status());
    }//if (response.status() == 200)  / else
  }
  else {
    Wt::log("error") << "Http::Client error: " << err.message();
    wxLogMessage("Http::Client error: %s", err.message().c_str());
  }

}

/**
  * On new window, we veto to stop extra windows appearing
  */
void InterSpecWebFrame::OnNewWindow(wxWebViewEvent& evt)
{
  // I guess we would handle downloads and stuff here?
  
  const wxString url_str = evt.GetURL();
  wxURI uri(url_str);

  const wxString& server = uri.GetServer(); //"http://<server>/mypath"
  const wxString& path = uri.GetPath();//"http://mysite.com<path>"
  const wxString& port = uri.GetPort();//"http://mysite.com:<port>"
  const wxString& query = uri.GetQuery();//"http://mysite.com/mypath?<query>"
  const wxString& scheme = uri.GetScheme(); //e.x., http, "<scheme>://mysite.com"

  wxLogMessage("OnNewWindow type=%i, url=%s", int(evt.GetNavigationAction()), url_str.c_str());

  if (evt.GetNavigationAction() == wxWEBVIEW_NAV_ACTION_USER)
  {
    // This logic is roughly the same as for the macOS app, see AppDelegate.mm
    if ((path.find("request=redirect&url=http") != wxString::npos)
      || (scheme == "mailto")
      || (uri.HasServer() && (server != "127.0.0.1") && (server != "localhost"))
      )
    {
      //external url or email
      const bool launched = wxLaunchDefaultBrowser(url_str);
      if (launched)
        wxLogMessage("Launched '%s' in other application.", url_str.utf8_string().c_str());
      else
        wxLogMessage("Failed to launch '%s' in other application.", url_str.utf8_string().c_str());
    }
    else if ( uri.HasServer() && (server == "127.0.0.1" || server == "localhost"))
    {
      //CSV, spectrum file, JSON file, etc
      wxLogMessage("%s", "File to save; url='" + evt.GetURL() + "', target='" + evt.GetTarget() + "'");
      
      auto server = Wt::WServer::instance();
      assert(server);
      Wt::Http::Client* client = new Wt::Http::Client( server->ioService() );
      client->setTimeout(10);
      client->setMaximumResponseSize(512 * 1024 * 1024);

      client->done().connect( boost::bind(&InterSpecWebFrame::handle_download_response, this, client, boost::placeholders::_1, boost::placeholders::_2) );

      if (client->get(url_str.utf8_string()))
      {
        wxLogMessage("Have started GET for download");
      }
      else
      {
        wxLogMessage("Error calling Get for download URL");
      }

      /*
      wxFileDialog
        saveFileDialog(this, _("Save XYZ file"), "", "",
          "XYZ files (*.xyz)|*.xyz", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
      if (saveFileDialog.ShowModal() == wxID_CANCEL)
        return;     // the user changed idea...

      // save the current contents in the file;
      // this can be done with e.g. wxWidgets output streams:
      wxFileOutputStream output_stream(saveFileDialog.GetPath());
      if (!output_stream.IsOk())
      {
        wxLogError("Cannot save current contents in file '%s'.", saveFileDialog.GetPath());
        return;
      }
      */
      
      /*
      wxHTTP get;
      get.SetHeader(_T("Content-type"), _T("text/html; charset=utf-8"));
      get.SetTimeout(10); // 10 seconds of timeout instead of 10 minutes ...
      while (!get.Connect( m_url ))
      {
        wxLogMessage("Error connecting." );
        wxSleep(5);
      }
   -  wxInputStream *httpStream = get.GetInputStream(path);

      if (get.GetError() == wxPROTO_NOERR)
      {
        wxString res;
        wxStringOutputStream out_stream(&res);
        httpStream->Read(out_stream);

        wxString summary = "File size: " + std::to_string(res.size()) + ", headers:\n";
       
        wxMessageBox(res);
      }
      else
      {
        wxMessageBox(_T("Unable to connect!"));
      }

      wxDELETE(httpStream);
      get.Close();
      */
      

      /*
      wxURL url(url_str);
      if (url.GetError() == wxURL_NOERR)
      {
        wxString htmldata;
        wxInputStream* in = url.GetInputStream();

        if (in && in->IsOk())
        {
          wxStringOutputStream html_stream(&htmldata);
          in->Read(html_stream);
          wxLogMessage(htmldata);
        }
        delete in;
      }
      */
    }
    else if ((!uri.HasServer() && (url_str.find("data:application/octet-stream") != wxString::npos))
      || (!uri.HasServer() && (url_str.find("data:image/svg+xml") != wxString::npos)))
    {
      // Edge seems to automatically download these to the Downloads directory. 
      wxLogMessage("Got a data:image/svg+xml or data:application/octet-stream URL we're not handling atm" );
    }

  }
  else
  {
    wxLogMessage("%s", "Unknown OnNewWindow event: url='" + evt.GetURL() + "', target='" + evt.GetTarget() + "', NavAction: "+ std::to_string(evt.GetNavigationAction()) );
  }// etc

  


  UpdateState();
}

void InterSpecWebFrame::OnTitleChanged(wxWebViewEvent& evt)
{
  //SetTitle(evt.GetString());
  wxLogMessage("%s", "Title changed; title='" + evt.GetString() + "'");
}

/*
void InterSpecWebFrame::OnFullScreenChanged(wxWebViewEvent& evt)
{
  // Doesnt seem to get called for us
  wxLogMessage("Full screen changed; status = %d", evt.GetInt());
  ShowFullScreen(evt.GetInt() != 0);

  const bool isMax = IsMaximized();
  if (isMax)
  {
    RunScript("Wt.WT.TitleBarChangeMaximized(true);");
  }
  else
  {
    RunScript("Wt.WT.TitleBarChangeMaximized(false);");
  }
}
*/


void InterSpecWebFrame::OnScriptMessage(wxWebViewEvent& evt)
{
  wxLogMessage("Script message received; value = %s, handler = %s", evt.GetString(), evt.GetMessageHandler());
  wxString msg = evt.GetString();
  if (msg == "MinimizeWindow")
  {
    Iconize();
  }
  else if (msg == "ToggleMaximizeWindow")
  {
    const bool wasMax = IsMaximized();
    Maximize(!wasMax);
  }
  else if (msg == "CloseWindow")
  {
    wxLogMessage("Closing window");
    Close();
  } else if(msg == "SessionFinishedLoading")
  {
    wxConfigBase* config = wxConfigBase::Get(true);
    config->Write("/NumLoadAttempts", 0);


#ifdef _WIN32
    check_url_association();
#endif
  }
  else if (msg == "OpenInExternalBrowser")
  {
    const wxString token = generate_token();

    const int add_token_status = InterSpecServer::add_allowed_session_token(token.c_str(), InterSpecServer::SessionType::ExternalBrowserInstance);
    assert(add_token_status == 0);

    const wxString url_str = m_url + "?apptoken=" + token + "&restore=no";

    const bool launched = wxLaunchDefaultBrowser(url_str);
    if (launched)
      wxLogMessage("Launched '%s' in default browser.", url_str.utf8_string().c_str());
    else
      wxLogMessage("Failed to launch '%s' in default browser.", url_str.utf8_string().c_str());
  }
  else if (msg == "NewAppWindow")
  {
    wxLogMessage("Opening a new app window");
    auto app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance());
    assert(app);
    if (app)
      app->new_app_window();
  }
  else if (msg == "MouseDownInTitleBar")
  {
    wxLogMessage("MouseDownInTitleBar");
    CaptureMouse();
    m_dragging_window = true;
    wxPoint pos = wxGetMousePosition();
    wxPoint origin = GetPosition();
    int dx = pos.x - origin.x;
    int dy = pos.y - origin.y;
    m_mouse_down_pos = wxPoint(dx, dy);
  }
  else
  {
    wxLogMessage("Unrocgnized message from JS: '%s'.", msg.utf8_string().c_str());
  }
}//void InterSpecWebFrame::OnScriptMessage(wxWebViewEvent& evt)


void InterSpecWebFrame::handle_mouse_up(wxMouseEvent& evt)
{
  wxLogMessage("handle_mouse_up" );
 
  if (m_dragging_window && HasCapture() )
  {
    wxLogMessage("handle_mouse_up: Had Capture");
    ReleaseMouse();
  }
  else
  {
    evt.Skip();
  }

  m_dragging_window = false;
}//handle_mouse_up

void InterSpecWebFrame::handle_mouse_move(wxMouseEvent& evt)
{
  if (!m_dragging_window)
  {
    evt.Skip();
    return;
  }

  if (!evt.LeftIsDown())
  {
    m_dragging_window = false;
    evt.Skip();
    return;
  }

  wxPoint pt = evt.GetPosition();
  if (evt.Dragging())
  {
    wxPoint pos = ClientToScreen(pt);
    Move(wxPoint(pos.x - m_mouse_down_pos.x, pos.y - m_mouse_down_pos.y));
  }
}//void handle_mouse_move(wxMouseEvent& evt)


void InterSpecWebFrame::OnScriptResult(wxWebViewEvent& evt)
{
  if (evt.IsError())
    wxLogError("Async script execution failed: %s", evt.GetString());
  else
    wxLogMessage("Async script result received; value = %s", evt.GetString());
}



void InterSpecWebFrame::RunScript(const wxString& javascript)
{
  wxLogMessage("Running JavaScript:\n%s\n", javascript);

  wxString result;
  if (m_browser->RunScript(javascript, &result))
  {
    wxLogMessage("RunScript() returned \"%s\"", result);
  }
  else
  {
    wxLogWarning("RunScript() failed");
  }
}


/**
  * Callback invoked when a loading error occurs
  */
void InterSpecWebFrame::OnError(wxWebViewEvent& evt)
{
#define WX_ERROR_CASE(type) \
    case type: \
        category = #type; \
        break;

  wxString category;
  switch (evt.GetInt())
  {
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_CONNECTION);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_CERTIFICATE);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_AUTH);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_SECURITY);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_NOT_FOUND);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_REQUEST);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_USER_CANCELLED);
    WX_ERROR_CASE(wxWEBVIEW_NAV_ERR_OTHER);
  }

  wxLogMessage("%s", "Error; url='" + evt.GetURL() + "', error='" + category + " (" + evt.GetString() + ")'");

  //Show the info bar with an error
  //m_info->ShowMessage(_("An error occurred loading ") + evt.GetURL() + "\n" +
  //  "'" + category + "'", wxICON_ERROR);

  UpdateState();
}


const wxString& InterSpecWebFrame::app_token() const
{
  return m_token;
}