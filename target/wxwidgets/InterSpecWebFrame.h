#ifndef InterSpec_WEB_FRAME_h
#define InterSpec_WEB_FRAME_h
/* SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.

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


#include <wx/frame.h>




class wxString;
class wxWebView;
class wxSizeEvent;
class wxIdleEvent;
class wxCloseEvent;
class wxFocusEvent;
class wxMouseEvent;
class wxWebViewEvent;
class wxMaximizeEvent;
class wxChildFocusEvent;


class InterSpecWebFrame : public wxFrame
{
public:
  InterSpecWebFrame(const wxString& url, const bool no_restore, const wxString &file_to_open );
  virtual ~InterSpecWebFrame();

  void UpdateState();
  void OnIdle(wxIdleEvent& evt);
  void OnNavigationRequest(wxWebViewEvent& evt);
  void OnNavigationComplete(wxWebViewEvent& evt);
  void OnDocumentLoaded(wxWebViewEvent& evt);
  void OnNewWindow(wxWebViewEvent& evt);
  void OnTitleChanged(wxWebViewEvent& evt);
  //void OnFullScreenChanged(wxWebViewEvent& evt);
  void OnScriptMessage(wxWebViewEvent& evt);
  void OnScriptResult(wxWebViewEvent& evt);
  void OnError(wxWebViewEvent& evt);

  void RunScript(const wxString& javascript);

  void handleOnClose(wxCloseEvent& evt);
  //void handleOnFocus(wxFocusEvent& evt);
  //void handleFocusLost(wxFocusEvent& evt);
  void handleChildFocus(wxChildFocusEvent& evt);
  void handleWindowMaximizeChange(wxMaximizeEvent& evt);
  void handleWinowSizeChange(wxSizeEvent& evt);

  const wxString& app_token() const;


  /* Stop hack to use HTML titlebar to move window around when #m_dragging_window is true. */
  void handle_mouse_up(wxMouseEvent& evt);

private:
  /** Dispatch a click on a `target="_blank"` link: download via the DLL for
   localhost URLs, hand off to the OS default app for mailto/external URLs.

   Called from both #OnNavigationRequest (which is needed on macOS WKWebView,
   where `wxEVT_WEBVIEW_NEWWINDOW` does not reliably fire for `target="_blank"`
   link activations) and #OnNewWindow (which is the path Windows Edge takes).
   #m_last_target_blank_url dedupes the two paths when both fire. */
  void handle_target_blank_url( const wxString &url_str );

public:

  /* Hack to use HTML titlebar to move window around, when #m_dragging_window is true. */
  void handle_mouse_move(wxMouseEvent& evt);


#ifdef _WIN32
  /* We need to overide MSWWindowProc to make a frameless window, otherwise there will be a few pixel
  white-area at the top of the window that looks terrible.
  */
  WXLRESULT MSWWindowProc(WXUINT nMsg, WXWPARAM wParam, WXLPARAM lParam) wxOVERRIDE;
#endif


protected:
  static wxString generate_token();

private:
  wxWebView* m_browser;

  wxString m_url;
  wxString m_token;

  /** The Edge WebView2 doesnt respect the "-webkit-app-region: drag" property, meaning we would need 
  the "native" header bar to allow dragging window around, taking up an extra 30 pixels or so.  So 
  intead we will detect when the user clicks on the titlebar, in JS, and then capture the mouse in
  c++, and move the window around until the user lets go.
  */
  bool m_dragging_window;
  wxPoint m_mouse_down_pos;

  bool m_currently_maximized;

  /** Last URL handled by #handle_target_blank_url, used to dedupe between
   wxEVT_WEBVIEW_NAVIGATING and wxEVT_WEBVIEW_NEWWINDOW (Windows Edge fires
   both for the same link click; macOS WKWebView only fires NAVIGATING).
   Cleared via CallAfter so the same link can be re-clicked. */
  wxString m_last_target_blank_url;
};//class InterSpecWebFrame




#endif //InterSpec_WEB_FRAME_h