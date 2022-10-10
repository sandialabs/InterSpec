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
class wxInfoBar;
class wxWebView;
class wxIdleEvent;
class wxCloseEvent;
class wxWebViewEvent;



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
  void OnFullScreenChanged(wxWebViewEvent& evt);
  void OnScriptMessage(wxWebViewEvent& evt);
  void OnScriptResult(wxWebViewEvent& evt);
  void OnError(wxWebViewEvent& evt);

  void RunScript(const wxString& javascript);

  void OnClose(wxCloseEvent& event);
private:
  wxWebView* m_browser;

  wxInfoBar* m_info;

  wxString m_url;
  wxString m_token;
};//class InterSpecWebFrame




#endif //InterSpec_WEB_FRAME_h