#ifndef INTERSPEC_WX_APP_h
#define INTERSPEC_WX_APP_h
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

#include <string>
#include <vector>

#include <wx/app.h>


class wxSingleInstanceChecker;


class IpcServer;
class InterSpecWebFrame;


#if !wxUSE_CMDLINE_PARSER
#error "ex command line parsing must be enabled"
#endif // !wxUSE_CMDLINE_PARSER

class InterSpecWxApp : public wxApp
{
public:
  InterSpecWxApp();

  virtual bool OnInit() override;
  virtual int OnExit() override;

  virtual void OnInitCmdLine(wxCmdLineParser& parser) override;
  virtual bool OnCmdLineParsed(wxCmdLineParser& parser) override;

  void handle_frame_closing(InterSpecWebFrame* frame);
  void handle_open_file_message(const std::string& message);

private:
  wxString m_url;
  std::vector<wxString> m_command_line_args;

  std::vector<InterSpecWebFrame*> m_frames;
  InterSpecWebFrame* m_active_frame;

  wxSingleInstanceChecker* m_checker;

  IpcServer* m_ipc_server;
};


#endif //INTERSPEC_WX_APP_h