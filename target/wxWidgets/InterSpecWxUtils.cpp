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

#include <wx/file.h>
#include <wx/log.h>
#include <wx/filedlg.h>
#include <wx/dirdlg.h>
#include <wx/msgdlg.h>
#include <wx/config.h>
#include <wx/filename.h>

#include <Wt/WServer>
#include <Wt/WApplication>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/DirectorySelector.h"

#include "InterSpecWxApp.h"
#include "InterSpecWxUtils.h"

namespace InterSpecWxUtils
{
  void handle_javascript_error( const std::string &error_msg, const std::string app_token )
  {
    InterSpecWxApp::handle_javascript_error( error_msg, app_token );
  }


  void save_file_data( std::string data, std::string suggested_name )
  {
    auto app = dynamic_cast<InterSpecWxApp *>( wxApp::GetInstance() );
    if( !app )
    {
      wxLogError( "InterSpecWxUtils::save_file_data: failed to get wxApp" );
      return;
    }

    wxWindow *topWindow = app->GetTopWindow();
    if( !topWindow )
    {
      wxLogError( "InterSpecWxUtils::save_file_data: failed to get top window" );
      return;
    }

    // Dispatch to the main UI thread for the file dialog
    topWindow->GetEventHandler()->CallAfter(
      [data = std::move( data ), suggested_name = std::move( suggested_name )]()
    {
      const std::string ext = SpecUtils::file_extension( suggested_name );

      wxString wildcard = "All files (*.*)|*.*";
      if( ext == ".png" )
        wildcard = "PNG images (*.png)|*.png|" + wildcard;
      else if( ext == ".svg" )
        wildcard = "SVG images (*.svg)|*.svg|" + wildcard;

      wxConfigBase *config = wxConfigBase::Get( true );
      wxString defaultDir = config->Read( "/LastSaveDir", wxString( "" ) );
      if( !defaultDir.empty() )
      {
        wxFileName defDirFile( defaultDir );
        if( !defDirFile.IsOk() || !defDirFile.DirExists() || !defDirFile.IsDirWritable() )
          defaultDir = "";
      }

      wxFileDialog dlg( nullptr, _("Save file"), defaultDir,
                        wxString::FromUTF8( suggested_name ),
                        wildcard,
                        wxFD_SAVE | wxFD_OVERWRITE_PROMPT );

      if( dlg.ShowModal() != wxID_OK )
        return;

      const wxString savePath = dlg.GetPath();
      wxFile file( savePath, wxFile::write );
      if( !file.IsOpened() )
      {
        wxLogError( "Cannot save file '%s'.", savePath );
        return;
      }

      file.Write( data.data(), data.size() );
      file.Close();

      wxFileName lastSaveName( savePath );
      const wxString lastSavePath = lastSaveName.GetPath();
      config->Write( "/LastSaveDir", lastSavePath );
      wxLogMessage( "Saved file: %s", savePath );
    } );
  }//void save_file_data(...)


  void browse_for_directory( const std::string &title,
                             const std::string &message,
                             std::function<void(const std::vector<std::string> &)> callback )
  {
    Wt::WApplication * const app = Wt::WApplication::instance();
    if( !app || !callback )
      return;

    const std::string session_id = app->sessionId();

    wxApp * const wxapp = dynamic_cast<wxApp *>( wxApp::GetInstance() );
    if( !wxapp )
    {
      wxLogError( "InterSpecWxUtils::browse_for_directory: failed to get wxApp" );
      return;
    }

    wxWindow * const topWindow = wxapp->GetTopWindow();
    if( !topWindow )
    {
      wxLogError( "InterSpecWxUtils::browse_for_directory: failed to get top window" );
      return;
    }

    // Dispatch to the wx main UI thread for the dialog; post the result back to the Wt session.
    topWindow->GetEventHandler()->CallAfter(
      [title, message, callback, session_id, topWindow]()
    {
      wxConfigBase * const config = wxConfigBase::Get( true );
      wxString defaultDir = config->Read( "/LastSaveDir", wxString( "" ) );
      if( !defaultDir.empty() )
      {
        wxFileName defDirFile = wxFileName::DirName( defaultDir );
        if( !defDirFile.IsOk() || !defDirFile.DirExists() )
          defaultDir = "";
      }

      const wxString prompt = message.empty() ? wxString::FromUTF8( title )
                                              : wxString::FromUTF8( message );

      wxDirDialog dlg( topWindow, prompt, defaultDir,
                       wxDD_DEFAULT_STYLE | wxDD_DIR_MUST_EXIST );
      if( !title.empty() )
        dlg.SetTitle( wxString::FromUTF8( title ) );

      std::vector<std::string> paths;
      if( dlg.ShowModal() == wxID_OK )
      {
        const wxString chosen = dlg.GetPath();
        paths.push_back( std::string( chosen.utf8_str() ) );
        config->Write( "/LastSaveDir", chosen );
      }

      // Post the result back to the Wt session thread.
      Wt::WServer * const server = Wt::WServer::instance();
      if( !server )
        return;

      server->post( session_id, [callback, paths](){
        Wt::WApplication * const app = Wt::WApplication::instance();
        if( !app )
          return;
        callback( paths );
        app->triggerUpdate();
      } );
    } );
  }//void browse_for_directory(...)


  void register_native_directory_picker()
  {
    set_wx_native_directory_picker( &browse_for_directory );
  }
}//namespace InterSpecWxUtils
