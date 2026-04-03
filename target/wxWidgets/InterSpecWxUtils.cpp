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
#include <wx/msgdlg.h>
#include <wx/config.h>
#include <wx/filename.h>

#include "SpecUtils/StringAlgo.h"

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
}//namespace InterSpecWxUtils
