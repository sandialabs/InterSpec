#ifndef WX_INTERSPEC_UTILS_h
#define WX_INTERSPEC_UTILS_h
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
#include <functional>

// This header / src could be eliminated now- but will leave in till we're sure we dont need something like this

namespace InterSpecWxUtils
{
  void handle_javascript_error( const std::string &error_msg, const std::string app_token );

  /** Saves binary file data to a user-selected path using a native wxFileDialog.
   Called from a Wt worker thread; dispatches to the wxWidgets main thread.
  */
  void save_file_data( std::string data, std::string suggested_name );

  /** Opens a native wxDirDialog and invokes `callback` on the Wt session
      thread with the selected path (single-element vector) or an empty
      vector if the user cancelled.
      Must be called from a Wt event-loop thread.
  */
  void browse_for_directory( const std::string &title,
                             const std::string &message,
                             std::function<void(const std::vector<std::string> &)> callback );

  /** Registers browse_for_directory as the wx native directory picker for
      LibInterSpec's DirectorySelector widget.  Must be called once, from the
      wx executable, after the wxApp instance has been constructed (e.g. from
      InterSpecWxApp::OnInit).  Cannot be done from inside LibInterSpec
      because wxWidgets is statically linked separately into the dylib and
      the executable, so wx static state does not survive that boundary.
  */
  void register_native_directory_picker();
}//namespace InterSpecWxUtils


#endif //WX_INTERSPEC_UTILS_h