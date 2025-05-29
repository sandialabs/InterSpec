#ifndef MacOsUtils_h
#define MacOsUtils_h
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

#include <string>

namespace macOsUtils
{
  /** When the macOS app starts, it checks the preference "DoNotResume" to see if
   it should block loading the apps previous state.  Once the app loads client
   side, InterSpecApp calls back to c++ from JS and sets "DoNotResume" to NO,
   which is what this function does.
   */
  void sessionSuccessfullyLoaded();
  
  std::string static_data_base_dir();
  std::string user_data_dir();
  
  bool openFinderToPath( const std::string &filepath );
  bool showFileInFinder( const std::string &filepath );
  
  /** Asynchronously, but modally, shows a native file/directory chooser, using obj-c code, calling the provided callback when
   selection is done.
   
   This function must be called from within the main Wt thread of your application, and the callback with results will be posted
   to this same application instance.
   
   @param title The utf-8 title of the selection dialog
   @param message The utf-8  message for the selection dialog
   @param canChooseFiles
   @param canChooseDirectories
   @param allowsMultipleSelection
   @param callback The callback called when user selection is done; the argument will have the filesystem paths
          of the selected files/directories, or if the user canceled, this vector will be empty.  This callback is made from
          within the main Wt thread of your the instance that called this function.  The WApplication::triggerUpdate() function
          will be called after calling your function.
   
   Upon failure, some messages is printed to stdout, but thats about it.
   */
  void showFilePicker( const std::string title, const std::string message,
                           const bool canChooseFiles,
                           const bool canChooseDirectories,
                           const bool allowsMultipleSelection,
                           const std::function<void(const std::vector<std::string> &)> callback );
}

#endif
