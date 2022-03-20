#ifndef ElectronUtils_h
#define ElectronUtils_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#if( BUILD_AS_ELECTRON_APP )

/* Since we are using node-addon-api to interface to LibInterSpec, we could use
 standard c++ (namespaced) functions, however to avoid any issues with compiler
 versions or whatever, we will only call extern "C" functions from N-API.
 
 Once you add the function you want into this header/src file, you need to add the
 wrapper function for it into InterSpecAddOn.cpp.
 */

#ifdef _WIN32
#define LIB_INTERFACE(type) __declspec(dllexport) type __cdecl
#else
#define LIB_INTERFACE(type) __attribute__ ((visibility ("default"))) type
#endif

extern "C"
{
  /** Returns port being served on.
   Will always serve on 127.0.0.1.
   
   All input paths should be UTF-8.
   \param basedir has only been tested to be a directory relative to CWD, i.e, "."
   \param xml_config_path has only been tested to be a directory relative to CWD, i.e, "./data/config/wt_config_electron.xml"
   
   ToDo: Make sure basedir and xml_config_path all work okay on Windows - maybe consider making them relative (in the JS probably)
   */
  LIB_INTERFACE(int) interspec_start_server( const char *process_name,
                                            const char *user_data_dir,
                                            const char *basedir,
                                            const char *xml_config_path );
  
  /** ToDo: also remove all session tokens
   */
  LIB_INTERFACE(void) interspec_kill_server();
  
  /** ToDo: allow adding of mutliple tokens; should also add a mechanism to
   mark the token as having been loaded, so it cant be reused.
   */
  LIB_INTERFACE(void) interspec_add_allowed_session_token( const char *session_token );
  
  /** Returns 0 if authorized but not alive, 1 if session is alive, -1 if dead, -2 if nor authorized. */
  LIB_INTERFACE(int) interspec_remove_allowed_session_token( const char *session_token );
  
  /** Returns 1 if session with specified token is alive, otherwise 0 */
  LIB_INTERFACE(int) interspec_session_is_alive( const char *session_token );
  
  /** Open one or more files from the filesystem.  For macOS this would be
   dropping a file to the app icon in the Finder.  Or on Windows it would be
   dropping multiple files onto app icon.
   
   ToDo: should add logic to try and figure out what files are background, foreground, etc.
   
   \param session_token session token you want to open the file.
   \param files_json: a JSON array of files to open.  Opened,
   ex. ["/path/to/some/file","/some/other/file"] .
   Input files are opened one at a time using #InterSpecApp::userOpenFromFileSystem
   \returns number of files opened.
   Will be zero if files were not valid spectrum files.
   -1 if files_json is invalid format.
   -2 if session_token is invalid (however, currently will ask ALL
   sessions to open the files... ToDo: decide on this behaviour).
   */
  LIB_INTERFACE(int) interspec_open_file( const char *session_token, const char *files_json );

  /** Returns if using the Electron-native menus.  If you are using the HTML-electron menu, will return false.
   */
  LIB_INTERFACE(bool) interspec_using_electron_menus();
}//extern "C"

#endif //#if( BUILD_AS_ELECTRON_APP )


namespace ElectronUtils
{
#if( USING_ELECTRON_NATIVE_MENU )
  /** Requests main.js to load a new clean session (i.e., don restore any state)
   This is a workaround to when the user requests a new session, the normal
   mechanism in c++ creates duplicate Electron menu items...
   
   Must be called from within a WApplication thread (e.g., wApp is valid).
   
   @returns whether message was succesfully sent or not.
   */
  bool requestNewCleanSession();
#endif
  
  /** Notify parent application (main.js, or objective-c) that the session has
   loaded.
   
   Must be called from within a WApplication thread (e.g., wApp is valid).
   
   @returns whether message was succesfully sent or not.
   
   Note: main.js will wait till recieveing this notification before asking the
   session to open any files the OS requested.
   */
  bool notifyNodeJsOfNewSessionLoad();

  /**
   
   Returns true if it thinks message was sent.
   */
  bool send_nodejs_message( const std::string &msg_name, const std::string &msg_data );
}//namespace ElectronUtils


#endif  //#ifndef ElectronUtils_h
