#ifndef InterSpecServer_hpp
#define InterSpecServer_hpp
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
#include <stdio.h>

#include "InterSpec_config.h"

#include <Wt/WApplication>


namespace InterSpecServer
{
  void startServer( int argc, char *argv[],
                                Wt::WApplication::ApplicationCreator createApplication );
  
  void startServerNodeAddon( std::string proccessname,
                             std::string basedir,
                             const std::string configpath );
  
  
  void killServer();
  
  
  
  //portBeingServedOn(): will only be valid if this instance of the app is
  //  serving the webpages.  Will be -1 if not serving.
  int portBeingServedOn();
  
  //urlBeingServedOn(): will only be valid if this instance of the app is
  //  serving the webpages.  Will be empty if not serving.
  //  Example value returned: "http://127.0.0.1:7234"
  std::string urlBeingServedOn();
  
  
  bool changeToBaseDir( int argc, char *argv[] );
  
  std::string getWtConfigXml( int argc, char *argv[] );
  

  /** Add a session token to be tracked or allowed.
   
   Session tokens are given as part of the url using the key 'apptoken' or 'externalid', and are used to interact with a specific session
   being displayed in a particular browser tab, or WebView, or whatever, programmatically through the c++.
   An example of adding a token to a url would be: "http://localhost:8080?apptoken=23324assd1&restore=no&..."
      
   By default any token, or no token at all are required, but you can change this by calling #set_require_tokened_sessions, but if you do,
   then you you also need to call #add_allowed_session_token with a particular token _before_ a session with that token in the URL
   can be loaded.  If you dont restrict requiring tokens, then the status of sessions with tokens that you have not called
   #add_allowed_session_token for will not be tracked, but the InterSpecApp class will still keep these tokens so that you can interact
   with individual sessions.
   
   @Returns:
   - 0 if hadnt been seen before
   - 1 if authorized but not seen yet
   - 2 if authorized and loaded
   -3 if deauthorized session
   - 4 if dead session
   
   if session had been seen before (e.g., non-zero return value), no changes will be made to that sessions allowed state.
   */
  int add_allowed_session_token( const char *session_id );

  /** Returns -1 if invalid token.  Returns +1 if valid token that had never been loaded.  Returns zero if was loaded.  */
  int remove_allowed_session_token( const char *session_token );

  /** Returns if we should allow sessions without, or without valid, tokens.
   Defaults to true unless you call #set_require_tokened_sessions.
   */
  bool allow_untokened_sessions();

  /** Set wether sessions without tokens, or without valid tokens are allowed.
      Default is all sessions allowed.
   */
  void set_require_tokened_sessions( const bool require );

  /** Returns the status for the specified session.
    Returns:
      - 0 if #add_allowed_session_token has not been called for the token
      - 1 if #add_allowed_session_token has been called, but session not yet loaded
      - 2 if session has been loaded
      - 3 If session is no longer authorized (e.g., #remove_allowed_session_token called for it)
      - 4 if session is dead
   */
  int session_status( const char *session_token );
  
  /** Function that should be called when a session is initially loaded.
   
   Returns 0 if had been authorized, 1 if had been seen, and 2 if dead or no longer authorized
   */
  int set_session_loaded( const char *session_token );

  /** Sets the session as dead. */
  void set_session_destructing( const char *session_token );
  

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
  int open_file_in_session( const char *session_token, const char *files_json );
}//namespace InterSpecServer


#endif /* InterSpecServer_hpp */
