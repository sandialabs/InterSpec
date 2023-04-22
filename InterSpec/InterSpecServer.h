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
  /** Used by iOS target to start the web-server.
   
   Only real difference from #startWebServer is that it allows specifying two arguments like:
     "--tempdir /some/.../path/", that set the "TMPDIR" and "WT_TMP_DIR" environment variables
     to the specified path for spooling uploads too.
   
   #TODO: Check if we can have the iOS target just use #startWebServer and set the environment variables itself (or using a dedicated function here)
   */
  void startServer( int argc, char *argv[],
                                Wt::WApplication::ApplicationCreator createApplication );
  
  /** Starts the web-server, with specified options; primarily called from #start_server.
   
   Does not set any of the InterSpec configuration variables (#InterSpecServer::start_server does
   all that), but does set the server-related static variables (i.e., port and url being served on).
   
   Also starts initializing DecayDataBaseServer as soon as the thread pool is available, so main
   GUI thread wont block waiting for it to initialize (hopefully).
   */
  void startWebServer( std::string proccessname,
                             std::string basedir,
                             const std::string configpath,
                             unsigned short int server_port_num = 0
#if( BUILD_FOR_WEB_DEPLOYMENT )
                             , std::string http_address = "127.0.0.1"
#endif
                              );
  
  /** Starts the server, given the various options.
   Sets user data database directory; creates user data directory if needed and sets writable
   directory; looks for serial_to_model in user data directory; upgrades user data database if
   needed; and if the static data directory hasnt already been explicitly set it will be set to
   `basedir + "/data"`.
   Then calls #startWebServer.

   @param process_name Name of the executable being ran - not sure if/why we actually need this 
          TODO: check if we can just remove this
   @param userdatadir The directory to store user data to.  
          On Windows this might be "C:\Users\<user>\AppData\Roaming\InterSpec"
   @param basedir The directory that contains the `InterSpec_resources` and `data` directories.
   @param xml_config_path The path for the Wt XML config file
   @param server_port_num The port to start the server on.  
          If zero, then a random port above 1024 will be chosen.  
          Ports below 1024 usually require admin privledges.
#if( BUILD_FOR_WEB_DEPLOYMENT )
   @param http_address The network adapter address to bind the web-server to.  "127.0.0.1" is
          localhost, while "0.0.0.0" will bind so it is visible on the external network
#endif
   @returns port being served on, or a negative value on error.

   This function is currently used to start the server for all versions of the app besides iOS.
   */
  int start_server( const char *process_name, const char *userdatadir,
                    const char *basedir, const char *xml_config_path,
                    unsigned short int server_port_num = 0
#if( BUILD_FOR_WEB_DEPLOYMENT )
                    , const char *http_address = "127.0.0.1"
#endif
                   );
  

  void killServer();
  
  /** Blocks current thread until the Wt server is stopped by some other means. */
  int wait_for_shutdown();
  
  //portBeingServedOn(): will only be valid if this instance of the app is
  //  serving the webpages.  Will be -1 if not serving.
  int portBeingServedOn();
  
  //urlBeingServedOn(): will only be valid if this instance of the app is
  //  serving the webpages.  Will be empty if not serving.
  //  Example value returned: "http://127.0.0.1:7234"
  std::string urlBeingServedOn();
  
  
  bool changeToBaseDir( int argc, char *argv[] );
  
  std::string getWtConfigXml( int argc, char *argv[] );
  
  enum class SessionType
  {
    PrimaryAppInstance, ExternalBrowserInstance
  };

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
  int add_allowed_session_token( const char *session_id, const SessionType type );

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

  /** Returns if the session was found, and if so the type set when #add_allowed_session_token
   was called.
   */
  std::pair<bool,SessionType> session_type( const char *session_token );

  /** Sets the session as dead. */
  void set_session_destructing( const char *session_token );
  

  /** Open one or more files from the filesystem.  For macOS this would be
   dropping a file to the app icon in the Finder.  Or on Windows it would be
   dropping multiple files onto app icon.
   
   URLs with the "interspec://" schema can also be included instead of a file path; 
   they will be processed after the files are opened.

   ToDo: should add logic to try and figure out what files are background, foreground, etc.
   
   \param session_token session token you want to open the file.
   \param files_json: a JSON array of files to open.  Opened,
   ex. ["/path/to/some/file","/some/other/file","interspec://drf/define?..."] .
   Input files are opened one at a time using #InterSpecApp::userOpenFromFileSystem
   \returns number of files opened.
   Will be zero if files were not valid spectrum files.
   -1 if files_json is invalid format.
   -2 if session_token is invalid (however, for everywhere besides iOS, currently asks ALL
   sessions to open the files... ToDo: decide on this behavior).

   TODO: currently if session_token is valid, file will be opened synchronously, meaning you 
         can delete the file immediately after this call.  But if token is invalid, it will
         open it asynchonously (in all open sessions), meaning you cant delete the file yet... 
         Not a very nice interface, should brobably remove this latter behaviour
   */
  int open_file_in_session( const char *session_token, const char *files_json );

  /** Passes the specified url to InterSpecApp::handleAppUrl(url), and returns
  if the URL was used.
  */
  bool pass_app_url_to_session( const char *session_token, const std::string &url );

  /** Sets a file to open upon session load.
   
   The \c session_token must have been set by #add_allowed_session_token, but must not have actually
   been loaded yet.
   
   The file path must either be a file-system path to a spectrum file, or it can be a App URL with
   the schema "interspec://..."
   
   Call #clear_file_to_open_on_load after you have succesfully handled the file, or else the "Clear Session..."
   action will cause the session to load this same file again.
   
   Throws exception if an invalid session token, or the session has already been loaded.
   */
  void set_file_to_open_on_load( const char *session_token, const std::string file_path );

  /** Returns the file path, or app url, of any files to be opened during creating of a new session.
   
   Returns empty string if none.
   */
  std::string file_to_open_on_load( const std::string &session_token );
  
  /** Once you have loaded the file/url from #file_to_open_on_load, call this function to avoid
   loading the spectrum again if "Clear Session..." action is taken by the user.
   */
  void clear_file_to_open_on_load( const std::string &session_token );


#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
  struct DesktopAppConfig
  {
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
    /**
      @brief The proxy settings to use; only relevant when using the maps tool, 
      and you are behind a proxy, and it is a Windows or Linux build. Defaults 
      to empty string.
     
       Valid values are: empty (default), 'direct', 'auto_detect', 'system'; any other 
       string will be interpreted as the 'proxyRules' setting documented at 
       https://www.electronjs.org/docs/latest/api/session#sessetproxyconfig 
       (a value might be 'http=proxy.somecompany.com:80,direct://')
     
       This setting is only applicable if you use the map tool (the only part of 
       InterSpec to use the internet); if map tiles wont load, and you are behind a proxy, 
       you _may_ need to set this setting; usually a value of 'system' will work.  
       If you set this settting, and you are not behind a proxy, InterSpec startup 
       time may hang for ~30 seconds, usually with a white screen, while a proxy 
       request times out.
    */
    std::string m_proxy;
#endif

    /**
     @brief The HTTP port to serve the application on; if zero (recomended!), 
     will choose random port on startup.  Defaults to zero (i.e., random high
     number port)/
     
     If non-zero, should be a port larger than 1024 (e.x., 8080) since ports 
     below this may require admin privledges (never run InterSpec as admin!).
    */
    unsigned short int m_http_port;

    /** 
    @brief Wether to require browser sessions to have a one-use-only token in 
    thier URL.  Defaults to true.

    Normally external sessions (i.e. 'View' -> 'Use in external browser') get 
    assigned a one-time-use token that is required to load InterSpec into the 
    browser.  Without a valid token, you cant load a session in the browser.  
    If you allow external sessions without tokens, then the token wont be 
    needed - and any application that can access your localhost network can 
    create a session and potentually access your data.  
    
    It is not recomended to to enable this setting.
    */
    bool m_require_token;

    /** Allow restoring previous session, if applicable. */
    bool m_allow_restore;

    /** If true, will open the WebView inspector console on launch. */
    bool m_open_dev_tools;

#if( USE_LEAFLET_MAP )
    /** The key to use to access the https://arcgis.com server to get map tiles from.
     This value overides the default one built into InterSpec, if specified.
     */
    std::string m_arcgis_key;
#endif
    
    /** Returns application config, as determined from defaults and then 
    InterSpec_app_settings.json files.

    First looks in "`app_data_dir`/config/InterSpec_app_settings.json"
    and modifies default values for settings found, then looks in 
    "`user_data_dir`/InterSpec_app_settings.json" and further modifies
    values.  
    If a file doesnt exist, it is ignored.
    IF a directory passed in is empty, it is skipped.

    If there is an error in the file (invalid JSON or invalid data type)
    then an exception is thrown.
    */
    static DesktopAppConfig init( const std::string &app_data_dir, 
                                  const std::string &user_data_dir );

  private:
    DesktopAppConfig();
  };//struct DesktopAppConfig
#endif

}//namespace InterSpecServer


#endif /* InterSpecServer_hpp */
