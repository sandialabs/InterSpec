#ifndef InterSpecApp_h
#define InterSpecApp_h
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

#include <set>
#include <mutex>
#include <chrono>
#include <string>

#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#if(  BUILD_AS_WX_WIDGETS_APP )
#include <functional>
#endif

class InterSpec;
namespace Wt
{
  class WText;
  class WGridLayout;
  class WPushButton;
}//namespace Wt

class PopupDivMenuItem;

#if( INCLUDE_ANALYSIS_TEST_SUITE )
class SpectrumViewerTester;
#endif

//Feature under development (20171203) to prompt user if they want to load their
//  previous state, before actually doing it.
// Could make this always true for anything but mobile devices...
// Still deciding if I like this or not... loading the state doesnt seem to be
//   a huge slowdown in the startup.
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
#define PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE 0
#else
#define PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE 0
#endif

class InterSpec_API InterSpecApp : public Wt::WApplication
{
public:
  InterSpecApp( const Wt::WEnvironment& env );
  virtual ~InterSpecApp();

  //isMobile(): Tests to see if a tablet or a phone platform
  bool isMobile() const;
  
  //isPhone(): Tests to see if on a iPhone, Android phone, or Blackberry device
  bool isPhone() const;
  
  //isTablet(): Tests to see if on a iPad or Android tablet
  bool isTablet() const;
  
  //isAndroid(): Searches for Android in the user agent string
  bool isAndroid() const;
  
#if( WT_VERSION>=0x3030400 )
  //ToDo: define a custom handler for javascript errors
//  virtual void handleJavaScriptError(const std::string& errorText)
//  {
//    std::cerr << "JS error: " << errorText << std::endl;
//  }
#endif
  
  /* svlog send the message to the InterSpec for display to the user. */
  void svlog( const Wt::WString& message, int priority = 1 );
  
  /* convenience overload */
  void svlog( const char *message, int priority = 1 );
  
  /* convenience overload */
  void svlog( const std::string& message, int priority = 1 );

  InterSpec *viewer();

  std::chrono::steady_clock::time_point::duration activeTimeInCurrentSession() const;
  
  /** Returns the amount of use-time since `InterSpecUser::totalTimInApp()` has been updated. */
  std::chrono::steady_clock::time_point::duration timeSinceTotalUseTimeUpdated() const;
  
  
  //userNameFromOS(): Caution, will return 'apache' if being served, from
  //  an apache server, 'mobile' if on a iOS device, or blank upon failure.
  static std::string userNameFromOS();

  //These functions functions should be moved to a differenct header/source file
  //  perhaps named ServerUtilityFunctions.h/.cpp
  std::string getUserNameFromEnvironment() const;  //returns CGI env variable "REMOTE_USER"
  
  //tempDirectory(): when the app is deployed as a FCGI application, the temp
  //  path may not be avaliable as a system environment variable, or if it is,
  //  it might not be the same as where Apache is configured to use.
  //  The answer returned by this function comes from first checking the CGI
  //  environemnt parameters, and then falling back to
  //  SpecUtils::temp_dir() if CGI doesnt specify temp directory
  //  locations.
  //  The actual answer returned is actually the result cached from the first
  //  call to this function where WApplication::instance() was non-null; the
  //  InterSpecApp constructor calls this function so the proper directory
  //  can be resolved at any future point.  Calls before a InterSpecApp
  //  is created may return just the result of SpecUtils::temp_dir().
  //This function is necassarry as essentually a workaround for FCGI deployment.
  //Note that currently this function requires locking a mutex, as well as some
  //  other not-super-fast actions.
  static std::string tempDirectory();

#if(  BUILD_AS_WX_WIDGETS_APP )
  static void setJavascriptErrorHandler( std::function<void(std::string, std::string)> fctn );
#endif
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  /** Returns the token passed as part of url using parameter 'externalid' or 'apptoken'.
   e.g., if url looked like "localhost:8080?externalid=blah", then this function would return "blah"
   
   Note that if you want to assign a session token but not have this be the primary window, add a URL argument of 'primary=no'.
  */
  std::string externalToken();
  
  /** Returns true if this app instance corresponds to the primary window instance. */
  static bool isPrimaryWindowInstance();
  
  /** Get a pointer to the app-instance, using the external token to identify the instance.
   
   Returns nullptr on failure.
   
   Note that you still need to take an WApplication::UpdateLock before accessing the app instance,
   and you may need to call WApplication::triggerUpdate() to cause changes to propagate to user.
   */
  static InterSpecApp *instanceFromExtenalToken( const std::string &externalToken );
  

  bool userOpenFromFileSystem( const std::string &path );
  
  static std::set<InterSpecApp *> runningInstances();
  
  //Need put put in a function that will minimize memmorry usage of all cuurently
  //  running sessions.
#endif

#if( defined(WIN32) && BUILD_AS_ELECTRON_APP )
  //When users drag files from Outlook on windows into the app
  //  you can call the following functions
  void dragEventWithFileContentsStarted();
  void dragEventWithFileContentsFinished();
#endif
  

#if( BUILD_AS_OSX_APP || IOS || BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
  static void osThemeChange( std::string name );
#endif
  
  Wt::JSignal<unsigned int> &hotkeySignal();
  
#if( IOS )
  /** Function called from objective-c to update the safe areas defined in iOS
      (e.g., the area for the notch).
      The orientation int cooresponds to UIDeviceOrientation,
      (e.g., UIDeviceOrientationLandscapeLeft==3, UIDeviceOrientationLandscapeRight==4)
      and the rest of the arguments is user pixels (i.e., not physical pixels).
   
      Initial safe are insets can be set as a command line argument "&SafeAreas=3,0,44,22,44"
   */
  void setSafeAreaInsets( const int orientation, const float top,
                          const float right, const float bottom,
                          const float left );
  
  enum class DeviceOrientation
  {
    Unknown = 0,
    Portrait = 1,
    PortraitUpsideDown = 2,
    LandscapeLeft = 3,   //Notch to the left on iPhoneX
    LandscapeRight = 4,  //Notch to the right on iPhoneX
    FaceUp = 5,
    FaceDown = 6
  };
    
  /** Get the safe area insets. */
  void getSafeAreaInsets( DeviceOrientation &orientation, float &top, float &right, float &bottom, float &left );
#endif
  
  //clearSession(): resets the session to a blank session (e.g. doesnt load
  //  any saved state).
  virtual void clearSession();
    
  void miscSignalHandler( const std::string &signal );
  
  /** Handles "deep" urls.
   
   Meant to handle URLs with the scheme "interspec://...", that would be passed to the application
   by the OS, like when a QR code is scanned.
   
   The prefix "interspec://" is optional, and may be omitted.
   
   Currently just passes through to InterSpec::handleAppUrl(...)
   
   Returns if the URL was successfully used.
   
   As of 20220405: currently being implemented and not fully fleshed out; in the future could
                   implement to handle spectra, display decay data for a isotope, DRFs, etc.  Could
                   also have an argument to the regular browsers URL that this function could
                   process.
   */
  bool handleAppUrl( const std::string &url );
  
  /** Loads the XML file for current locale, to use for localizing strings
   
   \param name The base-name of the XML file to load.  For this, you will usually call with an argument
          like "InterSpec", or "PeakEdit".  You do not need to give the path to the XML file, or the
          extension; this function will take care of adding those.
   
   If the XML file has already been loaded, no action will be taken.
   */
  void useMessageResourceBundle( const std::string &name );
  
  /** Looks in `InterSpec_resources/app_text/` for files like `InterSpec_(.+).xml`, and returns
   the list of the capture group.
   */
  static const std::set<std::string> &languagesAvailable();
  
protected:

  //notify(): over-riding WApplication::notify in order to catch any exceptions
  //  that may happen during event handling
  virtual void notify( const Wt::WEvent &event );

  //finalize(): called before destruction to take care of things that might
  //  involve calls to virtual functions.  We use it to record usage statistics
  virtual void finalize();
  
  /** In principle, called when user navigates away from the page or reloads; however, this appears to happen
   for Chrome and Firefox, but not Safari.
   */
  virtual void unload();
  
  //prepareForEndOfSession(): increments users total time in app and saves the
  //  state if their preferences ask for it
  virtual void prepareForEndOfSession();
  
  /** Adds `m_activeTimeSinceDbUpdate` to the database, saving the updated result. */
  void updateUsageTimeToDb();
  
#if(  BUILD_AS_WX_WIDGETS_APP )
  virtual void handleJavaScriptError( const std::string &errorText );
#endif

#if( !BUILD_FOR_WEB_DEPLOYMENT )
  /** Checks for URL argument "externalid", or equivalently "apptoken", and if found sets m_externalToken to it.
   
   Returns true if session should continue to be loaded; also notifies InterSpecServer we have loaded this token.
   */
  bool checkExternalTokenFromUrl();
#endif
  
  //setupDomEnvironment(): loads required JS and CSS resource, sets the
  //  appropriate meta tags.
  void setupDomEnvironment();
  
  //setupWidgets(): initializes the InterSpec widget, as well as the header
  //  and footer if applicable; if app has already been initialized, the
  //  existing widgets will be deleted and then re-created.
  //  If 'attemptStateLoad' is specified, the app will attempt to load the most
  //  recent saved state, unless the URL contains a restore=0 argument.
  void setupWidgets( const bool attemptStateLoad );
  
protected:
  InterSpec *m_viewer;
  Wt::WGridLayout *m_layout;

#if( !BUILD_FOR_WEB_DEPLOYMENT )
  Wt::Signal<const InterSpecApp *> m_destructing;
#endif
  
  std::chrono::steady_clock::time_point m_lastAccessTime;
  std::chrono::steady_clock::time_point::duration m_activeTimeInSession;
  /** We will occasionally update the use duration in the database - this field tracks how long of
   active use since the database was last updated.
   */
  std::chrono::steady_clock::time_point::duration m_activeTimeSinceDbUpdate;
  
#define OPTIMISTICALLY_SAVE_USER_STATE 0
  //If OPTIMISTICALLY_SAVE_USER_STATE is enabled, then the users state will
  //  attempt to be saved whenever a 'onbeforeunload' is recieved.  The downside
  //  is that right now for chrome this doesnt work for a simple refresh, and
  //  it may be different for different browsers.
  //For the case when the server doesnt get the signal, and the user has
  //  navigated away, we could decrease the timeout settings in wt_config.xml
  //  so the session gets destroyed quicker after refreshes
  
#if( OPTIMISTICALLY_SAVE_USER_STATE )
  //m_thinkingLeaveSignal: an experimental signal that tries to catch unload
  //  events a little better; still doesnt catch refresh events.
  std::unique_ptr<Wt::JSignal<> > m_thinkingLeaveSignal;
#endif
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  /** App token specified using the URL "externalid" or "apptoken" arguments. */
  std::string m_externalToken;
  
  /** Specifies wether this instance corresponds to the primary application window (e.g., uses the native menus, receives request to
   open files, etc).
   Defaults too true, if a 'apptoken'/'externalid' is specified in the url, but can be set to false by including a url argument of 'primary=no'.
   
   \TODO: should default this to false, even if there is an external session, and instead have the
   primary application (electron, obj-c, java, etc) instead add a 'primary=true' argument to the URL, and require the URL to have been
   allowed through #InterSpecServer::add_allowed_session_token before this can be set to true.
   */
  bool m_primaryApp;
#endif
  
  //iOS version of app handles not reloading app state through the iOS restore
  //  stuff.
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || ANDROID || BUILD_AS_WX_WIDGETS_APP )
  void loadSuccesfullCallback();
  std::unique_ptr<Wt::JSignal<> > m_sucessfullyLoadedSignal;
#endif
  
  /** A signal emitted from the JS with a string argument that specifies what caused the signal to be emitted.
   So far this signal is only used to allow manually specifying a button to be clicked on the notification messages.
   
   \sa #miscSignalHandler
   */
  std::unique_ptr<Wt::JSignal<std::string> > m_miscSignal;
  
  /** Signal triggered when user presses control + another button  */
  Wt::JSignal<unsigned int> m_hotkeySignal;
  
  
#if( IOS )
  /* Note that we could setup a JS based 'orientationchange' JSignal to let us know
     when there is an orientation change from the JS, however it still isnt clear how
     to get the 'safe-area-inset-left' type CSS values in JavaScript (not that
     it actually matters though I dont thing) - so for now will just rely on the objective-c
     notification (but maybe JS would be better to cover android to... - for a latter date!)
   */
  DeviceOrientation m_orientation;
  
  /** CSS px values for areas on iPhoneX's like the notch or bottom bar. */
  float m_safeAreas[4];
#endif
  
  friend class InterSpec;
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  friend class SpectrumViewerTester;
#endif
};//class InterSpecApp

#endif  //InterSpecApp_h


