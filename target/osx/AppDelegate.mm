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

#import "AppDelegate.h"
#import "Availability.h"
#import <AppKit/NSView.h>
#import <AppKit/NSSavePanel.h>

#include <set>
#include <string>
#include <thread>
#include <chrono>
#include <stdlib.h>

//We gotta fix some wierd errors...
#ifdef check
#undef check
#endif
#ifdef require
#undef require
#endif

#include <Wt/WString>
#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Serializer>

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

#if( USE_SPECRUM_FILE_QUERY_WIDGET )
#include "InterSpec/SpecFileQueryWidget.h"
#endif

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SerialToDetectorModel.h"

@implementation OurWebView
- (BOOL)performDragOperation:(id <NSDraggingInfo>)sender
{
  // This function lets files go through to WkWebView - we actually dont need this function,
  //  but leaving in for the moment as I experiment
  
  // We wont handle promises here, so this way WkWebView will keep handling them,
  //  But here is some code to intercept them, in case it turns out we should.
  //  (TODO: double check Mail and Outlook attachments can be dragged in)
  //NSPasteboard *pboard = [sender draggingPasteboard];
  //NSDictionary *options = @{};
  //NSArray<Class> *promises_class = @[[NSFilePromiseReceiver class]];
  //NSArray<NSURL*> *promises = [pboard readObjectsForClasses:promises_class options:options];
  //
  //for( NSFilePromiseReceiver *promise in promises )
  //{
  //  [promise receivePromisedFilesAtDestination: [NSURL fileURLWithPath: NSTemporaryDirectory()]
  //                                     options: @{}
  //                              operationQueue: [NSOperationQueue new]
  //                                      reader:^(NSURL * _Nonnull fileURL, NSError * _Nullable errorOrNil) {
  //    if (errorOrNil) {
  //      NSLog(@"Error: %@", errorOrNil);
  //      return;
  //    }
  //    NSLog(@"fileURL: %s", [fileURL fileSystemRepresentation]);
  //  }];
  //}//for( NSFilePromiseReceiver *promise in promises )
  
  return [super performDragOperation: sender];
}//performDragOperation - e.g., when file is dropped

- (NSDragOperation)draggingEntered:(id <NSDraggingInfo>)sender {
  
  NSLog( @"In draggingEntered!" );
  
  NSDragOperation sourceDragMask = [sender draggingSourceOperationMask];
  NSPasteboard *pboard = [sender draggingPasteboard];
  
  NSArray<NSPasteboardType> * types = [pboard types];
  if( !types )
    return [super draggingEntered: sender];
  
  
  NSArray<Class> *url_class = @[[NSURL class]];
  NSDictionary *options = @{};
  NSArray<NSURL*> *files = [pboard readObjectsForClasses:url_class options:options];
  
  // If we use obj-c, to create an array, and then print that to JSON, like is commented out,
  //  we cant use the NSJSONWritingWithoutEscapingSlashes option of NSJSONSerialization
  //  (requires macOS 10.15), and so we'll get slashes '/' escaped, which JSON.parse(), doesnt
  //  like, and we cant use as-is either.
  //  So we'll just bring in a little c++/Wt for the moment, so strings will be properly escaped.
  //NSMutableArray *array = [NSMutableArray arrayWithCapacity:[files count]];
  Wt::Json::Array json_array;
  
  // Form some JSON to set a variable in JS client-side, that can then make the call to server code
  for (NSURL *url in files)
  {
    const char *fs_path = [url fileSystemRepresentation];  //Will be UTF-8 encoded
    json_array.push_back( Wt::WString::fromUTF8(fs_path) );
    //NSString *fs_path_str = [[NSString alloc] initWithCString:fs_path encoding:NSUTF8StringEncoding];
    //[array addObject:fs_path_str];
    //NSLog( @"draggingEntered: git URL: %s", fs_path );
  }
  
  
  // We wont do anything about promises (like dragging an attachment from mail or Outlook),
  //  at the moment (since if we do, then WkWebView wont handle them) - so for these cases
  //  we will just let the upload handle as normal
  //NSArray<Class> *promises_class = @[[NSFilePromiseReceiver class]];
  //NSArray<NSURL*> *promises = [pboard readObjectsForClasses:promises_class options:options];
  //for( NSFilePromiseReceiver *promise in promises )
  //{
  //}
  
  
  if( json_array.empty() )
  {
    //Probably got file promise(s) here - we'll just fall-back to normal file upload, and
    //  let WkWebView take care of getting the files.
    NSString *js = @"(function(){$(document).data('dragOverFilePaths',null);})()";
    [self evaluateJavaScript: js completionHandler:nil];
  }else
  {
    //NSError *error = nil;
    //NSData *jsonData = [NSJSONSerialization dataWithJSONObject:array options:NSJSONWritingPrettyPrinted error:&error];
    //NSString *json_string = [[NSString alloc] initWithData:jsonData encoding:NSUTF8StringEncoding];
    
    const std::string json_data = Wt::Json::serialize(json_array);
    const std::string js_str =
    "(function(){\n\t"
    "let fns = {};\n\t"
    "fns.time = new Date();\n\t"
    "fns.filenames = " + json_data + ";\n\t"
    "$(document).data('dragOverFilePaths',fns);\n\t"
    "console.log(\"Set file paths\",fns);\n"
    "})();";
    
    NSString *js = [[NSString alloc] initWithCString:js_str.c_str() encoding:NSUTF8StringEncoding];
    
    //NSLog( @"draggingEntered: json_string: %@", js );
    
    [self evaluateJavaScript: js completionHandler:nil];
  }//if( !json_array.empty() )
  
  // Leaving in the below for development only.
  if ( [types containsObject:NSPasteboardTypeFileURL] ) {
    NSLog( @"draggingEntered: Is NSPasteboardTypeFileURL:" );
  }else if ( [types containsObject:NSPasteboardTypeURL] ) {
    NSLog( @"draggingEntered: Is NSPasteboardTypeURL:" );
  }else if ( [types containsObject:NSFileContentsPboardType] ) {
    NSLog( @"draggingEntered: Is NSFileContentsPboardType:" );
  }
  
  // We'll let the WKWebView handle the drag-n-drop as normal
  return [super draggingEntered: sender];
}
@end



@implementation AppDelegate

@synthesize persistentStoreCoordinator = _persistentStoreCoordinator;
@synthesize managedObjectModel = _managedObjectModel;
@synthesize managedObjectContext = _managedObjectContext;


+ (void)initialize
{
  
}


//See https://www.bignerdranch.com/blog/cocoa-ui-preservation-yall/ for how to actually use.
//See also : https://developer.apple.com/library/archive/documentation/General/Conceptual/MOSXAppProgrammingGuide/CoreAppDesign/CoreAppDesign.html

- (void)application:(NSApplication *)app didDecodeRestorableState:(NSCoder *)coder
{
  // This seems to get called only when the app was killed by like Xcode or force-kill of the app
  NSLog( @"didDecodeRestorableState" );
}


- (void)application:(NSApplication *)app willEncodeRestorableState:(NSCoder *)coder
{
  // Gets called when you change window size, or put app into background, etc.
  // (not when you quit the app normally)
  NSLog( @"willEncodeRestorableState" );
}


- (void)application:(NSApplication *)application openURLs:(NSArray<NSURL *> *)urls
{
  NSLog( @"In openURLs!" );
    
  // Note: Implementing this method causes `openFile` to not be called, so this function will
  //  be called anytime the user double-clicks on a spectrum file in the Finder, as well
  //  as whenever a user clicks on a URL of the scheme "interspec://..." in the operating system.


  // From the terminal, if we type: `open 'interspec://drf:v1:MyName:"My other Par":Eqn=323.32+5.2l(x)-3.5E-6*ln(y);someOtherPar'`
  //
  // Then here we get an NSArray with one entry, with url properties:
  //  [url scheme] == interspec
  //  [url absoluteString]: interspec://drf:v1:MyName:%22My%20other%20Par%22:Eqn=323.32+5.2l(x)-3.5E-6*ln(y);someOtherPar
  //  [url absoluteURL]: interspec://drf:v1:MyName:%22My%20other%20Par%22:Eqn=323.32+5.2l(x)-3.5E-6*ln(y);someOtherPar
  //  [url resourceSpecifier]: //drf:v1:MyName:%22My%20other%20Par%22:Eqn=323.32+5.2l(x)-3.5E-6*ln(y);someOtherPar
  //  [url standardizedURL]: interspec://drf
  //  [[url absoluteString] stringByRemovingPercentEncoding]: interspec://drf:v1:MyName:"My other Par":Eqn=323.32+5.2l(x)-3.5E-6*ln(y);someOtherPar
  
  std::unique_ptr<Wt::WApplication::UpdateLock> app_lock;
  InterSpecApp *specapp = nullptr;
  
  if( _UrlUniqueId )
    specapp = InterSpecApp::instanceFromExtenalToken( [_UrlUniqueId UTF8String] );
  
  if( specapp )
  {
    app_lock.reset( new Wt::WApplication::UpdateLock( specapp ) );
    if( !(*app_lock) )
    {
      NSLog( @"\topenURLs: _UrlUniqueId unable to get WApplication::UpdateLock for apptoken='%@' - not handling URL!", _UrlUniqueId );
      
      app_lock.reset();
    }
  }else
  {
    NSLog( @"\topenURLs: _UrlUniqueId unable to get app for apptoken='%@' - not handling URL! ", _UrlUniqueId );
  }
    
  for( NSURL *url in urls )
  {
    // NSLog( @"scheme: %@", [url scheme] );
    if( !url )
    {
      NSLog( @"openURLs: null url" );
      continue;
    }
    
    std::string urlcontent;
    
    if( [url isFileURL] )
    {
      urlcontent = [url fileSystemRepresentation];
      
      NSLog( @"Will open specfile '%s'", urlcontent.c_str() );
      
      if( app_lock )
      {
        specapp->userOpenFromFileSystem( urlcontent );
        specapp->triggerUpdate();
      }
    }else
    {
      NSString *absStr = [url absoluteString];
      if( !absStr )
      {
        NSLog( @"openURLs: null absoluteString" );
        continue;
      }
      
      urlcontent = [absStr UTF8String];
      
      if( app_lock )
      {
        specapp->handleAppUrl( urlcontent );
        specapp->triggerUpdate();
      }
    }// if( [url isFileURL] ) / else
    
    
    if( !app_lock && !urlcontent.empty() )
    {
      //The WebView may not have requested the URL yet when this function gets
      //  called, in the case a user double clicks on a file in the Finder to launch
      //  this app.  Therefore we will mark that we should open this file using the
      //  URL argument (that uses an intermediate database to store file location).
      NSLog( @"Will mark the file to open once a session is created" );
      
      // I dont quite understand why we need this next '[NSString alloc] initWithString', but we
      //  do or the app will crash when opened wit ha long app-url, but not with a file...
      //  This probably indicates something is wrong with my understanding of the obj-c the memory
      //  management somehow, but I guess for now I'll cross my fingers, or put my head in the sand.
      _fileNeedsOpening = [[NSString alloc] initWithString: [NSString stringWithUTF8String:urlcontent.c_str()]];
      
      NSLog( @"Initially fileNeedsOpening %@", _fileNeedsOpening );
    }//if( !app_lock )
    
  }//for( NSURL *url in urls )
}//openURLs


-(BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)sender
{
  return YES;
}

- (void) terminated: (NSNotification *)notification
{
  //Issue: if we remove the ovbserver, then the messages from the initial
  //       launch of the app no longer come through, if we dont, then a ton
  //       of cpu will eventually be used up
  [[NSNotificationCenter defaultCenter] removeObserver:self];
  
  [NSApp terminate:self];
}

Wt::WApplication *createApplication(const Wt::WEnvironment& env)
{
  return new InterSpecApp( env );
}

- (void)applicationWillFinishLaunching:(NSNotification *)aNotification
{
  NSLog(@"applicationWillFinishLaunching");
  _fileNeedsOpening = nil;
  _isServing = NO;
  _UrlServingOn = @"";
  _UrlUniqueId = nil;
}


- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
  NSLog(@"In applicationDidFinishLaunching");
  
  /*
   //Could maybe get rid of using XIB/NIB by manually creating a window like:
   NSRect frame = NSMakeRect(0, 0, 300, 300);
   NSWindow *window  = [[[NSWindow alloc] initWithContentRect:frame
                 styleMask:NSBorderlessWindowMask
                 backing:NSBackingStoreBuffered
                 defer:NO] autorelease];
   [window setBackgroundColor:[NSColor redColor]];
   [window makeKeyAndOrderFront:NSApp];
   [window setFrameUsingName: @"MainWindow"];
   [window setFrameAutosaveName: @"MainWindow"];
   */
  
  //ToDo: I think the OS provides a better mechanism than manually tracking if
  //      shutdown was clean or not... but first tries with didDe/willEn-codeRestorableState
  //      didnt work super well.
  NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
  //Will return false by defautl
  BOOL doResume = [defaults boolForKey:@"DoResume"];
  NSLog( @"DoResume %i", doResume );
  
  //Set to not resume; when we get confirmation all loaded okay, we will set
  // it back to resuming by default.
  [defaults setBool:NO forKey:@"DoResume"];
  
  /*
   //Check for file "DoNotResume" in data directory.  Probably not necassary since force-killing the app
   //  should achieve the same things.  But leaving commented out until I test stuff.
   NSFileManager *fileManager = [NSFileManager defaultManager];
   NSURL *testResumeFile = [[self applicationFilesDirectory] URLByAppendingPathComponent:@"DoNotResume"];
   NSString *testResumeFilePath = [testResumeFile path];
   if ([fileManager fileExistsAtPath:testResumeFilePath]) {
   NSLog( @"Found DoNotResume file at %@\n\tWill not resume previous state.", testResumeFilePath );
   doResume = NO;
   NSError *error = nil;
   [[NSFileManager defaultManager] removeItemAtPath:testResumeFilePath error:&error];
   }
   */
  
  
  if( doResume )
  {
    //ToDo: Make sure window is showing... or check if this is necassary even
  }else
  {
    //Place window at a reasonable size and location
    const int width = [[NSScreen mainScreen] frame].size.width;
    const int height = [[NSScreen mainScreen] frame].size.height;
    NSSize mySize;
    mySize.width = 0.85*width;
    mySize.height = 0.85*height;
    [_window setContentSize: mySize];
    [[self window] setFrameTopLeftPoint:NSMakePoint(0.025*width, 0.95*height)];
  }//if( doResume ) / else
  
  
  
  WKWebViewConfiguration *webConfig = [[WKWebViewConfiguration alloc] init];
  //Setting the config like bellow seems to slow down the rendering.
  
  if ([webConfig respondsToSelector:@selector(applicationNameForUserAgent)]) {
    webConfig.applicationNameForUserAgent = @"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_4) AppleWebKit/603.1.30 (KHTML, like Gecko) Version/10.1 Safari/603.1.30";
  }
  
  //webConfig.ignoresViewportScaleLimits = ;  //Not sure what this is for.
  if ([webConfig respondsToSelector:@selector(suppressesIncrementalRendering)]) {
    webConfig.suppressesIncrementalRendering = YES;
  }
  
  WKPreferences *prefs = [webConfig preferences];
  //prefs.javaEnabled = NO;    //default is false anyway, and deprecated in 10.15
  //prefs.plugInsEnabled = NO; //default is false anyway, and deprecated in 10.15
  //prefs.javaScriptEnabled = YES; //default is true anyway, and deprecated in 10.15
  
  if ([webConfig respondsToSelector:@selector(javaScriptCanOpenWindowsAutomatically)]) {
    //Does seem to get here for some reason
    prefs.javaScriptCanOpenWindowsAutomatically = YES;  //default is true in macOS anyway
  }
  
  //prefs.minimumFontSize = 6.0;
  if ([webConfig respondsToSelector:@selector(tabFocusesLinks)]) {
    prefs.tabFocusesLinks = NO;
  }
  
  // Make it so the obj-c didReceiveScriptMessage method will be called any time the javascript
  //  does something like:
  //      window.webkit.messageHandlers.interOp.postMessage({"action": "DoSomething"});
  [webConfig.userContentController addScriptMessageHandler: self name:@"interOp"];
  
  //Some additional settings we may want to set:
  //  see more at http://jonathanblog2000.blogspot.com/2016/11/understanding-ios-wkwebview.html
  //[prefs setValue:@YES forKey:@"allowFileAccessFromFileURLs"];
  //[prefs setValue:@YES forKey:@"acceleratedDrawingEnabled"];
  //[prefs setValue:@YES forKey:@"displayListDrawingEnabled"];
  //[prefs setValue:@YES forKey:@"visualViewportEnabled"];
  //[prefs setValue:@NO forKey:@"standalone"];
  //[prefs setValue:@YES forKey:@"fullScreenEnabled"];
  //[prefs setValue:@YES forKey:@"developerExtrasEnabled"];
  //[prefs setValue:@YES forKey:@"resourceUsageOverlayVisible"];
  //[prefs setValue:@NO forKey:@"peerConnectionEnabled"]; //WebRTC
  //[prefs setValue:@NO forKey:@"mediaDevicesEnabled"]; //cameras, microphones, etc
  //[prefs setValue:@NO forKey:@"screenCaptureEnabled"];
  //[prefs setValue:@YES forKey:@"javaScriptCanAccessClipboard"];  //doesnt seem to make a difference for flux tool (always HTML text)
  //hiddenPageDOMTimerThrottlingEnabled
  //hiddenPageDOMTimerThrottlingAutoIncreases
  //pageVisibilityBasedProcessSuppressionEnabled
  

  const std::string base_dir = [[[NSBundle mainBundle] resourcePath] UTF8String];
  const std::string data_dir = [[[[NSBundle mainBundle] resourcePath] stringByAppendingPathComponent: @"data"] UTF8String];
  const std::string user_data_dir = [[[self applicationFilesDirectory] path] UTF8String];
  
  bool require_token = true;
#if( PERFORM_DEVELOPER_CHECKS )
  bool allow_dev_tools = true;
#else
  bool allow_dev_tools = false;
#endif
  unsigned short int server_port = 0;
  try
  {
    const InterSpecServer::DesktopAppConfig app_config
                   = InterSpecServer::DesktopAppConfig::init( data_dir, user_data_dir );
    
    if( !app_config.m_allow_restore )
      doResume = NO;
    allow_dev_tools = (allow_dev_tools || app_config.m_open_dev_tools);
    server_port = app_config.m_http_port;
    require_token = app_config.m_require_token;
    UndoRedoManager::setMaxUndoRedoSteps( app_config.m_max_undo_steps );
  }catch( std::exception &e )
  {
    NSString *error = [NSString stringWithUTF8String:e.what()];
    
    NSAlert *alert = [[NSAlert alloc] init];
    [alert setMessageText:@"Error parsing application JSON options!"];
    [alert setInformativeText: error ];
    [alert addButtonWithTitle:@"Ok"];
    [alert setAlertStyle:NSAlertStyleCritical];
    [alert runModal];
  }//try / catch get confiugurations

  InterSpecServer::set_require_tokened_sessions( require_token );
  
  NSString *tempDir = NSTemporaryDirectory();
  if( tempDir == nil ) //shouldnt ever fail, right
    tempDir = @"/tmp";
  static const std::string tmpdr = [tempDir UTF8String];  //static since I'm not sure how long the location pointed to by setenv has to last
  
  setenv( "TMPDIR", tmpdr.c_str(), 1);
  setenv( "WT_TMP_DIR", tmpdr.c_str(), 1);
  
  
  if( allow_dev_tools )
  {
    [prefs setValue:@YES forKey:@"developerExtrasEnabled"];
    
    //Note: currently I disable right click in InterSpecApp using javascript if
    //      no PERFORM_DEVELOPER_CHECKS.  However, could also do:
    //https://stackoverflow.com/questions/28801032/how-can-the-context-menu-in-wkwebview-on-the-mac-be-modified-or-overridden#28981319
    //[_InterSpecWebView willOpenMenu:<#(nonnull NSMenu *)#> withEvent:<#(nonnull NSEvent *)#>];
  }//if( allow_dev_tools )
  
  
  //Create WKWebView manually, rather than in XIB to support macOS 10.10 and 10.11...
  self.InterSpecWebView = [[OurWebView alloc] initWithFrame: _window.contentView.frame configuration: webConfig];
  [_window.contentView addSubview:self.InterSpecWebView];
  
  //Make sure the WKWebView resizes.
  self.InterSpecWebView.autoresizingMask = (NSViewWidthSizable | NSViewHeightSizable);
  [_window.contentView setAutoresizesSubviews: YES];
  
  //Set the user agent string so
  NSString *userAgentStr = @"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_4) AppleWebKit/603.1.30 (KHTML, like Gecko) Version/10.1 Safari/603.1.30";
  if ([_InterSpecWebView respondsToSelector:@selector(setCustomUserAgent)]) {
    [_InterSpecWebView setCustomUserAgent: userAgentStr];
    _InterSpecWebView.allowsLinkPreview = NO;
  }
  
  _InterSpecWebView.allowsBackForwardNavigationGestures = NO;  //default is NO anyway
  
  
  if( allow_dev_tools )
  {
    [prefs setValue:@YES forKey:@"developerExtrasEnabled"];
    
    NSMenu *menu = [[[NSApp mainMenu] itemAtIndex: 1] submenu]; //Get "Edit" menu
    
    if( menu )
    {
      NSMenuItem *itemnow = [[NSMenuItem alloc]
                             initWithTitle:@"Enable Web Inspector"
                             action:@selector(enableWebInspector)
                             keyEquivalent:@""];
      [itemnow setTarget:self];
      [menu addItem:itemnow];
      [itemnow setEnabled:YES];
    }
    
    //Note: currently the right click in InterSpecApp is disabled using javascript if
    //      unles PERFORM_DEVELOPER_CHECKS is true.  However, could also do:
    //https://stackoverflow.com/questions/28801032/how-can-the-context-menu-in-wkwebview-on-the-mac-be-modified-or-overridden#28981319
    //[_InterSpecWebView willOpenMenu:<#(nonnull NSMenu *)#> withEvent:<#(nonnull NSEvent *)#>];
  }//if( allow_dev_tools )
  
  static const std::string argv0 = [[[NSBundle mainBundle] executablePath] UTF8String];
  const char *xml_config_path = "data/config/wt_config_osx.xml";
  
  InterSpecServer::start_server( argv0.c_str(), user_data_dir.c_str(),
                                base_dir.c_str(), xml_config_path, server_port );
  
  //now we'll wait for the server to start
  std::string url = InterSpecServer::urlBeingServedOn();
  
  if( url.empty() )
  {
    _isServing = NO;
     std::cerr << "Unable to start server!" << std::endl;
  }else
  {
    _isServing = YES;
    _UrlServingOn = [NSString stringWithUTF8String:url.c_str()];
    
    //NSLog( @"Serving at URL=%@", _UrlServingOn );
    
    //[_window setTitle: [NSString stringWithFormat:@"TRB InterSpec - %@", _UrlServingOn]];
    [_window setTitle: @"InterSpec"];
    
    _UrlUniqueId = [self generateSessionToken];
    
    InterSpecServer::add_allowed_session_token( [_UrlUniqueId UTF8String], InterSpecServer::SessionType::PrimaryAppInstance );
    
    NSString *actualURL = [NSString stringWithFormat:@"%@?apptoken=%@&primary=yes", _UrlServingOn, _UrlUniqueId];
    
    
    if( _fileNeedsOpening )
    {
      NSLog( @"There isa  file that needs opening" );
      NSLog( @"UrlUniqueId %@", _UrlUniqueId );
      NSLog( @"fileNeedsOpening %@", _fileNeedsOpening );
      
      const std::string session_token = [_UrlUniqueId UTF8String];
      const std::string file_path = [_fileNeedsOpening UTF8String];
      
      InterSpecServer::set_file_to_open_on_load( session_token.c_str(), file_path );
      
      NSLog( @"Have called InterSpecServer::set_file_to_open_on_load" );
      
      _fileNeedsOpening = nil;
    }//if( fileNeedsOpening )
    
    //Note: this bit of code does not seem to work from the themeChanged notification.
#ifdef AVAILABLE_MAC_OS_X_VERSION_10_14_AND_LATER   //preproccessor to compile on older macOS
    if( @available(macOS 10.14, *) )  //runtime check to run on older macOS
    {
      NSAppearanceName basicAppearance
        = [NSApp.mainWindow.effectiveAppearance
           bestMatchFromAppearancesWithNames:@[ NSAppearanceNameAqua,
                                               NSAppearanceNameDarkAqua ]];
      if( [basicAppearance isEqualToString:NSAppearanceNameDarkAqua] )
      {
        NSLog( @"requesting initial theme to dark" );
        actualURL = [NSString stringWithFormat:@"%@&colortheme=dark", actualURL];
        //The initial window is still white-ish untill the page loads - below
        // causes an even worse flash of white->black->white->theme color
        //NSString *js = @"document.body.style.background = \"black\"";
        //[_InterSpecWebView evaluateJavaScript:js completionHandler:nil];
      }
    }//if( >= macOS 10.14 )
#endif
  
    
    if( !doResume )
      actualURL = [NSString stringWithFormat:@"%@&restore=no", actualURL];
    
    // WkWebView registers _window for registerForDraggedTypes; if it didnt, we could use the following line
    //[_window registerForDraggedTypes: [NSArray arrayWithObjects:NSFilenamesPboardType, NSURLPboardType, nil]];
    
    //if( [_InterSpecWebView respondsToSelector:@selector(mainFrame)])
    //[[_InterSpecWebView mainFrame] loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]]];
    [_InterSpecWebView loadRequest: [NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]] ];
    
    [_InterSpecWebView setNavigationDelegate: self];
    [_InterSpecWebView setUIDelegate: self];
    
    
    if (@available(macOS 10.14, *)) {
      [NSDistributedNotificationCenter.defaultCenter addObserver:self selector:@selector(themeChanged:) name:@"AppleInterfaceThemeChangedNotification" object: nil];
    }
  }
}//applicationDidFinishLaunching:(NSNotification *)aNotification


-(void)themeChanged:(NSNotification *) notification
{
  NSString *appearance = [[NSUserDefaults standardUserDefaults] stringForKey:@"AppleInterfaceStyle"];
  
  if( !appearance )
  {
    InterSpecApp::osThemeChange( "light" );
  }else if( [appearance rangeOfString:@"Dark" options:NSCaseInsensitiveSearch].location != NSNotFound )
  {
    InterSpecApp::osThemeChange( "dark" );
  }
}//themeChanged

//WKUIDelegate
//webView:createWebViewWithConfiguration:forNavigationAction:windowFeatures:
//webViewDidClose:
- (void)webView:(WKWebView *)webView
        runOpenPanelWithParameters:(WKOpenPanelParameters *)parameters
        initiatedByFrame:(WKFrameInfo *)frame
        completionHandler:(void (^)(NSArray<NSURL *> *URLs))completionHandler
{
  // Create the File Open Dialog class.
  NSOpenPanel *openDlg = [NSOpenPanel openPanel];
  
  bool chooseDirectory = false;
  
#if( USE_SPECRUM_FILE_QUERY_WIDGET )
  chooseDirectory = SpecFileQuery::isSelectingDirectory();
#endif
  
  if( chooseDirectory )
  {
    [openDlg setCanChooseFiles:NO];
    [openDlg setCanChooseDirectories:YES];
  }else
  {
    [openDlg setCanChooseFiles:YES];
    [openDlg setCanChooseDirectories:NO];
  }
  
  if( [openDlg runModal] == NSModalResponseOK )
  {
    NSArray* files = [[openDlg URLs]valueForKey:@"relativePath"];
    
#if( USE_SPECRUM_FILE_QUERY_WIDGET )
    if( chooseDirectory && [files count] )
    {
      NSURL *u = [files objectAtIndex:0];
      if( u && [u respondsToSelector:@selector(UTF8String)]) {
        SpecFileQuery::setSearchDirectory( [u UTF8String] );
      }
    }
#endif
    
    completionHandler( [openDlg URLs] );
  }else
  {
    completionHandler( nil );
  }
}


//Implement WKNavigationDelegate
- (void)webView:(WKWebView *)webView
                decidePolicyForNavigationAction:(WKNavigationAction *)navigationAction
                decisionHandler:(void (^)(WKNavigationActionPolicy))decisionHandler
{
  NSLog(@"decidePolicyForNavigationAction" );
  
  switch( [navigationAction navigationType] )
  {
    case WKNavigationTypeLinkActivated:  //external URL, or CSV, or spectrum file download.
    {
      NSLog(@"WKNavigationTypeLinkActivated" );
      NSURLRequest *request = [navigationAction request];
      //NSEventModifierFlags modifierFlags = [navigationAction modifierFlags];
      //NSInteger buttonNumber = [navigationAction buttonNumber];
      
      if( request )
      {
        NSString *host = [[request URL] host];
        NSString *absurl = [[request URL] absoluteString];
        NSString *scheme = [[request URL] scheme];
        
        NSLog(@"scheme=%@, host=%@, absurl=%@", scheme, host, absurl );
        
        if ([absurl rangeOfString:@"request=redirect&url=http"].location != NSNotFound
            || [scheme isEqual:@"mailto"]
            || ( ([scheme isEqual:@"https"] || [scheme isEqual:@"http"])
                && ![host isEqual:@"127.0.0.1"]
                && ([host rangeOfString:@"localhost"].location == NSNotFound)) )
        {
          //external url or email
          [[NSWorkspace sharedWorkspace] openURL:[request URL]];
        }else if( (host && [host isEqualToString:@"127.0.0.1"])
                 || (!host && absurl && [absurl hasPrefix: @"data:application/octet-stream"])
                 || (!host && absurl && [absurl hasPrefix: @"data:image/svg+xml"]) )
        {
          //CSV, spectrum file, JSON file, etc
          NSLog(@"Will attempt to download spectrum, CSV, JSON, PNG, etc. file");
         
          NSURLRequest *theRequest = [NSURLRequest requestWithURL:[[request URL] absoluteURL]];
          NSURLDownload  *theDownload = [[NSURLDownload alloc] initWithRequest:theRequest delegate:self];
          
          if( !theDownload )
            NSLog(@"The download failed");  // Inform the user that the download failed.
        }
      }//if( request )
      
      decisionHandler(WKNavigationActionPolicyCancel);
      break;
    }//case WKNavigationTypeLinkActivated:
      
    case WKNavigationTypeBackForward:
      NSLog(@"WKNavigationTypeBackForward" );
      decisionHandler(WKNavigationActionPolicyCancel);
      break;
    
    case WKNavigationTypeFormSubmitted:
      NSLog(@"WKNavigationTypeFormSubmitted" );
      decisionHandler(WKNavigationActionPolicyAllow);  //file uploads
      break;
      
    case WKNavigationTypeReload:
    case WKNavigationTypeFormResubmitted:
      NSLog(@"WKNavigationTypeFormSubmitted, WKNavigationTypeReload, or WKNavigationTypeFormResubmitted recieved - canceling navigation" );
      decisionHandler(WKNavigationActionPolicyCancel);
      break;
      
    case WKNavigationTypeOther:
      //Initial page load comes to here
      decisionHandler(WKNavigationActionPolicyAllow);
      break;
      
    default:
      NSLog(@"Unknown WKNavigationType" );
      decisionHandler(WKNavigationActionPolicyAllow);
      break;
  }
}


/*
//For WKUIDelegate, if we want to open a new window...
- (WKWebView *)webView:(WKWebView *)webView
createWebViewWithConfiguration:(WKWebViewConfiguration *)configuration
   forNavigationAction:(WKNavigationAction *)navigationAction
        windowFeatures:(WKWindowFeatures *)windowFeatures
{
  
}
*/


//Implement NSURLDownloadDelegate
- (void)download:(NSURLDownload *)download decideDestinationWithSuggestedFilename:(NSString *)filename
{
  //NSLog(@"decideDestinationWithSuggestedFilename: %@", filename );
  
  //The PNG screenshot has a data encoded URI that the obj-c doesnt seem to grab
  //  the filename from, so lets do that manually.  I'm sure there is a much
  //  more elegant way to do this... but it seems to work.
  if( !filename || [filename isEqualToString: @"Unknown"] )
  {
    NSURLRequest *request = [download request];
    NSURL *url = request ? [request URL] : (NSURL *)nil;
    NSString *absurl = url ? [url absoluteString] : @"";
    NSRange pos = [absurl rangeOfString:@"filename="];
    if( pos.location != NSNotFound )
    {
      NSRange fnamerange;
      fnamerange.location = pos.location + pos.length;
      fnamerange.length = [absurl length] - fnamerange.location;
      
      NSRange semipos = [absurl rangeOfString: @";" options: 0 range: fnamerange locale: nil];
      
      if( semipos.location != NSNotFound )
      {
        fnamerange.length = semipos.location - fnamerange.location;
        filename = [absurl substringWithRange: fnamerange];
      }
    }
  }//if( filename is unknown )
  
  NSSavePanel *panel = [NSSavePanel savePanel];
  [panel setNameFieldStringValue:filename];
  
  NSInteger result = [panel runModal];
  if( result == NSFileHandlingPanelOKButton)
  {
    NSURL *path = [panel URL];
    NSLog(@"Saving to %@", [path path]);
    [download setDestination:[path path] allowOverwrite:NO];
  }else
  {
    [download cancel];
  }
}

//Implement NSURLDownloadDelegate
- (void)download:(NSURLDownload *)download didFailWithError:(NSError *)error
{
  NSLog(@"- (void)download:(NSURLDownload *)download didFailWithError:(NSError *)error");
  NSLog(@"Download failed! Error - %@ %@",
        [error localizedDescription],
        [[error userInfo] objectForKey:NSURLErrorFailingURLStringErrorKey]);
}

//Implement NSURLDownloadDelegate
- (void)downloadDidFinish:(NSURLDownload *)download
{
  //NSData *data = [[download request] HTTPBody];
  //NSUInteger length = [data length];
  //  long long expectedLength = [[download request] expectedContentLength];
  //std::cout << "length=" << length << std::endl;
  NSLog(@"- (void)downloadDidFinish:(NSURLDownload *)download: ");
}

//Implement NSURLDownloadDelegate
- (void)startDownloadingURL:sender
{
  NSLog(@"- (void)startDownloadingURL:sender");
}

//Implement NSURLDownloadDelegate
- (void)downloadDidBegin:(NSURLDownload *)download
{
  NSLog(@"downloadDidBegin");
}

//Implement NSURLDownloadDelegate
-(void)download:(NSURLDownload *)download didCreateDestination:(NSString *)path
{
  // path now contains the destination path
  // of the download, taking into account any
  // unique naming caused by -setDestination:allowOverwrite:
  NSLog(@"Final file destination: %@",path);
}




// Returns the directory the application uses to store the Core Data store file. This code uses a directory named "sandia.InterSpec" in the user's Application Support directory.
- (NSURL *)applicationFilesDirectory
{ 
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSURL *appSupportURL = [[fileManager URLsForDirectory:NSApplicationSupportDirectory inDomains:NSUserDomainMask] lastObject];
    return [appSupportURL URLByAppendingPathComponent:@"sandia.InterSpec"];
}

// Creates if necessary and returns the managed object model for the application.
- (NSManagedObjectModel *)managedObjectModel
{
  if( _managedObjectModel )
    return _managedObjectModel;
	
  NSURL *modelURL = [[NSBundle mainBundle] URLForResource:@"InterSpec" withExtension:@"momd"];
  _managedObjectModel = [[NSManagedObjectModel alloc] initWithContentsOfURL:modelURL];
  return _managedObjectModel;
}

// Returns the persistent store coordinator for the application. This implementation creates and return a coordinator, having added the store for the application to it. (The directory for the store is created, if necessary.)
- (NSPersistentStoreCoordinator *)persistentStoreCoordinator
{
    if( _persistentStoreCoordinator )
      return _persistentStoreCoordinator;
    
    NSManagedObjectModel *mom = [self managedObjectModel];
    if (!mom) {
        NSLog(@"%@:%@ No model to generate a store from", [self class], NSStringFromSelector(_cmd));
        return nil;
    }
    
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSURL *applicationFilesDirectory = [self applicationFilesDirectory];
    NSError *error = nil;
    
    NSDictionary *properties = [applicationFilesDirectory resourceValuesForKeys:@[NSURLIsDirectoryKey] error:&error];
    
    if (!properties) {
        BOOL ok = NO;
        if ([error code] == NSFileReadNoSuchFileError) {
            ok = [fileManager createDirectoryAtPath:[applicationFilesDirectory path] withIntermediateDirectories:YES attributes:nil error:&error];
        }
        if (!ok) {
            [[NSApplication sharedApplication] presentError:error];
            return nil;
        }
    } else {
//        if (![properties[NSURLIsDirectoryKey] boolValue]) {
            // Customize and localize this error.
//            NSString *failureDescription = [NSString stringWithFormat:@"Expected a folder to store application data, found a file (%@).", [applicationFilesDirectory path]];
      
//            NSMutableDictionary *dict = [NSMutableDictionary dictionary];
//            [dict setValue:failureDescription forKey:NSLocalizedDescriptionKey];
//            error = [NSError errorWithDomain:@"YOUR_ERROR_DOMAIN" code:101 userInfo:dict];
      
//            [[NSApplication sharedApplication] presentError:error];
//            return nil;
//        }
    }
    
    NSURL *url = [applicationFilesDirectory URLByAppendingPathComponent:@"InterSpec.storedata"];
    NSPersistentStoreCoordinator *coordinator = [[NSPersistentStoreCoordinator alloc] initWithManagedObjectModel:mom];
    if (![coordinator addPersistentStoreWithType:NSXMLStoreType configuration:nil URL:url options:nil error:&error]) {
        [[NSApplication sharedApplication] presentError:error];
        return nil;
    }
    _persistentStoreCoordinator = coordinator;
    
    return _persistentStoreCoordinator;
}

// Returns the managed object context for the application (which is already bound to the persistent store coordinator for the application.) 
- (NSManagedObjectContext *)managedObjectContext
{
  if( _managedObjectContext )
    return _managedObjectContext;
    
  NSPersistentStoreCoordinator *coordinator = [self persistentStoreCoordinator];
  if (!coordinator)
  {
      NSMutableDictionary *dict = [NSMutableDictionary dictionary];
      [dict setValue:@"Failed to initialize the store" forKey:NSLocalizedDescriptionKey];
      [dict setValue:@"There was an error building up the data file." forKey:NSLocalizedFailureReasonErrorKey];
      NSError *error = [NSError errorWithDomain:@"YOUR_ERROR_DOMAIN" code:9999 userInfo:dict];
      [[NSApplication sharedApplication] presentError:error];
      return nil;
  }
  _managedObjectContext = [[NSManagedObjectContext alloc] init];
  [_managedObjectContext setPersistentStoreCoordinator:coordinator];

  return _managedObjectContext;
}

// Returns the NSUndoManager for the application. In this case, the manager returned is that of the managed object context for the application.
- (NSUndoManager *)windowWillReturnUndoManager:(NSWindow *)window
{
  return [[self managedObjectContext] undoManager];
}


- (NSApplicationTerminateReply)applicationShouldTerminate:(NSApplication *)sender
{
  // Save changes in the application's managed object context before the application terminates.
  
  NSLog( @"Terminated %p", _window );
  
  if (!_managedObjectContext)
  {
    InterSpecServer::killServer();
    return NSTerminateNow;
  }
    
  if (![[self managedObjectContext] commitEditing])
  {
    InterSpecServer::killServer();
    NSLog(@"%@:%@ unable to commit editing to terminate", [self class], NSStringFromSelector(_cmd));
    return NSTerminateCancel;
  }
    
  if( ![[self managedObjectContext] hasChanges] )
  {
    InterSpecServer::killServer();
    return NSTerminateNow;
  }
    
  NSError *error = nil;
  if( ![[self managedObjectContext] save:&error] )
  {
    // Customize this code block to include application-specific recovery steps.
    BOOL result = [sender presentError:error];
    if (result)
      return NSTerminateCancel;

    NSString *question = NSLocalizedString(@"Could not save changes while quitting. Quit anyway?", @"Quit without saves error question message");
    NSString *info = NSLocalizedString(@"Quitting now will lose any changes you have made since the last successful save", @"Quit without saves error question info");
    NSString *quitButton = NSLocalizedString(@"Quit anyway", @"Quit anyway button title");
    NSString *cancelButton = NSLocalizedString(@"Cancel", @"Cancel button title");
    NSAlert *alert = [[NSAlert alloc] init];
    [alert setMessageText:question];
    [alert setInformativeText:info];
    [alert addButtonWithTitle:quitButton];
    [alert addButtonWithTitle:cancelButton];

    NSInteger answer = [alert runModal];
        
    if( answer == NSAlertFirstButtonReturn )
      return NSTerminateCancel;
  }

  InterSpecServer::killServer();
  return NSTerminateNow;
}

- (void)applicationWillTerminate:(NSNotification *)notification
{
  //This is the last function called before termination; not called if force-quit
  NSLog( @"applicationWillTerminate" );
  
  NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
  [defaults setBool: YES forKey: @"DoResume"];
  
  NSLog( @"applicationWillTerminate: setting to YES, DoResume=%i", [defaults boolForKey:@"DoResume"] );
}


// This method is for WKScriptMessageHandler, and is triggered each time 'interOp' is sent a
//    message from the JavaScript code with something like:
//      window.webkit.messageHandlers.interOp.postMessage({"val": 1});
- (void)userContentController:(WKUserContentController *)userContentController
      didReceiveScriptMessage:(WKScriptMessage *)message{
  NSDictionary *sentData = (NSDictionary*)message.body;
  
  if( sentData == nil )
  {
    NSLog( @"didReceiveScriptMessage: got nil dictionary" );
    return;
  }
  
  id val = [sentData objectForKey: @"action"];
  if( val != (id)[NSNull null] )
  {
    NSString *str = (NSString*)val;
    if( !str )
    {
      NSLog( @"didReceiveScriptMessage: no value for action" );
      return;
    }
    
    if( [str isEqualToString:@"ExternalInstance"] )
    {
      NSLog( @"didReceiveScriptMessage: request for external instance" );
      
      NSString *sessionToken = [self generateSessionToken];
      InterSpecServer::add_allowed_session_token( [sessionToken UTF8String], InterSpecServer::SessionType::ExternalBrowserInstance );
      
      //url = InterSpecServer::urlBeingServedOn(); //will be empty if not serving
      const NSInteger port = InterSpecServer::portBeingServedOn();
      NSString *url = [NSString stringWithFormat:@"http://localhost:%ld?apptoken=%@&restore=no", port, sessionToken];
      //NSString *url = [NSString stringWithFormat:@"%@?primary=no&restore=no", _UrlServingOn];
      [[NSWorkspace sharedWorkspace] openURL: [NSURL URLWithString:url] ];
    }else
    {
      NSLog( @"didReceiveScriptMessage: un-understood action: \"%@\"", str );
    }
  }else
  {
    NSLog( @"didReceiveScriptMessage: No \"action\" key" );
  }
  
  // Could call back into the WebView using
  //[_InterSpecWebView evaluateJavaScript: @"someJavascript" completionHandler:nil];
}//didReceiveScriptMessage


-(void)enableWebInspector
{
  NSLog( @"Will show WebInpector" );
  
  // For macOS build, the context menu is disabled on wApp->domRoot() using JavaScript, so we need
  //  to over-ride this
  NSString *js = @"document.querySelector('.Wt-domRoot').setAttribute('oncontextmenu','return true;');";
  [_InterSpecWebView evaluateJavaScript:js completionHandler:nil];
    
  NSAlert *alert = [[NSAlert alloc] init];
  [alert setMessageText:@"Web Inspector Enabled"];
  [alert setInformativeText: @"Right click somewhere other than the spectrum or other chart, "
                              "and select 'Inspect Element' from the contect menu" ];
  [alert addButtonWithTitle:@"Ok"];
  [alert setAlertStyle:NSAlertStyleInformational];
  [alert runModal];
  
  // TODO: disable, or remove the "Enable Web Inspector" menu item.
}//enableWebInspector


-(NSString *)generateSessionToken
{
#define TOKEN_LEN 14
  char data[TOKEN_LEN];
  for( int index = 0; index < TOKEN_LEN; ++index )
    data[index] = (char)('A' + (arc4random_uniform(26)));
  
  return [[NSString alloc] initWithBytes:data length:TOKEN_LEN encoding:NSUTF8StringEncoding];
}//NSString *generateSessionToken

@end
