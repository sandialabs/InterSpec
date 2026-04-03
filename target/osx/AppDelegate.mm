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
      
      // Use copy to ensure we have a strong reference that outlives this method
      _fileNeedsOpening = [[NSString alloc] initWithUTF8String:urlcontent.c_str()];
      
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
  webConfig.applicationNameForUserAgent = @"InterSpec";
  webConfig.suppressesIncrementalRendering = YES;
  
  WKPreferences *prefs = [webConfig preferences];
  prefs.javaScriptCanOpenWindowsAutomatically = YES;
  prefs.tabFocusesLinks = NO;
  
  // Make it so the obj-c didReceiveScriptMessage method will be called any time the javascript
  //  does something like:
  //      window.webkit.messageHandlers.interOp.postMessage({"action": "DoSomething"});
  [webConfig.userContentController addScriptMessageHandler: self name:@"interOp"];

  // Add a script message handler
  [webConfig.userContentController addScriptMessageHandler:self name:@"jsErrorHandler"];

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
  
  _InterSpecWebView.customUserAgent = @"InterSpec";
  _InterSpecWebView.allowsLinkPreview = NO;
  _InterSpecWebView.allowsBackForwardNavigationGestures = NO;
  
  
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
      
      // Add separator before session debug items
      [menu addItem:[NSMenuItem separatorItem]];
      
      // Add "Simulate WebView Termination" menu item
      NSMenuItem *simulateTerminate = [[NSMenuItem alloc]
                                       initWithTitle:@"Simulate WebView Termination"
                                       action:@selector(simulateWebViewTermination)
                                       keyEquivalent:@""];
      [simulateTerminate setTarget:self];
      [menu addItem:simulateTerminate];
      [simulateTerminate setEnabled:YES];
      
      // Add "Simulate Session Invalidation" menu item
      NSMenuItem *simulateInvalidate = [[NSMenuItem alloc]
                                        initWithTitle:@"Simulate Session Invalidation"
                                        action:@selector(simulateSessionInvalidation)
                                        keyEquivalent:@""];
      [simulateInvalidate setTarget:self];
      [menu addItem:simulateInvalidate];
      [simulateInvalidate setEnabled:YES];
      
      // Add "Force Session Validation" menu item (one-time check)
      NSMenuItem *forceValidation = [[NSMenuItem alloc]
                                     initWithTitle:@"Force Session Validation"
                                     action:@selector(periodicSessionHealthCheck)
                                     keyEquivalent:@""];
      [forceValidation setTarget:self];
      [menu addItem:forceValidation];
      [forceValidation setEnabled:YES];
      
      // Add "Start Periodic Health Check" menu item (toggles 60-second timer)
      NSMenuItem *togglePeriodicCheck = [[NSMenuItem alloc]
                                         initWithTitle:@"Start Periodic Health Check (60s)"
                                         action:@selector(togglePeriodicHealthCheck)
                                         keyEquivalent:@""];
      [togglePeriodicCheck setTarget:self];
      [menu addItem:togglePeriodicCheck];
      [togglePeriodicCheck setEnabled:YES];
    }
    
    //Note: currently the right click in InterSpecApp is disabled using javascript if
    //      unles PERFORM_DEVELOPER_CHECKS is true.  However, could also do:
    //https://stackoverflow.com/questions/28801032/how-can-the-context-menu-in-wkwebview-on-the-mac-be-modified-or-overridden#28981319
    //[_InterSpecWebView willOpenMenu:<#(nonnull NSMenu *)#> withEvent:<#(nonnull NSEvent *)#>];
  }//if( allow_dev_tools )
  
  InterSpecApp::setNativeFileSaveHandler( []( std::string data, std::string suggested_name ){
    NSData *nsData = [NSData dataWithBytes:data.data() length:data.size()];
    NSString *nsSuggestedName = [NSString stringWithUTF8String:suggested_name.c_str()];

    dispatch_async( dispatch_get_main_queue(), ^{
      NSSavePanel *savePanel = [NSSavePanel savePanel];
      savePanel.title = @"Save Chart Image";
      savePanel.prompt = @"Save";
      savePanel.nameFieldStringValue = nsSuggestedName;
      savePanel.canCreateDirectories = YES;

      [savePanel beginWithCompletionHandler:^( NSModalResponse result ){
        if( result != NSModalResponseOK )
          return;

        NSURL *destinationURL = savePanel.URL;
        if( !destinationURL )
          return;

        NSError *writeError = nil;
        const BOOL written = [nsData writeToURL:destinationURL
                                        options:NSDataWritingAtomic
                                          error:&writeError];
        if( !written )
        {
          NSLog( @"Failed to save chart image: %@", writeError.localizedDescription );
          dispatch_async( dispatch_get_main_queue(), ^{
            NSAlert *alert = [[NSAlert alloc] init];
            alert.messageText = @"Error Saving File";
            alert.informativeText = writeError.localizedDescription;
            alert.alertStyle = NSAlertStyleCritical;
            [alert addButtonWithTitle:@"OK"];
            [alert runModal];
          } );
        }//if( !written )
      }];//beginWithCompletionHandler
    } );//dispatch_async
  } );//setNativeFileSaveHandler

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
    
    /*
     //Right now drag-n-drop from Outlook does not work - should investigate this:
     // Probably need to sub-class WKWebView, and then do something like
     //https://developer.apple.com/library/archive/documentation/Cocoa/Conceptual/DragandDrop/Tasks/acceptingdrags.html#//apple_ref/doc/uid/20000993-BABHHIHC
     */
    
    
    
    

    //[_window registerForDraggedTypes: [NSArray arrayWithObjects:NSFilenamesPboardType, NSURLPboardType, nil]]; //requires macOS >=10.13, which is the lowest we target anyway
    [_window registerForDraggedTypes: [NSArray arrayWithObjects:NSPasteboardTypeFileURL, NSPasteboardTypeURL, NSFileContentsPboardType, nil]]; //requires macOS >=10.13, which is the lowest we target anyway
    
    
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

// This method is called when the web view finishes loading a page
- (void)webView:(WKWebView *)webView didFinishNavigation:(WKNavigation *)navigation {
  NSLog(@"WebView has finished loading (could be a refresh).");
  
  [self validateCurrentSession];
  
  if( _sessionIsValid )
  {
    _sessionRecoveryAttempts = 0; // Reset on successful load
    NSLog(@"Session validated successfully after navigation");
  }
  else
  {
    NSLog(@"Session invalid after navigation - will attempt recovery");
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Session Validation Failed After Navigation";
    alert.informativeText = [NSString stringWithFormat:
      @"After navigation finished, session token '%@' is no longer valid. "
      "Recovery will be attempted.", _UrlUniqueId];
    alert.alertStyle = NSAlertStyleWarning;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
    [self recoverSession];
  }
}

// Called when the WebView's content process is terminated by the system
- (void)webViewWebContentProcessDidTerminate:(WKWebView *)webView
{
  NSLog(@"WebView content process terminated by system!");
  
#if( PERFORM_DEVELOPER_CHECKS )
  NSAlert *alert = [[NSAlert alloc] init];
  alert.messageText = @"WebView Process Terminated";
  alert.informativeText = @"The WebView content process was terminated by the system. "
                          "This is likely the cause of session disconnection issues. "
                          "Attempting automatic recovery.";
  alert.alertStyle = NSAlertStyleWarning;
  [alert addButtonWithTitle:@"OK"];
  [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
  
  [self recoverSession];
}

// Called when navigation fails during the early stages (before response received)
- (void)webView:(WKWebView *)webView
    didFailProvisionalNavigation:(WKNavigation *)navigation
    withError:(NSError *)error
{
  NSLog(@"WebView provisional navigation failed: %@", error.localizedDescription);
  
#if( PERFORM_DEVELOPER_CHECKS )
  NSAlert *alert = [[NSAlert alloc] init];
  alert.messageText = @"Navigation Failed (Provisional)";
  alert.informativeText = [NSString stringWithFormat:
    @"Error Code: %ld\nDomain: %@\nDescription: %@",
    (long)error.code, error.domain, error.localizedDescription];
  alert.alertStyle = NSAlertStyleWarning;
  [alert addButtonWithTitle:@"OK"];
  [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
}

// Called when navigation fails after response has been received
- (void)webView:(WKWebView *)webView
    didFailNavigation:(WKNavigation *)navigation
    withError:(NSError *)error
{
  NSLog(@"WebView navigation failed: %@", error.localizedDescription);
  
#if( PERFORM_DEVELOPER_CHECKS )
  NSAlert *alert = [[NSAlert alloc] init];
  alert.messageText = @"Navigation Failed";
  alert.informativeText = [NSString stringWithFormat:
    @"Error Code: %ld\nDomain: %@\nDescription: %@",
    (long)error.code, error.domain, error.localizedDescription];
  alert.alertStyle = NSAlertStyleWarning;
  [alert addButtonWithTitle:@"OK"];
  [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
}

- (void)applicationDidEnterBackground:(NSNotification *)notification {
    NSLog(@"App entered background.");
}

- (void)applicationWillEnterForeground:(NSNotification *)notification {
    NSLog(@"App will enter foreground.");
}

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
  
  // You can allow the WKWebView to select directories, using the Wt::WFileUpload,
  //  by setting the 'webkitdirectory' attribute on the input element of the WFileUpload,
  //  through a call like:
  //    upload->doJavaScript( "document.querySelector('#" + upload->id() + " input').setAttribute('webkitdirectory', true);" );
  //    (and similarly to allow multiple files with the attribute 'multiple')
  //  The user can then select directories, and the following will will correctly get what we want,
  //  and below in the obj-c we can get the path of the directory selected, and pass it off to other
  //  parts of the program, or the DOM, or whatever.
  //  But the cleaner way of doing things is to call `macOsUtils::showFilePicker(...)`,
  //  Leaving in this mechanism, and note in the code for now because it is likely we will
  //  want a multi-file and/or directory selector, that uses a WFileUpload as the presentation
  //  to the user, so it would be reasonable to create a class that inherits from, or contains,
  //  WFileUpload and will just kinda take care of providing native filesystem paths normally,
  //  but will fallback to html upload.
  const BOOL chooseDirectory = [parameters allowsDirectories];
  const BOOL multiSelect = [parameters allowsMultipleSelection]; //Untested
  
  if( chooseDirectory )
  {
    [openDlg setCanChooseFiles:NO];
    [openDlg setCanChooseDirectories:YES];
  }else
  {
    [openDlg setCanChooseFiles:YES];
    [openDlg setCanChooseDirectories:NO];
  }
  
  if( multiSelect )
    [openDlg setAllowsMultipleSelection: YES];
  
  if( [openDlg runModal] == NSModalResponseOK )
  {
    NSArray<NSURL *> *urls = [openDlg URLs];
    
    /*
    {//Begin get user selected paths, and set to DOM, and wherever else in memory wanted
      Wt::Json::Array json_array;
      
      // Form some JSON to set a variable in JS client-side, that can then make the call to server code
      for (NSURL *url in urls)
      {
        const char *fs_path = [url fileSystemRepresentation];  //Will be UTF-8 encoded
        json_array.push_back( Wt::WString::fromUTF8(fs_path) );
      }
      
      if( json_array.empty() )
      {
        //Probably got file promise(s) here - we'll just fall-back to normal file upload, and
        //  let WkWebView take care of getting the files.
        NSString *js = @"(function(){$(document).data('SelectedPaths',null);})()";
        [_InterSpecWebView evaluateJavaScript: js completionHandler:nil];
        
        NSLog( @"Cleared path data to JS." );
      }else
      {
        const std::string json_data = Wt::Json::serialize(json_array);
        const std::string js_str =
        "(function(){\n\t"
        "let fns = {};\n\t"
        "fns.time = new Date();\n\t"
        "fns.isDir = " + std::string(chooseDirectory ? "1" : "0") + ";"
        "fns.filenames = " + json_data + ";\n\t"
        "$(document).data('SelectedPaths',fns);\n\t"
        "console.log(\"Set SelectedPaths paths\",fns);\n"
        "})();";
        
        NSString *js = [[NSString alloc] initWithCString:js_str.c_str() encoding:NSUTF8StringEncoding];
        [_InterSpecWebView evaluateJavaScript: js completionHandler:nil];
        
        NSLog( @"Set path data to JS=%s", js_str.c_str() );
      }//if( !json_array.empty() )
      
      // Set to wherever else in memory wanted
    }//end get user selected paths, and set to DOM, and wherever else in memory wanted
    */
    
    completionHandler( urls );
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

          //Using NSURLDownload is deprecated in 10.11, so we use NSURLSession instead.
          //The download goes to a temporary location, then we show a save panel and move it directly to the destination.
          //This avoids the file collision issues that occurred with the previous intermediate copy approach.

          NSURL *url = [[request URL] absoluteURL];
          NSURLSession *session = [NSURLSession sharedSession];

          // Create a download task
          NSURLSessionDownloadTask *downloadTask = [session downloadTaskWithURL:url
            completionHandler:^(NSURL *location, NSURLResponse *response, NSError *error) {
              if( error )
              {
                // Display an error alert to the user
                dispatch_async(dispatch_get_main_queue(), ^{
                  NSAlert *alert = [[NSAlert alloc] init];
                  alert.messageText = @"Download Failed";
                  alert.informativeText = error.localizedDescription;
                  alert.alertStyle = NSAlertStyleCritical;
                  [alert addButtonWithTitle:@"OK"];
                  [alert runModal];
                });
                return;
              }else
              {
                NSLog(@"Download succeeded. File is located at: %@", location.path);

                // The temp file is only guaranteed to exist during this completion handler.
                // We must move it to a persistent location immediately before showing the save panel.
                NSFileManager *fileManager = [NSFileManager defaultManager];
                
                // Create a unique filename to avoid conflicts with existing files
                NSString *originalFilename = response.suggestedFilename ?: @"downloadedFile";
                NSString *filenameWithoutExtension = [originalFilename stringByDeletingPathExtension];
                NSString *fileExtension = [originalFilename pathExtension];
                NSUUID *uuid = [NSUUID UUID];
                NSString *uniqueFilename = [NSString stringWithFormat:@"%@_%@", filenameWithoutExtension, [uuid UUIDString]];
                if( [fileExtension length] > 0 )
                {
                  uniqueFilename = [uniqueFilename stringByAppendingPathExtension:fileExtension];
                }
                
                NSURL *persistentTempURL = [NSURL fileURLWithPath:[NSTemporaryDirectory() stringByAppendingPathComponent:uniqueFilename]];

                // Move the file to a persistent location immediately (before completion handler returns)
                NSError *moveError = nil;
                [fileManager moveItemAtURL:location toURL:persistentTempURL error:&moveError];

                if( moveError )
                {
                  NSLog(@"Failed to move temporary file to persistent location: %@", moveError.localizedDescription);
                  dispatch_async(dispatch_get_main_queue(), ^{
                    NSAlert *alert = [[NSAlert alloc] init];
                    alert.messageText = @"Error Saving File";
                    alert.informativeText = moveError.localizedDescription;
                    alert.alertStyle = NSAlertStyleCritical;
                    [alert addButtonWithTitle:@"OK"];
                    [alert runModal];
                  });
                  return;
                }//if( moveError )

                // Now we can safely show the save panel - the file is in a persistent location
                dispatch_async(dispatch_get_main_queue(), ^{
                  NSSavePanel *savePanel = [NSSavePanel savePanel];
                  savePanel.title = @"Save Exported File";
                  savePanel.prompt = @"Save";
                  savePanel.nameFieldStringValue = response.suggestedFilename ?: @"downloadedFile";
                  savePanel.canCreateDirectories = YES;

                  [savePanel beginWithCompletionHandler:^(NSModalResponse result) {
                    if( result == NSModalResponseOK )
                    {
                      NSURL *destinationURL = savePanel.URL;
                      if( destinationURL )
                      {
                        // Ensure parent directory exists
                        NSURL *parentDir = [destinationURL URLByDeletingLastPathComponent];
                        NSError *dirError = nil;
                        [fileManager createDirectoryAtURL:parentDir
                                    withIntermediateDirectories:YES
                                    attributes:nil
                                    error:&dirError];

                        if( dirError )
                        {
                          NSLog(@"Failed to create parent directory: %@", dirError.localizedDescription);
                          // Clean up persistent temp file
                          [fileManager removeItemAtURL:persistentTempURL error:nil];
                          dispatch_async(dispatch_get_main_queue(), ^{
                            NSAlert *alert = [[NSAlert alloc] init];
                            alert.messageText = @"Error Creating Directory";
                            alert.informativeText = dirError.localizedDescription;
                            alert.alertStyle = NSAlertStyleCritical;
                            [alert addButtonWithTitle:@"OK"];
                            [alert runModal];
                          });
                          return;
                        }//if( dirError )

                        // Move from persistent temp location to final destination
                        // Use replaceItemAtURL to properly handle overwriting existing files (preserves permissions in sandboxed apps)
                        NSError *finalMoveError = nil;
                        NSURL *resultingURL = nil;
                        
                        if( [fileManager fileExistsAtPath:[destinationURL path]] )
                        {
                          // File exists - use replaceItemAtURL to overwrite it properly
                          [fileManager replaceItemAtURL:destinationURL
                                           withItemAtURL:persistentTempURL
                                          backupItemName:nil
                                                 options:NSFileManagerItemReplacementUsingNewMetadataOnly
                                        resultingItemURL:&resultingURL
                                                   error:&finalMoveError];
                        }else
                        {
                          // File doesn't exist - just move it
                          [fileManager moveItemAtURL:persistentTempURL toURL:destinationURL error:&finalMoveError];
                        }//if( file exists ) / else

                        if( finalMoveError )
                        {
                          NSLog(@"Error saving to: %@", destinationURL.path);
                          NSLog(@"Failed to save file: %@", finalMoveError.localizedDescription);
                          // Clean up persistent temp file on error
                          [fileManager removeItemAtURL:persistentTempURL error:nil];
                          
                          dispatch_async(dispatch_get_main_queue(), ^{
                            NSAlert *alert = [[NSAlert alloc] init];
                            alert.messageText = @"Error Saving File";
                            alert.informativeText = finalMoveError.localizedDescription;
                            alert.alertStyle = NSAlertStyleCritical;
                            [alert addButtonWithTitle:@"OK"];
                            [alert runModal];
                          });
                        }else
                        {
                          NSLog(@"File saved to: %@", destinationURL.path);
                        }//if( finalMoveError ) / else
                      }//if( destinationURL )
                    }else
                    {
                      // User canceled - clean up the persistent temp file
                      NSLog(@"User canceled saving the file.");
                      NSError *deleteError = nil;
                      [fileManager removeItemAtURL:persistentTempURL error:&deleteError];
                      if( deleteError )
                      {
                        NSLog(@"Failed to delete persistent temporary file: %@", deleteError.localizedDescription);
                      }
                    }//if( result == NSModalResponseOK ) / else
                  }]; //[savePanel beginWithCompletionHandler:^(NSModalResponse result) {...
                });//dispatch_async(dispatch_get_main_queue(), ^{...
              }//if( error ) / else
            }];//NSURLSessionDownloadTask *downloadTask....completionHandler:^(NSURL *location, NSURLResponse *response, NSError *error) {....

          // Start the download task
          [downloadTask resume];
        }//if( (host && [host isEqualToString:@"127.0.0.1"]) || ....
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



// Returns the directory the application uses to store user data.
// Uses a directory named "sandia.InterSpec" in the user's Application Support directory.
- (NSURL *)applicationFilesDirectory
{ 
  NSFileManager *fileManager = [NSFileManager defaultManager];
  NSURL *appSupportURL = [[fileManager URLsForDirectory:NSApplicationSupportDirectory inDomains:NSUserDomainMask] lastObject];
  NSURL *appFilesDir = [appSupportURL URLByAppendingPathComponent:@"sandia.InterSpec"];
  
  // Ensure the directory exists
  NSError *error = nil;
  NSDictionary *properties = [appFilesDir resourceValuesForKeys:@[NSURLIsDirectoryKey] error:&error];
  
  if( !properties )
  {
    // Directory doesn't exist, try to create it
    if( [error code] == NSFileReadNoSuchFileError )
    {
      NSError *createError = nil;
      BOOL created = [fileManager createDirectoryAtURL:appFilesDir
                           withIntermediateDirectories:YES
                                            attributes:nil
                                                 error:&createError];
      if( !created )
      {
        NSLog(@"Failed to create application files directory: %@", createError.localizedDescription);
      }
    }
  }else if( ![properties[NSURLIsDirectoryKey] boolValue] )
  {
    // A file exists at the path where we expected a directory
    NSLog(@"Expected a folder to store application data, found a file at: %@", [appFilesDir path]);
    
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Application Data Error";
    alert.informativeText = [NSString stringWithFormat:
      @"Expected a folder to store application data, but found a file at:\n%@\n\n"
      "Please remove this file and restart the application.",
      [appFilesDir path]];
    alert.alertStyle = NSAlertStyleCritical;
    [alert addButtonWithTitle:@"OK"];
    [alert runModal];
  }
  
  return appFilesDir;
}


// Returns nil since InterSpec uses its own undo/redo system via UndoRedoManager
- (NSUndoManager *)windowWillReturnUndoManager:(NSWindow *)window
{
  return nil;
}


- (NSApplicationTerminateReply)applicationShouldTerminate:(NSApplication *)sender
{
  NSLog( @"applicationShouldTerminate called" );
  
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
//      window.webkit.messageHandlers.jsErrorHandler.postMessage({"message": message, "source": source, ...} );
- (void)userContentController:(WKUserContentController *)userContentController
      didReceiveScriptMessage:(WKScriptMessage *)message
{
  // Validate that message.body is a dictionary
  if( ![message.body isKindOfClass:[NSDictionary class]] )
  {
    NSLog( @"didReceiveScriptMessage: expected NSDictionary, got %@", [message.body class] );
    return;
  }
  
  NSDictionary *sentData = (NSDictionary *)message.body;

  if( [message.name isEqualToString:@"jsErrorHandler"] )
  {
    NSString *errorMessage = sentData[@"message"];
    NSString *source = sentData[@"source"];
    NSNumber *lineNumber = sentData[@"lineno"];
    NSNumber *columnNumber = sentData[@"colno"];
    NSString *errorDetails = sentData[@"error"];

    // Create the alert as a sheet (non-blocking) so multiple errors don't queue up modal dialogs
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"JavaScript Error";
    alert.informativeText = [NSString stringWithFormat:
      @"Unexpected JavaScript Error - please consider reporting to InterSpec@sandia.gov and restarting app:\n\n"
      "Message: %@\nSource: %@\nLine: %@\nColumn: %@\nDetails: %@",
      errorMessage ?: @"Unknown",
      source ?: @"Unknown",
      lineNumber ?: @0,
      columnNumber ?: @0,
      errorDetails ?: @"No additional details"];
    alert.alertStyle = NSAlertStyleWarning;
    [alert addButtonWithTitle:@"OK"];
    
    // Use sheet instead of modal to avoid blocking
    [alert beginSheetModalForWindow:_window completionHandler:nil];

    return;
  }//if( a JS error )

  id val = sentData[@"action"];
  if( val && val != (id)[NSNull null] && [val isKindOfClass:[NSString class]] )
  {
    NSString *action = (NSString *)val;
    
    if( [action isEqualToString:@"ExternalInstance"] )
    {
      NSLog( @"didReceiveScriptMessage: request for external instance" );
      
      NSString *sessionToken = [self generateSessionToken];
      InterSpecServer::add_allowed_session_token( [sessionToken UTF8String], InterSpecServer::SessionType::ExternalBrowserInstance );
      
      const NSInteger port = InterSpecServer::portBeingServedOn();
      NSString *url = [NSString stringWithFormat:@"http://localhost:%ld?apptoken=%@&restore=no", port, sessionToken];
      [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:url]];
    }else
    {
      NSLog( @"didReceiveScriptMessage: unrecognized action: \"%@\"", action );
    }
  }else
  {
    NSLog( @"didReceiveScriptMessage: missing or invalid \"action\" key" );
  }
}//didReceiveScriptMessage


#pragma mark - Session Validation and Recovery

- (void)validateCurrentSession
{
  if( !_UrlUniqueId || !_isServing )
  {
    _sessionIsValid = NO;
    _lastSessionValidation = [NSDate date];
    NSLog(@"validateCurrentSession: No token or not serving, sessionIsValid=NO");
    return;
  }
  
  const int status = InterSpecServer::session_status( [_UrlUniqueId UTF8String] );
  
  // Status: 0=unknown, 1=authorized not loaded, 2=loaded, 3=no longer auth, 4=dead
  _sessionIsValid = (status == 2);
  _lastSessionValidation = [NSDate date];
  
  NSLog(@"validateCurrentSession: token='%@', status=%d, valid=%@",
        _UrlUniqueId, status, _sessionIsValid ? @"YES" : @"NO");
  
  if( !_sessionIsValid && status >= 3 )
  {
    NSLog(@"Session is dead (status=%d) or no longer authorized - recovery needed", status);
  }
}


- (void)recoverSession
{
  NSLog(@"Attempting session recovery (attempt %ld)", (long)++_sessionRecoveryAttempts);
  
  // Always get the current server URL fresh from InterSpecServer
  // This is more reliable than using _UrlServingOn which might be stale or corrupted
  const int serverPort = InterSpecServer::portBeingServedOn();
  if( serverPort <= 0 )
  {
    NSLog(@"Cannot recover session: server not running (port=%d)", serverPort);
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Session Recovery Failed";
    alert.informativeText = [NSString stringWithFormat:
      @"Cannot recover session because server is not running.\n\n"
      "Server port: %d", serverPort];
    alert.alertStyle = NSAlertStyleCritical;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
    return;
  }
  
  // Construct the base URL from the port
  _UrlServingOn = [NSString stringWithFormat:@"http://127.0.0.1:%d", serverPort];
  _isServing = YES;
  NSLog(@"Using server URL: %@", _UrlServingOn);
  
  if( _sessionRecoveryAttempts > 3 )
  {
    NSLog(@"Too many recovery attempts (%ld), giving up", (long)_sessionRecoveryAttempts);
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Session Recovery Failed";
    alert.informativeText = @"Multiple recovery attempts failed. Please restart the application.";
    alert.alertStyle = NSAlertStyleCritical;
    [alert addButtonWithTitle:@"OK"];
    [alert runModal];
#endif
    return;
  }
  
  // Mark old session as dead if it exists
  NSString *oldToken = _UrlUniqueId;
  if( oldToken )
  {
    NSLog(@"Marking old session '%@' as destructing", oldToken);
    InterSpecServer::set_session_destructing( [oldToken UTF8String] );
  }
  
  // Generate new session token
  _UrlUniqueId = [self generateSessionToken];
  
  NSLog(@"Session recovery: old token='%@', new token='%@'", oldToken, _UrlUniqueId);
  
  // Register the new token as a primary app instance
  InterSpecServer::add_allowed_session_token(
    [_UrlUniqueId UTF8String],
    InterSpecServer::SessionType::PrimaryAppInstance
  );
  
  // Build the new URL with restore=yes to restore previous state
  // Extra safety check - should not happen due to guard at top of function
  if( !_UrlServingOn || !_UrlUniqueId )
  {
    NSLog(@"CRITICAL: recoverSession has nil UrlServingOn='%@' or UrlUniqueId='%@'",
          _UrlServingOn ?: @"(nil)", _UrlUniqueId ?: @"(nil)");
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Session Recovery Error";
    alert.informativeText = [NSString stringWithFormat:
      @"Unexpected nil value during recovery.\n\n"
      "UrlServingOn: %@\nUrlUniqueId: %@",
      _UrlServingOn ?: @"(nil)", _UrlUniqueId ?: @"(nil)"];
    alert.alertStyle = NSAlertStyleCritical;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
    return;
  }
  
  NSString *newURL = [NSString stringWithFormat:@"%@?apptoken=%@&primary=yes&restore=yes",
                      _UrlServingOn, _UrlUniqueId];
  
  // Check for dark mode
#ifdef AVAILABLE_MAC_OS_X_VERSION_10_14_AND_LATER
  if( @available(macOS 10.14, *) )
  {
    NSAppearanceName basicAppearance =
      [NSApp.mainWindow.effectiveAppearance
       bestMatchFromAppearancesWithNames:@[NSAppearanceNameAqua, NSAppearanceNameDarkAqua]];
    if( [basicAppearance isEqualToString:NSAppearanceNameDarkAqua] )
    {
      newURL = [NSString stringWithFormat:@"%@&colortheme=dark", newURL];
    }
  }
#endif
  
  NSLog(@"Loading recovered session URL: %@", newURL);
  NSLog(@"  _UrlServingOn: '%@'", _UrlServingOn);
  NSLog(@"  _UrlUniqueId: '%@'", _UrlUniqueId);
  
  // Create and validate the URL
  NSURL *url = [NSURL URLWithString:newURL];
  if( !url )
  {
    NSLog(@"ERROR: Failed to create NSURL from string: %@", newURL);
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Invalid Recovery URL";
    alert.informativeText = [NSString stringWithFormat:
      @"Failed to create URL from string.\n\n"
      "URL String: %@\n\nUrlServingOn: %@\nUrlUniqueId: %@",
      newURL, _UrlServingOn, _UrlUniqueId];
    alert.alertStyle = NSAlertStyleCritical;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
    return;
  }
  
  // Load the new session
  [_InterSpecWebView loadRequest:[NSURLRequest requestWithURL:url]];
  
  _sessionIsValid = NO; // Will be set to YES after successful load in didFinishNavigation
  
#if( PERFORM_DEVELOPER_CHECKS )
  NSAlert *alert = [[NSAlert alloc] init];
  alert.messageText = @"Session Recovery Attempted";
  alert.informativeText = [NSString stringWithFormat:
    @"Old token: %@\nNew token: %@\nRecovery attempt: %ld\n\n"
    "The session should restore automatically with previous state.",
    oldToken ?: @"(none)", _UrlUniqueId, (long)_sessionRecoveryAttempts];
  alert.alertStyle = NSAlertStyleInformational;
  [alert addButtonWithTitle:@"OK"];
  [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
}


- (void)periodicSessionHealthCheck
{
  if( !_isServing || !_UrlUniqueId )
    return;
  
  [self validateCurrentSession];
  
  if( !_sessionIsValid )
  {
    NSLog(@"Periodic health check detected invalid session - triggering recovery");
    
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Session Health Check Failed";
    alert.informativeText = [NSString stringWithFormat:
      @"Periodic health check detected the session '%@' is no longer valid. "
      "Attempting automatic recovery.", _UrlUniqueId];
    alert.alertStyle = NSAlertStyleWarning;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
    
    [self recoverSession];
  }
}


#pragma mark - Debug Simulation Methods

- (void)simulateWebViewTermination
{
  NSLog(@"Simulating WebView termination...");
  
#if( PERFORM_DEVELOPER_CHECKS )
  NSAlert *alert = [[NSAlert alloc] init];
  alert.messageText = @"Simulating WebView Termination";
  alert.informativeText = @"This will call the webViewWebContentProcessDidTerminate: handler "
                          "as if the system terminated the WebView process.";
  alert.alertStyle = NSAlertStyleInformational;
  [alert addButtonWithTitle:@"Continue"];
  [alert addButtonWithTitle:@"Cancel"];
  
  if( [alert runModal] == NSAlertFirstButtonReturn )
  {
    [self webViewWebContentProcessDidTerminate:_InterSpecWebView];
  }
#else
  [self webViewWebContentProcessDidTerminate:_InterSpecWebView];
#endif
}


- (void)simulateSessionInvalidation
{
  NSLog(@"Simulating session invalidation...");
  
  if( _UrlUniqueId )
  {
    // Mark the current session as dead
    InterSpecServer::set_session_destructing( [_UrlUniqueId UTF8String] );
    
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Session Invalidated";
    alert.informativeText = [NSString stringWithFormat:
      @"Session '%@' has been marked as dead.\n\n"
      "Use 'Force Session Validation' to check immediately and trigger recovery.\n\n"
      "Or start the periodic health check timer to detect this automatically.",
      _UrlUniqueId];
    alert.alertStyle = NSAlertStyleInformational;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
  }
  else
  {
    NSLog(@"simulateSessionInvalidation: No session token to invalidate");
  }
}


- (void)togglePeriodicHealthCheck
{
  if( _sessionHealthCheckTimer )
  {
    // Stop the timer
    [_sessionHealthCheckTimer invalidate];
    _sessionHealthCheckTimer = nil;
    NSLog(@"Stopped periodic session health check timer");
    
    // Update menu item title
    NSMenu *menu = [[[NSApp mainMenu] itemAtIndex: 1] submenu];
    for( NSMenuItem *item in [menu itemArray] )
    {
      if( [item action] == @selector(togglePeriodicHealthCheck) )
      {
        [item setTitle:@"Start Periodic Health Check (60s)"];
        break;
      }
    }
    
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Periodic Health Check Stopped";
    alert.informativeText = @"The 60-second periodic session health check timer has been stopped.";
    alert.alertStyle = NSAlertStyleInformational;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
  }
  else
  {
    // Start the timer
    _sessionHealthCheckTimer = [NSTimer scheduledTimerWithTimeInterval:60.0
                                                                target:self
                                                              selector:@selector(periodicSessionHealthCheck)
                                                              userInfo:nil
                                                               repeats:YES];
    NSLog(@"Started periodic session health check timer (60 second interval)");
    
    // Update menu item title
    NSMenu *menu = [[[NSApp mainMenu] itemAtIndex: 1] submenu];
    for( NSMenuItem *item in [menu itemArray] )
    {
      if( [item action] == @selector(togglePeriodicHealthCheck) )
      {
        [item setTitle:@"Stop Periodic Health Check"];
        break;
      }
    }
    
#if( PERFORM_DEVELOPER_CHECKS )
    NSAlert *alert = [[NSAlert alloc] init];
    alert.messageText = @"Periodic Health Check Started";
    alert.informativeText = @"The session will be checked every 60 seconds. "
                            "If the session becomes invalid, automatic recovery will be attempted.";
    alert.alertStyle = NSAlertStyleInformational;
    [alert addButtonWithTitle:@"OK"];
    [alert beginSheetModalForWindow:_window completionHandler:nil];
#endif
  }
}


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
