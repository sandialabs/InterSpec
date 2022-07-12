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
#include <boost/filesystem.hpp>  //toso: get rid of using boost in this file

//We gotta fix some wierd errors...
#ifdef check
#undef check
#endif
#ifdef require
#undef require
#endif

#include <Wt/WString>

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

#if( USE_SPECRUM_FILE_QUERY_WIDGET )
#include "InterSpec/SpecFileQueryWidget.h"
#endif

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SerialToDetectorModel.h"


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


-(BOOL)application:(NSApplication *)theApplication openFile:(NSString *)filename
{
  BOOL loaded = NO;
  
  if( !filename )
    return NO;
  
  NSLog( @"openFile" );
  [self setDbDirectory];
  
  NSLog( @"have setDbDirectory" );
  filename=[NSString stringWithFormat:@"\"%@\"",filename];
  
  std::string specfile = [filename UTF8String];
  if( specfile.length() && specfile[0] == '\"' )
    specfile = specfile.substr(1);
  if( specfile.length() && specfile[specfile.length()-1] == '\"' )
    specfile = specfile.substr(0,specfile.size()-1);

  NSLog( @"Will open specfile %s", specfile.c_str() );
  
  std::string name_ending;
  const std::string::size_type lastperiod = specfile.find_last_of( '.' );
  if( lastperiod != std::string::npos )
    name_ending = specfile.substr( lastperiod+1 );
  
  //There is a race condition here, the InterSpecApp might destruct before
  //  we call it
  InterSpecApp *specapp = 0;
  const std::set<InterSpecApp *> apps = InterSpecApp::runningInstances();
  
  if( !apps.empty() && !!_UrlUniqueId )
  {
    //for some reasons accessing the _UrlUniqueId was causing a crash bellow,
    //  so will instead cop out and just use the first app instance in the set,
    //  which really isnt the correct thing to do
//    const char *uniqueid = [_UrlUniqueId UTF8String];
//    const char *uniqueid = [_UrlUniqueId cStringUsingEncoding:NSASCIIStringEncoding];
//    if( uniqueid )
//      specapp = InterSpecApp::instanceFromExtenalIdString( uniqueid );
//    else
      specapp = *apps.begin();
    
    if( specapp )
    {
      NSLog( @"Will use existing session to open file" );
      
      Wt::WApplication::UpdateLock lock( specapp );
      
      if( lock )
      {
        loaded = specapp->userOpenFromFileSystem( specfile );
        specapp->triggerUpdate();
      }else
      {
        NSLog( @"Couldnt get WApplication::UpdateLock when opening a file" );
      }
    }//if( specapp )
    
    [[NSApplication sharedApplication] activateIgnoringOtherApps:YES];
  }else if( apps.empty() )
  {
    //The WebView may not have requested the URL yet when this function gets
    //  called, in the case a user double clicks on a file in the Finder to launch
    //  this app.  Therefore we will mark that we should open this file using the
    //  URL argument (that uses an intermediate database to store file location).
    NSLog( @"Will mark the file to open once a session is created" );
    
    const std::string dbpath = [[[self applicationFilesDirectory] path] UTF8String];
    DbToFilesystemLink::setFileNumToFilePathDBNameBasePath( dbpath );
    
    _fileNeedsOpening = DbToFilesystemLink::addFileToOpenToDatabase( specfile );
  }//if( instances.empty() )
  
  return loaded;
}//openFile


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
  _PreferenceDbPath = nil;
  _fileNeedsOpening = -1;
  _isServing = NO;
  _UrlServingOn = @"";
  _UrlUniqueId = nil;
}


- (void)setDbDirectory
{
  NSLog( @"setDbDirectory" );
  if( _PreferenceDbPath )
    return;
  
  NSLog( @"will get path" );
  boost::filesystem::path datadir = [[[self applicationFilesDirectory] path] UTF8String];
  
  if( !boost::filesystem::exists( datadir ) )
  {
    try
    {
      boost::filesystem::create_directories( datadir );
      NSLog( @"Created directory directory %s", datadir.c_str() );
    }catch(...){}
  }
  
  if( boost::filesystem::exists( datadir ) )
  {
    InterSpec::setWritableDataDirectory( datadir.string<std::string>() );
    
    const std::vector<std::string> serial_db = SpecUtils::ls_files_in_directory( datadir.string<std::string>(), "serial_to_model.csv" );
    if( !serial_db.empty() )
      SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
    
    
    
    datadir /= "InterSpecUserData.db";
   
    try
    {
      DataBaseUtils::setPreferenceDatabaseFile( datadir.string<std::string>() );
    }catch( std::exception &e )
    {
      NSLog( @"Error: %s", e.what() );
    }
    
    _PreferenceDbPath = [NSString stringWithFormat:@"%s", datadir.c_str()];
    NSLog( @"Datadir=%s", datadir.c_str() );
  }else
    NSLog( @"Failed to creade directory %s", datadir.c_str() );
}//- (void)setDbDirectory:(void);


- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
  NSLog(@"Finished Launching");

  
  Wt::WString::setDefaultEncoding( Wt::UTF8 );
  
  /*
   //Could maybe get rid of using XIB/NIB by manueally creating a window like:
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
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  [prefs setValue:@YES forKey:@"developerExtrasEnabled"];
  
  //Note: currently I disable right click in InterSpecApp using javascript if
  //      no PERFORM_DEVELOPER_CHECKS.  However, could also do:
  //https://stackoverflow.com/questions/28801032/how-can-the-context-menu-in-wkwebview-on-the-mac-be-modified-or-overridden#28981319
  //[_InterSpecWebView willOpenMenu:<#(nonnull NSMenu *)#> withEvent:<#(nonnull NSEvent *)#>];
#endif

  //To allow deep integration, could
  //[_webConfig setURLSchemeHandler:<#(nullable id<WKURLSchemeHandler>)#> forURLScheme: @"helloworld://"];
  
  //Create WKWebView manually, rather than in XIB to support macOS 10.10 and 10.11...
  self.InterSpecWebView = [[WKWebView alloc] initWithFrame: _window.contentView.frame configuration: webConfig];
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
  
  
  [self setDbDirectory];
  DataBaseVersionUpgrade::checkAndUpgradeVersion();

  
  static const std::string basedir = std::string("--basedir=")
                            + [[[NSBundle mainBundle] resourcePath] UTF8String];
  static const std::string argv0 = [[[NSBundle mainBundle] executablePath] UTF8String];
  
  NSString *tempDir = NSTemporaryDirectory();
  if( tempDir == nil ) //shouldnt ever fail, right
    tempDir = @"/tmp";
  static const std::string tmpdr = [tempDir UTF8String];  //static since I'm not sure how long the location pointed to by setenv has to last
  
  const char *argv[] = { argv0.c_str(), "--forceserve", "--nobrowsertab", basedir.c_str(), "-c", "data/config/wt_config_osx.xml", "--tempdir", tmpdr.c_str() };
  int argc = sizeof(argv) / sizeof(argv[0]);
  
  InterSpecServer::startServer( argc, (char **)argv, &createApplication );
  
  //now we'll wait for the server to start
  std::string url;
  int running = -1;
  int numtries = 0;
  while( url.empty() && running<0 && numtries < 110 )  //110 sections of at least 10 ms, is 1.1 seconds
  {
    url = InterSpecServer::urlBeingServedOn(); //will be empty if not serving
    if( url.empty() )
      std::this_thread::sleep_for( std::chrono::milliseconds(10) );
  }//while( numtries < 11 )
  
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
    
    const int randint = arc4random();
    _UrlUniqueId = [NSString stringWithFormat:@"%i", randint]; //@"123456789";
    
    InterSpecServer::add_allowed_session_token( [_UrlUniqueId UTF8String] );
    
    NSString *actualURL = [NSString stringWithFormat:@"%@?apptoken=%@&primary=yes", _UrlServingOn, _UrlUniqueId];
    
    if( _fileNeedsOpening > -1 )
    {
      actualURL = [NSString stringWithFormat:@"%@&specfile=%i", actualURL, _fileNeedsOpening ];
      _fileNeedsOpening = -1;
    }//if( fileNeedsOpening.size() )
    
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
    
    
    //if( [_InterSpecWebView respondsToSelector:@selector(mainFrame)])
    //[[_InterSpecWebView mainFrame] loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]]];
    [_InterSpecWebView loadRequest: [NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]] ];
    
    
    /*
     //Right now drag-n-drop from Outlook does not work - should investigate this:
     // Probably need to sub-class WKWebView, and then do something like 
     //https://developer.apple.com/library/archive/documentation/Cocoa/Conceptual/DragandDrop/Tasks/acceptingdrags.html#//apple_ref/doc/uid/20000993-BABHHIHC
    [_InterSpecWebView registerForDraggedTypes:[NSArray arrayWithObjects:
                      NSStringPboardType, NSFilenamesPboardType, NSTIFFPboardType,
                      NSRTFPboardType, NSTabularTextPboardType, NSFontPboardType,
                      NSRulerPboardType, NSColorPboardType, NSRTFDPboardType,
                      NSHTMLPboardType, NSURLPboardType, NSPDFPboardType,
                      NSMultipleTextSelectionPboardType, NSPostScriptPboardType,
                      NSVCardPboardType, NSInkTextPboardType, NSFilesPromisePboardType,
                      NSPasteboardTypeFindPanelSearchOptions, nil]];  // NSFilenamesPboardType is probably mostly what is needed

    - (NSDragOperation)draggingEntered:(id <NSDraggingInfo>)sender {
      NSPasteboard *pboard;
      NSDragOperation sourceDragMask;
      NSLog( @"In draggingEntered!" );
      sourceDragMask = [sender draggingSourceOperationMask];
      pboard = [sender draggingPasteboard];
      
      if ( [[pboard types] containsObject:NSColorPboardType] ) {
        if (sourceDragMask & NSDragOperationGeneric) {
          return NSDragOperationGeneric;
        }
      }
      if ( [[pboard types] containsObject:NSFilenamesPboardType] ) {
        if (sourceDragMask & NSDragOperationLink) {
          return NSDragOperationLink;
        } else if (sourceDragMask & NSDragOperationCopy) {
          return NSDragOperationCopy;
        }
      }
      return NSDragOperationNone;
    }
    */
    
    
    [_InterSpecWebView setNavigationDelegate: self];
    [_InterSpecWebView setUIDelegate: self];
    
    
    if (@available(macOS 10.14, *)) {
      [NSDistributedNotificationCenter.defaultCenter addObserver:self selector:@selector(themeChanged:) name:@"AppleInterfaceThemeChangedNotification" object: nil];
    }
  }
  
  // Try adding a popover, to eventually display a map
  //  See https://developer.apple.com/documentation/appkit/nspopover
  /*
   //  Need to implement https://developer.apple.com/documentation/appkit/nsviewcontroller?language=objc
   //  note need to override -loadView method
  NSRect bounds = [[self.window contentView] bounds];
  
  self.myPopover = [[NSPopover alloc] init];
  
  [self.myPopover setContentSize:NSMakeSize(0.85*bounds.size.width, 0.85*bounds.size.height)];
  [self.myPopover setBehavior: NSPopoverBehaviorTransient]; //change to NSPopoverBehaviorApplicationDefined
  [self.myPopover setAnimates:YES];
  [self.myPopover setContentViewController: self.popoverViewController];
  
  [self.myPopover showRelativeToRect: bounds
                            ofView:[[NSApp mainWindow] contentView]
                     preferredEdge:NSMinYEdge];
   */
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
        
        NSLog(@"host=%@, absurl=%@", host, absurl );
        
        if ([absurl rangeOfString:@"request=redirect&url=http"].location != NSNotFound
            || [[[request URL] scheme] isEqual:@"mailto"])
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
      
      //url = InterSpecServer::urlBeingServedOn(); //will be empty if not serving
      const NSInteger port = InterSpecServer::portBeingServedOn();
      NSString *url = [NSString stringWithFormat:@"http://localhost:%ld?primary=no&restore=no", port];
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

@end
