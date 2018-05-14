//
//  AppDelegate.m
//  InterSpec OSX
//
//  Created by Johnson, William C on 7/3/13.
//  Copyright (c) 2013 Johnson, William C. All rights reserved.
//

#import "AppDelegate.h"
#import <WebKit/WebView.h>
#import <WebKit/WebFrame.h>
#import <AppKit/NSSavePanel.h>

//A hack to speed up canvas rendering for smooth zoom out and pans
#include "target/osx/WebKit/WebPreferencesPrivate.h"

#include <set>
#include <string>
#include <stdlib.h>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

//We gotta fix some wierd errors...
#ifdef check
#undef check
#endif
#ifdef require
#undef require
#endif

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "InterSpec/DataBaseVersionUpgrade.h"
//#include "target/DesktopWithNativeBrowser/DesktopWithNativeBrowser.h"

#define SHOW_WEBVIEW_INSPECTOR 0

//#if( !BUILD_DESKTOP_STYLE_SERVER )
//#error You must enable the BUILD_DESKTOP_STYLE_SERVER preproccessor switch to make the OSX InterSpec app
//#endif
/*
//http://stackoverflow.com/questions/3344157/setting-default-application-for-given-file-extension-on-mac-os-x-from-code
 
#import <Foundation/Foundation.h>

@interface LaunchServicesWrapper : NSObject

+ (BOOL)setMyselfAsDefaultApplicationForFileExtension:
(NSString *)fileExtension;

@end


#import <ApplicationServices/ApplicationServices.h>
#import "LaunchServicesWrapper.h"

@implementation LaunchServicesWrapper

+ (NSString *)UTIforFileExtension:(NSString *)extension
{
  return (NSString *)CFBridgingRelease(
              UTTypeCreatePreferredIdentifierForTag( kUTTagClassFilenameExtension, (__bridge CFStringRef)extension, NULL ) );
}

+ (BOOL)setMyselfAsDefaultApplicationForFileExtension: (NSString *)fileExtension
{
  return LSSetDefaultRoleHandlerForContentType(
                                               (__bridge CFStringRef) [LaunchServicesWrapper
                                                                       UTIforFileExtension:fileExtension], kLSRolesAll,
                                               (__bridge CFStringRef) [[NSBundle mainBundle]
                                                                       bundleIdentifier]
                                               );
}

@end
*/


@implementation MyDownloadDelegate
- (void)download:(NSURLDownload *)download decideDestinationWithSuggestedFilename:(NSString *)filename
{
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

- (void)download:(NSURLDownload *)download didFailWithError:(NSError *)error
{
  NSLog(@"- (void)download:(NSURLDownload *)download didFailWithError:(NSError *)error");
  NSLog(@"Download failed! Error - %@ %@",
        [error localizedDescription],
        [[error userInfo] objectForKey:NSURLErrorFailingURLStringErrorKey]);
}

- (void)downloadDidFinish:(NSURLDownload *)download
{
  NSData *data = [[download request] HTTPBody];
  NSUInteger length = [data length];
//  long long expectedLength = [[download request] expectedContentLength];
  std::cout << "length=" << length << std::endl;
  NSLog(@"- (void)downloadDidFinish:(NSURLDownload *)download: ");
}

- (void)startDownloadingURL:sender
{
  NSLog(@"- (void)startDownloadingURL:sender");
}

- (void)downloadDidBegin:(NSURLDownload *)download
{
  NSLog(@"downloadDidBegin");

}


-(void)download:(NSURLDownload *)download didCreateDestination:(NSString *)path
{
  // path now contains the destination path
  // of the download, taking into account any
  // unique naming caused by -setDestination:allowOverwrite:
  NSLog(@"Final file destination: %@",path);
}
@end


@implementation InterSpecWebPolicyDecisionListener

- (NSArray *)webView:(WebView *)sender contextMenuItemsForElement:(NSDictionary *)element defaultMenuItems:(NSArray *)defaultMenuItems
{
#if( PERFORM_DEVELOPER_CHECKS )
  return defaultMenuItems;
#else
  //disable the right-click context menu
  return nil;
#endif
}

//Open dialog
- (void)webView:(WebView *)sender runOpenPanelForFileButtonWithResultListener:(id < WebOpenPanelResultListener >)resultListener
{
    // Create the File Open Dialog class.
    NSOpenPanel* openDlg = [NSOpenPanel openPanel];
    
    // Enable the selection of files in the dialog.
    [openDlg setCanChooseFiles:YES];
    
    // Enable the selection of directories in the dialog.
    [openDlg setCanChooseDirectories:NO];
    
    if ( [openDlg runModal] == NSOKButton )
    {
        NSArray* files = [[openDlg URLs]valueForKey:@"relativePath"];
        [resultListener chooseFilenames:files];
    }
    
}

// Probably won't call this function because everything should load in new window, so go to decidePolicyForNewWindowAction
-  (void)webView:(WebView *)webView decidePolicyForNavigationAction:(NSDictionary *)actionInformation
         request:(NSURLRequest *)request
           frame:(WebFrame *)frame
decisionListener:(id<WebPolicyDecisionListener>)listener
{
  [listener use];
}

//If link was clicked and opening in new window (URL/mailto links)
-  (void)webView:(WebView *)webView decidePolicyForNewWindowAction:(NSDictionary *)actionInformation
         request:(NSURLRequest *)request
    newFrameName:(NSString *)frameName
decisionListener:(id <WebPolicyDecisionListener>)listener
{
    if (WebNavigationTypeLinkClicked == [[actionInformation objectForKey:WebActionNavigationTypeKey] intValue])
    {
        // link was clicked and webview want to open it in new window do something with it...
        NSString *host = [[request URL] host];
        NSString *absurl = [[request URL] absoluteString];
        if ([absurl rangeOfString:@"request=redirect&url=http"].location != NSNotFound || [[[request URL] scheme] isEqual:@"mailto"]) {
            //external url or email
            [[NSWorkspace sharedWorkspace] openURL:[request URL]];
        } else if ([host isEqualToString:@"127.0.0.1"]) {
            //Email or URL open in native OSX client
            [listener download];
        }
        else {
            //Should not get here...
            [listener use];
        }
    }
    else
    {
        //Should not get here
        [listener ignore];
    }
}
@end

@implementation AppDelegate

@synthesize persistentStoreCoordinator = _persistentStoreCoordinator;
@synthesize managedObjectModel = _managedObjectModel;
@synthesize managedObjectContext = _managedObjectContext;


+ (void)initialize
{
#if( SHOW_WEBVIEW_INSPECTOR )
  [[NSUserDefaults standardUserDefaults] registerDefaults:@{@"WebKitDeveloperExtras": @YES,
                                                            @"WebKitScriptDebuggerEnabled": @YES,
                                                            @"WebKitScriptProfilerEnabled": @YES}];
#endif
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

#if( PERFORM_DEVELOPER_CHECKS )
  [[NSUserDefaults standardUserDefaults] setBool:TRUE forKey:@"WebKitDeveloperExtras"];
  [[NSUserDefaults standardUserDefaults] synchronize];
#endif
  
  //Lets adjust the size and position of the main window to be reasonable
  const int width = [[NSScreen mainScreen] frame].size.width;
  const int height = [[NSScreen mainScreen] frame].size.height;
  NSSize mySize;
  mySize.width = 0.85*width;
  mySize.height = 0.85*height;
  [_window setContentSize: mySize];
  [[self window] setFrameTopLeftPoint:NSMakePoint(0.025*width, 0.95*height)];
  
  [self setDbDirectory];
  DataBaseVersionUpgrade::checkAndUpgradeVersion();
  
  
  static const std::string basedir = std::string("--basedir=")
                            + [[[NSBundle mainBundle] resourcePath] UTF8String];
  static const std::string argv0 = [[[NSBundle mainBundle] executablePath] UTF8String];
  
  
  const char *argv[] = { argv0.c_str(), "--forceserve", "--nobrowsertab", basedir.c_str(), "-c", "data/config/wt_config_osx.xml" };
  int argc = sizeof(argv) / sizeof(argv[0]);
  
  InterSpecServer::startServer( argc, (char **)argv, &createApplication );
  
//  boost::function<void(void)> worker =
//          boost::bind( DesktopWithNativeBrowser::DesktopRunInDefaultBrowser,
//                      argc, (char **)argv, &createApplication );
//  boost::thread slave( worker );
//  slave.detach();
  
  //now well wait for the server to start
  std::string url;
  int running = -1;
  int numtries = 0;
  while( url.empty() && running<0 && numtries < 110 )  //110 sections of at least 10 ms, is 1.1 seconds
  {
    url = InterSpecServer::urlBeingServedOn(); //will be empty if not serving
    if( url.empty() )
      boost::this_thread::sleep( boost::posix_time::milliseconds(10) );
  }//while( numtries < 11 )
  
  if( url.empty() )
  {
    _isServing = NO;
     std::cerr << "Unable to start server!" << std::endl;
  }else
  {
    _isServing = YES;
    _UrlServingOn = [NSString stringWithUTF8String:url.c_str()];
    [_window setTitle: [NSString stringWithFormat:@"TRB InterSpec - %@", _UrlServingOn]];
    
    const int randint = arc4random();
    _UrlUniqueId = [NSString stringWithFormat:@"%i", randint]; //@"123456789";
    NSString *actualURL = [NSString stringWithFormat:@"%@?externalid=%@", _UrlServingOn, _UrlUniqueId];
    
    if( _fileNeedsOpening > -1 )
    {
      actualURL = [NSString stringWithFormat:@"%@&specfile=%i", actualURL, _fileNeedsOpening ];
      _fileNeedsOpening = -1;
    }//if( fileNeedsOpening.size() )
    
    //A hack to speed up canvas rendering for smooth zoom out and pans
    //  (Note, this is playing with private prefernces Apple doesnt want us to,
    //   so there is a user option way to do the bellow that should be good,
    //   but is untested)
    WebPreferences * prefs = [_InterSpecWebView preferences];
    
    //this is the only setting actually needed to speed things up
    if ([prefs respondsToSelector:@selector(setCanvasUsesAcceleratedDrawing:)]) {
      [prefs setCanvasUsesAcceleratedDrawing:YES];
    }
    
    if ([prefs respondsToSelector:@selector(setAcceleratedDrawingEnabled:)]) {
      [prefs setAcceleratedDrawingEnabled:YES];
    }
    
    if ([prefs respondsToSelector:@selector(setAcceleratedCompositingEnabled:)]) {
      [prefs setAcceleratedCompositingEnabled:YES];
    }
    
    if ([prefs respondsToSelector:@selector(setAccelerated2dCanvasEnabled:)]) {
      [prefs setAccelerated2dCanvasEnabled:YES];
    }
    
    if ([prefs respondsToSelector:@selector(setWebGLEnabled:)]) {
      [prefs setWebGLEnabled:YES];
    }
     
    
    /*
    //The bellow should work instead of the above
    [[NSUserDefaults standardUserDefaults] setObject:@YES
                                              forKey:@"WebKitWebGLEnabled"];
    [[NSUserDefaults standardUserDefaults] setObject:@YES
                                              forKey:@"WebKitAccelerated2dCanvasEnabled"];
    [[NSUserDefaults standardUserDefaults] setObject:@YES
                                              forKey:@"WebKitAcceleratedCompositingEnabled"];
    [[NSUserDefaults standardUserDefaults] setObject:@YES
                                              forKey:@"WebKitAcceleratedDrawingEnabled"];
    [[NSUserDefaults standardUserDefaults] setObject:@YES
                                              forKey:@"WebKitCanvasUsesAcceleratedDrawing"];
    */
    
    //Set the user agent string so
    NSString *userAgentStr = @"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_4) AppleWebKit/603.1.30 (KHTML, like Gecko) Version/10.1 Safari/603.1.30";
    [_InterSpecWebView setCustomUserAgent: userAgentStr];
    
    
    [[_InterSpecWebView mainFrame] loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]]];
    
    _UrlDownload = [[MyDownloadDelegate alloc] init];
    InterSpecWebPolicyDecisionListener* listener = [[InterSpecWebPolicyDecisionListener alloc] init];

#if( SHOW_WEBVIEW_INSPECTOR )
    [[_InterSpecWebView inspector] show:nil];
#endif
    
    //Download delegate
    [_InterSpecWebView setDownloadDelegate:_UrlDownload];
    //Link delegate
    [_InterSpecWebView setPolicyDelegate:listener];
    //Open dialog delegate
    [_InterSpecWebView setUIDelegate:listener];
  }
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

// Performs the save action for the application, which is to send the save: message to the application's managed object context. Any encountered errors are presented to the user.
- (IBAction)saveAction:(id)sender
{
  NSError *error = nil;
    
  if (![[self managedObjectContext] commitEditing])
    NSLog(@"%@:%@ unable to commit editing before saving", [self class], NSStringFromSelector(_cmd));
    
  if( ![[self managedObjectContext] save:&error] )
    [[NSApplication sharedApplication] presentError:error];
}

- (NSApplicationTerminateReply)applicationShouldTerminate:(NSApplication *)sender
{
  // Save changes in the application's managed object context before the application terminates.
  
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
        
    if( answer == NSAlertAlternateReturn )
      return NSTerminateCancel;
  }

  InterSpecServer::killServer();
  return NSTerminateNow;
}

- (IBAction)ShowInBrowserClicked:(id)sender
{
  std::string url = ([_UrlServingOn UTF8String]);
  system( ("open "+ url).c_str() );
  //DesktopWithNativeBrowser::openUrl( [_UrlServingOn UTF8String] );
}



- (void)webView:(WebView *)webView
  decidePolicyForNewWindowAction:(NSDictionary *)actionInformation
                         request:(NSURLRequest *)request
                    newFrameName:(NSString *)frameName
                decisionListener:(id < WebPolicyDecisionListener >)listener
{
  //This function gets called when the user "Downloads" a spectrum, or a CSV
  //  file (e.g. when a new window is tried to be openend).
  
  //The next line causes MyDownloadDelegate decideDestinationWithSuggestedFilename
  //  to be called
  [listener download];
}//decidePolicyForNewWindowAction


@end
