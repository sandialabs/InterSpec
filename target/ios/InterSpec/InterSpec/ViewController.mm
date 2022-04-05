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
//
//  Created by Johnson, William C on 9/30/13.
//

#import <Foundation/NSData.h>

#import <WebKit/WebKit.h>

#import "ViewController.h"

#include <stdlib.h>


#include <string>
#include <fstream>
#include <iostream>

#include <Wt/WServer>

#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "AppDelegate.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/DbToFilesystemLink.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

#define IS_IPAD() (UI_USER_INTERFACE_IDIOM() == UIUserInterfaceIdiomPad)

@interface ViewController () <UIWebViewDelegate,UIDocumentInteractionControllerDelegate>

@end


#pragma mark - 


@implementation ViewController


Wt::WApplication *createApplication(const Wt::WEnvironment& env)
{
  return new InterSpecApp( env );
}



-(void) fixPermissions:(NSString *)path
{
  //Starting with iOS 10.3 (APFS) default file creation mode is '4', which is read only, so any files we create 
  NSFileManager *filemgr;
  NSDictionary *attribs;
  
  filemgr = [NSFileManager defaultManager];
  
  attribs = [filemgr attributesOfItemAtPath: path error: NULL];
  
  NSLog (@"File type %@", [attribs objectForKey: NSFileType]);
  NSLog (@"POSIX Permissions %@", [attribs objectForKey: NSFilePosixPermissions]);
  
  
  if (![filemgr fileExistsAtPath:path]) {
    NSLog(@"Not setting permissions for %@; does not exist", path);
    return;
  }   
  
  
  NSError *error = nil;
  // Get the current permissions
  NSDictionary *currentPerms = [filemgr attributesOfFileSystemForPath:path error:&error];
  if (currentPerms) {
    // Update the permissions with the new permission
    NSMutableDictionary *attributes = [currentPerms mutableCopy];
    attributes[NSFilePosixPermissions] = @(0666);
    if (![filemgr setAttributes:attributes ofItemAtPath: path error:&error]) {
      NSLog(@"Unable to make %@ read-only: %@", path, error);
    }
  } else {
    NSLog(@"Unable to read permissions for %@: %@", path, error);
  }
  
}

-(NSString *)generateSessionToken
{
#define TOKEN_LEN 14
  char data[TOKEN_LEN];
  for( int index = 0; index < TOKEN_LEN; ++index )
    data[index] = (char)('A' + (arc4random_uniform(26)));
  
  return [[NSString alloc] initWithBytes:data length:TOKEN_LEN encoding:NSUTF8StringEncoding];
}//NSString *generateSessionToken



-(BOOL)startServer
{
  //XXX - could look at if we could call Wt::WServer::resume() instead, when returning fmor the background
  if( _isServing )
    return NO;
  
  //We set temp dir in +initialize, but I guise I dont really know if it can
  //  change in iOS every once in a while, so do it here to, JIC
  NSString *tempDir = NSTemporaryDirectory();
  if( tempDir == nil ) //shouldnt ever fail, right
    tempDir = @"/tmp";
  const std::string tmpdirstr = [tempDir UTF8String];
  
  
  static std::string argv0 = "--docroot";
  static std::string argv1 = ".";
  static std::string argv4 = "-c";
  static std::string argv5 = "./data/config/wt_config_ios.xml";
  static std::string argv6 = "--tempdir";
  static std::string argv7 = [tempDir UTF8String];
  static char *argv[] = { &argv0[0], &argv1[0], &argv4[0], &argv5[0], &argv6[0], &argv7[0] };
  int argc = sizeof(argv) / sizeof(argv[0]);
  InterSpecServer::startServer( argc, argv, &createApplication );
  
  _isServing = YES;
  const std::string url = InterSpecServer::urlBeingServedOn();
  _UrlServingOn = [NSString stringWithUTF8String:url.c_str()];
  
  return YES;
}//-(BOOL)startServer


-(BOOL)enteredBackground
{
  //Use this method to release shared resources, save user data, invalidate
  //  timers, and store enough application state information to restore your
  //  application to its current state in case it is terminated later.
  //If your application supports background execution, this method is called
  //  instead of applicationWillTerminate: when the user quits.

  //For some reason the WServer stops responding after the app is placed in the
  //  background, and then brought back into foreground.  So we'll instead kill
  //  off the WServer (which will save user state if option is enabled).  And
  //  then upon re-entering the foreground, we'll load the saved user state.
  
  //We have about 5 seconds in this function before the OS kills the app
  //  - should consider using a beginBackgroundTaskWithExpirationHandler call
  
  NSLog( @"In enteredBackground" );
  
  void (^killStuffBlock)(void) = ^{
    NSLog( @"In enteredBackground ---> killStuffBlock" );
    _appHasGoneIntoBackground = YES;
    _appComminFromBackground = YES;
    
    InterSpecServer::killServer();
    NSLog( @"Killed server" );
    
    _isServing = NO;
    _UrlServingOn = @"";

  
    if( _UrlUniqueId )
    {
      const int status = InterSpecServer::remove_allowed_session_token( [_UrlUniqueId UTF8String] );
      if( status != 0 )
      {
        NSLog( @"InterSpecServer::remove_allowed_session_token gave unexpected status %ld for token '%@'.", _UrlUniqueId );
      }
    }
    _UrlUniqueId = @"";
    
    //We could actually delete the webview here and recreate when bringing to
    //  the foreground, but whatever.  For now lets just clear the webview to
    //  save a couple megabytes of ram.
    [_webView loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:@"about:blank"]]];
    
    //The following doesnt seem to make a big difference, so lets not bother with it
    //see: http://twobitlabs.com/2012/01/ios-ipad-iphone-nsurlcache-uiwebview-memory-utilization/
    [[NSURLCache sharedURLCache] removeAllCachedResponses];
  };
  
  if( [NSThread isMainThread] )
  {
    killStuffBlock();
  }else
  {
    dispatch_group_t group = dispatch_group_create();
    dispatch_group_async(group, dispatch_get_main_queue(), killStuffBlock );
    dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
    NSLog(@"Syncing up after killStuffBlock" );
  }
  
  NSLog( @"Done in enteredBackground" );
  
  return YES;
}//-(BOOL)enteredBackground


-(BOOL)willEnterForeground
{
  //Called as part of the transition from the background to the inactive state;
  //  here you can undo many of the changes made on entering the background.
  //This function is called _before_ openSpectrumFile
  
  NSLog( @"In willEnterForeground" );
  
//  [self startServer];
  //XXX - could look at if we could call Wt::WServer::resume() instead
  
  return YES;
}//-(BOOL)willEnterForeground


-(BOOL)wakeupFromBackground
{
  //Restart any tasks that were paused (or not yet started) while the
  //application was inactive. If the application was previously in the background,
  //  optionally refresh the user interface.
  //This function is called _after_ openSpectrumFile
  
  NSLog( @"In wakeupFromBackground" );
  
  boost::filesystem::path apppath = [[[NSBundle mainBundle] executablePath] UTF8String];
  apppath = apppath.parent_path();
  
//  NSLog( @"Setting CWD=%s", apppath.string<std::string>().c_str() );
//  boost::filesystem::current_path( apppath );

  
  if( !_appComminFromBackground )
  {
    NSLog( @"App never went in background, skipping from wakeupFromBackground" );
    if( _dbIndexOfFileToOpen < 0 )
      return YES;

    if( _dbIndexOfFileToOpen >= 0 )
    {
      NSLog( @"We should open up a file though" );
      InterSpecApp *app = InterSpecApp::instanceFromExtenalToken( [_UrlUniqueId UTF8String] );

      if( app )
      {
        Wt::WApplication::UpdateLock applock( app );
        Wt::WApplication *wtapp = Wt::WApplication::instance();
        
        if( applock && wtapp )
        {
          const bool success = app->openFileFromDbFileSystemLink( _dbIndexOfFileToOpen );
          app->triggerUpdate();
          NSLog( @"Success opening file = %i", int(success) );
          return YES;
        }//if( applock && wtapp )
      }//if( app )
    }//if( _dbIndexOfFileToOpen >= 0 )
    
    //If were here, we failed to get the WApplication or something, so lets kill
    //  off the server (to save any active sessions), and then restore the state
    [self enteredBackground];
  }//if( !_appComminFromBackground )
  
  if( !_isServing )
    [self startServer];
  
  _appComminFromBackground = NO;
  
  
  _UrlUniqueId = [this generateSessionToken];
  InterSpecServer::add_allowed_session_token( [_UrlUniqueId UTF8String], InterSpecServer::SessionType::PrimaryAppInstance );
  
  NSString *actualURL = [NSString stringWithFormat:@"%@?apptoken=%@", _UrlServingOn, _UrlUniqueId];
  
  //Check to see if we should open a spectrum file; if so append which on to URL
  if( _dbIndexOfFileToOpen >= 0 )
    actualURL = [NSString stringWithFormat:@"%@&specfile=%i", actualURL, _dbIndexOfFileToOpen];
  _dbIndexOfFileToOpen = -1;
  
  //See if we should specify orientation and safe area in the URL, so it will be
  //  known before the html/JS gets loaded.  We could do the same thing for initial
  //  positioning/sizing of things as well...
  if( @available(iOS 11, *) )
  {
    UIEdgeInsets insets = [[[UIApplication sharedApplication] delegate] window].safeAreaInsets;
    UIInterfaceOrientation orientation = [[UIApplication sharedApplication] statusBarOrientation];
    
    NSLog( @"wakeupFromBackground Orientation=%i, SafeAreas = {t=%f, r=%f, b=%f, l=%f}",
          int(orientation), insets.top, insets.right, insets.bottom, insets.right );
    
    if( insets.top > 0 || insets.right > 0 || insets.bottom > 0 || insets.left > 0 )
    {
      actualURL = [NSString stringWithFormat:@"%@&SafeAreas=%i,%i,%i,%i,%i", actualURL,
                   int(orientation), int(insets.top), int(insets.right), int(insets.bottom), int(insets.left)];
    }
  }//if( @available(iOS 11, *) )
  
  bool restore = false;
  if( [[NSUserDefaults standardUserDefaults] boolForKey:@"ShouldRestore"] )
    restore = true;
  
  //make sure next time the app starts up, we default to not trying to
  //  restore (see note in shouldRestoreApplicationState about possible moving
  //  not using UserDefaults)
  [[NSUserDefaults standardUserDefaults] setBool:NO forKey:@"ShouldRestore"];
  [[NSUserDefaults standardUserDefaults] synchronize];
  
  //If the user has killed the app, or the last session crashed, then restore
  //  will be false.
  if( !restore && !_appHasGoneIntoBackground )
  {
    NSLog(@"Not attempting to restore state");
    actualURL = [NSString stringWithFormat:@"%@&restore=0", actualURL ];
  }else
  {
    NSLog(@"Will allow state restoration");
  }
  
  NSLog(@"\n\nwillEnterForeground: started server at URL=%@", actualURL);
  
  [_webView loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]]];
  
  return YES;
}//-(BOOL)wakeupFromBackground




-(BOOL)openSpectrumFile:(NSURL *)url
{
  std::string urlstr = [[url path] UTF8String];
  //std::string urlstr = [url fileSystemRepresentation];
  
  // I dont think we need this try /catch any more, but leaving in, because well, why not
  try
  {
    NSLog( @"Will check to try and open %s", urlstr.c_str() );
    
    NSFileManager *filemgr = [NSFileManager defaultManager];
    if( ![filemgr isReadableFileAtPath: [url path]] ){
      NSLog(@"Can NOT access file at path '%@'.", [url path]);
      return NO;
    }
    
    NSLog( @"Can access, and will try to open '%s'", urlstr.c_str() );
    
    std::string path_json = "[\"" + urlstr + "\"]";
    int status = InterSpecServer::open_file_in_session( [_UrlUniqueId UTF8String], path_json.c_str() );
    if( status > 0 ){
      NSLog( @"Opened file at '%@' directly via InterSpecServer.", [url path] );
      return YES;
    }
    
    NSLog( @"Will try to open file via DbToFilesystemLink from '%@'.", [url path] );
  
    // We will copy the file to a temporary location, since by the time we try to read it into
    //  InterSpec, we may not have permissions to access the file an further.
    
    NSString *tempDir = NSTemporaryDirectory();
    NSString *fname = [url lastPathComponent];
    // TODO: should create unique directory so there wont be a filename clash
    NSString *tmpfile = [tempDir stringByAppendingPathComponent:fname];
    
    NSError *error;
    BOOL copied = [filemgr copyItemAtPath: [url path] toPath: tmpfile error:&error];
    if( copied ){
      NSLog( @"Copied file from '%@' t0 '%@'.", [url path], tmpfile );
      urlstr = [tmpfile UTF8String];
    }else{
      NSLog( @"Failed to copy file from '%@' t0 '%@', with error %@.", [url path], tmpfile, error );
    }
    
    // TODO: Try parsing the file here, and see if it is a valid spectrum file, so this way we can return a proper YES/NO.  Also, we could just hold it in memorry and pass it off to the InterSpec session once it starts (but would need to implement this in InterSpec class)
    // TODO: cleanup temporary file if we copied it over
    // TODO: skip this whole database thing, and just have a member variable to hold the URL, and open it once the app loads
    
    //TODO - could try to see if file is a background file to the current
    //       spectrum, and if so load it as the background spectrum `
    DbToFilesystemLink::FileIdToLocation requestinfo;
    requestinfo.m_foregroundFilePath = urlstr;
    //    requestinfo.m_backgroundFilePath = "/path/to/background";
    const int dbid = DbToFilesystemLink::addFileToOpenToDatabase( requestinfo );
    if( dbid < 0 )
    {
      NSLog( @"openSpectrumFile: Error saving file path to database" );
      return NO;
    }
      
    _dbIndexOfFileToOpen = dbid;
  }catch( std::exception &e )
  {
    NSLog( @"Caught exception trying to access file: %s", e.what() );
    return NO;
  }

  return YES;
}//openSpectrumFile




- (void)viewDidLoad
{
  NSLog( @"Starting viewDidLoad" );
  
  _isServing = NO;
  _UrlUniqueId = @"";
  _UrlServingOn = @"";
  _dbIndexOfFileToOpen = -1;
  _appHasGoneIntoBackground = NO;
  _appComminFromBackground = YES;
  
  //self.edgesForExtendedLayout = UIRectEdgeNone;
  [super viewDidLoad];
  
  //For some reason we still dont fill the space on iPhoneXs... maybe we need
  //  to do something in JavaScript... maybe try css-tricks.com/the-notch-and-css and then probably convert to WkWebView...
  WKWebViewConfiguration *wvconfig = [[WKWebViewConfiguration alloc] init];
  //[wvconfig setApplicationNameForUserAgent:<#(NSString * _Nullable)#>]
  _webView = [[WKWebView alloc] initWithFrame: CGRectZero configuration: wvconfig];
  [self.view addSubview: _webView];
  
  [_webView setTranslatesAutoresizingMaskIntoConstraints: NO];
  [[[_webView leadingAnchor] constraintEqualToAnchor: self.view.leadingAnchor constant: 0.0] setActive: YES];
  [[[_webView trailingAnchor] constraintEqualToAnchor: self.view.trailingAnchor constant: 0.0]setActive: YES];
  [[[_webView topAnchor] constraintEqualToAnchor: self.view.topAnchor constant: 0.0] setActive: YES];
  [[[_webView bottomAnchor] constraintEqualToAnchor: self.view.bottomAnchor constant: 0.0] setActive: YES];

  UIScrollView *view = [_webView scrollView];
  [view setScrollEnabled: NO];
  
  if (@available(iOS 11.0, *)) {
    view.contentInsetAdjustmentBehavior = UIScrollViewContentInsetAdjustmentNever;
  }
  
  self.webView.UIDelegate = self;
  [_webView setAllowsLinkPreview: NO];
  //[_webView setAllowsMagnification: NO]; //default is not
  [_webView setAllowsBackForwardNavigationGestures: NO];
  
  if( UI_USER_INTERFACE_IDIOM() == UIUserInterfaceIdiomPad )
  {
    //[_webView evaluateJavaScript:@"navigator.userAgent" completionHandler:^(id __nullable userAgent, NSError * __nullable error) {
    //  NSLog(@"The actual user agent: %@", userAgent);
    //  On my iPad prints "Mozilla/5.0 (iPad; CPU OS 14_4 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Mobile/15E148"
    //}];
    
    // For some reason I get 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_6) AppleWebKit/605.1.15 (KHTML, like Gecko)' in the C++ on my iPad,
    //  which means the c++ of InterSpec doesnt detect it as an iPad - so we'll set a custom user agent that includes "iPad"
    [_webView setCustomUserAgent: @"Mozilla/5.0 (iPad; CPU OS 12_2 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Mobile/15E148"];
  }else if( UI_USER_INTERFACE_IDIOM() == UIUserInterfaceIdiomPhone )
  {
    // This looks to always be fine, so we wont set a custom user agent
  }
  
  [_webView setNavigationDelegate: self];
  
  NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
  NSString *documentsDirectory = [paths objectAtIndex:0];
  
  //The database needs the sqlite ending for some reason.  If I use the name
  //  SRBPreferences.db, then it doesnt work
  NSString *writableDBPath = [documentsDirectory stringByAppendingPathComponent:@"data.sqlite"];
  DataBaseUtils::setPreferenceDatabaseFile( [writableDBPath UTF8String] );
  InterSpec::setWritableDataDirectory( [documentsDirectory UTF8String] );
  
  NSLog( @"Writing Database to '%@'", writableDBPath );


  [self fixPermissions: writableDBPath];
  
  DataBaseVersionUpgrade::checkAndUpgradeVersion();

  [self fixPermissions: writableDBPath];
  
   NSString *fileNumToFilePathToDBPath = [documentsDirectory stringByAppendingPathComponent:@"FileNumToFilePath.sqlite"];
  [self fixPermissions: fileNumToFilePathToDBPath];
  
  const boost::filesystem::path dbpath = [documentsDirectory UTF8String];
  const bool validdir = DbToFilesystemLink::setFileNumToFilePathDBNameBasePath( dbpath.string<std::string>() );
  if( !validdir )
    NSLog( @"viewDidLoad: Error setting FileNumToFilePathDBNameBasePath: %@", documentsDirectory );
  
  [self fixPermissions: fileNumToFilePathToDBPath];
  
  NSLog( @"Done in viewDidLoad" );
}


- (void)onKeyboardHide
{
  //This function is vestigiual from debugging an issue with the keyboard hiding
  //  NSLog(@"Triggering screen resize" );
  //  NSString *jsstring = [NSString stringWithFormat:@"setTimeout( function(){$('.Wt-domRoot').height(window.innerHeight); window.scrollTo(0,0);}, 0 );" ];
  //  [_webView stringByEvaluatingJavaScriptFromString: jsstring];
}

- (void)webView:(WKWebView *)webView didFailNavigation:(WKNavigation *)navigation withError:(NSError *)error
{
  NSLog(@"didFailNavigation" );
}//didFailNavigation

- (void)webView:(WKWebView *)webView didCommitNavigation:(WKNavigation *)navigation
{
  NSLog(@"didCommitNavigation" );
}//didCommitNavigation

- (void)webView:(WKWebView *)webView didFinishNavigation:(WKNavigation *)navigation
{
  NSLog(@"didFinishNavigation" );
  
  //NSLog( @"Doing horrible hack for iOS 9, to force a re-size of the contents." ) ;
  //[_webView evaluateJavaScript:@"setTimeout( function(){ $('.Wt-domRoot').get(0).wtResize(); }, 500 );" completionHandler: nil];
  
  
  //[_webView evaluateJavaScript:@"document.body.style.background = \"background: red;\"" completionHandler: nil];
  //See https://stackoverflow.com/questions/4629969/ios-return-bad-value-for-window-innerheight-width
  //  for description of issue
  //[_webView evaluateJavaScript:@"window.innerHeight" completionHandler:
  // ^(id _Nullable result, NSError * _Nullable error) {
  //   NSString *heightString = nil;
  //   if (error == nil) {
  //     if (result != nil) {
  //       heightString = [NSString stringWithFormat:@"%@", result];
  //     }
  //   } else {
  //     NSLog(@"evaluateJavaScript error : %@", error.localizedDescription);
  //   }
  //   if( heightString )
  //     NSLog( @"innerHeight=%@\n", heightString ) ;  //window.innerHeight=552.000000, window.innerWidth=980
  // }];
}//didFinishNavigation


//Opens URL in Safari, and email links in native Mail
- (void)webView:(WKWebView *)webView decidePolicyForNavigationAction:(WKNavigationAction *)navigationAction
                decisionHandler:(void (^)(WKNavigationActionPolicy))decisionHandler
{
  NSLog( @"decidePolicyForNavigationAction" );
  
  if( navigationAction == nil )
  {
    NSLog( @"nil navigationAction" );
    decisionHandler( WKNavigationActionPolicyCancel );
    return;
  }
  
  if( navigationAction.navigationType == WKNavigationTypeLinkActivated )
  {
    //CSVs from decay widget get here
    NSLog( @"WKNavigationTypeLinkActivated" );
    
    //To open the URL in safari, you can do:
    //UIApplication *application = [UIApplication sharedApplication];
    //[application openURL:navigationAction.request.URL options:@{} completionHandler:nil];
    
    //However, to try and save data to a file, we will download the data, and
    //  let the user save it
    void (^completionHandlerBlock)(NSData *, NSURLResponse *, NSError *)
    = ^(NSData *data, NSURLResponse *response, NSError *error){
       NSLog( @"completionHandlerBlock" );
      
      if( error || (data==nil) )
      {
        NSLog( @"completionHandlerBlock Error" );
        return;
      }
      
      NSString *name = nil, *mimetype = nil;
      if( response )
      {
        mimetype = [response MIMEType];
        name = [response suggestedFilename];
      }
      
      if( name == nil )
        name = @"download.txt";
      
      NSString *tempDir = NSTemporaryDirectory();
      if( tempDir == nil ) //shouldnt ever fail, right
        tempDir = @"/tmp";
      NSString *tmpfile = [tempDir stringByAppendingPathComponent:name];
      if( [data writeToFile: tmpfile  atomically: YES] )
      {
        dispatch_async( dispatch_get_main_queue(), ^{
          AppDelegate *appDelegate = (AppDelegate *)[[UIApplication sharedApplication] delegate];
          [appDelegate sendSpectrumFileToOtherApp: tmpfile];
        } );
      }else
      {
        NSLog( @"Failed to write file '%@'", tmpfile );
      }
    };//completionHandlerBlock
    
    NSURLSession *session = [NSURLSession sharedSession];
    NSURLSessionDataTask *task = [session dataTaskWithURL:navigationAction.request.URL completionHandler:completionHandlerBlock];
    [task resume];
  
    
    decisionHandler( WKNavigationActionPolicyCancel );
    return;
  }//
  
  if( navigationAction.navigationType == WKNavigationTypeBackForward )
  {
    NSLog( @"WKNavigationTypeBackForward" );
    decisionHandler( WKNavigationActionPolicyCancel );
    return;
  }
  
  if( navigationAction.navigationType == WKNavigationTypeReload )
    NSLog( @"WKNavigationTypeReload" );
  
  if( navigationAction.navigationType == WKNavigationTypeFormResubmitted )
    NSLog( @"WKNavigationTypeFormResubmitted" );
  
  if( navigationAction.navigationType == WKNavigationTypeFormResubmitted )
    NSLog( @"WKNavigationTypeFormSubmitted" );
  
  if( navigationAction.navigationType == WKNavigationTypeOther )
  {
    //Get here for the initial load, or opening AuxWindow
    NSLog( @"WKNavigationTypeOther" );
  }
  
  
  decisionHandler( WKNavigationActionPolicyAllow );
}//decidePolicyForNavigationAction


/*
- (void)webView:(UIWebView *)webView didFailLoadWithError:(NSError *)error
{
  NSLog(@"\n\n\ndidFailLoadWithError: error => %@ ", [error userInfo] );
}

- (void)webViewDidFinishLoad:(UIWebView *)webView
{
  NSLog(@"WebViewDidFinishLoad" );
  
  UIWindow *topWindow = [[UIApplication sharedApplication] keyWindow];
  
  [topWindow setInsetsLayoutMarginsFromSafeArea: NO];
  [webView setInsetsLayoutMarginsFromSafeArea: NO];
  
  UILayoutGuide *marginGuid = [topWindow layoutMarginsGuide];
  CGRect marginrect = [marginGuid layoutFrame];
  NSLog( @"marginrect = {%f, %f}", marginrect.size.width, marginrect.size.height ); //552,304
  
  UILayoutGuide *layoutGuid = [topWindow safeAreaLayoutGuide];
  CGRect layoutrect = [layoutGuid layoutFrame];
  NSLog( @"layoutrect = {%f, %f}", layoutrect.size.width, layoutrect.size.height ); //568x320
  
  UIEdgeInsets insets = [topWindow safeAreaInsets];
  NSLog( @"TopInsets = {%f, %f, %f, %f}", insets.right, insets.left, insets.top, insets.bottom );
  
  CGRect webViewRect = [webView bounds];
  NSLog( @"webViewRect = {%f, %f, %f, %f}", webViewRect.origin.x, webViewRect.origin.y, webViewRect.size.width, webViewRect.size.height );
  
  UIEdgeInsets webVieEdge = [webView layoutMargins];
  NSLog( @"webVieEdge = {%f, %f, %f, %f}", webVieEdge.left, webVieEdge.right, webVieEdge.bottom, webVieEdge.top );
  
  //For some reason on iOS 9, there is an issue where Wt thinks the window height
  //  is the same as the window width (some talk on internet of something
  //  asyncrounous), so here we will do some JavaScript to force a resize of height
  //  to the proper height.  Note that window.innerHeight does give correct value here.
  NSLog( @"Doing horrible hack for iOS 9, to force a re-size of the contents." ) ;
  NSString * heightString = [webView stringByEvaluatingJavaScriptFromString:@"window.innerHeight"];
  CGFloat contentHeight = [ heightString doubleValue ] ;
  NSLog( @"contentHeight=%f\n", contentHeight ) ;
  
//  NSString *jsstring = [NSString stringWithFormat:@"setTimeout( function(){$('.Wt-domRoot').height(window.innerHeight); window.scrollTo(0,0);}, 1000 );" ];
//  [webView stringByEvaluatingJavaScriptFromString: jsstring];
  
//  webView.scrollView.scrollEnabled = NO;    // Property available in iOS 5.0 and later 
//  CGRect frame = webView.frame;
//  frame.size.height = 200; //webView.scrollView.contentSize.height;
//  webView.frame = frame;
//  CGSize fittingSize = [webView sizeThatFits:CGSizeZero];
//  frame.size = fittingSize;
//  webView.frame = frame;
//  NSLog(@"size: %f, %f", frame.size.width, frame.size.height);
}

- (void)webViewDidStartLoad:(UIWebView *)webView
{
  NSLog(@"\n\n\nwebViewDidStartLoad" );
}

//Opens URL in Safari, and email links in native Mail
-(BOOL) webView:(UIWebView *)inWeb shouldStartLoadWithRequest:(NSURLRequest *)inRequest navigationType:(UIWebViewNavigationType)inType {
    if ( inType == UIWebViewNavigationTypeLinkClicked ) {
        [[UIApplication sharedApplication] openURL:[inRequest URL]];
        return NO;
    }
    
    return YES;
}
 */

//- (void)webView:(WKWebView *)webView
//        runOpenPanelWithParameters:(WKOpenPanelParameters *)parameters
//        initiatedByFrame:(WKFrameInfo *)frame
//        completionHandler:(void (^)(NSArray<NSURL *> *URLs))completionHandler
//{
//}//

- (BOOL)webView:(WKWebView *)webView shouldPreviewElement:(WKPreviewElementInfo *)elementInfo
{
  return NO;
}

- (void)webViewDidClose:(WKWebView *)webView
{
  NSLog(@"webViewDidClose!" );
}


- (void)didReceiveMemoryWarning
{
  using namespace std;
  NSLog(@"didReceiveMemoryWarning!" );
  //  vector<const InterSpecApp *> viewers = InterSpecServer::runningInstances();
  
  [super didReceiveMemoryWarning];
  // Dispose of any resources that can be recreated.
}

//===================================================================
- (UIViewController *)documentInteractionControllerViewControllerForPreview:(UIDocumentInteractionController *)controller
{
  //required to implement UIDocumentInteractionControllerDelegate
  NSLog(@"documentInteractionControllerViewControllerForPreview" );
  return self;
}

- (UIView *)documentInteractionControllerViewForPreview:(UIDocumentInteractionController *)controller
{
  //required to implement UIDocumentInteractionControllerDelegate
  NSLog(@"documentInteractionControllerViewForPreview" );
  return self.view;
}

- (CGRect)documentInteractionControllerRectForPreview:(UIDocumentInteractionController *)controller
{
  //required to implement UIDocumentInteractionControllerDelegate
  NSLog(@"documentInteractionControllerRectForPreview" );
  return self.view.frame;
}

/*
 
#pragma mark - Delegate Methods

- (void)documentInteractionController:(UIDocumentInteractionController *)controller willBeginSendingToApplication:(NSString *)application {
  //Doesnt seem to make it here when sending files to other apps
  NSLog(@"willBeginSendingToApplication" );
}

- (void)documentInteractionController:(UIDocumentInteractionController *)controller didEndSendingToApplication:(NSString *)application{
  //Doesnt seem to make it here when sending files to other apps
  NSLog(@"didEndSendingToApplication" );
}

- (void)documentInteractionControllerDidDismissOpenInMenu:(UIDocumentInteractionController *)controller {
  //Doesnt seem to make it here when sending files to other apps
  NSLog(@"documentInteractionControllerDidDismissOpenInMenu" );
}
*/

- (void)viewSafeAreaInsetsDidChange
{
  NSLog(@"viewSafeAreaInsetsDidChange" );
  
  if( @available(iOS 11, *) )
  {
    [super viewSafeAreaInsetsDidChange];
    [self setSafeAreasToClient];
  }//if( @available(iOS 11, *) )
}

- (void)setSafeAreasToClient
{
  if( @available(iOS 11.0, *) )
  {
    UIEdgeInsets insets = [[[UIApplication sharedApplication] delegate] window].safeAreaInsets;
    UIInterfaceOrientation orient = [[UIApplication sharedApplication] statusBarOrientation];
    
    Wt::WServer *server = Wt::WServer::instance();
    if( !server )
      return;
    
    const float top = insets.top;
    const float right = insets.right;
    const float bottom = insets.bottom;
    const float left = insets.left;
    const int orientation = orient;
    
    NSLog(@"setSafeAreasToClient: Orient=%i, {%f, %f, %f, %f}", orientation, top, right, bottom, left );
    
    server->postAll( std::bind([orientation,top,right,bottom,left](){
      InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
      if( !app ) //Shouldnt ever happen
        return;
      app->setSafeAreaInsets( orientation, top, right, bottom, left );
      app->triggerUpdate();
    }));
  }
}//setSafeAreasToClient


- (void)viewDidLayoutSubviews
{
  [super viewDidLayoutSubviews];
  
  NSLog(@"viewDidLayoutSubviews" );
}

@end
