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

#import "ViewController.h"

#include <stdlib.h>


#include <string>
#include <fstream>
#include <iostream>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

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
  
  _appHasGoneIntoBackground = YES;
  _appComminFromBackground = YES;
  
  InterSpecServer::killServer();
  NSLog( @"Killed server" );
  
  _isServing = NO;
  _UrlServingOn = @"";
  _UrlUniqueId = @"";

  //We could actually delete the webview here and recreate when bringing to
  //  the foreground, but whatever.  For now lets just clear the webview to
  //  save a couple megabytes of ram.
  [_webView loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:@"about:blank"]]];

//The following doesnt seem to make a big difference, so lets not bother with it
//see: http://twobitlabs.com/2012/01/ios-ipad-iphone-nsurlcache-uiwebview-memory-utilization/
  [[NSURLCache sharedURLCache] removeAllCachedResponses];
  
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
      InterSpecApp *app = InterSpecApp::instanceFromExtenalIdString( [_UrlUniqueId UTF8String] );

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
  
  //
  _UrlUniqueId = [NSString stringWithFormat:@"%i", arc4random_uniform(10000000)];
  NSString *actualURL = [NSString stringWithFormat:@"%@?externalid=%@", _UrlServingOn, _UrlUniqueId];
  
  //Check to see if we should open a spectrum file; if so append which on to URL
  if( _dbIndexOfFileToOpen >= 0 )
    actualURL = [NSString stringWithFormat:@"%@&specfile=%i", actualURL, _dbIndexOfFileToOpen];
  _dbIndexOfFileToOpen = -1;
  
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
  [_webView setScalesPageToFit: NO];
  [_webView loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:actualURL]]];
  
  return YES;
}//-(BOOL)wakeupFromBackground




-(BOOL)openSpectrumFile:(NSURL *)url
{
  const std::string urlstr = [[url path] UTF8String];
  
  if( boost::filesystem::exists( urlstr ) )
  {
    NSLog( @"Will try to open %s", urlstr.c_str() );
    
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
    
//We cant use this next block of code because wakeupFromBackground: hasnt been
//  called yet, so a valid InterSpecApp wont exist, since we dont instruct
//  the webview to load the URL until wakeupFromBackground:
//  if( _UrlUniqueId.length > 0 && _dbIndexOfFileToOpen >= 0 )
//  {
//    InterSpecApp *app = InterSpecApp::instanceFromExtenalIdString( [_UrlUniqueId UTF8String] );
//    if( app )
//    {
//      Wt::WApplication::UpdateLock applock( app );
//      Wt::WApplication *wtapp = Wt::WApplication::instance();
//      if( wtapp && applock )
//      {
//        app->openFileFromDbFileSystemLink( _dbIndexOfFileToOpen ) )
//        app->triggerUpdate();
//      }
//    }
//  }//if( _dbIndexOfFileToOpen >= 0 )
  }else
  {
    NSLog( @"The file %@ could not be accessed", [url path] );
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
  
  self.edgesForExtendedLayout = UIRectEdgeNone;
  [super viewDidLoad];
  
  if( IS_IPAD() )
    _webView = _padWebView;
  else
    _webView = _phoneWebView;
  
  NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
  NSString *documentsDirectory = [paths objectAtIndex:0];
  
  //The database needs the sqlite ending for some reason.  If I use the name
  //  SRBPreferences.db, then it doesnt work
  NSString *writableDBPath = [documentsDirectory stringByAppendingPathComponent:@"data.sqlite"];
  DataBaseUtils::setPreferenceDatabaseFile( [writableDBPath UTF8String] );
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
  
  [_webView setScalesPageToFit: NO];
  UIScrollView *view = [_webView scrollView];
  [view setScrollEnabled: NO];
  self.webView.delegate = self;
  
  NSLog( @"Done in viewDidLoad" );
}

- (void)webView:(UIWebView *)webView didFailLoadWithError:(NSError *)error
{
  NSLog(@"\n\n\ndidFailLoadWithError: error => %@ ", [error userInfo] );
}

- (void)onKeyboardHide
{
  //This function is vestigiual from debugging an issue with the keyboard hiding
//  NSLog(@"Triggering screen resize" );
//  NSString *jsstring = [NSString stringWithFormat:@"setTimeout( function(){$('.Wt-domRoot').height(window.innerHeight); window.scrollTo(0,0);}, 0 );" ];
//  [_webView stringByEvaluatingJavaScriptFromString: jsstring];
}


- (void)webViewDidFinishLoad:(UIWebView *)webView
{
  NSLog(@"WebViewDidFinishLoad" );
  
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
  return self;
}

- (UIView *)documentInteractionControllerViewForPreview:(UIDocumentInteractionController *)controller
{
  //required to implement UIDocumentInteractionControllerDelegate
  return self.view;
}

- (CGRect)documentInteractionControllerRectForPreview:(UIDocumentInteractionController *)controller
{
  //required to implement UIDocumentInteractionControllerDelegate
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

@end
