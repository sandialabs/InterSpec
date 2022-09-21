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

#import "AppDelegate.h"

#import "ViewController.h"

#define USE_SHOW_FINGERTIP 0

#if( USE_SHOW_FINGERTIP )
#import "MBFingerTipWindow.h"
#endif

#include <string>

#include <boost/filesystem.hpp>

#include "InterSpec/InterSpecServer.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"


@interface InterSpecAppView : UIWindow
{
}
@end

@implementation InterSpecAppView
- (id) initWithFrame: (CGRect)frame
{
  self = [super initWithFrame:frame];
  if( self )
  {
    NSLog( @"Created InterSpecAppView" );
  }
  return self;
}

- (void)safeAreaInsetsDidChange
{
  NSLog( @"safeAreaInsetsDidChange" );
  
  if( @available(iOS 11, *) )
  {
    [super safeAreaInsetsDidChange];
    //blah blah blah
    if( [self insetsLayoutMarginsFromSafeArea] )
    {
      NSLog( @"safeAreaInsetsDidChange" );
    }
  }
}//safeAreaInsetsDidChange

@end  //InterSpecAppView



@implementation AppDelegate
 
#if( USE_SHOW_FINGERTIP )
- (MBFingerTipWindow *)window
{
  static MBFingerTipWindow *visualFeedbackWindow = nil;
  if (!visualFeedbackWindow) visualFeedbackWindow = [[MBFingerTipWindow alloc] initWithFrame:[[UIScreen mainScreen] bounds]];
  [visualFeedbackWindow setAlwaysShowTouches: true];
  return visualFeedbackWindow;
}
#endif



Wt::WApplication *createThisApplication(const Wt::WEnvironment& env)
{
  return new InterSpecApp( env );
}

+ (void)initialize
{
  static const std::string basedir = std::string("--basedir=")
                             + [[[NSBundle mainBundle] resourcePath] UTF8String];
  
  boost::filesystem::path apppath = [[[NSBundle mainBundle] executablePath] UTF8String];
  apppath = apppath.parent_path();
  NSLog( @"Setting CWD=%s", apppath.string<std::string>().c_str() );
  boost::filesystem::current_path( apppath );
  
  //Since app is run in a container, we need to set the correct temp directory
  //  for Wt, especially for "uploading" files.
  NSString *tempDir = NSTemporaryDirectory();
  if( tempDir == nil ) //shouldnt ever fail, right
    tempDir = @"/tmp";
  static const std::string tmpdirstr = [tempDir UTF8String];
  
  {
    static const char *argv0 = "--docroot";
    static const char *argv1 = ".";
    static const char *argv4 = "-c";
    static const char *argv5 = "./data/config/wt_config_ios.xml";
    static const char *argv6 = "--accesslog=-";
    static const char *argv7 = "--tempdir";
    static const char *argv8 = tmpdirstr.c_str();
    static const char *argv[] = { argv0, argv1, argv4, argv5, argv6, argv7, argv8 };
    int argc = sizeof(argv) / sizeof(argv[0]);
    InterSpecServer::startServer( argc, (char **)argv, &createThisApplication );
  }
}

-(BOOL)application:(UIApplication *)application
           openURL:(NSURL *)url
           options:(NSDictionary<UIApplicationOpenURLOptionsKey, id> *)options
{
  //Will see if we should open
  if( ![url isFileURL] ){
    NSLog(@"Try to handle non-file URL with scheme '%@'.", [url scheme] );
    
    if ([[url scheme] caseInsensitiveCompare:@"interspec"] == NSOrderedSame) {
      return [_viewController handleURL: url];
    }
    
    NSLog(@"URL passed in is not a file and not 'interspec' scheme" );
    return NO;
  }//if( not a file URL )
  
  NSLog(@"Will see if we should open file in place for URL %@", url );
  BOOL openInPlace = false;
  
  if( options && ([options objectForKey: UIApplicationOpenURLOptionsOpenInPlaceKey] != nil) )
    openInPlace = [options[UIApplicationOpenURLOptionsOpenInPlaceKey] boolValue];
  
  if( ![url isFileURL] ){
    NSLog(@"URL passed in is not a file" );
    return NO;
  }
  
  BOOL didStartAccessing = [url startAccessingSecurityScopedResource];
  if( !didStartAccessing ) {
    NSLog(@"Error getting access to %@", url );
  }
  
  BOOL success = NO;
  if( !openInPlace )
  {
    NSLog(@"openInPlace=False" );
    // TODO: we should copy the file over...
    success = [_viewController openSpectrumFile: url];
  }else
  {
    NSLog(@"openInPlace=True" );
    success = [_viewController openSpectrumFile: url];
  }
  
  if( didStartAccessing )
    [url stopAccessingSecurityScopedResource];
  
  return success;
}


-(BOOL)application:(UIApplication *)application
           openURL:(NSURL *)url
 sourceApplication:(NSString *)sourceApplication
        annotation:(id)annotation
{
  BOOL success = NO;
  
  // Make sure url indicates a file (as opposed to, e.g., http://)
  if (url != nil && [url isFileURL])
  {
    BOOL didStartAccessing = [url startAccessingSecurityScopedResource];
    if( !didStartAccessing ) {
      NSLog(@"Error getting access to %@", url );
    }
    
    NSLog(@"File to open %@ from %@", url, sourceApplication );
    
    success = [_viewController openSpectrumFile: url];
    
    if( didStartAccessing )
      [url stopAccessingSecurityScopedResource];
  }
  
  // Indicate if we have successfully opened the URL
  return success;
}


- (BOOL)application:(UIApplication *)application didFinishLaunchingWithOptions:(NSDictionary *)launchOptions
{
  // Override point for customization after application launch.
  NSLog(@"didFinishLaunchingWithOptions");
  [[UIApplication sharedApplication] setMinimumBackgroundFetchInterval:UIApplicationBackgroundFetchIntervalMinimum];
  
  self.window = [[InterSpecAppView alloc] initWithFrame:[[UIScreen mainScreen] bounds]];
  self.viewController = [[ViewController alloc] initWithNibName: nil bundle: nil];
  self.window.rootViewController = self.viewController;
  [self.window makeKeyAndVisible];
  self.isInBackground = NO;
  self.backgroundTask = UIBackgroundTaskInvalid;
  
  //Nesary for iOS 9, when keyboard hides, need to re-align the app
  NSNotificationCenter *center = [NSNotificationCenter defaultCenter];
  //  [center addObserver:self selector:@selector(didShowKeyboard) name:UIKeyboardDidShowNotification object:nil];
  [center addObserver:_viewController selector:@selector(onKeyboardHide) name:UIKeyboardDidHideNotification object:nil]; //UIKeyboardWillHideNotification
  
  
  //Register to get notified of device orientation changes
  UIDevice *device = [UIDevice currentDevice];
  [device beginGeneratingDeviceOrientationNotifications];
  NSNotificationCenter *nc = [NSNotificationCenter defaultCenter];
  [nc addObserver:self
         selector:@selector(orientationChanged:)
             name:UIDeviceOrientationDidChangeNotification
           object:device];
  
  
  return YES;
}



- (BOOL)application:(UIApplication *)application shouldSaveApplicationState:(NSCoder *)coder
{
  //Called when going into the background, even if will continue to run for a
  //  limitied time in the background
  NSLog(@"shouldSaveApplicationState");
  return YES;
}

- (BOOL)application:(UIApplication *)application shouldRestoreApplicationState:(NSCoder *)coder
{
  //Called on start up if there is a state to restore
  NSLog(@"shouldRestoreApplicationState");

  //We could probably use a variable to communicate to the ViewController that
  //  we can try to restore the app state, but for now just use the UserDefaults
  [[NSUserDefaults standardUserDefaults] setBool:YES forKey:@"ShouldRestore"];
  [[NSUserDefaults standardUserDefaults] synchronize];
  
  return NO;
}

- (void)applicationWillResignActive:(UIApplication *)application
{
  //Sent when the application is about to move from active to inactive state.
  //  This can occur for certain types of temporary interruptions (such as an
  //  incoming phone call or SMS message) or when the user quits the application
  //  and it begins the transition to the background state.
  //Use this method to pause ongoing tasks, disable timers, and throttle down
  //  OpenGL ES frame rates. Games should use this method to pause the game.
  NSLog(@"applicationWillResignActive");
}

- (void)applicationDidEnterBackground:(UIApplication *)application
{
  //Use this method to release shared resources, save user data, invalidate
  //  timers, and store enough application state information to restore your
  //  application to its current state in case it is terminated later.
  //If your application supports background execution, this method is called
  //  instead of applicationWillTerminate: when the user quits.
  NSLog(@"applicationDidEnterBackground");
  
  self.isInBackground = YES;
  
  //If we can run a background thread (typically for at most 180 seconds), lets
  //  not kill the Wt server immediately, and instead waite untill we're almost
  //  out of time to do this
  //Code modified from http://hayageek.com/ios-background-task/
  if([[UIDevice currentDevice] respondsToSelector:@selector(isMultitaskingSupported)])
  {
    NSLog(@"Multitasking Supported");
    
    self.backgroundTask = [application beginBackgroundTaskWithExpirationHandler:^ {
      
      //Clean up code. Tell the system that we are done.
      [application endBackgroundTask: self.backgroundTask];
      self.backgroundTask = UIBackgroundTaskInvalid;
    }];
    
    //To make the code block asynchronous
    //dispatch_async(dispatch_get_main_queue(), ^{
    
    dispatch_async( dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
      //### background task starts
      NSLog(@"Running in the background\n");
      
      double (^getTimeRemainingBlock)(void) = ^{
        __block double timeRemain = 0.0;
      
        void (^getTimeRemainingWorker)(void) = ^{
          timeRemain = [[UIApplication sharedApplication] backgroundTimeRemaining];
          //NSLog(@"From main thread time remaining: %f", timeRemain);
        };
      
        if( [NSThread isMainThread] )
        {
          getTimeRemainingWorker();
        }else
        {
          dispatch_group_t group = dispatch_group_create();
          dispatch_group_async(group, dispatch_get_main_queue(), getTimeRemainingWorker );
          dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
          //NSLog(@"Should be done getting time remaining: %f", timeRemain);
        }
        
        //NSLog(@"Returning time remaining: %f", timeRemain);
        return timeRemain;
      };//getTimeRemainingBlock
      
      double timeleft = getTimeRemainingBlock();
      
      NSLog(@"Initial background time: %f", timeleft);
      while(timeleft > 15)  //well give ourselves 15 seconds to save the current InterSpec apps states
      {
        if( !self.isInBackground )
        {
          NSLog(@"isInBackground set");
          break;
        }
        timeleft = getTimeRemainingBlock();
        //NSLog(@"Background time Remaining: %f", timeleft);
        [NSThread sleepForTimeInterval: 0.25]; //wait for 0.25 sec
      }
       
      //#### background task ends
      if( self.isInBackground )
      {
        NSLog(@"Will Notify controller to kill server" );
//        [[NSUserDefaults standardUserDefaults] setBool:YES forKey:@"ShouldRestore"];
//        [[NSUserDefaults standardUserDefaults] synchronize];
        [self.viewController enteredBackground];
      }else
      {
        NSLog(@"Didnt have to kill server" );
      }
      
      //Clean up code. Tell the system that we are done.
      [application endBackgroundTask: self.backgroundTask];
      self.backgroundTask = UIBackgroundTaskInvalid;
    });
  }
  else
  {
    NSLog(@"Multitasking Not Supported");
    [self.viewController enteredBackground];
  }
}

- (void)applicationWillEnterForeground:(UIApplication *)application
{
  //Called as part of the transition from the background to the inactive state;
  //  here you can undo many of the changes made on entering the background.
  NSLog(@"applicationWillEnterForeground");
  [self.viewController willEnterForeground];
}

- (void)applicationDidBecomeActive:(UIApplication *)application
{
  //Restart any tasks that were paused (or not yet started) while the
  //  application was inactive. If the application was previously in the
  //  background, optionally refresh the user interface.
  NSLog(@"applicationDidBecomeActive");
  
  if( self.backgroundTask != UIBackgroundTaskInvalid )
  {
    NSLog(@"Will end background task");
    self.isInBackground = NO;
  }
  
  [self.viewController wakeupFromBackground];
}

- (void)applicationWillTerminate:(UIApplication *)application
{
  //Called when the application is about to terminate. Save data if appropriate.
  //  See also applicationDidEnterBackground:.
  NSLog(@"applicationWillTerminate");
}


-(void)sendSpectrumFileToOtherApp: (NSString *) filename;
{
  if( !filename )
  {
    NSLog(@"sendSpectrumFileToOtherApp: No file-name passed in." );
    return;
  }
  
  NSLog(@"sendSpectrumFileToOtherApp: Starting for '%@'.", filename );
  
  NSString *ext = [filename pathExtension];
  NSLog(@"sendSpectrumFileToOtherApp: fileextension '%@'.", ext );
  
  _documentController = [UIDocumentInteractionController interactionControllerWithURL:[NSURL fileURLWithPath:filename]];
  _documentController.delegate = _viewController;
//  [_documentController retain];
  if( _documentController.UTI == nil )
  {
    NSLog(@"The UTI type is unknown, we'll try to fix this up." );
    
    // To get the UTI, on macOS use command `mdls -name kMDItemContentType myfile.txt`
    if ([ext caseInsensitiveCompare:@"CSV"] == NSOrderedSame)
      _documentController.UTI = @"public.comma-separated-values-text"; //kUTTypeCommaSeparatedText;
    else if ([ext caseInsensitiveCompare:@"XML"] == NSOrderedSame)
      _documentController.UTI = @"public.xml";
    else if ([ext caseInsensitiveCompare:@"PNG"] == NSOrderedSame)
      _documentController.UTI = @"public.png";
    else if ([ext caseInsensitiveCompare:@"SVG"] == NSOrderedSame)
      _documentController.UTI = @"public.svg-image";
    else if ([ext caseInsensitiveCompare:@"HTML"] == NSOrderedSame)
      _documentController.UTI = @"public.html";
    else if ([ext caseInsensitiveCompare:@"CALp"] == NSOrderedSame)
      _documentController.UTI = @"sandia.InterSpec.spectrum";
    else // Assume a spectrum file
      _documentController.UTI = @"sandia.InterSpec.spectrum"; // on macOS this is "gov.sandia.interspec.gamma-spectrum"
  }else
  {
    NSLog(@"The UTI type is '%@'.", _documentController.UTI );
  }
  
  //CGRect navRect = _viewController.view.frame;  //CGRectZero
  //On the iPad the dialog doesnt show up, so lets just hack it for the moment.
  CGRect navRect = CGRectMake(0.0, 0.0, 0.0, 0.0);

  [_documentController presentOptionsMenuFromRect:navRect inView:_viewController.view animated:YES];
 
  //I get the message: Attempting to load the view of a view controller while it
  //                   is deallocating is not allowed and may result in
  //                   undefined behavior
  //at this point...
  //see maybye: http://stackoverflow.com/questions/32282401/attempting-to-load-the-view-of-a-view-controller-while-it-is-deallocating-uis
}


//********** ORIENTATION CHANGED **********
- (void)orientationChanged:(NSNotification *)note
{
  //NSLog(@"Orientation  has changed: %ld", [[note object] orientation]);
  
  if( _viewController )
    [_viewController setSafeAreasToClient];
}//orientationChanged


- (void)openURL:(NSURL *)url
          options:(NSDictionary<UIApplicationOpenExternalURLOptionsKey, id> *)options
          completionHandler:(void (^)(BOOL success))completion
//-(BOOL)application:(UIApplication *)app openURL:(NSURL *)url options:(NSDictionary<UIApplicationOpenURLOptionsKey, id> *)options;
{
  NSLog(@"url received: %@", url);
  
  const BOOL status = [_viewController handleURL: url];
  completion( status );
}

@end
