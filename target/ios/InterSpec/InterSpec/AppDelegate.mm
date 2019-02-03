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

#include <string>

#include <boost/filesystem.hpp>

#include "InterSpec/InterSpecServer.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"



@implementation AppDelegate


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
    static const char *argv6 = "--tempdir";
    static const char *argv7 = tmpdirstr.c_str();
    static const char *argv[] = { argv0, argv1, argv4, argv5, argv6, argv7 };
    int argc = sizeof(argv) / sizeof(argv[0]);
    InterSpecServer::startServer( argc, (char **)argv, &createThisApplication );
  }
  
  
  NSString *userAgentStr;
  
  if( [[UIDevice currentDevice] userInterfaceIdiom] == UIUserInterfaceIdiomPhone )
  {
    userAgentStr = @"Mozilla/5.0 (iPhone; CPU iPhone OS 6_1 like Mac OS X) AppleWebKit/536.26 (KHTML, like Gecko) Version/6.0 Mobile/10B141 Safari/8536.25";
  }else
  {
    userAgentStr = @"Mozilla/5.0 (iPad; CPU OS 6_1 like Mac OS X) AppleWebKit/536.26 (KHTML, like Gecko) Version/6.0 Mobile/10B141 Safari/8536.25";
  }
  
  NSLog( @"Setting UserAgent to %@", userAgentStr );
  
  // Set user agent (the only problem is that we can't modify the User-Agent later in the program)
  NSDictionary *dictionnary = [[NSDictionary alloc] initWithObjectsAndKeys:userAgentStr, @"UserAgent", nil];
  [[NSUserDefaults standardUserDefaults] registerDefaults:dictionnary];
}


-(BOOL)application:(UIApplication *)application
           openURL:(NSURL *)url
 sourceApplication:(NSString *)sourceApplication
        annotation:(id)annotation
{
  // Make sure url indicates a file (as opposed to, e.g., http://)
  if (url != nil && [url isFileURL])
  {
    NSLog(@"File to open %@ from %@", url, sourceApplication );
    
    [_viewController openSpectrumFile: url];
  }
  
  // Indicate that we have successfully opened the URL
  return YES;
}

- (BOOL)application:(UIApplication *)application didFinishLaunchingWithOptions:(NSDictionary *)launchOptions
{
  // Override point for customization after application launch.
  NSLog(@"didFinishLaunchingWithOptions");
  [[UIApplication sharedApplication] setMinimumBackgroundFetchInterval:UIApplicationBackgroundFetchIntervalMinimum];
  self.window = [[UIWindow alloc] initWithFrame:[[UIScreen mainScreen] bounds]];
  
  if ([[UIDevice currentDevice] userInterfaceIdiom] == UIUserInterfaceIdiomPhone) {
    self.viewController = [[ViewController alloc] initWithNibName:@"ViewController_iPhone" bundle:nil];
  } else {
    self.viewController = [[ViewController alloc] initWithNibName:@"ViewController_iPad" bundle:nil];
  }
  self.window.rootViewController = self.viewController;
  [self.window makeKeyAndVisible];
  self.isInBackground = NO;
  self.backgroundTask = UIBackgroundTaskInvalid;
  
  
  //Nesary for iOS 9, when keyboard hides, need to re-align the app
  NSNotificationCenter *center = [NSNotificationCenter defaultCenter];
  //  [center addObserver:self selector:@selector(didShowKeyboard) name:UIKeyboardDidShowNotification object:nil];
  [center addObserver:_viewController selector:@selector(onKeyboardHide) name:UIKeyboardDidHideNotification object:nil]; //UIKeyboardWillHideNotification
  
  
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
    dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
      
      //### background task starts
      NSLog(@"Running in the background\n");
      double timeRemain = [[UIApplication sharedApplication] backgroundTimeRemaining];
      NSLog(@"Initial background time: %f", timeRemain);
      while(timeRemain > 25)  //well give ourselves 25 seconds to save the current InterSpec apps states
      {
        if( !self.isInBackground )
        {
          NSLog(@"isInBackground set");
          break;
        }
        timeRemain = [[UIApplication sharedApplication] backgroundTimeRemaining];
        NSLog(@"Background time Remaining: %f", timeRemain);
        [NSThread sleepForTimeInterval:5]; //wait for 5 sec
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
  NSLog(@"Begining sendSpectrumFileToOtherApp" );

  
  _documentController = [UIDocumentInteractionController interactionControllerWithURL:[NSURL fileURLWithPath:filename]];
  _documentController.delegate = _viewController;
//  [_documentController retain];
  _documentController.UTI = @"sandia.InterSpec.spectrum";
  
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




@end
