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

#import <Cocoa/Cocoa.h>
#import <WebKit/WebKit.h>

@interface AppDelegate : NSObject <NSApplicationDelegate,WKNavigationDelegate,WKUIDelegate,NSURLDownloadDelegate>
@property (nonatomic,strong) WKWebView *InterSpecWebView;
//@property (nonatomic, strong) WKWebViewConfiguration *webConfig;

@property (nonatomic) BOOL isServing;
@property (nonatomic) int fileNeedsOpening;
@property (strong, nonatomic) NSString *UrlServingOn;
@property (strong, nonatomic) NSString *UrlUniqueId;
@property (strong, nonatomic) NSString *PreferenceDbPath;

@property (assign) IBOutlet NSWindow *window;

@property (readonly, strong, nonatomic) NSPersistentStoreCoordinator *persistentStoreCoordinator;
@property (readonly, strong, nonatomic) NSManagedObjectModel *managedObjectModel;
@property (readonly, strong, nonatomic) NSManagedObjectContext *managedObjectContext;

-(void)setDbDirectory;

-(void) terminated: (NSNotification *)notification;
-(BOOL) application:(NSApplication *)theApplication openFile:(NSString *)filename;
-(BOOL) applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)sender;
-(void)themeChanged:(NSNotification *) notification;
@end
