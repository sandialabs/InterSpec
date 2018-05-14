//
//  AppDelegate.h
//  InterSpec OSX
//
//  Created by Johnson, William C on 7/3/13.
//  Copyright (c) 2013 Johnson, William C. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <WebKit/WebKit.h>
#import <WebKit/WebView.h>
#import <WebKit/WebUIDelegate.h>
#import <WebKit/WebPolicyDelegate.h>

//Need to create the ability to pick location and name of download files, so we
//  need to implement a download delegate.
@interface MyDownloadDelegate : NSURLDownload{
}
- (void)download:(NSURLDownload *)download decideDestinationWithSuggestedFilename:(NSString *)filename;
- (void)download:(NSURLDownload *)download didFailWithError:(NSError *)error;
- (void)downloadDidFinish:(NSURLDownload *)download;
- (void)startDownloadingURL:sender;
- (void)downloadDidBegin:(NSURLDownload *)download;
- (void)download:(NSURLDownload *)download didCreateDestination:(NSString *)path;
@end


@interface InterSpecWebPolicyDecisionListener: NSObject {
}
- (void)webView:(WebView *)sender runOpenPanelForFileButtonWithResultListener:(id < WebOpenPanelResultListener >)resultListener;
- (void)webView:(WebView *)webView decidePolicyForNavigationAction:(NSDictionary *)actionInformation
    request:(NSURLRequest *)request
    frame:(WebFrame *)frame
    decisionListener:(id<WebPolicyDecisionListener>)listener;

- (void)webView:(WebView *)webView decidePolicyForNewWindowAction:(NSDictionary *)actionInformation
    request:(NSURLRequest *)request
    newFrameName:(NSString *)frameName
    decisionListener:(id <WebPolicyDecisionListener>)listener;

  //disable the right-click context menu
- (NSArray *)webView:(WebView *)sender contextMenuItemsForElement:(NSDictionary *)element defaultMenuItems:(NSArray *)defaultMenuItems;
@end

@interface AppDelegate : NSObject <NSApplicationDelegate>
@property IBOutlet WebView *InterSpecWebView;
@property (strong) MyDownloadDelegate *UrlDownload;

@property (nonatomic) BOOL isServing;
@property (nonatomic) int fileNeedsOpening;
@property (strong, nonatomic) NSString *UrlServingOn;
@property (strong, nonatomic) NSString *UrlUniqueId;
@property (strong, nonatomic) NSString *PreferenceDbPath;

@property (assign) IBOutlet NSWindow *window;

@property (readonly, strong, nonatomic) NSPersistentStoreCoordinator *persistentStoreCoordinator;
@property (readonly, strong, nonatomic) NSManagedObjectModel *managedObjectModel;
@property (readonly, strong, nonatomic) NSManagedObjectContext *managedObjectContext;
@property (unsafe_unretained) IBOutlet NSTextView *textView;

- (void)setDbDirectory;

- (IBAction)saveAction:(id)sender;
- (IBAction)ShowInBrowserClicked:(id)sender;

-(void) terminated: (NSNotification *)notification;
-(BOOL)application:(NSApplication *)theApplication openFile:(NSString *)filename;
-(BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)sender;


//Implement methods necassary for WebPolicyDelegate that controls how to deal
//  with opening new windows, or downloading files.  For now we'll just
//  modify how we open new windows since these are CSV or file downloads.
- (void)webView:(WebView *)webView
  decidePolicyForNewWindowAction:(NSDictionary *)actionInformation
                         request:(NSURLRequest *)request
                    newFrameName:(NSString *)frameName
                decisionListener:(id < WebPolicyDecisionListener >)listener;
/*
- (void)webView:(WebView *)webView
  decidePolicyForNavigationAction:(NSDictionary *)actionInformation
                          request:(NSURLRequest *)request
                            frame:(WebFrame *)frame
                 decisionListener:(id < WebPolicyDecisionListener >)listener;

- (void)webView:(WebView *)webView
  decidePolicyForMIMEType:(NSString *)type
                  request:(NSURLRequest *)request
                    frame:(WebFrame *)frame
         decisionListener:(id < WebPolicyDecisionListener >)listener;
 */
@end
