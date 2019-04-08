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

#import <UIKit/UIKit.h>
#import <WebKit/WebKit.h>

@interface ViewController : UIViewController <WKUIDelegate,WKNavigationDelegate,UIDocumentInteractionControllerDelegate>
@property (nonatomic) WKWebView *webView;

@property (nonatomic) BOOL isServing;
@property (nonatomic) BOOL appHasGoneIntoBackground;
@property (nonatomic) BOOL appComminFromBackground;
@property (nonatomic) int dbIndexOfFileToOpen;
@property (nonatomic) NSString *UrlServingOn;
@property (nonatomic) NSString *UrlUniqueId;

/*
- (void)webView:(UIWebView *)webView didFailLoadWithError:(NSError *)error;
- (void)webViewDidFinishLoad:(UIWebView *)webView;
- (void)webViewDidStartLoad:(UIWebView *)webView;
//- (BOOL)webView:(UIWebView *)webView shouldStartLoadWithRequest:(NSURLRequest *)request navigationType:(UIWebViewNavigationType)navigationType;
-(BOOL) webView:(UIWebView *)inWeb shouldStartLoadWithRequest:(NSURLRequest *)inRequest navigationType:(UIWebViewNavigationType)inType;
*/

//WKNavigationDelegate items
- (void)webView:(WKWebView *)webView didFailNavigation:(WKNavigation *)navigation withError:(NSError *)error;
- (void)webView:(WKWebView *)webView didCommitNavigation:(WKNavigation *)navigation;
- (void)webView:(WKWebView *)webView didFinishNavigation:(WKNavigation *)navigation;
- (void)webView:(WKWebView *)webView decidePolicyForNavigationAction:(WKNavigationAction *)navigationAction
                 decisionHandler:(void (^)(WKNavigationActionPolicy))decisionHandler;

// WKUIDelegate items
//- (void)webView:(WKWebView *)webView
//                runOpenPanelWithParameters:(WKOpenPanelParameters *)parameters
//                initiatedByFrame:(WKFrameInfo *)frame
//                completionHandler:(void (^)(NSArray<NSURL *> *URLs))completionHandler;
- (BOOL)webView:(WKWebView *)webView
        shouldPreviewElement:(WKPreviewElementInfo *)elementInfo;
- (void)webViewDidClose:(WKWebView *)webView;

//- (WKWebView *)webView:(WKWebView *)webView
//createWebViewWithConfiguration:(WKWebViewConfiguration *)configuration
//   forNavigationAction:(WKNavigationAction *)navigationAction
//        windowFeatures:(WKWindowFeatures *)windowFeatures;

//- (void)webView:(WKWebView *)webView
//runJavaScriptAlertPanelWithMessage:(NSString *)message
//initiatedByFrame:(WKFrameInfo *)frame
//completionHandler:(void (^)(void))completionHandler;



-(BOOL)startServer;
-(BOOL)openSpectrumFile:(NSURL *)url;
-(BOOL)enteredBackground;
-(BOOL)willEnterForeground;
-(BOOL)wakeupFromBackground;
-(void)onKeyboardHide;
-(void)fixPermissions:(NSString *)path;


- (void)viewDidLayoutSubviews;
/*
//implemented UIDocumentInteractionControllerDelegate methods
//  Commented out as they dont seem to ever be called...
- (void)documentInteractionController:(UIDocumentInteractionController *)controller willBeginSendingToApplication:(NSString *)application;
- (void)documentInteractionController:(UIDocumentInteractionController *)controller didEndSendingToApplication:(NSString *)application;
- (void)documentInteractionControllerDidDismissOpenInMenu:(UIDocumentInteractionController *)controller;
*/
@end
