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

@interface ViewController : UIViewController <UIWebViewDelegate,UIDocumentInteractionControllerDelegate>
@property (weak, nonatomic) IBOutlet UIWebView *webView;
@property (weak, nonatomic) IBOutlet UIWebView *padWebView;
@property (weak, nonatomic) IBOutlet UIWebView *phoneWebView;

@property (nonatomic) BOOL isServing;
@property (nonatomic) BOOL appHasGoneIntoBackground;
@property (nonatomic) BOOL appComminFromBackground;
@property (nonatomic) int dbIndexOfFileToOpen;
@property (nonatomic) NSString *UrlServingOn;
@property (nonatomic) NSString *UrlUniqueId;

- (void)webView:(UIWebView *)webView didFailLoadWithError:(NSError *)error;
- (void)webViewDidFinishLoad:(UIWebView *)webView;
- (void)webViewDidStartLoad:(UIWebView *)webView;
//- (BOOL)webView:(UIWebView *)webView shouldStartLoadWithRequest:(NSURLRequest *)request navigationType:(UIWebViewNavigationType)navigationType;
-(BOOL) webView:(UIWebView *)inWeb shouldStartLoadWithRequest:(NSURLRequest *)inRequest navigationType:(UIWebViewNavigationType)inType;
-(BOOL)startServer;
-(BOOL)openSpectrumFile:(NSURL *)url;
-(BOOL)enteredBackground;
-(BOOL)willEnterForeground;
-(BOOL)wakeupFromBackground;
-(void)onKeyboardHide;
-(void)fixPermissions:(NSString *)path;

/*
//implemented UIDocumentInteractionControllerDelegate methods
//  Commented out as they dont seem to ever be called...
- (void)documentInteractionController:(UIDocumentInteractionController *)controller willBeginSendingToApplication:(NSString *)application;
- (void)documentInteractionController:(UIDocumentInteractionController *)controller didEndSendingToApplication:(NSString *)application;
- (void)documentInteractionControllerDidDismissOpenInMenu:(UIDocumentInteractionController *)controller;
*/
@end
