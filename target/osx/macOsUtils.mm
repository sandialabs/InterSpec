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

#include "InterSpec_config.h"

#import <dispatch/dispatch.h>
#import <Foundation/NSUserDefaults.h>

#import <AppKit/AppKit.h>
#import <AppKit/NSWorkspace.h>    //NSWorkspace
#import <Foundation/Foundation.h> //NSFileManager, NSBundle, NSURL

#include <Wt/WServer>
#include <Wt/WApplication>

#include "target/osx/macOsUtils.h"

namespace macOsUtils
{
  
void sessionSuccessfullyLoaded()
{
  dispatch_async( dispatch_get_main_queue(), ^{
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    [defaults setBool:NO forKey:@"DoNotResume"];
  });
}//void macOsUtils_sessionSuccessfullyLoaded()
  

std::string static_data_base_dir()
{
  return [[[NSBundle mainBundle] resourcePath] UTF8String];
}
  
  
std::string user_data_dir()
{
  NSFileManager *fileManager = [NSFileManager defaultManager];
  NSURL *appSupportURL = [[fileManager URLsForDirectory:NSApplicationSupportDirectory inDomains:NSUserDomainMask] lastObject];
  NSURL *appUrl = [appSupportURL URLByAppendingPathComponent:@"sandia.InterSpec"];
  return [[appUrl path] UTF8String];
}
  
bool openFinderToPath( const std::string &filepath )
{
  NSString *nsstr_filepath = [[NSString alloc] initWithCString:filepath.c_str() encoding:NSUTF8StringEncoding];
  
  return [[NSWorkspace sharedWorkspace] selectFile: nil inFileViewerRootedAtPath: nsstr_filepath];
}
  
bool showFileInFinder( const std::string &filepath )
{
  //  See https://chromium.googlesource.com/chromium/src/+/refs/heads/main/chrome/browser/platform_util_mac.mm for how this could/should be implemented, at least for BUILD_AS_OSX_APP
  NSFileManager *fileManager = [NSFileManager defaultManager];
  
  NSString *nsstr_filepath = [[NSString alloc] initWithCString:filepath.c_str() encoding:NSUTF8StringEncoding];
  BOOL valid = [fileManager fileExistsAtPath:nsstr_filepath];
  
  NSURL *nsurl_filepath = [NSURL fileURLWithPath:nsstr_filepath];
  [[NSWorkspace sharedWorkspace] activateFileViewerSelectingURLs:@[ nsurl_filepath ]];
  
  return valid;
}
  
  
void showFilePicker( const std::string title, const std::string message,
                           const bool canChooseFiles,
                           const bool canChooseDirectories,
                           const bool allowsMultipleSelection,
                           const std::function<void(const std::vector<std::string> &)> callback )
{
  Wt::WApplication * const app = wApp;
  assert( app );
  if( !app )
    throw std::runtime_error( "showDirectoryPicker must be called from WApplication main thread" );
  const std::string session_id = app->sessionId();
  
  // Perform UI operations on the main thread.
  dispatch_async(dispatch_get_main_queue(), ^{
    @autoreleasepool {
      NSOpenPanel *panel = [NSOpenPanel openPanel];

      if( !title.empty() )
        panel.title = [NSString stringWithUTF8String:title.c_str()];
      
      if( !message.empty() )
        panel.message = [NSString stringWithUTF8String:message.c_str()];

      panel.canChooseFiles = canChooseFiles;
      panel.canChooseDirectories = canChooseDirectories;
      panel.allowsMultipleSelection = allowsMultipleSelection;
      if( canChooseDirectories )
        panel.canCreateDirectories = YES;

      // Run the panel modally
      [panel beginWithCompletionHandler:^(NSModalResponse result) {
        std::vector<std::string> selectedPaths; //will be empty for cancel

        if( result == NSModalResponseOK )
        {
          for( NSURL *selectedURL in panel.URLs )
          {
            if( !selectedURL.isFileURL )
            {
              NSLog(@"Warning: non-file URL returned: %@", selectedURL);
              continue;
            }
            
            const char *fsRep = [selectedURL fileSystemRepresentation];
            if( fsRep )
            {
              selectedPaths.push_back( std::string(fsRep) );
            }else
            {
              NSLog(@"Warning: Could not get fileSystemRepresentation for URL: %@", selectedURL);
              // We could try [selectedURL path] here...
              // selectedPath = std::string([[selectedURL path] UTF8String]);
            }
          }//for (NSURL *selectedURL in panel.URLs)
         }//if( result == NSModalResponseOK )

        // Call the callback from the apps main Wt thread
        try
        {
          Wt::WServer * const server = Wt::WServer::instance();
          assert( server );
          if( !server )
            throw std::runtime_error( "showDirectoryPicker: could get WServer" );
            
          server->post(session_id, [selectedPaths,callback](){
            callback( selectedPaths );
            wApp->triggerUpdate();
          });
        }catch (const std::exception& e)
        {
          NSLog(@"Exception in C++ callback: %s", e.what());
        }//try / catch
      }];//[panel beginWithCompletion
    } // @autoreleasepool
  }); // dispatch_async
}//void showDirectoryPicker(...)


bool createAndStoreSecurityScopedBookmark(const std::string &path, const std::string &bookmarkKey)
{
  @autoreleasepool {
    NSString *nsPath = [NSString stringWithUTF8String:path.c_str()];
    NSURL *url = [NSURL fileURLWithPath:nsPath];
    
    if( !url )
    {
      NSLog(@"createAndStoreSecurityScopedBookmark: Could not create NSURL for path: %s", path.c_str());
      return false;
    }
    
    NSError *error = nil;
    NSData *bookmarkData = [url bookmarkDataWithOptions:NSURLBookmarkCreationWithSecurityScope
                                 includingResourceValuesForKeys:nil
                                                  relativeToURL:nil
                                                          error:&error];
    
    if( !bookmarkData || error )
    {
      NSLog(@"createAndStoreSecurityScopedBookmark: Failed to create bookmark for path %s: %s", 
            path.c_str(), error ? [[error localizedDescription] UTF8String] : "Unknown error");
      return false;
    }
    
    // Store in NSUserDefaults
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSString *nsKey = [NSString stringWithUTF8String:bookmarkKey.c_str()];
    [defaults setObject:bookmarkData forKey:nsKey];
    [defaults synchronize];
    
    NSLog(@"createAndStoreSecurityScopedBookmark: Successfully created and stored bookmark for path: %s", path.c_str());
    return true;
  }
}

bool resolveAndStartAccessingSecurityScopedBookmark(const std::string &bookmarkKey, std::string &resolvedPath)
{
  @autoreleasepool {
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSString *nsKey = [NSString stringWithUTF8String:bookmarkKey.c_str()];
    NSData *bookmarkData = [defaults objectForKey:nsKey];
    
    if( !bookmarkData )
    {
      NSLog(@"resolveAndStartAccessingSecurityScopedBookmark: No bookmark data found for key: %s", bookmarkKey.c_str());
      return false;
    }
    
    BOOL isStale = NO;
    NSError *error = nil;
    NSURL *url = [NSURL URLByResolvingBookmarkData:bookmarkData
                                           options:NSURLBookmarkResolutionWithSecurityScope
                                     relativeToURL:nil
                               bookmarkDataIsStale:&isStale
                                              error:&error];
    
    if( !url || error )
    {
      NSLog(@"resolveAndStartAccessingSecurityScopedBookmark: Failed to resolve bookmark for key %s: %s", 
            bookmarkKey.c_str(), error ? [[error localizedDescription] UTF8String] : "Unknown error");
      return false;
    }
    
    if( isStale )
    {
      NSLog(@"resolveAndStartAccessingSecurityScopedBookmark: Bookmark is stale for key: %s", bookmarkKey.c_str());
      // Could recreate the bookmark here if needed, but for now just warn
    }
    
    // Start accessing the security scoped resource
    BOOL accessStarted = [url startAccessingSecurityScopedResource];
    if( !accessStarted )
    {
      NSLog(@"resolveAndStartAccessingSecurityScopedBookmark: Failed to start accessing security scoped resource for key: %s", bookmarkKey.c_str());
      return false;
    }
    
    // Get the resolved path
    const char *fsRep = [url fileSystemRepresentation];
    if( fsRep )
    {
      resolvedPath = std::string(fsRep);
      NSLog(@"resolveAndStartAccessingSecurityScopedBookmark: Successfully resolved bookmark for key %s to path: %s", 
            bookmarkKey.c_str(), resolvedPath.c_str());
      return true;
    }
    else
    {
      NSLog(@"resolveAndStartAccessingSecurityScopedBookmark: Could not get filesystem representation for resolved URL");
      [url stopAccessingSecurityScopedResource]; // Clean up
      return false;
    }
  }
}

void stopAccessingSecurityScopedBookmark(const std::string &bookmarkKey)
{
  @autoreleasepool {
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSString *nsKey = [NSString stringWithUTF8String:bookmarkKey.c_str()];
    NSData *bookmarkData = [defaults objectForKey:nsKey];
    
    if( !bookmarkData )
    {
      NSLog(@"stopAccessingSecurityScopedBookmark: No bookmark data found for key: %s", bookmarkKey.c_str());
      return;
    }
    
    BOOL isStale = NO;
    NSError *error = nil;
    NSURL *url = [NSURL URLByResolvingBookmarkData:bookmarkData
                                           options:NSURLBookmarkResolutionWithSecurityScope
                                     relativeToURL:nil
                               bookmarkDataIsStale:&isStale
                                              error:&error];
    
    if( url && !error )
    {
      [url stopAccessingSecurityScopedResource];
      NSLog(@"stopAccessingSecurityScopedBookmark: Stopped accessing security scoped resource for key: %s", bookmarkKey.c_str());
    }
    else
    {
      NSLog(@"stopAccessingSecurityScopedBookmark: Could not resolve bookmark to stop accessing for key: %s", bookmarkKey.c_str());
    }
  }
}

void removeSecurityScopedBookmark(const std::string &bookmarkKey)
{
  @autoreleasepool {
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSString *nsKey = [NSString stringWithUTF8String:bookmarkKey.c_str()];
    
    // Stop accessing the resource first
    stopAccessingSecurityScopedBookmark(bookmarkKey);
    
    // Remove from defaults
    [defaults removeObjectForKey:nsKey];
    [defaults synchronize];
    
    NSLog(@"removeSecurityScopedBookmark: Removed bookmark for key: %s", bookmarkKey.c_str());
  }
}

std::vector<std::string> getSecurityScopedBookmarkKeys(const std::string &keyPrefix)
{
  std::vector<std::string> keys;
  
  @autoreleasepool {
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSDictionary *allDefaults = [defaults dictionaryRepresentation];
    NSString *nsPrefix = [NSString stringWithUTF8String:keyPrefix.c_str()];
    
    for( NSString *key in allDefaults )
    {
      if( [key hasPrefix:nsPrefix] )
      {
        keys.push_back([key UTF8String]);
      }
    }
  }
  
  return keys;
}

}//namespace macOsUtils
