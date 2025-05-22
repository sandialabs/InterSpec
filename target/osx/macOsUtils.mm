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
}//namespace macOsUtils
