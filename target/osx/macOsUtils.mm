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

#import <AppKit/NSWorkspace.h>    //NSWorkspace
#import <Foundation/Foundation.h> //NSFileManager, NSBundle, NSURL

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
}//namespace macOsUtils
