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

#import <Cocoa/Cocoa.h>
#import <QuickLookUI/QuickLookUI.h>

#include "SpecPreviewCommon.h"


@interface PreviewViewController : NSViewController <QLPreviewingController>
@end

@implementation PreviewViewController

- (NSString *)nibName
{
  return nil;
}

- (void)loadView
{
  self.view = [[NSView alloc] initWithFrame:NSMakeRect( 0, 0, 800, 600 )];
}

- (void)preparePreviewOfFileAtURL:(NSURL *)url
                completionHandler:(void (^)(NSError * _Nullable))handler
{
  const char *filePath = [[url path] fileSystemRepresentation];
  if( !filePath )
  {
    handler( [NSError errorWithDomain:@"gov.sandia.SpecFilePreview"
                                code:1
                            userInfo:@{NSLocalizedDescriptionKey: @"Invalid file URL"}] );
    return;
  }

  CGImageRef cgImage = render_spec_file_to_cgimage( filePath, 800, 600, SpectrumPreview, NULL );

  if( !cgImage )
  {
    handler( [NSError errorWithDomain:@"gov.sandia.SpecFilePreview"
                                code:2
                            userInfo:@{NSLocalizedDescriptionKey: @"Failed to render spectrum"}] );
    return;
  }

  NSImage *nsImage = [[NSImage alloc] initWithCGImage:cgImage
                                                 size:NSMakeSize( 800, 600 )];
  CGImageRelease( cgImage );

  NSImageView *imageView = [NSImageView imageViewWithImage:nsImage];
  imageView.frame = self.view.bounds;
  imageView.autoresizingMask = NSViewWidthSizable | NSViewHeightSizable;
  imageView.imageScaling = NSImageScaleProportionallyUpOrDown;
  [self.view addSubview:imageView];

  handler( nil );
}

@end
