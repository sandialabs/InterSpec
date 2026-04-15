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
#import <Quartz/Quartz.h>

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

  // Render spectrum to PDF using the shared rendering code
  uint8_t *pdfData = NULL;
  size_t pdfLen = 0;
  render_spec_file_to_pdf( &pdfData, &pdfLen, filePath, 800, 600, SpectrumPreview, NULL );

  if( !pdfData || pdfLen == 0 )
  {
    handler( [NSError errorWithDomain:@"gov.sandia.SpecFilePreview"
                                code:2
                            userInfo:@{NSLocalizedDescriptionKey: @"Failed to render spectrum"}] );
    return;
  }

  // Create PDFDocument from rendered data
  NSData *nsData = [NSData dataWithBytesNoCopy:pdfData length:pdfLen freeWhenDone:YES];
  PDFDocument *pdfDoc = [[PDFDocument alloc] initWithData:nsData];

  if( !pdfDoc )
  {
    handler( [NSError errorWithDomain:@"gov.sandia.SpecFilePreview"
                                code:3
                            userInfo:@{NSLocalizedDescriptionKey: @"Failed to create PDF document"}] );
    return;
  }

  // Display in a PDFView (native, works in Quick Look view bridge)
  PDFView *pdfView = [[PDFView alloc] initWithFrame:self.view.bounds];
  pdfView.autoresizingMask = NSViewWidthSizable | NSViewHeightSizable;
  pdfView.autoScales = YES;
  pdfView.document = pdfDoc;
  pdfView.displayMode = kPDFDisplaySinglePage;
  [self.view addSubview:pdfView];

  handler( nil );
}

@end
