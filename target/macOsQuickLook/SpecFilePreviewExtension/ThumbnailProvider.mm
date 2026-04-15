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
#import <QuickLookThumbnailing/QuickLookThumbnailing.h>

#include "SpecPreviewCommon.h"


@interface ThumbnailProvider : QLThumbnailProvider
@end

@implementation ThumbnailProvider

- (void)provideThumbnailForFileRequest:(QLFileThumbnailRequest *)request
                     completionHandler:(void (^)(QLThumbnailReply * _Nullable, NSError * _Nullable))handler
{
  const CGSize maxSize = request.maximumSize;
  const CGFloat scale = request.scale;

  // Use points for rendering
  const float width_pt = (float)maxSize.width;
  const float height_pt = (float)maxSize.height;

  const char *filePath = [[request.fileURL path] fileSystemRepresentation];
  if( !filePath )
  {
    handler( nil, [NSError errorWithDomain:@"gov.sandia.SpecFileThumbnail"
                                     code:1
                                 userInfo:@{NSLocalizedDescriptionKey: @"Invalid file URL"}] );
    return;
  }

  // Skip thumbnail for large files (slow to parse) or very small sizes (not useful)
  // Returning handler(nil, nil) causes the system to use the default file icon.
  NSNumber *fileSize = nil;
  [request.fileURL getResourceValue:&fileSize forKey:NSURLFileSizeKey error:nil];
  const unsigned long long maxFileSize = 15 * 1024 * 1024; // 15 MB

  if( (fileSize && [fileSize unsignedLongLongValue] > maxFileSize)
     || width_pt < 32 || height_pt < 32 )
  {
    handler( nil, nil );
    return;
  }

  NSString *logoPath = [[NSBundle mainBundle] pathForResource:@"InterSpec128" ofType:@"png"];
  const char *logoCStr = logoPath ? [logoPath fileSystemRepresentation] : NULL;


#if MACOS_QUICK_LOOK_USE_PDF

  // Render spectrum to PDF
  uint8_t *pdfData = NULL;
  size_t pdfLen = 0;
  render_spec_file_to_pdf( &pdfData, &pdfLen, filePath, width_pt, height_pt, SpectrumThumbnail, logoCStr );

  if( !pdfData || pdfLen == 0 )
  {
    handler( nil, [NSError errorWithDomain:@"gov.sandia.SpecFileThumbnail"
                                     code:2
                                 userInfo:@{NSLocalizedDescriptionKey: @"Failed to render spectrum"}] );
    return;
  }

  // Create a CGPDFDocument from the PDF data
  NSData *nsData = [NSData dataWithBytesNoCopy:pdfData length:pdfLen freeWhenDone:YES];
  CGDataProviderRef dataProvider = CGDataProviderCreateWithCFData( (__bridge CFDataRef)nsData );

  if( !dataProvider )
  {
    handler( nil, [NSError errorWithDomain:@"gov.sandia.SpecFileThumbnail"
                                     code:3
                                 userInfo:@{NSLocalizedDescriptionKey: @"Failed to create PDF data provider"}] );
    return;
  }

  CGPDFDocumentRef pdfDoc = CGPDFDocumentCreateWithProvider( dataProvider );
  CGDataProviderRelease( dataProvider );

  if( !pdfDoc )
  {
    handler( nil, [NSError errorWithDomain:@"gov.sandia.SpecFileThumbnail"
                                     code:4
                                 userInfo:@{NSLocalizedDescriptionKey: @"Failed to create PDF document"}] );
    return;
  }

  CGPDFPageRef pdfPage = CGPDFDocumentGetPage( pdfDoc, 1 );
  if( !pdfPage )
  {
    CGPDFDocumentRelease( pdfDoc );
    handler( nil, [NSError errorWithDomain:@"gov.sandia.SpecFileThumbnail"
                                     code:5
                                 userInfo:@{NSLocalizedDescriptionKey: @"No pages in PDF"}] );
    return;
  }

  // Get the PDF page dimensions to compute the drawing rect
  const CGRect mediaBox = CGPDFPageGetBoxRect( pdfPage, kCGPDFMediaBox );
  const CGFloat pdfWidth = CGRectGetWidth( mediaBox );
  const CGFloat pdfHeight = CGRectGetHeight( mediaBox );

  // Compute the context size (maintain aspect ratio within maxSize)
  CGFloat contextWidth = maxSize.width;
  CGFloat contextHeight = maxSize.height;
  if( pdfWidth > 0 && pdfHeight > 0 )
  {
    const CGFloat aspect = pdfWidth / pdfHeight;
    if( contextWidth / contextHeight > aspect )
      contextWidth = contextHeight * aspect;
    else
      contextHeight = contextWidth / aspect;
  }

  // Retain pdfDoc for use in the drawing block; release it when done
  CGPDFDocumentRetain( pdfDoc );

  QLThumbnailReply *reply =
    [QLThumbnailReply replyWithContextSize:CGSizeMake( contextWidth, contextHeight )
                                 currentContextDrawingBlock:^BOOL {
      CGContextRef ctx = [[NSGraphicsContext currentContext] CGContext];
      if( !ctx )
      {
        CGPDFDocumentRelease( pdfDoc );
        return NO;
      }

      // Fill with white background
      CGContextSetRGBFillColor( ctx, 1.0, 1.0, 1.0, 1.0 );
      CGContextFillRect( ctx, CGRectMake( 0, 0, contextWidth, contextHeight ) );

      // Use CGPDFPageGetDrawingTransform to properly handle all PDF coordinate transforms
      const CGRect drawRect = CGRectMake( 0, 0, contextWidth, contextHeight );
      const CGAffineTransform pdfTransform = CGPDFPageGetDrawingTransform( pdfPage, kCGPDFMediaBox, drawRect, 0, true );
      CGContextConcatCTM( ctx, pdfTransform );

      CGContextDrawPDFPage( ctx, pdfPage );

      CGPDFDocumentRelease( pdfDoc );
      return YES;
    }];

  handler( reply, nil );


#else /* !MACOS_QUICK_LOOK_USE_PDF */

  CGImageRef cgImage = render_spec_file_to_cgimage( filePath, width_pt, height_pt,
                                                     SpectrumThumbnail, logoCStr );
  if( !cgImage )
  {
    handler( nil, [NSError errorWithDomain:@"gov.sandia.SpecFileThumbnail"
                                     code:2
                                 userInfo:@{NSLocalizedDescriptionKey: @"Failed to render spectrum"}] );
    return;
  }

  QLThumbnailReply *reply =
    [QLThumbnailReply replyWithContextSize:CGSizeMake( width_pt, height_pt )
                    currentContextDrawingBlock:^BOOL {
      CGContextRef ctx = [[NSGraphicsContext currentContext] CGContext];
      if( !ctx )
      {
        CGImageRelease( cgImage );
        return NO;
      }

      CGContextDrawImage( ctx, CGRectMake( 0, 0, width_pt, height_pt ), cgImage );
      CGImageRelease( cgImage );
      return YES;
    }];

  handler( reply, nil );

#endif /* MACOS_QUICK_LOOK_USE_PDF */
}

@end
