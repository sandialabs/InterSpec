#include <stdio.h>
#include <stdlib.h>

#include <QuickLook/QuickLook.h>
#include <CoreServices/CoreServices.h>
#include <CoreFoundation/CoreFoundation.h>

#import <Foundation/Foundation.h>
#import <AppKit/AppKit.h>

#include "SpecPreviewCommon.h"

#define MIN_SIZE(a,b) ((a)<(b)?(a):(b))

OSStatus GenerateThumbnailForURL(void *thisInterface, QLThumbnailRequestRef thumbnail, CFURLRef url, CFStringRef contentTypeUTI, CFDictionaryRef options, CGSize maxSize);
void CancelThumbnailGeneration(void *thisInterface, QLThumbnailRequestRef thumbnail);

/* -----------------------------------------------------------------------------
    Generate a thumbnail for file

   This function's job is to create thumbnail for designated file as fast as possible
   ----------------------------------------------------------------------------- */

OSStatus GenerateThumbnailForURL(void *thisInterface, QLThumbnailRequestRef thumbnail, CFURLRef url, CFStringRef contentTypeUTI, CFDictionaryRef options, CGSize maxSize)
{
  //To Run from the command line:
  //qlmanage -d1 -t /Users/wcjohns/rad_ana/InterSpec/example_spectra/ba133_source_640s_20100317.n42 -g ./SpecFilePreview.qlgenerator -c gov.sandia.gamma-spectrum
  
  printf( "\n\nIn SpecFileThumbnail\n" );
  @autoreleasepool {
  // To complete your generator please implement the function GenerateThumbnailForURL in GenerateThumbnailForURL.c
  CFStringRef filepath = CFURLCopyFileSystemPath(url, kCFURLPOSIXPathStyle);
  if( !filepath )
  {
    printf( "SpecFileThumbnail: Failed to get filepath!\n" );
    return 1;
  }
  
  char filename_buffer[2*1024];
  if( !CFStringGetCString(filepath, filename_buffer, sizeof(filename_buffer), kCFStringEncodingUTF8) )
  {
    printf( "SpecFileThumbnail: Failed to get filepath c-string!\n" );
    return 2;
  }
  
  printf( "SpecFileThumbnail: Got file %s\n", filename_buffer );
  
  //@autoreleasepool {
  // The above might have taken some time, so before proceeding make sure the user didn't cancel the request
  if( QLThumbnailRequestIsCancelled(thumbnail) )
  {
    printf( "SpecFileThumbnail: preview canceled 0\n" );
    return noErr;
  }
  
  //FILE *f = fopen ( "/Users/wcjohns/Downloads/len.txt", "a" );
  //fprintf(f, "Thumb: %i, %i\n", (int)maxSize.width, (int)maxSize.height );
  //fclose(f);
  
  //maxSize looks to be 128x128 points.
  //  Things like [[NSScreen mainScreen] convertRectToBacking: rect] seem to be pointless
  const int width = (int)maxSize.width;
  const int height = (int)maxSize.height;
  
  uint8_t *preview_data = NULL;
  size_t preview_data_len = 0;
  render_spec_file_to_preview( &preview_data, &preview_data_len, filename_buffer, width, height, SpectrumThumbnail );
  
  if( QLThumbnailRequestIsCancelled(thumbnail) )
  {
    printf( "SpecFileThumbnail: preview canceled 1\n" );
    free( preview_data );
    return noErr;
  }
  
  if( !preview_data )
  {
    printf( "SpecFileThumbnail: failed to grab spectrum data\n" );
    return noErr;
  }
  
  CFDataRef hmtl_dataref = CFDataCreate( kCFAllocatorDefault, (const UInt8 *)preview_data, preview_data_len );
  
  free( preview_data );
  
  printf( "SpecFileThumbnail: rendered chart\n" );
  
  NSDictionary *unused_properties = @{};
  
  // Put metadata and attachment in a dictionary
  NSDictionary *preview_properties = @{ // properties for the HTML data
    (__bridge NSString *)kQLPreviewPropertyTextEncodingNameKey : @"UTF-8",
    (__bridge NSString *)kQLPreviewPropertyMIMETypeKey : @"text/html"
    
    // properties for attaching the CSS stylesheet
    //(__bridge NSString *)kQLPreviewPropertyAttachmentsKey : @{
    //@"plistStylesheet.css" : @{
    //    (__bridge NSString *)kQLPreviewPropertyMIMETypeKey : @"text/css",
    //    (__bridge NSString *)kQLPreviewPropertyAttachmentDataKey: cssData,
    //    },
    //},
  };
  
  // Pass preview data and metadata/attachment dictionary to QuickLook
  //QL_EXPORT void QLThumbnailRequestSetImage(QLThumbnailRequestRef thumbnail, CGImageRef image, CFDictionaryRef properties);
  //QL_EXPORT void QLThumbnailRequestSetImageWithData(QLThumbnailRequestRef thumbnail, CFDataRef data, CFDictionaryRef properties);
  //I dont things this next call is valid: @discussion Currently supported UTIs are: none. This call only works if your generator is set to be run in the main thread
  //kUTTypePDF
  QLThumbnailRequestSetThumbnailWithDataRepresentation(thumbnail,
                                        hmtl_dataref,
                                        kUTTypeHTML,
                                        (__bridge CFDictionaryRef)preview_properties,
                                        (__bridge CFDictionaryRef)unused_properties
                                                       );
  }
  
  printf( "Done in SpecFileThumbnail\n\n\n" );
  
  return noErr;
}

void CancelThumbnailGeneration(void *thisInterface, QLThumbnailRequestRef thumbnail)
{
    // Implement only if supported
}
