#include <CoreFoundation/CoreFoundation.h>
#include <CoreServices/CoreServices.h>
#import <Foundation/Foundation.h>
#include <QuickLook/QuickLook.h>


#include "SpecPreviewCommon.h"

OSStatus GeneratePreviewForURL(void *thisInterface, QLPreviewRequestRef preview, CFURLRef url, CFStringRef contentTypeUTI, CFDictionaryRef options);
void CancelPreviewGeneration(void *thisInterface, QLPreviewRequestRef preview);

/* -----------------------------------------------------------------------------
   Generate a preview for file

   This function's job is to create preview for designated file
   ----------------------------------------------------------------------------- */

OSStatus GeneratePreviewForURL(void *thisInterface, QLPreviewRequestRef preview, CFURLRef url, CFStringRef contentTypeUTI, CFDictionaryRef options)
{
  //sudo qlmanage -d4 -g ./SpecFilePreview.qlgenerator -c gov.sandia.gamma-spectrum -t /Users/wcjohns/rad_ana/InterSpec/example_spectra/ba133_source_640s_20100317.n42
  
  //https://developer.apple.com/library/content/documentation/UserExperience/Conceptual/Quicklook_Programming_Guide/Articles/QLDynamicGeneration.html#//apple_ref/doc/uid/TP40005020-CH15-SW5
  
  //Could use libRSVG (https://wiki.gnome.org/action/show/Projects/LibRsvg) to render SVG to pdf, to support multiple pages
  //Could just create html with svg embeded.
  
  printf( "\n\nIn SpecFilePreview\n" );
  @autoreleasepool {
  // To complete your generator please implement the function GenerateThumbnailForURL in GenerateThumbnailForURL.c
  CFStringRef filepath = CFURLCopyFileSystemPath(url, kCFURLPOSIXPathStyle);
  if( !filepath )
  {
    printf( "SpecFile Preview: Failed to get filepath!\n" );
    return 1;
  }
  
  char filename_buffer[2*1024];
  if( !CFStringGetCString(filepath, filename_buffer, sizeof(filename_buffer), kCFStringEncodingUTF8) )
  {
    printf( "SpecFilePreview: Failed to get filepath c-string!\n" );
    return 2;
  }
  
  printf( "SpecFilePreview: Got file %s\n", filename_buffer );
  
  //NSURL *nsurl = (__bridge NSURL *)url;
  //NSData *data = [NSData dataWithContentsOfURL:nsurl];
  //if (!data) return noErr;
  
    // Load the property list from the URL
  
  printf( "SpecFilePreview: will check if cancelled\n" );
    // The above might have taken some time, so before proceeding make sure the user didn't cancel the request
  if( QLPreviewRequestIsCancelled(preview) )
  {
    printf( "SpecFilePreview: preview canceled 0\n" );
    return noErr;
  }
   
  //const char *html = "Root";
  //CFDataRef hmtl_dataref = CFDataCreate( kCFAllocatorDefault, html, strlen(html) );
  //CFDataRef hmtl_dataref = CFDataCreateWithBytesNoCopy( NULL, html, strlen(html), kCFAllocatorNull );
  
    // Load a CSS stylesheet to attach to the HTML
    //NSBundle *bundle = [NSBundle bundleForClass:[HTMLPreviewBuilder class]];
    //NSURL *cssFile = [bundle URLForResource:@"plistStylesheet" withExtension:@"css"];
    //NSData *cssData = [NSData dataWithContentsOfURL:cssFile];
  //NSURL *furl = [NSURL fileURLWithPath: @"/Users/wcjohns/Downloads/modular_RPM_diagram_ver89.svg" ];
  //NSURL *furl = [NSURL fileURLWithPath: @"/Users/wcjohns/Library/Developer/Xcode/DerivedData/SpecFilePreview-aktaslzgzzavbrdhhxiymujbzepb/Build/Products/Release/example.html" ];
  //NSData *hmtl_ptr = [NSData dataWithContentsOfURL: furl ];
  //{
    //NSString* newStr = [NSString stringWithUTF8String:[hmtl_ptr bytes]];
    //const char *cString = [newStr UTF8String];
    //printf( "Contents: %s\n\n", cString);
  //}
  
    //FILE *f = fopen ( "/Users/wcjohns/Downloads/len.txt", "a" );
    //fprintf(f, "Preview\n" );
    //fclose(f);
    
    printf( "SpecFilePreview: about to render to preview\n" );

    uint8_t *preview_data = NULL;
    size_t preview_data_length = 0;
    render_spec_file_to_preview( &preview_data, &preview_data_length, filename_buffer, 800, 600, SpectrumPreview );
  
    if( !preview_data || preview_data_length==0 )
    {
      printf( "SpecFilePreview: failed to grab spectrum data\n" );
      return noErr;
    }
    
    if( QLPreviewRequestIsCancelled(preview) )
    {
      printf( "SpecFilePreview: preview canceled 1\n" );
      free( preview_data );
      return noErr;
    }
  
    printf( "SpecFilePreview: done rendering to preview\n" );
  
    CFDataRef hmtl_dataref = CFDataCreate( kCFAllocatorDefault, (const UInt8 *)preview_data, preview_data_length );
  
    free( preview_data );
  
  printf( "SpecFilePreview: rendered chart" );
#if( RENDER_PREVIEWS_AS_PDF )
    NSDictionary *properties = @{ };
    QLPreviewRequestSetDataRepresentation(preview, hmtl_dataref, kUTTypePDF, (__bridge CFDictionaryRef)properties);
#else
  // Put metadata and attachment in a dictionary
  NSDictionary *properties = @{ // properties for the HTML data
                                 (__bridge NSString *)kQLPreviewPropertyTextEncodingNameKey : @"UTF-8",
                                 (__bridge NSString *)kQLPreviewPropertyMIMETypeKey : @"text/html",
                                 };
  
  //Currently supported UTIs are: kUTTypeImage, kUTTypePDF, kUTTypeHTML,
  //            kUTTypeXML, kUTTypePlainText, kUTTypeRTF, kUTTypeMovie, kUTTypeAudio,
  
  // Pass preview data and metadata/attachment dictionary to QuickLook
  QLPreviewRequestSetDataRepresentation(preview,
                                          hmtl_dataref,
                                          //(__bridge CFDataRef)hmtl_ptr,
                                          kUTTypeHTML,
                                          (__bridge CFDictionaryRef)properties);
#endif
  }
  
  //Could do multiple pages in a PDF context:
  //CGContextRef pdfcontext = QLPreviewRequestCreatePDFContext(QLPreviewRequestRef preview, const CGRect *mediaBox, CFDictionaryRef auxiliaryInfo, CFDictionaryRef properties);
  //CGPDFContextBeginPage(pdfcontext, CFDictionaryRef pageInfo);
  //...
  //CGPDFContextEndPage(pdfcontext);
  //CGPDFContextBeginPage(pdfcontext, CFDictionaryRef pageInfo);
  //...
  //CGPDFContextEndPage(pdfcontext);
  //...

  
  printf( "Done in SpecFilePreview\n\n\n" );
  
  return noErr;
}

void CancelPreviewGeneration(void *thisInterface, QLPreviewRequestRef preview)
{
    // Implement only if supported
}
