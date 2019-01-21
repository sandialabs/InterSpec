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
//  Created by Johnson, William C on 12/25/17.
//

#ifndef SpecPreviewCommon_h
#define SpecPreviewCommon_h

#include <stddef.h>
#include <stdint.h>
/*
 Instructions for compiling:
   cd InterSpec
   mkdir build_quick_look
   cd build_quick_look
   cmake -DBUILD_AS_QUICK_LOOK_LIBRARY=ON ..
   make -j8
   cd ../target/macOsQuickLook/SpecFilePreview
   open SpecFilePreview.xcodeproj
   and then build
 
 TODO:
   -Create different UTTypeIdentifier types for each spectrum file format
     -Make “Conforms to UITs” indicate XML or text or whatever, as apropriate,
      so the Finder "Open With" menu will present other reasonable programs for opening
   -The Thumbnail preview sizing could be improved, as could the padding and
    text sizes individually for both preview and thumbnail
   -Should consider getting rid of y-axis numbers/title for true icon renders.
    -Need to verify, or figure out how to tell if a true icon render
   -Files with multiple spectra, but arent passthrough, currenly display as a
    single spectrum; it would be nice to allow changing spectra in the preview,
    like a PDF or Word document.  Implementing this would require generating
    the preview as a PDF... this would take some decent work.
     -See https://www.cocoanetics.com/2010/06/rendering-pdf-is-easier-than-you-thought/
      for a potential method of rendering SVG (as a webpage) to PDF.
     -Also see WPdfImage to render to PDF, and then push to apples stuff.
      Looks like WPdfImage can render to a specified page in a document.
 */

#include <Wt/WConfig.h>

/* If rendering as PDF, then Wt must be compiled with PDF (libharu) enabled, and
   you must link to libharu, libpng, and libz (in addition to Wt and boost).
   When using PDF, I am slowly working towards enabing multiple pages...
 */
#if( defined(WT_HAS_WPDFIMAGE) )
#define RENDER_PREVIEWS_AS_PDF 1
#else
#define RENDER_PREVIEWS_AS_PDF 0
#endif


#ifdef __cplusplus
extern "C" {
#endif
  enum SpecPreviewType
  {
    SpectrumThumbnail,
    SpectrumPreview
  };
  
  //Renders either the preview or thumbnail to a buffer.
  //  Thumbnails are always returned in HTML format (kUTTypeHTML).
  //  Previews are returned as PDFs if RENDER_PREVIEWS_AS_PDF is set to 1, or
  //    else they are returned as HTML.
  //
  //  (*result) will be set to buffer containing results.
  // (*result_size) will be set to result length.
  //  size will be zero on failure, and (*result) null.
  void render_spec_file_to_preview( uint8_t **result, size_t *result_size,
                                    const char * const filename, const float width_pt, const float height_pt, const enum SpecPreviewType type );
#ifdef __cplusplus
}
#endif

#endif /* SpecPreviewCommon_h */
