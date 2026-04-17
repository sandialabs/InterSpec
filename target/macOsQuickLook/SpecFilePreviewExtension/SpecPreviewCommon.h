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

#ifndef SpecPreviewCommon_h
#define SpecPreviewCommon_h

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

  enum SpecPreviewType
  {
    SpectrumPreview,
    SpectrumThumbnail
  };


#if MACOS_QUICK_LOOK_USE_PDF

  /** Renders a spectrum file to PDF for use as a Quick Look preview or thumbnail.

   @param result Will be set to a malloc'd buffer containing PDF data. Caller must free().
   @param result_size Will be set to the byte length of the result buffer. Zero on failure.
   @param filename Path to the spectrum file to render.
   @param width_pt Width in points for the rendered image.
   @param height_pt Height in points for the rendered image.
   @param type SpectrumPreview for full preview (axes, titles, time series for passthrough),
          or SpectrumThumbnail for compact rendering (no axes/titles, spectrum only).
   @param logo_path Path to logo PNG for thumbnail watermark, or NULL to skip.
   */
  void render_spec_file_to_pdf( uint8_t **result, size_t *result_size,
                                const char * const filename,
                                const float width_pt, const float height_pt,
                                const enum SpecPreviewType type,
                                const char * const logo_path );

#else /* !MACOS_QUICK_LOOK_USE_PDF */

#include <CoreGraphics/CoreGraphics.h>

  /** Renders a spectrum file to a CGImage for Quick Look preview or thumbnail.

   Uses native Core Graphics rendering (no libharu dependency).

   @param filename Path to the spectrum file to render.
   @param width Width in pixels for the rendered image.
   @param height Height in pixels for the rendered image.
   @param type SpectrumPreview for full preview (axes, titles, time series for passthrough),
          or SpectrumThumbnail for compact rendering (no axes/titles, spectrum only).
   @param logo_path Path to logo PNG for thumbnail watermark, or NULL to skip.
   @returns A CGImageRef on success (caller must CGImageRelease), or NULL on failure.
   */
  CGImageRef render_spec_file_to_cgimage( const char * const filename,
                                          const float width, const float height,
                                          const enum SpecPreviewType type,
                                          const char * const logo_path );

#endif /* MACOS_QUICK_LOOK_USE_PDF */

#ifdef __cplusplus
}
#endif

#endif /* SpecPreviewCommon_h */
