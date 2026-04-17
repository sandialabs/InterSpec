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

#ifndef CGPaintDevice_h
#define CGPaintDevice_h

#include <CoreGraphics/CoreGraphics.h>
#include <CoreText/CoreText.h>

#include <Wt/WFont>
#include <Wt/WLength>
#include <Wt/WPaintDevice>

/** A Wt::WPaintDevice backed by a macOS CGBitmapContext and Core Text.

 Provides raster rendering of Wt painter output using native macOS graphics,
 eliminating the need for libharu/WPdfImage for bitmap output.

 Two construction modes:
 - Primary: creates and owns a CGBitmapContext of the given dimensions.
 - Sub-region: shares a parent device's context, painting into a translated sub-region.
   Used for passthrough spectra where the spectrum and time series occupy different
   vertical portions of the same image.
 */
class CGPaintDevice : public Wt::WPaintDevice
{
public:
  /** Creates a device that owns a CGBitmapContext of the given pixel dimensions. */
  CGPaintDevice( double width, double height );

  /** Creates a sub-region device that paints into a portion of the parent's context.
   The sub-region has local coordinates (0,0)-(width,height) and is placed at
   offset (x,y) within the parent's coordinate space.
   The parent must outlive this device.
   */
  CGPaintDevice( CGPaintDevice &parent, double x, double y, double width, double height );

  ~CGPaintDevice();

  // Non-copyable
  CGPaintDevice( const CGPaintDevice & ) = delete;
  CGPaintDevice &operator=( const CGPaintDevice & ) = delete;

  /** Returns a CGImage from the bitmap context. Caller must CGImageRelease().
   Only valid for the primary (owning) constructor.
   */
  CGImageRef createCGImage() const;

  /** Returns the underlying CGContextRef. */
  CGContextRef cgContext() const;

  // -- WPaintDevice interface --
  Wt::WFlags<FeatureFlag> features() const override;
  Wt::WLength width() const override;
  Wt::WLength height() const override;
  void setChanged( Wt::WFlags<ChangeFlag> flags ) override;
  void drawArc( const Wt::WRectF &rect, double startAngle, double spanAngle ) override;
  void drawImage( const Wt::WRectF &rect, const std::string &imageUri,
                  int imgWidth, int imgHeight, const Wt::WRectF &sourceRect ) override;
  void drawLine( double x1, double y1, double x2, double y2 ) override;
  void drawPath( const Wt::WPainterPath &path ) override;
  void drawText( const Wt::WRectF &rect, Wt::WFlags<Wt::AlignmentFlag> flags,
                 Wt::TextFlag textFlag, const Wt::WString &text,
                 const Wt::WPointF *clipPoint ) override;
  Wt::WTextItem measureText( const Wt::WString &text, double maxWidth = -1,
                              bool wordWrap = false ) override;
  Wt::WFontMetrics fontMetrics() override;
  void init() override;
  void done() override;
  bool paintActive() const override;

protected:
  Wt::WPainter *painter() const override;
  void setPainter( Wt::WPainter *painter ) override;

private:
  Wt::WLength m_width;
  Wt::WLength m_height;
  Wt::WPainter *m_painter;
  CGContextRef m_ctx;
  bool m_ownsContext;
  double m_originX;
  double m_originY;
  double m_bitmapHeight; // Total bitmap height (for computing Y-flip in sub-regions)

  // Current Core Text font
  CTFontRef m_ctFont;

  void applyTransform( const Wt::WTransform &t );
  void drawPlainPath( const Wt::WPainterPath &path );
  void strokeAndFillPath();
  CTFontRef createCTFont( const Wt::WFont &font );
  void updateFont( const Wt::WFont &font );
};//class CGPaintDevice

#endif /* CGPaintDevice_h */
