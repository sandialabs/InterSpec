#ifndef PlutoPaintDevice_h
#define PlutoPaintDevice_h
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

#include <string>

#include <Wt/WFont.h>
#include <Wt/WColor.h>
#include <Wt/WLength.h>
#include <Wt/WPaintDevice.h>

struct plutovg_canvas;
struct plutovg_surface;

namespace Wt
{
  class WAbstractDataInfo;
}

/** A `Wt::WPaintDevice` that draws into an in-memory image with plutovg, for rendering charts to
 PNG or JPEG with no browser and no platform graphics library.

 Text is drawn with one embedded font (DejaVu Sans), whatever the requested family or weight;
 plutovg does not shape complex scripts, so labels should be in a Latin script.

 Ported from the CoreGraphics device of the macOS Quick Look extension (CGPaintDevice).
 */
class PlutoPaintDevice : public Wt::WPaintDevice
{
public:
  /** Creates an image of the given size, filled with `background`. */
  PlutoPaintDevice( int width, int height, const Wt::WColor &background = Wt::WColor( 255, 255, 255 ) );
  ~PlutoPaintDevice();

  PlutoPaintDevice( const PlutoPaintDevice & ) = delete;
  PlutoPaintDevice &operator=( const PlutoPaintDevice & ) = delete;

  /** Returns the image encoded as "png" or "jpeg".
   @throws std::exception on failure, or for any other format.
   */
  std::string encode( const std::string &format ) const;

  // WPaintDevice interface
  Wt::WFlags<Wt::PaintDeviceFeatureFlag> features() const override;
  Wt::WLength width() const override;
  Wt::WLength height() const override;
  void setChanged( Wt::WFlags<Wt::PainterChangeFlag> flags ) override;
  void drawArc( const Wt::WRectF &rect, double startAngle, double spanAngle ) override;
  void drawImage( const Wt::WRectF &rect, const std::string &imageUri,
                  int imgWidth, int imgHeight, const Wt::WRectF &sourceRect ) override;
  void drawImage( const Wt::WRectF &rect, const Wt::WAbstractDataInfo *imageInfo,
                  int imgWidth, int imgHeight, const Wt::WRectF &sourceRect ) override;
  void drawLine( double x1, double y1, double x2, double y2 ) override;
  void drawPath( const Wt::WPainterPath &path ) override;
  void drawRect( const Wt::WRectF &rectangle ) override;
  void drawText( const Wt::WRectF &rect, Wt::WFlags<Wt::AlignmentFlag> flags,
                 Wt::TextFlag textFlag, const Wt::WTextF &text,
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
  /** Adds `path` to the canvas's current path. */
  void addPath( const Wt::WPainterPath &path );

  /** Fills (with the brush) and strokes (with the pen) the current path, then clears it. */
  void fillAndStroke();

  /** Sets the canvas paint to `color`. */
  void setColor( const Wt::WColor &color );

  /** Sets the canvas transform to `t`. */
  void setTransform( const Wt::WTransform &t );

  /** The font size, in pixels, for the painter's current font. */
  double fontSize() const;

  const int m_width;
  const int m_height;
  Wt::WPainter *m_painter;
  plutovg_surface *m_surface;
  plutovg_canvas *m_canvas;
};//class PlutoPaintDevice

#endif //PlutoPaintDevice_h
