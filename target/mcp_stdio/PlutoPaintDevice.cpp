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

#include "PlutoPaintDevice.h"

#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <stdexcept>

#include <plutovg.h>

#include <Wt/WPen.h>
#include <Wt/WBrush.h>
#include <Wt/WRectF.h>
#include <Wt/WPointF.h>
#include <Wt/WPainter.h>
#include <Wt/WTransform.h>
#include <Wt/WPainterPath.h>
#include <Wt/WFontMetrics.h>

using namespace Wt;

// The TTF compiled into the executable (see embed_file.cmake).
extern const unsigned char g_embedded_font_ttf[];
extern const size_t g_embedded_font_ttf_size;


namespace
{
  const double ns_pi = 3.14159265358979323846;

  /** The embedded font, loaded once; plutovg font faces may be shared between canvases. */
  plutovg_font_face_t *font_face()
  {
    static plutovg_font_face_t * const face
          = plutovg_font_face_load_from_data( g_embedded_font_ttf,
                                              static_cast<unsigned int>( g_embedded_font_ttf_size ),
                                              0, nullptr, nullptr );
    return face;
  }

  /** Adds an elliptic arc to the canvas path, continuing from its current point.  Angles are in
   degrees, counter-clockwise, with y down (as WSvgImage maps them); drawn as short line segments,
   which handles ellipses and either direction.
   */
  void add_elliptic_arc( plutovg_canvas_t *canvas, const double cx, const double cy,
                         const double rx, const double ry, const double start_deg, const double sweep_deg )
  {
    const double start = -start_deg * ns_pi / 180.0;
    double sweep = -sweep_deg * ns_pi / 180.0;
    if( std::fabs( sweep ) > 2.0 * ns_pi )
      sweep = std::copysign( 2.0 * ns_pi, sweep );

    // One line segment per two degrees.
    const int nsteps = std::max( 8, static_cast<int>( std::ceil( std::fabs( sweep ) / (ns_pi / 90.0) ) ) );
    for( int step = 0; step <= nsteps; ++step )
    {
      const double angle = start + sweep * step / nsteps;
      // (A line_to with no current point is a move_to.)
      plutovg_canvas_line_to( canvas, static_cast<float>( cx + rx * std::cos( angle ) ),
                              static_cast<float>( cy + ry * std::sin( angle ) ) );
    }
  }//add_elliptic_arc(...)


  /** Appends what plutovg's encoders write to a std::string. */
  void append_to_string( void *closure, void *data, int size )
  {
    static_cast<std::string *>( closure )->append( static_cast<const char *>( data ), size );
  }
}//namespace


PlutoPaintDevice::PlutoPaintDevice( const int width, const int height, const WColor &background )
  : m_width( width ),
    m_height( height ),
    m_painter( nullptr ),
    m_surface( nullptr ),
    m_canvas( nullptr )
{
  if( (width <= 0) || (height <= 0) )
    throw std::runtime_error( "PlutoPaintDevice: invalid image size" );

  m_surface = plutovg_surface_create( width, height );
  if( !m_surface )
    throw std::runtime_error( "PlutoPaintDevice: could not create the image" );

  plutovg_color_t color;
  plutovg_color_init_rgba8( &color, background.red(), background.green(), background.blue(), background.alpha() );
  plutovg_surface_clear( m_surface, &color );

  m_canvas = plutovg_canvas_create( m_surface );
  if( !m_canvas )
  {
    plutovg_surface_destroy( m_surface );
    throw std::runtime_error( "PlutoPaintDevice: could not create the canvas" );
  }
}//PlutoPaintDevice constructor


PlutoPaintDevice::~PlutoPaintDevice()
{
  plutovg_canvas_destroy( m_canvas );
  plutovg_surface_destroy( m_surface );
}


std::string PlutoPaintDevice::encode( const std::string &format ) const
{
  std::string data;
  bool ok = false;
  if( format == "png" )
    ok = plutovg_surface_write_to_png_stream( m_surface, &append_to_string, &data );
  else if( (format == "jpeg") || (format == "jpg") )
    ok = plutovg_surface_write_to_jpg_stream( m_surface, &append_to_string, &data, 90 );
  else
    throw std::runtime_error( "PlutoPaintDevice: unsupported image format '" + format + "'" );

  if( !ok || data.empty() )
    throw std::runtime_error( "PlutoPaintDevice: failed to encode the image as " + format );

  return data;
}//encode(...)


WFlags<PaintDeviceFeatureFlag> PlutoPaintDevice::features() const
{
  return PaintDeviceFeatureFlag::FontMetrics;
}


WLength PlutoPaintDevice::width() const
{
  return WLength( m_width, LengthUnit::Pixel );
}


WLength PlutoPaintDevice::height() const
{
  return WLength( m_height, LengthUnit::Pixel );
}


bool PlutoPaintDevice::paintActive() const
{
  return (m_painter != nullptr);
}


WPainter *PlutoPaintDevice::painter() const
{
  return m_painter;
}


void PlutoPaintDevice::setPainter( WPainter *painter )
{
  m_painter = painter;
}


void PlutoPaintDevice::init()
{
  // The base state; setChanged() restores to it, and re-applies the clip and transform.
  plutovg_canvas_save( m_canvas );
}


void PlutoPaintDevice::done()
{
  plutovg_canvas_restore( m_canvas );
}


void PlutoPaintDevice::setTransform( const WTransform &t )
{
  plutovg_matrix_t matrix;
  plutovg_matrix_init( &matrix, static_cast<float>( t.m11() ), static_cast<float>( t.m12() ),
                       static_cast<float>( t.m21() ), static_cast<float>( t.m22() ),
                       static_cast<float>( t.dx() ), static_cast<float>( t.dy() ) );
  plutovg_canvas_set_matrix( m_canvas, &matrix );
}


void PlutoPaintDevice::setColor( const WColor &c )
{
  plutovg_canvas_set_rgba( m_canvas, c.red() / 255.0f, c.green() / 255.0f, c.blue() / 255.0f, c.alpha() / 255.0f );
}


double PlutoPaintDevice::fontSize() const
{
  return m_painter ? m_painter->font().sizeLength( 16 ).toPixels() : 16.0;
}


void PlutoPaintDevice::setChanged( WFlags<PainterChangeFlag> flags )
{
  // Colors are set before each drawing operation (plutovg has a single paint for filling and
  //  stroking), and the font when drawing text; only the clip, transform and pen geometry are
  //  canvas state.
  if( flags.test( PainterChangeFlag::Transform ) || flags.test( PainterChangeFlag::Clipping ) )
  {
    plutovg_canvas_restore( m_canvas );
    plutovg_canvas_save( m_canvas );

    if( m_painter->hasClipping() && !m_painter->clipPath().isEmpty() )
    {
      setTransform( m_painter->clipPathTransform() );
      plutovg_canvas_new_path( m_canvas );
      addPath( m_painter->clipPath() );
      plutovg_canvas_clip( m_canvas );
    }

    setTransform( m_painter->combinedTransform() );
    flags = PainterChangeFlag::Pen;
  }//if( transform or clipping changed )

  if( flags.test( PainterChangeFlag::Pen ) )
  {
    const WPen &pen = m_painter->pen();
    const double pen_width = pen.width().toPixels();
    plutovg_canvas_set_line_width( m_canvas, static_cast<float>( (pen_width == 0.0) ? 1.0 : pen_width ) );

    switch( pen.capStyle() )
    {
      case PenCapStyle::Flat:   plutovg_canvas_set_line_cap( m_canvas, PLUTOVG_LINE_CAP_BUTT );   break;
      case PenCapStyle::Square: plutovg_canvas_set_line_cap( m_canvas, PLUTOVG_LINE_CAP_SQUARE ); break;
      case PenCapStyle::Round:  plutovg_canvas_set_line_cap( m_canvas, PLUTOVG_LINE_CAP_ROUND );  break;
    }

    switch( pen.joinStyle() )
    {
      case PenJoinStyle::Miter: plutovg_canvas_set_line_join( m_canvas, PLUTOVG_LINE_JOIN_MITER ); break;
      case PenJoinStyle::Bevel: plutovg_canvas_set_line_join( m_canvas, PLUTOVG_LINE_JOIN_BEVEL ); break;
      case PenJoinStyle::Round: plutovg_canvas_set_line_join( m_canvas, PLUTOVG_LINE_JOIN_ROUND ); break;
    }

    // The same dash patterns as the other InterSpec paint devices.
    std::vector<float> dashes;
    switch( pen.style() )
    {
      case PenStyle::None:
      case PenStyle::SolidLine:
        break;
      case PenStyle::DashLine:       dashes = { 4.0f, 2.0f };                         break;
      case PenStyle::DotLine:        dashes = { 1.0f, 2.0f };                         break;
      case PenStyle::DashDotLine:    dashes = { 4.0f, 2.0f, 1.0f, 2.0f };             break;
      case PenStyle::DashDotDotLine: dashes = { 4.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f }; break;
    }
    plutovg_canvas_set_dash( m_canvas, 0.0f, dashes.empty() ? nullptr : dashes.data(),
                             static_cast<int>( dashes.size() ) );
  }//if( pen changed )
}//setChanged(...)


void PlutoPaintDevice::addPath( const WPainterPath &path )
{
  const std::vector<WPainterPath::Segment> &segments = path.segments();
  if( segments.empty() )
    return;

  if( segments[0].type() != MoveTo )
    plutovg_canvas_move_to( m_canvas, 0.0f, 0.0f );

  for( size_t i = 0; i < segments.size(); ++i )
  {
    const WPainterPath::Segment &s = segments[i];

    switch( s.type() )
    {
      case MoveTo:
        plutovg_canvas_move_to( m_canvas, static_cast<float>( s.x() ), static_cast<float>( s.y() ) );
        break;

      case LineTo:
        plutovg_canvas_line_to( m_canvas, static_cast<float>( s.x() ), static_cast<float>( s.y() ) );
        break;

      case CubicC1:
        assert( (i + 2) < segments.size() );
        plutovg_canvas_cubic_to( m_canvas,
                                 static_cast<float>( s.x() ), static_cast<float>( s.y() ),
                                 static_cast<float>( segments[i+1].x() ), static_cast<float>( segments[i+1].y() ),
                                 static_cast<float>( segments[i+2].x() ), static_cast<float>( segments[i+2].y() ) );
        i += 2;
        break;

      case QuadC:
        assert( (i + 1) < segments.size() );
        plutovg_canvas_quad_to( m_canvas,
                                static_cast<float>( s.x() ), static_cast<float>( s.y() ),
                                static_cast<float>( segments[i+1].x() ), static_cast<float>( segments[i+1].y() ) );
        i += 1;
        break;

      case ArcC:
        assert( (i + 2) < segments.size() );
        add_elliptic_arc( m_canvas, s.x(), s.y(), segments[i+1].x(), segments[i+1].y(),
                          segments[i+2].x(), segments[i+2].y() );
        i += 2;
        break;

      case CubicC2:
      case CubicEnd:
      case QuadEnd:
      case ArcR:
      case ArcAngleSweep:
        assert( 0 );  // consumed with the segment that starts them
        break;
    }//switch( s.type() )
  }//for( loop over segments )
}//addPath(...)


void PlutoPaintDevice::fillAndStroke()
{
  const bool fill = (m_painter->brush().style() != BrushStyle::None);
  const bool stroke = (m_painter->pen().style() != PenStyle::None);

  if( fill )
  {
    setColor( m_painter->brush().color() );
    plutovg_canvas_fill_preserve( m_canvas );
  }

  if( stroke )
  {
    setColor( m_painter->pen().color() );
    plutovg_canvas_stroke_preserve( m_canvas );
  }

  plutovg_canvas_new_path( m_canvas );
}//fillAndStroke()


void PlutoPaintDevice::drawPath( const WPainterPath &path )
{
  plutovg_canvas_new_path( m_canvas );
  addPath( path );
  fillAndStroke();
}


void PlutoPaintDevice::drawLine( double x1, double y1, double x2, double y2 )
{
  if( m_painter->pen().style() == PenStyle::None )
    return;

  plutovg_canvas_new_path( m_canvas );
  plutovg_canvas_move_to( m_canvas, static_cast<float>( x1 ), static_cast<float>( y1 ) );
  plutovg_canvas_line_to( m_canvas, static_cast<float>( x2 ), static_cast<float>( y2 ) );
  setColor( m_painter->pen().color() );
  plutovg_canvas_stroke( m_canvas );
}


void PlutoPaintDevice::drawRect( const WRectF &rectangle )
{
  WPainterPath path;
  path.addRect( rectangle );
  drawPath( path );
}


void PlutoPaintDevice::drawArc( const WRectF &rect, double startAngle, double spanAngle )
{
  plutovg_canvas_new_path( m_canvas );
  add_elliptic_arc( m_canvas, rect.center().x(), rect.center().y(), rect.width() / 2.0, rect.height() / 2.0,
                    startAngle, spanAngle );
  fillAndStroke();
}


void PlutoPaintDevice::drawImage( const WRectF &rect, const std::string &imageUri,
                                  int imgWidth, int imgHeight, const WRectF &sourceRect )
{
  // Only local files (PNG/JPEG) are supported.
  plutovg_surface_t * const image = plutovg_surface_load_from_image_file( imageUri.c_str() );
  if( !image || (sourceRect.width() <= 0.0) || (sourceRect.height() <= 0.0) )
  {
    if( image )
      plutovg_surface_destroy( image );
    return;
  }

  // Map the source rectangle of the image onto `rect`.
  const double sx = rect.width() / sourceRect.width();
  const double sy = rect.height() / sourceRect.height();
  plutovg_matrix_t texture_matrix;
  plutovg_matrix_init( &texture_matrix, static_cast<float>( sx ), 0.0f, 0.0f, static_cast<float>( sy ),
                       static_cast<float>( rect.x() - sourceRect.x() * sx ),
                       static_cast<float>( rect.y() - sourceRect.y() * sy ) );

  plutovg_canvas_set_texture( m_canvas, image, PLUTOVG_TEXTURE_TYPE_PLAIN, 1.0f, &texture_matrix );
  plutovg_canvas_new_path( m_canvas );
  plutovg_canvas_rect( m_canvas, static_cast<float>( rect.x() ), static_cast<float>( rect.y() ),
                       static_cast<float>( rect.width() ), static_cast<float>( rect.height() ) );
  plutovg_canvas_fill( m_canvas );
  plutovg_surface_destroy( image );
}//drawImage(...)


void PlutoPaintDevice::drawImage( const WRectF &, const WAbstractDataInfo *, int, int, const WRectF & )
{
  // Not used by InterSpec's charts.
}


void PlutoPaintDevice::drawText( const WRectF &rect, WFlags<AlignmentFlag> flags, TextFlag,
                                 const WTextF &text, const WPointF *clipPoint )
{
  if( clipPoint && m_painter && !m_painter->clipPath().isEmpty() )
  {
    if( !m_painter->clipPathTransform().map( m_painter->clipPath() )
          .isPointInPath( m_painter->worldTransform().map( *clipPoint ) ) )
      return;
  }

  plutovg_font_face_t * const face = font_face();
  const std::string utf8 = text.text().toUTF8();
  if( !face || utf8.empty() )
    return;

  const float size = static_cast<float>( fontSize() );
  float ascent = 0.0f, descent = 0.0f, line_gap = 0.0f;
  plutovg_font_face_get_metrics( face, size, &ascent, &descent, &line_gap, nullptr );
  descent = std::fabs( descent );  // negative, as in the font file

  const float text_width = plutovg_font_face_text_extents( face, size, utf8.data(), static_cast<int>( utf8.size() ),
                                                           PLUTOVG_TEXT_ENCODING_UTF8, nullptr );

  const WFlags<AlignmentFlag> horizontal = flags & AlignHorizontalMask;
  double x = rect.left();
  if( horizontal.test( AlignmentFlag::Right ) )
    x = rect.right() - text_width;
  else if( horizontal.test( AlignmentFlag::Center ) )
    x = rect.center().x() - text_width / 2.0;

  const WFlags<AlignmentFlag> vertical = flags & AlignVerticalMask;
  double baseline = rect.top() + ascent;
  if( vertical.test( AlignmentFlag::Middle ) )
    baseline = rect.center().y() + (ascent - descent) / 2.0;
  else if( vertical.test( AlignmentFlag::Bottom ) )
    baseline = rect.bottom() - descent;

  setColor( m_painter->pen().color() );
  plutovg_canvas_set_font( m_canvas, face, size );
  plutovg_canvas_fill_text( m_canvas, utf8.data(), static_cast<int>( utf8.size() ), PLUTOVG_TEXT_ENCODING_UTF8,
                            static_cast<float>( x ), static_cast<float>( baseline ) );
}//drawText(...)


WTextItem PlutoPaintDevice::measureText( const WString &text, double maxWidth, bool wordWrap )
{
  plutovg_font_face_t * const face = font_face();
  const std::string utf8 = text.toUTF8();
  if( !face || utf8.empty() )
    return WTextItem( text, 0 );

  const float size = static_cast<float>( fontSize() );
  const auto width_of = [face, size]( const std::string &s ) -> double {
    return plutovg_font_face_text_extents( face, size, s.data(), static_cast<int>( s.size() ),
                                           PLUTOVG_TEXT_ENCODING_UTF8, nullptr );
  };

  const double full_width = width_of( utf8 );
  if( !wordWrap || (maxWidth < 0.0) || (full_width <= maxWidth) )
    return WTextItem( text, full_width );

  // Break at the last space that fits.
  size_t break_pos = std::string::npos;
  for( size_t pos = utf8.find( ' ' ); pos != std::string::npos; pos = utf8.find( ' ', pos + 1 ) )
  {
    if( width_of( utf8.substr( 0, pos ) ) > maxWidth )
      break;
    break_pos = pos;
  }

  if( break_pos == std::string::npos )
    return WTextItem( text, full_width );

  const std::string first = utf8.substr( 0, break_pos );
  return WTextItem( WString::fromUTF8( first ), width_of( first ), width_of( utf8.substr( break_pos + 1 ) ) );
}//measureText(...)


WFontMetrics PlutoPaintDevice::fontMetrics()
{
  const WFont font = m_painter ? m_painter->font() : WFont();
  plutovg_font_face_t * const face = font_face();
  if( !face )
    return WFontMetrics( font, 0, 0, 0 );

  float ascent = 0.0f, descent = 0.0f, line_gap = 0.0f;
  plutovg_font_face_get_metrics( face, static_cast<float>( fontSize() ), &ascent, &descent, &line_gap, nullptr );
  return WFontMetrics( font, line_gap, ascent, std::fabs( descent ) );
}//fontMetrics()
