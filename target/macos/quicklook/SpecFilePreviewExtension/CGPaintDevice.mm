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

#include "CGPaintDevice.h"

#include <cmath>
#include <cassert>
#include <string>
#include <iostream>

#include <ImageIO/ImageIO.h>

#include <Wt/WPen.h>
#include <Wt/WBrush.h>
#include <Wt/WColor.h>
#include <Wt/WRectF.h>
#include <Wt/WPointF.h>
#include <Wt/WPainter.h>
#include <Wt/WTransform.h>
#include <Wt/WPainterPath.h>
#include <Wt/WFontMetrics.h>
#include <Wt/WAbstractDataInfo.h>

using namespace Wt;


CGPaintDevice::CGPaintDevice( double width, double height )
  : m_width( WLength( width, LengthUnit::Pixel ) )
  , m_height( WLength( height, LengthUnit::Pixel ) )
  , m_painter( nullptr )
  , m_ctx( nullptr )
  , m_ownsContext( true )
  , m_originX( 0.0 )
  , m_originY( 0.0 )
  , m_bitmapHeight( height )
  , m_ctFont( nullptr )
{
  const size_t w = static_cast<size_t>( std::ceil( width ) );
  const size_t h = static_cast<size_t>( std::ceil( height ) );

  CGColorSpaceRef colorSpace = CGColorSpaceCreateWithName( kCGColorSpaceSRGB );
  m_ctx = CGBitmapContextCreate( nullptr, w, h, 8, 4 * w,
                                  colorSpace,
                                  kCGImageAlphaPremultipliedLast );
  CGColorSpaceRelease( colorSpace );

  if( !m_ctx )
    throw std::runtime_error( "CGPaintDevice: failed to create bitmap context" );
}//CGPaintDevice( width, height )


CGPaintDevice::CGPaintDevice( CGPaintDevice &parent,
                              double x, double y,
                              double width, double height )
  : m_width( WLength( width, LengthUnit::Pixel ) )
  , m_height( WLength( height, LengthUnit::Pixel ) )
  , m_painter( nullptr )
  , m_ctx( parent.m_ctx )
  , m_ownsContext( false )
  , m_originX( x )
  , m_originY( y )
  , m_bitmapHeight( parent.m_bitmapHeight )
  , m_ctFont( nullptr )
{
}//CGPaintDevice( parent, x, y, width, height )


CGPaintDevice::~CGPaintDevice()
{
  if( m_ctFont )
    CFRelease( m_ctFont );

  if( m_ownsContext && m_ctx )
    CGContextRelease( m_ctx );
}//~CGPaintDevice()


CGImageRef CGPaintDevice::createCGImage() const
{
  assert( m_ownsContext );
  if( !m_ownsContext || !m_ctx )
    return nullptr;
  return CGBitmapContextCreateImage( m_ctx );
}//createCGImage()


CGContextRef CGPaintDevice::cgContext() const
{
  return m_ctx;
}//cgContext()


WFlags<PaintDeviceFeatureFlag> CGPaintDevice::features() const
{
  return PaintDeviceFeatureFlag::FontMetrics;
}


WLength CGPaintDevice::width() const
{
  return m_width;
}


WLength CGPaintDevice::height() const
{
  return m_height;
}


bool CGPaintDevice::paintActive() const
{
  return (m_painter != nullptr);
}


WPainter *CGPaintDevice::painter() const
{
  return m_painter;
}


void CGPaintDevice::setPainter( WPainter *painter )
{
  m_painter = painter;
}


void CGPaintDevice::init()
{
  // Save the base context state
  CGContextSaveGState( m_ctx );

  // CGBitmapContext has origin at bottom-left (Y-up); Wt uses top-left (Y-down).
  // Translate so that local (0,0) maps to the sub-region's top-left in bitmap coords,
  // then flip Y so Y increases downward.
  //
  // In bitmap coords (Y-up), the sub-region's top-left is at:
  //   bitmap_x = m_originX
  //   bitmap_y = m_bitmapHeight - m_originY  (top of image = max Y in bitmap)
  //
  // After Scale(1,-1), local y maps to (bitmap_y - y), which correctly makes
  // y=0 at the sub-region top and y=m_height at the sub-region bottom.
  CGContextTranslateCTM( m_ctx, m_originX, m_bitmapHeight - m_originY );
  CGContextScaleCTM( m_ctx, 1.0, -1.0 );

  // Save a second state for transform/clipping changes (the per-frame state)
  CGContextSaveGState( m_ctx );
}//init()


void CGPaintDevice::done()
{
  CGContextRestoreGState( m_ctx );  // per-frame state
  CGContextRestoreGState( m_ctx );  // base coordinate flip
}//done()


void CGPaintDevice::applyTransform( const WTransform &t )
{
  // WTransform stores [m11, m12, m21, m22, dx, dy] in the same order
  // as CGAffineTransform [a, b, c, d, tx, ty]
  const CGAffineTransform cg = CGAffineTransformMake( t.m11(), t.m12(),
                                                       t.m21(), t.m22(),
                                                       t.dx(), t.dy() );
  CGContextConcatCTM( m_ctx, cg );
}//applyTransform(...)


void CGPaintDevice::setChanged( WFlags<PainterChangeFlag> flags )
{
  if( flags.test(PainterChangeFlag::Transform) || flags.test(PainterChangeFlag::Clipping) )
  {
    // Restore to the clean per-frame state, then re-save for next change
    CGContextRestoreGState( m_ctx );
    CGContextSaveGState( m_ctx );

    // Apply clipping if active
    if( m_painter->hasClipping() && !m_painter->clipPath().isEmpty() )
    {
      const WTransform &clipT = m_painter->clipPathTransform();
      applyTransform( clipT );
      drawPlainPath( m_painter->clipPath() );
      CGContextClip( m_ctx );
      applyTransform( clipT.inverted() );
    }

    // Apply the combined model+view transform
    applyTransform( m_painter->combinedTransform() );

    // After transform change, must re-apply pen/brush/font
    flags = PainterChangeFlag::Pen | PainterChangeFlag::Brush | PainterChangeFlag::Font;
  }

  if( flags.test(PainterChangeFlag::Pen) )
  {
    const WPen &pen = m_painter->pen();

    if( pen.style() != PenStyle::None )
    {
      const WColor &c = pen.color();
      CGContextSetRGBStrokeColor( m_ctx,
                                   c.red() / 255.0, c.green() / 255.0,
                                   c.blue() / 255.0, c.alpha() / 255.0 );

      // CG treats width 0 as invisible; Wt convention (matching WRasterImage-gm) is 1px
      const double penWidth = pen.width().toPixels();
      CGContextSetLineWidth( m_ctx, (penWidth == 0.0) ? 1.0 : penWidth );

      switch( pen.capStyle() )
      {
        case PenCapStyle::Flat:   CGContextSetLineCap( m_ctx, kCGLineCapButt );   break;
        case PenCapStyle::Square: CGContextSetLineCap( m_ctx, kCGLineCapSquare ); break;
        case PenCapStyle::Round:  CGContextSetLineCap( m_ctx, kCGLineCapRound );  break;
      }

      switch( pen.joinStyle() )
      {
        case PenJoinStyle::Miter: CGContextSetLineJoin( m_ctx, kCGLineJoinMiter ); break;
        case PenJoinStyle::Bevel: CGContextSetLineJoin( m_ctx, kCGLineJoinBevel ); break;
        case PenJoinStyle::Round: CGContextSetLineJoin( m_ctx, kCGLineJoinRound ); break;
      }

      switch( pen.style() )
      {
        case PenStyle::None:
          break;
        case PenStyle::SolidLine:
          CGContextSetLineDash( m_ctx, 0, nullptr, 0 );
          break;
        case PenStyle::DashLine:
        {
          const CGFloat pattern[] = { 4.0, 2.0 };
          CGContextSetLineDash( m_ctx, 0, pattern, 2 );
          break;
        }
        case PenStyle::DotLine:
        {
          const CGFloat pattern[] = { 1.0, 2.0 };
          CGContextSetLineDash( m_ctx, 0, pattern, 2 );
          break;
        }
        case PenStyle::DashDotLine:
        {
          const CGFloat pattern[] = { 4.0, 2.0, 1.0, 2.0 };
          CGContextSetLineDash( m_ctx, 0, pattern, 4 );
          break;
        }
        case PenStyle::DashDotDotLine:
        {
          const CGFloat pattern[] = { 4.0, 2.0, 1.0, 2.0, 1.0, 2.0 };
          CGContextSetLineDash( m_ctx, 0, pattern, 6 );
          break;
        }
      }//switch( pen.style() )
    }
  }//if( flags & PainterChangeFlag::Pen )

  if( flags.test(PainterChangeFlag::Brush) )
  {
    const WBrush &brush = m_painter->brush();
    if( brush.style() != BrushStyle::None )
    {
      const WColor &c = brush.color();
      CGContextSetRGBFillColor( m_ctx,
                                 c.red() / 255.0, c.green() / 255.0,
                                 c.blue() / 255.0, c.alpha() / 255.0 );
    }
  }

  if( flags.test(PainterChangeFlag::Font) )
  {
    updateFont( m_painter->font() );
  }
}//setChanged(...)


void CGPaintDevice::drawLine( double x1, double y1, double x2, double y2 )
{
  if( m_painter->pen().style() == PenStyle::None )
    return;

  CGContextBeginPath( m_ctx );
  CGContextMoveToPoint( m_ctx, x1, y1 );
  CGContextAddLineToPoint( m_ctx, x2, y2 );
  CGContextStrokePath( m_ctx );
}//drawLine(...)


void CGPaintDevice::drawPlainPath( const WPainterPath &path )
{
  if( path.isEmpty() )
    return;

  const std::vector<WPainterPath::Segment> &segments = path.segments();

  if( !segments.empty() && segments[0].type() != MoveTo )
    CGContextMoveToPoint( m_ctx, 0, 0 );

  for( size_t i = 0; i < segments.size(); ++i )
  {
    const WPainterPath::Segment &s = segments[i];

    switch( s.type() )
    {
      case MoveTo:
        CGContextMoveToPoint( m_ctx, s.x(), s.y() );
        break;

      case LineTo:
        CGContextAddLineToPoint( m_ctx, s.x(), s.y() );
        break;

      case CubicC1:
      {
        const double x1 = s.x();
        const double y1 = s.y();
        const double x2 = segments[i + 1].x();
        const double y2 = segments[i + 1].y();
        const double x3 = segments[i + 2].x();
        const double y3 = segments[i + 2].y();
        CGContextAddCurveToPoint( m_ctx, x1, y1, x2, y2, x3, y3 );
        i += 2;
        break;
      }

      case CubicC2:
      case CubicEnd:
        assert( false );
        break;

      case ArcC:
      {
        const double cx = s.x();
        const double cy = s.y();
        const double radius = segments[i + 1].x();
        const double startAngleDeg = segments[i + 2].x();
        double spanAngleDeg = segments[i + 2].y();

        static const double EPSILON = 1e-4;
        if( std::fabs( spanAngleDeg ) >= (360.0 - EPSILON) )
          spanAngleDeg = 360.0;

        // Wt angles: degrees, counter-clockwise from 3 o'clock
        // CG angles: radians, counter-clockwise from 3 o'clock in Y-up
        // Since we flipped Y, CG clockwise=1 gives visual counter-clockwise
        const double startRad = startAngleDeg * M_PI / 180.0;
        const double endRad = (startAngleDeg + spanAngleDeg) * M_PI / 180.0;

        // In flipped coords, clockwise=1 renders as CCW visually (matching Wt)
        const int clockwise = (spanAngleDeg > 0) ? 1 : 0;
        CGContextAddArc( m_ctx, cx, cy, radius, startRad, endRad, clockwise );

        i += 2;
        break;
      }

      case ArcR:
      case ArcAngleSweep:
        assert( false );
        break;

      case QuadC:
      {
        const double cpx = s.x();
        const double cpy = s.y();
        const double endx = segments[i + 1].x();
        const double endy = segments[i + 1].y();
        CGContextAddQuadCurveToPoint( m_ctx, cpx, cpy, endx, endy );
        i += 1;
        break;
      }

      case QuadEnd:
        assert( false );
        break;
    }//switch( s.type() )
  }//for( size_t i = 0; ... )
}//drawPlainPath(...)


void CGPaintDevice::strokeAndFillPath()
{
  const bool hasStroke = (m_painter->pen().style() != PenStyle::None);
  const bool hasFill = (m_painter->brush().style() != BrushStyle::None);

  if( hasStroke && hasFill )
    CGContextDrawPath( m_ctx, kCGPathFillStroke );
  else if( hasStroke )
    CGContextStrokePath( m_ctx );
  else if( hasFill )
    CGContextFillPath( m_ctx );
  else
    CGContextBeginPath( m_ctx );  // discard the path
}//strokeAndFillPath()


void CGPaintDevice::drawPath( const WPainterPath &path )
{
  CGContextBeginPath( m_ctx );
  drawPlainPath( path );
  strokeAndFillPath();
}//drawPath(...)


void CGPaintDevice::drawRect( const WRectF &rectangle )
{
  WPainterPath path;
  path.addRect( rectangle );
  drawPath( path );
}//drawRect(...)


void CGPaintDevice::drawArc( const WRectF &rect, double startAngle, double spanAngle )
{
  CGContextSaveGState( m_ctx );

  // Transform to unit circle centered at the ellipse center
  CGContextTranslateCTM( m_ctx, rect.center().x(), rect.center().y() );
  CGContextScaleCTM( m_ctx, rect.width() / 2.0, rect.height() / 2.0 );

  const double startRad = startAngle * M_PI / 180.0;
  const double endRad = (startAngle + spanAngle) * M_PI / 180.0;

  // In flipped coords, clockwise=1 renders as CCW visually (matching Wt positive span)
  const int clockwise = (spanAngle > 0) ? 1 : 0;

  CGContextBeginPath( m_ctx );
  CGContextAddArc( m_ctx, 0, 0, 1.0, startRad, endRad, clockwise );

  // Stroke/fill while the ellipse scaling CTM is still active
  strokeAndFillPath();
  CGContextRestoreGState( m_ctx );
}//drawArc(...)


void CGPaintDevice::drawImage( const WRectF &rect, const std::string &imageUri,
                                int imgWidth, int imgHeight,
                                const WRectF &sourceRect )
{
  // Load image from file path
  CFStringRef cfPath = CFStringCreateWithCString( nullptr, imageUri.c_str(),
                                                   kCFStringEncodingUTF8 );
  if( !cfPath )
    return;

  CFURLRef url = CFURLCreateWithFileSystemPath( nullptr, cfPath,
                                                 kCFURLPOSIXPathStyle, false );
  CFRelease( cfPath );
  if( !url )
    return;

  CGImageSourceRef imgSource = CGImageSourceCreateWithURL( url, nullptr );
  CFRelease( url );
  if( !imgSource )
    return;

  CGImageRef image = CGImageSourceCreateImageAtIndex( imgSource, 0, nullptr );
  CFRelease( imgSource );
  if( !image )
    return;

  // Crop to sourceRect if it doesn't cover the full image
  CGImageRef drawImage = image;
  const bool needsCrop = (sourceRect.x() != 0 || sourceRect.y() != 0
                          || sourceRect.width() != imgWidth
                          || sourceRect.height() != imgHeight);
  if( needsCrop )
  {
    const CGRect cropRect = CGRectMake( sourceRect.x(), sourceRect.y(),
                                         sourceRect.width(), sourceRect.height() );
    drawImage = CGImageCreateWithImageInRect( image, cropRect );
    CGImageRelease( image );
    if( !drawImage )
      return;
  }

  // CGContextDrawImage draws images upside-down in our flipped coordinate system.
  // Temporarily un-flip for the image.
  CGContextSaveGState( m_ctx );
  CGContextTranslateCTM( m_ctx, rect.x(), rect.y() + rect.height() );
  CGContextScaleCTM( m_ctx, 1.0, -1.0 );
  CGContextDrawImage( m_ctx, CGRectMake( 0, 0, rect.width(), rect.height() ), drawImage );
  CGContextRestoreGState( m_ctx );

  CGImageRelease( drawImage );
}//drawImage(...)


void CGPaintDevice::drawImage( const WRectF &rect, const WAbstractDataInfo *imageInfo,
                                int imgWidth, int imgHeight,
                                const WRectF &sourceRect )
{
  // The QuickLook code paths only construct WPainter::Image via the file-path
  // ctor, so the painter dispatches to the string-URI overload above; this
  // overload exists only to satisfy the Wt 4 pure-virtual contract.
  assert( false );
}//drawImage(WAbstractDataInfo*)


CTFontRef CGPaintDevice::createCTFont( const WFont &font )
{
  CFStringRef familyName = nullptr;

  // Use specific family name if set
  const std::string specific = font.specificFamilies().toUTF8();
  if( !specific.empty() )
  {
    // Strip quotes and take the first family name (before any comma)
    std::string name = specific;
    const size_t comma = name.find( ',' );
    if( comma != std::string::npos )
      name = name.substr( 0, comma );

    // Remove leading/trailing quotes and whitespace
    while( !name.empty() && (name.front() == '\'' || name.front() == '"' || name.front() == ' ') )
      name.erase( name.begin() );
    while( !name.empty() && (name.back() == '\'' || name.back() == '"' || name.back() == ' ') )
      name.pop_back();

    if( !name.empty() )
      familyName = CFStringCreateWithCString( nullptr, name.c_str(), kCFStringEncodingUTF8 );
  }

  // Fall back to generic family
  if( !familyName )
  {
    switch( font.genericFamily() )
    {
      case FontFamily::Serif:     familyName = CFSTR( "Times New Roman" ); CFRetain( familyName ); break;
      case FontFamily::Monospace: familyName = CFSTR( "Courier" );         CFRetain( familyName ); break;
      case FontFamily::Cursive:   familyName = CFSTR( "Apple Chancery" );  CFRetain( familyName ); break;
      case FontFamily::Fantasy:   familyName = CFSTR( "Papyrus" );         CFRetain( familyName ); break;
      default:                    familyName = CFSTR( "Helvetica" );       CFRetain( familyName ); break;
    }
  }

  const double size = font.sizeLength( 16 ).toPixels();
  CTFontRef baseFont = CTFontCreateWithName( familyName, size, nullptr );
  CFRelease( familyName );

  if( !baseFont )
    return nullptr;

  // Apply bold/italic traits
  CTFontSymbolicTraits desiredTraits = 0;
  if( font.style() == FontStyle::Italic || font.style() == FontStyle::Oblique )
    desiredTraits |= kCTFontItalicTrait;
  if( font.weight() == FontWeight::Bold || font.weight() == FontWeight::Bolder )
    desiredTraits |= kCTFontBoldTrait;

  if( desiredTraits != 0 )
  {
    CTFontRef styledFont = CTFontCreateCopyWithSymbolicTraits( baseFont, size, nullptr,
                                                                desiredTraits, desiredTraits );
    if( styledFont )
    {
      CFRelease( baseFont );
      return styledFont;
    }
    // If the trait variant doesn't exist, fall through to the base font
  }

  return baseFont;
}//createCTFont(...)


void CGPaintDevice::updateFont( const WFont &font )
{
  if( m_ctFont )
  {
    CFRelease( m_ctFont );
    m_ctFont = nullptr;
  }

  m_ctFont = createCTFont( font );
}//updateFont(...)


void CGPaintDevice::drawText( const WRectF &rect, WFlags<AlignmentFlag> flags,
                               TextFlag textFlag, const WString &text,
                               const WPointF *clipPoint )
{
  // Clip-point check: if a clip point is provided and is outside the clip path, skip
  if( clipPoint && m_painter && !m_painter->clipPath().isEmpty() )
  {
    if( !m_painter->clipPathTransform().map( m_painter->clipPath() )
          .isPointInPath( m_painter->worldTransform().map( *clipPoint ) ) )
      return;
  }

  if( !m_ctFont )
    return;

  const std::string utf8 = text.toUTF8();
  if( utf8.empty() )
    return;

  // Build attributed string with current font and stroke color for text
  const WColor &c = m_painter->pen().color();
  CGFloat colorComponents[4] = {
    c.red() / 255.0, c.green() / 255.0,
    c.blue() / 255.0, c.alpha() / 255.0
  };
  CGColorSpaceRef cs = CGColorSpaceCreateWithName( kCGColorSpaceSRGB );
  CGColorRef cgColor = CGColorCreate( cs, colorComponents );
  CGColorSpaceRelease( cs );

  CFStringRef cfStr = CFStringCreateWithCString( nullptr, utf8.c_str(), kCFStringEncodingUTF8 );
  if( !cfStr )
  {
    CGColorRelease( cgColor );
    return;
  }

  const void *keys[] = { kCTFontAttributeName, kCTForegroundColorAttributeName };
  const void *vals[] = { m_ctFont, cgColor };
  CFDictionaryRef attrs = CFDictionaryCreate( nullptr, keys, vals, 2,
                                               &kCFTypeDictionaryKeyCallBacks,
                                               &kCFTypeDictionaryValueCallBacks );
  CFAttributedStringRef attrStr = CFAttributedStringCreate( nullptr, cfStr, attrs );
  CFRelease( cfStr );
  CFRelease( attrs );
  CGColorRelease( cgColor );

  CTLineRef line = CTLineCreateWithAttributedString( attrStr );
  CFRelease( attrStr );

  if( !line )
    return;

  // Measure the text to compute alignment
  double ascent = 0, descent = 0, leading = 0;
  const double textWidth = CTLineGetTypographicBounds( line, &ascent, &descent, &leading );

  // Horizontal alignment
  const WFlags<AlignmentFlag> hAlign = flags & AlignHorizontalMask;
  double textX;
  if( hAlign.test( AlignmentFlag::Right ) )
    textX = rect.right() - textWidth;
  else if( hAlign.test( AlignmentFlag::Center ) )
    textX = rect.center().x() - textWidth / 2;
  else
    textX = rect.left();

  // Vertical alignment
  const WFlags<AlignmentFlag> vAlign = flags & AlignVerticalMask;
  double textY;
  if( vAlign.test( AlignmentFlag::Middle ) )
    textY = rect.center().y() + (ascent - descent) / 2;
  else if( vAlign.test( AlignmentFlag::Bottom ) )
    textY = rect.bottom() - descent;
  else
    textY = rect.top() + ascent;

  // Core Text draws text with Y-up baseline coordinates.
  // In our flipped coordinate system, we must un-flip locally for text rendering.
  CGContextSaveGState( m_ctx );
  CGContextTranslateCTM( m_ctx, 0, textY );
  CGContextScaleCTM( m_ctx, 1.0, -1.0 );
  CGContextSetTextPosition( m_ctx, textX, 0 );
  CTLineDraw( line, m_ctx );
  CGContextRestoreGState( m_ctx );

  CFRelease( line );
}//drawText(...)


WTextItem CGPaintDevice::measureText( const WString &text, double maxWidth, bool wordWrap )
{
  if( !m_ctFont )
    return WTextItem( text, 0 );

  const std::string utf8 = text.toUTF8();
  if( utf8.empty() )
    return WTextItem( WString::Empty, 0 );

  CFStringRef cfStr = CFStringCreateWithCString( nullptr, utf8.c_str(), kCFStringEncodingUTF8 );
  if( !cfStr )
    return WTextItem( text, 0 );

  const void *keys[] = { kCTFontAttributeName };
  const void *vals[] = { m_ctFont };
  CFDictionaryRef attrs = CFDictionaryCreate( nullptr, keys, vals, 1,
                                               &kCFTypeDictionaryKeyCallBacks,
                                               &kCFTypeDictionaryValueCallBacks );
  CFAttributedStringRef attrStr = CFAttributedStringCreate( nullptr, cfStr, attrs );
  CFRelease( attrs );

  if( maxWidth >= 0 && wordWrap )
  {
    // Use typesetter to find word-break boundary within maxWidth
    CTTypesetterRef typesetter = CTTypesetterCreateWithAttributedString( attrStr );
    const CFIndex breakIndex = CTTypesetterSuggestLineBreak( typesetter, 0, maxWidth );
    CFRelease( typesetter );

    if( breakIndex < CFStringGetLength( cfStr ) )
    {
      // Only measure up to the break point
      CFAttributedStringRef subStr = CFAttributedStringCreateWithSubstring(
        nullptr, attrStr, CFRangeMake( 0, breakIndex ) );
      CTLineRef subLine = CTLineCreateWithAttributedString( subStr );
      const double w = CTLineGetTypographicBounds( subLine, nullptr, nullptr, nullptr );
      CFRelease( subLine );
      CFRelease( subStr );

      // Get next-line width
      CFAttributedStringRef restStr = CFAttributedStringCreateWithSubstring(
        nullptr, attrStr, CFRangeMake( breakIndex, CFStringGetLength( cfStr ) - breakIndex ) );
      CTLineRef restLine = CTLineCreateWithAttributedString( restStr );
      const double nextW = CTLineGetTypographicBounds( restLine, nullptr, nullptr, nullptr );
      CFRelease( restLine );
      CFRelease( restStr );

      // Extract the broken text
      const CFIndex bufSize = breakIndex * 4 + 1;  // generous UTF-8 buffer
      char *buf = new char[bufSize];
      CFStringGetCString( cfStr, buf, bufSize, kCFStringEncodingUTF8 );
      // Truncate to breakIndex characters (approximate -- CFIndex is UTF-16)
      std::string broken( buf, static_cast<size_t>( breakIndex ) );
      delete[] buf;

      CFRelease( cfStr );
      CFRelease( attrStr );
      return WTextItem( WString::fromUTF8( broken ), w, nextW );
    }
  }

  // Measure the full text
  CTLineRef line = CTLineCreateWithAttributedString( attrStr );
  const double w = CTLineGetTypographicBounds( line, nullptr, nullptr, nullptr );
  CFRelease( line );
  CFRelease( cfStr );
  CFRelease( attrStr );

  return WTextItem( text, w );
}//measureText(...)


WFontMetrics CGPaintDevice::fontMetrics()
{
  if( !m_ctFont )
  {
    // Ensure we have a font
    if( m_painter )
      updateFont( m_painter->font() );
    if( !m_ctFont )
      return WFontMetrics( WFont(), 0, 0, 0 );
  }

  const double ascent = CTFontGetAscent( m_ctFont );
  const double descent = CTFontGetDescent( m_ctFont );
  const double leading = CTFontGetLeading( m_ctFont );

  return WFontMetrics( m_painter ? m_painter->font() : WFont(),
                       leading, ascent, descent );
}//fontMetrics()
