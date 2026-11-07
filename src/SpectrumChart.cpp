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

#include "InterSpec_config.h"

#include <map>
#include <set>
#include <cmath>
#include <regex>
#include <tuple>
#include <limits>
#include <string>
#include <cstring>
#include <optional>
#include <algorithm>
#include <functional>

#include <boost/functional/hash.hpp>

#include <Wt/WAny.h>
#include <Wt/WFont.h>
#include <Wt/WRectF.h>
#include <Wt/WLength.h>
#include <Wt/WPainter.h>
#include <Wt/Chart/WAxis.h>
#include <Wt/WApplication.h>
#include <Wt/WPaintDevice.h>
#include <Wt/WPainterPath.h>
#include <Wt/Chart/WDataSeries.h>
#include <Wt/Chart/WCartesianChart.h>
#include <Wt/Chart/WStandardChartProxyModel.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/SpectrumDataModel.h"

#include "src/js_inline/SpectrumChart.js"

using namespace std;
using namespace Wt;

using SpecUtils::Measurement;

namespace
{
  const WColor ns_defaultTimeHighlightColors[3] = { {255,255,0,155}, {0,255,255,75}, {0,128,0,75} };
  const WColor ns_default_peakFillColor( 0, 51, 255 );
  const WColor ns_default_occupiedTimeColor( 128, 128, 128, 75 );

  // Appearance constants of SpectrumChartD3.js (and its CSS)
  const int ns_major_tick_length = 7;
  const int ns_minor_tick_length = 4;
  const double ns_peak_fill_opacity = 0.6;
  const double ns_ref_line_top_pad = 30.0;
  const WColor ns_major_grid_color( 179, 179, 179 );
  const WColor ns_minor_grid_color( 230, 230, 230 );
  const WColor ns_legend_fill_color( 245, 245, 245, 128 );
  const WColor ns_legend_stroke_color( 185, 185, 185, 204 );
  
struct MyTickLabel
{
  enum TickLength { Zero, Short, Long };
    
  double u;
  TickLength tickLength;
  WString    label;
    
  MyTickLabel(double v, TickLength length, const WString& l = WString())
  : u(v), tickLength(length), label(l)
  { }
};//struct MyTickLabel


double estimated_text_width( const size_t nchars, const double font_px );


/** Returns if an axis' labels have been hidden, by setting a blank label format (e.g., " "). */
bool labels_hidden( const Chart::WAxis &axis )
{
  return SpecUtils::trim_copy( axis.labelFormat().toUTF8() ).empty();
}


/** JavaScript's `Number.toPrecision(precision)`; e.g., (1500,3) gives "1.50e+3", and (0.05,3) gives "0.0500". */
std::string js_to_precision( const double v, const int precision )
{
  if( IsNan(v) || IsInf(v) )
    return IsNan(v) ? "NaN" : ((v > 0) ? "Infinity" : "-Infinity");

  char buffer[64] = { '\0' };
  if( v == 0.0 )
  {
    snprintf( buffer, sizeof(buffer), "%.*f", std::max(precision - 1, 0), 0.0 );
    return buffer;
  }

  // The exponent, after rounding to `precision` significant digits
  snprintf( buffer, sizeof(buffer), "%.*e", std::max(precision - 1, 0), v );
  const char * const epos = strchr( buffer, 'e' );
  const int exponent = epos ? atoi( epos + 1 ) : 0;

  if( epos && ((exponent < -6) || (exponent >= precision)) )
    return string( buffer, static_cast<size_t>(epos - buffer) ) + "e" + ((exponent < 0) ? "-" : "+") + std::to_string( std::abs(exponent) );

  snprintf( buffer, sizeof(buffer), "%.*f", std::max( precision - 1 - exponent, 0 ), v );
  return buffer;
}//std::string js_to_precision( const double v, const int precision )


/** Formats a y-axis number as SpectrumChartD3.js's `formatYNumber` does; e.g., 1000 gives "1e+3",
 1500 gives "1.5e+3", 0.5 gives "0.5", and 500 gives "500".
 */
std::string format_y_number( const double v )
{
  if( (v >= 1000.0) || (v < 0.1) )
  {
    static const std::regex zeros_before_e( "\\.0+e" ), trailing_zeros( "\\.0+$" ), zeros_e( "0+e" );

    string t = js_to_precision( v, 3 );
    t = std::regex_replace( t, zeros_before_e, "e" );
    t = std::regex_replace( t, trailing_zeros, "" );
    if( t.find('.') != string::npos )
      t = std::regex_replace( t, zeros_e, "e" );
    return t;
  }//if( v >= 1000.0 || v < 0.1 )

  char buffer[64] = { '\0' };
  const int ndecimal = std::max( 0, 2 - static_cast<int>( std::floor( std::log10(v) ) ) );
  snprintf( buffer, sizeof(buffer), "%.*f", ndecimal, v );
  string t = buffer;
  if( t.find('.') != string::npos )
  {
    while( !t.empty() && (t.back() == '0') )
      t.pop_back();
  }
  if( !t.empty() && (t.back() == '.') )
    t.pop_back();
  return t;
}//format_y_number(...)


/** Formats a rate to about `ndig` significant digits, as SpectrumChartD3.js's legend does. */
std::string legend_rate_str( const double val, const int ndig )
{
  char buffer[64] = { '\0' };
  const double powTen = std::floor( std::log10( std::fabs(val) ) );

  if( (val == std::floor(val)) || IsNan(powTen) || IsInf(powTen) )
    snprintf( buffer, sizeof(buffer), "%.0f", val );
  else if( (powTen < -4) || (powTen > ndig) )
  {
    // As JavaScript's `toExponential(ndig)`; e.g., "6.096e-6", not "6.096e-06"
    static const std::regex exponent_zeros( "e([+-])0+([0-9])" );
    snprintf( buffer, sizeof(buffer), "%.*e", ndig, val );
    return std::regex_replace( string(buffer), exponent_zeros, "e$1$2" );
  }
  else if( powTen < 3 )
    snprintf( buffer, sizeof(buffer), "%.*f", static_cast<int>( ndig - powTen ), val );
  else
    snprintf( buffer, sizeof(buffer), "%.0f", val );

  return buffer;
}//std::string legend_rate_str( const double val, const int ndig )


/** The number of decimals d3's `scale.tickFormat(20, "f")` uses for an axis spanning `span`, which is
 how SpectrumChartD3.js formats x-axis labels.
 */
int d3_tick_precision( const double span )
{
  if( !(span > 0.0) || IsInf(span) )
    return 0;

  const double m = 20.0;
  double step = std::pow( 10.0, std::floor( std::log10( span / m ) ) );
  const double err = m / span * step;
  if( err <= 0.15 )
    step *= 10.0;
  else if( err <= 0.35 )
    step *= 5.0;
  else if( err <= 0.75 )
    step *= 2.0;

  const double p = -std::floor( std::log10( step ) + 0.01 );
  return (IsNan(p) || IsInf(p)) ? 0 : static_cast<int>( std::clamp( p, 0.0, 20.0 ) );
}//int d3_tick_precision( const double span )


/** "Nice" linear ticks from `rmin` to `rmax`, about `nlabel` of them labeled; a port of
 `_linearTicks` in SpectrumChartD3.js.  `format(value, interval)` gives the label of a major tick,
 where `interval` is the spacing of major ticks.
 */
void linear_ticks( const double rmin, const double rmax, const double nlabel, const int max_iter,
                   const std::function<WString(double,double)> &format,
                   std::vector<MyTickLabel> &ticks )
{
  const double EPSILON = 1.0E-3;
  const double range = rmax - rmin;
  const double n = std::pow( 10.0, std::floor( std::log10( range/nlabel ) ) );
  const double msd = range / (n * nlabel);

  if( IsNan(n) || IsInf(n) || !(nlabel > 0.0) || (n <= 0.0) || IsNan(msd) )
    return;

  int subdashes = 0;
  double interval = 0.0;
  if( msd < 1.5 )
  {
    subdashes = 2;
    interval = 0.5*n;
  }else if( msd < 3.3 )
  {
    subdashes = 5;
    interval = 0.5*n;
  }else if( msd < 7 )
  {
    subdashes = 5;
    interval = n;
  }else
  {
    subdashes = 10;
    interval = n;
  }

  const double biginterval = subdashes * interval;
  double start = biginterval * std::floor( (rmin + 1.0E-15) / biginterval );
  if( start < (rmin - EPSILON*interval) )
    start += biginterval;

  for( int i = -static_cast<int>( std::floor( std::floor(start - rmin) / interval ) ); ; ++i )
  {
    if( i > max_iter )  //JIC
      break;

    const double v = start + interval * i;
    if( (v - rmax) > EPSILON * interval )
      break;

    const bool major = ((i % subdashes) == 0);
    const WString label = (major && (i >= 0)) ? format( v, biginterval ) : WString();
    ticks.push_back( MyTickLabel( v, major ? MyTickLabel::Long : MyTickLabel::Short, label ) );
  }//for( intervals to draw ticks for )
}//void linear_ticks(...)


void getXAxisLabelTicks( const SpectrumChart *chart,
                         const Chart::WAxis& axis,
                         std::vector<MyTickLabel>& ticks )
{
  // The ticks fall on "nice" values (multiples of 5, 10, 25, 50, 100, etc.), with labels formatted
  //  like SpectrumChartD3.js does - but always with enough decimals to tell labels apart.
  if( axis.scale() != Chart::AxisScale::Linear )
    throw runtime_error( "getXAxisLabelTicks called for wrong type of graphs");

  const double rendermin = axis.minimum();
  const double rendermax = axis.maximum();

  const double minxpx = chart->mapToDevice(rendermin,0).x();
  const double maxxpx = chart->mapToDevice(rendermax,0).x();
  const double widthpx = maxxpx - minxpx;

  const bool hide_labels = labels_hidden( axis );
  const int d3_precision = d3_tick_precision( rendermax - rendermin );

  const auto format = [=]( const double value, const double interval ) -> WString {
    if( hide_labels )
      return WString();

    int ndecimal = 0;
    while( (ndecimal < 6)
           && (std::fabs( interval*std::pow(10.0,ndecimal) - std::round(interval*std::pow(10.0,ndecimal)) ) > 1.0E-6) )
      ++ndecimal;

    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.*f", std::max( ndecimal, d3_precision ), value );
    if( std::string(buffer).find_first_not_of( "-0." ) == string::npos )
      return WString::fromUTF8( "0" );  //Avoid "-0"
    return WString::fromUTF8( buffer );
  };//format

  // Labels are at least `xAxisMinimumLabelSpacing()` apart, and further if they are wide enough
  //  to otherwise run into each other (as SpectrumChartD3.js allows more room for more digits).
  const double font_px = axis.labelFont().sizeLength( 12 ).toPixels();
  double nlabel = std::floor( widthpx / chart->xAxisMinimumLabelSpacing() );
  for( int iteration = 0; iteration < 5; ++iteration )
  {
    ticks.clear();
    linear_ticks( rendermin, rendermax, nlabel, 5000, format, ticks );

    size_t max_chars = 0;
    double label_spacing_px = widthpx;
    std::optional<double> prev_label_px;
    for( const MyTickLabel &tick : ticks )
    {
      if( tick.label.empty() )
        continue;
      max_chars = std::max( max_chars, tick.label.value().size() );
      const double px = chart->mapToDevice( tick.u, 0 ).x();
      if( prev_label_px )
        label_spacing_px = std::min( label_spacing_px, px - *prev_label_px );
      prev_label_px = px;
    }//for( const MyTickLabel &tick : ticks )

    const double needed_px = estimated_text_width( max_chars, font_px ) + font_px;
    if( (max_chars == 0) || (label_spacing_px >= needed_px) || (nlabel <= 1.0) )
      break;

    nlabel = std::max( 1.0, std::floor( nlabel * label_spacing_px / needed_px ) );
  }//for( iterations to get labels spaced far enough apart )
}//getXAxisLabelTicks(...)
  
  
void getYAxisLabelTicks( const SpectrumChart *chart,
                         const Chart::WAxis& axis,
                         std::vector<MyTickLabel> &ticks )
{
  // A port of SpectrumChartD3.js's `yticks()`: ticks fall on "nice" values, and on a log axis, at
  //  each decade, with minor ticks at 2 through 9 times the decade, when there is room.
  const double EPSILON = 1.0E-3;

  const double renderymin = axis.minimum();
  const double renderymax = axis.maximum();
  const double minypx = chart->mapToDevice(0.0,renderymax,axis.id()).y();
  const double maxypx = chart->mapToDevice(0.0,renderymin,axis.id()).y();
  const double heightpx = maxypx - minypx;

  const bool hide_labels = labels_hidden( axis );
  const auto label = [hide_labels]( const double v ) -> WString {
    return hide_labels ? WString() : WString::fromUTF8( format_y_number( v ) );
  };

  switch( axis.scale() )
  {
    case Chart::AxisScale::Discrete: case Chart::AxisScale::Date: case Chart::AxisScale::DateTime:
      throw runtime_error( "getYAxisLabelTicks() for SpectrumChart can only handle linear and log scales");

    case Chart::AxisScale::Linear:
    {
      linear_ticks( renderymin, renderymax, heightpx / 50.0, 500,
                    [&label]( const double v, const double ){ return label( v ); }, ticks );
      break;
    }//case Chart::AxisScale::Linear:

    case Chart::AxisScale::Log:
    {
      const int yMinPrefPow = -1;  // i.e., floor( log10(logYAxisMin=0.1) )

      //Get the power of 10 just below or equal to rendermin.
      int minpower = (renderymin > 0.0) ? static_cast<int>( std::floor( log10(renderymin) ) ) : yMinPrefPow;

      //Get the power of 10 just above or equal to renderymax.  If renderymax
      //  is less than or equal to 0, set power to be 0.
      int maxpower = (renderymax > 0.0) ? static_cast<int>( std::ceil( log10(renderymax) ) ) : 0;

      //Adjust minpower and maxpower
      if( maxpower == minpower )
      {
        //Happens when renderymin==renderymax which is a power of 10
        ++maxpower;
        --minpower;
      }else if( (maxpower > 2) && (minpower < yMinPrefPow) )
      {
        //We had a tiny value (possibly a fraction of a count), as well as a
        //  large value (>1000).
        minpower = yMinPrefPow;
      }else if( (maxpower >= 0) && (minpower < yMinPrefPow) && ((maxpower - minpower) > 6) )
      {
        //we had a tiny power (1.0E-5), as well as one between 1 and 999,
        //  so we will only show the most significant decades
        minpower = maxpower - 5;
      }//if( minpower == maxpower ) / else / else

      //numdecades: number of decades the data covers, including the decade
      //  above and below the data.
      const int numdecades = maxpower - minpower + 1;

      //minpxdecade: minimum number of pixels we need per decade.
      const double minpxdecade = 25.0;

      //labeldelta: the number of decades between successive labeled large ticks;
      //  each decade will have a large tick regardless
      int labeldelta = 1;

      //nticksperdecade: number of small+large ticks per decade
      int nticksperdecade = 10;

      if( std::floor( heightpx / minpxdecade ) < numdecades )
      {
        labeldelta = static_cast<int>( std::round( numdecades / (heightpx / minpxdecade) ) );
        labeldelta = std::max( labeldelta, 1 );
        nticksperdecade = 1;
      }

      int nmajorticks = 0;

      for( int decade = minpower; decade <= maxpower; ++decade )
      {
        const double startcounts = std::pow( 10.0, decade );
        const double deltacounts = 10.0 * startcounts / nticksperdecade;
        const double eps = deltacounts * EPSILON;

        if( ((startcounts - renderymin) > -eps) && ((startcounts - renderymax) < eps) )
        {
          ++nmajorticks;
          const WString t = ((decade % labeldelta) == 0) ? label( startcounts ) : WString();
          ticks.push_back( MyTickLabel( startcounts, MyTickLabel::Long, t ) );
        }

        for( int i = 1; i < (nticksperdecade - 1); ++i )
        {
          const double y = startcounts + i*deltacounts;
          if( ((y - renderymin) > -eps) && ((y - renderymax) < eps) )
            ticks.push_back( MyTickLabel( y, MyTickLabel::Short ) );
        }
      }//for( int decade = minpower; decade <= maxpower; ++decade )

      // Only decades are labeled; if none are in range, we'll label the first and last ticks.
      if( (nmajorticks < 1) && !ticks.empty() )
      {
        ticks.front().label = label( ticks.front().u );
        ticks.back().label = label( ticks.back().u );
      }

      if( ticks.empty() )
      {
        const double y = 0.5*(renderymin + renderymax);
        ticks.push_back( MyTickLabel( y, MyTickLabel::Long, label( y ) ) );
      }

      break;
    }//case Chart::AxisScale::Log:
  }//switch( axis.scale() )
}//getYAxisLabelTicks(...)


/** Appends the pixel points `i0` through `i1` (in reverse, if `rev`) to `path`, leaving out points
 within 0.1 px (vertically) of the line drawn past them; a port of `coords` in SpectrumChartD3.js.
 The first point is moved to, if `move_to`, or else a line is drawn to it.
 */
void append_thinned( WPainterPath &path, const std::vector<double> &xpx,
                     const std::function<double(size_t)> &ypx,
                     const size_t i0, const size_t i1, const bool rev, const bool move_to )
{
  double ax = 0.0, ay = 0.0, px = 0.0, py = 0.0;
  double lo = 0.0, hi = 0.0;  // The range of slopes, from the last point drawn, of a line within 0.1 px of all points since
  bool first = true;

  const auto draw = [&]( const double x, const double y ){
    if( first && move_to )
      path.moveTo( x, y );
    else
      path.lineTo( x, y );
    first = false;
    ax = x;
    ay = y;
    lo = -std::numeric_limits<double>::infinity();
    hi = std::numeric_limits<double>::infinity();
  };//draw

  for( size_t j = i0; j <= i1; ++j )
  {
    const size_t i = rev ? (i0 + i1 - j) : j;
    const double x = xpx[i], y = ypx(i);

    if( j == i0 )
    {
      draw( x, y );
    }else
    {
      double u = std::fabs( x - ax );
      const double slope = (y - ay) / u;
      if( !(slope >= lo) || !(slope <= hi) )
      {
        draw( px, py );
        u = std::fabs( x - ax );
      }

      if( u > 0.0 )
      {
        lo = std::max( lo, (y - ay - 0.1) / u );
        hi = std::min( hi, (y - ay + 0.1) / u );
      }
    }//if( j == i0 ) / else

    px = x;
    py = y;
  }//for( size_t j = i0; j <= i1; ++j )

  if( i1 > i0 )
    draw( px, py );
}//void append_thinned(...)


/** Adds a vertical line, from `y0` to `y1` at `x`, to `path` - dashed by the given pattern of
 alternating dash and gap lengths (as SVG `stroke-dasharray`), or solid if it is empty.

 Dashes are drawn as separate segments, so any pattern works with every paint device.
 */
void add_dashed_vline( WPainterPath &path, const double x, const double y0, const double y1,
                       const std::vector<double> &pattern )
{
  const double len = std::fabs( y1 - y0 );
  const double dir = (y1 >= y0) ? 1.0 : -1.0;

  double total = 0.0;
  for( const double d : pattern )
    total += d;

  if( pattern.empty() || !(total > 0.0) )
  {
    path.moveTo( x, y0 );
    path.lineTo( x, y1 );
    return;
  }

  // An odd number of entries is repeated to make it even, as SVG does.
  std::vector<double> dashes = pattern;
  if( dashes.size() % 2 )
    dashes.insert( end(dashes), begin(pattern), end(pattern) );

  double pos = 0.0;
  for( size_t i = 0; pos < len; i = (i + 1) % dashes.size() )
  {
    const double next = std::min( len, pos + dashes[i] );
    if( (i % 2) == 0 )
    {
      path.moveTo( x, y0 + dir*pos );
      path.lineTo( x, y0 + dir*next );
    }
    pos = next;
  }//for( loop over dashes )
}//void add_dashed_vline(...)


/** An estimate of the width of `nchars` characters of a font `font_px` pixels tall, for when a
 paint device cant measure text (e.g., WSvgImage); generous, as digits are about 0.56 em wide.
 */
double estimated_text_width( const size_t nchars, const double font_px )
{
  return 0.6 * font_px * nchars;
}


/** The width of `text` in the painters current font: measured if the paint device can, or else
 estimated.
 */
double text_width( WPainter &painter, const WString &text )
{
  WPaintDevice * const device = painter.device();
  if( device && device->features().test( PaintDeviceFeatureFlag::FontMetrics ) )
    return device->measureText( text ).width();
  return estimated_text_width( text.value().size(), painter.font().sizeLength( 12 ).toPixels() );
}


/** Returns if two colors have the same RGBA values. */
bool same_color( const WColor &lhs, const WColor &rhs )
{
  return (lhs.red() == rhs.red()) && (lhs.green() == rhs.green())
         && (lhs.blue() == rhs.blue()) && (lhs.alpha() == rhs.alpha());
}


/** Returns `color` with its alpha multiplied by `frac`. */
WColor with_alpha( const WColor &color, const double frac )
{
  const int alpha = static_cast<int>( std::round( frac * color.alpha() ) );
  return WColor( color.red(), color.green(), color.blue(), std::clamp( alpha, 0, 255 ) );
}


/** Returns the SpectrumDataModel powering `chart`, or nullptr. */
SpectrumDataModel *spectrum_model( const Chart::WCartesianChart &chart )
{
  auto proxyModel = dynamic_cast<Wt::Chart::WStandardChartProxyModel *>( chart.model().get() );
  return dynamic_cast<SpectrumDataModel *>( proxyModel ? proxyModel->sourceModel().get() : nullptr );
}
}//namespace



class SpectrumRenderIterator;

class SeriesRenderer
{
public:
  virtual ~SeriesRenderer() { }
  virtual void addValue(double x, double y, double stacky,
                        const WModelIndex &xIndex, const WModelIndex &yIndex) = 0;
  virtual void paint() = 0;
  
protected:
  const SpectrumChart &chart_;
  WPainter            &painter_;
  const Chart::WDataSeries &series_;
  
  SeriesRenderer( const SpectrumChart &chart,
                 WPainter &painter,
                 const Chart::WDataSeries& series,
                 SpectrumRenderIterator& it)
  : chart_(chart),
  painter_(painter),
  series_(series),
  it_(it)
  { }
  
  static double crisp(double u) {
    return std::floor(u) + 0.5;
  }
  
  WPointF hv(const WPointF& p) {
    return chart_.hv(p);
  }
  
  WPointF hv(double x, double y) {
    return chart_.hv(x, y);
  }
  
protected:
  SpectrumRenderIterator& it_;
}; //class SeriesRenderer


class SpectrumRenderIterator : public Chart::SeriesIterator
{
public:
  SpectrumRenderIterator( const SpectrumChart &chart, Wt::WPainter &painter);
  virtual void startSegment(int currentXSegment, int currentYSegment,
                            const WRectF& currentSegmentArea);
  virtual void endSegment();
  virtual bool startSeries( const Chart::WDataSeries& series, double groupWidth,
                           int , int);
  virtual void endSeries();
  virtual void newValue(const Chart::WDataSeries& series, double x, double y,
                        double stackY,
                        const WModelIndex& xIndex,
                        const WModelIndex& yIndex);
  double breakY(double y);
  
private:
  const SpectrumChart &chart_;
  Wt::WPainter &painter_;
  const Chart::WDataSeries *series_;
  
  SpectrumSeriesRenderer    *seriesRenderer_;
  double             minY_, maxY_;
};//class SpectrumRenderIterator



class SpectrumSeriesRenderer : public SeriesRenderer
{
public:
  SpectrumSeriesRenderer( const SpectrumChart &chart,
                         WPainter &painter,
                         const Chart::WDataSeries &series,
                         SpectrumRenderIterator& it)
  : SeriesRenderer( chart, painter, series, it ),
  curveLength_(0)
  { }
  
  void addValue( double x, double y, double stacky,
                const WModelIndex &xIndex, const WModelIndex &yIndex)
  {
    WPointF p = chart_.map( x, y, series_.axis(),
                           it_.currentXSegment(), it_.currentYSegment());
    
    if( curveLength_ == 0 )
    {
      curve_.moveTo(hv(p));
      
      if (series_.fillRange() != Chart::FillRangeType::None
          && series_.brush().style() != BrushStyle::None) {
        fill_.moveTo(hv(fillOtherPoint(x)));
        fill_.lineTo(hv(p));
      }
    }else
    {
      if (series_.type() == Chart::SeriesType::Line) {
        curve_.lineTo(hv(p));
        fill_.lineTo(hv(p));
      } else {
        if (curveLength_ == 1) {
          computeC(p0, p, c_);
        } else {
          WPointF c1, c2;
          computeC(p_1, p0, p, c1, c2);
          curve_.cubicTo(hv(c_), hv(c1), hv(p0));
          fill_.cubicTo(hv(c_), hv(c1), hv(p0));
          c_ = c2;
        }
      }
    }
    
    p_1 = p0;
    p0 = p;
    lastX_ = x;
    ++curveLength_;
  }//addValue(...)
  
  void paint()
  {
    //On my 2.6 GHz Intel Core i7 macbook pro, during rendering many times, I
    //  got cpu usage (20150316):
    //  WPainter::strokePath:             10ms
    //  SpectrumDataModel::index:         3ms
    //  SpectrumDataModel::data:          2ms
    //  SpectrumSeriesRenderer::addValue: 2ms
    //  SpectrumDataModel::rowWidth:      1ms
    
    if( curveLength_ > 1 )
    {
      if( series_.type() == Wt::Chart::SeriesType::Curve )
      {
        WPointF c1;
        computeC(p0, p_1, c1);
        curve_.cubicTo(hv(c_), hv(c1), hv(p0));
        fill_.cubicTo(hv(c_), hv(c1), hv(p0));
      }
      
      if (series_.fillRange() != Chart::FillRangeType::None
          && series_.brush().style() != BrushStyle::None) {
        fill_.lineTo(hv(fillOtherPoint(lastX_)));
        fill_.closeSubPath();
        painter_.setShadow(series_.shadow());
        painter_.fillPath(fill_, series_.brush());
      }
      
      if (series_.fillRange() == Chart::FillRangeType::None)
        painter_.setShadow(series_.shadow());
      else
        painter_.setShadow(WShadow());
      
      painter_.strokePath(curve_, series_.pen());
    }//if( curveLength_ > 1 )
    
    curveLength_ = 0;
    curve_ = WPainterPath();
    fill_ = WPainterPath();
  }//void paint()
  
private:
  int curveLength_;
  WPainterPath curve_;
  WPainterPath fill_;
  
  double  lastX_;
  WPointF p_1, p0, c_;
  
  static double dist(const WPointF& p1, const WPointF& p2)
  {
    double dx = p2.x() - p1.x();
    double dy = p2.y() - p1.y();
    return std::sqrt( dx*dx + dy*dy);
  }
  
  static void computeC( const WPointF &p, const WPointF &p1, WPointF &c )
  {
    c.setX(p.x() + 0.3 * (p1.x() - p.x()));
    c.setY(p.y() + 0.3 * (p1.y() - p.y()));
  }
  
  static void computeC(const WPointF& p_1, const WPointF& p0,
                       const WPointF& p1,
                       WPointF& c1, WPointF& c2)
  {
    double m1x = (p_1.x() + p0.x())/2.0;
    double m1y = (p_1.y() + p0.y())/2.0;
    
    double m2x = (p0.x() + p1.x())/2.0;
    double m2y = (p0.y() + p1.y())/2.0;
    
    double L1 = dist(p_1, p0);
    double L2 = dist(p0, p1);
    double r = L1/(L1 + L2);
    
    c1.setX(p0.x() - r * (m2x - m1x));
    c1.setY(p0.y() - r * (m2y - m1y));
    
    r = 1-r;
    
    c2.setX(p0.x() - r * (m1x - m2x));
    c2.setY(p0.y() - r * (m1y - m2y));
  }//void computeC(...)
  
  WPointF fillOtherPoint(double x) const
  {
    Chart::FillRangeType fr = series_.fillRange();
    
    switch( fr )
    {
      case Chart::FillRangeType::MinimumValue:
        return WPointF(chart_.map(x, 0, series_.axis(),
                                  it_.currentXSegment(),
                                  it_.currentYSegment()).x(),
                       chart_.chartArea().bottom());
      case Chart::FillRangeType::MaximumValue:
        return WPointF(chart_.map(x, 0, series_.axis(),
                                  it_.currentXSegment(),
                                  it_.currentYSegment()).x(),
                       chart_.chartArea().top());
      case Chart::FillRangeType::ZeroValue:
        return WPointF(chart_.map(x, 0, series_.axis(),
                                  it_.currentXSegment(),
                                  it_.currentYSegment()));
      default:
        return WPointF();
    }//switch( fr )
  }//WPointF fillOtherPoint(double x) const
};//class SpectrumSeriesRenderer : public SeriesRenderer



SpectrumRenderIterator::SpectrumRenderIterator( const SpectrumChart &chart,
                                                Wt::WPainter &painter )
  : chart_( chart ),
  painter_( painter ),
  series_(0)
  { }
  
void SpectrumRenderIterator::startSegment(int currentXSegment, int currentYSegment,
                            const WRectF& currentSegmentArea)
  {
    Chart::SeriesIterator::startSegment(currentXSegment, currentYSegment,
                                        currentSegmentArea);
    
    const Chart::WAxis& yAxis = chart_.axis(series_->axis());
    
    if (currentYSegment == 0)
      maxY_ = DBL_MAX;
    else
      maxY_ = currentSegmentArea.bottom();
    
    if (currentYSegment == yAxis.segmentCount() - 1)
      minY_ = -DBL_MAX;
    else
      minY_ = currentSegmentArea.top();
  }
  
  
  
void SpectrumRenderIterator::endSegment()
{
  Chart::SeriesIterator::endSegment();
  seriesRenderer_->paint();
}
  
  
bool SpectrumRenderIterator::startSeries( const Chart::WDataSeries& series, double groupWidth,
                           int , int )
{
  seriesRenderer_ = new SpectrumSeriesRenderer( chart_, painter_, series, *this);
  series_ = &series;
  painter_.save();

  return !series.isHidden();
}
  
void SpectrumRenderIterator::endSeries()
{
  seriesRenderer_->paint();
  painter_.restore();
    
  delete seriesRenderer_;
  series_ = 0;
}
  
  
void SpectrumRenderIterator::newValue(const Chart::WDataSeries& series, double x, double y,
                        double stackY,
                        const WModelIndex& xIndex,
                        const WModelIndex& yIndex)
{
  if( IsNan(x) || IsNan(y) )
    seriesRenderer_->paint();
  else
    seriesRenderer_->addValue( x, y, stackY, xIndex, yIndex);
}


double SpectrumRenderIterator::breakY(double y)
{
  if (y < minY_)
    return minY_;
  else if (y > maxY_)
    return maxY_;
  else
    return y;
}



void SpectrumChart::calcRenderRectangle( const Wt::WRectF &rectangle ) const
{
  int w, h;
  if( rectangle.isNull() || rectangle.isEmpty() )
  {
    w = static_cast<int>( width().toPixels());
    h = static_cast<int>( height().toPixels());
  }else
  {
    w = static_cast<int>( rectangle.width());
    h = static_cast<int>( rectangle.height());
  }
  
  const int padLeft = plotAreaPadding(Side::Left);
  const int padRight = plotAreaPadding(Side::Right);
  const int padTop = plotAreaPadding(Side::Top);
  const int padBottom = plotAreaPadding(Side::Bottom);
  
  if( orientation() == Orientation::Vertical )
    m_renderRectangle = WRectF( padLeft, padTop,
                  std::max(10, w - padLeft - padRight),
                  std::max(10, h - padTop - padBottom));
  else
    m_renderRectangle = WRectF( padTop, padRight,
                  std::max(10, w - padTop - padBottom),
                  std::max(10, h - padRight - padLeft));
  
}//void calcRenderRectangle( const Wt::WRectF &rectangle )


void SpectrumChart::render(	Wt::WPainter &painter, const Wt::WRectF &rectangle ) const
{
  calcRenderRectangle( rectangle);
  
  // Grid lines, if enabled, are drawn by `paint(...)` over the spectrum and peaks, as in
  //  SpectrumChartD3.js
  if( initLayout( rectangle, painter.device() ) )
  {
    renderChartBackground(painter, rectangle);
    renderSeries(painter);
    renderAxes(painter, Chart::AxisProperty::Line | Chart::AxisProperty::Labels);
    renderLegend(painter);
  }//if( initLayout(rectangle) )
}//void render(...)


const Wt::WRectF &SpectrumChart::chartArea() const
{
  return m_renderRectangle;
}//WRectF chartArea() const


void SpectrumChart::renderChartBackground( Wt::WPainter &painter, const Wt::WRectF &rectangle ) const
{
  const WFont oldFont = painter.font();
  
  auto paintMargins = [&](){
    const double rx = m_renderRectangle.x();
    const double ry = m_renderRectangle.y();
    const double rw = m_renderRectangle.width();
    const double rh = m_renderRectangle.height();
    if( ry > 0 )  //Top strip
      painter.fillRect( hv(WRectF(rx,0,rw,ry)), m_chartMarginBrush);
    if( (ry+rh) < rectangle.height() ) //Bottom strip
      painter.fillRect( hv(WRectF(rx,ry+rh,rw,rectangle.height()-ry)), m_chartMarginBrush);
    if( rx > 0 )
      painter.fillRect( hv(WRectF(0,0,rx,rectangle.height())), m_chartMarginBrush);
    if( (rx+rw) < rectangle.width() )
      painter.fillRect( hv(WRectF(rx+rw,0,rectangle.width()-rw-rx,rectangle.height())), m_chartMarginBrush);
  };
  
  if( background().style() != BrushStyle::None && m_chartMarginBrush.style() != BrushStyle::None )
  {
    paintMargins();
    painter.fillRect( hv(chartArea()), background());
  }else if( background().style() != BrushStyle::None )
  {
    painter.fillRect( hv(rectangle), background());
  }else if( m_chartMarginBrush.style() != BrushStyle::None )
  {
    paintMargins();
  }
  
  
  SpectrumDataModel * const specModel = spectrum_model( *this );
  assert( specModel );

  
  if( textInMiddleOfChart().empty()
     && specModel && !specModel->histUsedForXAxis() )
  {
    const double renderxmin = axis(Chart::Axis::X).minimum();
    const double renderxmax = axis(Chart::Axis::X).maximum();
    const double renderymin = axis(Chart::Axis::Y1).minimum();
    const double renderymax = axis(Chart::Axis::Y1).maximum();
    const WPointF upperLeft = mapToDevice(renderxmin,renderymax);
    const WPointF lowerRight = mapToDevice(renderxmax,renderymin);
    
    const double minxpx = upperLeft.x();
    const double maxxpx = lowerRight.x();
    const double minypx = upperLeft.y();
    const double maxypx = lowerRight.y();
    const double xrange = maxxpx - minxpx;
    const double yrange = maxypx - minypx;
    
    //ToDo: Get SVG version of logo
    WPainter::Image snllogo( "InterSpec_resources/images/SNL_Stacked_Black_Blue.png", 688, 265);  //17 kb
    
    //Lets draw the Sandia logo to take up about 1/4th the height of the screen
    const double logo_scale = 0.25 * yrange / snllogo.height();
    const double logo_width = logo_scale * snllogo.width();
    const double logo_height = logo_scale * snllogo.height();
    
    const double middle_px_x = 0.5*(minxpx + maxxpx);
    const double middle_px_y = 0.5*(minypx + maxypx);
    
    //Try to Make "InterSpec" be about the same width as the SNL logo
    const int font_size = static_cast<int>( logo_width / 4);
    
    WRectF snlrect( middle_px_x - 0.5*logo_width, middle_px_y - 0.5*logo_height - font_size, logo_width, logo_height);
    painter.drawImage( snlrect, snllogo);
    
    painter.save();
    WFont font( FontFamily::Monospace);
    font.setSize( font_size);
    font.setWeight( FontWeight::Normal);
    font.setVariant( FontVariant::SmallCaps);
    painter.setFont( font);
    painter.setPen( m_textPen);
    double texty = middle_px_y + 0.5*logo_height - 1.05*font_size;
    painter.drawText( minxpx, texty, xrange, 1.05*font_size, AlignmentFlag::Center, "InterSpec");
    painter.restore();
  }//if( specModel && !specModel->histUsedForXAxis() )
  
  painter.setFont( oldFont);
}//void renderChartBackground( Wt::WPainter &painter ) const;


void SpectrumChart::renderGridLines( Wt::WPainter &painter,
                                const Wt::Chart::Axis axistype ) const
{
  const Chart::WAxis &axis = this->axis(axistype);
  
  if( !axis.isGridLinesEnabled() )
    return;
  
  bool vertical = true;
  vector<MyTickLabel> ticks;
  
  switch( axistype )
  {
    case Wt::Chart::Axis::X:
      vertical = true;
      getXAxisLabelTicks( this, axis, ticks);
    break;
    
    case Wt::Chart::Axis::Y1: case Wt::Chart::Axis::Y2:
//    case Wt::Chart::Axis::Y1: case Wt::Chart::Axis::Ordinate:
      vertical = false;
      getYAxisLabelTicks( this, axis, ticks);
    break;
      
    default:
      return;
      break;
  }//switch( axistype )
  
  // As SpectrumChartD3.js: a line at every tick, lighter for minor ticks.  If the axis' grid pen
  //  was set (i.e., isnt Wt's default), it is used; otherwise the D3 charts colors are (which it
  //  uses for every color theme).
  const bool default_pen = (axis.gridLinesPen() == WPen( StandardColor::Gray ));
  const WPen major_pen = default_pen ? WPen( ns_major_grid_color ) : axis.gridLinesPen();
  WPen minor_pen = major_pen;
  minor_pen.setColor( default_pen ? ns_minor_grid_color : with_alpha( major_pen.color(), 0.45 ) );
  const Chart::Axis yaxis_id = vertical ? Chart::Axis::Y1 : axistype;
  const WRectF &area = chartArea();

  WPainterPath majorGridPath, minorGridPath;
  majorGridPath.setOpenSubPathsEnabled( true );
  minorGridPath.setOpenSubPathsEnabled( true );

  for( const MyTickLabel &tick : ticks )
  {
    if( tick.tickLength == MyTickLabel::Zero )
      continue;

    WPainterPath &path = (tick.tickLength == MyTickLabel::Long) ? majorGridPath : minorGridPath;
    if( vertical )
    {
      const double x = std::floor( mapToDevice( tick.u, 0.0 ).x() ) + 0.5;
      path.moveTo( x, area.top() );
      path.lineTo( x, area.bottom() );
    }else
    {
      const double y = std::floor( mapToDevice( 0.0, tick.u, yaxis_id ).y() ) + 0.5;
      path.moveTo( area.left(), y );
      path.lineTo( area.right(), y );
    }
  }//for( const MyTickLabel &tick : ticks )

  if( !minorGridPath.isEmpty() )
    painter.strokePath( minorGridPath, minor_pen );

  if( !majorGridPath.isEmpty() )
    painter.strokePath( majorGridPath, major_pen );
}//void renderGridLines( Wt::WPainter &painter, const Wt::WAxis &axis ) const



void SpectrumChart::renderSeries( Wt::WPainter &painter ) const
{
  if( !isLargeEnough( painter ) )
  {
    const double height = static_cast<double>( painter.viewPort().height());
    const double width = static_cast<double>( painter.viewPort().width());
    painter.drawText( 0.0, 0.5*height - 8.0, width, 20.0,
                       AlignmentFlag::Center, "Enlarge to view chart");
    return;
  }//if( !isLargeEnough() )
  
  {
    SpectrumRenderIterator iterator( *this, painter);
    iterateSpectrum( &iterator, painter);
  }
  
  //Label and marker iterator code removed 20150224 since they were causing a
  // slowdown (the marker iterator was accessing all rows of the data model
  // no matter the x-range), and we arent currently using them anyway.
//  {
//    LabelRenderIterator iterator(*this, painter);
//    iterateSeries( &iterator, &painter);
//  }
//  {
//    MarkerRenderIterator iterator(*this, painter);
//    iterateSeries( &iterator, &painter);
//  }
}//void renderSeries( Wt::WPainter &painter ) const;


void SpectrumChart::renderAxes( Wt::WPainter &painter,
                      Wt::WFlags<  Wt::Chart::AxisProperty > properties ) const
{
  if( !isLargeEnough(painter) )
    return;
  
  SpectrumDataModel * const m = spectrum_model( *this );
  assert( m );
  
  if( m )
  {
    renderXAxis( painter, axis(Chart::Axis::X), properties);
    renderYAxis( painter, axis(Chart::Axis::Y1), properties);
    renderYAxis( painter, axis(Chart::Axis::Y2), properties);
  }else
  {
    renderAxis( painter, axis(Chart::Axis::X), properties);
    renderAxis( painter, axis(Chart::Axis::Y1), properties);
    renderAxis( painter, axis(Chart::Axis::Y2), properties);
  }//if( m ) / else
}//void renderAxes(...) const;


void SpectrumChart::renderYAxis( Wt::WPainter &painter,
                                 const Wt::Chart::WAxis &yaxis,
                                Wt::WFlags<Wt::Chart::AxisProperty> properties ) const
{
  //Adapted from Wt 3.3.4, WCartesianChart::renderAxis(...).
  //  wcjohns specialized so that y-min/y-max can be arbitrary, but the axis
  //  markers will still be multiples of 10.
  if( !yaxis.isVisible() )
    return;

  WPointF axisStart, axisEnd;

  axisStart.setY(chartArea().bottom() + 0.5);
  axisEnd.setY(chartArea().top() + 0.5);

  enum { Left = 0x1, Right = 0x2, Both = 0x3 } tickPos = Left;

  Chart::AxisValue location = yaxis.location();
  if( yaxis.id() == Chart::Axis::Y2 )
    location = Chart::AxisValue::Maximum;

  switch( location )
  {
    case Chart::AxisValue::Minimum:
    {
      const double x = chartArea().left() - yaxis.margin() + 0.5;
      axisStart.setX( x);
      axisEnd.setX( x);
      tickPos = Left;

      break;
    }//case Chart::AxisValue::Minimum:

    case Chart::AxisValue::Maximum:
    {
      double x = chartArea().right() + yaxis.margin() + 0.5;
      axisStart.setX(x);
      axisEnd.setX(x);
      tickPos = Right;

      break;
    }//case Chart::AxisValue::Maximum:

    case Chart::AxisValue::Zero:
    {
      double x = std::floor(map(0, 0, Chart::Axis::Y1).x()) + 0.5;
      axisStart.setX(x);
      axisEnd.setX(x);
      tickPos = Both;

      break;
    }//case Chart::AxisValue::Zero:


    case Chart::AxisValue::Both:
    {
      throw runtime_error( "SpectrumChart::renderYAxis(...): Chart::AxisValue::Both not implemented.");
      break;
    }//case Chart::AxisValue::Both:
  }//switch( location )


  //Render Title, rotated, along the left edge of the image (or the right edge, for Y2)
  if( properties.test(Chart::AxisProperty::Labels) && !yaxis.title().empty() )
  {
    const double h = axisStart.y() - axisEnd.y();
    const double fontheight = yaxis.titleFont().sizeLength().toPixels();

    painter.save();
    painter.setPen( m_textPen);
    painter.setBrush( WBrush());
    painter.setFont( yaxis.titleFont());

    if( yaxis.id() == Chart::Axis::Y2 )
    {
      // Reading downwards, against the right edge
      const double right = painter.window().width();
      painter.translate( right, axisEnd.y() );
      painter.rotate( 90 );
      painter.drawText( 0.0, 1.0, h, fontheight + 6, AlignmentFlag::Center | AlignmentFlag::Top, yaxis.title() );
    }else
    {
      painter.rotate( -90);
      painter.drawText( -axisStart.y(), 1.0, h, fontheight + 6,
                        AlignmentFlag::Center | AlignmentFlag::Top, yaxis.title() );
    }
    painter.restore();
  }//if( draw labels and there is a title )
  
  
  //If paint the chart line
  if( properties.test(Chart::AxisProperty::Line) )
  {
    const WPen oldPen = painter.pen();
    painter.setPen( yaxis.pen());
    painter.drawLine( hv(axisStart), hv(axisEnd));
    painter.setPen( oldPen);
  }//if( properties.test(Chart::AxisProperty::Line) )
  
  
  //Draw the labels
  std::vector<MyTickLabel> ticks;
  getYAxisLabelTicks( this, yaxis, ticks);
  const size_t nticks = std::min( int(ticks.size()), 1000);
  
  WPainterPath ticksPath;
  ticksPath.setOpenSubPathsEnabled( true );
  
  painter.save();
  const WFont oldFont = painter.font();
  const WPen oldPen = painter.pen();
  WFont labelFont = yaxis.labelFont();
  painter.setFont( labelFont);
  
  for( size_t i = 0; i < nticks; ++i )
  {
    const double yval = ticks[i].u;
    const double ypx = std::floor(mapToDevice(1.0,yval,yaxis.id()).y()) + 0.5;
    
    const int tickLength = (ticks[i].tickLength == MyTickLabel::Long)
                            ? ns_major_tick_length : ns_minor_tick_length;

    if( ticks[i].tickLength != MyTickLabel::Zero )
    {
      double x = axisStart.x() + (tickPos & Right ? +tickLength : 0);
      const double y = hv(x,ypx).y();
      ticksPath.moveTo( x, y);
      x = axisStart.x() + (tickPos & Left ? -tickLength : 0);
      ticksPath.lineTo( x, y);
    }//if( non zero tick length )
    
    const int txtwidth = 200;
    const int txtheight = 20;
    const int txtPadding = 3;
    
    double labelxpx, labelypx;
    WFlags<AlignmentFlag> labelAlignFlags;
    
    switch( location )
    {
      case Chart::AxisValue::Minimum:
      case Chart::AxisValue::Zero:
        labelxpx = axisStart.x() + (tickPos & Left ? -tickLength : 0) - txtwidth - txtPadding;
        labelypx = ypx - txtheight/2;
        labelAlignFlags = AlignmentFlag::Right | AlignmentFlag::Middle;
      break;
        
      case Chart::AxisValue::Maximum:
        labelxpx = axisStart.x() + (tickPos & Right ? +tickLength : 0);
        labelypx = ypx - txtheight/2;
        labelAlignFlags = AlignmentFlag::Left | AlignmentFlag::Middle;
      break;
        
      case Chart::AxisValue::Both:
      {
        throw runtime_error( "SpectrumChart::renderYAxis(...): Chart::AxisValue::Both not implemented.");
        break;
      }//case Chart::AxisValue::Both:
    }//switch( location )
    
    if( properties.test(Chart::AxisProperty::Labels) && !ticks[i].label.empty() )
    {   
      painter.setPen( m_textPen);
      painter.setFont( labelFont);
      painter.drawText( WRectF(labelxpx, labelypx, txtwidth, txtheight),
                        labelAlignFlags, ticks[i].label);
      painter.setPen( oldPen);
    }
  }//for( size_t i = 0; i < nticks; ++i )
  
  painter.restore();
  painter.setFont( oldFont);
  painter.setPen( oldPen);
  
  if ((properties.test(Chart::AxisProperty::Line)) && !ticksPath.isEmpty() )
    painter.strokePath( ticksPath, yaxis.pen());
}//void renderYAxis(...)



void SpectrumChart::renderXAxis( Wt::WPainter &painter,
                                const Wt::Chart::WAxis& axis,
                                Wt::WFlags<Wt::Chart::AxisProperty> properties ) const
{
  if( !axis.isVisible() )
    return;
  
  //This function adapted from WChart2DRenderer::renderAxis(...) in order to
  //  make x axis labels always line up nice
  const WPen oldPen = painter.pen();
  WFont oldFont1 = painter.font();
  WFont labelFont = axis.labelFont();
  painter.setFont(labelFont);
  
  double lineypx = 0;
  enum { Left = 0x1, Right = 0x2, Both = 0x3 } tickPos = Left;
  AlignmentFlag labelHFlag = AlignmentFlag::Left;
  
  Chart::AxisValue location = axis.location();
  
  switch( location )
  {
    case Chart::AxisValue::Minimum:
      tickPos = Left;
      labelHFlag = AlignmentFlag::Top;
      lineypx = chartArea().bottom() + 0.5 + axis.margin();
      break;
      
    case Chart::AxisValue::Maximum:
      tickPos = Right;
      labelHFlag = AlignmentFlag::Bottom;
      lineypx = chartArea().top() - 0.5 - axis.margin();
      break;
      
    case Chart::AxisValue::Zero:
      tickPos = Both;
      labelHFlag = AlignmentFlag::Top;
      lineypx = std::floor(map(0, 0, Chart::Axis::Y1).y()) + 0.5;
      break;
      
    case Chart::AxisValue::Both:
      throw runtime_error( "SpectrumChart::renderXAxis(...): Chart::AxisValue::Both not implemented.");
      break;
  }//switch( location_[axis.id()] )
  
  
  //  WPainterPath tildeStartMarker = WPainterPath();
  //  tildeStartMarker.moveTo(0, 0);
  //  tildeStartMarker.lineTo(0, axis.segmentMargin() - 25);
  //  tildeStartMarker.moveTo(-15, axis.segmentMargin() - 10);
  //  tildeStartMarker.lineTo(15, axis.segmentMargin() - 20);
  
  //  WPainterPath tildeEndMarker = WPainterPath();
  //  tildeEndMarker.moveTo(0, 0);
  //  tildeEndMarker.lineTo(0, -(axis.segmentMargin() - 25));
  //  tildeEndMarker.moveTo(-15, -(axis.segmentMargin() - 20));
  //  tildeEndMarker.lineTo(15, -(axis.segmentMargin() - 10));
  
  //  static const int CLIP_MARGIN = 5;
  
  double minxpx = mapToDevice(axis.minimum(),0).x();
  double maxxpx = mapToDevice(axis.maximum(),0).x();
  
  {
    if( properties.test(Chart::AxisProperty::Line) )
    {
      painter.setPen( axis.pen());
      const WPointF begin = hv( minxpx, lineypx);
      const WPointF end = hv( maxxpx, lineypx);
      painter.drawLine(begin, end);
    }
    
    const double titlewidthpx = 10.0*axis.title().narrow().size();

    WPainterPath ticksPath;
    ticksPath.setOpenSubPathsEnabled( true );

    std::vector<MyTickLabel> ticks;
    getXAxisLabelTicks( this, axis, ticks);

    const size_t nticks = std::min( int(ticks.size()), 1000);

    // The last label is left off if it would run past the right of the chart, as in SpectrumChartD3.js
    size_t last_label_index = nticks;
    for( size_t i = 0; i < nticks; ++i )
      last_label_index = ticks[i].label.empty() ? last_label_index : i;

    for( size_t i = 0; i < nticks; ++i )
    {
      const double d = ticks[i].u;
      const double dd = std::floor(mapToDevice(d,0).x()) + 0.5;

      const int tickLength = (ticks[i].tickLength == MyTickLabel::Long)
                              ? ns_major_tick_length : ns_minor_tick_length;

      const Chart::AxisValue location = axis.location();
      
      if( ticks[i].tickLength != MyTickLabel::Zero )
      {
        ticksPath.moveTo(hv(dd, lineypx + (tickPos & Right ? -tickLength : 0)));
        ticksPath.lineTo(hv(dd, lineypx + (tickPos & Left ? +tickLength : 0)));
      }
      
      WPointF labelPos;
      switch( location )
      {
        case Chart::AxisValue::Minimum:
        case Chart::AxisValue::Zero:
          labelPos = WPointF(dd, lineypx + tickLength);
          break;
          
        case Chart::AxisValue::Maximum:
          labelPos = WPointF(dd, lineypx - tickLength);
          break;
          
        case Chart::AxisValue::Both:
          throw runtime_error( "SpectrumChart::renderXAxis(...): Chart::AxisValue::Both not implemented.");
          break;
      }//switch( location )
      
      if( properties.test(Chart::AxisProperty::Labels) && !ticks[i].label.empty() )
      {
        WFlags<AlignmentFlag> labelFlags = labelHFlag;
        
        if (axis.labelAngle() == 0)
          labelFlags |= AlignmentFlag::Center;
        else if (axis.labelAngle() > 0)
          labelFlags |= AlignmentFlag::Right;
        else
          labelFlags |= AlignmentFlag::Left;
        
        bool drawText = true;
        
        
        if( m_compactAxis )
        {
          const double titlestartpx = (chartArea().right() - titlewidthpx);
          const double labelwidthpx = 10.0*ticks[i].label.toUTF8().size();
          const double labelendpx = labelPos.x() + 0.5*labelwidthpx;
          drawText = (labelendpx < titlestartpx);
        }else if( i == last_label_index )
        {
          const double labelendpx = labelPos.x() + 0.5*text_width( painter, ticks[i].label );
          drawText = ((labelendpx - 10.0) <= chartArea().right());
        }//if( m_compactAxis ) / else
        
        if( drawText )
        {
          const int margin = m_compactAxis ? 0 : 3;
          renderLabel( painter, ticks[i].label, labelPos, labelFlags, axis.labelAngle(), margin );
        }//if( drawText )
      }//if( draw label text )
    }//for( size_t i = 0; i < nticks; ++i )
    
    if( properties.test(Chart::AxisProperty::Line) )
    {
      painter.strokePath(ticksPath, axis.pen());

      // An arrow at the right end of the axis, when there is more data to the right, as
      //  SpectrumChartD3.js does
      const SpectrumDataModel * const m = spectrum_model( *this );
      const shared_ptr<const Measurement> xhist = m ? m->histUsedForXAxis() : nullptr;
      const double data_max = (xhist && xhist->num_gamma_channels()) ? xhist->gamma_energy_max() : 0.0;
      if( (location == Chart::AxisValue::Minimum) && (plotAreaPadding(Side::Right) >= 9)
          && (data_max > (axis.maximum() + 1.0E-6*fabs(data_max))) )
      {
        WPainterPath arrow;
        arrow.moveTo( maxxpx + 2.0, lineypx - 3.0 );
        arrow.lineTo( maxxpx + 2.0, lineypx + 8.0 );
        arrow.lineTo( maxxpx + 8.0, lineypx + 2.0 );
        arrow.closeSubPath();
        painter.save();
        painter.setPen( axis.pen() );
        painter.setBrush( WBrush( axis.pen().color() ) );
        painter.drawPath( arrow );
        painter.restore();
      }//if( there is data past the right of the chart )
    }//if( properties.test(Chart::AxisProperty::Line) )

    if( properties.test(Chart::AxisProperty::Labels) && !axis.title().empty() )
    {
      WFont oldFont2 = painter.font();
      WFont titleFont = axis.titleFont();
      painter.setFont(titleFont);

      //Normal mode is vertical
      const bool vertical = orientation() == Orientation::Vertical;

      if( m_compactAxis && vertical )
      {
        const WPen oldpen = painter.pen();
        const WBrush oldbrush = painter.brush();
        painter.setBrush( m_chartMarginBrush);
        painter.setPen( WPen(Wt::StandardColor::Transparent));
        const double fontheight = titleFont.sizeLength().toPixels();
        const WString &title = axis.title();
        const WPointF pos(chartArea().right() - titlewidthpx, lineypx + ns_major_tick_length);

        const double x = pos.x(), y = lineypx + ns_major_tick_length, h = fontheight + 1;
        painter.drawRect( x, y, titlewidthpx, h);
        painter.setPen( oldpen);
        painter.setBrush( oldbrush);
        
        renderLabel( painter, title, pos, AlignmentFlag::Top | AlignmentFlag::Left, 0, 0 );
      }else
      {
        // Just below the tick labels
        const double label_height = axis.labelFont().sizeLength( 12 ).toPixels();
        const double title_y = lineypx + ns_major_tick_length + 3 + label_height + 2;
        if( vertical )
          renderLabel( painter, axis.title(), WPointF(chartArea().center().x(), title_y),
                       AlignmentFlag::Top | AlignmentFlag::Center, 0, 0 );
        else
          renderLabel( painter, axis.title(), WPointF(chartArea().right(), lineypx),
                       AlignmentFlag::Top | AlignmentFlag::Left, 0, 8 );
      }
      
      painter.setFont(oldFont2);
    }
  }
  
  painter.setFont(oldFont1);
  painter.setPen( oldPen );
}//void renderXAxis(...)


void SpectrumChart::renderLegend( Wt::WPainter &painter ) const
{
  return;
//#error SpectrumChart needs to implement renderLegend for Wt>=3.3.2
}//void renderLegend( Wt::WPainter &painter ) const


bool SpectrumChart::isLargeEnough( Wt::WPainter &painter ) const
{
  const double width  = painter.viewPort().width();
  const double height = painter.viewPort().height();
  const int wpadding  = plotAreaPadding(Side::Left)
                        + plotAreaPadding(Side::Right);
  const int hpadding  = plotAreaPadding(Side::Top)
                        + plotAreaPadding(Side::Bottom);
 
  if( (height-hpadding) < 50 )
    return false;
  
  if( (width-wpadding) < 50 )
    return false;

  return true;
}//bool isLargeEnough() const


void SpectrumChart::iterateSpectrum( SpectrumRenderIterator *iterator,
                                     WPainter &painter ) const
{
  using namespace Wt::Chart;
  const auto &seriesv = series();
  //  unsigned rows = m_model ? m_model->rowCount() : 0;
  
  
  SpectrumDataModel * const m = spectrum_model( *this );
  assert( m);


  if( !m )
    return;

  const int nrow = m->rowCount();

  for( unsigned g = 0; g < seriesv.size(); ++g )
  {
    if( seriesv[g]->isHidden() )
      continue;

    const int endSeries = g, startSeries = g;
    int i = startSeries;

    for( ; ; )
    {
      WDataSeries &series = *seriesv[i];

      const bool doSeries = iterator->startSeries( series, 0.0, 1, 0);

      if( doSeries )
      {
        for( int currentXSegment = 0;
            currentXSegment < axis(Chart::Axis::X).segmentCount();
            ++currentXSegment)
        {
          for( int currentYSegment = 0;
              currentYSegment < axis(series.axis()).segmentCount();
              ++currentYSegment)
          {
            WRectF csa = chartSegmentArea( axis(series.axis()),
                                           currentXSegment, currentYSegment);
            
            const double minXPx = csa.left();
            const double maxXPx = csa.right();
            const double minx = mapFromDevice( WPointF(minXPx,0.0) ).x();
            const double maxx = mapFromDevice( WPointF(maxXPx,0.0) ).x();
            
            // The rows in view (even if just part of one row is)
            const int minRow = std::max( m->findRow( minx ), 0 );
            const int maxRow = std::min( m->findRow( maxx ), nrow - 1 );
            if( (nrow < 1) || (minRow > maxRow) )
              continue;

            iterator->startSegment( currentXSegment, currentYSegment, csa);
            
            painter.save();
            WPainterPath clipPath;
            clipPath.addRect(hv(csa));
            painter.setClipPath(clipPath);
            painter.setClipping(true);
            
            for( int row = minRow; row <= maxRow; ++row )
            {
              int c = series.XSeriesColumn();
              if( c == -1 )
                c = XSeriesColumn();
              const int column = series.modelColumn();
              const WModelIndex xIndex = m->index( row, c);
              const WModelIndex yIndex = m->index( row, column);
              
              if( !yIndex.isValid() || !xIndex.isValid() )
                continue;
              
              const double y = asNumber( m->data(yIndex));

              // Drawn as a histogram (stair steps), as SpectrumChartD3.js does; the rebin
              //  factor should keep channels at least a pixel wide (see displayRebinFactor(...)).
              const double x_start = m->rowLowEdge( xIndex.row());
              const double x_end = x_start + m->rowWidth( xIndex.row());
              iterator->newValue( series, x_start, y, 0, xIndex, yIndex);
              iterator->newValue( series, x_end, y, 0, xIndex, yIndex);
            }//for( unsigned row = 0; row < rows; ++row )
            
            iterator->endSegment();
            painter.restore();
          }//for( int currentYSegment ...
        }//for( int currentXSegment = 0;...)
        
        iterator->endSeries();
      }//if( doSeries )
      
      if( i == endSeries )
        break;
      else
      {
        if( endSeries < startSeries )
          --i;
        else
          ++i;
      }//if( i == endSeries ) / lese
    }//for (;;)
  }//for( unsigned g = 0; g < series.size(); ++g )
}

////////////////////////////////////////////////////////////////////////////////











SpectrumChart::SpectrumChart()
: Wt::Chart::WCartesianChart(),
  m_widthInPixels( 0 ),
  m_heightInPixels ( 0 ),
  m_peakModel( nullptr ),
  m_defaultPeakColor( ns_default_peakFillColor ),
  m_occupiedMarkerColor( ns_default_occupiedTimeColor ),
  m_timeHighlightColors{ ns_defaultTimeHighlightColors[0], ns_defaultTimeHighlightColors[1], ns_defaultTimeHighlightColors[2] },
  m_legendType( SpectrumChart::NoLegend ),
  m_xAxisUnits( SpectrumChart::kUndefinedUnits ),
  m_compactAxis( false ),
  m_xAxisMinLabelSpacing( 50.0 ),
  m_textInMiddleOfChart(),
  m_verticalLinesShowing( false ),
  m_horizontalLinesShowing( false )
{
  addStyleClass( "SpectrumChart");
  
  //Using ScatterPlot instead of CategoryChart since the value of the templates
  //  are reported as the sum of the templates of lower series index than,
  //  and including itself (Chart::ScatterPlot does not support stackings)
  setType( Chart::ChartType::Scatter);
  //m_chart->axis(Chart::Axis::X).setLabelFormat( "%.2f");
  axis(Chart::Axis::X).setScale( Chart::AxisScale::Linear);
  
  for( int label = 0; label < static_cast<int>(PeakLabelType::NumPeakLabels); ++label )
  {
    m_peakLabelsToShow[label] = false;
    m_peakLabelsColors[label] = 0;  //i.e., use the text color
  }//for( loop over all labels )
  
//  setLayoutSizeAware( false);
  
  axis(Chart::Axis::X).setMinimum( 0.0);
  axis(Chart::Axis::X).setMaximum( 3000.0);
}//SpectrumChart( constructor )


SpectrumChart::~SpectrumChart()
{
}//~SpectrumChart()

void SpectrumChart::setPeakLabelColor( PeakLabelType label,
                                       uint32_t red, uint32_t green,
                                       uint32_t blue, uint32_t alpha )
{
  if( red > 255 )   red = 255;
  if( green > 255 ) green = 255;
  if( blue > 255 )  blue = 255;
  if( alpha > 255 ) alpha = 255;
  
  green = (green << 8);
  blue  = (blue << 16);
  alpha = (alpha << 24);
  
  m_peakLabelsColors[static_cast<int>(label)] = (red | green | blue | alpha);
}//void setPeakLabelColor(...)


void SpectrumChart::setShowPeakLabel( PeakLabelType label, bool show )
{
  m_peakLabelsToShow[static_cast<int>(label)] = show;
}//void setShowPeakLabel(...)


bool SpectrumChart::showingPeakLabel( PeakLabelType label ) const
{
  return m_peakLabelsToShow[static_cast<int>(label)];
}

void SpectrumChart::peakModelDataChanged( const WModelIndex &topLeft,
                                          const WModelIndex &bottomRight )
{
  const int left = topLeft.column();
  const int right = bottomRight.column();
  
  if( topLeft.row() != bottomRight.row()
      || left <= PeakModel::kAmplitude
      || left >= PeakModel::kUserLabel
      || right >= PeakModel::kUserLabel )
  {
    update( Wt::WFlags<PaintFlag>{});
  }else if( left <= PeakModel::kPhotoPeakEnergy )
  {
  }
  else if( (left == PeakModel::kUserLabel && showingPeakLabel( PeakLabelType::ShowPeakUserLabel ))
          || (left == PeakModel::kIsotope && showingPeakLabel( PeakLabelType::ShowPeakNuclideLabel ))
          || (left == PeakModel::kPhotoPeakEnergy && showingPeakLabel( PeakLabelType::ShowPeakNuclideEnergies ))  )
  {
    //This next statment is poorly formed, and only catches if a single data thingy
    //  is changed - I'll fix this at some point when I'm not so tired
    update( Wt::WFlags<PaintFlag>{});
  }
}//void peakModelDataChanged(...)


void SpectrumChart::setPeakModel( std::shared_ptr<PeakModel> model )
{
  if( m_peakModel )
    throw runtime_error( "SpectrumChart::setPeakModel(...): a model has "
                         "already been set");

  if( !model )
    throw runtime_error( "SpectrumChart::setPeakModel(...): invalid input model");

  // Co-own the model for as long as this chart references it (no dangling `PeakModel*`).
  m_peakModel = std::move( model );

  //Note that if you first remove rows from m_peakModel, then add some, then
  //  modify some data, all within one server callback, this will all only
  //  cause one update to be pushed to the user
  m_peakModel->dataChanged().connect( this, [this]( const Wt::WModelIndex &a1, const Wt::WModelIndex &a2 ){
    peakModelDataChanged( a1, a2);
  });
  m_peakModel->rowsRemoved().connect( this, [this]( const Wt::WModelIndex &, int, int ){ update(); } );
  m_peakModel->rowsInserted().connect( this, [this]( const Wt::WModelIndex &, int, int ){ update(); } );
}//void setPeakModel( std::shared_ptr<PeakModel> model )


void SpectrumChart::setPeaks( std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks )
{
  m_peaks = std::move( peaks );
  update();
}//void setPeaks(...)


void SpectrumChart::setXAxisUnits( SpectrumChart::XAxisUnits units )
{
  m_xAxisUnits = units;
}//setXAxisUnits( SpectrumChart::XAxisUnits units )


SpectrumChart::XAxisUnits SpectrumChart::xAxisUnits() const
{
  return m_xAxisUnits;
}//XAxisUnits xAxisUnits() const

void SpectrumChart::showGridLines( const bool draw )
{
  Chart::WCartesianChart::axis(Chart::Axis::X).setGridLinesEnabled( draw);
  Chart::WCartesianChart::axis(Chart::Axis::Y1).setGridLinesEnabled( draw);
  update(); //force re-render
}//void showGridLines( const bool draw )

void SpectrumChart::showVerticalLines( const bool draw )
{
  Chart::WCartesianChart::axis(Chart::Axis::X).setGridLinesEnabled( draw);
  update(); //force re-render
  m_verticalLinesShowing = draw;
}

void SpectrumChart::showHorizontalLines( const bool draw )
{
  Chart::WCartesianChart::axis(Chart::Axis::Y1).setGridLinesEnabled( draw);
  update(); //force re-render
  m_horizontalLinesShowing = draw;
}

bool SpectrumChart::verticalLinesShowing() const
{
  return m_verticalLinesShowing;
}

bool SpectrumChart::horizontalLinesShowing() const
{
  return m_horizontalLinesShowing;
}


Wt::WLength SpectrumChart::paintedHeight() const
{
  return WLength( m_heightInPixels, WLength::Unit::Pixel);
}//Wt::WLength paintedHeight() const



Wt::WLength SpectrumChart::paintedWidth() const
{
  return WLength( m_widthInPixels, WLength::Unit::Pixel);
}//Wt::WLength paintedWidth() const


void SpectrumChart::setTextInMiddleOfChart( const Wt::WString &s )
{
  m_textInMiddleOfChart = s;
}


const Wt::WString &SpectrumChart::textInMiddleOfChart() const
{
  return m_textInMiddleOfChart;
}

void SpectrumChart::setCompactAxis( const bool compact )
{
  if( m_compactAxis == compact )
    return;
  
  const int old_pad = plotAreaPadding(Side::Bottom);
  const int new_pad = std::max(old_pad + (compact ? -1 : 1)*17, 0);
  
  m_compactAxis = compact;
  setPlotAreaPadding( new_pad, Side::Bottom);
}

void SpectrumChart::setXAxisMinimumLabelSpacing( const double pixels )
{
  m_xAxisMinLabelSpacing = std::max( pixels, 1.0 );
}


double SpectrumChart::xAxisMinimumLabelSpacing() const
{
  return m_xAxisMinLabelSpacing;
}


bool SpectrumChart::isAxisCompacted() const
{
  return m_compactAxis;
}



void SpectrumChart::setTimeHighLightRegions( const vector< pair<double,double> > &p,
                                         const HighlightRegionType type )
{
  WColor color;
  switch( type )
  {
    case ForegroundHighlight:       color = m_timeHighlightColors[0]; break;
    case BackgroundHighlight:       color = m_timeHighlightColors[1];  break;
    case SecondForegroundHighlight: color = m_timeHighlightColors[2];    break;
  }//switch( type )
  
  vector< HighlightRegion > tokeep;
  for( size_t i = 0; i < m_highlights.size(); ++i )
    if( m_highlights[i].hash != size_t(type) )
      tokeep.push_back( m_highlights[i]);
  
  for( size_t i = 0; i < p.size(); ++i )
  {
    HighlightRegion region;
    region.lowerx = static_cast<float>( p[i].first);
    region.upperx = static_cast<float>( p[i].second);
    region.color = color;
    region.hash = size_t(type);
    tokeep.push_back( region);
  }
  
  tokeep.swap( m_highlights);
}//void setTimeHighLightRegions( vector< pair<double,double> > )


void SpectrumChart::clearOccupancyRegions()
{
  setOccupancyRegions( vector< pair<double,double> >());
}//void clearOccupancyRegions( const HighlightRegionType region )


void SpectrumChart::setOccupancyRegions( const std::vector< std::pair<double,double> > &p )
{
  m_occupancy_regions.clear();
  
  for( const auto &pe : p )
  {
    HighlightRegion region;
    region.lowerx = static_cast<float>( pe.first);
    region.upperx = static_cast<float>( pe.second);
    region.color = m_occupiedMarkerColor;
    region.hash = 0;//size_t(type);
    m_occupancy_regions.push_back( region);
  }
}//void setOccupancyRegions(...)


bool SpectrumChart::removeDecorativeHighlightRegion( size_t hash )
{
  if( hash < 3 )
    return false;
  
  const size_t nprev = m_highlights.size();
  for( size_t i = 0; i < nprev; ++i )
  {
    if( m_highlights[i].hash == hash )
    {
      m_highlights.erase( m_highlights.begin() + i);
      return true;
    }
  }
  
  return false;
}//bool removeDecorativeHighlightRegion( size_t hash )


size_t SpectrumChart::addDecorativeHighlightRegion( const float lowerx,
                                                  const float upperx,
                                                  const Wt::WColor &color )
{
  HighlightRegion region;
  region.lowerx = lowerx;
  region.upperx = upperx;
  region.color = color;
  region.hash = 0;
  boost::hash_combine( region.hash, lowerx);
  boost::hash_combine( region.hash, upperx);
  boost::hash_combine( region.hash, color.red());
  boost::hash_combine( region.hash, color.green());
  boost::hash_combine( region.hash, color.blue());
  boost::hash_combine( region.hash, color.alpha());
  
  if( region.hash <= 2 )
    region.hash += 3;
  
  //should in principle check for collision, but whatever
  
  m_highlights.push_back( region);
  
  return region.hash;
}//void addDecorativeHighlightRegion(...)



void SpectrumChart::setXAxisRange( const double x1, const double x2 )
{
  // (labels are formatted according to the range; see getXAxisLabelTicks(...))
  axis(Chart::Axis::X).setRange( x1, x2);
}//void setXAxisRange()


void SpectrumChart::setYAxisRange( const double y1, const double y2 )
{
  axis(Chart::Axis::Ordinate).setRange( y1, y2);
}//void setYAxisRange()


//TODO: note that setting autolimits on a zoomed in section doesnt
//      really do what we want (eg rescale y-axis according to visible
//      range), therefore you should call setAutoYAxisRange()
void SpectrumChart::setAutoXAxisRange()
{
//  axis(Chart::Axis::X).setAutoLimits( Chart::AxisValue::Minimum | Chart::AxisValue::Maximum);
  //XXX
  //  WAxis doesnt immediately recalculate X and Y axis ranges when
  //  setAutoXAxisRange() and setAutoYAxisRange() are called, therefore
  //  we will set our own limits as an approximation, so
  //  that any function that calls axis(Chart::Axis::X).minimum()/maximum()
  //  before the range is computer by Wt, will at least approximately get
  //  the correct range...  you probably want to call setAutoYAxisRange()
  //  as well.
  //  Doing this means it's neccasary to re-call setAutoXAxisRange()
  //  whenever you load new data
  SpectrumDataModel * const theModel = spectrum_model( *this );
  assert( theModel );
  
  if( theModel )
  {
    std::shared_ptr<const Measurement> axisH = theModel->histUsedForXAxis();
    if( axisH )
    {
      float xmin = min( 0.0f, axisH->gamma_energy_min());
      float xmax = axisH->gamma_energy_max();
      setXAxisRange( xmin, xmax);
    }//if( axisH )
  }else
  {
    axis(Chart::Axis::X).setAutoLimits( Chart::AxisValue::Minimum | Chart::AxisValue::Maximum);
  }//if( theModel is a SpectrumDataModel ) / else
}//void setAutoXAxisRange()


void SpectrumChart::yAxisRangeFromXRange( const double x0, const double x1,
                                          double &y0, double &y1 )
{
  y0 = y1 = 0.0;
  SpectrumDataModel * const theModel = spectrum_model( *this );
  assert( theModel );
  
  if( !theModel )
    throw runtime_error( "SpectrumChart::setAutoXAxisRange(): warning, you "
                         "should only use the SpectrumChart class with a "
                         "SpectrumChart!");
  
  double min_non_zero = 0.0;
  theModel->yRangeInXRange( x0, x1, y0, y1, &min_non_zero );
    
  if( axis(Chart::Axis::Ordinate).scale() == Chart::AxisScale::Log )
  {
    //Specify the (approx) fraction of the chart that the scale should extend
    //  past where the data where hit.
    const double yfractop = 0.05, yfracbottom = 0.025;

    // As SpectrumChartD3.js: if the data never reaches 1.0, base the minimum on the smallest
    //  non-zero value, rather than always going down to 0.1.
    if( (y0 <= 0.0) && (min_non_zero > 0.0) && (min_non_zero < 0.1) && (y1 < 1.0) )
      y0 = 0.5*min_non_zero;
    
    const double y0Intitial = ((y0<=0.0) ? 0.1 : y0);
    double y1Intitial = ((y1<=0.0) ? 1.0 : y1);
    y1Intitial = ((y1Intitial<=y0Intitial) ? 1.1*y0Intitial : y1Intitial);
    
    
    const double logY0 = log10(y0Intitial);
    const double logY1 = log10(y1Intitial);
    
    const double logLowerY = ((y0<=0.0) ? -1.0 : (logY0 - yfracbottom*(logY1-logY0)));
    const double logUpperY = logY1 + yfractop*(logY1-logY0);
    
    const double ylower = std::pow( 10.0, logLowerY);
    const double yupper = std::pow( 10.0, logUpperY);
    
    y0 = ((y0<=0.0) ? 0.1 : ylower);
    y1 = ((y1<=0.0) ? 1.0 : yupper);
  }else //if( axis(Chart::Axis::Y1).scale() == Chart::AxisScale::Linear )
  {
    y0 = ((y0 <= 0.0) ? 1.1*y0 : 0.9*y0);
    y1 = 1.1*y1;
  }//if( Y is LogScale ) / else
}//void yAxisRangeFromXRange(...)


void SpectrumChart::setAutoYAxisRange()
{
  double y0=0.0, y1=0.0;
  const double x0 = axis(Chart::Axis::X).minimum();
  const double x1 = axis(Chart::Axis::X).maximum();

  yAxisRangeFromXRange( x0, x1, y0, y1);

  if( y0 == y1 )
    axis(Chart::Axis::Ordinate).setAutoLimits( Chart::AxisValue::Minimum | Chart::AxisValue::Maximum);
  else
    setYAxisRange( y0, y1);
}//void setAutoYAxisRange()


void SpectrumChart::paintTextInMiddleOfChart( WPainter &painter ) const
{
  if( m_textInMiddleOfChart.empty() )
    return;
  
  painter.save();
  const WFont oldFont = painter.font();
  WFont txtFont = painter.font();
  //    txtFont.setSize( WFont::XLarge);
  txtFont.setSize( WLength(32));
  painter.setFont( txtFont);
  const int left_padding = plotAreaPadding(Side::Left);
  const int right_padding = plotAreaPadding(Side::Right);
  const int top_padding = plotAreaPadding(Side::Top);
  const int width = painter.window().width() - left_padding - right_padding;
  const int height = painter.window().height()-plotAreaPadding(Side::Top)-plotAreaPadding(Side::Bottom);
  
  try
  {
    //TextWordWrap not supported by all painters
    WRectF textArea( left_padding, top_padding, width, height);
    painter.drawText( textArea, AlignmentFlag::Center | AlignmentFlag::Middle,
                     TextFlag::WordWrap, m_textInMiddleOfChart);
  }catch(...)
  {
    const double lineHeight = 40.0;
    const double fontlen = txtFont.sizeLength().toPixels();
    const size_t ncharacters = static_cast<size_t>( 2.0 * width / fontlen);
    string text = m_textInMiddleOfChart.narrow();
    
    if( text.size() < ncharacters )
    {
      const double fromTop = top_padding + 0.5*height - 0.5*lineHeight;
      WRectF textArea( left_padding, fromTop, width, lineHeight);
      painter.drawText( textArea, AlignmentFlag::Center | AlignmentFlag::Middle, m_textInMiddleOfChart);
    }else
    {
      text = m_textInMiddleOfChart.toUTF8();
      vector<string> words, lines;
      SpecUtils::split( words, text, " \t");
      string line;
      for( const string &word : words )
      {
        if( line.empty() || ((line.size()+1+word.size()) < ncharacters ) )
        {
          if( !line.empty() )
            line += " ";
          line += word;
        }else
        {
          lines.push_back( line);
          line = word;
        }//if( theres room for this word ) / else
      }//for( const string &word : words )
      
      if( !line.empty() )
        lines.push_back( line);
      
      const double initialHeight = top_padding + 0.5*height
      - 0.5*lineHeight*lines.size();
      for( size_t i = 0; i < lines.size(); ++i )
      {
        const double fromTop = initialHeight + lineHeight*i;
        //          WRectF textArea( left_padding + 20.0, fromTop, width, lineHeight);
        //          painter.drawText( textArea, AlignmentFlag::Left | AlignmentFlag::Middle, lines[i]);
        WRectF textArea( left_padding, fromTop, width, lineHeight);
        painter.drawText( textArea, AlignmentFlag::Center | AlignmentFlag::Middle, lines[i]);
      }//for( size_t i = 0; i < lines.size(); ++i )
    }//if( text.size() < ncharacters ) / else
  }//try / catch
  
  painter.restore();
  painter.setFont( oldFont);
}//void paintTextInMiddleOfChart(...)



void SpectrumChart::paint( WPainter &painter, const WRectF &rectangle ) const
{
  // WCartesianChart only sets this from paintEvent(...); without it, SVGs get jagged curves.
  painter.setRenderHint( RenderHint::Antialiasing, true );

  // Draws the background, spectra, and axes (see render(...))
  Chart::WCartesianChart::paint( painter, rectangle);
  
  paintTextInMiddleOfChart( painter);
  
  //Check if data is both positive and negative, and if so, draw a grey line at
  //  y==0; this is mostly for looking at second derivatives in comparison to
  //  to the spectrum (a development feature).
  if( axis(Chart::Axis::Y1).scale() == Chart::AxisScale::Linear )
  {
    const double ymin = axis(Chart::Axis::Y1).minimum();
    const double ymax = axis(Chart::Axis::Y1).maximum();
    if( ymax > 0.0 && ymin < -0.05*ymax )
    {
      const double x_min = axis(Chart::Axis::X).minimum();
      const double x_max = axis(Chart::Axis::X).maximum();
      const WPointF leftpoint = mapToDevice( x_min, 0.0);
      const WPointF rightpoint = mapToDevice( x_max, 0.0);
      
      painter.save();
      painter.setPen( WPen(Wt::StandardColor::Gray));
      painter.drawLine( leftpoint, rightpoint);
      painter.restore();
    }//if( ymax > 0.0 && ymin < -0.01*ymax )
  }//if( axis(Chart::Axis::Y1).scale() == Chart::AxisScale::Linear )
  
  paintHighlightRegions( painter);

  paintPeaks( painter);

  // Same stacking order as SpectrumChartD3.js: grid over the spectrum and peaks, then reference
  //  lines, and the legend on top.
  renderGridLines( painter, Chart::Axis::X );
  renderGridLines( painter, Chart::Axis::Y1 );
  renderGridLines( painter, Chart::Axis::Y2 );

  renderReferncePhotoPeakLines( painter);

  paintOnChartLegend( painter);

  paintTitle( painter );
}//paint( ... )


void SpectrumChart::paintTitle( Wt::WPainter &painter ) const
{
  // Wt draws the title from `renderLegend(...)`, which we override, so we draw it here.
  if( title().empty() )
    return;

  const WFont font = titleFont();
  const double fontheight = font.sizeLength( 16 ).toPixels();
  const double height = std::max( chartArea().top(), fontheight );

  painter.save();
  painter.setFont( font );
  painter.setPen( m_textPen );
  painter.setBrush( WBrush() );
  painter.drawText( WRectF( 0.0, 0.0, painter.window().width(), height ),
                    AlignmentFlag::Center | AlignmentFlag::Middle, title() );
  painter.restore();
}//void paintTitle( Wt::WPainter &painter ) const


void SpectrumChart::paintHighlightRegions( Wt::WPainter &painter ) const
{
  const WFont  oldFont  = painter.font();
  const WPen   oldPen   = painter.pen();
  const WBrush oldBrush = painter.brush();
  
  for( size_t i = 0; i < m_highlights.size(); ++i )
  {
    const HighlightRegion &reg = m_highlights[i];
    const double y_min = axis(Chart::Axis::Y1).minimum();
    const double y_max = axis(Chart::Axis::Y1).maximum();
    
    const WPointF topLeft = mapToDevice( reg.lowerx, y_max);
    WPointF bottomRight = mapToDevice( reg.upperx, y_min);
    
    if( fabs(bottomRight.x()-topLeft.x()) < 2.0 )
    {
      const double mult = (reg.lowerx>reg.upperx) ? -1.0 : 1.0;
      bottomRight.setX( bottomRight.x() + mult*2.0);
    }
    
    painter.fillRect( WRectF(topLeft, bottomRight), WBrush(reg.color));
  }//for( size_t i = 0; i < m_highlights.size(); ++i )
  
  for( const auto &reg : m_occupancy_regions )
  {
    const double y_min_counts = axis(Chart::Axis::Y1).minimum();
    WPointF bottomLeft = mapToDevice( reg.lowerx, y_min_counts);
    WPointF bottomRight = mapToDevice( reg.upperx, y_min_counts);
    
    const double bottomPx = chartArea().bottom() + 0.5;// + axis.margin();
    const double topPx = chartArea().top();// + axis.margin();
    
    //If we wanted to shade in the x-axis
    //bottomLeft.setY( bottomPx);
    //bottomRight.setY( bottomPx + 5);
    //painter.fillRect( WRectF(bottomLeft, bottomRight), WBrush(reg.color));
    
    painter.setPen( WPen(reg.color));
    painter.setBrush( WBrush(reg.color));
    
    const float leftPx = bottomLeft.x();
    const float rightPx = bottomRight.x();
    painter.drawLine( WPointF(leftPx,bottomPx), WPointF(leftPx,topPx));
    painter.drawLine( WPointF(rightPx,bottomPx), WPointF(rightPx,topPx));
    
    if( fabs(rightPx - leftPx) > 50 )
    {
      WFont font( FontFamily::Default);
      font.setSize( 8);
      painter.setFont( font);
      painter.drawText(leftPx+2, topPx, 30, 10.0, AlignmentFlag::Left, "occ. start");
      painter.drawText(rightPx-2-30, topPx, 30, 10.0, AlignmentFlag::Right,  "occ. end");
    }
  }//for( const auto region : m_underlines )
  
  painter.setPen( oldPen);
  painter.setFont( oldFont);
  painter.setBrush( oldBrush);
}//void paintHighlightRegions( Wt::WPainter& painter ) const


bool SpectrumChart::roiCurves( const std::vector<std::shared_ptr<const PeakDef>> &roi_peaks,
                               const std::shared_ptr<const SpecUtils::Measurement> &data,
                               const std::vector<double> &edges,
                               const std::vector<double> &counts,
                               const std::vector<double> *background,
                               const double xmin, const double xmax, const double px_per_kev,
                               RoiCurves &curves )
{
  // A port of `roiGeometry` in SpectrumChartD3.js; see the comments there too.
  curves = RoiCurves();

  const size_t npts = counts.size();
  if( roi_peaks.empty() || !roi_peaks.front() || (npts < 2) || (edges.size() != (npts + 1))
     || !(px_per_kev > 0.0) || (background && (background->size() != npts)) )
    return false;

  const shared_ptr<const PeakContinuum> continuum = roi_peaks.front()->continuum();
  if( !continuum )
    return false;

  const double xa = std::max( roi_peaks.front()->lowerX(), xmin );
  const double xb = std::min( roi_peaks.back()->upperX(), xmax );
  if( !(xb > xa) )
    return false;

  // The energy range each peak is drawn over: Gaussians out to an amplitude-dependent number of
  //  sigma, and skewed peaks until their tails are negligible.
  const size_t npeaks = roi_peaks.size();
  const double inf = std::numeric_limits<double>::infinity();
  vector<pair<double,double>> extents( npeaks, {-inf, inf} );
  double min_sigma = inf;
  for( size_t j = 0; j < npeaks; ++j )
  {
    const PeakDef &peak = *roi_peaks[j];
    if( !peak.gausPeak() )
      continue;

    if( peak.skewType() == PeakDef::NoSkew )
    {
      const double amp = std::max( 1.0, std::fabs( peak.amplitude() ) );
      const double nsigma = std::max( 4.0, std::sqrt( 2.0*std::log( 4.0E3*amp ) ) );
      extents[j] = { peak.mean() - nsigma*peak.sigma(), peak.mean() + nsigma*peak.sigma() };
    }else if( const optional<pair<double,double>> range = peak.visRange() )
    {
      extents[j] = *range;
    }

    if( (extents[j].first <= xb) && (extents[j].second >= xa) )
      min_sigma = std::min( min_sigma, peak.sigma() );
  }//for( size_t j = 0; j < npeaks; ++j )

  // Samples are spaced about a pixel apart, but no closer than sigma/16 of the narrowest peak.
  const double step_px = std::min( 20.0, std::max( 1.0, min_sigma*px_per_kev/16.0 ) );

  const auto bin_width = [&edges]( const size_t k ) -> double { return edges[k+1] - edges[k]; };
  const auto bin_at = [&]( const double x ) -> size_t {
    const ptrdiff_t k = (std::upper_bound( begin(edges), begin(edges) + npts, x ) - begin(edges)) - 1;
    return static_cast<size_t>( std::clamp( k, ptrdiff_t(0), static_cast<ptrdiff_t>(npts - 1) ) );
  };

  // Sample energies `xs`, and the displayed channel width at each, `ws`.  Samples are anchored to
  //  channels (an odd number per channel, so its center is included), with the last sample at or
  //  before xa, and the first at or after xb, kept so the curves can be clipped to [xa, xb].
  vector<double> &xs = curves.energies;
  vector<double> ws;
  xs.push_back( xa );
  ws.push_back( bin_width( bin_at(xa) ) );
  bool done = false;
  for( size_t k = (bin_at(xa) > 0) ? (bin_at(xa) - 1) : 0; !done && (k < npts); ++k )
  {
    const double lo = edges[k], w = bin_width(k);
    if( !(w > 0.0) )
      continue;

    // (just the center, for channels too wide to resolve the peaks)
    const size_t m = (w > 2.0*min_sigma) ? 1 : (1 + 2*static_cast<size_t>( std::floor( 0.5*w*px_per_kev/step_px ) ));
    const double dx = w / m;
    const double s0 = std::min( m - 1.0, std::max( 0.0, std::floor( (xa - lo)/dx - 0.5 ) ) );
    for( size_t s = static_cast<size_t>( s0 ); !done && (s < m); ++s )
    {
      const double x = lo + (s + 0.5)*dx;
      if( x <= xa )
      {
        xs[0] = x;
        ws[0] = w;
      }else
      {
        xs.push_back( x );
        ws.push_back( w );
        done = (x >= xb);
      }
    }//for( loop over samples in this channel )
  }//for( loop over displayed channels )

  if( !done )
  {
    xs.push_back( xb );
    ws.push_back( bin_width( bin_at(xb) ) );
  }

  const size_t n = xs.size();
  assert( n >= 2 );

  // Each samples window, [x - w/2, x + w/2]
  vector<double> win_lo( n ), win_hi( n );
  for( size_t i = 0; i < n; ++i )
  {
    win_lo[i] = xs[i] - 0.5*ws[i];
    win_hi[i] = xs[i] + 0.5*ws[i];
  }

  // The end samples get moved onto xa and xb, along the line to their neighbor.
  const double fa = (xa - xs[0]) / (xs[1] - xs[0]);
  const double fb = (xb - xs[n-2]) / (xs[n-1] - xs[n-2]);
  xs[0] = xa;
  xs[n-1] = xb;

  // Sets out[d0..d1] to the counts in each samples window, for samples i0 through i1 (zero for
  //  others), with the end samples moved onto xa and xb.
  const auto window_counts = [&]( const std::function<double(double,double)> &integral, vector<double> &out,
                                  const size_t d0, const size_t d1, const size_t i0, const size_t i1 ){
    for( size_t i = d0; i <= d1; ++i )
      out[i] = ((i < i0) || (i > i1)) ? 0.0 : std::max( 0.0, integral( win_lo[i], win_hi[i] ) );
    const double v0 = out[0], v1 = out[1], vm = out[n-2], vn = out[n-1];
    if( d0 == 0 )
      out[0] = v0 + fa*(v1 - v0);
    if( d1 == (n - 1) )
      out[n-1] = vm + fb*(vn - vm);
  };//window_counts

  // Content of displayed channels at energy x, linearly interpolated between channel centers.
  const auto hist_at = [&]( const double x, const vector<double> &values ) -> double {
    const auto center = [&edges]( const size_t k ){ return 0.5*(edges[k] + edges[k+1]); };
    size_t k = bin_at( x );
    if( (k > 0) && (x < center(k)) )
      --k;
    if( k == (npts - 1) )
      return values[k];
    const double frac = std::clamp( (x - center(k)) / (center(k+1) - center(k)), 0.0, 1.0 );
    return values[k] + frac*(values[k+1] - values[k]);
  };//hist_at

  // The background is taken at channel centers, the same way the background-subtracted data is.
  const auto bg_at = [&]( const double x ) -> double {
    return background ? hist_at( x, *background ) : 0.0;
  };

  // The continuum: step and External continua are defined channel-by-channel, so are drawn from
  //  the cumulative sum of their channel counts; polynomials are integrated directly.
  std::function<double(double,double)> continuum_integral;
  vector<float> cont_energies;
  vector<double> cont_cumulative;
  if( PeakDef::continuum_channel_counts( roi_peaks, data, cont_energies, cont_cumulative ) )
  {
    const size_t nchan = cont_cumulative.size();
    cont_cumulative.insert( begin(cont_cumulative), 0.0 );
    for( size_t i = 1; i <= nchan; ++i )
      cont_cumulative[i] += cont_cumulative[i-1];

    const auto cdf = [nchan, &cont_energies, &cont_cumulative]( const double x ) -> double {
      const ptrdiff_t k = (std::upper_bound( begin(cont_energies), end(cont_energies), x ) - begin(cont_energies)) - 1;
      const size_t i = static_cast<size_t>( std::clamp( k, ptrdiff_t(0), static_cast<ptrdiff_t>(nchan - 1) ) );
      const double density = (cont_cumulative[i+1] - cont_cumulative[i]) / (cont_energies[i+1] - cont_energies[i]);
      return cont_cumulative[i] + (x - cont_energies[i])*density;
    };
    continuum_integral = [cdf]( const double x0, const double x1 ){ return cdf(x1) - cdf(x0); };
  }else if( continuum->isPolynomial() && !PeakContinuum::is_step_continuum( continuum->type() ) )
  {
    continuum_integral = [continuum]( const double x0, const double x1 ){
      return continuum->offset_integral( x0, x1, nullptr, nullptr, 0 );
    };
  }//if( step or External continuum ) / else if( polynomial )

  // The continuum `c`, minus the background
  vector<double> &c = curves.continuum;
  c.assign( n, 0.0 );
  if( continuum_integral )
    window_counts( continuum_integral, c, 0, n - 1, 0, n - 1 );
  for( size_t i = 0; i < n; ++i )
    c[i] -= bg_at( xs[i] );

  // First sample index for which `pred(i)` is true, for a `pred` that is false and then true.
  const auto first_true = [n]( const std::function<bool(size_t)> &pred ) -> size_t {
    size_t lo = 0, hi = n;
    while( lo < hi )
    {
      const size_t mid = (lo + hi) / 2;
      if( pred(mid) )
        hi = mid;
      else
        lo = mid + 1;
    }
    return lo;
  };//first_true

  curves.peaks.resize( npeaks );
  for( size_t j = 0; j < npeaks; ++j )
  {
    const PeakDef &peak = *roi_peaks[j];
    RoiCurves::PeakCurve &curve = curves.peaks[j];

    // The samples whose windows the peak is in (the only peak in a ROI is drawn across the whole
    //  ROI); its fill and outline span one more sample each side, to start and end on what it
    //  sits on.
    const bool whole = (!peak.gausPeak() || (npeaks == 1));
    const double lo = extents[j].first, hi = extents[j].second;
    const size_t i0 = whole ? 0 : first_true( [&]( size_t i ){ return win_hi[i] >= lo; } );
    const size_t i1_end = whole ? n : first_true( [&]( size_t i ){ return win_lo[i] > hi; } );
    if( (i1_end == 0) || (i0 > (i1_end - 1)) )
      continue;
    const size_t i1 = i1_end - 1;

    curve.first = (i0 > 0) ? (i0 - 1) : size_t(0);
    curve.last = std::min( n - 1, i1 + 1 );
    curve.counts.assign( n, 0.0 );

    if( peak.gausPeak() )
    {
      window_counts( [&peak]( const double x0, const double x1 ){ return peak.gauss_integral( x0, x1 ); },
                     curve.counts, curve.first, curve.last, i0, i1 );
    }else
    {
      // A data-defined peak is the data above the continuum
      for( size_t i = 0; i < n; ++i )
        curve.counts[i] = std::max( 0.0, hist_at( xs[i], counts ) - bg_at( xs[i] ) - c[i] );
    }
  }//for( size_t j = 0; j < npeaks; ++j )

  return true;
}//bool roiCurves(...)


int SpectrumChart::displayRebinFactor( const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                       const double xmin, const double xmax,
                                       const double plot_width_px, const double line_width )
{
  const size_t nchannel = meas ? meas->num_gamma_channels() : size_t(0);
  if( (nchannel < 2) || !(plot_width_px > 2.0) || !(xmax > xmin) )
    return 1;

  const size_t first = meas->find_gamma_channel( static_cast<float>( xmin ) );
  const size_t last = meas->find_gamma_channel( static_cast<float>( xmax ) );
  const double nvisible = (last > first) ? static_cast<double>( last - first ) : 1.0;
  const double factor = std::ceil( std::max( line_width, 0.0 ) * nvisible / plot_width_px );

  return std::max( 1, static_cast<int>( std::min( factor, 0.5*nchannel ) ) );
}//int displayRebinFactor(...)


void SpectrumChart::paintPeaks( Wt::WPainter& painter ) const
{
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = m_peakModel ? m_peakModel->peaks() : m_peaks;
  const SpectrumDataModel * const model = spectrum_model( *this );
  if( !peaks || peaks->empty() || !model )
    return;

  const shared_ptr<const Measurement> data = model->getData();
  const int nrow = model->rowCount();
  if( !data || (nrow < 2) )
    return;

  const Chart::WAxis &yaxis = axis( Chart::Axis::Y1 );
  const double xmin = axis( Chart::Axis::X ).minimum();
  const double xmax = axis( Chart::Axis::X ).maximum();
  const double ymin = yaxis.minimum();
  if( !(xmax > xmin) || !(yaxis.maximum() > ymin) )
    return;

  // Group the peaks into ROIs (by continuum), ordered by energy
  vector<vector<shared_ptr<const PeakDef>>> rois;
  {
    vector<shared_ptr<const PeakDef>> sorted;
    for( const shared_ptr<const PeakDef> &peak : *peaks )
    {
      if( peak && peak->continuum() )
        sorted.push_back( peak );
    }
    std::sort( begin(sorted), end(sorted), &PeakDef::lessThanByMeanShrdPtr );

    std::map<const PeakContinuum *, size_t> roi_index;
    for( const shared_ptr<const PeakDef> &peak : sorted )
    {
      const auto pos = roi_index.emplace( peak->continuum().get(), rois.size() );
      if( pos.second )
        rois.emplace_back();
      rois[pos.first->second].push_back( peak );
    }

    std::sort( begin(rois), end(rois), []( const auto &lhs, const auto &rhs ){
      return lhs.front()->lowerX() < rhs.front()->lowerX();
    } );
  }

  // The displayed channels, from one before the chart to one past it
  const int row0 = std::clamp( model->findRow( xmin ) - 1, 0, nrow - 1 );
  const int row1 = std::clamp( model->findRow( xmax ) + 1, 0, nrow - 1 );
  const bool subtract = (model->backgroundSubtract() && model->getBackground());
  const auto value_of = [model]( const int row, const SpectrumDataModel::ColumnType column ) -> double {
    const double value = asNumber( model->displayBinValue( row, column ) );
    return (IsNan(value) || IsInf(value)) ? 0.0 : value;
  };

  vector<double> edges, counts, background;
  for( int row = row0; row <= row1; ++row )
  {
    edges.push_back( model->rowLowEdge( row ) );
    counts.push_back( value_of( row, SpectrumDataModel::DATA_COLUMN ) );
    if( subtract )
      background.push_back( value_of( row, SpectrumDataModel::BACKGROUND_COLUMN ) );
  }
  edges.push_back( model->rowLowEdge( row1 ) + model->rowWidth( row1 ) );

  // Mapping to pixels: on a log axis, values off the chart go to its edges (as d3's clamped log
  //  scale does), and on a linear axis values are kept within reach of the chart (for the paint
  //  devices' sake; it is clipped anyway).
  const WRectF &area = chartArea();
  const double x0px = mapToDevice( xmin, ymin ).x();
  const double px_per_kev = (mapToDevice( xmax, ymin ).x() - x0px) / (xmax - xmin);
  const double ymargin = (yaxis.scale() == Chart::AxisScale::Log) ? 0.0 : area.height();
  const auto ypx = [&]( const double y ) -> double {
    const double py = mapToDevice( xmin, y ).y();
    if( IsNan(py) || IsInf(py) )
      return area.bottom();
    return std::clamp( py, area.top() - ymargin, area.bottom() + ymargin );
  };

  const auto peak_pen = []( const WColor &color ) -> WPen {
    WPen pen( color );
    pen.setWidth( WLength(1.0) );
    pen.setCapStyle( PenCapStyle::Flat );
    pen.setJoinStyle( PenJoinStyle::Miter );
    return pen;
  };

  const bool showlabels = (showingPeakLabel( PeakLabelType::ShowPeakUserLabel )
                           || showingPeakLabel( PeakLabelType::ShowPeakEnergyLabel )
                           || showingPeakLabel( PeakLabelType::ShowPeakNuclideLabel ));
  vector<pair<shared_ptr<const PeakDef>,WPointF>> labels;

  painter.save();
  WPainterPath clip;
  clip.addRect( area );
  painter.setClipPath( clip );
  painter.setClipping( true );

  for( const vector<shared_ptr<const PeakDef>> &roi : rois )
  {
    RoiCurves curves;
    if( !roiCurves( roi, data, edges, counts, subtract ? &background : nullptr,
                    xmin, xmax, px_per_kev, curves ) )
      continue;

    const size_t n = curves.energies.size();
    const vector<double> &c = curves.continuum;

    vector<double> xpx( n ), cpx( n );
    for( size_t i = 0; i < n; ++i )
    {
      xpx[i] = x0px + px_per_kev*(curves.energies[i] - xmin);
      cpx[i] = ypx( c[i] );
    }

    // Peaks are stacked in ROI order: each peaks fill sits on top of the continuum plus the earlier
    //  peaks, while (in multi-peak ROIs) each peaks outline sits on just the continuum.
    vector<double> stack = c, stack_px = cpx, top_px( n );
    vector<pair<WPainterPath,WColor>> fills, outlines;
    const bool multi = (roi.size() > 1);

    for( size_t j = 0; j < roi.size(); ++j )
    {
      const RoiCurves::PeakCurve &curve = curves.peaks[j];
      if( curve.first > curve.last )
        continue;

      const PeakDef &peak = *roi[j];
      const vector<double> &p = curve.counts;
      const size_t d0 = curve.first, d1 = curve.last;
      const WColor color = peak.lineColor().isDefault() ? m_defaultPeakColor : peak.lineColor();

      for( size_t i = d0; i <= d1; ++i )
        top_px[i] = ypx( stack[i] + p[i] );

      WPainterPath fill;
      append_thinned( fill, xpx, [&stack_px]( size_t i ){ return stack_px[i]; }, d0, d1, false, true );
      append_thinned( fill, xpx, [&top_px]( size_t i ){ return top_px[i]; }, d0, d1, true, false );
      fill.closeSubPath();
      fills.emplace_back( fill, color );

      if( multi )
      {
        WPainterPath outline;
        append_thinned( outline, xpx, [&cpx]( size_t i ){ return cpx[i]; }, d0, d0, false, true );
        append_thinned( outline, xpx, [&]( size_t i ){ return ypx( c[i] + p[i] ); }, d0, d1, false, false );
        append_thinned( outline, xpx, [&cpx]( size_t i ){ return cpx[i]; }, d1, d1, false, false );
        outlines.emplace_back( outline, color );
      }//if( multi )

      for( size_t i = d0; i <= d1; ++i )
      {
        stack[i] += p[i];
        stack_px[i] = top_px[i];
      }

      // Labels are placed above the continuum plus this peak, at its centroid
      const double mean = peak.mean();
      if( showlabels && (mean >= curves.energies.front()) && (mean <= curves.energies.back()) )
      {
        const size_t i = std::max( size_t(1), static_cast<size_t>(
                            std::lower_bound( begin(curves.energies), end(curves.energies), mean ) - begin(curves.energies) ) );
        const double f = (mean - curves.energies[i-1]) / (curves.energies[i] - curves.energies[i-1]);
        const double cont = c[i-1] + f*(c[i] - c[i-1]);
        const double pk = p[i-1] + f*(p[i] - p[i-1]);
        labels.emplace_back( roi[j], WPointF( x0px + px_per_kev*(mean - xmin), ypx( cont + pk ) ) );
      }
    }//for( size_t j = 0; j < roi.size(); ++j )

    for( const pair<WPainterPath,WColor> &fill : fills )
    {
      painter.setPen( peak_pen( fill.second ) );
      painter.setBrush( WBrush( with_alpha( fill.second, ns_peak_fill_opacity ) ) );
      painter.drawPath( fill.first );
    }

    for( const pair<WPainterPath,WColor> &outline : outlines )
      painter.strokePath( outline.first, peak_pen( outline.second ) );

    // The continuum line - unless a single peak's fill outline already draws it.
    if( multi || (fills.size() != 1) || !same_color( fills.front().second, m_defaultPeakColor ) )
    {
      WPainterPath line;
      append_thinned( line, xpx, [&cpx]( size_t i ){ return cpx[i]; }, 0, n - 1, false, true );
      painter.strokePath( line, peak_pen( m_defaultPeakColor ) );
    }
  }//for( const vector<shared_ptr<const PeakDef>> &roi : rois )

  for( const auto &label : labels )
    drawPeakText( *label.first, painter, label.second );

  painter.restore();
}//void paintPeaks( Wt::WPainter& painter ) const


void SpectrumChart::drawPeakText( const PeakDef &peak, Wt::WPainter &painter,
                                  Wt::WPointF peakTip ) const
{
  // The rows of text, as SpectrumChartD3.js has them
  vector<pair<string,PeakLabelType>> rows;
  char buffer[64] = { '\0' };

  if( showingPeakLabel( PeakLabelType::ShowPeakUserLabel ) && !peak.userLabel().empty() )
    rows.emplace_back( peak.userLabel(), PeakLabelType::ShowPeakUserLabel );

  if( showingPeakLabel( PeakLabelType::ShowPeakEnergyLabel ) )
  {
    snprintf( buffer, sizeof(buffer), "%.2f keV", peak.mean() );
    rows.emplace_back( buffer, PeakLabelType::ShowPeakEnergyLabel );
  }

  if( showingPeakLabel( PeakLabelType::ShowPeakNuclideLabel ) )
  {
    string txt;
    double energy = 0.0;
    // The energies are of the source gamma (e.g., not 511 keV less, for a single escape peak), as
    //  `PeakDef::gaus_peaks_to_json(...)` gives them to the D3 chart.
    bool have_energy = true;
    if( peak.parentNuclide() )
    {
      txt = peak.parentNuclide()->symbol;
      const SandiaDecay::RadParticle * const particle = peak.decayParticle();
      if( peak.sourceGammaType() == PeakDef::AnnihilationGamma )
        energy = 510.9989;
      else if( particle )
        energy = particle->energy;
      else
        have_energy = false;
    }else if( peak.xrayElement() )
    {
      txt = peak.xrayElement()->name;
      energy = peak.xrayEnergy();
    }else if( peak.reaction() )
    {
      txt = peak.reaction()->name();
      energy = peak.reactionEnergy();
    }

    if( !txt.empty() && have_energy && showingPeakLabel( PeakLabelType::ShowPeakNuclideEnergies ) )
    {
      const char *type = "";
      switch( peak.sourceGammaType() )
      {
        case PeakDef::NormalGamma:       type = "";        break;
        case PeakDef::AnnihilationGamma: type = " annih."; break;
        case PeakDef::XrayGamma:         type = " x-ray";  break;
        case PeakDef::DoubleEscapeGamma: type = " D.E.";   break;
        case PeakDef::SingleEscapeGamma: type = " S.E.";   break;
      }//switch( peak.sourceGammaType() )

      snprintf( buffer, sizeof(buffer), ", %.2f keV%s", energy, type );
      txt += buffer;
    }//if( show nuclide energies )

    if( !txt.empty() )
      rows.emplace_back( txt, PeakLabelType::ShowPeakNuclideLabel );
  }//if( showing nuclide labels )

  if( rows.empty() )
    return;

  // The rows are centered on the peak, the last ending 10 px above it
  const WFont font = axis( Chart::Axis::Y1 ).labelFont();
  const double fontheight = font.sizeLength( 12 ).toPixels();

  painter.save();
  painter.setFont( font );
  painter.setBrush( WBrush() );
  for( size_t i = 0; i < rows.size(); ++i )
  {
    const uint32_t c = m_peakLabelsColors[static_cast<int>(rows[i].second)];
    painter.setPen( c ? WPen( WColor( c & 0xFF, (c >> 8) & 0xFF, (c >> 16) & 0xFF, (c >> 24) & 0xFF ) ) : m_textPen );

    const double bottom = peakTip.y() - 10.0 - (rows.size() - 1 - i)*fontheight;
    painter.drawText( WRectF( peakTip.x() - 500.0, bottom - 1.5*fontheight, 1000.0, 1.5*fontheight ),
                      AlignmentFlag::Center | AlignmentFlag::Bottom, WString::fromUTF8( rows[i].first ) );
  }//for( size_t i = 0; i < rows.size(); ++i )
  painter.restore();
}//void drawPeakText(...)



void SpectrumChart::paintOnChartLegend( Wt::WPainter &painter ) const
{
  // A port of SpectrumChartD3.js's legend (`updateLegend()`): a translucent box in the upper right,
  //  with each spectrum's line, title and counts, live and real times, scale factor, and neutrons.
  if( m_legendType != OnChartLegend )
    return;

  const SpectrumDataModel * const mdl = spectrum_model( *this );
  if( !mdl )
    return;

  const int nhists = (!!mdl->getData() + !!mdl->getSecondData() + !!mdl->getBackground());

  //For time series data, dont show the legend for gamma only data
  if( (nhists < 1) || ((m_xAxisUnits == XAxisUnits::kSeconds) && (nhists < 2)) )
    return;

  struct Entry
  {
    WPen pen;
    vector<WString> rows;  // The first row is the title
  };//struct Entry

  vector<Entry> entries;

  const auto add_entry = [&]( const shared_ptr<const Measurement> &hist, const int column,
                              const double sf, const char *default_title ){
    if( !hist )
      return;

    Entry entry;
    try
    {
      entry.pen = series( column ).pen();
    }catch( std::exception & )
    {
    }

    const double nsum = sf * hist->gamma_count_sum();
    string title = hist->title().empty() ? string( default_title ) : hist->title();
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), (nsum > 1000.0) ? " (%.0f counts)" : " (%.1f counts)", nsum );
    entry.rows.push_back( WString::fromUTF8( title + buffer ) );

    const double lt = hist->live_time(), rt = hist->real_time();
    if( lt > 0.0 )
      entry.rows.push_back( WString::fromUTF8( "Live Time: " + js_to_precision( sf*lt, 4 ) + " s" ) );
    if( rt > 0.0 )
      entry.rows.push_back( WString::fromUTF8( "Real Time: " + js_to_precision( sf*rt, 4 ) + " s" ) );
    if( std::fabs( sf - 1.0 ) > 1.0E-6 )
      entry.rows.push_back( WString::fromUTF8( "Scaled by " + js_to_precision( sf, 4 ) ) );

    if( hist->contained_neutron() )
    {
      const double nrt = (hist->neutron_live_time() > 0.0f) ? hist->neutron_live_time()
                                                             : ((rt > 1.0E-6) ? rt : lt);
      const double neutrons = hist->neutron_counts_sum();
      if( nrt > 0.0 )
        entry.rows.push_back( WString::fromUTF8( "Neutrons: " + legend_rate_str( neutrons/nrt, 3 ) + " cps" ) );
      else
        entry.rows.push_back( WString::fromUTF8( "Neutrons: " + legend_rate_str( sf*neutrons, 3 ) ) );
    }//if( hist->contained_neutron() )

    entries.push_back( entry );
  };//add_entry

  // In the same order as SpectrumChartD3.js
  add_entry( mdl->getData(), mdl->dataColumn(), 1.0, "Foreground" );
  add_entry( mdl->getSecondData(), mdl->secondDataColumn(), mdl->secondDataScaledBy(), "Secondary" );
  add_entry( mdl->getBackground(), mdl->backgroundColumn(), mdl->backgroundScaledBy(), "Background" );

  if( entries.empty() )
    return;

  painter.save();

  const WFont font = axis( Chart::Axis::Y1 ).labelFont();
  painter.setFont( font );

  const double fontheight = font.sizeLength( 12 ).toPixels();
  const double rowheight = 1.2*fontheight;
  const double entry_spacing = 5.0;

  // The body: title rows start 15 px in, the rest 20 px
  double body_width = 0.0, body_height = 0.0;
  for( const Entry &entry : entries )
  {
    for( size_t i = 0; i < entry.rows.size(); ++i )
      body_width = std::max( body_width, (i ? 20.0 : 15.0) + text_width( painter, entry.rows[i] ) );
    body_height += rowheight*entry.rows.size();
  }
  body_height += entry_spacing*(entries.size() - 1);

  // The box is 30 px in from the right of the image, and 10 px below the top of the chart
  const double box_width = body_width + 15.0, box_height = body_height + 10.0;
  const double image_width = painter.window().width();
  const double box_x = std::max( chartArea().left(), image_width - 30.0 - box_width );
  const double box_y = chartArea().top() + 10.0;

  // A rounded rectangle (using quadratic curves, which every paint device supports)
  const double r = 5.0;
  WPainterPath box;
  box.moveTo( box_x + r, box_y );
  box.lineTo( box_x + box_width - r, box_y );
  box.quadTo( box_x + box_width, box_y, box_x + box_width, box_y + r );
  box.lineTo( box_x + box_width, box_y + box_height - r );
  box.quadTo( box_x + box_width, box_y + box_height, box_x + box_width - r, box_y + box_height );
  box.lineTo( box_x + r, box_y + box_height );
  box.quadTo( box_x, box_y + box_height, box_x, box_y + box_height - r );
  box.lineTo( box_x, box_y + r );
  box.quadTo( box_x, box_y, box_x + r, box_y );
  box.closeSubPath();

  const WBrush &chart_brush = background();
  const WColor fill = ((chart_brush.style() == BrushStyle::Solid) && !chart_brush.color().isDefault())
                      ? with_alpha( chart_brush.color(), 0.5 )
                      : ns_legend_fill_color;
  painter.setPen( WPen( ns_legend_stroke_color ) );
  painter.setBrush( WBrush( fill ) );
  painter.drawPath( box );

  // The entries
  painter.setBrush( WBrush() );
  double y = box_y + 6.0;
  const double x = box_x + 8.0;
  for( const Entry &entry : entries )
  {
    WPen swatch_pen = entry.pen;
    swatch_pen.setCapStyle( PenCapStyle::Flat );
    painter.setPen( swatch_pen );
    painter.drawLine( x, y + 0.5*rowheight - 1.0, x + 12.0, y + 0.5*rowheight - 1.0 );

    painter.setPen( m_textPen );
    for( size_t i = 0; i < entry.rows.size(); ++i )
    {
      const double left = x + (i ? 20.0 : 15.0);
      painter.drawText( WRectF( left, y, image_width, rowheight ),
                        AlignmentFlag::Left | AlignmentFlag::Middle, entry.rows[i] );
      y += rowheight;
    }
    y += entry_spacing;
  }//for( const Entry &entry : entries )

  painter.restore();
}//void paintOnChartLegend()


void SpectrumChart::paintEvent( WPaintDevice *paintDevice )
{
  //This function is called before SpectrumChart::paint(...)
  m_widthInPixels  = paintDevice->width().toPixels();
  m_heightInPixels = paintDevice->height().toPixels();

#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
  {
    const double w = paintDevice->width().toPixels();
    const double h = paintDevice->height().toPixels();
    const int status = setLeftYAxisPadding( w, h, paintDevice);
    if( status < 0 )
      cerr << "SpectrumChart::paintEvent() warning: failed to adjust left axis"
           << endl;
  }
#endif

  //paint the chart
  Wt::Chart::WCartesianChart::paintEvent( paintDevice);
}//void paintEvent( WPaintDevice *paintDevice )


void SpectrumChart::clearTimeHighlightRegions( const HighlightRegionType region )
{
  setTimeHighLightRegions( vector< pair<double,double> >(), region);
}//void clearTimeHighlightRegions()


void SpectrumChart::setYAxisScale( Wt::Chart::AxisScale scale )
{
  axis( Chart::Axis::Y1 ).setScale( scale);
}//void setYAxisScale( Wt::Chart::AxisScale scale )


void SpectrumChart::showHistogramIntegralsInLegend( const bool show )
{
  SpectrumDataModel * const the_model = spectrum_model( *this );
  assert( the_model );
  
  if( the_model )
    the_model->addIntegralOfHistogramToLegend( show);
}//void showHistogramIntegralsInLegend( const bool show = true);


bool SpectrumChart::legendIsEnabled() const
{
  return (m_legendType != NoLegend);
}//bool legendIdEnabled() const


void SpectrumChart::disableLegend()
{
  if( m_legendType == NoLegend )
    return;

  m_legendType = NoLegend;
  update();
}//void disableLegend()


void SpectrumChart::enableLegend()
{
  if( m_legendType != NoLegend )
    return;

  m_legendType = OnChartLegend;
  update();
}//void enableLegend()


#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
int SpectrumChart::setLeftYAxisPadding( double width, double height, Wt::WPaintDevice *paintDevice )
{
  const int oldleft = plotAreaPadding(Side::Left);
  
//  if( chartArea().width() > 5 && chartArea().height() > 5 )
//  {
//    width = chartArea().width() + plotAreaPadding(Side::Left) + plotAreaPadding(Side::Right);
//    height = chartArea().height() + plotAreaPadding(Side::Top) + plotAreaPadding(Side::Bottom);
//  }
  
  const bool inited = initLayout( WRectF(0.0, 0.0, width, height ), paintDevice);
  if( !inited )
    return -1;
  
  Chart::WAxis &yaxis = axis(Chart::Axis::Ordinate);
  
  vector<MyTickLabel> ticks;
  getYAxisLabelTicks( this, yaxis, ticks);
  
/*
  if( axis(Chart::Axis::Y1).autoLimits() || ticks.empty() )
  {
    double ymin = yaxis.minimum();
    double ymax = yaxis.maximum();

    //I dont thing this section of the code will get called anymore.
 #if( WT_VERSION < 0x3030600 )
    SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model().get());
 #else
    auto proxyModel = dynamic_cast<Wt::Chart::WStandardChartProxyModel *>( model().get());
    SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( proxyModel ? proxyModel->sourceModel().get() : nullptr);
    assert( theModel);
 #endif

    if( !theModel )
      return -1;
    
    const int nbin = theModel->rowCount();
    
    double x0 = axis(Chart::Axis::X).minimum();
    double x1 = axis(Chart::Axis::X).maximum();
    
    if( axis(Chart::Axis::X).autoLimits() )
    {      
      x0 = theModel->rowLowEdge(0);
      x1 = theModel->rowLowEdge(nbin-1) + theModel->rowWidth(nbin-1);
    }//if( axis(Chart::Axis::X).autoLimits() )
  
    yAxisRangeFromXRange( x0, x1, ymin, ymax);
    
    getYAxisLabelTicks( this, yaxis, ticks, ymin, ymax);
  }
*/
  
  size_t nchars = 2;
  for( size_t i = 0; i < ticks.size(); ++i )
    nchars = std::max( nchars, ticks[i].label.narrow().size());
  
  //A better solution would be to use a WPdfImage paint device to measure the
  //  label width...
//  WPdfImage *img = new WPdfImage(100, 100);
//  WPainter *painter = new WPainter(img);
//  painter->setFont( yaxis.labelFont());
//  const double maxWidth = -1;
//  const bool wordWrap = false;
//  WTextItem t = painter->device()->measureText("MyLabel", maxWidth, wordWrap);
//  const double REL_ERROR = 1.02;
//  double left = t.width() * REL_ERROR;
//  delete painter;
//  delete img;
  
  
  const int left = 29 + static_cast<int>(6.9*nchars);
//  cout << "Got left=" << left << ", oldleft=" << oldleft << ", nchars=" << nchars << " from " << ticks.size() << " ticks" << endl;
  

  if( left == oldleft )
    return 0;
    
  setPlotAreaPadding( left, Side::Left);
  
  return 1;
}//setLeftYAxisPadding()
#endif


void SpectrumChart::setDefaultPeakColor( const Wt::WColor &color )
{
  m_defaultPeakColor = color.isDefault() ? ns_default_peakFillColor : color;
}

void SpectrumChart::setOccupiedTimeSamplesColor( const Wt::WColor &color )
{
  m_occupiedMarkerColor = color.isDefault() ? ns_default_occupiedTimeColor : color;
  m_occupiedMarkerColor.setRgb( m_occupiedMarkerColor.red(), m_occupiedMarkerColor.green(), m_occupiedMarkerColor.blue(), 75);
}

void SpectrumChart::setForegroundHighlightColor( const Wt::WColor &color )
{
  auto &c = m_timeHighlightColors[0];
  c = color.isDefault() ? ns_defaultTimeHighlightColors[0] : color;
  c.setRgb( c.red(), c.green(), c.blue(), 155);
}

void SpectrumChart::setBackgroundHighlightColor( const Wt::WColor &color )
{
  auto &c = m_timeHighlightColors[1];
  c = color.isDefault() ? ns_defaultTimeHighlightColors[1] : color;
  c.setRgb( c.red(), c.green(), c.blue(), 75);
}

void SpectrumChart::setSecondaryHighlightColor( const Wt::WColor &color )
{
  auto &c = m_timeHighlightColors[2];
  c = color.isDefault() ? ns_defaultTimeHighlightColors[2] : color;
  c.setRgb( c.red(), c.green(), c.blue(), 75);
}//void setSecondaryHighlightColor( const Wt::WColor &color )

void SpectrumChart::setAxisLineColor( const Wt::WColor &color )
{
  WPen pen( color);
  axis(Chart::Axis::X).setPen( pen);
  axis(Chart::Axis::Y1).setPen( pen);
}//void setAxisLineColor( const Wt::WColor &color )

void SpectrumChart::setChartMarginColor( const Wt::WColor &color )
{
  if( color.isDefault() )
    m_chartMarginBrush = WBrush();
  else
    m_chartMarginBrush = WBrush( color);
}//void setChartMarginColor( const Wt::WColor &color )

void SpectrumChart::setChartBackgroundColor( const Wt::WColor &color )
{
  if( color.isDefault() )
    setBackground( WBrush());
  else
    setBackground( WBrush(color));
}//void setChartBackgroundColor( const Wt::WColor &color )

void SpectrumChart::setTextColor( const Wt::WColor &color )
{
  if( color.isDefault() )
    m_textPen = WPen( WColor(StandardColor::Black));
  else
    m_textPen = WPen( color);
  setTextPen( m_textPen);
}//void setTextColor( const Wt::WColor &color )


void SpectrumChart::saveChartToPng( const std::string &filename )
{
  LOAD_JAVASCRIPT(wApp, "src/js_inline/SpectrumChart.js", "SpectrumChart", wtjsCanvasToPngDownload);
  
  doJavaScript( "Wt.WT.CanvasToPngDownload('" + id() + "','" + filename +"');" );
}//void saveChartToPng( const std::string &name )


void SpectrumChart::setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )
{
  m_referencePhotoPeakLines = nuc;
  update( Wt::WFlags<PaintFlag>{});
}//void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )


void SpectrumChart::persistCurrentReferncePhotoPeakLines()
{
  if( m_referencePhotoPeakLines.m_validity != ReferenceLineInfo::InputValidity::Valid )
    return;
  
  const string current_input = m_referencePhotoPeakLines.m_input.m_input_txt;
  
  ReferenceLineInfo *prev_line = nullptr;
  for( size_t i = 0; !prev_line && (i < m_persistedPhotoPeakLines.size()); ++i )
  {
    if( m_persistedPhotoPeakLines[i].m_input.m_input_txt == current_input )
      prev_line = &(m_persistedPhotoPeakLines[i]);
  }
  
  if( prev_line )
    *prev_line = m_referencePhotoPeakLines;
  else if( m_referencePhotoPeakLines.m_validity == ReferenceLineInfo::InputValidity::Valid )
    m_persistedPhotoPeakLines.push_back( m_referencePhotoPeakLines);
  
  m_referencePhotoPeakLines.reset();
  
  update( Wt::WFlags<PaintFlag>{});
}//void persistCurrentReferncePhotoPeakLines()


void SpectrumChart::clearAllReferncePhotoPeakLines()
{
  m_referencePhotoPeakLines.reset();
  m_persistedPhotoPeakLines.clear();
  
  update( Wt::WFlags<PaintFlag>{});
}//void clearAllReferncePhotoPeakLines()


void SpectrumChart::renderReferncePhotoPeakLines( Wt::WPainter &painter,
                            const ReferenceLineInfo &nuc ) const
{
  // A port of SpectrumChartD3.js's `drawRefGammaLines()`, with lines at the same energy (to 10 eV)
  //  combined, as `ReferenceLineInfo::toJson(...)` does for the JS.
  const double minx = axis(Chart::Axis::X).minimum();
  const double maxx = axis(Chart::Axis::X).maximum();
  const double miny = axis(Chart::Axis::Y1).minimum();

  struct DrawnLine
  {
    double energy = 0.0, h = 0.0, hu = 0.0;  // hu is the height without cascade summing
    std::set<std::string> particles;
    bool major = false;
    WColor color;
  };//struct DrawnLine

  vector<const ReferenceLineInfo::RefLine *> sorted;
  for( const ReferenceLineInfo::RefLine &line : nuc.m_ref_lines )
  {
    if( line.m_normalized_intensity > std::numeric_limits<float>::min() )
      sorted.push_back( &line );
  }
  std::stable_sort( begin(sorted), end(sorted), []( const auto *lhs, const auto *rhs ){
    return lhs->m_energy < rhs->m_energy;
  } );

  const auto round_energy = []( const double e ) -> double { return std::round(100.0*e)/100.0; };

  vector<DrawnLine> lines;
  for( const ReferenceLineInfo::RefLine *line : sorted )
  {
    const double energy = round_energy( line->m_energy );
    if( (energy < minx) || (energy > maxx) )
      continue;

    if( lines.empty() || (lines.back().energy != energy) )
    {
      lines.emplace_back();
      lines.back().energy = energy;
      lines.back().color = line->m_color;
    }

    DrawnLine &drawn = lines.back();
    drawn.h += line->m_normalized_intensity;
    drawn.hu += line->m_uncorrected_intensity;
    drawn.particles.insert( line->particlestr() );
    drawn.major = (drawn.major || line->m_major_line);
  }//for( const ReferenceLineInfo::RefLine *line : sorted )

  // The uncorrected height is only used if cascade summing changed the height
  for( DrawnLine &line : lines )
  {
    if( IsNan(line.hu) || IsInf(line.hu) || (std::fabs(line.hu - line.h) <= 0.0001*line.h) )
      line.hu = line.h;
  }

  double max_amp = 0.0;
  for( const DrawnLine &line : lines )
    max_amp = std::max( max_amp, std::max( line.h, line.hu ) );

  if( !(max_amp > 0.0) )
    return;

  const WColor source_color = nuc.m_input.m_color.isDefault() ? WColor( "#0000FF" ) : nuc.m_input.m_color;

  // Heights are linear in intensity, with the largest line in view reaching 30 px below the top
  const WRectF &area = chartArea();
  const double h = area.height();
  const double top_pad = std::min( h, ns_ref_line_top_pad );
  const auto ypx = [&]( const double amp ) -> double {
    return area.top() + std::min( h - (h - top_pad)*amp/max_amp, h - 2.0 );
  };

  const auto dashes_for = []( const std::set<std::string> &particles ) -> vector<double> {
    if( particles.size() != 1 )
      return {};
    const std::string &particle = *particles.begin();
    if( particle == "xray" )
      return { 3.0, 3.0 };
    if( particle == "beta" )
      return { 1.0, 1.0 };
    if( particle == "alpha" )
      return { 3.0, 2.0, 1.0 };
    if( particle == "positron" )
      return { 3.0, 1.0 };
    if( (particle == "electronCapture") || (particle == "cascade-sum") )
      return { 6.0, 6.0 };
    if( (particle == "S.E.") || (particle == "D.E.") )
      return { 4.0, 1.0 };
    return {};
  };//dashes_for

  // Lines are collected into a path for each color and opacity, and then drawn
  vector<tuple<WColor,double,WPainterPath>> paths;
  const auto path_for = [&paths]( const WColor &color, const double opacity ) -> WPainterPath & {
    for( auto &entry : paths )
    {
      if( same_color( std::get<0>(entry), color ) && (std::get<1>(entry) == opacity) )
        return std::get<2>(entry);
    }
    paths.emplace_back( color, opacity, WPainterPath() );
    std::get<2>(paths.back()).setOpenSubPathsEnabled( true );
    return std::get<2>(paths.back());
  };//path_for

  const double threshold = 1.0E-4 * max_amp;
  for( const DrawnLine &line : lines )
  {
    if( (std::max( line.h, line.hu ) < threshold) && !line.major )
      continue;

    const WColor color = line.color.isDefault() ? source_color : line.color;
    const double x = mapToDevice( line.energy, miny ).x();
    const bool gain = (line.h > 1.0001*line.hu);  // summing-in made the line taller
    const bool loss = (line.hu > 1.0001*line.h);  // summing-out made the line shorter

    add_dashed_vline( path_for( color, 1.0 ), x, area.bottom(), ypx( gain ? line.hu : line.h ),
                      dashes_for( line.particles ) );

    // Summing-out: a translucent segment shows the amplitude lost, up to the no-summing height.
    // Summing-in: a dashed segment shows the amplitude gained, above the no-summing height.
    if( loss )
      add_dashed_vline( path_for( color, 0.35 ), x, ypx( line.hu ), ypx( line.h ), {} );
    else if( gain )
      add_dashed_vline( path_for( color, 1.0 ), x, ypx( line.hu ), ypx( line.h ), { 3.0, 2.0 } );
  }//for( const DrawnLine &line : lines )

  for( const auto &entry : paths )
  {
    WPen pen( with_alpha( std::get<0>(entry), std::get<1>(entry) ) );
    pen.setWidth( WLength(1.0) );
    pen.setCapStyle( PenCapStyle::Flat );
    painter.strokePath( std::get<2>(entry), pen );
  }
}//void renderReferncePhotoPeakLines(...)


void SpectrumChart::renderReferncePhotoPeakLines( Wt::WPainter &painter ) const
{
  for( auto iter = m_persistedPhotoPeakLines.rbegin(); iter != m_persistedPhotoPeakLines.rend(); ++iter )
    renderReferncePhotoPeakLines( painter, *iter);
  
  if( m_referencePhotoPeakLines.m_validity == ReferenceLineInfo::InputValidity::Valid )
    renderReferncePhotoPeakLines( painter, m_referencePhotoPeakLines);
}//void renderReferncePhotoPeakLines( Wt::WPainter &painter);


void SpectrumChart::modelChanged()
{
  Wt::Chart::WCartesianChart::modelChanged();

  SpectrumDataModel * const m = spectrum_model( *this );
  if( !m && model() )
    throw runtime_error( "SpectrumChart can only be used with SpectrumDataModel models");
}//void modelChanged()



