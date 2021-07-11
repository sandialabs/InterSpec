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

#include <set>
#include <map>
#include <string>
#include <algorithm>

#include <boost/any.hpp>
#include <boost/functional/hash.hpp>

#include <Wt/WFont>
#include <Wt/Utils>
#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLineF>
#include <Wt/WPoint>
#include <Wt/WRectF>
#include <Wt/WEvent>
#include <Wt/WLength>
#include <Wt/WPainter>
#include <Wt/WIconPair>
#include <Wt/WPopupMenu>
#include <Wt/WPushButton>
#include <Wt/WCircleArea>
#include <Wt/WApplication>
#include <Wt/WPainterPath>
#include <Wt/Chart/WAxis>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WPaintDevice>
#include <Wt/WFontMetrics>
#include <Wt/Chart/WDataSeries>
#include <Wt/Chart/WCartesianChart>


#if( WT_VERSION<=0x3030100 )
#include <Wt/Chart/WChart2DRenderer>
typedef Wt::Chart::WChart2DRenderer ChartRenderHelper_t;
#else
typedef Wt::Chart::WCartesianChart ChartRenderHelper_t;
#endif

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/ReferenceLineInfo.h"

#include "InterSpec/CanvasForDragging.h"
#include "js/SpectrumChart.js"

using namespace std;
using namespace Wt;

using SpecUtils::Measurement;

namespace
{
  const WColor ns_defaultTimeHighlightColors[3] = { {255,255,0,155}, {0,255,255,75}, {0,128,0,75} };
  const WColor ns_default_peakFillColor( 0, 51, 255, 155 );
  const WColor ns_default_occupiedTimeColor( 128, 128, 128, 75 );
  
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
  
void getXAxisLabelTicks( const SpectrumChart *chart,
                         const Chart::WAxis& axis,
                         std::vector<MyTickLabel>& ticks )
{
  //This function equivalentish to WAxis::getLabelTicks(...) but makes it so
  //  the x axis labels (hopefully) always line up nicely where we want them
  //  e.g. kinda like multiple of 5, 10, 25, 50, 100, etc.
  static double EPSILON = 1E-3;
    
  //  int rc = 0;
  //  if (axis.chart()->model())
  //    rc = axis.chart()->model()->rowCount();
    
  if( axis.scale() != Chart::LinearScale )
    throw runtime_error( "getXAxisLabelTicks called for wrong type of graphs" );
    
  const double rendermin = axis.minimum();
  const double rendermax = axis.maximum();
    
  const double minxpx = chart->mapToDevice(rendermin,0).x();
  const double maxxpx = chart->mapToDevice(rendermax,0).x();
  const double widthpx = maxxpx - minxpx;
  const int nlabel = widthpx / 50;
    
  const double range = rendermax - rendermin;
  double renderInterval;// = range / 10.0;
  const double n = std::pow(10, std::floor(std::log10(range/nlabel)));
  const double msd = range / nlabel / n;
    
    
  if( IsNan(n) || IsInf(n) || nlabel<=0 || n<=0.0 )  //JIC
  {
    ticks.clear();
    return;
  }
    
  int subdashes = 0;
    
  if (msd < 1.5)
  {
    subdashes = 2;
    renderInterval = 0.5*n;
  }else if (msd < 3.3)
  {
    subdashes = 5;
    renderInterval = 0.5*n;
  }else if (msd < 7)
  {
    subdashes = 5;
    renderInterval = n;
  }else
  {
    subdashes = 10;
    renderInterval = n;
  }
  
  
  //Should set axis label format here based on range and number of labels.
  //Right now it is based on only range (see setXAxisRange(...)), and in
  //  situations where there is only one digit after the deicmal place, multiple
  //  labels in a row can be exactly the same.
  
  const double biginterval = subdashes * renderInterval;
  double startEnergy = biginterval * std::floor((rendermin + std::numeric_limits<double>::epsilon()) / biginterval);
    
  if( startEnergy < (rendermin-EPSILON*renderInterval) )
    startEnergy += biginterval;
    
  for( int i = -int(floor(startEnergy-rendermin)/renderInterval); ; ++i)
  {
    if( i > 5000 )  //JIC
      break;
      
    double v = startEnergy + renderInterval * i;
      
    if( (v - rendermax) > EPSILON * renderInterval )
      break;
      
    WString t;
    
    if( i>=0 && ((i % subdashes) == 0) )
    {
      if( renderInterval < 1.0 )
      {
        //If x-range is more than 10, but we have more than 10 labels, then labelFormat flag is
        //  "%.0f", meaining no decimals, but we need decimals or else multiple ticks will have same
        //  label - so here we will .
        //  The intervals should almost always be something like 0.1, 0.2, etc, but because I didnt
        //  go through the logic totally yet to make sure, I'll leave in a default catchall that
        //  will print to four places, and then just remove trailing zeros.
        const char *fmtflag = "%.4f";
        bool removetrail = false;
        if( renderInterval == 0.1 || renderInterval == 0.2 || renderInterval == 0.5 )
          fmtflag = "%.1f";
        else if( renderInterval == 0.25 )
          fmtflag = "%.2f";
        else
          removetrail = true;
        
        char buffer[32] = { '\0' };
        int nchar = snprintf( buffer, sizeof(buffer), fmtflag, v );
        
        //Remove trailing zeros
        if( removetrail && (--nchar > 0) && (nchar < static_cast<int>(sizeof(buffer))) )
        {
          while( (nchar > 0) && ((buffer[nchar] == '0') || (buffer[nchar] == '.')) )
            buffer[nchar--] = '\0';
        }
        
        t = buffer;
      }else
      {
        t = axis.label(v);
      }
    }//If( a major tick )
    
    ticks.push_back(MyTickLabel(v, i % subdashes == 0 ? MyTickLabel::Long : MyTickLabel::Short, t));
  }//for( intervals to draw ticks for )
}//getXAxisLabelTicks(...)
  
  
void getYAxisLabelTicks( const SpectrumChart *chart,
                         const Chart::WAxis& axis,
                         std::vector<MyTickLabel> &ticks )
{
  //This function equivalentish to WAxis::getLabelTicks(...) but makes it so
  //  the y axis labels (hopefully) always line up nicely where we want them
  //  e.g. kinda like multiple of 5, 10, 25, 50, 100, etc.
  const double EPSILON = 1.0E-3;
    
  const double renderymin = axis.minimum();
  const double renderymax = axis.maximum();
  const double minypx = chart->mapToDevice(0.0,renderymax,axis.id()).y();
  const double maxypx = chart->mapToDevice(0.0,renderymin,axis.id()).y();
  const double heightpx = maxypx - minypx;
  
  //note that maxypx actually cooresponds to the smalles y-values
  const double range = renderymax - renderymin;
    
  switch( axis.scale() )
  {
    case Chart::CategoryScale: case Chart::DateScale: case Chart::DateTimeScale:
    {
      throw runtime_error( "getYAxisLabelTicks() for SpectrumChart can only"
                          " handle linear and log scales");
      break;
    }//case: unsupported scale type
        
    case Chart::LinearScale:
    {
      //px_per_div: pixels between major and/or minor labels.
      const int px_per_div = 50;
        
      //nlabel: approx number of major + minor labels we would like to have.
      const int nlabel = heightpx / px_per_div;
        
      //renderInterval: Inverse of how many large labels to place between powers
      //  of 10 (1 is none, 0.5 is is one).
      double renderInterval;
        
      //n: approx how many labels will be used
      const double n = std::pow(10, std::floor(std::log10(range/nlabel)));
        
      //msd: approx how many sub dashes we would expect there to have to be
      //     to satisfy the spacing of labels we want with the given range.
      const double msd = range / nlabel / n;
        
        
      if( IsNan(n) || IsInf(n) || nlabel<=0 || n<=0.0 )  //JIC
      {
        ticks.clear();
        return;
      }
        
      int subdashes = 0;
        
      if( msd < 1.5 )
      {
        subdashes = 2;
        renderInterval = 0.5*n;
      }else if (msd < 3.3)
      {
        subdashes = 5;
        renderInterval = 0.5*n;
      }else if (msd < 7)
      {
        subdashes = 5;
        renderInterval = n;
      }else
      {
        subdashes = 10;
        renderInterval = n;
      }
        
      const double biginterval = subdashes * renderInterval;
      double starty = biginterval * std::floor((renderymin + std::numeric_limits<double>::epsilon()) / biginterval);
        
      if( starty < (renderymin-EPSILON*renderInterval) )
        starty += biginterval;
        
      for( int i = -int(floor(starty-renderymin)/renderInterval); ; ++i)
      {
        if( i > 500 )  //JIC
          break;
          
        const double v = starty + renderInterval * i;
          
        if( (v - renderymax) > EPSILON * renderInterval )
          break;
          
        WString t;
        if( i>=0 && ((i % subdashes) == 0) )
          t = axis.label(v);
        MyTickLabel::TickLength len = ((i % subdashes) == 0) ? MyTickLabel::Long : MyTickLabel::Short;
          
        ticks.push_back( MyTickLabel(v, len, t) );
      }//for( intervals to draw ticks for )
        
      break;
    }//case Chart::LinearScale:
        
    case Chart::LogScale:
    {
      //Get the power of 10 just below or equal to rendermin.
      int minpower = (renderymin > 0.0) ? (int)floor( log10(renderymin) ) : -1;
        
      //Get the power of 10 just above or equal to renderymax.  If renderymax
      //  is less than or equal to 0, set power to be 0.
      int maxpower = (renderymax > 0.0) ? (int)ceil( log10(renderymax) ): 0;
        
      //Adjust minpower and maxpower
      if( maxpower == minpower )
      {
        //Happens when renderymin==renderymax which is a power of 10
        ++maxpower;
        --minpower;
      }else if( maxpower > 2 && minpower < -1)
      {
        //We had a tiny value (possibly a fraction of a count), as well as a
        //  large value (>1000).
        minpower = -1;
      }else if( maxpower >= 0 && minpower < -1 && (maxpower-minpower) > 6 )
      {
        //we had a tiny power (1.0E-5), as well as one between 1 and 999,
        //  so we will only show the most significant decades
        minpower = maxpower - 5;
      }//if( minpower == maxpower ) / else / else
        
        
      //numdecades: number of decades the data covers, including the decade
      //  above and below the data.
      const int numdecades = maxpower - minpower + 1;
        
      //minpxdecade: minimum number of pixels we need per decade.
      const int minpxdecade = 25;
        
      //labeldelta: the number of decades between successive labeled large ticks
      //  each decade will have a large tick regardless
      int labeldelta = 1;
        
      //nticksperdecade: number of small+large ticks per decade
      int nticksperdecade = 10;
        
        
      if( (heightpx / minpxdecade) < numdecades )
      {
        labeldelta = numdecades / (heightpx / minpxdecade);
        nticksperdecade = 1;
//        if( labeldelta < 3 )
//          nticksperdecade = 10;
//        else if( labeldelta < 4 )
//          nticksperdecade = 5;
//        else
//          nticksperdecade = 3;
      }//if( (heightpx / minpxdecade) < numdecades )
        
        
      int nticks = 0;
      int nmajorticks = 0;
        
      for( int decade = minpower; decade <= maxpower; ++decade )
      {
        const double startcounts = std::pow( 10.0, decade );
        const double deltacounts = 10.0 * startcounts / nticksperdecade;
        const double eps = deltacounts * EPSILON;
          
        if( (startcounts - renderymin) > -eps && (startcounts - renderymax) < eps )
        {
          const WString t = ((decade%labeldelta)==0) ? axis.label(startcounts)
          : WString("");
          ++nticks;
          ++nmajorticks;
          ticks.push_back( MyTickLabel(startcounts, MyTickLabel::Long, t) );
        }//if( startcounts >= renderymin && startcounts <= renderymax )
          
          
        for( int i = 1; i < (nticksperdecade-1); ++i )
        {
          const double y = startcounts + i*deltacounts;
          const WString t = axis.label(y);
            
          if( (y - renderymin) > -eps && (y - renderymax) < eps )
          {
            ++nticks;
            ticks.push_back( MyTickLabel(y, MyTickLabel::Short, t) );
          }
        }//for( int i = 1; i < nticksperdecade; ++i )
      }//for( int decade = minpower; decade <= maxpower; ++decade )
        
      //If we have a decent number of (sub) labels, the user can orient
      //  themselves okay, so well get rid of the minor labels.
      if( (nticks > 8 || (heightpx/nticks) < 25 || nmajorticks > 1) && nmajorticks > 0 )
      {
        for( size_t i = 0; i < ticks.size(); ++i )
          if( ticks[i].tickLength == MyTickLabel::Short )
            ticks[i].label = "";
      }
        
      if( ticks.empty() )
      {
        cerr << "Forcing a single axis point in" << endl;
        const double y = 0.5*(renderymin+renderymax);
        const WString t = axis.label(y);
        ticks.push_back( MyTickLabel(y, MyTickLabel::Long, t) );
      }
    }//case Chart::LogScale:
  }//switch( axis.scale() )
    
}//getYAxisLabelTicks(...)
  
  
/** Places all peaks that share a ROI into contiguous order, where the contigous
   peaks are sorted by the mean of their first peak, and within the ROI the
   peaks are sorted as well.
 */
void prepare_gaus_peaks_for_iterating( std::vector< PeakModel::PeakShrdPtr > &gauspeaks )
{
  if( gauspeaks.size() < 2 )
    return;
  
  std::sort( begin(gauspeaks), end(gauspeaks), &PeakDef::lessThanByMeanShrdPtr );
  
  //We will almost never need to re-group the peaks (only if there is a peak
  //  in a ROI energy range, that doesnt share the ROI with the other peaks - a
  //  very rare thing to happen), so we will do a check that is hopefully a
  //  little faster than doing the grouping to detect this condition.
  bool need_to_group = false;
  set<const PeakContinuum *> roiptrsset;
  const PeakContinuum *last_cont = nullptr;
  for( size_t i = 0; !need_to_group && (i < gauspeaks.size()); ++i )
  {
    const PeakContinuum * const cont = gauspeaks[i]->continuum().get();
    if( cont != last_cont )
    {
      need_to_group = roiptrsset.count( cont );
      roiptrsset.insert( cont );
      last_cont = cont;
    }
  }
  if( !need_to_group )
    return;
  
  //We need to do the grouping.
  std::vector< std::pair<const PeakContinuum *,std::vector<uint16_t> > > roi_to_index;
  for( size_t i = 0; i < gauspeaks.size(); ++i )
  {
    const uint16_t peak_index = static_cast<uint16_t>(i);
    const PeakContinuum * const cont = gauspeaks[i]->continuum().get();
    
    //Find index to this PeakContinuum in the vector
    const size_t num_index = roi_to_index.size();
    size_t index = roi_to_index.size();
    for( size_t i = num_index; (index==num_index) && (i > 0); --i )
    {
      if( roi_to_index[i-1].first == cont )
        index = i-1;
    }
    
    //If we didnt find an index, add an entry into the vector
    if( index == num_index )
    {
      roi_to_index.resize(num_index + 1);
      roi_to_index.back().first = cont;
    }
    
    roi_to_index[index].second.push_back( peak_index );
  }//for( size_t i = 0; i < gauspeaks.size(); ++i )
  
  /*
   //If we hadnt already detected if we needed needed to do the grouping, we
   //  could with this block of code (leaving in to test against a little more).
  bool all_sequential = true;
  for( const auto &roi_p : roi_to_index )
  {
    for( size_t i = 1; all_sequential && (i < roi_p.second.size()); ++i )
      all_sequential = ((roi_p.second[i-1] + 1) == roi_p.second[i]);
    if( !all_sequential )
      break;
  }//for( const auto &roi_p : roi_to_index )
  
  if( all_sequential )
    return;
   */
  
  std::vector< PeakModel::PeakShrdPtr > newpeaks;
  newpeaks.reserve( gauspeaks.size() );
  for( const auto &roi_pind : roi_to_index )
  {
    for( auto ind : roi_pind.second )
      newpeaks.push_back( gauspeaks[ind] );
  }
  
  assert( newpeaks.size() == gauspeaks.size() );
  
  newpeaks.swap( gauspeaks );
}//void prepare_gaus_peaks_for_iterating(...)
  
  
std::vector<std::shared_ptr<const PeakDef> >::const_iterator next_gaus_peaks(
   std::vector<std::shared_ptr<const PeakDef> >::const_iterator peak_start,
   std::vector<std::shared_ptr<const PeakDef> >::const_iterator peak_end,
   std::vector<std::shared_ptr<const PeakDef> > &peaks,
   std::shared_ptr<const Measurement> data )
{
  peaks.clear();
  
  if( peak_start == peak_end )
    return peak_end;
  
  std::shared_ptr<const PeakDef> first_peak = *peak_start;

  peaks.push_back( first_peak );
  
  for( auto iter = peak_start + 1; iter != peak_end; ++iter )
  {
    std::shared_ptr<const PeakDef> peak = *iter;
    if( peak->continuum() == first_peak->continuum() )
    {
      peaks.push_back( peak );
    }else
    {
      return iter;
    }
  }//for( iter = peak_start; iter != peak_end; ++iter )
  
  return peak_end;
}//next_gaus_peaks(...)
  
}//namespace


#if( WT_VERSION<=0x3030100 )
class SpectrumRenderer : public Chart::WChart2DRenderer
{
public:
  SpectrumRenderer( SpectrumChart *chart,
                    WPainter &painter, const WRectF &rectangle )
    : Chart::WChart2DRenderer( chart, painter, rectangle ),
      m_rect( rectangle ),
      m_chart( chart ),
      m_model( dynamic_cast<SpectrumDataModel *>( chart->model() ))
  {
    if( !m_model )
      throw runtime_error( "SpectrumRenderer: you must pass in a SpectrumChart"
                           " with a valid SpectrumDataModel" );
  }
  
  virtual ~SpectrumRenderer()
  {
  }

  
protected:
  virtual void render()
  {
    m_chart->customRender( painter(), m_rect );
  }
  

  WRectF m_rect;
  SpectrumChart *m_chart;
  SpectrumDataModel *m_model;
};//class SpectrumRenderer

#endif

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
  SpectrumRenderIterator( const SpectrumChart &chart, Wt::WPainter &painter );
  virtual void startSegment(int currentXSegment, int currentYSegment,
                            const WRectF& currentSegmentArea);
  virtual void endSegment();
  virtual bool startSeries( const Chart::WDataSeries& series, double groupWidth,
                           int , int );
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
                           it_.currentXSegment(), it_.currentYSegment() );
    
    if( curveLength_ == 0 )
    {
      curve_.moveTo(hv(p));
      
      if (series_.fillRange() != Chart::NoFill
          && series_.brush() != NoBrush) {
        fill_.moveTo(hv(fillOtherPoint(x)));
        fill_.lineTo(hv(p));
      }
    }else
    {
      if (series_.type() == Chart::LineSeries) {
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
      if( series_.type() == Wt::Chart::CurveSeries )
      {
        WPointF c1;
        computeC(p0, p_1, c1);
        curve_.cubicTo(hv(c_), hv(c1), hv(p0));
        fill_.cubicTo(hv(c_), hv(c1), hv(p0));
      }
      
      if (series_.fillRange() != Chart::NoFill
          && series_.brush() != NoBrush) {
        fill_.lineTo(hv(fillOtherPoint(lastX_)));
        fill_.closeSubPath();
        painter_.setShadow(series_.shadow());
        painter_.fillPath(fill_, series_.brush());
      }
      
      if (series_.fillRange() == Chart::NoFill)
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
    return std::sqrt( dx*dx + dy*dy );
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
      case Chart::MinimumValueFill:
        return WPointF(chart_.map(x, 0, series_.axis(),
                                  it_.currentXSegment(),
                                  it_.currentYSegment()).x(),
                       chart_.chartArea().bottom());
      case Chart::MaximumValueFill:
        return WPointF(chart_.map(x, 0, series_.axis(),
                                  it_.currentXSegment(),
                                  it_.currentYSegment()).x(),
                       chart_.chartArea().top());
      case Chart::ZeroValueFill:
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
  seriesRenderer_ = new SpectrumSeriesRenderer( chart_, painter_, series, *this );
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
    seriesRenderer_->addValue( x, y, stackY, xIndex, yIndex );
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
    w = static_cast<int>( rectangle.width() );
    h = static_cast<int>( rectangle.height() );
  }
  
  const int padLeft = plotAreaPadding(Left);
  const int padRight = plotAreaPadding(Right);
  const int padTop = plotAreaPadding(Top);
  const int padBottom = plotAreaPadding(Bottom);
  
  if( orientation() == Vertical )
    m_renderRectangle = WRectF( padLeft, padTop,
                  std::max(10, w - padLeft - padRight),
                  std::max(10, h - padTop - padBottom) );
  else
    m_renderRectangle = WRectF( padTop, padRight,
                  std::max(10, w - padTop - padBottom),
                  std::max(10, h - padRight - padLeft) );
  
}//void calcRenderRectangle( const Wt::WRectF &rectangle )


#if(WT_VERSION<=0x3030100)
Wt::WPointF SpectrumChart::map( double xValue, double yValue,
                Wt::Chart::Axis yAxis,
                int XSegment,
                int YSegment ) const
{
  return mapToDevice( xValue, yValue, yAxis, XSegment, YSegment );
}//map

Wt::WPointF SpectrumChart::hv( double x, double y ) const
{
  return hv( WPointF(x,y) );
}

Wt::WPointF SpectrumChart::hv( const WPointF &r ) const
{
  if( orientation() == Vertical )
    return r;
  throw runtime_error( "SpectrumChart does not support non-Vertical orientation for Wt 3.3.1 or earlier" );
}

Wt::WRectF SpectrumChart::hv( const Wt::WRectF &p ) const
{
  if( orientation() == Vertical )
    return p;
  throw runtime_error( "SpectrumChart does not support non-Vertical orientation for Wt 3.3.1 or earlier" );
}


void SpectrumChart::renderLabel(WPainter& painter, const WString& text,
                                 const WPointF& p, const WColor& color,
                                 WFlags<AlignmentFlag> flags,
                                 double angle, int margin) const
{
  AlignmentFlag horizontalAlign = flags & AlignHorizontalMask;
  AlignmentFlag verticalAlign = flags & AlignVerticalMask;
  
  AlignmentFlag rHorizontalAlign = horizontalAlign;
  AlignmentFlag rVerticalAlign = verticalAlign;
  
  double width = 1000;
  double height = 20;
  
  WPointF pos = hv(p);
  
  if (orientation() == Horizontal) {
    switch (horizontalAlign) {
      case AlignLeft:
        rVerticalAlign = AlignTop; break;
      case AlignCenter:
        rVerticalAlign = AlignMiddle; break;
      case AlignRight:
        rVerticalAlign = AlignBottom; break;
      default:
        break;
    }
    
    switch (verticalAlign) {
      case AlignTop:
        rHorizontalAlign = AlignRight; break;
      case AlignMiddle:
        rHorizontalAlign = AlignCenter; break;
      case AlignBottom:
        rHorizontalAlign = AlignLeft; break;
      default:
        break;
    }
  }
  
  double left = pos.x();
  double top = pos.y();
  
  switch (rHorizontalAlign) {
    case AlignLeft:
      left += margin; break;
    case AlignCenter:
      left -= width/2; break;
    case AlignRight:
      left -= width + margin;
    default:
      break;
  }
  
  switch (rVerticalAlign) {
    case AlignTop:
      top += margin; break;
    case AlignMiddle:
      top -= height/2; break;
    case AlignBottom:
      top -= height + margin; break;
    default:
      break;
  }
  
  WPen pen(color);
  WPen oldPen = painter.pen();
  painter.setPen(pen);
  
  if (angle == 0)
    painter.drawText(WRectF(left, top, width, height),
                     rHorizontalAlign | rVerticalAlign, text);
  else {
    painter.save();
    painter.translate(pos);
    painter.rotate(-angle);
    painter.drawText(WRectF(left - pos.x(), top - pos.y(), width, height),
                     rHorizontalAlign | rVerticalAlign, text);
    painter.restore();
  }
  
  painter.setPen(oldPen);
}//#if(WT_VERSION<=0x3030100)



Wt::WRectF SpectrumChart::chartSegmentArea( Wt::Chart::WAxis yAxis, int xSegment,
                                              int ySegment) const
{
  cerr << "SpectrumChart::chartSegmentArea(...) for Wt<=3.3.1 isnt quite correct" << endl;
  return m_renderRectangle;
  
//  const WAxis& xAxis = axis(Chart::XAxis);
//  
//  const Chart::WAxis::Segment& xs = xAxis.segments_[xSegment];
//  const Chart::WAxis::Segment& ys = yAxis.segments_[ySegment];
//  
//  // margin used when clipping, see also WAxis::prepareRender(),
//  // when the renderMinimum/maximum is 0, clipping is done exact
//  
//  double x1 = xs.renderStart
//  + (xSegment == 0
//     ? (xs.renderMinimum == 0 ? 0 : -chart_->axisPadding())
//     : -segmentMargin_/2);
//  double x2 = xs.renderStart + xs.renderLength
//  + (xSegment == xAxis.segmentCount() - 1
//     ? (xs.renderMaximum == 0 ? 0 : chart_->axisPadding())
//     : segmentMargin_/2);
//  
//  double y1 = ys.renderStart - ys.renderLength
//  - (ySegment == yAxis.segmentCount() - 1
//     ? (ys.renderMaximum == 0 ? 0 : chart_->axisPadding())
//     : segmentMargin_/2);
//  double y2 = ys.renderStart
//  + (ySegment == 0
//     ? (ys.renderMinimum == 0 ? 0 : chart_->axisPadding())
//     : segmentMargin_/2);
//  
//  return WRectF(std::floor(x1 + 0.5), std::floor(y1 + 0.5),
//                std::floor(x2 - x1 + 0.5), std::floor(y2 - y1 + 0.5));
}
#endif


void SpectrumChart::labelRender( Wt::WPainter& painter, const Wt::WString& text,
                                const Wt::WPointF& p, const Wt::WColor& color,
                                Wt::WFlags<Wt::AlignmentFlag> flags,
                                double angle, int margin) const
{
#if(WT_VERSION>=0x3030400)
  const WPen oldpen = painter.pen();
  painter.setPen( WPen(color) );
  renderLabel( painter, text, p, /*color,*/ flags, angle, margin );
  painter.setPen( oldpen );
#else
  renderLabel( painter, text, p, color, flags, angle, margin );
#endif
}//void labelRender(...)


#if(WT_VERSION<=0x3030100)
void SpectrumChart::customRender(	Wt::WPainter &painter, const Wt::WRectF &rectangle )
#else
void SpectrumChart::render(	Wt::WPainter &painter, const Wt::WRectF &rectangle ) const
#endif
{
  calcRenderRectangle( rectangle );
  
#if(WT_VERSION<=0x3030100)
  initLayout(rectangle);
#else
  if( initLayout(rectangle) )
#endif
  {
    renderChartBackground(painter, rectangle);
    renderGridLines(painter, Chart::XAxis);
    renderGridLines(painter, Chart::Y1Axis);
    renderGridLines(painter, Chart::Y2Axis);
    renderSeries(painter);
    
//    if( axis(Chart::Y2Axis).isVisible() )
//      axis(Chart::Y2Axis).setTitle( "Neutron CPS" ); //
    
    renderAxes(painter, Chart::Line | Chart::Labels);
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
      painter.fillRect( hv(WRectF(rx,0,rw,ry)), m_chartMarginBrush );
    if( (ry+rh) < rectangle.height() ) //Bottom strip
      painter.fillRect( hv(WRectF(rx,ry+rh,rw,rectangle.height()-ry)), m_chartMarginBrush );
    if( rx > 0 )
      painter.fillRect( hv(WRectF(0,0,rx,rectangle.height())), m_chartMarginBrush );
    if( (rx+rw) < rectangle.width() )
      painter.fillRect( hv(WRectF(rx+rw,0,rectangle.width()-rw-rx,rectangle.height())), m_chartMarginBrush );
  };
  
  if( background().style() != NoBrush && m_chartMarginBrush.style() != NoBrush )
  {
    paintMargins();
    painter.fillRect( hv(chartArea()), background() );
  }else if( background().style() != NoBrush )
  {
    painter.fillRect( hv(rectangle), background() );
  }else if( m_chartMarginBrush.style() != NoBrush )
  {
    paintMargins();
  }
  
  SpectrumDataModel *specModel = dynamic_cast<SpectrumDataModel *>( model() );
  if( textInMiddleOfChart().empty()
     && specModel && !specModel->histUsedForXAxis() )
  {
    const double renderxmin = axis(Chart::XAxis).minimum();
    const double renderxmax = axis(Chart::XAxis).maximum();
    const double renderymin = axis(Chart::YAxis).minimum();
    const double renderymax = axis(Chart::YAxis).maximum();
    const WPointF upperLeft = mapToDevice(renderxmin,renderymax);
    const WPointF lowerRight = mapToDevice(renderxmax,renderymin);
    
    const double minxpx = upperLeft.x();
    const double maxxpx = lowerRight.x();
    const double minypx = upperLeft.y();
    const double maxypx = lowerRight.y();
    const double xrange = maxxpx - minxpx;
    const double yrange = maxypx - minypx;
    
    //ToDo: Get SVG version of logo
    WPainter::Image snllogo( "InterSpec_resources/images/SNL_Stacked_Black_Blue.png", 688, 265 );  //17 kb
    
    //Lets draw the Sandia logo to take up about 1/4th the height of the screen
    const double logo_scale = 0.25 * yrange / snllogo.height();
    const double logo_width = logo_scale * snllogo.width();
    const double logo_height = logo_scale * snllogo.height();
    
    const double middle_px_x = 0.5*(minxpx + maxxpx);
    const double middle_px_y = 0.5*(minypx + maxypx);
    
    //Try to Make "InterSpec" be about the same width as the SNL logo
    const int font_size = static_cast<int>( logo_width / 4 );
    
    WRectF snlrect( middle_px_x - 0.5*logo_width, middle_px_y - 0.5*logo_height - font_size, logo_width, logo_height );
    painter.drawImage( snlrect, snllogo );
    
    painter.save();
    WFont font( WFont::Monospace );
    font.setSize( font_size );
    font.setWeight( WFont::Weight::NormalWeight );
    font.setVariant( WFont::SmallCaps );
    painter.setFont( font );
    painter.setPen( m_textPen );
    double texty = middle_px_y + 0.5*logo_height - 1.05*font_size;
    painter.drawText( minxpx, texty, xrange, 1.05*font_size, AlignCenter, "InterSpec" );
    painter.restore();
  }//if( specModel && !specModel->histUsedForXAxis() )
  
  painter.setFont( oldFont );
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
    case Wt::Chart::XAxis:
      vertical = true;
      getXAxisLabelTicks( this, axis, ticks );
    break;
    
    case Wt::Chart::YAxis: case Wt::Chart::Y2Axis:
//    case Wt::Chart::Y1Axis: case Wt::Chart::OrdinateAxis:
      vertical = false;
      getYAxisLabelTicks( this, axis, ticks );
    break;
      
    default:
      return;
      break;
  }//switch( axistype )
  
  const Chart::WAxis &yaxis = WCartesianChart::axis(Chart::Y1Axis);
  const Chart::WAxis &xaxis = WCartesianChart::axis(Chart::XAxis);
  
  const double ymin = yaxis.minimum();
  const double ymax = yaxis.maximum();
  const double xmin = xaxis.minimum();
  const double xmax = xaxis.maximum();
  
  
  WPainterPath majorGridPath, minorGridPath;
  
  for( size_t i = 0; i < ticks.size(); ++i )
  {
    const double d = ticks[i].u;
    
    const double x0 = vertical ? d : xmin;
    const double x1 = vertical ? d : xmax;
    const double y0 = vertical ? ymin : d;
    const double y1 = vertical ? ymax : d;
    
    if( ticks[i].tickLength == MyTickLabel::Long )
    {
      majorGridPath.moveTo( hv(mapToDevice(x0,y0)) );
      majorGridPath.lineTo( hv(mapToDevice(x1,y1)) );
    }else if( ticks[i].tickLength == MyTickLabel::Short )
    {
      minorGridPath.moveTo( hv(mapToDevice(x0,y0)) );
      minorGridPath.lineTo( hv(mapToDevice(x1,y1)) );
    }
  }//for( size_t i = 0; i < ticks.size(); ++i )
  
  if( !majorGridPath.isEmpty() )
    painter.strokePath(majorGridPath, axis.gridLinesPen());
  
  WPen minorpen = axis.gridLinesPen();
  
  WColor minorColor = minorpen.color();
  if( minorColor.name().empty() )
  {
    minorColor.setRgb( minorColor.red(), minorColor.green(), minorColor.blue(), 40 );
    minorpen.setColor( minorColor );
  }//if( minorColor.name().empty() )
  
  minorpen.setWidth( 1 );
  
  if( !minorGridPath.isEmpty() )
    painter.strokePath(minorGridPath, minorpen );
}//void renderXGrid( Wt::WPainter &painter, const Wt::WAxis &axis ) const



void SpectrumChart::renderSeries( Wt::WPainter &painter ) const
{
  if( !isLargeEnough( painter ) )
  {
    const double height = static_cast<double>( painter.viewPort().height() );
    const double width = static_cast<double>( painter.viewPort().width() );
    painter.drawText( 0.0, 0.5*height - 8.0, width, 20.0,
                       AlignCenter, "Enlarge to view chart" );
    return;
  }//if( !isLargeEnough() )
  
  {
    SpectrumRenderIterator iterator( *this, painter );
    iterateSpectrum( &iterator, painter );
  }
  
  //Label and marker iterator code removed 20150224 since they were causing a
  // slowdown (the marker iterator was accessing all rows of the data model
  // no matter the x-range), and we arent currently using them anyway.
//  {
//    LabelRenderIterator iterator(*this, painter);
//    iterateSeries( &iterator, &painter );
//  }
//  {
//    MarkerRenderIterator iterator(*this, painter);
//    iterateSeries( &iterator, &painter );
//  }
}//void renderSeries( Wt::WPainter &painter ) const;


void SpectrumChart::renderAxes( Wt::WPainter &painter,
                      Wt::WFlags<  Wt::Chart::AxisProperty > properties ) const
{
  if( !isLargeEnough(painter) )
    return;
  
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  
  if( m )
  {
    renderXAxis( painter, axis(Chart::XAxis), properties);
    renderYAxis( painter, axis(Chart::Y1Axis), properties );
    renderYAxis( painter, axis(Chart::Y2Axis), properties );
  }else
  {
#if(WT_VERSION<=0x3030100)
    throw runtime_error( "SpectrumChart only supports SpectrumDataModel for Wt 3.3.1 or before" );
#else
    renderAxis( painter, axis(Chart::XAxis), properties);
    renderAxis( painter, axis(Chart::Y1Axis), properties );
    renderAxis( painter, axis(Chart::Y2Axis), properties );
#endif
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
  AlignmentFlag labelHFlag = AlignCenter;
  
  axisStart.setY(chartArea().bottom() + 0.5);
  axisEnd.setY(chartArea().top() + 0.5);
  
  enum { Left = 0x1, Right = 0x2, Both = 0x3 } tickPos = Left;
  
  const int TICK_LENGTH = 5;
  Chart::AxisValue location = yaxis.location();
  if( yaxis.id() == Chart::Y2Axis )
    location = Chart::MaximumValue;
  
  switch( location )
  {
    case Chart::MinimumValue:
    {
      const double x = chartArea().left() - yaxis.margin() + 0.5;
      axisStart.setX( x );
      axisEnd.setX( x );
      labelHFlag = AlignRight;
      tickPos = Left;
      
      break;
    }//case Chart::MinimumValue:
      
    case Chart::MaximumValue:
    {
      double x = chartArea().right() + yaxis.margin() + 0.5;
      axisStart.setX(x);
      axisEnd.setX(x);
      labelHFlag = AlignLeft;
      tickPos = Right;
      
      break;
    }//case Chart::MaximumValue:
    
    case Chart::ZeroValue:
    {
      double x = std::floor(map(0, 0, Chart::YAxis).x()) + 0.5;
      axisStart.setX(x);
      axisEnd.setX(x);
      labelHFlag = AlignRight;
      tickPos = Both;
      
      break;
    }//case Chart::ZeroValue:
  }//switch( location )
  
  
  //Render Title
  if( (properties & Chart::Labels) && !yaxis.title().empty() )
  {
/*
     WFont oldFont2 = painter.font();
     painter.setFont( axis.titleFont() );
     
     const double u = axisStart.x();
     const double x = u + (labelHFlag == AlignRight ? 15 : -15);
     const double y = chartArea().top() - 8;
     const WFlags<AlignmentFlag> flag = labelHFlag | AlignBottom;
     labelRender( painter, axis.title(), WPointF(x,y), black, flag, 0, 0);
     painter.setFont(oldFont2);
*/
    
    //XXX - wont work for Y2Axis
    const WPen oldPen = painter.pen();
    const WFont oldFont = painter.font();
    const WBrush oldBrush = painter.brush();
    
    painter.save();
    painter.rotate( -90 );
    painter.setPen( m_textPen ); //axis(Chart::YAxis).titleFont()
    painter.setBrush( WBrush() );
    painter.setFont( yaxis.titleFont() );

    const double h = axisStart.y() - axisEnd.y();
    
    const WString ytitle = axis(Chart::YAxis).title();
    if( !ytitle.empty() )
      painter.drawText( -axisStart.y(), 1.0, h, 20,
                        AlignCenter | AlignTop, ytitle );
    painter.restore();
    painter.setPen( oldPen );
    painter.setFont( oldFont );
    painter.setBrush( oldBrush );
  }//if( draw labels and there is a title )
  
  
  //If paint the chart line
  if( properties & Chart::Line )
  {
    const WPen oldPen = painter.pen();
    painter.setPen( yaxis.pen() );
    painter.drawLine( hv(axisStart), hv(axisEnd) );
    painter.setPen( oldPen );
  }//if( properties & Chart::Line )
  
  
  //Draw the labels
  std::vector<MyTickLabel> ticks;
  getYAxisLabelTicks( this, yaxis, ticks );
  const size_t nticks = std::min( int(ticks.size()), 1000 );
  
  WPainterPath ticksPath;
  
  const WFont oldFont = painter.font();
  const WPen oldPen = painter.pen();
  
  painter.setFont( yaxis.labelFont() );
  
  for( size_t i = 0; i < nticks; ++i )
  {
    const double yval = ticks[i].u;
    const double ypx = std::floor(mapToDevice(1.0,yval,yaxis.id()).y()) + 0.5;
    
    const int tickLength = ticks[i].tickLength == MyTickLabel::Long
                            ? TICK_LENGTH : TICK_LENGTH / 2;
    
    if( ticks[i].tickLength != MyTickLabel::Zero )
    {
      double x = axisStart.x() + (tickPos & Right ? +tickLength : 0);
      const double y = hv(x,ypx).y();
      ticksPath.moveTo( x, y );
      x = axisStart.x() + (tickPos & Left ? -tickLength : 0);
      ticksPath.lineTo( x, y );
    }//if( non zero tick length )
    
    const int txtwidth = 200;
    const int txtheight = 20;
    const int txtPadding = 3;
    
    double labelxpx, labelypx;
    WFlags<AlignmentFlag> labelAlignFlags;
    
    switch( location )
    {
      case Wt::Chart::MinimumValue:
      case Wt::Chart::ZeroValue:
        labelxpx = axisStart.x() + (tickPos & Left ? -tickLength : 0) - txtwidth - txtPadding;
        labelypx = ypx - txtheight/2;
        labelAlignFlags = AlignRight | AlignMiddle;
      break;
        
      case Wt::Chart::MaximumValue:
        labelxpx = axisStart.x() + (tickPos & Right ? +tickLength : 0);
        labelypx = ypx - txtheight/2;
        labelAlignFlags = AlignLeft | AlignMiddle;
      break;
    }//switch( location )
    
    if( (properties & Chart::Labels) && !ticks[i].label.empty()
       && (!m_avoidMobileMenu || labelypx > 42) ) //Check that we dont overlap with the hamburger menu.
    {   
      painter.setPen( m_textPen );
      painter.drawText( WRectF(labelxpx, labelypx, txtwidth, txtheight),
                        labelAlignFlags, ticks[i].label);
      painter.setPen( oldPen );
    }
  }//for( size_t i = 0; i < nticks; ++i )
  
  painter.setFont( oldFont );
  
  if ((properties & Chart::Line) && !ticksPath.isEmpty() )
    painter.strokePath( ticksPath, yaxis.pen() );
}//void renderYAxis(...)



void SpectrumChart::renderXAxis( Wt::WPainter &painter,
                                const Wt::Chart::WAxis& axis,
                                Wt::WFlags<Wt::Chart::AxisProperty> properties ) const
{
  if( !axis.isVisible() )
    return;
  
  //This function adapted from WChart2DRenderer::renderAxis(...) in order to
  //  make x axis labels always line up nice
  WFont oldFont1 = painter.font();
  WFont labelFont = axis.labelFont();
  painter.setFont(labelFont);
  
  double lineypx = 0;
  enum { Left = 0x1, Right = 0x2, Both = 0x3 } tickPos = Left;
  AlignmentFlag labelHFlag = AlignLeft;
  
  Chart::AxisValue location = axis.location();
  
  switch( location )
  {
    case Chart::MinimumValue:
      tickPos = Left;
      labelHFlag = AlignTop;
      lineypx = chartArea().bottom() + 0.5 + axis.margin();
      break;
      
    case Chart::MaximumValue:
      tickPos = Right;
      labelHFlag = AlignBottom;
      lineypx = chartArea().top() - 0.5 - axis.margin();
      break;
      
    case Chart::ZeroValue:
      tickPos = Both;
      labelHFlag = AlignTop;
      lineypx = std::floor(map(0, 0, Chart::YAxis).y()) + 0.5;
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
    if( properties & Chart::Line )
    {
      painter.setPen( axis.pen() );
      const WPointF begin = hv( minxpx, lineypx );
      const WPointF end = hv( maxxpx, lineypx );
      painter.drawLine(begin, end);
    }
    
    const double titlewidthpx = 10.0*axis.title().narrow().size();
    
    WPainterPath ticksPath;
    
    std::vector<MyTickLabel> ticks;
    getXAxisLabelTicks( this, axis, ticks );
    
    const int TICK_LENGTH = 5;
    
    const size_t nticks = std::min( int(ticks.size()), 1000 );
    
    for( size_t i = 0; i < nticks; ++i )
    {
      const double d = ticks[i].u;
      const double dd = std::floor(mapToDevice(d,0).x()) + 0.5;
      
      const int tickLength = ticks[i].tickLength == MyTickLabel::Long
      ? TICK_LENGTH : TICK_LENGTH / 2;
      
      const Chart::AxisValue location = axis.location();
      
      if( ticks[i].tickLength != MyTickLabel::Zero )
      {
        ticksPath.moveTo(hv(dd, lineypx + (tickPos & Right ? -tickLength : 0)));
        ticksPath.lineTo(hv(dd, lineypx + (tickPos & Left ? +tickLength : 0)));
      }
      
      WPointF labelPos;
      switch( location )
      {
        case Wt::Chart::MinimumValue:
        case Wt::Chart::ZeroValue:
          labelPos = WPointF(dd, lineypx + tickLength);
          break;
          
        case Wt::Chart::MaximumValue:
          labelPos = WPointF(dd, lineypx - tickLength);
          break;
      }//switch( location )
      
      if( (properties & Chart::Labels) && !ticks[i].label.empty() )
      {
        WFlags<AlignmentFlag> labelFlags = labelHFlag;
        
        if (axis.labelAngle() == 0)
          labelFlags |= AlignCenter;
        else if (axis.labelAngle() > 0)
          labelFlags |= AlignRight;
        else
          labelFlags |= AlignLeft;
        
        bool drawText = true;
        
        
        if( m_compactAxis )
        {
          const double titlestartpx = (chartArea().right() - titlewidthpx);
          const double labelwidthpx = 10.0*ticks[i].label.toUTF8().size();
          const double labelendpx = labelPos.x() + 0.5*labelwidthpx;
          drawText = (labelendpx < titlestartpx);
        }//if( m_compactAxis )
        
        if( drawText )
        {
          const int margin = m_compactAxis ? 0 : 3;
          labelRender( painter, ticks[i].label,
                      labelPos, black, labelFlags, axis.labelAngle(), margin );
        }//if( labelRender )
      }//if( draw label text )
    }//for( size_t i = 0; i < nticks; ++i )
    
    if( properties & Chart::Line )
      painter.strokePath(ticksPath, axis.pen());
    
    if( (properties & Chart::Labels) && !axis.title().empty() )
    {
      WFont oldFont2 = painter.font();
      WFont titleFont = axis.titleFont();
      painter.setFont(titleFont);
      
      //Normal mode is vertical
      const bool vertical = orientation() == Vertical;
      
      if( m_compactAxis && vertical )
      {
        const WPen oldpen = painter.pen();
        const WBrush oldbrush = painter.brush();
        painter.setBrush( m_chartMarginBrush );
        painter.setPen( WPen(Wt::GlobalColor::transparent) );
        const double fontheight = titleFont.sizeLength().toPixels();
        const WString &title = axis.title();
        const WPointF pos(chartArea().right() - titlewidthpx, lineypx + TICK_LENGTH );
        
        const double x = pos.x(), y = lineypx + TICK_LENGTH, h = fontheight + 1;
        painter.drawRect( x, y, titlewidthpx, h );
        painter.setPen( oldpen );
        painter.setBrush( oldbrush );
        
        labelRender( painter, title, pos, m_textPen.color(), AlignTop | AlignLeft, 0, 0 );
      }else
      {
        if( vertical )
          labelRender( painter, axis.title(),
                      WPointF(chartArea().center().x(), lineypx + 22),
                      m_textPen.color(), AlignTop | AlignCenter, 0, 0 );
        else
          labelRender( painter, axis.title(),
                      WPointF(chartArea().right(), lineypx),
                      m_textPen.color(), AlignTop | AlignLeft, 0, 8 );
      }
      
      painter.setFont(oldFont2);
    }
  }
  
  painter.setFont(oldFont1);
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
  const int wpadding  = plotAreaPadding(Left)
                        + plotAreaPadding(Right);
  const int hpadding  = plotAreaPadding(Top)
                        + plotAreaPadding(Bottom);
 
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
#if( WT_VERSION >= 0x3030800 )
  const vector<WDataSeries*> seriesv = series();
  series();
#else
  const std::vector<WDataSeries> &seriesv = series();
#endif
  //  unsigned rows = m_model ? m_model->rowCount() : 0;
  
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  
  const int nrow = (!!m ? m->rowCount() : 0);
  
  for( unsigned g = 0; g < seriesv.size(); ++g )
  {
#if( WT_VERSION >= 0x3030800 )
    if( seriesv[g]->isHidden() )
      continue;
#else
    if( seriesv[g].isHidden() )
      continue;
#endif
    
    const int endSeries = g, startSeries = g;
    int i = startSeries;
    
    for( ; ; )
    {
#if( WT_VERSION >= 0x3030800 )
      WDataSeries &series = *seriesv[i];
#else
      const WDataSeries &series = seriesv[i];
#endif
     
      bool doSeries = iterator->startSeries( series, 0.0, 1, 0 );
      
      
      
      //should consider if (pxPerBin > 3*penWidth)
      
      if( doSeries )
      {
        for( int currentXSegment = 0;
            currentXSegment < axis(XAxis).segmentCount();
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
            
            int minRow = (!!m ? m->findRow( minx ) : 0);
            int maxRow = (!!m ? m->findRow( maxx ) : 0);
            
            if( minRow == maxRow )
              continue;
            
            minRow = std::max( minRow, 0 );
            maxRow = std::min( maxRow, nrow-1 );
            
            const double pxPerBin = ((maxXPx-minXPx) / (maxRow - minRow));
            const double penWidth = series.pen().width().toPixels();
#if( IOS || ANDROID )
            const bool drawHist = (pxPerBin > 2.0*penWidth);
#else
            const bool drawHist = (pxPerBin > penWidth);
#endif
//            cerr << "pxPerBin=" << pxPerBin << ", penWidth=" << penWidth << endl;
            
            iterator->startSegment( currentXSegment, currentYSegment, csa );
            
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
              const WModelIndex xIndex = m->index( row, c );
              const WModelIndex yIndex = m->index( row, column );
              
              if( !yIndex.isValid() || !xIndex.isValid() )
                continue;
              
              const double y = asNumber( m->data(yIndex) );
              
              if( !drawHist )
              {
                const double x = asNumber( m->data(xIndex) );
                iterator->newValue( series, x, y, 0, xIndex, yIndex );
              }else
              {
                //should consider if (pxPerBin > 3*penWidth)
                const double x_start = m->rowLowEdge( xIndex.row() );
                const double x_end = x_start + m->rowWidth( xIndex.row() );
                iterator->newValue( series, x_start, y, 0, xIndex, yIndex );
                iterator->newValue( series, x_end, y, 0, xIndex, yIndex );
              }//if( pxPerBin < 3.0 )
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











/*
 * See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
 */
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

SpectrumChart::SpectrumChart( Wt::WContainerWidget *parent )
: Wt::Chart::WCartesianChart( parent ),
  m_widthInPixels( 0 ),
  m_heightInPixels ( 0 ),
  m_isDragging( false ),
  m_dragAction( SpectrumChart::ZoomIn ),
  m_allowMultipleHighlightRegions( true ),
  m_allowSingleClickHighlight( false ),
  m_peakModel( NULL ),
  m_defaultPeakColor( ns_default_peakFillColor ),
  m_occupiedMarkerColor( ns_default_occupiedTimeColor ),
  m_timeHighlightColors{ ns_defaultTimeHighlightColors[0], ns_defaultTimeHighlightColors[1], ns_defaultTimeHighlightColors[2] },
  m_legend( NULL ),
  m_legendType( SpectrumChart::NoLegend ),
  m_legendNeedsRender( true ),
  m_paintedLegendWidth(-100.0f),
  m_overlayCanvas( NULL ),
  m_scrollY( 0 ),
  m_xAxisUnits( SpectrumChart::kUndefinedUnits ),
  m_compactAxis( false ),
  m_textInMiddleOfChart(),
  m_verticalLinesShowing( false ),
  m_horizontalLinesShowing( false ),
  m_avoidMobileMenu( false )
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  , m_mouseDownX( -999 ),
  m_mouseDownY( -999 ),
  m_currentDownX( -999 ),
  m_currentDownY( -999 ),
  m_mouseDownEnergy( -999.9f ),
  m_mouseDownCounts( -999.9f ),
  m_mouseDownXMin( -999.9f ),
  m_mouseDownXMax( -999.9f ),
  m_mouseDownYMin( -999.9f ),
  m_mouseDownYMax( -999.9f )
#endif
{
  addStyleClass( "SpectrumChart" );
  
  //Using ScatterPlot instead of CategoryChart since the value of the templates
  //  are reported as the sum of the templates of lower series index than,
  //  and including itself (Chart::ScatterPlot does not support stackings)
  setType( Chart::ScatterPlot );
  //m_chart->axis(Chart::XAxis).setLabelFormat( "%.2f" );
  axis(Chart::XAxis).setScale( Chart::LinearScale );
  axis(Chart::YAxis).setLabelFormat( "%.3g" );
  
  for( PeakLabels label = PeakLabels(0);
       label < kNumPeakLabels; label = PeakLabels(label+1) )
  {
    m_peakLabelsToShow[label] = false;
    m_peakLabelsColors[label] = (0xFF << 24);  //black with alpha=255
  }//for( loop over all labels )
  
//  setLayoutSizeAware( false );
  
  axis(Chart::XAxis).setMinimum( 0.0 );
  axis(Chart::XAxis).setMaximum( 3000.0 );
  
  //make it so we can use control-drag and right click on the cavas without
  //  the browsers context menu poping up
  setAttributeValue( "oncontextmenu", "return false;" );
}//SpectrumChart( constructor )


//void SpectrumChart::handleClick( Wt::WMouseEvent e )
//{
//  cerr << "Was clicked" << endl;
  //Wt-popup PopupDivMenu wt-no-reparent Wt-popupmenu Wt-outset
  //Wt-popup Wt-popupmenu Wt-outset
//  WPoint p( e.window().x, e.window().y );
//  m_rightClickMenu->setHidden( false );
//  m_rightClickMenu->popup( p );
  
//  string js = INLINE_JAVASCRIPT(
//  var fcn = function(id,x,y){
//    var e = Wt.WT.getElement(id);
//    e.style['left'] = x + 'px';
//    e.style['top'] = y + 'px';
//  }; );
//  js += " fcn('" + m_rightClickMenu->id() + "',"
//        + std::to_string(e.window().x) + ","
//        + std::to_string(e.window().y) + ");";
//  doJavaScript( js );
//}

SpectrumChart::~SpectrumChart()
{
  if( m_legend )
  {
    delete m_legend;
    m_legend = 0;
  }

  if( m_overlayCanvas )
  {
    delete m_overlayCanvas;
    m_overlayCanvas = 0;
  }
}//~SpectrumChart()

void SpectrumChart::setPeakLabelColor( SpectrumChart::PeakLabels label,
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
  
  m_peakLabelsColors[label] = (red | green | blue | alpha);
}//void setPeakLabelColor(...)


void SpectrumChart::setShowPeakLabel( SpectrumChart::PeakLabels label, bool show )
{
  m_peakLabelsToShow[label] = show;
}//void setShowPeakLabel(...)


bool SpectrumChart::showingPeakLabel( SpectrumChart::PeakLabels label ) const
{
  return m_peakLabelsToShow[label];
}


#if( WT_VERSION<=0x3030100 )
Wt::Chart::WChart2DRenderer *SpectrumChart::createRenderer( WPainter& painter,
                                                const WRectF &rectangle ) const
{
  SpectrumDataModel *specModel = dynamic_cast<SpectrumDataModel *>( model() );
  
  if( !specModel )
    throw runtime_error( "SpectrumChart::createRenderer: invalid model" );
    
    
  return new SpectrumRenderer( const_cast<SpectrumChart *>(this),
                                        painter, rectangle);
  
//  return createRenderer( painter, rectangle );
//  return Chart::WCartesianChart::createRenderer( painter, rectangle );
}//createRenderer(...)
#endif

void SpectrumChart::peakModelDataChanged( const WModelIndex &topLeft,
                                          const WModelIndex &bottomRight )
{
  //right now peakModelDataChanged(...) assumes loadPeakInfoToClient()
  //  only loads nuclide symbol and energy to client, and
  //  not kUseForShieldingSourceFit, kCandidateIsotopes, kUseForCalibration;
  //  if this changes in the future, make sure to update
  //  peakModelDataChanged(...)
  const int left = topLeft.column();
  const int right = bottomRight.column();
  
  if( topLeft.row() != bottomRight.row()
      || left <= PeakModel::kAmplitude
      || left >= PeakModel::kUserLabel
      || right >= PeakModel::kUserLabel )
    update( WFlags<PaintFlag>(0) );
  else if( left <= PeakModel::kPhotoPeakEnergy )
      loadPeakInfoToClient();
  //This next statment is poorly formed, and only catches if a single data thingy
  //  is changed - I'll fix this at some point when I'm not so tired
  else if( (left == PeakModel::kUserLabel && m_peakLabelsToShow[kShowPeakUserLabel])
          || (left == PeakModel::kIsotope && m_peakLabelsToShow[kShowPeakNuclideLabel])
          || (left == PeakModel::kPhotoPeakEnergy && m_peakLabelsToShow[kShowPeakNuclideEnergies])  )
    update( WFlags<PaintFlag>(0) );
}//void peakModelDataChanged(...)


void SpectrumChart::setPeakModel( PeakModel *model )
{
  if( m_peakModel )
    throw runtime_error( "SpectrumChart::setPeakModel(...): a model has "
                         "already been set" );
  
  if( !model )
    throw runtime_error( "SpectrumChart::setPeakModel(...): invalid input model" );
  
  m_peakModel = model;
  
  //Note that if you first remove rows from m_peakModel, then add some, then
  //  modify some data, all within one server callback, this will all only
  //  cause one update to be pushed to the user
  m_peakModel->dataChanged().connect( boost::bind( &SpectrumChart::peakModelDataChanged, this, _1, _2 ) );
  m_peakModel->rowsRemoved().connect( boost::bind( &SpectrumChart::update, this, WFlags<PaintFlag>(0) ) );
  m_peakModel->rowsInserted().connect( boost::bind( &SpectrumChart::update, this, WFlags<PaintFlag>(0) ) );
}//void setPeakModel( PeakModel *model )



//int SpectrumChart::graphableAreaWidth() const
//{
//  return m_widthInPixels - plotAreaPadding(Left) - plotAreaPadding(Right);
//}//int SpectrumChart::graphableAreaWidth() const


Wt::Signal<double,double> &SpectrumChart::xRangeChanged()
{
  return m_xRangeChanged;
}//xRangeChanged()


Wt::Signal<double,double> &SpectrumChart::rightMouseDragg()
{
  return m_rightMouseDragg;
}//rightMouseDragg()

Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> &SpectrumChart::leftClick()
{
  return m_leftClick;
}//Wt::Signal<double> &leftClick()

Wt::Signal<double,double> &SpectrumChart::doubleLeftClick()
{
  return m_doubleLeftClick;
}

Wt::Signal<double,double,int,int> &SpectrumChart::rightClick()
{
  return m_rightClick;
}//Wt::Signal<double> &rightClick()


Wt::Signal<double,double> &SpectrumChart::controlMouseMoved()
{
  return m_controlMouseMoved;
}

Wt::Signal<double,double,int,int> &SpectrumChart::controlKeyDragged()
{
  return m_controlKeyDragg;
}//controlKeyDragged();


Wt::Signal<double,double> &SpectrumChart::shiftKeyDragged()
{
  return m_shiftKeyDragg;
}//controlKeyDragged();


Wt::Signal<double,double> &SpectrumChart::shiftAltKeyDragged()
{
  return m_shiftAltKeyDragg;
}

Wt::Signal<double,double> &SpectrumChart::altKeyDragged()
{
  return m_altKeyDragg;
}

void SpectrumChart::setXAxisUnits( SpectrumChart::XAxisUnits units )
{
  m_xAxisUnits = units;
  string unitstr;
  switch( units )
  {
    case kkeV:           unitstr = "keV";     break;
    case kSeconds:       unitstr = "seconds"; break;
    case kUndefinedUnits: break;
  }//switch( units )

  if( unitstr.length() && m_overlayCanvas )
    wApp->doJavaScript( "$('#c" + m_overlayCanvas->id() + "').data('xunit','" + unitstr + "');" );
  else if( m_overlayCanvas )
    wApp->doJavaScript( "$('#c" + m_overlayCanvas->id() + "').data('xunit',null);" );
  else
    cerr << "SpectrumChart::setXAxisUnits(...)\n\tSpectrumChart::setXAxisUnits(...): overlayCanvas not set"
         << endl;
}//setXAxisUnits( SpectrumChart::XAxisUnits units )


SpectrumChart::XAxisUnits SpectrumChart::xAxisUnits() const
{
  return m_xAxisUnits;
}//XAxisUnits xAxisUnits() const


void SpectrumChart::setDragAction( const SpectrumChart::DragAction action )
{
  m_dragAction = action;
}//void setDragAction( const DragAction action )


void SpectrumChart::allowMultipleRegionHighlight( bool allow )
{
  m_allowMultipleHighlightRegions = allow;
}//void allowMultipleRegionHighlight( bool allow )



void SpectrumChart::allowSingleClickHighlight( bool allow )
{
  m_allowSingleClickHighlight = allow;
}//void allowSingleClickHighlight( bool allow )


void SpectrumChart::showGridLines( const bool draw )
{
  Chart::WCartesianChart::axis(Chart::XAxis).setGridLinesEnabled( draw );
  Chart::WCartesianChart::axis(Chart::YAxis).setGridLinesEnabled( draw );
  update(); //force re-render
}//void showGridLines( const bool draw )

void SpectrumChart::showVerticalLines( const bool draw )
{
  Chart::WCartesianChart::axis(Chart::XAxis).setGridLinesEnabled( draw );
  update(); //force re-render
  m_verticalLinesShowing = draw;
}

void SpectrumChart::showHorizontalLines( const bool draw )
{
  Chart::WCartesianChart::axis(Chart::YAxis).setGridLinesEnabled( draw );
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
  return WLength( m_heightInPixels, WLength::Pixel );
}//Wt::WLength paintedHeight() const



Wt::WLength SpectrumChart::paintedWidth() const
{
  return WLength( m_widthInPixels, WLength::Pixel );
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
  
  const int old_pad = plotAreaPadding(Wt::Bottom);
  const int new_pad = std::max(old_pad + (compact ? -1 : 1)*17, 0);
  
  m_compactAxis = compact;
  setPlotAreaPadding( new_pad, Bottom );
  if( m_overlayCanvas )
    m_overlayCanvas->setChartPadding();
}

bool SpectrumChart::isAxisCompacted() const
{
  return m_compactAxis;
}

void SpectrumChart::wtDragged( Wt::WMouseEvent /*event*/ )
{
  //This function handles the standard Wt mouseMoved() signal
  m_isDragging = true;
}//void wtDragged( int x0, int y0, int dx, int dy )


void SpectrumChart::wtMouseUp( Wt::WMouseEvent event )
{
  //This function handles the standard Wt mouseWentUp() signal
  const int x1 = event.widget().x;
  const int y1 = event.widget().y;
  const int dx = event.dragDelta().x;
  const int dy = event.dragDelta().y;
  const int x0 = x1 - dx;
  const int y0 = y1 - dy;


  const bool wasDragging = m_isDragging;
  m_isDragging = false;

  if( !wasDragging && (event.button() == WMouseEvent::RightButton) )
  {
    handleRightClick( x1, y1, event.modifiers(),
                     event.document().x, event.document().y );
    return;
  }
  
  if( !wasDragging && (!m_allowSingleClickHighlight || (m_dragAction!=HighLight)) )
    return;

  if( !wasDragging )
    handleLeftClick( x1, y1, event.modifiers(),
                     event.document().x, event.document().y );
  else
    handleDrag( x0, y0, x1, y0,
                event.window().x-event.widget().x,
                event.window().y-event.widget().y,
                event.modifiers(), event.button(), -1 );
}//void wtMouseUp( Wt::WMouseEvent event )



#if( USE_OverlayDragEvent )
void SpectrumChart::userDragged( const OverlayDragEvent &e )
{
  handleDrag( e.x0, e.y0, e.x1, e.y1,
              e.offsetLeft, e.offsetTop,
              e.keyModifiers, WMouseEvent::Button(e.button), e.dt );
}//void userDragged( OverlayDragEvent event )
#else
void SpectrumChart::userDragged( int x0, int y0, Wt::WMouseEvent event )
{
  const int x1 = event.widget().x;
  const int y1 = event.widget().y;

  handleDrag( x0, y0, x1, y1,
              event.window().x-event.widget().x,
              event.window().y-event.widget().y,
              event.modifiers(), event.button(), -1 );
}//void userDragged( int x0, int x1, Wt::WMouseEvent event )
#endif  //#if( USE_OverlayDragEvent )



Wt::JSignal<Wt::WKeyEvent> &SpectrumChart::keyPressedWhileMousedOver()
{
  if( !m_overlayCanvas )
  {
    throw std::runtime_error( "SpectrumChart::keyPressedWhileMousedOver(): "
                              "you may only call this function if you have "
                              "first called "
                              "SpectrumChart::enableOverlayCanvas(bool,bool)!" );
  }//if( !m_overlayCanvas )

  return m_overlayCanvas->keyPressWhileMousedOver();
}//Wt::JSignal<Wt::WKeyEvent> &keyPressedWhileMousedOver()


void SpectrumChart::allowArrowToMoveSingleClickRegion( bool allow )
{
  if( !allow )
  {
    if( !m_overlayCanvas )
    {
      cerr << "SpectrumChart::allowArrowToMoveSingleClickRegion( bool allow ): "
           << "some funny error: signal is connected, but m_overlayCanvas==NULL"
           << " - you probably willl never see this message" << endl;
      return;
    }//if( !m_overlayCanvas )

    if( m_arrowRespondConnection.connected() )
      m_overlayCanvas->keyPressWhileMousedOver().disconnect( m_arrowRespondConnection );
    return;
  }//if( !allow )

  if( !m_overlayCanvas )
  {
    throw std::runtime_error( "SpectrumChart::allowArrowToMoveSingleClickRegion(bool): "
                              "you may only call this function if you have "
                              "first called "
                              "SpectrumChart::enableOverlayCanvas(bool,bool)!" );
  }//if( !m_overlayCanvas )

  m_arrowRespondConnection = m_overlayCanvas->keyPressWhileMousedOver().connect( boost::bind( &SpectrumChart::handleArrowPress, this, _1 ) );
}//void allowArrowToMoveSingleClickRegion( bool allow )


void SpectrumChart::handleArrowPress( const Wt::WKeyEvent &event )
{
  if( !m_allowSingleClickHighlight )
    return;

  if( (event.key() != Key_Left) && (event.key() != Key_Right) )
    return;

  vector< HighlightRegion * > highlights;
  for( size_t i = 0; i < m_highlights.size(); ++i )
    if( m_highlights[i].hash == size_t(ForegroundHighlight) )
      highlights.push_back( &(m_highlights[i]) );
  
  if( highlights.size() != 1 )
    return;
  
  if( event.modifiers()!=0x0 && event.modifiers() != ShiftModifier )
    return;
  
  const double lowerx = highlights[0]->lowerx;
  const double upperx = highlights[0]->upperx;
  SpectrumDataModel *th1model = dynamic_cast<SpectrumDataModel *>( model() );

  if( !th1model )
    return;

  double newlowerx = 0.0, newupperx = 0.0;
  
  const int nrowmove = ((event.key()==Key_Left) ? -1 : 1);
  int lowrow = th1model->findRow( lowerx+(1.0E-6) );
  int highrow = th1model->findRow( (float)upperx-(1.0E-6) );

  if( event.modifiers() != ShiftModifier )
  {
    lowrow += nrowmove;
    highrow += nrowmove;
  }else
  {
    if( event.key() == Key_Left )
      lowrow += nrowmove;
    if(event.key() == Key_Right )
      highrow += nrowmove;
  }//if( event.modifiers() != ShiftModifier ) / else
  
  lowrow = std::max( lowrow, 0 );
  highrow = std::min( highrow, th1model->rowCount()-1 );
  newlowerx = th1model->rowLowEdge( lowrow );
  newupperx = th1model->rowLowEdge( highrow ) + th1model->rowWidth( highrow );
  
  highlights[0]->lowerx = static_cast<float>( newlowerx );
  highlights[0]->upperx = static_cast<float>( newupperx );

  update();//force re-rendering

  //XXX
  //   see XXX note in handleLeftClick(...)!
  //    m_xRangeChanged.emit( lowerx, upperx );
  m_xRangeChanged.emit( newlowerx, newupperx );
}//void SpectrumChart::handleArrowPress( Wt::WKeyEvent &event )


void SpectrumChart::handleDoubleTap( int x, int y )
{
  const WPointF point1( x, y );
  const WPointF modelPoint1 = mapFromDevice( point1, Chart::OrdinateAxis );
  m_doubleLeftClick.emit( modelPoint1.x(), modelPoint1.y() );
}//void handleDoubleTap( int x, int y )


void SpectrumChart::handleDoubleLeftClick( Wt::WMouseEvent event )
{
  if( event.button()!=WMouseEvent::LeftButton )
    return;

  const int x1 = event.widget().x;
  const int y1 = event.widget().y;
  const WPointF point1( x1, y1 );
  const WPointF modelPoint1 = mapFromDevice( point1, Chart::OrdinateAxis );
  
  m_doubleLeftClick.emit( modelPoint1.x(), modelPoint1.y() );
}//void handleDoubleLeftClick( Wt::WMouseEvent event )


void SpectrumChart::handleLeftClick( const int x1, const int y1, const int modifier,
                                     const int pageX, const int pageY )
{
  //This function assumes the user actually single clicked, and didnt drag
  //  or anything else like that (and not the standard Wt mouseWentUp() signal).  
  const WPointF point1( x1, y1 );
  const WPointF modelPoint1 = mapFromDevice( point1, Chart::OrdinateAxis );

  if( modifier == 0x0 )
    m_leftClick.emit( modelPoint1.x(), modelPoint1.y(), pageX, pageY );
  
  if( !m_allowSingleClickHighlight )
    return;

  SpectrumDataModel *th1model = dynamic_cast<SpectrumDataModel *>( model() );

  if( !th1model )
  {
    cerr << "Super lame programmer logic error SpectrumChart::wtMouseUp(...)" << endl;
    return;
  }//if( !th1model )
  
  const int row = th1model->findRow( modelPoint1.x() );

  if( row < 0 || row >= th1model->rowCount() )
    return;
  
  const double lowerx = th1model->rowLowEdge( row );
  const double upperx = lowerx + th1model->rowWidth( row );
  
  if( !m_allowMultipleHighlightRegions || (modifier != ShiftModifier) )
    m_xRangeChanged.emit( lowerx, upperx );
  else
    m_shiftKeyDragg.emit( (lowerx+upperx)/2.0, (lowerx+upperx)/2.0 );
}//void handleLeftClick( int x, int y )


void SpectrumChart::handleRightClick( const int x1, const int y1,
                                      const int ,//modifier
                                      const int pageX, const int pageY )
{
  const WPointF point1( x1, y1 );
  const WPointF modelPoint = mapFromDevice( point1, Chart::OrdinateAxis );

  m_rightClick.emit( modelPoint.x(), modelPoint.y(), pageX, pageY );
}//void handleRightClick( Wt::WMouseEvent event )


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
      tokeep.push_back( m_highlights[i] );
  
  for( size_t i = 0; i < p.size(); ++i )
  {
    HighlightRegion region;
    region.lowerx = static_cast<float>( p[i].first );
    region.upperx = static_cast<float>( p[i].second );
    region.color = color;
    region.hash = size_t(type);
    tokeep.push_back( region );
  }
  
  tokeep.swap( m_highlights );
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
    region.lowerx = static_cast<float>( pe.first );
    region.upperx = static_cast<float>( pe.second );
    region.color = m_occupiedMarkerColor;
    region.hash = 0;//size_t(type);
    m_occupancy_regions.push_back( region );
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
      m_highlights.erase( m_highlights.begin() + i );
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
  boost::hash_combine( region.hash, lowerx );
  boost::hash_combine( region.hash, upperx );
  boost::hash_combine( region.hash, color.red() );
  boost::hash_combine( region.hash, color.green() );
  boost::hash_combine( region.hash, color.blue() );
  boost::hash_combine( region.hash, color.alpha() );
  
  if( region.hash <= 2 )
    region.hash += 3;
  
  //should in principle check for collision, but whatever
  
  m_highlights.push_back( region );
  
  return region.hash;
}//void addDecorativeHighlightRegion(...)


void SpectrumChart::handleDrag( int _x0, int _y0, int _x1, int _y1,
                                int offsetLeft, int offsetTop,
                               Wt::WFlags<Wt::KeyboardModifier> modifiers,
                               Wt::WMouseEvent::Button button,
                               int dt )
{
  //What this function does:
  // 1) If user dragged mouse less than 'mouse_epsilon' in either x or y
  //    direction, than do nothing
  // 2) If user drags to the right more than up/down, then zoom in the 'x' axis
  // 3) If user drags down more 'y' than 'x', than zoom in on the 'y' axis
  // 4) If user drags to the left, more than up/down, then zoom-out; for
  //    highlighting, add region to highlighted regions

  // This compensates for the canvas not moving when the user scrolls up or down
  // In the Anthony app, the user may vertically scroll the chart.  Although
  // the chart vertically moves, we don't vertically scroll the overlay
  // canvas, because that might overlap with the tabs, making the tabs
  // non-clickable.  The workaround is to detect how much the user has
  // scrolled in the y axis (in number of pixels), and compensate for that
  // in this function, which is a simple y-axis translation, in pixels.
  _y0 += m_scrollY;
  _y1 += m_scrollY;

   
  const double mouse_epsilon = std::max( 5.0, m_widthInPixels/model()->rowCount() );
  const double left_mouse_epsilon = 10.0;

  const double x0 = static_cast<double>( _x0 );
  const double y0 = static_cast<double>( _y0 );
  const double x1 = static_cast<double>( _x1 );
  const double y1 = static_cast<double>( _y1 );
  const double dx = x1 - x0;
  const double dy = y1 - y0;

  const WPointF point0( x0, y0 );
  const WPointF point1( x1, y1 );
  const WPointF modelPoint0 = mapFromDevice( point0, Chart::OrdinateAxis );
  const WPointF modelPoint1 = mapFromDevice( point1, Chart::OrdinateAxis );

  //Make it so if there are multiple modifiers, than nothing is done
  const int nmodifiers = ((modifiers & Wt::ControlModifier)==Wt::ControlModifier)
                         + ((modifiers & Wt::ShiftModifier)==Wt::ShiftModifier)
                         + ((modifiers & Wt::AltModifier)==Wt::AltModifier)
//                         + (button==WMouseEvent::RightButton)
                         ;
  const unsigned int shiftAltMod = (Wt::ShiftModifier | Wt::AltModifier);
  const bool isShiftAlt = ((modifiers & shiftAltMod) == shiftAltMod);

  if( m_dragAction == HighLight )
  {
    if( button != WMouseEvent::LeftButton )
      return;
    if( (modifiers != 0x0) && (modifiers != Wt::ShiftModifier)
        && !isShiftAlt && (modifiers!=Wt::AltModifier) )
      return;
  }//if( m_dragAction == HighLight )

  
  if( (nmodifiers > 1) && !isShiftAlt )
    return;

  //case 1
  if( ((dx>0.0 && dx<mouse_epsilon)
       || (dx<0.0 && fabs(dx)<left_mouse_epsilon))
     && (fabs(dy) < mouse_epsilon) )
  {
    cerr << "\nDrag epsilon too small; doing nothing." << endl;
    return;
  }//case 1
  
  
  if( dt > 0 && dt < 100 ) //This 100 should mach same value in 'OverlayOnMouseUp' JS function as well
  {
    cerr << "\nDrag to short of time (" << dt << "); not doing anything" << endl;
    return;
  }
  
  if( button==WMouseEvent::RightButton )
  {
    if( nmodifiers == 0 )
      m_rightMouseDragg.emit( modelPoint0.x(), modelPoint1.x() );
    return;
  }//if( button==WMouseEvent::RightButton )

  if( modifiers == Wt::ControlModifier )
  {
    m_controlKeyDragg.emit( modelPoint0.x(), modelPoint1.x(),
                            offsetLeft + _x1,
                            offsetTop + _y1 );
    return;
  }//if( modifiers | Wt::ControlModifier )

  if( modifiers == Wt::AltModifier
     || button == WMouseEvent::RightButton )
  {
    if( m_dragAction != HighLight )
      handleAltDrag( modelPoint1.x() - modelPoint0.x() );
    else
      m_altKeyDragg.emit( modelPoint0.x(), modelPoint1.x() );
    return;
  }//if( modifiers == Wt::AltModifier )

  
  if( modifiers == Wt::ShiftModifier )
  {
    m_shiftKeyDragg.emit( modelPoint0.x(), modelPoint1.x() );
    
    if( m_dragAction != HighLight )
      return;
  }//if( modifiers == Wt::ShiftModifier )

  if( isShiftAlt )
  {
    m_shiftAltKeyDragg.emit( modelPoint0.x(), modelPoint1.x() );
    return;
  }//if( isShiftAlt )

  //case 2
  if( (dx > 0.0) && (fabs(dx) >= fabs(dy)) )
  {
    switch( m_dragAction )
    {
      case ZoomIn:
      {
        //In principle, we should probably on zoom into regions that correspond
        //  to the edges of bins.
        const SpectrumDataModel * const theModel = dynamic_cast<const SpectrumDataModel *>( model() );

        if( !theModel ) // This should never actually happen.
          break;

        const std::shared_ptr<const Measurement> xHist = theModel->histUsedForXAxis();

        if( !xHist )
          break;

        //We want to keep the user from zooming in past were it is useful;
        //  e.g. to an area less than the width of a bin
        //  XXX - This method of controlling this could use improvment as well
        //        as some kind of notification to the user upon failure
        //
        
        const size_t channel_X0 = xHist->find_gamma_channel( (float)modelPoint0.x() );
        const size_t channel_X1 = xHist->find_gamma_channel( (float)modelPoint1.x() );
        const size_t dbin = ((channel_X1 > channel_X0) ? (channel_X1 - channel_X0) : size_t(0));
        if( dbin > 3 )
        {
          setXAxisRange( modelPoint0.x(), modelPoint1.x() );
        }else if( dbin > 0 )
        {
          const double delta = modelPoint1.x() - modelPoint0.x();
          const double mindelta = 0.5*(xHist->gamma_channel_width(channel_X0)
                                       + xHist->gamma_channel_width(channel_X1));
          if( delta > mindelta )
            setXAxisRange( modelPoint0.x(), modelPoint1.x() );
        }//if( binX0 != binX1 )

//        setAutoYAxisRange();
        m_xRangeChanged.emit( modelPoint0.x(), modelPoint1.x() );
        break;
      }//case ZoomIn:


      case HighLight:
        if( !m_allowMultipleHighlightRegions || (modifiers != ShiftModifier) )
          m_xRangeChanged.emit( modelPoint0.x(), modelPoint1.x() );
//Next lines commented out because this signal has already been emited
//        else
//          m_shiftKeyDragg.emit( modelPoint0.x(), modelPoint1.x() );
      break;
      case NoAction:
      break;
    };//switch( m_dragAction )

    return;
  }//case 2

  //case 3
  if( fabs(dx) < fabs(dy) )
  {
    switch( m_dragAction )
    {
      case ZoomIn:
        if( dy > 0.0 )
        {
          //If the y-axis range gets too small, we can get a JS error on the
          //  client side - lets protect against this.  Note that this could
          //  cause problems on extremely scalled spectrums but I expect this
          //  to never actually happen to anyone
          const double ydiff = fabs( modelPoint1.y()-modelPoint0.y() );
          if( ydiff > 0.005
          || (axis(Chart::Y1Axis).scale()==Chart::LinearScale && ydiff>1.0E-8) )
            setYAxisRange( modelPoint1.y(), modelPoint0.y() );
        }
        else
        {
          const int can_height = (int)paintedHeight().toPixels();

          if( fabs(dy) >= 0.075*can_height )  //XXX the 0.075 and 0.05 must be same as in javascript for this app
          {
            setAutoYAxisRange();
          }else
          {
            SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model() );

            if( !theModel )
              return;

            const double x0 = axis(Chart::XAxis).minimum();
            const double x1 = axis(Chart::XAxis).maximum();
            const double old_y0 = axis(Chart::YAxis).minimum();
            const double old_y1 = axis(Chart::YAxis).maximum();
            const double old_range = old_y1 - old_y0;

            //XXX the below (and above to some extent) shares enough logic
            //    with setAutoYAxisRange() that it should be refactored to
            //    ensure consistency in future code upgrades

            double data_y0 = 0.0, data_y1 = 0.0;
            double min_disp = 0.0, max_disp = 0.0;

            theModel->yRangeInXRange( x0, x1, data_y0, data_y1 );

            if( axis(Chart::OrdinateAxis).scale() == Chart::LogScale )
            {
              if( data_y0 <= 0.0 ) min_disp = 0.1;
              else                 min_disp = 0.1*data_y0;
              max_disp = 2.5*data_y1;
            }else //if( axis(Chart::YAxis).scale() == Chart::LinearScale )
            {
              if( data_y0 <= 0.0 ) min_disp = 1.1*data_y0;
              else                 min_disp = 0.9*data_y0;
              max_disp = 1.1*data_y1;
            }//if( Y is ogScale ) / else

            const double centroid = old_y0 + 0.5*old_range;
            const bool zoomX2 = (fabs(dy) < 0.05*can_height);
            const double mult = (zoomX2 ? 1.0 : 2.0);
            const double new_y0 = max(centroid - mult*old_range, min_disp);
            const double new_y1 = min(new_y0 + 2.0*mult*old_range, max_disp);

            setYAxisRange( new_y0, new_y1 );
          }//if( do full zoom out ) / else due to or 4
        }//if( zoom-in ) / else ( zoom-out )
      break;

      case HighLight:
      case NoAction:
      break;
    };//switch( m_dragAction )

    return;
  }//case 3

  //case 4
  if( (dx < 0.0) && (fabs(dx) >= fabs(dy)) )
  {
    switch( m_dragAction )
    {
      case ZoomIn:
      {
        const int can_width = (int)paintedWidth().toPixels();

        if( -dx >= 0.04*can_width ) //XXX this 0.04 and 0.02 need to be same as with the JavaScript
        {
          setAutoXAxisRange();
        }else
        {
          SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model() );
          std::shared_ptr<const Measurement> axisH = (theModel ? theModel->histUsedForXAxis()
                                       : std::shared_ptr<const Measurement>());
          const bool hasSpectrum = (axisH && axisH->num_gamma_channels() > 3);
          const double hist_min = (hasSpectrum ? axisH->gamma_energy_min() : 0.0f );
          const double hist_max = (hasSpectrum ? axisH->gamma_energy_max() : 3000.0f);

          const double old_xmin = axis(Chart::XAxis).minimum();
          const double old_xmax = axis(Chart::XAxis).maximum();

          const double old_range = old_xmax - old_xmin;
          const double centroid = old_xmin + 0.5*old_range;

          const bool is2xZoom = (-dx < 0.02*can_width); //else 4x xoom

          const double mult = (is2xZoom ? 1.0 : 2.0);
          const double new_xmin = max(centroid - mult*old_range, hist_min);
          const double new_xmax = min(new_xmin + 2.0*mult*old_range, hist_max);

          axis(Chart::XAxis).setMinimum( new_xmin );
          axis(Chart::XAxis).setMaximum( new_xmax );
        }//if( full zoom out ) / else

        setAutoYAxisRange();
        m_xRangeChanged.emit( axis(Chart::XAxis).minimum(),
                             axis(Chart::XAxis).maximum() );
        break;
      }//case ZoomIn:

      case HighLight:
        if( !m_allowMultipleHighlightRegions || (modifiers != ShiftModifier) )
          m_xRangeChanged.emit( modelPoint0.x(), modelPoint1.x() );
//Next couple lines commented out because this signal should have been emmitted
//        else
//          m_shiftKeyDragg.emit( modelPoint0.x(), modelPoint1.x() );
      break;

      case NoAction:
      break;
    };//switch( m_dragAction )

    return;
  }//case 4

  //case 5
  //...
}//void handleDrag( Wt::WMouseEvent event )


void SpectrumChart::handleAltDrag( const double dx )
{
  const double oldmin = axis(Chart::XAxis).minimum();
  const double oldmax = axis(Chart::XAxis).maximum();
  double newmin = oldmin - dx;
  double newmax = oldmax - dx;

  //Check and make sure the user isnt dragging past where the axis is defined
  SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model() );
  if( theModel )
  {
    std::shared_ptr<const Measurement> axisH = theModel->histUsedForXAxis();

    if( axisH
        && (axisH->energy_calibration_model() != SpecUtils::EnergyCalType::InvalidEquationType) )
    {
      if( dx < 0 )
      {
        const double max_energy = axisH->gamma_energy_max();
        if( newmax > max_energy )
        {
          newmax = max_energy;
          newmin = max_energy - (oldmax-oldmin);
        }//if( newmax > max )
      }else
      {
        const double min_energy = axisH->gamma_energy_min();
        if( newmin < min_energy )
        {
          newmin = min_energy;
          newmax = min_energy + (oldmax-oldmin);
        }//if( newmin < min )
      }//if( dragging to the left ) / else (dragging to the right )
    }//if( axisH )
  }//if( theModel )

  setXAxisRange( newmin, newmax );
  m_xRangeChanged.emit( newmin, newmax );
}//void handleAltDrag( const double x1, const double x2 )


#if( USE_HIGH_BANDWIDTH_INTERACTIONS )

void SpectrumChart::setRenderWaitingStatus( const bool waiting ) const
{
  if( !m_overlayCanvas )
    return;
    
  const string js = "$('#c" + m_overlayCanvas->id() + "').data('RenderWaiting',"
                    + string( waiting ? "true" : "false" ) + ");";
  
  wApp->doJavaScript( js, true );
//  wApp->doJavaScript( "setTimeout(function(){" + js + "},10);",true);
}//void setRenderWaitingStatus( const bool waiting ) const;


bool SpectrumChart::checkHighBandWidthInterActionApplicable() const
{
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  const bool app = ((m_dragAction==ZoomIn) && m && !!m->histUsedForXAxis());
  if( !app )
    setRenderWaitingStatus( false );
  
  return app;
}//bool checkHighBandWidthInterActionApplicable() const


void SpectrumChart::handleMouseLeft()
{
  m_mouseDownX = m_mouseDownY = -999;
  m_currentDownX = m_currentDownY = -999;
  m_mouseDownEnergy = m_mouseDownCounts = -999.9f;
  m_mouseDownXMin = m_mouseDownXMax = -999.9f;
  m_mouseDownYMin = m_mouseDownYMax = -999.9f;
  
  setRenderWaitingStatus( false );
}//void handleMouseLeft()


void SpectrumChart::handlePinchZoomChange( int x0_t0, int x_t0,
                                           int x0_t1, int x_t1, int y )
{
  if( !checkHighBandWidthInterActionApplicable() )
    return;
  
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  if( !m )
    return;
  
  std::shared_ptr<const Measurement> axisH = m->histUsedForXAxis();
  if( !axisH || y != 0 )
  {
    setRenderWaitingStatus( false );
    return;
  }
  
  //The below assumes the x-axis is a continuous, linear scale (which is
  //  true for InterSpec, but not Wt::Chart::WCartesianChart in general).
  const double factor = fabs((double)(x_t0 - x_t1)) / fabs((double)(x0_t0 - x0_t1));
  
  const double xaxismin = axisH->gamma_energy_min();
  const double xaxismax = axisH->gamma_energy_max();
  const double initialDE = fabs(m_mouseDownXMax - m_mouseDownXMin);
  
  const float currentXMin = static_cast<float>( axis(Chart::XAxis).minimum() );
  const float currentXMax = static_cast<float>( axis(Chart::XAxis).maximum() );
  const float currentYMin = static_cast<float>( axis(Chart::OrdinateAxis).minimum() );
  const float currentYMax = static_cast<float>( axis(Chart::OrdinateAxis).maximum() );

  const WPointF lowerleft = mapToDevice( currentXMin, currentYMin );
  const WPointF upperright = mapToDevice( currentXMax, currentYMax );
  const double dpx = upperright.x() - lowerleft.x();
  if( dpx==0.0 || IsInf(dpx) || IsNan(dpx) || m_mouseDownXMax<0.0f )
  {
    setRenderWaitingStatus( false );
    return;
  }
  
  const double origDE = m_mouseDownXMax - m_mouseDownXMin;
  const double t0E0 = m_mouseDownXMin + origDE*((x0_t0-lowerleft.x())/dpx);
  const double t1E0 = m_mouseDownXMin + origDE*((x0_t1-lowerleft.x())/dpx);
  const double origMidE = 0.5*(t0E0 + t1E0);
  const double midPxFrac = (0.5*(x_t0+x_t1) - lowerleft.x())/dpx;
  
  double newxmin = xaxismin, newxmax = xaxismax;
  
  if( factor > 1.0 )
  {
    //zooming-in - fingers are further apart then when they startet
    //0.5 below is arbitrary, just to speed up the zooming
    const double newDE = 0.5 * initialDE / factor;
    
    //origMidE should be at midPx, and the total x-range should span newDE
    newxmin = std::max( origMidE - midPxFrac*newDE, xaxismin );
    newxmax = std::min( origMidE + (1.0-midPxFrac)*newDE, xaxismax );
  }else
  {
    //zooming-out - fingers are closer together then when they started
    const double maxDE = xaxismax - xaxismin;
    const double originalzoom = (m_mouseDownXMax - m_mouseDownXMin) / maxDE;
    
    //we want to get to a zoom of 1.0 by the time fingers are 25px apart
    const double delta = 1.0 - originalzoom;
    const double originalDPx = fabs((double)(x0_t0 - x0_t1));
    const double currentDPx = fabs((double)(x_t0 - x_t1));
    const double frac = 1.0 - (currentDPx-25.0)/(originalDPx-25.0);
    
    double zoom = originalzoom + delta*frac;
    
    zoom = std::max( zoom, originalzoom/factor );
    zoom = std::min( zoom, 1.0 );
    const double newDE = zoom * maxDE;
    
    double shift = 0.0;
    newxmin = origMidE - midPxFrac*newDE;
    newxmax = origMidE + (1.0-midPxFrac)*newDE;
    
    if( newxmin < xaxismin )
      shift = xaxismin - newxmin;
    else if( newxmax > xaxismax )
      shift = xaxismax - newxmax;
    
    newxmin += shift;
    newxmax += shift;
  }//if( factor < 1.0 )
  
  
  if( IsInf(newxmin) || IsNan(newxmin) || IsInf(newxmax) || IsNan(newxmax)
      || newxmin==newxmax || newxmax<newxmin )
  {
    newxmin = m_mouseDownXMin;
    newxmax = m_mouseDownXMax;
  }
  
  setXAxisRange( newxmin, newxmax );
  m_xRangeChanged.emit( newxmin, newxmax );
}//handlePinchZoomChange(...)


void SpectrumChart::handleMouseEnter()
{
  handleMouseLeft();
}


void SpectrumChart::handleMouseDown( int x, int y, int modifiers )
{
  if( !checkHighBandWidthInterActionApplicable() )
    return;
  
  y += m_scrollY;
  
  const WPointF point = mapFromDevice( WPointF(x,y), Chart::OrdinateAxis );

//  const int ntouches = (modifiers / 1000);
  
  m_mouseDownX = x;
  m_mouseDownY = y;
  m_currentDownX = x;
  m_currentDownY = y;
  m_mouseDownEnergy = point.x();
  m_mouseDownCounts = point.y();

  m_mouseDownXMin = static_cast<float>( axis(Chart::XAxis).minimum() );
  m_mouseDownXMax = static_cast<float>( axis(Chart::XAxis).maximum() );
  m_mouseDownYMin = static_cast<float>( axis(Chart::OrdinateAxis).minimum() );
  m_mouseDownYMax = static_cast<float>( axis(Chart::OrdinateAxis).maximum() );
  
  setRenderWaitingStatus( false );
}//void handleMouseDown( int x, int y, int modifiers )


void SpectrumChart::handleMouseUp( int x, int y, int modifiers, int dt )
{
  if( !checkHighBandWidthInterActionApplicable() )
    return;
  
//  if( dt < 120 )
//    return;
  
  handleMouseLeft();
}//void handleMouseUp( int x, int y, int modifiers )


void SpectrumChart::handleLeftMouseMove( int x, int y, int dt )
{
  if( !checkHighBandWidthInterActionApplicable() )
    return;
    
  if( x > m_mouseDownX )
    return;
  
  if( (m_mouseDownX == -999)
     || (m_currentDownX==x && m_currentDownY==y) )
  {
    setRenderWaitingStatus( false );
    return;
  }

  //The zooming out on the y-axis is a bit messed up
  //Should implement handling the user hit the escape key while zooming out
  const int dx = m_mouseDownX - x;
  const int dy = m_mouseDownY - y;
  
  if( (abs(dy) > abs(dx)) )
  {
    if( abs(dy) > 10 )
    {
      const float xmin = static_cast<float>( axis(Chart::XAxis).minimum() );
      const float xmax = static_cast<float>( axis(Chart::XAxis).maximum() );
  
      if( fabs(m_mouseDownXMin-xmin)>0.001 || fabs(m_mouseDownXMax-xmax)>0.001 )
      {
        setXAxisRange( m_mouseDownXMin, m_mouseDownXMax );
        m_xRangeChanged.emit( m_mouseDownXMin, m_mouseDownXMax );
        setYAxisRange( m_mouseDownYMin, m_mouseDownYMax ); //doesnt seem to be working
      }
    }//if( abs(dy) > 10 )
    
    setRenderWaitingStatus( false );
    
    return;
  }//if( (abs(dy) > abs(dx)) && (abs(dy) > 20) )
  
  
  m_currentDownX = x;
  m_currentDownY = y;
  
  y += m_scrollY;
  
  double xaxismin = 0.0, xaxismax = 3000.0;
  
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  if( !m )
    return;
  
  std::shared_ptr<const Measurement> axisH = m->histUsedForXAxis();
  if( !!axisH )
  {
    xaxismin = axisH->gamma_energy_min();
    xaxismax = axisH->gamma_energy_max();
  }//if( !!axisH )
  
  double frac = 4.0 * (m_mouseDownX - x) / double(m_mouseDownX);
  if( IsInf(frac) || IsNan(frac) || frac < 0.0 )
    frac = 0.0;
  
  const double origdx = m_mouseDownXMax - m_mouseDownXMin;
  const double maxdx = xaxismax - xaxismin;
  const double newdx = origdx + frac*(maxdx - origdx);
  
  const double deltadx = newdx - origdx;
  const double newxmin = std::max( m_mouseDownXMin - 0.5*deltadx, xaxismin );
  const double newxmax = std::min( m_mouseDownXMax + 0.5*deltadx, xaxismax );
    
  setXAxisRange( newxmin, newxmax );
  m_xRangeChanged.emit( newxmin, newxmax );
}//void handleLeftMouseMove( int x, int y )


void SpectrumChart::handleAltLeftMouseMove( int x, int y, int dt )
{
  if( !checkHighBandWidthInterActionApplicable() )
    return;

  if( (m_mouseDownX == -999)
      || (m_currentDownX==x && m_currentDownY==y) )
  {
    setRenderWaitingStatus( false );
    return;
  }
  
  m_currentDownX = x;
  m_currentDownY = y;
  
  y += m_scrollY;
  
  const WPointF point = mapFromDevice( WPointF(x,y), Chart::OrdinateAxis );
  
  const double delta = m_mouseDownEnergy - point.x();
  const double oldmin = axis(Chart::XAxis).minimum();
  const double oldmax = axis(Chart::XAxis).maximum();

  double newmin = oldmin + delta;
  double newmax = oldmax + delta;
  
  double xaxismin = 0.0, xaxismax = 3000.0;
  
  //checked in checkHighBandWidthInterActionApplicable() already
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  if( !m )
    return;
  
  std::shared_ptr<const Measurement> axisH = m->histUsedForXAxis();
  if( !!axisH )
  {
    xaxismin = axisH->gamma_energy_min();
    xaxismax = axisH->gamma_energy_max();
  }//if( !!axisH )
  
  if( newmin < xaxismin )
  {
    newmin = xaxismin;
    newmax = newmin + (oldmax-oldmin);
  }
  
  if( newmax > xaxismax )
  {
    newmax = xaxismax;
    newmin = newmax - (oldmax-oldmin);
  }
 
  setXAxisRange( newmin, newmax );
  m_xRangeChanged.emit( newmin, newmax );
}//void handleAltLeftMouseMove( int x, int y );
#endif  //#if( USE_HIGH_BANDWIDTH_INTERACTIONS )


void SpectrumChart::handleControlMouseMove( int x_start, int x_finish )
{
  try
  {
    const WPointF point0( static_cast<double>( x_start ), 0.0 );
    const WPointF point1( static_cast<double>( x_finish ), 0.0 );
    const WPointF modelPoint0 = mapFromDevice( point0, Chart::OrdinateAxis );
    const WPointF modelPoint1 = mapFromDevice( point1, Chart::OrdinateAxis );

    m_controlMouseMoved.emit( modelPoint0.x(), modelPoint1.x() );
  }catch(...)
  {
    cerr << "Caught exception in handleControlMouseMove(...)" << endl;
  }
}//void handleControlMouseMove( int x_start )




void SpectrumChart::setXAxisRange( const double x1, const double x2 )
{
  const double xrange = x2 - x1;
  
  axis(Chart::XAxis).setRange( x1, x2 );
  
  //These labels dont
  if( xrange < 1.0 )
    axis(Chart::XAxis).setLabelFormat( "%.2f" );
  else if( xrange < 10.0 )
    axis(Chart::XAxis).setLabelFormat( "%.1f" );
  else
    axis(Chart::XAxis).setLabelFormat( "%.0f" );
}//void setXAxisRange()


void SpectrumChart::setYAxisRange( const double y1, const double y2 )
{
  axis(Chart::OrdinateAxis).setRange( y1, y2 );
}//void setYAxisRange()


//TODO: note that setting autolimits on a zoomed in section doesnt
//      really do what we want (eg rescale y-axis according to visible
//      range), therefore you should call setAutoYAxisRange()
void SpectrumChart::setAutoXAxisRange()
{
//  axis(Chart::XAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  //XXX
  //  WAxis doesnt immediately recalculate X and Y axis ranges when
  //  setAutoXAxisRange() and setAutoYAxisRange() are called, therefore
  //  we will set our own limits as an approximation, so
  //  that any function that calls axis(Chart::XAxis).minimum()/maximum()
  //  before the range is computer by Wt, will at least approximately get
  //  the correct range...  you probably want to call setAutoYAxisRange()
  //  as well.
  //  Doing this means it's neccasary to re-call setAutoXAxisRange()
  //  whenever you load new data
  SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model() );
  if( theModel )
  {
    std::shared_ptr<const Measurement> axisH = theModel->histUsedForXAxis();
    if( axisH )
    {
      float xmin = min( 0.0f, axisH->gamma_energy_min() );
      float xmax = axisH->gamma_energy_max();
      setXAxisRange( xmin, xmax );
    }//if( axisH )
  }else
  {
    axis(Chart::XAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  }//if( theModel is a SpectrumDataModel ) / else
}//void setAutoXAxisRange()


void SpectrumChart::yAxisRangeFromXRange( const double x0, const double x1,
                                          double &y0, double &y1 )
{
  y0 = y1 = 0.0;
  SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model() );
  if( !theModel )
    throw runtime_error( "SpectrumChart::setAutoXAxisRange(): warning, you "
                         "should only use the SpectrumChart class with a "
                         "SpectrumChart!" );
  
  theModel->yRangeInXRange( x0, x1, y0, y1 );
    
  if( axis(Chart::OrdinateAxis).scale() == Chart::LogScale )
  {
    //Specify the (approx) fraction of the chart that the scale should extend
    //  past where the data where hit.
    const double yfractop = 0.05, yfracbottom = 0.025;
    
    const double y0Intitial = ((y0<=0.0) ? 0.1 : y0);
    double y1Intitial = ((y1<=0.0) ? 1.0 : y1);
    y1Intitial = ((y1Intitial<=y0Intitial) ? 1.1*y0Intitial : y1Intitial);
    
    
    const double logY0 = log10(y0Intitial);
    const double logY1 = log10(y1Intitial);
    
    const double logLowerY = ((y0<=0.0) ? -1.0 : (logY0 - yfracbottom*(logY1-logY0)));
    const double logUpperY = logY1 + yfractop*(logY1-logY0);
    
    const double ylower = std::pow( 10.0, logLowerY );
    const double yupper = std::pow( 10.0, logUpperY );
    
    y0 = ((y0<=0.0) ? 0.1 : ylower);
    y1 = ((y1<=0.0) ? 1.0 : yupper);
  }else //if( axis(Chart::YAxis).scale() == Chart::LinearScale )
  {
    y0 = ((y0 <= 0.0) ? 1.1*y0 : 0.9*y0);
    y1 = 1.1*y1;
  }//if( Y is LogScale ) / else
}//void yAxisRangeFromXRange(...)


void SpectrumChart::setAutoYAxisRange()
{
  double y0=0.0, y1=0.0;
  const double x0 = axis(Chart::XAxis).minimum();
  const double x1 = axis(Chart::XAxis).maximum();

  yAxisRangeFromXRange( x0, x1, y0, y1 );

  if( y0 == y1 )
    axis(Chart::OrdinateAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  else
    setYAxisRange( y0, y1 );
}//void setAutoYAxisRange()


void SpectrumChart::paintTextInMiddleOfChart( WPainter &painter ) const
{
  if( m_textInMiddleOfChart.empty() )
    return;
  
  painter.save();
  const WFont oldFont = painter.font();
  WFont txtFont = painter.font();
  //    txtFont.setSize( WFont::XLarge );
  txtFont.setSize( WLength(32) );
  painter.setFont( txtFont );
  const int left_padding = plotAreaPadding(Left);
  const int right_padding = plotAreaPadding(Right);
  const int top_padding = plotAreaPadding(Top);
  const int width = painter.window().width() - left_padding - right_padding;
  const int height = painter.window().height()-plotAreaPadding(Top)-plotAreaPadding(Bottom);
  
  try
  {
    //TextWordWrap not supported by all painters
    WRectF textArea( left_padding, top_padding, width, height );
    painter.drawText( textArea, AlignCenter|AlignMiddle,
                     TextWordWrap, m_textInMiddleOfChart );
  }catch(...)
  {
    const double lineHeight = 40.0;
    const double fontlen = txtFont.sizeLength().toPixels();
    const size_t ncharacters = static_cast<size_t>( 2.0 * width / fontlen );
    string text = m_textInMiddleOfChart.narrow();
    
    if( text.size() < ncharacters )
    {
      const double fromTop = top_padding + 0.5*height - 0.5*lineHeight;
      WRectF textArea( left_padding, fromTop, width, lineHeight );
      painter.drawText( textArea, AlignCenter|AlignMiddle, m_textInMiddleOfChart );
    }else
    {
      text = m_textInMiddleOfChart.toUTF8();
      vector<string> words, lines;
      SpecUtils::split( words, text, " \t" );
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
          lines.push_back( line );
          line = word;
        }//if( theres room for this word ) / else
      }//for( const string &word : words )
      
      if( !line.empty() )
        lines.push_back( line );
      
      const double initialHeight = top_padding + 0.5*height
      - 0.5*lineHeight*lines.size();
      for( size_t i = 0; i < lines.size(); ++i )
      {
        const double fromTop = initialHeight + lineHeight*i;
        //          WRectF textArea( left_padding + 20.0, fromTop, width, lineHeight );
        //          painter.drawText( textArea, AlignLeft|AlignMiddle, lines[i] );
        WRectF textArea( left_padding, fromTop, width, lineHeight );
        painter.drawText( textArea, AlignCenter|AlignMiddle, lines[i] );
      }//for( size_t i = 0; i < lines.size(); ++i )
    }//if( text.size() < ncharacters ) / else
  }//try / catch
  
  painter.restore();
  painter.setFont( oldFont );
}//void paintTextInMiddleOfChart(...)



void SpectrumChart::paint( WPainter &painter, const WRectF &rectangle ) const
{
  //This function is called after (or maybye from within) SpectrumChart::paintEvent(...)
  Chart::WCartesianChart::paint( painter, rectangle );
  
  paintTextInMiddleOfChart( painter );
  
  //Check if data is both positive and negative, and if so, draw a grey line at
  //  y==0; this is mostly for looking at second derivatives in comparison to
  //  to the spectrum (a development feature).
  if( axis(Chart::YAxis).scale() == Chart::LinearScale )
  {
    const double ymin = axis(Chart::YAxis).minimum();
    const double ymax = axis(Chart::YAxis).maximum();
    if( ymax > 0.0 && ymin < -0.05*ymax )
    {
      const WPen oldPen = painter.pen();
      
      const double x_min = axis(Chart::XAxis).minimum();
      const double x_max = axis(Chart::XAxis).maximum();
      const WPointF leftpoint = mapToDevice( x_min, 0.0 );
      const WPointF rightpoint = mapToDevice( x_max, 0.0 );
      
      painter.save();
      painter.setPen( WPen(Wt::gray) );
      painter.drawLine( leftpoint, rightpoint );
      painter.restore();
      painter.setPen( oldPen );
    }//if( ymax > 0.0 && ymin < -0.01*ymax )
  }//if( axis(Chart::YAxis).scale() == Chart::LinearScale )
  
  
  paintHighlightRegions( painter );

  paintPeaks( painter );
  
  paintOnChartLegend( painter );
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  renderReferncePhotoPeakLines( painter );
#endif
  
  loadPeakInfoToClient();
  loadXAndYRangeEquationToClient();
  
  if( m_overlayCanvas )
  {
#if( IOS || ANDROID )
    wApp->doJavaScript( "Wt.WT.OverlayTouchChange(" + m_overlayCanvas->jsRef() + ",null);" ); 
#else
    wApp->doJavaScript( "Wt.WT.DrawGammaLines('c" + m_overlayCanvas->id() + "',true);" );
#endif
  }
}//paint( ... )


void SpectrumChart::paintHighlightRegions( Wt::WPainter &painter ) const
{
  const WFont  oldFont  = painter.font();
  const WPen   oldPen   = painter.pen();
  const WBrush oldBrush = painter.brush();
  
  for( size_t i = 0; i < m_highlights.size(); ++i )
  {
    const HighlightRegion &reg = m_highlights[i];
    const double y_min = axis(Chart::YAxis).minimum();
    const double y_max = axis(Chart::YAxis).maximum();
    
    const WPointF topLeft = mapToDevice( reg.lowerx, y_max );
    WPointF bottomRight = mapToDevice( reg.upperx, y_min );
    
    if( fabs(bottomRight.x()-topLeft.x()) < 2.0 )
    {
      const double mult = (reg.lowerx>reg.upperx) ? -1.0 : 1.0;
      bottomRight.setX( bottomRight.x() + mult*2.0 );
    }
    
    painter.fillRect( WRectF(topLeft, bottomRight), WBrush(reg.color) );
  }//for( size_t i = 0; i < m_highlights.size(); ++i )
  
  for( const auto &reg : m_occupancy_regions )
  {
    const double y_min_counts = axis(Chart::YAxis).minimum();
    WPointF bottomLeft = mapToDevice( reg.lowerx, y_min_counts );
    WPointF bottomRight = mapToDevice( reg.upperx, y_min_counts );
    
    const double bottomPx = chartArea().bottom() + 0.5;// + axis.margin();
    const double topPx = chartArea().top();// + axis.margin();
    
    //If we wanted to shade in the x-axis
    //bottomLeft.setY( bottomPx );
    //bottomRight.setY( bottomPx + 5 );
    //painter.fillRect( WRectF(bottomLeft, bottomRight), WBrush(reg.color) );
    
    painter.setPen( WPen(reg.color) );
    painter.setBrush( WBrush(reg.color) );
    
    const float leftPx = bottomLeft.x();
    const float rightPx = bottomRight.x();
    painter.drawLine( WPointF(leftPx,bottomPx), WPointF(leftPx,topPx) );
    painter.drawLine( WPointF(rightPx,bottomPx), WPointF(rightPx,topPx) );
    
    if( fabs(rightPx - leftPx) > 50 )
    {
      WFont font( WFont::Default );
      font.setSize( 8 );
      painter.setFont( font );
      painter.drawText(leftPx+2, topPx, 30, 10.0, Wt::AlignmentFlag::AlignLeft, "occ. start");
      painter.drawText(rightPx-2-30, topPx, 30, 10.0, Wt::AlignmentFlag::AlignRight,  "occ. end");
    }
  }//for( const auto region : m_underlines )
  
  painter.setPen( oldPen );
  painter.setFont( oldFont );
  painter.setBrush( oldBrush );
}//void paintHighlightRegions( Wt::WPainter& painter ) const


void SpectrumChart::loadPeakInfoToClient() const
{
  //right now peakModelDataChanged(...) assumes loadPeakInfoToClient()
  //  only loads nuclide symbol and energy to client, and
  //  not kUseForShieldingSourceFit, kCandidateIsotopes, kUseForCalibration;
  //  if this changes in the future, make sure to update
  //  peakModelDataChanged(...)
  
  if( !m_overlayCanvas || !m_peakModel )
    return;

#if( WT_VERSION >= 0x3030800 )
  auto absModel = model();
#else
  const WAbstractItemModel *absModel = model();
#endif
  const SpectrumDataModel *th1Model = dynamic_cast<const SpectrumDataModel *>( absModel );
  
  if( !th1Model )
  {
    cerr << "SpectrumChart::loadPeakInfoToClient(): no SpectrumDataModel" << endl;
    return;
  }

  std::shared_ptr<const Measurement> data;
  int dataCol = -1, backgroundCol = -1;

  if( th1Model )
  {
    data = th1Model->getData();

    dataCol = th1Model->dataColumn();
    backgroundCol = th1Model->backgroundColumn();
  }//if( th1Model )


  stringstream js;
  js << ""
     << "" "{"
     << ""   "try{"
     << ""      "$('#c" << m_overlayCanvas->id() << "').data('peaks',"
     << ""        "[";

  const size_t npeaks = m_peakModel->npeaks();
  for( size_t peakn = 0; peakn < npeaks; ++peakn )
  {
    const PeakDef &peak = m_peakModel->peak(peakn);
    if( peakn )
      js << ",";

    try
    {
      double ymin = 0.0;
      
      double lowx(0.0), upperx(0.0);
      findROIEnergyLimits( lowx, upperx, peak, data );
      
      const int peakrow = th1Model->findRow( peak.mean() );
//      const int lowerrow = th1Model->findRow( lowx );
//      const int upperrow = th1Model->findRow( upperx );

      double ymax = th1Model->data( peakrow, dataCol );

      if( data && th1Model->backgroundSubtract() && (backgroundCol >= 0) )
      {
        ymin = th1Model->data( peakrow, backgroundCol );
        ymax -= ymin;
      }//if( data && th1Model )

      double contarea = 0;
      switch( peak.continuum()->type() )
      {
        case PeakContinuum::NoOffset:     case PeakContinuum::Constant:
        case PeakContinuum::Linear:       case PeakContinuum::Quadratic:
        case PeakContinuum::Cubic:        case PeakContinuum::External:
          contarea = peak.offset_integral( lowx, upperx, data );
          break;
          
        case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
        case PeakContinuum::BiLinearStep:
          if( data )
            contarea = peak.offset_integral( lowx, upperx, data );
          break;
      }//switch( peak.continuum()->type() )
      

      double lowe = 0.0, highe = 0.0;
      findROIEnergyLimits( lowe, highe, peak, data );
      
      const WPointF lower = mapToDevice( lowe, ymin );
      const WPointF upper = mapToDevice( highe, ymax );
      const WPointF meanpoint = mapToDevice( peak.mean(), ymax );

//      cout << "testpeaks.push_back( EnergyIntensityPair( " << fixed << setprecision(2) << peak.mean
//           << ", " << fixed << setprecision(1) << peak.amplitude << " );" << endl;

      if( IsInf(peak.mean()) || IsNan(peak.mean()) )
        continue;

      js << "{mean:"   << fixed << setprecision(2) << peak.mean();
      const double ampuncert = peak.amplitudeUncert();
      
      switch ( peak.type() )
      {
        case PeakDef::GaussianDefined:
          js << ",sigma:" << fixed << setprecision(2) << peak.sigma()
             << ",amp:"   << fixed << setprecision(1) << peak.peakArea();
          if( ampuncert > 0.0 && !IsNan( ampuncert ) && !IsInf( ampuncert ) )
            js << ",ampuncert:" << fixed << setprecision(2) << ampuncert;
        break;
          
        default:
          js << ",sigma:'na',amp:'na'";
        break;
      }//switch ( peak.type() )


      float energy = 0.0f;
      try
      {
        energy = peak.gammaParticleEnergy();
      }catch(...)
      {
      }
      
      if( peak.parentNuclide() && peak.decayParticle() )
      {
        js << ",nuclide:'" << peak.parentNuclide()->symbol << " (";
        
        switch( peak.sourceGammaType() )
        {
          case PeakDef::NormalGamma:
          case PeakDef::AnnihilationGamma:
          {
            const SandiaDecay::Transition *trans = peak.nuclearTransition();
            if( trans && trans->parent && trans->parent != peak.parentNuclide() )
              js << trans->parent->symbol << ", ";
            break;
          }
            
          case PeakDef::SingleEscapeGamma:
          {
            energy += 510.99891f;
            js << "S.E. from ";
            break;
          }
            
          case PeakDef::DoubleEscapeGamma:
          {
            energy += 2.0f*510.99891f;
            js << "D.E. from ";
            break;
          }
            
          case PeakDef::XrayGamma:
          {
            const SandiaDecay::Transition *trans = peak.nuclearTransition();
            if( trans && trans->parent && trans->parent != peak.parentNuclide() )
              js << trans->parent->symbol << ", ";
            js << "xray ";
            break;
          }
        }//switch( peak.sourceGammaType() )
        
        js << std::setprecision(2) << energy << " keV)'";
      }else if( peak.parentNuclide() && energy>0.0f)
      {
        js << ",nuclide:'" << peak.parentNuclide()->symbol
           << " (" << std::setprecision(2) << energy << " keV ann.)'";
      }else if( peak.xrayElement() )
      {
        js << ",nuclide:'" << peak.xrayElement()->symbol << " x-ray ("
           << std::setprecision(2) << energy << " keV)'";
      }else if( peak.reaction() )
      {
        js << ",nuclide:'" << peak.reaction()->name() << " ("
        << std::setprecision(2) << energy << " keV)'";
      }//if( parent nuclide ) / else ( xray ) / else if( reaction )

      if( peak.chi2Defined() )
        js << ",chi2:'" << std::setprecision(2) << peak.chi2dof() << "'";
      
      js <<   ",xminp:" << static_cast<int>(floor(lower.x()+0.5))
         <<   ",xmaxp:" << static_cast<int>(floor(upper.x()+0.5))
         <<   ",yminp:" << static_cast<int>(floor(lower.y()+0.5))
         <<   ",ymaxp:" << static_cast<int>(floor(upper.y()+0.5))
         <<   ",meanp:" << static_cast<int>(floor(meanpoint.x()+0.5))
         <<   ",cont:" << fixed << setprecision(1) << contarea
         << "}";
    }catch(...)
    {
      //
    }//try / catch
  }//for( loop over peaks )

  js << ""        "]);"
     << ""     "}catch(e){}"
     << ""   "}"
        ;

  wApp->doJavaScript( js.str() );
}//void loadPeakInfoToClient()


std::string SpectrumChart::pixelToCoordinateMappingJs( const Wt::Chart::Axis axisType ) const
{
  const Chart::WAxis &yaxis = axis(axisType);
  const Chart::AxisScale yscale = yaxis.scale();
  const bool logy = (yscale==Chart::LogScale); //CategoryScale LinearScale LogScale DateScale DateTimeScale

  if( logy )
  {
    const double ymin = max(0.1,yaxis.minimum());
    const double ymax = yaxis.maximum();

    WPointF lowerLog, upperLog;
    if( axisType == Chart::XAxis )
    {
      lowerLog = mapToDevice( ymin, 0.0, axisType );
      upperLog = mapToDevice( ymax, 0.0, axisType );
    }else
    {
      lowerLog = mapToDevice( 0.0, ymin, axisType );
      upperLog = mapToDevice( 0.0, ymax, axisType );
    }//if( axisType == Chart::XAxis ) / else

/*
    double h, y;
    const double lowerYCoord = lowerLog.y();
    const double upperYCoord = upperLog.y();
    const double midy_pixel = (lowerLog.y() - upperLog.y())/2.0;
    const WPointF midLog = mapFromDevice( WPointF(0.0, midy_pixel) );
    h = midy_pixel;
    y = std::exp(std::log(ymin) + h * (std::log(ymax) - std::log(ymin)) / (lowerLog.y() - upperLog.y()) );
    WPointF answer = mapFromDevice( WPointF(0.0, lowerYCoord-h) );
    cerr << "trying h=" << h << ", y=" << y << ", ymax=" << ymax << ", answer.y()=" << answer.y() << endl;
*/
    const double y_pixel_range = lowerLog.y() - upperLog.y();
    const double log_coord_range_ratio = std::log(ymax) - std::log(ymin);
    const double pixel_multiplier = log_coord_range_ratio / y_pixel_range;
    WStringStream js;
    js << "function(y){"
       << "return Math.exp(" << std::log(ymin)
                             << "+(" << lowerLog.y() << "-y)*"
                             << pixel_multiplier
       << ");}";
    return js.str();
  }else
  {
    const WPointF lowerLeft = mapFromDevice( WPointF(0.0,m_heightInPixels), axisType );
    const WPointF upperRight = mapFromDevice( WPointF(m_widthInPixels,0.0), axisType );

    WStringStream js;
    if( axisType == Chart::XAxis )
    {
      const double x0 = lowerLeft.x();
      const double x1 = (upperRight.x() - x0) / m_widthInPixels;
      js << "function(x){return " << x0 << "+" << x1 << "*x;}";
    }else
    {
      const double y0 = lowerLeft.y();
      const double y1 = (upperRight.y()-y0) / m_heightInPixels;
      js << "function(y){return " << y0 << "-" << y1 << "*(y-" << m_heightInPixels << ");}";
    }//if( axisType == Chart::XAxis ) / else

    return js.str();
  }//if( logy ) / else

  return "null"; //wont ever get here
}//std::string pixelToCoordinateMappingJs( Wt::Chart::Axis axis )

//#if(DRAW_GAMMA_LINES_LOG_AND_LIN)
std::string SpectrumChart::countsToPixelMappingJs() const
{
  const Chart::WAxis &yaxis = axis(Chart::YAxis);
  const Chart::AxisScale yscale = yaxis.scale();
  const bool logy = (yscale==Chart::LogScale); //CategoryScale LinearScale LogScale DateScale DateTimeScale
  
  if( logy )
  {
    const double ymin = max(0.1,yaxis.minimum());
    const double ymax = yaxis.maximum();
    
    WPointF lowerLog, upperLog;
    lowerLog = mapToDevice( 0.0, ymin, Chart::YAxis );
    upperLog = mapToDevice( 0.0, ymax, Chart::YAxis );
    
    const double y_pixel_range = lowerLog.y() - upperLog.y();
    const double log_coord_range_ratio = std::log(ymax) - std::log(ymin);
    const double pixel_multiplier = log_coord_range_ratio / y_pixel_range;
    WStringStream js;
    
    js << "function(counts){ return " << lowerLog.y() << "-(Math.log(counts)-("
       << std::log(ymin) << "))/" << pixel_multiplier << ";}";
    return js.str();
  }else
  {
    const WPointF lowerLeft = mapFromDevice( WPointF(0.0,m_heightInPixels), Chart::YAxis );
    const WPointF upperRight = mapFromDevice( WPointF(m_widthInPixels,0.0), Chart::YAxis );
    
    WStringStream js;
    const double y0 = lowerLeft.y();
    const double y1 = (upperRight.y()-y0) / m_heightInPixels;
    const double high = 1.0*m_heightInPixels + y0/y1;
    js << "function(counts){return " << high << "-counts/" << y1 << ";}";
    return js.str();
  }//if( logy ) / else
  
  return "null"; //wont ever get here
}//std::string SpectrumChart::countsToPixelMappingJs() const
//#endif

std::string SpectrumChart::energyToPixelMappingJs() const
{
  const WPointF lowerLeft = mapFromDevice( WPointF(0.0,m_heightInPixels), Wt::Chart::XAxis );
  const WPointF upperRight = mapFromDevice( WPointF(m_widthInPixels,0.0), Wt::Chart::XAxis );

  WStringStream js;
  const double x0 = lowerLeft.x();
  const double x1 = (upperRight.x() - x0) / m_widthInPixels;

  if( x0 >= 0.0 )
    js << "function(x){return (x-" << x0 << ")/"<< x1 << ";}";
  else
    js << "function(x){return (x+" << -x0 << ")/"<< x1 << ";}";

  return js.str();
}//std::string SpectrumChart::energyToPixelMappingJs() const


void SpectrumChart::loadXAndYRangeEquationToClient() const
{
  if( m_xAxisUnits == kUndefinedUnits )
    return;
  
  WStringStream js;
  js << ""
     << "" "{"
     << ""   "try{"
     << ""      "var can = $('#c" << m_overlayCanvas->id() << "');"
     << ""      "can.data('xeqn'," << pixelToCoordinateMappingJs( Wt::Chart::XAxis ) << ");"
     << ""      "can.data('yeqn'," << pixelToCoordinateMappingJs( Wt::Chart::YAxis ) << ");"
     << ""      "can.data('exeqn'," << energyToPixelMappingJs() << ");"
//#if(DRAW_GAMMA_LINES_LOG_AND_LIN)
     << ""      "can.data('cyeqn'," << countsToPixelMappingJs() << ");"
//#endif
        ;
  if( axis(Chart::Y2Axis).isVisible() )
    js << ""    "can.data('y2eqn'," << pixelToCoordinateMappingJs( Wt::Chart::Y2Axis ) << ");";
  else
    js << ""    "can.data('y2eqn',null);";
  js << ""     "}catch(e){}"
     << ""   "}";

  wApp->doJavaScript( js.str() );
}//void loadXAndYRangeEquationToClient() const




void SpectrumChart::setShowRefLineInfoForMouseOver( const bool show )
{
  if( m_overlayCanvas )
    m_overlayCanvas->setShowRefLineInfoForMouseOver( show );
}//void setShowRefLineInfoForMouseOver( const bool show )


void SpectrumChart::paintNonGausPeak( const PeakDef &peak, Wt::WPainter& painter ) const
{
  using boost::any_cast;

  if( peak.type() != PeakDef::DataDefined )
    return;

#if( WT_VERSION >= 0x3030800 )
  auto absModel = model();
#else
  const WAbstractItemModel *absModel = model();
#endif
  const SpectrumDataModel *th1Model = dynamic_cast<const SpectrumDataModel *>( absModel );
  if( !th1Model )
    throw runtime_error( "SpectrumChart::paintNonGausPeak(...): stupid programmer error"
                         ", SpectrumChart may only be powered by a SpectrumDataModel." );

  const std::shared_ptr<const Measurement> data = th1Model->getData();
  if( !data )
  {
    cerr << "SpectrumChart::paint: No data histogram defined" << endl;
    return;
  }//if( !data )

  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  
  const double minDispX = axis(Chart::XAxis).minimum();
  const double maxDispX = axis(Chart::XAxis).maximum();
  const double minDispY = axis(Chart::YAxis).minimum();

  const WPointF minXDev = mapToDevice( minDispX, minDispY, Chart::XAxis );
  const WPointF maxXDev = mapToDevice( maxDispX, minDispY, Chart::XAxis );
  
  const double axisMinX = mapFromDevice( minXDev ).x();
  const double axisMaxX = mapFromDevice( maxXDev ).x();
  
  const float lowerx = static_cast<float>( std::max( peak.lowerX(), axisMinX ) );
  const float upperx = static_cast<float>( std::min( peak.upperX(), axisMaxX ) );
  
  if( lowerx > upperx )
    return;
  
  const int lowerrow = th1Model->findRow( lowerx );
  const int upperrow = th1Model->findRow( upperx );
  
  //First, draw the continuum
  
  WPen pen( peak.lineColor().isDefault() ? m_defaultPeakColor : peak.lineColor() );
  pen.setWidth( 2 );
  painter.setPen( pen );
  switch( continuum->type() )
  {
    case PeakContinuum::NoOffset:
    break;
      
    case PeakContinuum::Constant:     case PeakContinuum::Linear:
    case PeakContinuum::Quadratic:    case PeakContinuum::Cubic:
    case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep: case PeakContinuum::External:
    {
      vector<WPointF> path;
      for( int row = lowerrow; row <= upperrow; ++row )
      {
        const double rowLowX = th1Model->rowLowEdge( row );
        const double rowUpperX = rowLowX + th1Model->rowWidth( row );
        const double lowerx = std::max( rowLowX, axisMinX );
        const double upperx = std::min( rowUpperX, axisMaxX );
        const double contheight = continuum->offset_integral( rowLowX, rowUpperX, data );
      
        const WPointF startPoint = mapToDevice( lowerx, contheight );
        const WPointF endPoint = mapToDevice( upperx, contheight );
      
        if( fabs(startPoint.x() - endPoint.x()) < 2.0 )
        {
          WPointF point( 0.5*(startPoint.x() + endPoint.x()), endPoint.y() );
          path.push_back( point );
        }else
        {
          path.push_back( startPoint );
          path.push_back( endPoint );
        }
      }//for( int row = lowerrow; row <= upperrow; ++row )
      
      painter.drawLines( path );
      break;
    }//case external or polynomial continuum
  }//switch( continuum->type() )
  
  //Now draw the areas of data above the continuum
  WPointF startpoint;
  std::unique_ptr<WPainterPath> drawpath;
  vector<WPointF> continuumvalues;
  
  for( int row = lowerrow; row <= upperrow; ++row )
  {
    const double rowLowX = th1Model->rowLowEdge( row );
    const double rowUpperX = rowLowX + th1Model->rowWidth( row );
    
    const double lowerx = std::max( rowLowX, axisMinX );
    const double upperx = std::min( rowUpperX, axisMaxX );
    const boost::any dataheightany = th1Model->displayBinValue( row, SpectrumDataModel::DATA_COLUMN );
    const double contheight = continuum->offset_integral( rowLowX, rowUpperX, data );
    const double dataheight = asNumber( dataheightany );
    
    if( dataheightany.empty() )
      continue;
    
    if( dataheight > contheight )
    {
      WPointF startPoint = mapToDevice( lowerx, contheight );
      WPointF endPoint = mapToDevice( upperx, contheight );
      
      const bool isnarrow = (fabs(startPoint.x() - endPoint.x()) < 2.0);
      
      if( isnarrow )
        startPoint = WPointF( 0.5*(startPoint.x() + endPoint.x()), endPoint.y() );
      continuumvalues.push_back( startPoint );
      if( !isnarrow )
        continuumvalues.push_back( endPoint );
      
      startPoint = mapToDevice( lowerx, dataheight );
      endPoint = mapToDevice( upperx, dataheight );
     
      if( isnarrow )
        startPoint = WPointF( 0.5*(startPoint.x() + endPoint.x()), endPoint.y() );
      
      if( !drawpath )
      {
        startpoint = startPoint;
        drawpath.reset( new WPainterPath() );
      }
      
      drawpath->lineTo( startPoint );
      if( !isnarrow )
        drawpath->lineTo( endPoint );
    }//if( dataheight > contheight )
    
    if( (contheight >= dataheight) || (row == upperrow) )
    {
      if( !!drawpath )
      {
        for( auto iter = continuumvalues.rbegin(); iter != continuumvalues.rend(); ++iter )
          drawpath->lineTo( *iter );
        
        drawpath->lineTo( startpoint );
        
        painter.fillPath( *drawpath, WBrush( (peak.lineColor().isDefault() ? m_defaultPeakColor : peak.lineColor()) ) );
        //painter.fillPath( *drawpath, WBrush(m_defaultPeakColor) );
        
        drawpath.reset();
        continuumvalues.clear();
      }//if( !!drawpath )
    }//if( (contheight >= dataheight) || (row == upperrow) )
  }//for( int row = lowerrow; row <= upperrow; ++row )
}//void paintNonGausPeak( const Peak &peak )


double SpectrumChart::peakBackgroundVal( const int row, const PeakDef &peak,
                    const SpectrumDataModel *th1Model,
                    std::shared_ptr<const PeakDef> prevPeak,
                    std::shared_ptr<const PeakDef> nextPeak )
{
  const double bin_x_min = th1Model->rowLowEdge( row );
  const double bin_x_max = bin_x_min + th1Model->rowWidth( row );
  
  
  double yval = 0.0;
  
  switch( peak.continuum()->type() )
  {
    case PeakContinuum::NoOffset:     case PeakContinuum::Constant:
    case PeakContinuum::Linear:       case PeakContinuum::Quadratic:
    case PeakContinuum::Cubic:        case PeakContinuum::External:
      yval = peak.offset_integral( bin_x_min, bin_x_max, nullptr );
      break;
    
    case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    {
      shared_ptr<const SpecUtils::Measurement> data = th1Model->getData();
      if( data )
        yval = peak.offset_integral( bin_x_min, bin_x_max, data );
      break;
    }
  }//switch( peak.continuum()->type() )
  
  if( th1Model->backgroundSubtract() && (th1Model->backgroundColumn() >= 0) )
    yval -= th1Model->data( row, th1Model->backgroundColumn() );

  return yval;
  
//  double distToPrev = std::numeric_limits<double>::infinity();
//  double distToNext = std::numeric_limits<double>::infinity();
//  const double distToCurrent = fabs( peak.mean() - bin_x_min );
//  if( prevPeak )
//    distToPrev = fabs( prevPeak->mean() - bin_x_min );
//  if( nextPeak )
//    distToNext = fabs( nextPeak->mean() - bin_x_min );
//  
//  if( (distToCurrent<=distToPrev) && (distToCurrent<=distToNext) )
//  {
//    if( peak.offsetType() >= PeakDef::Constant )
//      ystart += peak.offset_integral( bin_x_min, bin_x_max );
//  }else if( prevPeak && (distToPrev<=distToNext) )
//  {
//    if( peak.offsetType() >= PeakDef::Constant )
//      ystart += prevPeak->offset_integral( bin_x_min, bin_x_max );
//  }else if( nextPeak )  //the if( nextPeak ) shouldnt be necessary
//  {
//    ystart += nextPeak->offset_integral( bin_x_min, bin_x_max );
//  }

//  return ystart;
}//double peakBackgroundVal(...)


double SpectrumChart::peakYVal( const int row, const PeakDef &peak,
                                const SpectrumDataModel *th1Model,
                               std::shared_ptr<const PeakDef> prevPeak,
                               std::shared_ptr<const PeakDef> nextPeak )
{
  double bin_x_min = th1Model->rowLowEdge( row );
  double bin_width = th1Model->rowWidth( row );
  double bin_x_max = bin_x_min + bin_width;

  double ystart = 0.0;
  ystart += peak.gauss_integral( bin_x_min, bin_x_max );
  ystart += peakBackgroundVal( row, peak, th1Model, prevPeak, nextPeak );
  
  return ystart;
}//double peakYVal(...)


void SpectrumChart::visibleRange( double &axisMinX, double &axisMaxX,
                                  double &axisMinY, double &axisMaxY ) const
{
  axisMinX = axis(Chart::XAxis).minimum();
  axisMaxX = axis(Chart::XAxis).maximum();
  axisMinY = axis(Chart::YAxis).minimum();
  axisMaxY = axis(Chart::YAxis).maximum();
  
  
  //Wt::Chart draws about 5 pixels to the left of axisMinX, so lets
  //  compensate for this.
  //Note that for newer versions of Wt, a function axisPadding() may be
  //  available rather than hard-coding in this 5.0
  // XXX - below might not be completely working
  WPointF minPoint = mapToDevice( axisMinX, axisMinY, Chart::XAxis );
  minPoint.setX( minPoint.x() - axisPadding() );
  minPoint = mapFromDevice( minPoint, Chart::XAxis );
  axisMinX = minPoint.x();
  
  WPointF maxPoint = mapToDevice( axisMaxX, axisMaxY, Chart::XAxis );
  maxPoint.setX( maxPoint.x() + axisPadding() );
  maxPoint = mapFromDevice( maxPoint, Chart::XAxis );
  axisMaxX = maxPoint.x();
}//void visibleXRange(...) const


float SpectrumChart::xUnitsPerPixel() const
{
  const double yValMin = axis(Chart::YAxis).minimum();
  const double xValMin = axis(Chart::XAxis).minimum();
  const double xValMax = axis(Chart::XAxis).maximum();
  
  const WPointF minpx = mapToDevice( xValMin, yValMin, Chart::XAxis );
  const WPointF maxpx = mapToDevice( xValMax, yValMin, Chart::XAxis );
  
  return static_cast<float>( (xValMax - xValMin) / (maxpx.x() - minpx.x()) );
}//float xUnitsPerPixel() const


bool SpectrumChart::gausPeakDrawRange( int &minrow, int &maxrow,
                                       const PeakDef &peak,
                                       const double axisMinX,
                                       const double axisMaxX ) const
{
  if( !peak.gausPeak() )
    return false;
  
  const SpectrumDataModel *th1Model
                           = dynamic_cast<const SpectrumDataModel *>( model() );
  assert( th1Model );

  double minx = peak.lowerX();
  double maxx = peak.upperX();
  
  std::shared_ptr<const Measurement> dataH = th1Model->getData();
  if( dataH )
  {
    size_t lower_channel, upper_channel;
    estimatePeakFitRange( peak, dataH, lower_channel, upper_channel );
    minx = dataH->gamma_channel_lower(lower_channel) + 0.001;
    maxx = dataH->gamma_channel_upper(upper_channel);
  }//if( dataH )

  
  try
  {
    minrow = th1Model->findRow( minx );
    maxrow = th1Model->findRow( maxx );
    if( maxrow >= th1Model->rowCount() )
      maxrow = th1Model->rowCount() - 1;
    
    if( minx>maxx )
      return false;

    double bin_x_min = th1Model->rowLowEdge( minrow );
    double bin_width = th1Model->rowWidth( minrow );
    double bin_x_max = bin_x_min + bin_width;
  
    //Lets not draw the part of the peak than might less than the minimum X of
    //  the x-axis
    while( (bin_x_max < axisMinX) && (minrow < maxrow) )
    {
      ++minrow;
      bin_x_min = th1Model->rowLowEdge( minrow );
      bin_width = th1Model->rowWidth( minrow );
      bin_x_max = bin_x_min + bin_width;
    }//while( (xstart < axisMinX) && (minrow < maxrow) )
    
    bin_x_min = th1Model->rowLowEdge( maxrow );
    while( (bin_x_min > axisMaxX) && (maxrow > minrow) )
    {
      --maxrow;
      bin_x_min = th1Model->rowLowEdge( maxrow );
    }//while( (xstart < axisMinX) && (minrow < maxrow) )
    
    return (minrow < maxrow);
  }catch(...)
  {
    cerr << "Shouldnt have caught errror in SpectrumChart::gausPeakDrawRange(...)" << endl;
    return false;
  }
    
  return true;
}//bool gausPeakDrawRange(...)





void SpectrumChart::gausPeakBinValue( double &bin_center, double &y,
                       const PeakDef &peak, const int row,
                       const SpectrumDataModel *th1Model,
                       const double axisMinX, const double axisMaxX,
                       const double axisMinY, const double axisMaxY,
                       std::shared_ptr<const PeakDef> prevPeak,
                       std::shared_ptr<const PeakDef> nextPeak ) const 
{
  y = 0.0;
  
  bin_center = th1Model->rowCenter( row );
    
  y = peakYVal( row, peak, th1Model, prevPeak, nextPeak );
  
  if( bin_center < axisMinX && row>1 )
  {
    double prevMid = th1Model->rowCenter( row - 1 );
    const double prevYVal = peakYVal( row-1, peak,
                                      th1Model, prevPeak, nextPeak );
    const double dist = (axisMinX - prevMid) / (bin_center - prevMid);
    y = prevYVal + dist * (y - prevYVal);
    bin_center = axisMinX;
  }//if( xstart < axisMinX )
    
  if( bin_center > axisMaxX && row>1 )
  {
    double prevMid = th1Model->rowCenter( row - 1 );
    const double prevYVal = peakYVal( row-1, peak, th1Model,
                                      prevPeak, nextPeak );
    const double dist = (axisMaxX - prevMid) / (bin_center - prevMid);
    y = prevYVal + dist * (y - prevYVal);
    bin_center = axisMaxX;
  }//if( xstart < axisMinX )
    
  y = std::max( y, axisMinY );
  y = std::min( y, axisMaxY );
  
  bin_center = std::max( bin_center, axisMinX );
  bin_center = std::min( bin_center, axisMaxX );
}//gausPeakBinValue



void SpectrumChart::drawIndependantGausPeak( const PeakDef &peak,
                                         Wt::WPainter& painter,
                                         const bool doFill,
                                         const double xstart,
                                         const double xend ) const
{
  double axisMinX, axisMaxX, axisMinY, axisMaxY;
  visibleRange( axisMinX, axisMaxX, axisMinY, axisMaxY );
  
  if( xstart != xend )
  {
    if( !IsInf(xstart) && !IsNan(xstart) && xstart>axisMinX )
      axisMinX = xstart;
    if( !IsInf(xend) && !IsNan(xend) && xend<axisMaxX )
      axisMaxX = xend;
  }//if( xstart != xend )
  
  const SpectrumDataModel *th1Model
                          = dynamic_cast<const SpectrumDataModel *>( model() );
  assert( th1Model );
  
  
  std::shared_ptr<const PeakDef> prevPeak, nextPeak;
  
  int minrow, maxrow;
  bool visible = gausPeakDrawRange( minrow, maxrow, peak, axisMinX, axisMaxX );
  
  if( !visible )
    return;
  
  //Decide if we should draw peaks to the center of the bins on either end of the ROI, or
  //  if we should explicitly go to the min/max of the left/right-most bins.
  //  Drawing to the extremes is probably more correct, but for some reason it
  //  still has some very minor artifacts
#define DRAW_PEAKS_TO_BIN_EXTREMES 0
  
  
  WPointF start_point(0.0, 0.0), peakMaxInPx( 999999.9, 999999.9 );

  try
  {
    double bin_center, y;
    gausPeakBinValue( bin_center, y, peak, minrow, th1Model,
                     axisMinX, axisMaxX, axisMinY, axisMaxY,
                     prevPeak, nextPeak );
#if( DRAW_PEAKS_TO_BIN_EXTREMES )
    double lower_y = y;
    if( minrow > 1 )
      gausPeakBinValue( bin_center, lower_y, peak, minrow-1, th1Model,
                       axisMinX, axisMaxX, axisMinY, axisMaxY,
                       prevPeak, nextPeak );
    const double lowx = std::max( axisMinX, std::min(axisMaxX,th1Model->rowLowEdge(minrow)) );
    lower_y = std::max( axisMinY, std::min(axisMaxY,0.5*(y+lower_y)) );
    start_point = mapToDevice( lowx, lower_y );
#else
    start_point = mapToDevice( bin_center, y );
#endif
  }catch( std::exception &e )
  {
    cerr << "Caught exception: SpectrumChart::drawIndependantGausPeak(...)" << "\n\t" << e.what() << endl;
  }
  
  
  
  WPainterPath path( start_point );
  
  for( int row = minrow; row <= maxrow; ++row )
  {
    try
    {
      double bin_center, y;
      gausPeakBinValue( bin_center, y, peak, row, th1Model,
                       axisMinX, axisMaxX, axisMinY, axisMaxY,
                       prevPeak, nextPeak );
      
      const WPointF pointInPx = mapToDevice( bin_center, y );
      path.lineTo( pointInPx );
      
#if( DRAW_PEAKS_TO_BIN_EXTREMES )
      //For right-most point, draw a line to the very right edge of bin
      if( row == maxrow )
      {
        double upper_y = y;
        if( (row+1) < th1Model->rowCount() )
          gausPeakBinValue( bin_center, upper_y, peak, row+1, th1Model,
                            axisMinX, axisMaxX, axisMinY, axisMaxY,
                             prevPeak, nextPeak );
        const double bin_x_min = th1Model->rowLowEdge( row );
        const double bin_width = th1Model->rowWidth( row );
        const double bin_x_max = bin_x_min + bin_width;
        const double upperx = std::max( axisMinX, std::min(axisMaxX,bin_x_max) );
        upper_y = std::max( axisMinY, std::min(axisMaxY,0.5*(y+upper_y)) );
        
        path.lineTo( mapToDevice(upperx, upper_y) );
      }//if( left most bin ) else if (right most bin )
#endif
      
      if( pointInPx.y() < peakMaxInPx.y() )
        peakMaxInPx = pointInPx;
    }catch(...){}
  }//for( int row = minrow; row <= maxrow; ++row )
  
  
  for( int row = maxrow; row >= minrow; --row )
  {
    try
    {
      const double bin_x_min = th1Model->rowLowEdge( row );
      const double bin_width = th1Model->rowWidth( row );
      const double bin_center = std::max(axisMinX, std::min(axisMaxX, bin_x_min+0.5*bin_width) );
      
      double val = peakBackgroundVal( row, peak, th1Model, prevPeak, nextPeak );
      val = std::min( axisMaxY, std::max(axisMinY, val) );
      
#if( DRAW_PEAKS_TO_BIN_EXTREMES )
      //For right most bin, need to extend to right side of bin
      if( row == maxrow )
      {
        try
        {
          double upper_y = val;
          if( (row+1) < th1Model->rowCount() )
            upper_y = peakBackgroundVal( row+1, peak, th1Model, prevPeak, nextPeak );
          const double upperx = std::max( axisMinX, std::min(axisMaxX,bin_x_min + bin_width) );
          upper_y = std::max( axisMinY, std::min(axisMaxY,0.5*(val+upper_y)) );
          path.lineTo( mapToDevice(upperx, upper_y) );
        }catch(...){}
      }//if( right most bin )
#endif
      
      const WPointF pointInPx = mapToDevice( bin_center, val );
      path.lineTo( pointInPx );
      
      if( pointInPx.y() < peakMaxInPx.y() )
        peakMaxInPx = pointInPx;
    
#if( DRAW_PEAKS_TO_BIN_EXTREMES )
      //For left-most bins,  we actually need to start at very left edge of bin,
      //  to match what we did in the loop above
      if( row == minrow )
      {
        try
        {
          double lower_y = val;
          if( row > 1 )
            lower_y = peakBackgroundVal( row-1, peak, th1Model, prevPeak, nextPeak );
          const double lowx = std::max( axisMinX, std::min(axisMaxX,bin_x_min) );
          lower_y = std::max( axisMinY, std::min(axisMaxY,0.5*(val+lower_y)) );
          path.lineTo( mapToDevice(lowx, lower_y) );
        }catch(...){}
      }//if( left most bin )
#endif
    }catch(...){}
  }//for( int row = maxrow; row >= minrow; --row )
  
  path.lineTo( start_point );
  
  const WPen oldPen = painter.pen();
  const WBrush oldBrush = painter.brush();
  
  WColor color = peak.lineColor().isDefault() ? m_defaultPeakColor : peak.lineColor();
  
  if( doFill )
  {
    color.setRgb( color.red(), color.green(), color.blue(), 155 );
    painter.setBrush( WBrush(color) );
    painter.setPen( WPen(WColor(color.red(),color.green(),color.blue())) );
    painter.drawPath( path );
    //painter.fillPath( path, WBrush(color) );
  }else
  {
    painter.setPen( WPen(color) );
    painter.drawPath( path );
  }
  
  painter.setPen( oldPen );
  painter.setBrush( oldBrush );
  
  drawPeakText( peak, painter, peakMaxInPx );
}//void drawIndependantGausPeak(...)


void SpectrumChart::drawPeakText( const PeakDef &peak, Wt::WPainter &painter,
                                  Wt::WPointF peakTip ) const
{
  int numLabels = 0, labelnum = 0;
  for( PeakLabels i = PeakLabels(0); i < kNumPeakLabels; i = PeakLabels(i+1) )
  {
    if( !m_peakLabelsToShow[i] )
      continue;
    
    switch ( i )
    {
      case kShowPeakUserLabel: numLabels += !peak.userLabel().empty(); break;
      case kShowPeakEnergyLabel: ++numLabels; break;
      case kShowPeakNuclideLabel:
        numLabels += (peak.parentNuclide() || peak.reaction() || peak.xrayElement());
      break;
      case kShowPeakNuclideEnergies:
      case kNumPeakLabels:
      break;
    }//switch ( i )
  }//for( PeakLabels i = PeakLabels(0); i < kNumPeakLabels; i = PeakLabels(i+1) )
  
  if( numLabels )
  {
    for( PeakLabels i = PeakLabels(0); i < kNumPeakLabels; i=PeakLabels(i+1) )
    {
      if( !m_peakLabelsToShow[i] || i==kShowPeakNuclideEnergies )
        continue;
   
      const int red   = (m_peakLabelsColors[i] & 0xFF);
      const int green = ((m_peakLabelsColors[i]>>8) & 0xFF);
      const int blue  = ((m_peakLabelsColors[i]>>16) & 0xFF);
      const int alpha = ((m_peakLabelsColors[i]>>24) & 0xFF);

      WString txt;
      char buffer[32];
      switch ( i )
      {
        case kShowPeakUserLabel:
          txt = WString::fromUTF8( peak.userLabel() );
          break;
        case kShowPeakEnergyLabel:
          snprintf( buffer, sizeof(buffer), "%.2f keV", (peak.mean()/SandiaDecay::keV) );
          txt = buffer;
          break;
        case kShowPeakNuclideLabel:
          if( peak.parentNuclide() )
          {
            txt = peak.parentNuclide()->symbol;
            if( m_peakLabelsToShow[kShowPeakNuclideEnergies] )
            {
              try
              {
                const float energy = peak.gammaParticleEnergy();
                snprintf( buffer, sizeof(buffer), ", %.2f keV", (energy/SandiaDecay::keV) );
                txt += buffer;
              }catch(...){}
            }
          }else if( peak.reaction() )
          {
            txt = peak.reaction()->name();
            if( m_peakLabelsToShow[kShowPeakNuclideEnergies] )
            {
              snprintf( buffer, sizeof(buffer), ", %.2f keV", (peak.reactionEnergy()/SandiaDecay::keV) );
              txt += buffer;
            }
          }else if( peak.xrayElement() )
          {
            txt = peak.xrayElement()->symbol + " xray";
            if( m_peakLabelsToShow[kShowPeakNuclideEnergies] )
            {
              snprintf( buffer, sizeof(buffer), ", %.2f keV", (peak.xrayEnergy()/SandiaDecay::keV) );
              txt += buffer;
            }
          }//if( peak.parentNuclide() ) / if else ...
          break;
          
        case kShowPeakNuclideEnergies:
        case kNumPeakLabels:
          break;
      }//switch ( i )
      
      if( txt.empty() )
        continue;
      
      const WFont  oldFont  = painter.font();
      const WPen   oldPen   = painter.pen();
      const WBrush oldBrush = painter.brush();
      
      const int npixelx = 12 * static_cast<int>(txt.toUTF8().length());
      const double txt_yval = peakTip.y() - 15.0*(numLabels-labelnum) - 5;
      const double txt_xval = peakTip.x() - 0.5*npixelx;
      WFont font( WFont::Default );
      font.setSize( 12 );
      painter.setPen( WPen(WColor(red,green,blue,alpha)) );
      painter.setFont( font );
      painter.drawText(txt_xval, txt_yval, npixelx, 20.0, AlignCenter, txt);
      ++labelnum;
      
      painter.setPen( oldPen );
      painter.setBrush( oldBrush );
      painter.setFont( oldFont );
    }//for( loop over PeakLabels )
  }//if( we have to write at least one label )
}//void drawPeakText(...)



void SpectrumChart::paintGausPeaks( const vector<std::shared_ptr<const PeakDef> > &peaks,
                                    Wt::WPainter& painter ) const
{
  //XXX - this function is a mess!  Probably because I dont know what I want.
  
  //ToDo 20190421: Changed it so that all `peaks` should have the
  //  same continuum (before this, this function was called with all peaks that
  //  had ROIs that overlapped in energy, but I think parts of this function
  //  assumed shared ROI, but other parts didnt assume this), which means we can
  //  simplify and cleanup this function a bit; notably get rid of
  //  `continuumToPeaks` and logic around it.
  
  const bool outline = (peaks.size() > 1);

  typedef std::map<int,pair<double,double> > IntToDPairMap;
  std::map<int,pair<double,double> > peaksSum, continuumVals;

  double axisMinX, axisMaxX, axisMinY, axisMaxY;
  visibleRange( axisMinX, axisMaxX, axisMinY, axisMaxY );
  
  const SpectrumDataModel *th1Model  = dynamic_cast<const SpectrumDataModel *>( model() );
  
  assert( th1Model );
  assert( std::is_sorted( begin(peaks), end(peaks), &PeakDef::lessThanByMeanShrdPtr ) );
  
  std::shared_ptr<const Measurement> data = th1Model->getData();
  
  bool hadGlobalCont = false;
  const bool sharedContinuum = (peaks.size()>1);
  
  ////////////////////////////////////////////////////////////////////////////
  //(Begin experimental code for multicolored peaks)
  //  If peaks are all the same color, this code results in something that
  //  looks about the same as the old way of doing things, but probably
  //  slightly larger JS since there are more paths to draw
  std::set<string> lineColors;
  for( const auto &p : peaks )
    lineColors.insert( p->lineColor().cssText(false) );
  
  if( lineColors.size() > 1 )
  {
    //Find bin limits we are workging with
    int minrow = std::numeric_limits<int>::max(), maxrow = -std::numeric_limits<int>::max();
    vector<pair<int,int>> peak_bin_ranges;
    //Looping throw all the peaks may be a bit overkill; most of the time just
    //  looking at the left and right most peaks would be sufficient
    for( const auto &p : peaks )
    {
      int minrowlocal, maxrowlocal;
      if( gausPeakDrawRange(minrowlocal, maxrowlocal, *p, axisMinX, axisMaxX) )
      {
        minrow = std::min(minrow,minrowlocal);
        maxrow = std::max(maxrow,maxrowlocal);
        peak_bin_ranges.push_back( {minrowlocal,maxrowlocal} );
      }else
      {
        peak_bin_ranges.push_back( {-999,-999} );
      }
    }//for( const auto &p : peaks )
    
    //Sanity check
    if( maxrow < 0 || maxrow < minrow || ((maxrow-minrow) > 16384) )
      return;
    
    const int nrows = 1 + maxrow - minrow;
    
    vector<double> xvalues(nrows, 0.0);
    //peak_yval: Zeroeth entry is continuum_y
    //  Each entry, j, after that is the sum of the continuum, and all peaks
    //  up to j.  The idea is that to draw the peaks left to right stacking on
    //  top of eachother, we can walk through j, and have peak_yval[j]
    //  give the height of all peaks, up to and including j, and then j-1
    //  will give the lower level of that peak.
    vector<vector<double>> peak_yval( 1 + peaks.size(), vector<double>(nrows,0.0) ); //zeroth entry is continuum_y
    
    //Go through and get xvalues and continuum y values.
    std::shared_ptr<const PeakContinuum> continuum = peaks[0]->continuum();
    
    assert( continuum );  //should always be true
    
    for( int row = minrow; row <= maxrow; ++row )
    {
      const double bin_x_min = th1Model->rowLowEdge( row );
      const double bin_width = th1Model->rowWidth( row );
      const double bin_x_max = bin_x_min + bin_width;
      const double bin_center = bin_x_min + 0.5*bin_width;
      
      double cont_y = 0.0;
      
      switch( peaks[0]->continuum()->type() )
      {
        case PeakContinuum::NoOffset:     case PeakContinuum::Constant:
        case PeakContinuum::Linear:       case PeakContinuum::Quadratic:
        case PeakContinuum::Cubic:        case PeakContinuum::External:
          cont_y = peaks[0]->offset_integral( bin_x_min, bin_x_max, data );
          break;
          
        case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
        case PeakContinuum::BiLinearStep:
        {
          if( data )
            cont_y = peaks[0]->offset_integral( bin_x_min, bin_x_max, data );
          break;
        }
      }//switch( peak.continuum()->type() )
      
      
      if( th1Model->backgroundSubtract() && (th1Model->backgroundColumn() >= 0) )
        cont_y -= th1Model->data( row, th1Model->backgroundColumn() );
      cont_y = std::min( std::max(cont_y, axisMinY), axisMaxY );
      
      for( size_t i = 0; i < peak_yval.size(); ++i )
        peak_yval[i][row-minrow] = cont_y;
      
      xvalues[row-minrow] = std::max( std::min(bin_center, axisMaxX), axisMinX );
      
      assert( peak_bin_ranges.size() == peaks.size() );
      
      for( size_t i = 0; i < peaks.size(); ++i )
      {
        const PeakDef &peak = *(peaks[i]);
        const auto &peak_bin_range = peak_bin_ranges[i];
        if( row < peak_bin_range.first || row > peak_bin_range.second )
          continue;
        try
        {
          double y = peak.gauss_integral( bin_x_min, bin_x_max );
          for( size_t j = i+1; j < peak_yval.size(); ++j )
            peak_yval[j][row-minrow] += y;
        }catch(...)
        {
          cerr << "Caught exception: SpectrumChart::paintGausPeaks(...)" << endl;  //I dont think this will happen, but JIC
        }
      }//for( size_t i = 0; i < peaks.size(); ++i )
    }//for( int row = minrow; row <= maxrow; ++row )
    
    for( size_t j = 1; j < peak_yval.size(); ++j )
    {
      const auto &peak = peaks[j-1];
      int firstbin = peak_bin_ranges[j-1].first;
      int lastbin = peak_bin_ranges[j-1].second;
      if( lastbin == -999 )
        continue;
      
      //Try to constrain how much will need to get drawn, by requiring at
      //  least half a pixel height from the peak fill area before drawing
      //  (for example set of peaks in Ba133 spectrum at 80keV: go from 53->18
      //   bins, and 53->21 bins)
      const double px_threshold = 0.5;
      double diff = 0.0;
      do
      {
        const auto l = mapToDevice( xvalues[firstbin-minrow], peak_yval[j-1][firstbin-minrow] );
        const auto u = mapToDevice( xvalues[firstbin-minrow], peak_yval[j][firstbin-minrow] );
        diff = (l.y() - u.y());
        if( diff < px_threshold )
          ++firstbin;
      }while( (diff < px_threshold) && (firstbin < lastbin) );
      
      do
      {
        const auto l = mapToDevice( xvalues[lastbin-minrow], peak_yval[j-1][lastbin-minrow] );
        const auto u = mapToDevice( xvalues[lastbin-minrow], peak_yval[j][lastbin-minrow] );
        diff = (l.y() - u.y());
        if( diff < px_threshold )
          --lastbin;
      }while( (diff < px_threshold) && (firstbin < lastbin) );
      
      //Note: we are going to bin center here, and not the left/right side
      //      of bins on the left and right side of drawing range; this is
      //      inconsistent with how drawIndependantGausPeak(...) does it,
      //      if the DRAW_PEAKS_TO_BIN_EXTREMES preprocessor variable is set
      //      to true.
      
      const double y_start = std::min(std::max( peak_yval[j][firstbin-minrow],axisMinY), axisMaxY);
      WPointF start_point( mapToDevice( xvalues[firstbin-minrow], y_start ) );
      
      WPainterPath path( start_point );
      for( int bin = firstbin+1; bin <= lastbin; ++bin )
      {
        const double y = std::min(std::max( peak_yval[j][bin-minrow],axisMinY), axisMaxY);
        try{ path.lineTo( mapToDevice( xvalues[bin-minrow], y ) ); }catch(...){}
      }
      
      for( int bin = lastbin; bin >= firstbin; --bin )
      {
        const double y = std::min(std::max( peak_yval[j-1][bin-minrow],axisMinY), axisMaxY);
        try{ path.lineTo( mapToDevice( xvalues[bin-minrow], y ) ); }catch(...){}
      }
      
      path.lineTo( start_point );
      
      WColor color = peak->lineColor().isDefault() ? m_defaultPeakColor : peak->lineColor();
      color.setRgb( color.red(), color.green(), color.blue(), 155 );
      
      painter.fillPath( path, WBrush(color) );
    }//for( size_t j = 1; j < peak_yval.size(); ++j )
    
    //Now draw the outlines for the peaks.
    //  Right now this results in the continuum being colored by the peak
    //  farthest to the right - eh, should implement a custom solution here.
    for( const auto &p : peaks  )
      drawIndependantGausPeak( *p, painter, false );
    
    //(End experimental code for multicolored peaks)
    ////////////////////////////////////////////////////////////////////////////
  }else //if( lineColors.size() > 1 )
  {
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      const PeakDef peak = *(peaks[i]);
      std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
      
      hadGlobalCont |= (continuum->type()==PeakContinuum::External);
      
      std::shared_ptr<const PeakDef> prevPeak, nextPeak;
      if( i )
        prevPeak = peaks[i-1];
      if( i < (peaks.size()-1) )
        nextPeak = peaks[i+1];
      
      //      if( independantContinuum && continuum->isPolynomial() )
      if( !sharedContinuum && continuum->isPolynomial() )
      {
        //drawIndependantGausPeak( peak, painter, false );
        drawIndependantGausPeak( peak, painter, true );
        continue;
      }//if( independantContinuum && peak.continuum().isPolynomial() )
      
      
      //If we are here, we are drawing multiple peaks that share a continuum
      int minrow, maxrow;
      bool visible = gausPeakDrawRange(minrow, maxrow, peak, axisMinX,axisMaxX);
      
      if( !visible )
        continue;
      
      //now limit minrow and maxrow for visible region.
      
      if( outline )
        drawIndependantGausPeak( peak, painter, false );
      
      WPointF start_point(0.0, 0.0);
      try
      {
        double xstart, ystart;
        gausPeakBinValue( xstart, ystart, peak, minrow, th1Model,
                         axisMinX, axisMaxX, axisMinY, axisMaxY,
                         prevPeak, nextPeak );
        xstart = max( th1Model->rowLowEdge( minrow ), axisMinX );
        start_point = mapToDevice( xstart, ystart );
      }catch(...)
      {
        cerr << "Caught exception: in SpectrumChart::paintGausPeaks() from gausPeakBinValue" << endl;
      }
      
      WPainterPath path( start_point );
      
      for( int row = minrow; row <= maxrow; ++row )
      {
        try
        {
          double bin_center, y;
          gausPeakBinValue( bin_center, y, peak, row, th1Model,
                           axisMinX, axisMaxX, axisMinY, axisMaxY,
                           prevPeak, nextPeak );
          if( outline )
          {
            if( peaksSum.find(row) == peaksSum.end() )
              peaksSum[row] = make_pair(0.0,0.0);
            
            pair<double,double> &peakval = peaksSum[row];
            
            const double backHeight = peakBackgroundVal( row, peak,
                                                        th1Model, prevPeak, nextPeak );
            peakval.first = bin_center;
            peakval.second += (y - backHeight);
          }//if( outline )
          
          path.lineTo( mapToDevice( bin_center, y ) );
        }catch(...){}
      }//for( int row = minrow; row <= maxrow; ++row )
      
      
      bool hasDrawnOne = false;
      WPointF peakTip( 999999.9, 999999.9 );
      
      for( int row = maxrow; row >= minrow; --row )
      {
        try
        {
          const double bin_x_min = th1Model->rowLowEdge( row );
          const double bin_width = th1Model->rowWidth( row );
          const double bin_x_max = bin_x_min + bin_width;
          
          bool shouldBeLast = (row==minrow);
          if( bin_x_max < axisMinX )
            shouldBeLast = true;
          
          double bin_center = th1Model->rowCenter( row );
          
          if( !hasDrawnOne )
          {
            hasDrawnOne = true;
            bin_center = bin_x_max;
          }//if( hasDrawnOne )
          
          if( shouldBeLast )
            bin_center = bin_x_min;
          
          bin_center = std::min( bin_center, axisMaxX );
          bin_center = std::max( bin_center, axisMinX );
          
          double val = peakBackgroundVal( row, peak, th1Model,
                                         std::shared_ptr<const PeakDef>(),
                                         std::shared_ptr<const PeakDef>() );
          
          val = std::max( val, axisMinY );
          val = std::min( val, axisMaxY );
          
          const WPointF pointInPx = mapToDevice( bin_center, val );
          path.lineTo( pointInPx );
          
          if( pointInPx.y() < peakTip.y() )
            peakTip = pointInPx;
          
          //          if( outline )
          {
            val = peakBackgroundVal( row, peak, th1Model,
                                    prevPeak, nextPeak );
            val = std::max( val, axisMinY );
            val = std::min( val, axisMaxY );
            continuumVals[row] = make_pair(bin_center,val);
          }
          
          if( shouldBeLast )
            row = minrow;  //same as a break;
        }catch(...){}
      }//for( int row = maxrow; row >= minrow; --row )
      
      path.lineTo( start_point );
      
      if( !outline )
      {
        WColor color = peak.lineColor().isDefault() ? m_defaultPeakColor : peak.lineColor();
        color.setRgb( color.red(), color.green(), color.blue(), 155 );
        
        //if( lineColors.size() > 1 )
        //{
        //}else
        //{
        //}//if( lineColors.size() > 1 ) / else
        
        painter.fillPath( path, WBrush(color) );
      }else
      {
        const WPen oldPen = painter.pen();
        WColor color( peak.lineColor().isDefault() ? m_defaultPeakColor : peak.lineColor() );
        
        painter.setPen( WPen(color) );
        painter.drawPath( path );
        painter.setPen( oldPen );
      }//if( outline )
      
      drawPeakText( peak, painter, peakTip );
    }//for( size_t i = 0; i < peaks.size(); ++i )
    
  }//if( lineColors.size() > 1 ) / else
    

  if( outline )
  {
    //If we're here, all we did was draw the outlines of the individual peaks
    //  Now we need to draw the filled region of all the sums of paths
    std::unique_ptr<WPainterPath> path;

    //first draw the peaks from left to right
    for( const IntToDPairMap::value_type &vt : peaksSum )
    {
      const int bin = vt.first;
      const double bin_center = vt.second.first;
      double val = vt.second.second;

      if( continuumVals.find(bin) != continuumVals.end() )
        val += continuumVals[bin].second;

      const WPointF point = mapToDevice( bin_center, val );

      if( !path )
        path.reset( new WPainterPath( point ) );
      else
        path->lineTo( point );
    }//for( int bin : bins )

    //assert( peaksSum.empty() == !path );
    
    if( path )
    {
      //ToDo 20190421: empirically (because the mess of logic that is this
      //  function is not straightforward to follow) it seems we only get here
      //  if all the peaks in the ROI are the same color - should probably
      //  verify this
      
      //draw the continuum, from right to left
      for( auto iter = continuumVals.rbegin(); iter != continuumVals.rend(); ++iter )
        path->lineTo( mapToDevice( iter->second.first, iter->second.second ) );
     
      WColor color = peaks.front()->lineColor().isDefault() ? m_defaultPeakColor : peaks.front()->lineColor();
      color.setRgb( color.red(), color.green(), color.blue(), 155 );
      painter.fillPath( *path, WBrush(color) );
    }//if( path )
  }else if( hadGlobalCont && continuumVals.size()>2 )
  {
    WPainterPath path;
    std::map<int,pair<double,double> >::const_iterator vt = continuumVals.begin();
    const WPointF start = mapToDevice( vt->second.first, vt->second.second );
    path.moveTo( start );
    for( vt++; vt != continuumVals.end(); ++vt )
      path.lineTo( mapToDevice( vt->second.first, vt->second.second ) );

    WPen pen( m_defaultPeakColor );
    //WPen pen( peak.lineColor().empty() ? m_defaultPeakColor : WColor(WString(peak.lineColor())) );
    
    pen.setWidth( 2 );
    painter.strokePath( path, pen );
  }//if( outline ) / else if( draw global continuum )
}//void paintGausPeak( const PeakDef &peak, Wt::WPainter& painter ) const




void SpectrumChart::paintPeaks( Wt::WPainter& painter ) const
{
  using boost::any_cast;
  
  if( !m_peakModel )
    return;

#if( WT_VERSION >= 0x3030800 )
  auto absModel = model();
#else
  const WAbstractItemModel *absModel = model();
#endif
  
  const SpectrumDataModel *th1Model = dynamic_cast<const SpectrumDataModel *>( absModel );
  if( !th1Model )
    throw runtime_error( "SpectrumChart::paintPeaks(...): stupid programmer"
                         " error, SpectrumChart may only be powered by a"
                         " SpectrumDataModel." );

  const double minx = axis(Chart::XAxis).minimum();
  const double maxx = axis(Chart::XAxis).maximum();

  std::shared_ptr<const Measurement> dataH = th1Model->getData();
  
  std::vector< PeakModel::PeakShrdPtr > gauspeaks, nongauspeaks;
  std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > peaksDeque = m_peakModel->peaks();

  if( peaksDeque )
  {
    for( const PeakModel::PeakShrdPtr &ptr : *peaksDeque )
    {
      double lowx(0.0), upperx(0.0);
      findROIEnergyLimits( lowx, upperx, *ptr, dataH );
      
      if( !ptr || lowx>maxx || upperx<minx )
        continue;

      if( ptr->gausPeak() )
        gauspeaks.push_back( ptr );
      else
        nongauspeaks.push_back( ptr );
    }//for( const PeakModel::PeakShrdPtr &ptr : *peaksDeque )
  }//if( peaksDeque )

  prepare_gaus_peaks_for_iterating( gauspeaks );
  std::vector<std::shared_ptr<const PeakDef> >::const_iterator iter = gauspeaks.begin();
  while( iter != gauspeaks.end() )
  {
    vector<std::shared_ptr<const PeakDef> > peaks_to_draw;
    iter = next_gaus_peaks( iter, gauspeaks.end(), peaks_to_draw, dataH );
    paintGausPeaks( peaks_to_draw, painter );
  }//while( iter != gauspeaks.end() )

  for( const PeakModel::PeakShrdPtr &ptr : nongauspeaks )
    paintNonGausPeak( *ptr, painter );
}//paint( ... )


void SpectrumChart::renderFloatingLegend()
{
  if( m_legendType != FloatingLegend )
    return;
  
  char buffer[64];
  
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  if( !m )
    return;
  
  int nhists = (!!m->getData() + !!m->getSecondData() + !!m->getBackground());

  //Dont draw legend text if theres no spectrum on the chart, or for
  //  time series data, dont show the legend for gamma only data
  if( (nhists < 1) || ((m_dragAction == HighLight) && (nhists < 2)) )
  {
    if( m_legend )
    {
      m_legendMovedSlot.reset();
      delete m_legend;
      m_legend = 0;
    }//if( m_legend )
    
    return;
  }//if( dont show legend )
  
  bool needRender = m_legendNeedsRender;
  
  if( !m_legend )
  {
    m_legend = new AuxWindow( "" );
    m_legend->setResizable( false );
    m_legend->show();
    m_legend->disableCollapse();
    m_legend->addStyleClass( "legend" );
    m_legend->finished().connect( this, &SpectrumChart::disableLegend );
    m_legend->expanded().connect( this, &SpectrumChart::legendExpandedCallback );
    
    //Set the offsets so
    const string jsthis = "$('#" + id() + "')";
    stringstream moveJs;
    
    //WRT == Widget Relative To
    //LTM == Legend Top Margin
    //LRM == Legend Right Margin
    moveJs << "var s=" << jsthis << ";"
    << "if(!s.data('LTM'))s.data('LTM'," << plotAreaPadding(Top)+5 << ");"
    << "if(!s.data('LRM'))s.data('LRM'," << plotAreaPadding(Right)+5 << ");"
    << jsthis  << ".data('legid','" << m_legend->id() << "');"
    << "Wt.WT.AlignLegend('" << id() << "');";
    doJavaScript( moveJs.str() );
    
    WContainerWidget *legendContents = m_legend->contents();
    WContainerWidget *titleBar = m_legend->titleBar();
    
    titleBar->addStyleClass( "legend-title" );
    legendContents->addStyleClass( "legend-body" );
    
    string js = "function(e,s){this.WT.LegendMoved('" + id() + "');}";
    m_legendMovedSlot.reset( new JSlot( js, m_legend ) );
    titleBar->clicked().connect( *m_legendMovedSlot );
    m_legend->show();
    needRender = true;
  }//if( !m_legend )

  if( !needRender )
    return;

  m_legendNeedsRender = false;
  
  WContainerWidget *content = m_legend->contents();
  content->clear();

  
  for( int index = 1; index < model()->columnCount(); ++index )
  {
    // Skip this one if it's not an active data field.
    if( !m->columnHasData( index ) )
      continue;

    WWidget *legendItem = createLegendItemWidget( index );
    content->addWidget( legendItem );

    
    std::shared_ptr<const Measurement> hist;
    
    if( index == m->dataColumn() )
    {
      hist = m->getData();
      
      // place a box to contain live time, dead time, and neutron counts
      WContainerWidget *statsBox = new WContainerWidget( content );
      statsBox->setAttributeValue(
          "style", "margin-left:5ex;margin-top:-0.5em;margin-bottom:0.5em;" );

      const float dataLiveTime = m->dataLiveTime();
      if( dataLiveTime > 0.0f )
      {
        snprintf( buffer, sizeof( buffer ), "Live Time %.1f", dataLiveTime );
        WText *txt = new WText( buffer, statsBox );
        txt->setInline( false );
      }//if( dataLiveTime > 0.0 )

      const float dataRealTime = m->dataRealTime();
      if( dataRealTime > 0.0f )
      {
        snprintf( buffer, sizeof( buffer ), "Real Time %.1f", dataRealTime );
        WText *txt = new WText( buffer, statsBox );
        txt->setInline( false );
      }//if( dataRealTime > 0.0 )

      const float dataNeutronCounts = m->dataNeutronCounts();
      if( (hist && hist->contained_neutron()) || (dataNeutronCounts > 0.0f) )
      {
        snprintf( buffer, sizeof( buffer ), "Neutron Count %.1f",
                  dataNeutronCounts );
        WText *txt = new WText( buffer, statsBox );
        txt->setInline( false );
      }//if( dataNeutronCounts > 0.0 )
    }else if( index == m->secondDataColumn() ||
              index == m->backgroundColumn() )
    {
      float scaledBy, neutronCount;
      if( index == m->secondDataColumn() )
      {
        hist = m->getSecondData();
        scaledBy = m->secondDataScaledBy();
        neutronCount = m->dataRealTime() *
                       m->secondDataNeutronCounts() /
                       m->secondDataRealTime();
      }else
      {
        hist = m->getBackground();
        scaledBy = m->backgroundScaledBy();
        neutronCount = m->dataRealTime() *
                       m->backgroundNeutronCounts() /
                       m->backgroundRealTime();
      }// if( index == m->secondDataColumn() ) / else( index ==
        // m->backgroundColumn() )

      const bool isScaled = (scaledBy > 0.0 && scaledBy != 1.0);
      const bool hasNeutron =
          ( !IsNan( neutronCount ) && !IsInf( neutronCount ) &&
             (neutronCount > 1.0E-3 || (hist && hist->contained_neutron())) );

      WContainerWidget *statsBox = new WContainerWidget( content );
      statsBox->setAttributeValue(
          "style", "margin-left:5ex;margin-top:-0.5em;margin-bottom:0.5em;" );

      const float thislt = (index == m->secondDataColumn())
                               ? m->secondDataLiveTime()
                               : m->backgroundLiveTime();

      if( thislt > 0.0 )
      {
        snprintf( buffer, sizeof( buffer ), "Live Time %.1f", thislt );
        WText *txt = new WText( buffer, statsBox );
        txt->setInline( false );

        if( hasNeutron )
        {
          snprintf( buffer, sizeof( buffer ), "Neutron Count %.1f",
                    neutronCount );
          txt = new WText( buffer, statsBox );
          txt->setToolTip( "Neutron count is scaled for real time to be "
                           "equivalent cps to main data" );
          txt->setInline( false );
        }//if( hasNeutron )


        if( isScaled )
        {
          const float datalt = m->dataLiveTime();
          const float ltscale = datalt / thislt;
          const bool isltscale = ( fabs( scaledBy - ltscale ) < 0.001 );

          string frmt = "Scaled by %.";
          frmt += ( scaledBy >= 0.1 ) ? "2f" : "1E";
          if( isltscale )
            frmt += " for LT";
          snprintf( buffer, sizeof( buffer ), frmt.c_str(), scaledBy );
          txt = new WText( buffer, statsBox );
          txt->setInline( false );
        }// if( isScaled )
      }//if( thislt > 0.0 )
    }//if( index == m->dataColumn() ) / else background or second
  }//for( int index = 0; index < ; ++index )

  wApp->doJavaScript( "Wt.WT.AlignLegend('" + id() + "');" );
}//void renderFloatingLegend()


void SpectrumChart::paintOnChartLegend( Wt::WPainter &painter ) const
{
  //I havent perfected how I'de like the legend to be rendered, but its pretty
  //  close.
  if( m_legendType != OnChartLegend )
    return;
  
  SpectrumDataModel *mdl = dynamic_cast<SpectrumDataModel *>( model() );
  if( !mdl )
    return;
  
  const int nhists = (!!mdl->getData() + !!mdl->getSecondData()
                      + !!mdl->getBackground());
  
  //Dont draw legend text if theres no spectrum on the chart
  if( nhists < 1 )
    return;
  
  //For time series data, dont show the legend for gamma only data
  if( (m_dragAction == HighLight) && (nhists < 2) )
    return;
  
  char buffer[64];
  const double y_min = axis(Chart::YAxis).minimum();
  const double y_max = axis(Chart::YAxis).maximum();
  const double x_min = axis(Chart::XAxis).minimum();
  const double x_max = axis(Chart::XAxis).maximum();
  
  const WPointF topLeft = mapToDevice( x_min, y_max );
  const WPointF bottomRight = mapToDevice( x_max, y_min );
  
  painter.save();
  
  WFont font( WFont::Default );
  font.setSize( 12 );

  const double rowheight = 20.0;
  double width = m_paintedLegendWidth + (nhists>1 ? 10.0 : 0.0);
  
  if( width <= 10.1 )
  {
//    if( painter.device()->features().testFlag(WPaintDevice::HasFontMetrics) )
//    {
//      WFontMetrics meteric = painter.device()->fontMetrics();
//      meteric...
//    }
    width = 13.0*font.sizeLength().toPixels();
  }//if we dont have font metrix
  
  double ystart = topLeft.y();
  const double xstart = bottomRight.x() - width;
  
  for( int index = 1; index < model()->columnCount(); ++index )
  {
    if( !mdl->columnHasData( index ) )
      continue;
    
    std::shared_ptr<const Measurement> hist;
    const Chart::WDataSeries &serie = series(index);
    
    WString title = asString(model()->headerData(index));
    float lt = 0.0f, rt = 0.0f, neutron = 0.0f, sf = 1.0f;
    
    if( index == mdl->dataColumn() )
    {
      hist = mdl->getData();
      if( title.empty() )
        title = "Foreground";
      lt = mdl->dataLiveTime();
      rt = mdl->dataRealTime();
      neutron = mdl->dataNeutronCounts();
    }else if( index == mdl->backgroundColumn() )
    {
      hist = mdl->getBackground();
      if( title.empty() )
        title = "Background";
      lt = mdl->backgroundLiveTime();
      rt = mdl->backgroundRealTime();
      sf = mdl->backgroundScaledBy();
      neutron = mdl->backgroundNeutronCounts();
      neutron *= mdl->dataRealTime() / rt;
    }else if( index == mdl->secondDataColumn() )
    {
      hist = mdl->getSecondData();
      if( title.empty() )
        title = "Second";
      lt = mdl->secondDataLiveTime();
      rt = mdl->secondDataRealTime();
      sf = mdl->secondDataScaledBy();
      neutron = mdl->secondDataNeutronCounts();
      neutron *= mdl->dataRealTime() / rt;
    }else
    {
      continue;
    }
    
    WRectF textArea;
    
    const WFont oldFont = painter.font();
    const WPen oldPen = painter.pen();
    const WBrush oldBrush = painter.brush();
    
    if( nhists > 1 )
    {
      const double x = xstart - 20.0;
      const double y = ystart+0.25*rowheight;
    
      painter.setPen( serie.pen() );
      painter.drawLine( x, y, x+16.0, y );
      
//      textArea = WRectF( xstart, ystart, width, rowheight );
//      painter.drawText( textArea, AlignLeft|AlignTop, title );
//      ystart += rowheight;
    }//if( nhists > 1 )
    
    painter.setPen( m_textPen );
    painter.setFont( font );
    
    bool wroteText = false;
    
    //display live time
    if( lt > 0.0 )
    {
      snprintf( buffer, sizeof(buffer), "Live Time %.1f s", lt );
      textArea = WRectF( xstart, ystart, width, rowheight );
      painter.drawText( textArea, AlignLeft|AlignTop, buffer );
      ystart += rowheight;
      wroteText = true;
    }//if( lt > 0.0 )
    
    
    if( index == mdl->dataColumn() )
    {
      //display live time
      if( rt > 0.0 )
      {
        snprintf( buffer, sizeof(buffer), "Real Time %.1f s", rt );
        textArea = WRectF( xstart, ystart, width, rowheight );
        painter.drawText( textArea, AlignLeft|AlignTop, buffer );
        ystart += rowheight;
        wroteText = true;
      }//if( rt > 0.0 )
    }//if( data column )
    
    //display neutrons
    const bool hasNeut = (!IsNan(neutron) && !IsInf(neutron)
                           && (neutron>1.0E-4 || (hist && hist->contained_neutron())));
    if( hasNeut )
    {
      snprintf( buffer, sizeof(buffer), "%.1f Neutrons", neutron );
      textArea = WRectF( xstart, ystart, width, rowheight );
      painter.drawText( textArea, AlignLeft|AlignTop, buffer );
      ystart += rowheight;
      wroteText = true;
    }//if( hasNeut )
    
    if( index != mdl->dataColumn() )
    {
      //display scale factor
      const bool isScaled = (sf>0.0 && fabs(sf-1.0)>1.0e-3 );
      
      if( isScaled )
      {
        const float datalt = mdl->dataLiveTime();
        const float ltscale = datalt / lt;
        const bool isltscale = (fabs(sf-ltscale) < 1.0e-3);
        
        string frmt = "Scaled %.";
        frmt += (sf >= 0.1) ? "2f" : "1E";
        if( isltscale )
          frmt += " for LT";
        snprintf( buffer, sizeof( buffer ), frmt.c_str(), sf );
        textArea = WRectF( xstart, ystart, width, rowheight );
        painter.drawText( textArea, AlignLeft|AlignTop, buffer );
        ystart += rowheight;
        wroteText = true;
      } // if( isScaled )
    }//if( not data column )
    
    InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
    std::shared_ptr<SpecMeas> meas = ((app && app->viewer()) ? app->viewer()->measurment( SpecUtils::SpectrumType::Foreground ) : std::shared_ptr<SpecMeas>());
      
    if( meas && meas->sample_numbers().size() > 1 )
    {
        const std::set<int> total_sample_nums = meas->sample_numbers();
        const std::set<int> &currentSamples = app->viewer()->displayedSamples( SpecUtils::SpectrumType::Foreground );
        int currentSample = -1;
        
        if( currentSamples.empty() && total_sample_nums.size() )
            currentSample = *(total_sample_nums.begin());
        else if( currentSamples.size() > 1 )
            currentSample = *(currentSamples.rbegin());
        else if( currentSamples.size() )
            currentSample = *(currentSamples.begin());
        
        snprintf( buffer, sizeof(buffer), "Sample: %d/%zu", currentSample , total_sample_nums.size());
        textArea = WRectF( xstart, ystart, width, rowheight );
        painter.drawText( textArea, AlignLeft|AlignTop, buffer );
        ystart += rowheight;
        wroteText = true;
    }
      
      
    if( !wroteText )
    {
      if( !title.empty() )
      {
        textArea = WRectF( xstart, ystart, width, rowheight );
        painter.drawText( textArea, AlignLeft|AlignTop, title );
      }//if( !title.empty() )
      
      ystart += rowheight;
    }//if( !wroteText )
    
    painter.setPen( oldPen );
    painter.setBrush( oldBrush );
    painter.setFont( oldFont );
    
    ystart += 5.0;
  }//for( int index = 1; index < model()->columnCount(); ++index )
  
  painter.restore();
}//void paintOnChartLegend()


void SpectrumChart::legendTextSizeCallback( const float legendwidth )
{
  m_paintedLegendWidth = legendwidth;
  m_legendTextMetric.reset();
}//void legendTextSizeCallback(...)


void SpectrumChart::paintEvent( WPaintDevice *paintDevice )
{
  //This function is called before SpectrumChart::paint(...)
  m_widthInPixels  = paintDevice->width().toPixels();
  m_heightInPixels = paintDevice->height().toPixels();

#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
  {
    const double w = paintDevice->width().toPixels();
    const double h = paintDevice->height().toPixels();
    const int status = setLeftYAxisPadding( w, h );
    if( status < 0 )
      cerr << "SpectrumChart::paintEvent() warning: failed to adjust left axis"
           << endl;
  }
#endif
  
  renderFloatingLegend();
  
  //When we are drawinging the legend on the chart, and its a <cavnas> element,
  //  we dont have server side font metrics, so we will call to the client side
  //  to get them, but only if we havent done this before, and the chart is
  //  actually rendered.
  //Has to be done from this function since the chart wont be rendered to the
  //  client side by the time enableLegend() has been called, and this all must
  //  be done from a non const member function.
  if( wApp
     && (m_legendType == OnChartLegend)
     && (preferredMethod() == WPaintedWidget::HtmlCanvas)
     && (m_widthInPixels>1.0)
     && (m_heightInPixels>1.0)
     && (m_paintedLegendWidth<=-10.0f) )
  {
    //m_paintedLegendWidth is -100 if we have never made it here, so lets
    //  indicate we have made it here
    m_paintedLegendWidth = -1.0f;
    LOAD_JAVASCRIPT(wApp, "js/SpectrumChart.js", "SpectrumChart", wtjsDrawnLegendTextMetric);
    
    m_legendTextMetric.reset( new JSignal<float>( this, "LegendTextMetric" ) );
    m_legendTextMetric->connect( boost::bind(&SpectrumChart::legendTextSizeCallback, this, _1) );
    
    wApp->doJavaScript( "Wt.WT.DrawnLegendTextMetric('"+id()+"');" );
  }//if( preferredMethod() == WPaintedWidget::HtmlCanvas )

  
  //paint the chart
  Wt::Chart::WCartesianChart::paintEvent( paintDevice );
  
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  setRenderWaitingStatus( false );
#endif
}//void paintEvent( WPaintDevice *paintDevice )


void SpectrumChart::clearTimeHighlightRegions( const HighlightRegionType region )
{
  setTimeHighLightRegions( vector< pair<double,double> >(), region );
}//void clearTimeHighlightRegions()


void SpectrumChart::setControlDragDebouncePeriod( int milliseconds )
{
  if( m_overlayCanvas )
    m_overlayCanvas->setControlDragDebounceTimePeriod( milliseconds );
}


void SpectrumChart::setYAxisScale( Wt::Chart::AxisScale scale )
{
  axis( Chart::YAxis ).setScale( scale );
}//void setYAxisScale( Wt::Chart::AxisScale scale )


void SpectrumChart::showHistogramIntegralsInLegend( const bool show )
{
  SpectrumDataModel *the_model = dynamic_cast<SpectrumDataModel *>( model() );
  if( the_model )
    the_model->addIntegralOfHistogramToLegend( show );
}//void showHistogramIntegralsInLegend( const bool show = true );


bool SpectrumChart::legendIsEnabled() const
{
  return (m_legendType != NoLegend);
}//bool legendIdEnabled() const


Signal<> &SpectrumChart::legendEnabled()
{
  return m_legendEnabled;
}//Signal<> &legendEnabled()


Signal<> &SpectrumChart::legendDisabled()
{
  return m_legendDisabled;
}//Signal<> &legendDisabled()


void SpectrumChart::connectWtMouseConnections()
{
  if( !m_doubleClickConnection.connected() )
    m_doubleClickConnection = doubleClicked().connect( this, &SpectrumChart::handleDoubleLeftClick );
  if( !m_wtMouseWentUpConnection.connected() )
    m_wtMouseWentUpConnection = mouseWentUp().connect(  boost::bind( &SpectrumChart::wtMouseUp, this, _1 ) );
  if( !m_wtMouseDraggedConnection.connected() )
    m_wtMouseDraggedConnection = mouseDragged().connect( boost::bind( &SpectrumChart::wtDragged, this, _1 ) );
}//void connectWtMouseConnections()


void SpectrumChart::disconnectWtMouseConnections()
{
  if( m_doubleClickConnection.connected() )
    doubleClicked().disconnect( m_doubleClickConnection );
  if( m_wtMouseWentUpConnection.connected() )
    mouseWentUp().disconnect( m_wtMouseWentUpConnection );
  if( m_wtMouseDraggedConnection.connected() )
    mouseDragged().disconnect( m_wtMouseDraggedConnection );
}//void disconnectWtMouseConnections()


bool SpectrumChart::overlayCanvasEnabled() const
{
  return (m_overlayCanvas && m_overlayCanvas->isEnabled() && !m_overlayCanvas->isHidden());
}//bool overlayCanvasEnabled() const


void SpectrumChart::disableOverlayCanvas()
{
  if( m_overlayCanvas )
  {
    delete m_overlayCanvas;
    m_overlayCanvas = NULL;
  }//if( m_overlayCanvas )

  // XXX - if you have previously hardcoded alignOverlayCanvas()->execJs() into your javascript,
  // I'm not sure what will happen
  if( m_alignOverlayCanvas )
    m_alignOverlayCanvas.reset();

  connectWtMouseConnections();
}//void disableOverlayCanvas()



Wt::JSlot *SpectrumChart::alignOverlayCanvas() //returns NULL if overlay canvas not enabled
{
  return m_alignOverlayCanvas.get();
}



Wt::JSignal<std::string> *SpectrumChart::overlayCanvasJsException()  //returns NULL if not available
{
  if( !m_overlayCanvas )
    return NULL;
  return m_overlayCanvas->jsException();
}//JSignal<string> overlayCanvasJsException()


CanvasForDragging *SpectrumChart::overlayCanvas()
{
  return m_overlayCanvas;
}//CanvasForDragging *overlayCanvas()


void SpectrumChart::enableOverlayCanvas( bool outline, bool highlight, bool altShiftHighlight )
{
  if( !m_overlayCanvas )
  {
    disconnectWtMouseConnections();
    
    if( !m_doubleClickConnection.connected() )
      m_doubleClickConnection = doubleClicked().connect( this, &SpectrumChart::handleDoubleLeftClick );
    
    WApplication *app = WApplication::instance();
    LOAD_JAVASCRIPT(app, "js/SpectrumChart.js", "SpectrumChart", wtjsAlignOverlay);
    
    m_overlayCanvas = new CanvasForDragging( this, outline, highlight, altShiftHighlight );
    m_overlayCanvas->userDragged().connect( this, &SpectrumChart::userDragged );
    m_overlayCanvas->userSingleClicked().connect( this, &SpectrumChart::handleLeftClick );
//    m_overlayCanvas->controlMouseDown().connect( boost::bind( &SpectrumChart::handleControlMouseDown, this, _1 ) );
    m_overlayCanvas->controlMouseMove().connect( boost::bind( &SpectrumChart::handleControlMouseMove, this, _1, _2 ) );
    m_overlayCanvas->rightClickSignal().connect( this, &SpectrumChart::handleRightClick );
    if( m_overlayCanvas->dblTap() )
      m_overlayCanvas->dblTap()->connect( boost::bind( &SpectrumChart::handleDoubleTap, this, _1, _2 ) );
    
//    const string js = "function(e,s){Wt.WT.AlignOverlay('" + id() + "Cover', '" + id() + "' );};";
    const string js = "function(e,s){this.WT.AlignOverlay('" + m_overlayCanvas->id() + "', '" + id() + "' );}";
    m_alignOverlayCanvas.reset( new JSlot( js, this ) );
    
#if(DRAW_GAMMA_LINES_LOG_AND_LIN)
    string gamloglinoptjs = "$('#c" + m_overlayCanvas->id() + "').data('GammaLinesLogAndLin',true);";
#else
    string gamloglinoptjs = "$('#c" + m_overlayCanvas->id() + "').data('GammaLinesLogAndLin',false);";
#endif
    
   app->doJavaScript( gamloglinoptjs );
    
    
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
    m_overlayCanvas->mousedown().connect( this, &SpectrumChart::handleMouseDown );
    m_overlayCanvas->mouseup().connect( this, &SpectrumChart::handleMouseUp );
    m_overlayCanvas->leftDownMouseMove().connect( this, &SpectrumChart::handleLeftMouseMove );
    m_overlayCanvas->altLeftDownMouseMove().connect( this, &SpectrumChart::handleAltLeftMouseMove );
    m_overlayCanvas->mouseWentOut().connect( this, &SpectrumChart::handleMouseLeft );
    m_overlayCanvas->mouseWentOver().connect( this, &SpectrumChart::handleMouseEnter );
    m_overlayCanvas->pinchZoomChange().connect( this, &SpectrumChart::handlePinchZoomChange );
#endif

#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
    const string renderjs = "$('#c" + m_overlayCanvas->id() + "').data('ServerReferencePhotopeaks',true);";
    wApp->doJavaScript( renderjs );
#endif
  }//if( !m_overlayCanvas )
}//void SpectrumChart::enableOverlayCanvas()



void SpectrumChart::disableLegend()
{
  if( m_legendType == NoLegend )
    return;
  
  const bool rerender = (m_legendType == OnChartLegend);
  
  m_legendType = NoLegend;
  
  if( m_legend )
  {
    delete m_legend;
    m_legend = 0;
    
    if( wApp )
      doJavaScript( "$('#" + id() + "').data('legid',null);" );
  }//if( m_legend )
  
  if( rerender )
    update();
  
  m_legendDisabled.emit();
}//void disableLegend()


void SpectrumChart::legendExpandedCallback()
{
  if( !m_legend )
    return;
  
  wApp->doJavaScript( "Wt.WT.AlignLegend('" + id() + "');" );
}//void SpectrumChart::legendExpandedCallback()


void SpectrumChart::enableLegend( const bool forceMobileStyle )
{  
  if( m_legendType != NoLegend )
    return;

  auto wtapp = wApp;
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wtapp );
  if( forceMobileStyle || (app && app->isMobile()) || (BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT) )
  {
    m_legendType = OnChartLegend;
    m_legendEnabled.emit();  //notify observers (mostly to togle menu item)
    
    update(); //force a re-draw
    return;
  }
  
  if( wtapp )
  {
    LOAD_JAVASCRIPT(wtapp, "js/SpectrumChart.js", "SpectrumChart", wtjsAlignLegend);
    LOAD_JAVASCRIPT(wtapp, "js/SpectrumChart.js", "SpectrumChart", wtjsLegendMoved);
  }
  
  m_legendType = FloatingLegend;
  renderFloatingLegend();
  m_legendEnabled.emit();
}//void enableLegend()



void SpectrumChart::setHidden( bool hidden, const Wt::WAnimation &animation )
{
  WCartesianChart::setHidden( hidden, animation );

  
  if( m_legend )
  {
    m_legendMovedSlot.reset();
    delete m_legend;
    m_legend = 0;
  }//if( m_legend )

  if( m_overlayCanvas )
    m_overlayCanvas->setHidden( hidden, animation );
}//void setHidden( bool hidden, const Wt::WAnimation &animation )


#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
int SpectrumChart::setLeftYAxisPadding( double width, double height )
{
  const int oldleft = plotAreaPadding(Wt::Left);
  
//  if( chartArea().width() > 5 && chartArea().height() > 5 )
//  {
//    width = chartArea().width() + plotAreaPadding(Left) + plotAreaPadding(Right);
//    height = chartArea().height() + plotAreaPadding(Top) + plotAreaPadding(Bottom);
//  }
  
#if( WT_VERSION<=0x3030100 )
  initLayout( WRectF(0.0, 0.0, width, height ) );
#else
  const bool inited = initLayout( WRectF(0.0, 0.0, width, height ) );
  if( !inited )
    return -1;
#endif
  
  Chart::WAxis &yaxis = axis(Chart::OrdinateAxis);
  
  vector<MyTickLabel> ticks;
  getYAxisLabelTicks( this, yaxis, ticks );
  
/*
  if( axis(Chart::YAxis).autoLimits() || ticks.empty() )
  {
    double ymin = yaxis.minimum();
    double ymax = yaxis.maximum();

    //I dont thing this section of the code will get called anymore.
    SpectrumDataModel *theModel = dynamic_cast<SpectrumDataModel *>( model() );
    if( !theModel )
      return -1;
    
    const int nbin = theModel->rowCount();
    
    double x0 = axis(Chart::XAxis).minimum();
    double x1 = axis(Chart::XAxis).maximum();
    
    if( axis(Chart::XAxis).autoLimits() )
    {      
      x0 = theModel->rowLowEdge(0);
      x1 = theModel->rowLowEdge(nbin-1) + theModel->rowWidth(nbin-1);
    }//if( axis(Chart::XAxis).autoLimits() )
  
    yAxisRangeFromXRange( x0, x1, ymin, ymax );
    
    getYAxisLabelTicks( this, yaxis, ticks, ymin, ymax );
  }
*/
  
  size_t nchars = 2;
  for( size_t i = 0; i < ticks.size(); ++i )
    nchars = std::max( nchars, ticks[i].label.narrow().size() );
  
  //A better solution would be to use a WPdfImage paint device to measure the
  //  label width...
//  WPdfImage *img = new WPdfImage(100, 100);
//  WPainter *painter = new WPainter(img);
//  painter->setFont( yaxis.labelFont() );
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
    
  setPlotAreaPadding( left, Wt::Left );

  if( m_overlayCanvas )
    m_overlayCanvas->setChartPadding();
  
  return 1;
}//setLeftYAxisPadding()
#endif


void SpectrumChart::setScrollY( int scrollY )
{
  m_scrollY = scrollY;
}


void SpectrumChart::setScrollingParent( Wt::WContainerWidget *parent )
{
  if( m_overlayCanvas )
    m_overlayCanvas->setScrollingParent( parent );
  else
    cerr << "SpectrumChart::setScrollingParent(...)\n\tWarning, calling this funtion has no effect\n";
}//void setScrollingParent( Wt::WContainerWidget *parent )


void SpectrumChart::setOverlayCanvasVisible( bool visible )
{
  if( m_overlayCanvas )
    m_overlayCanvas->setHidden( !visible );
}//void setOverlayCanvasVisible( bool visible )


void SpectrumChart::setDefaultPeakColor( const Wt::WColor &color )
{
  m_defaultPeakColor = color.isDefault() ? ns_default_peakFillColor : color;
}

void SpectrumChart::setOccupiedTimeSamplesColor( const Wt::WColor &color )
{
  m_occupiedMarkerColor = color.isDefault() ? ns_default_occupiedTimeColor : color;
  m_occupiedMarkerColor.setRgb( m_occupiedMarkerColor.red(), m_occupiedMarkerColor.green(), m_occupiedMarkerColor.blue(), 75 );
}

void SpectrumChart::setForegroundHighlightColor( const Wt::WColor &color )
{
  auto &c = m_timeHighlightColors[0];
  c = color.isDefault() ? ns_defaultTimeHighlightColors[0] : color;
  c.setRgb( c.red(), c.green(), c.blue(), 155 );
}

void SpectrumChart::setBackgroundHighlightColor( const Wt::WColor &color )
{
  auto &c = m_timeHighlightColors[1];
  c = color.isDefault() ? ns_defaultTimeHighlightColors[1] : color;
  c.setRgb( c.red(), c.green(), c.blue(), 75 );
}

void SpectrumChart::setSecondaryHighlightColor( const Wt::WColor &color )
{
  auto &c = m_timeHighlightColors[2];
  c = color.isDefault() ? ns_defaultTimeHighlightColors[2] : color;
  c.setRgb( c.red(), c.green(), c.blue(), 75 );
}//void setSecondaryHighlightColor( const Wt::WColor &color )

void SpectrumChart::setAxisLineColor( const Wt::WColor &color )
{
  WPen pen( color );
  axis(Chart::XAxis).setPen( pen );
  axis(Chart::YAxis).setPen( pen );
}//void setAxisLineColor( const Wt::WColor &color )

void SpectrumChart::setChartMarginColor( const Wt::WColor &color )
{
  if( color.isDefault() )
    m_chartMarginBrush = WBrush();
  else
    m_chartMarginBrush = WBrush( color );
}//void setChartMarginColor( const Wt::WColor &color )

void SpectrumChart::setChartBackgroundColor( const Wt::WColor &color )
{
  if( color.isDefault() )
    setBackground( WBrush() );
  else
    setBackground( WBrush(color) );
}//void setChartBackgroundColor( const Wt::WColor &color )

void SpectrumChart::setTextColor( const Wt::WColor &color )
{
  if( color.isDefault() )
    m_textPen = WPen( WColor(GlobalColor::black) );
  else
    m_textPen = WPen( color );
  setTextPen( m_textPen );
}//void setTextColor( const Wt::WColor &color )


void SpectrumChart::saveChartToPng( const std::string &filename )
{
  LOAD_JAVASCRIPT(wApp, "js/SpectrumChart.js", "SpectrumChart", wtjsCanvasToPngDownload);
  
  doJavaScript( "Wt.WT.CanvasToPngDownload('" + id() + "','" + filename +"');" );
}//void saveChartToPng( const std::string &name )


void SpectrumChart::setAvoidMobileMenu( const bool avoid )
{
  m_avoidMobileMenu = avoid;
}//void setAvoidMobileMenu( const bool avoid )


#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
void SpectrumChart::setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )
{
  m_referencePhotoPeakLines = nuc;
  update( WFlags<PaintFlag>(0) );
}//void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )


void SpectrumChart::persistCurrentReferncePhotoPeakLines()
{
  if( m_referencePhotoPeakLines.energies.empty() )
    return;
  
  vector<ReferenceLineInfo>::iterator pos;
  pos = std::find( m_persistedPhotoPeakLines.begin(),
                   m_persistedPhotoPeakLines.end(),
                   m_referencePhotoPeakLines );
  
  if( pos == m_persistedPhotoPeakLines.end()
      && m_referencePhotoPeakLines.displayLines )
    m_persistedPhotoPeakLines.push_back( m_referencePhotoPeakLines );
  
  m_referencePhotoPeakLines.reset();
  
  update( WFlags<PaintFlag>(0) );
}//void persistCurrentReferncePhotoPeakLines()


void SpectrumChart::clearAllReferncePhotoPeakLines()
{
  m_referencePhotoPeakLines.reset();
  m_persistedPhotoPeakLines.clear();
  
  update( WFlags<PaintFlag>(0) );
}//void clearAllReferncePhotoPeakLines()


void SpectrumChart::renderReferncePhotoPeakLines( Wt::WPainter &painter,
                            const ReferenceLineInfo &nuc ) const
{
  painter.save();
  
  const double minx = axis(Chart::XAxis).minimum();
  const double maxx = axis(Chart::XAxis).maximum();
  const double miny = axis(Chart::YAxis).minimum();
  const double maxy = axis(Chart::YAxis).maximum();
  
  const vector<double> &x = nuc.energies;
  const vector<double> &y = nuc.intensities;
  
  vector<double>::const_iterator startiter, enditer;
  startiter = std::lower_bound( x.begin(), x.end(), minx );
  enditer = std::upper_bound( x.begin(), x.end(), maxx );

  if( startiter == enditer )
    return;

  const size_t startindex = startiter - x.begin();
  const size_t endindex = enditer - x.begin();
  
  double max_amp = -1.0;
  
  for( size_t i = startindex; i < endindex; ++i )
  {
#if(DRAW_GAMMA_LINES_LOG_AND_LIN)
    max_amp = std::max( max_amp, miny + yrange*y[i] )
#else
    max_amp = std::max( y[i], max_amp );
#endif
  }//for( size_t i = startindex; i < endindex; ++i )
  
  if( max_amp <= 0.0 )
    return;
  
  const bool isValidColor = !nuc.lineColor.isDefault();
  const WColor color( isValidColor ? nuc.lineColor : WColor("#0000FF") );
  
  for( size_t i = startindex; i < endindex; ++i )
  {
    if( y[i] <= 0.0 )
      continue;
    
#if( DRAW_GAMMA_LINES_LOG_AND_LIN )
    cerr << "renderReferncePhotoPeakLines() Not tested for DRAW_GAMMA_LINES_LOG_AND_LIN=true!" << endl;
    const double yrange = maxy - miny;
    const double yval = miny + yrange*(miny + yrange*y[i])/max_amp;
    WPointF lower = mapToDevice( x[i], miny );
    WPointF upper = mapToDevice( x[i], yval );
#else
    WPointF lower = mapToDevice( x[i], miny );
    WPointF upper = mapToDevice( x[i], maxy );
    lower.setY( lower.y() + axisPadding() );
    const double yvalpx = upper.y() + (1.0 - y[i]/max_amp)*(lower.y()-upper.y());
    upper.setY( yvalpx );
#endif

    //Lets make sure the line is at least visible
    if( (lower.y() - upper.y()) < 2.0 )
      upper.setY( lower.y() - 2.0 );
    
    WPen pen;
    pen.setColor( color );
    pen.setWidth( 1 );
    if( nuc.particlestrs[i] != "gamma" )
      pen.setStyle( Wt::DashLine );
    painter.setPen( pen );
        
    painter.drawLine( lower, upper );
  }//for( size_t i = startindex; i < endindex; ++i )
  
  painter.restore();
}//void renderReferncePhotoPeakLines(...)


void SpectrumChart::renderReferncePhotoPeakLines( Wt::WPainter &painter ) const
{
  for( auto iter = m_persistedPhotoPeakLines.rbegin(); iter != m_persistedPhotoPeakLines.rend(); ++iter )
    renderReferncePhotoPeakLines( painter, *iter );
  
  if( m_referencePhotoPeakLines.energies.size() )
    renderReferncePhotoPeakLines( painter, m_referencePhotoPeakLines );
}//void renderReferncePhotoPeakLines( Wt::WPainter &painter );

#endif


void SpectrumChart::modelChanged()
{
  m_legendNeedsRender = true;
  Wt::Chart::WCartesianChart::modelChanged();
  
  SpectrumDataModel *m = dynamic_cast<SpectrumDataModel *>( model() );
  if( !m && model() )
    throw runtime_error( "SpectrumChart can only be used with SpectrumDataModel models" );
  
  if( m )
    m->dataSet().connect( this, &SpectrumChart::setLegendNeedsRendered );
}//void modelChanged()


void SpectrumChart::setLegendNeedsRendered()
{
  m_legendNeedsRender = true;
}//void setLegendNeedsRenderd();



