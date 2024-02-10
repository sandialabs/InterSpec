#ifndef SpectrumChart_h
#define SpectrumChart_h
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

#include <memory>
#include <vector>
#include <utility>


#include <Wt/WColor>
#include <Wt/WRectF>
#include <Wt/WLength>
#include <Wt/WEvent>
#include <Wt/WPainter>
#include <Wt/WAnimation>
#include <Wt/WPaintDevice>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/ReferenceLineInfo.h"

//Forward declarations
namespace Wt
{
  class WWidget;
  class WSpinBox;
  class WCheckBox;
  class WPopupMenu;
  class WButtonGroup;
  class WDoubleSpinBox;
  class WContainerWidget;
}//namespace Wt

class AuxWindow;
class PeakModel;

//See comments about SpectrumChart::setLeftYAxisPadding() for what
//  DYNAMICALLY_ADJUST_LEFT_CHART_PADDING controls.  This feature has not
//  been tested well enough to fully use yet, although it seems to work well.
#define DYNAMICALLY_ADJUST_LEFT_CHART_PADDING 1

class SeriesRenderer;
class LabelRenderIterator;
class MarkerRenderIterator;
class SpectrumSeriesRenderer;
class SpectrumRenderIterator;


static_assert( WT_VERSION > 0x3030100, "Wt should be 3.3.2 or newer (really 3.3.4 at least)" );

class SpectrumChart : public Wt::Chart::WCartesianChart
{
//XXX
// Whenever you load a new data histogram into the model which powers this
// chart, please call SpectrumChart::setAutoXAxisRange() as well, this is a
// workaround, for more info see XXX note in SpectrumChart::handleDrag(...)
 
public:
  enum XAxisUnits
  {
    kkeV,
    kSeconds,
    kUndefinedUnits
  };//enum XAxisUnits

  enum HighlightRegionType
  {
    ForegroundHighlight,
    BackgroundHighlight,
    SecondForegroundHighlight
  };//enum HighlightRegionType
  
  
  enum class OccupancyRegionType
  {
    Occupied,
    NotOccupied
  };//enum HighlightRegionType
  
  //HighlightRegion is what is used to overlay a color over the chart, and are
  //  used for roughly two different usages:
  //    1) On the time history chart, show which time segments are used to sum
  //       for the foreground, background, and secondary spectrums.
  //       These regions are set and cleared by calling:
  //       setTimeHighLightRegions(...), and clearTimeHighlightRegions(...).
  //    2) Used by other tools to highlight an energy range on the spectrum
  //       chart.  These regions are added and removed by
  //       addDecorativeHighlightRegion(...), and
  //       removeDecorativeHighlightRegion(...).  When a region is added, a
  //       unique ID is returned, so that region can be removed, without
  //       effecting other decorative regions.
  struct HighlightRegion
  {
    float lowerx, upperx;
    size_t hash;  /* Hash values of 0, 1, and 2 are reserved for to indicate
                     ForegroundHighlight, BackgroundHighlight, or 
                     SecondForegroundHighlight, respectively*/
    Wt::WColor color;
  };//HighlightRegion
  
  
  enum PeakLabels
  {
    kShowPeakUserLabel,
    kShowPeakEnergyLabel,
    kShowPeakNuclideLabel,
    kShowPeakNuclideEnergies,
    kNumPeakLabels
  };//enum PeakLabels
  
public:
  SpectrumChart( Wt::WContainerWidget *parent = NULL );
  virtual ~SpectrumChart();

  void setPeakModel( PeakModel *model );
  
  //Note: do not overide height() and width() or else when used in a stretching
  //      layout side will get really messed up after the first painting
  Wt::WLength paintedWidth() const;
  Wt::WLength paintedHeight() const;
  
  //setTextInMiddleOfChart(...): draws some large text over the middle of the
  //  chart - used int the spectrum quizzer for text based questions
  void setTextInMiddleOfChart( const Wt::WString &txt );
  
  //textInMiddleOfChart(): returns current m_textInMiddleOfChart
  const Wt::WString &textInMiddleOfChart() const;

  /** Set whether to slim down axis for small displays (e.g. phone).
      If setting to compact from non-compact, 17px of bottom margin will be
      removed, or if goign the other way 17px added.
      It seems 42px bottom padding for non-compact, and 25px good for compact
      is a good amount of padding to start with if you would like to customize
      padding after calling this nuction.
   
      Default is to not have compact.
      Currently only effects x-axis.
   */
  void setCompactAxis( const bool compact );
  bool isAxisCompacted() const;

  
  void showGridLines( const bool draw ); //shows horizantal and vertical
  void showVerticalLines( const bool draw );
  void showHorizontalLines( const bool draw );
  bool verticalLinesShowing() const;  // Added by christian (20170425)
  bool horizontalLinesShowing() const;
  
  void showHistogramIntegralsInLegend( const bool show = true );

  void setPeakLabelColor( PeakLabels label, uint32_t red, uint32_t green,
                          uint32_t blue, uint32_t alpha = 255 );
  void setShowPeakLabel( PeakLabels label, bool show );
  bool showingPeakLabel( PeakLabels label ) const;
  
  virtual void paint( Wt::WPainter& painter,
                     const Wt::WRectF& rectangle = Wt::WRectF() ) const;
  virtual void paintEvent( Wt::WPaintDevice *paintDevice );
  
  //renderFloatingLegend(): only renders if valid m_legend and
  //  m_legentType==FloatingLegend
  virtual void renderFloatingLegend();
  
  //paintOnChartLegend(): only paints if m_legentType==OnChartLegen
  virtual void paintOnChartLegend( Wt::WPainter &painter ) const;
  void paintTextInMiddleOfChart( Wt::WPainter &painter ) const;

  virtual void paintPeaks( Wt::WPainter& painter ) const;
  virtual void paintNonGausPeak( const PeakDef &peak, Wt::WPainter &painter ) const;

  virtual void drawPeakText( const PeakDef &peak, Wt::WPainter &painter,
                             Wt::WPointF peakMaxInPx ) const;
  
  //peakModelDataChanged(...): decides if the whole chart should be re-drawn
  //  if the peak mode data is changed.
  virtual void peakModelDataChanged( const Wt::WModelIndex &topLeft,
                                     const Wt::WModelIndex &bottomRight );

  //peakYVal(...): returns the y-value to paint for the center of the bin, for
  //  given the peak; includes the background value
  static double peakYVal( const int bin, const PeakDef &peak,
                          const SpectrumDataModel *th1Model,
                          std::shared_ptr<const PeakDef> prevPeak,
                          std::shared_ptr<const PeakDef> nextPeak );
  
  //gausPeakBinValue(...): gets the x and y values for the marker, for bin
  //  of data specified.
  //May throw exception
  void gausPeakBinValue( double &xvalue, double &yvalue,
                         const PeakDef &peak, const int bin,
                         const SpectrumDataModel *th1Model,
                         const double axisMinX, const double axisMaxX,
                         const double axisMinY, const double axisMaxY,
                         std::shared_ptr<const PeakDef> prevPeak,
                         std::shared_ptr<const PeakDef> nextPeak ) const;

  //peakBackgroundVal(...): gives the bottom of the peak (e.g. the continuum)
  //  to paint.  Acounts for data background subtraction as well.
  static double peakBackgroundVal( const int bin,
                                   const PeakDef &peak,
                                   const SpectrumDataModel *th1Model,
                                   std::shared_ptr<const PeakDef> prevPeak,
                                   std::shared_ptr<const PeakDef> nextPeak );
  

  //paintGausPeaks(): Paints peaks that share a ROI.  All peaks passed into
  //  this function must share a PeakContinuum, be sorted by peak mean, and be
  //  valid pointers in the vector (e.g., no nullptrs).
  virtual void paintGausPeaks( const std::vector<std::shared_ptr<const PeakDef> > &peaks,
                              Wt::WPainter &painter ) const;

  //drawIndependantGausPeak(...): draw a peak, not taking into account any other
  //  peaks.  If doFill==false then only the outline for the peak is drawn, or
  //  else the fill area and outline of the peak will be drawn.
  void drawIndependantGausPeak( const PeakDef &peak,
                                Wt::WPainter& painter,
                                const bool doFill = false,
                                const double xstart = -999.9,
                                const double xend = -999.9 ) const;
  
  void visibleRange( double &xmin, double &xmax,
                     double &ymin, double &ymax ) const;

  float xUnitsPerPixel() const;
  
  //gausPeakDrawRange(...): returns if the peak is visible at all
  bool gausPeakDrawRange( int &minbin, int &maxbin,
                          const PeakDef &peak,
                          const double axisMinX,
                          const double axisMaxX ) const;
  
  virtual void paintHighlightRegions( Wt::WPainter& painter ) const;

  virtual void setYAxisScale( Wt::Chart::AxisScale scale );
  virtual void setXAxisRange( const double x1, const double x2 );
  virtual void setYAxisRange( const double y1, const double y2 );
  virtual void yAxisRangeFromXRange( const double x0, const double x1,
                                     double &y0, double &y1 );
  virtual void setAutoXAxisRange();  //sets x-range to full range
  virtual void setAutoYAxisRange();  //sets y-range, from current x-range,
                                     //  note that if the rebin factor of the
                                     //  SpectrumDataModel is updated, or the data is
                                     //  changed, you must re-call this funtion

  virtual void setXAxisUnits( XAxisUnits units );
  virtual XAxisUnits xAxisUnits() const;

  //Adding or clearing the highlight regions does not cause an update to the
  //  chart; you need to manually call update() to cause a rerender.
  void clearTimeHighlightRegions( const HighlightRegionType region );
  void setTimeHighLightRegions( const std::vector< std::pair<double,double> > &p,
                            const HighlightRegionType type ); //does not refresh the chart
  
  void clearOccupancyRegions();
  void setOccupancyRegions( const std::vector< std::pair<double,double> > &p ); //does not refresh the chart
  
  
  //Decorative highlights get assigned an ID, so this way you can only remove
  //  a decorative highlight, if you know its ID.  This is an attempt to allow
  //  multiple widgets to have multiple regions that they can manage.
  //
  //removeDecorativeHighlightRegion(): returns if the removal was successful.
  //  If regionid is less than 3, no action will be taken.
  bool removeDecorativeHighlightRegion( size_t regionid );
  
  //addDecorativeHighlightRegion(): returns the ID for the highlihgted region
  //  added; you will need to retain this in order to remove the region in the
  //  future.  Returned ID will alway be greater than 2.
  size_t addDecorativeHighlightRegion( const float lowerx,
                                       const float upperx,
                                       const Wt::WColor &color );

  
  //enableLegend():  Makes it so the legend will be rendered, with some caviots.
  //  If no spectrums are present, will not be rendered.
  //  If a time series chart, only rendered if more than
  //  one sereies is to be plotted.
  //  For non-phone devices, creates a globablly floating legend (AuxWindow)
  //  that will be kept at a consistent distance from the top and right hand
  //  side of chart.
  //  For phone devices, the legend will be rendered directly on the chart.
  void enableLegend( const bool forceMobileStyle );
  
  //disableLegend(): causes legend to stop being drawn.  Deletes m_legend if
  //  it is non-null.
  void disableLegend();
  
  bool legendIsEnabled() const;
  
  Wt::Signal<> &legendEnabled();
  Wt::Signal<> &legendDisabled();

  //legendExpandedCallback(): rerenders the legend after its been expanded
  void legendExpandedCallback();
  
  void legendTextSizeCallback( const float legendwidth );

  void setDefaultPeakColor( const Wt::WColor &color );
  
  /** Set the color of the vertical lines (and text) that mark the beggining and
      end of an occupancy on a time chart.
      Currently alpha is ignored and always set to 75.
   */
  void setOccupiedTimeSamplesColor( const Wt::WColor &color );
  
  /** */
  void setForegroundHighlightColor( const Wt::WColor &color );
  void setBackgroundHighlightColor( const Wt::WColor &color );
  void setSecondaryHighlightColor( const Wt::WColor &color );

  /** Set the color of axis label, axis title, and other text in the chart. */
  void setTextColor( const Wt::WColor &color );
  
  /** Set the color of the axis lines, and ticks of the chart. */
  void setAxisLineColor( const Wt::WColor &color );
  
  /** Set the color of the chart margin area (the area outside the axis lines).
      If you set this to an uninitaialized color (default constructed WColor),
      then the margins will be painted to the same color as the chart background
      (if the background color is set).
   */
  void setChartMarginColor( const Wt::WColor &color );
  
  /** Set the background color of the chart area, and if the margin color isnt
      set, then that too.
   */
  void setChartBackgroundColor( const Wt::WColor &color );
  
  /** Executes appropriate javascript to generate and download a PNG based on
   the currently showing spectrum.  PNG generation is done client side.
   */
  void saveChartToPng( const std::string &name );

  virtual void setHidden( bool hidden,
                          const Wt::WAnimation &animation = Wt::WAnimation() );

#if(DYNAMICALLY_ADJUST_LEFT_CHART_PADDING)
  //setLeftYAxisPadding(): tries to guess how many characters will be used for
  //  the y-axis label, and then sets the left padding of the chart accordingly,
  //  updating the client side mouse coordinate equations as well.
  //  The width and height are as given by the WPaintDevice the chart will be
  //  painted which is usually equivalent to
  //    width = chartArea().width() + plotAreaPadding(Left) + plotAreaPadding(Right);
  //    height = chartArea().height() + plotAreaPadding(Top) + plotAreaPadding(Bottom);
  //  Returns -1 on error, 0 if no adjustment was necassary, and 1 if adjustment
  //  was made.
  int setLeftYAxisPadding( double width, double height, Wt::WPaintDevice *paintDevice );
#endif
  
  void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc );
  void persistCurrentReferncePhotoPeakLines();
  void clearAllReferncePhotoPeakLines();
  void renderReferncePhotoPeakLines( Wt::WPainter &painter ) const;
  void renderReferncePhotoPeakLines( Wt::WPainter &painter,
                          const ReferenceLineInfo &line ) const;

protected:
  //modelChanged(): hooks up the signals to catch when the legend needs to be
  //  re-rendered.  Also enforces that only SpectrumDataModel models can be
  //  used with this class, by throwing an exception otherwise.
  virtual void modelChanged();
  
  
  void setLegendNeedsRendered();
  
protected:
  void labelRender( Wt::WPainter& painter, const Wt::WString& text,
                   const Wt::WPointF& p, const Wt::WColor& color,
                   Wt::WFlags<Wt::AlignmentFlag> flags,
                   double angle, int margin) const;
  
  virtual void render( Wt::WPainter &painter, const Wt::WRectF &rectangle ) const;
  
  virtual void renderChartBackground( Wt::WPainter &painter, const Wt::WRectF &rectangle ) const;  //Take places of WCartesianChart::renderBackground(...)
  virtual void renderGridLines( Wt::WPainter &painter, const Wt::Chart::Axis axis ) const;
  virtual void renderSeries( Wt::WPainter &painter ) const;
  virtual void renderAxes( Wt::WPainter &painter, Wt::WFlags<  Wt::Chart::AxisProperty > properties ) const;
  virtual void renderXAxis( Wt::WPainter &painter, const Wt::Chart::WAxis& axis,
                            Wt::WFlags<Wt::Chart::AxisProperty> properties ) const;
  virtual void renderYAxis( Wt::WPainter &painter, const Wt::Chart::WAxis &axis,
                           Wt::WFlags<Wt::Chart::AxisProperty> properties ) const;
  virtual void renderLegend( Wt::WPainter &painter ) const;
  virtual bool isLargeEnough( Wt::WPainter &painter ) const;
  
  virtual void iterateSpectrum( SpectrumRenderIterator *iterator,
                                Wt::WPainter &painter ) const;
  
  //since we cant access WCartesianChar::chartArea_ ... (stupid Wt not allowing access)
  //  Note that providing a default WRectF() is a hack and will lead to
  //  rendering
  mutable Wt::WRectF m_renderRectangle;
  void calcRenderRectangle( const Wt::WRectF &rectangle ) const;
  const Wt::WRectF &chartArea() const;
  
  friend class SeriesRenderer;
  friend class LabelRenderIterator;
  friend class MarkerRenderIterator;
  friend class SpectrumRenderIterator;
  friend class SpectrumSeriesRenderer;
  
protected:
  double m_widthInPixels;        //The width of the chart in pixels
  double m_heightInPixels;       //The height of the chart in pixels
  
  PeakModel *m_peakModel;
  Wt::WColor m_defaultPeakColor;
  
  /** The vertical lines (and text) that mark the beggining and end of an
      occupancy on a time chart.
   */
  Wt::WColor m_occupiedMarkerColor;
  
  Wt::WColor m_timeHighlightColors[3];
  
  /** If a chart margin brush is not specifified but WAbstractChart::background()
       is, then the entire graphic area (e.g. <canvas>) will be filled according
       to WAbstractChart::background().
      If both chart margin and background are specified, then
       WAbstractChart::background() will be used to fill the area inside the chart
       axis, and m_chartMarginBrush used to fill area outside the axis.
      If only m_chartMarginBrush is specified, then only area outside axis
       will be painted.
   */
  Wt::WBrush m_chartMarginBrush;
  
  /** Since we dont have acess to WCartesianChart::textPen_, we will track
      text color our selves.
   */
  Wt::WPen m_textPen;
  
  //Potential ToDo: Right now m_highlights is mainly modified and added
  //  to by this class, however the effects of what is highlighted are also seen
  //  in the InterSpec class, which is a possibility for the displayed time
  //  range to get out of sync with the highlighted tiem ranges.  It might make
  //  sense for this class to never modify m_highlights, but instead
  //  just have the InterSpec class say which regions should be highlighted

  std::vector< HighlightRegion > m_highlights;
  
  std::vector< HighlightRegion > m_occupancy_regions;
  
  
  enum LegendType
  {
    NoLegend,
    FloatingLegend,
    OnChartLegend
  };//enum LegendType
  
  //m_legend: valid pointer only when m_legendType==FloatingLegend
  AuxWindow *m_legend;
  
  //m_legendType: indicates the type of legend that should be rendered.  Decided
  //  in enableLegend() based on the device type.
  LegendType m_legendType;
  
  //m_legendNeedsRender: the legend doesnt need to be rendered everytime the
  //  chart is updated, but only when the data changes, so this variable tracks
  //  that.
  bool m_legendNeedsRender;
  
  //m_paintedLegendWidth: for m_legendType==OnChartLegend, with the chart being
  //  rendered to a <canvas> element, we dont have server side font metrics, so
  //  will execute some javascript to measure the font metrics, and then call
  //  back to legendTextSizeCallback(...) to set the following variable.
  float m_paintedLegendWidth;
  std::unique_ptr< Wt::JSignal<float> > m_legendTextMetric;
  

  //
  bool m_peakLabelsToShow[kNumPeakLabels];
  uint32_t m_peakLabelsColors[kNumPeakLabels]; //{red, green, blue, alpha}
  
  //if m_xAxisUnits==kUndefinedUnits, then mouse coordinates wont be displayed
  XAxisUnits m_xAxisUnits;
  
  bool m_compactAxis;
  Wt::WString m_textInMiddleOfChart;
  
  // Added by christian (20170425)
  bool m_verticalLinesShowing;
  bool m_horizontalLinesShowing;
    
  ReferenceLineInfo m_referencePhotoPeakLines;
  std::vector<ReferenceLineInfo> m_persistedPhotoPeakLines;
};//class SpectrumChart




 #endif
