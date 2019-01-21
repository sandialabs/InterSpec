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
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/SpectrumDataModel.h"
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
#include "InterSpec/ReferenceLineInfo.h"
#endif

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
#if( WT_VERSION<=0x3030100 )
  namespace Chart
  {
    class WChart2DRenderer;
    enum AxisProperty{ Labels = 0x1, Line = 0x2 };
    W_DECLARE_OPERATORS_FOR_FLAGS(AxisProperty)
  }//namespace Chart
#endif
}//namespace Wt

class AuxWindow;
class PeakModel;
class CanvasForDragging;

//See comments about SpectrumChart::setLeftYAxisPadding() for what
//  DYNAMICALLY_ADJUST_LEFT_CHART_PADDING controls.  This feature has not
//  been tested well enough to fully use yet, although it seems to work well.
#define DYNAMICALLY_ADJUST_LEFT_CHART_PADDING 1

#if(WT_VERSION<=0x3030100)
class SpectrumRenderer;
#endif

//#if(WT_VERSION>0x3030100)
class SeriesRenderer;
class LabelRenderIterator;
class MarkerRenderIterator;
class SpectrumSeriesRenderer;
class SpectrumRenderIterator;
//#endif

class SpectrumChart : public Wt::Chart::WCartesianChart
{
  /****************************************************************************\
  | This class provides the basic viewing functionality for the interactive
  | analysis of gamma spectra. User input (zooming, clicking, etc.) is handled
  | in one of two ways:
  |   1. Via the CanvasForDragging overlay canvas (primary method)
  |   2. As a fallback for browsers that don't support the <canvas> tag (just IE?)
  |      the dragged(...) and mouseUp(...) functions handle the mouseWentUp() and
  |      mouseDragged() signals. This has been tested and is current with #1.
  |
  | Will likely eventually be replaced by D3.js based implementation in
  | D3SpectrumDisplayDiv.cpp.
  \****************************************************************************/

//XXX
// Whenever you load a new data histogram into the model which powers this
// chart, please call SpectrumChart::setAutoXAxisRange() as well, this is a
// workaround, for more info see XXX note in SpectrumChart::handleDrag(...)
 
public:
  enum DragAction
  {
    ZoomIn,
    HighLight,
    NoAction
  };//enum DragAction

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

  //graphableAreaWidth(): the width of the actual chart area, in pixels
//  int graphableAreaWidth() const;
  
  //Note: do not overide height() and width() or else when used in a stretching
  //      layout side will get really messed up after the first painting
  Wt::WLength paintedWidth() const;
  Wt::WLength paintedHeight() const;
  
  //setTextInMiddleOfChart(...): draws some large text over the middle of the
  //  chart - used int the spectrum quizzer for text based questions
  void setTextInMiddleOfChart( const Wt::WString &txt );
  
  //textInMiddleOfChart(): returns current m_textInMiddleOfChart
  const Wt::WString &textInMiddleOfChart() const;

  //setCompactAxis(): whether to slim down axis for small displays (e.g. on
  //  phone).  Note that effects wont be seen until next time chart is rendered.
  //  You should also adjust padding axis title text appropriately; x-axis
  //  padding of 23px seems to be a reasonable value.
  //  Default is to not have compact axis.
  //Currently only effects x-axis.
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
  //  given the peak; includes the backgorund value
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

  
  typedef std::shared_ptr<const PeakDef> PeadDefPtr;
  typedef std::vector<PeadDefPtr > PeadDefPtrVec;
  typedef PeadDefPtrVec::const_iterator PeakDefPtrVecIter;
  //next_peaks(...): hmmm, assumes all pointers are valid
  static PeakDefPtrVecIter next_peaks( PeakDefPtrVecIter peak_start,
                                       PeakDefPtrVecIter peak_end,
                                       PeadDefPtrVec &peaks,
                                       std::shared_ptr<const Measurement> data );

  //paintGausPeak(): Paints peaks which are in a connected region.
  virtual void paintGausPeak( const std::vector<std::shared_ptr<const PeakDef> > &peaks,
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

  virtual void loadPeakInfoToClient() const;
  
  virtual void loadXAndYRangeEquationToClient() const;
  std::string pixelToCoordinateMappingJs( const Wt::Chart::Axis axis ) const;
  std::string energyToPixelMappingJs() const;
//#if(DRAW_GAMMA_LINES_LOG_AND_LIN)
  std::string countsToPixelMappingJs() const;
//#endif
  
  //setShowRefLineInfoForMouseOver(): set wether or not the text information
  //  should be shown for the line that the mouse is currently over.  Default is
  //  to show the information.
  void setShowRefLineInfoForMouseOver( const bool show );

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

  //The following are <x-begin,x-end> coordinates in units of energy for the
  //  dragging action.
  Wt::Signal<double,double> &xRangeChanged();
  Wt::Signal<double,double,int/*last pageX*/,int/*last pageY*/> &controlKeyDragged();
  Wt::Signal<double,double> &shiftKeyDragged();
  Wt::Signal<double,double> &altKeyDragged();
  Wt::Signal<double,double> &shiftAltKeyDragged();
  Wt::Signal<double,double> &rightMouseDragg();
  
  //For the following, the doubles are the <x,y> coordinates of the action, in
  //  units of energy and counts
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> &leftClick();
  Wt::Signal<double,double> &doubleLeftClick();
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> &rightClick();

  //The following emits the enrgy of the start and current energy of current
  //  action
  Wt::Signal<double,double> &controlMouseMoved();

  //The keyPressedWhileMousedOver() function/signal may only be called
  //  after enableOverlayCanvas(bool,bool) has been called, otherwise will
  //  throw an exception.
  Wt::JSignal<Wt::WKeyEvent> &keyPressedWhileMousedOver();

  //setControlDragDebouncePeriod(...): set how often the client will phone home
  //  to the server during a control-drag event.  Periods of <= 0 will result
  //  in it happening whenever the mouse is moved.
  void setControlDragDebouncePeriod( int milliseconds );

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

  void setDragAction( const DragAction action );

  //below only applicable when setDragAction(HighLight) has been called
  void allowMultipleRegionHighlight( bool allow = true );

  //below only applicable when setDragAction(HighLight) has been called
  void allowSingleClickHighlight( bool allow = true );

  //below may only be called after enableOverlayCanvas(...), and will only
  //  have an effect if allowSingleClickHighlight(true) has been called, and
  //  there is exactly one singly-clicked highlighted region
  void allowArrowToMoveSingleClickRegion( bool allow = true );

  //enableLegend():  Makes it so the legend will be rendered, with some caviots.
  //  If no spectrums are present, will not be rendered.
  //  If a time series chart (DragAction==Highlight), only rendered if more than
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

  
  void enableOverlayCanvas( bool outline, bool highlight, bool altShiftHighlight );
  CanvasForDragging *overlayCanvas();
  void disableOverlayCanvas();
  bool overlayCanvasEnabled() const;
  void setOverlayCanvasVisible( bool visible );

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
  
  /** When the mobile (hamburger) menu is showing, we shouldnt place text (the
   y-axis labels) underneath it; it looks bad.
   */
  void setAvoidMobileMenu( const bool avoid );
  
  //setScrollingParent(...): sets the scrolling frame which contains the chart
  //  from which the overlay canvas should not extend beyond.  Calling this
  //  function when the overlay canvas is not enabled has no effect. Calling
  //  this function with a NULL argument removes this parent.
  void setScrollingParent( Wt::WContainerWidget *parent );
  Wt::JSlot *alignOverlayCanvas();  //returns NULL if overlay canvas not enabled

  //overlayCanvasJsException() is mostly for debugging and will
  //  probably be romived in the future
  Wt::JSignal<std::string> *overlayCanvasJsException();  //returns NULL if not available

  //functions to connect/disconnect the Wt::Chart::WCartesianChart based
  //  mouse signals.  If you are using the overlay canvas, you can (but dont
  //  have to) disconnect these signals.  If you are not using the overlay
  //  canvas, than you should connect these signals.
  void connectWtMouseConnections();
  void disconnectWtMouseConnections();

  virtual void setHidden( bool hidden,
                          const Wt::WAnimation &animation = Wt::WAnimation() );

  //setScrollY(...): sets the number of pixels the page has been scrooled, to
  //  help interpret clicked loacations - note that this is a hack for the
  //  anthony NM app.
  void setScrollY( int scrollY );

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
  int setLeftYAxisPadding( double width, double height );
#endif
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc );
  void persistCurrentReferncePhotoPeakLines();
  void clearAllReferncePhotoPeakLines();
  void renderReferncePhotoPeakLines( Wt::WPainter &painter ) const;
  void renderReferncePhotoPeakLines( Wt::WPainter &painter,
                          const ReferenceLineInfo &line ) const;
#endif

protected:
  //modelChanged(): hooks up the signals to catch when the legend needs to be
  //  re-rendered.  Also enforces that only SpectrumDataModel models can be
  //  used with this class, by throwing an exception otherwise.
  virtual void modelChanged();
  
  
  void setLegendNeedsRendered();
  
protected:
#if(WT_VERSION<=0x3030100)
  friend class SpectrumRenderer;
  Wt::Chart::WChart2DRenderer *createRenderer( Wt::WPainter &painter,
                                            const Wt::WRectF &rectangle ) const;

  //The functions bellow here were added to WCartesianChart in newer versions
  //  of Wt, so to be able to use the same rendering code, we had reimplement.
  Wt::WRectF chartSegmentArea( Wt::Chart::WAxis yAxis, int xSegment,
                                            int ySegment) const;
  Wt::WPointF map( double xValue, double yValue,
                   Wt::Chart::Axis yAxis = Wt::Chart::OrdinateAxis,
                   int 	currentXSegment = 0,
                  int 	currentYSegment = 0 ) const;
  
  Wt::WPointF hv( double x, double y ) const;
  Wt::WPointF hv( const Wt::WPointF &p ) const;
  Wt::WRectF hv( const Wt::WRectF &p ) const;
  void renderLabel( Wt::WPainter& painter, const Wt::WString& text,
                   const Wt::WPointF& p, const Wt::WColor& color,
                                    Wt::WFlags<Wt::AlignmentFlag> flags,
                   double angle, int margin) const;
#endif

  void labelRender( Wt::WPainter& painter, const Wt::WString& text,
                   const Wt::WPointF& p, const Wt::WColor& color,
                   Wt::WFlags<Wt::AlignmentFlag> flags,
                   double angle, int margin) const;
  
#if(WT_VERSION<=0x3030100)
  virtual void customRender( Wt::WPainter &painter, const Wt::WRectF &rectangle );
#else
  virtual void render( Wt::WPainter &painter, const Wt::WRectF &rectangle ) const;
#endif
  
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
  
  
  //wtDragged and wtMouseUp are backup functions incase there is now overlay
  //  canvas, and are only called in this case.  These functions then call the
  //  appropriate functions to provide similar functionality to what wede expect
  //  with the overlay canvas.
  void wtDragged( Wt::WMouseEvent event );
  void wtMouseUp( Wt::WMouseEvent event );

  
#if( USE_OverlayDragEvent )
  void userDragged( const OverlayDragEvent &event );
#else
  void userDragged( int x0, int x1, Wt::WMouseEvent event );
#endif

  
  void handleLeftClick( const int x1, const int y1, const int modifier,
                        const int pageX, const int pageY );
  void handleRightClick( const int x1, const int y1, const int modifier,
                         const int pageX, const int pageY );
  void handleDoubleTap( int x, int y );
  void handleDoubleLeftClick( Wt::WMouseEvent event );
  void handleDrag( int widgetStartX, int widgetStartY,
                   int widgetFinalX, int widgetFinalY,
                   int widgetOffsetLeft, int widgetOffsetTop,
                   Wt::WFlags<Wt::KeyboardModifier> modifiers,
                   Wt::WMouseEvent::Button button,
                   int dt );
  
  //handleArrowPress(...) only reacts to left and right arrow presses, and only
  //  if m_allowSingleClickHighlight==true.  If the shift key is pressed then
  //  bins on the left or right side of the region are added, otherwise the
  //  highlighted region is moved one bin to the left or right.  If multiple
  //  regions are highlighted, no action is taken.
  void handleArrowPress( const Wt::WKeyEvent &event );

  void handleControlMouseMove( int x_low, int x_high );
  void handleAltDrag( const double dx );
  
  
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  void setRenderWaitingStatus( const bool waiting ) const;
  bool checkHighBandWidthInterActionApplicable() const;
  void handleMouseDown( int x, int y, int modifiers );
  void handleMouseUp( int x, int y, int modifiers, int dt );
  void handleLeftMouseMove( int x, int y, int dt );
  void handleAltLeftMouseMove( int x, int y, int dt );
  void handleMouseLeft();
  void handleMouseEnter();
  void handlePinchZoomChange( int x0_t0, int x_t0, int x0_t1, int x_t1, int y );
#endif
  
  
protected:
  double m_widthInPixels;        //The width of the chart in pixels
  double m_heightInPixels;       //The height of the chart in pixels
  
  //Based just off the mouseWentUp() signal we dont know if the user is clicking
  //  or dragging.  So we will mark m_isDragging is true using the
  //  mouseDragged() signal
  bool m_isDragging;
  
  //xRangeChanged() signal will be emitted whenever the user does a 'drag'
  //  action when m_dragAction==ZoomIn.
  //If m_dragAction==HighLight then the xRangeChanged() will only be emmitted
  //  if it was with the left mouse button, and there where no modifier keys.
  DragAction m_dragAction;
  bool m_allowMultipleHighlightRegions; //only applicable when m_dragAction==HighLight
  bool m_allowSingleClickHighlight;     //only applicable when m_dragAction==HighLight
  
  //for all the bellow, the doubles are all the <x,y> coordinated of the action
  //  where x is in energy, and y is in counts.
  Wt::Signal<double,double> m_xRangeChanged;
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> m_controlKeyDragg;
  Wt::Signal<double,double> m_shiftKeyDragg;
  Wt::Signal<double,double> m_altKeyDragg;
  Wt::Signal<double,double> m_shiftAltKeyDragg;
  Wt::Signal<double,double> m_rightMouseDragg;
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> m_leftClick;
  Wt::Signal<double,double> m_doubleLeftClick;
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> m_rightClick;
  
  Wt::Signal<double/*start energy*/,double /*current energy*/> m_controlMouseMoved;
  
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
  
  Wt::Signal<> m_legendEnabled;
  Wt::Signal<> m_legendDisabled;
  
  CanvasForDragging *m_overlayCanvas;
  //If we use the OverlayCanvas, we might as well disconnect the
  //  mouseWentUp() and mouseDragged() signals of *this, so we have to keep
  //  track of these connections
#ifdef WT_USE_BOOST_SIGNALS2
  boost::signals2::connection m_arrowRespondConnection;
  boost::signals2::connection m_doubleClickConnection;
  boost::signals2::connection m_wtMouseWentUpConnection;
  boost::signals2::connection m_wtMouseDraggedConnection;
#else
  boost::signals::connection m_arrowRespondConnection;
  boost::signals::connection m_doubleClickConnection;
  boost::signals::connection m_wtMouseWentUpConnection;
  boost::signals::connection m_wtMouseDraggedConnection;
#endif
  
  std::unique_ptr<Wt::JSlot> m_alignOverlayCanvas;
  std::unique_ptr<Wt::JSlot> m_legendMovedSlot;
  
  // Used in case the user scrolls in the Anthony app
  // avoids overlapping the tabs
  int m_scrollY;
  
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
  
  /** On mobil avoid writing y-axis label text where the hamburger menu is. */
  bool m_avoidMobileMenu;
  
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  //Not all these variables are needed
  int m_mouseDownX, m_mouseDownY;
  int m_currentDownX, m_currentDownY;
  float m_mouseDownEnergy, m_mouseDownCounts;
  float m_mouseDownXMin, m_mouseDownXMax;
  float m_mouseDownYMin, m_mouseDownYMax;
#endif
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  ReferenceLineInfo m_referencePhotoPeakLines;
  std::vector<ReferenceLineInfo> m_persistedPhotoPeakLines;
#endif
};//class SpectrumChart




 #endif
