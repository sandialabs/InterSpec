#ifndef SpectrumDisplayDiv_h
#define SpectrumDisplayDiv_h
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

#include <Wt/WEvent>
#include <Wt/WSignal>
#include <Wt/WContainerWidget>

#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
#include "InterSpec/ReferenceLineInfo.h"
#endif

//Forward declarations
class SpecMeas;
class PeakModel;
class InterSpec;

class SpectrumChart;
class SpectrumDataModel;
class CanvasForDragging;
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }

namespace Wt
{
  class WColor;
  class WGridLayout;
  class WApplication;
  class WButtonGroup;
  class WStackedWidget;
  class WContainerWidget;
}//namespace Wt


//
//2012-03-14: Instead of using a WBorderLayout to strech the charts <canvas>
//  element, which adds a whole bunch of tables and divs to the page, we can use
//  a simple customization of the javascript function wtResize to do this.
//  Note: I think in Newer Wt, this is no longer the case, so we could remove
//  this feature.
#define SpectrumDisplayDiv_USE_CUSTOM_LAYOUT 0


class SpectrumDisplayDiv : public Wt::WContainerWidget
{
public:
  SpectrumDisplayDiv( Wt::WContainerWidget *parent = 0 );
  virtual ~SpectrumDisplayDiv();

  //setTextInMiddleOfChart(...): draws some large text over the middle of the
  //  chart - used int the spectrum quizzer for text based questions.
  void setTextInMiddleOfChart( const Wt::WString &s );
  
  //setCompactAxis(): whether to slim down axis for small displays (e.g. on
  //  phone).  Note that effects wont be seen until next time chart is rendered.
  //  You should also adjust padding axis title text appropriately; x-axis
  //  padding of 23px seems to be a reasonable value.
  //Currently only effects x-axis.
  void setCompactAxis( const bool compact );
  bool isAxisCompacted() const;
  
  void setAvoidMobileMenu( const bool avoid );
  
  Wt::Signal<double/*keV*/,double/*counts*/,int/*pageX*/,int/*pageY*/> &chartClicked();
  Wt::Signal<double/*kev*/,double/*counts*/,int/*pageX*/,int/*pageY*/> &rightClicked();
  Wt::Signal<double/*keV*/,double/*counts*/> &doubleLeftClick();
  Wt::Signal<double/*keV start*/,double/*keV end*/,int/*last pageX*/,int/*last pageY*/> &controlKeyDragged();
  Wt::Signal<double/*keV start*/,double/*keV end*/> &shiftKeyDragged();

  
  void setPeakModel( PeakModel *model );
  
  void setData( std::shared_ptr<SpecUtils::Measurement> data_hist, const bool keep_curent_xrange );
  void setSecondData( std::shared_ptr<SpecUtils::Measurement> hist, const bool ownAxis );
  void setBackground( std::shared_ptr<SpecUtils::Measurement> background );


  // These 8 functions retrieve the corresponding info from the model.
  //std::shared_ptr<SpecUtils::Measurement> data();
  std::shared_ptr<const SpecUtils::Measurement> data()       const;
  //std::shared_ptr<SpecUtils::Measurement> secondData();
  std::shared_ptr<const SpecUtils::Measurement> secondData() const;
  //std::shared_ptr<SpecUtils::Measurement> background();
  std::shared_ptr<const SpecUtils::Measurement> background() const;

  float foregroundLiveTime() const;
  float foregroundRealTime() const;
  
  float backgroundLiveTime() const;
  float backgroundRealTime() const;

  float secondForegroundLiveTime() const;
  float secondForegroundRealTime() const;
  
  std::shared_ptr<const SpecUtils::Measurement> histUsedForXAxis() const;

  //displayScaleFactor():  This is the multiple 
  float displayScaleFactor( const SpecUtils::SpectrumType spectrum_type ) const;
  
  //setDisplayScaleFactor(): set the effective live time of 'spectrum_type'
  //  to be 'sf' timess the live time of 'spectrum_type'.
  void setDisplayScaleFactor( const float sf,
                              const SpecUtils::SpectrumType spectrum_type );
  
  
  void visibleRange( double &xmin, double &xmax,
                    double &ymin, double &ymax ) const;

  virtual void setPlotAreaPadding( int left, int top, int right, int bottom );
  virtual int plotAreaPadding( const Wt::Side side ) const;
  virtual void setXAxisTitle( const std::string &title );
  virtual void setYAxisTitle( const std::string &title );
  virtual void setY2AxisTitle( const std::string &title );
  const std::string xAxisTitle() const;
  const std::string yAxisTitle() const;
  const std::string y2AxisTitle() const;

  float xUnitsPerPixel() const;
  
  void setForegroundSpectrumColor( const Wt::WColor &color );
  void setBackgroundSpectrumColor( const Wt::WColor &color );
  void setSecondarySpectrumColor( const Wt::WColor &color );
  
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  
  void setChartMarginColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  
  
  
  void setDefaultPeakColor( const Wt::WColor &color );
  void setOccupiedTimeSamplesColor( const Wt::WColor &color );
  
  void setForegroundHighlightColor( const Wt::WColor &color );
  void setBackgroundHighlightColor( const Wt::WColor &color );
  void setSecondaryHighlightColor( const Wt::WColor &color );
  
  
  void clearTimeHighlightRegions( const SpecUtils::SpectrumType type );
  void setTimeHighLightRegions( const std::vector< std::pair<double,double> > &p,
                            const SpecUtils::SpectrumType type );
  bool removeDecorativeHighlightRegion( size_t regionid );
  size_t addDecorativeHighlightRegion( const float lowerx,
                                    const float upperx,
                                    const Wt::WColor &color );
  void clearOccupancyRegions();
  void setOccupancyRegions( const std::vector< std::pair<double,double> > &p );
  
  void setAutoAxisRange();

  void enableLegend( const bool forceMobileStyle );
  void disableLegend();
  bool legendIsEnabled() const;
  
  Wt::Signal<> &legendEnabled();
  Wt::Signal<> &legendDisabled();

  //Some functions to tell the SpectrumChart what kind of data it is displaying.
  //  This is currently only used for telling the units of the mouse on the
  //  overlay canvas.  If you do not set that this display is for either
  //  energy or time, then mouse coordanents will not be displayed.
  void setIsEnergyDisplay();
  void setIsTimeDisplay();
  bool isEnergyDisplay() const;
  bool isTimeDisplay() const;

  /** Executes appropriate javascript to generate and download a PNG based on
      the currently showing spectrum.  PNG generation is done client side.
   */
  void saveChartToPng( const std::string &name );
  
  //setControlDragDebouncePeriod(...): set how often the client will phone home
  //  to the server during a control-drag event.  Periods of <= 0 will result
  //  in it happening whenever the mouse is moved.
  void setControlDragDebouncePeriod( int milliseconds );
  
  //setControlDragContinuumPreview(...): draws a straight line from the bin
  //  contents at energy_start to energy_finish, but only if the overlay
  //  canvas is being used
  void setControlDragContinuumPreview( double energy_start, double x_finish );
  
  void showHistogramIntegralsInLegend( const bool show );

  void enableOverlayCanvas( bool outline, bool highlight, bool enalbeAltShiftHighlight );
  CanvasForDragging *overlayCanvas();
  void disableOverlayCanvas();
  bool overlayCanvasEnabled() const;
  void setOverlayCanvasVisible( bool visible );
  Wt::JSlot *alignOverlayCanvas();  //returns NULL if overlay canvas not enabled

  //setScrollingParent(...): sets the scrolling frame which contains the chart
  //  from which the overlay canvas should not extend beyond.  Calling this
  //  function when the overlay canvas is not enabled has no effect. Calling
  //  this function with a NULL argument removes this parent.
  void setScrollingParent( Wt::WContainerWidget *parent );

  void setScrollY( int scrollY );
  
  //overlayCanvasJsException() is mostly for debugging and will
  //  probably be removed in the future
  Wt::JSignal<std::string> *overlayCanvasJsException();  //returns NULL if not available

  //functions to connect/disconnect the Wt::Chart::WCartesianChart based
  //  mouse signals.  If you are using the overlay canvas, you can (but dont
  //  have to) disconnect these signals.  If you are not using the overlay
  //  canvas, than you should connect these signals.
  void connectWtMouseConnections();
  void disconnectWtMouseConnections();

  virtual void setHidden( bool hidden,
                         const Wt::WAnimation &animation = Wt::WAnimation() );

  void setMouseDragZooms(); // default action
  void setMouseDragHighlights( const bool allowMultiple = false,
                               const bool allowSingleClick = true );
  void allowArrowToMoveSingleClickRegion( bool allow = true );
  void disableMouseDragActions(); // signal will still be emitted
  Wt::Signal<double,double> &xRangeChanged();
  Wt::Signal<double,double> &rightMouseDragg();

  Wt::Signal<double,double> &altKeyDragged();
  Wt::Signal<double,double> &shiftAltKeyDragged();
  
  Wt::Signal<double,double> &controlMouseMoved();
  
  void prefferPngRenderForChart(); //default is <canvas>

  //By default SpectrumDisplayDiv has setLayoutSizeAware(true) set, so if the
  //  widget is being sized by a Wt layout manager, layoutWidth() and
  //  layoutHeight() will return this widget width and height respectively
  int layoutWidth() const;
  int layoutHeight() const;

  //For the case of auto-ranging x-axis, the below _may_ return 0 when auto
  //  range is set, but chart hasnt been rendered  (although maybe +-DBL_MAX)
  double xAxisMinimum() const;
  double xAxisMaximum() const;

  double yAxisMinimum() const;
  double yAxisMaximum() const;
  
  bool yAxisIsLog() const;
  void setYAxisLog( bool log );
  
  void showGridLines( bool show );  //shows horizantal and vertical
  void showVerticalLines( const bool draw );
  void showHorizontalLines( const bool draw );
  bool verticalLinesShowing() const;  // Added by christian (20170425)
  bool horizontalLinesShowing() const;
  
  bool backgroundSubtract() const;
  void setBackgroundSubtract( bool subtract );
  
  void setXAxisMinimum( const double minimum );
  void setXAxisMaximum( const double maximum );
  void setXAxisRange( const double minimum, const double maximum );

  void setYAxisMinimum( const double minimum );
  void setYAxisMaximum( const double maximum );
  void setYAxisRange( const double minimum, const double maximum );
  
  
  Wt::WPointF energyCountsToPixels( double energy, double counts ) const;
  

  double mapEnergyToXPixel( double energy ) const;
  double mapCountsToYPixel( double energy ) const;
  
  //peakLabel should be of type SpectrumChart::PeakLabels, but I didnt want
  //  to include the SpectrumChart header for just this
  void setShowPeakLabel( int peakLabel, bool show );
  bool showingPeakLabel( int peakLabel ) const;
  
#if( BUILD_AS_UNIT_TEST_SUITE || BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE || BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
  SpectrumDataModel *model(){ return m_model; }
  SpectrumChart *chart(){ return m_chart; }
#endif

  void setDisplayRebinFactor( const int factor );
  //relies on setLayoutSizeAware(true) being set.
  void guessAndUpdateDisplayRebinFactor();
  void setAutoAdjustDisplayRebinFactor( bool auto_rebin = true );
  
  int displayRebinFactor() const;
  
#if( IOS || ANDROID )
  //forceOverlayAlign(): does the JS for force an align of the overlay
  void forceOverlayAlign();
#endif
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc );
  void persistCurrentReferncePhotoPeakLines();
  void clearAllReferncePhotoPeakLines();
#endif

  //setShowRefLineInfoForMouseOver(): set wether or not the text information
  //  should be shown for the line that the mouse is currently over.  Default is
  //  to show the information.
  void setShowRefLineInfoForMouseOver( const bool show );
  
protected:

  void initUserTools();

  //chartXRangeChangedCallback(...): rebins the displayed data, and sets the
  //  y-axis to be auto-range
  void chartXRangeChangedCallback( double x, double y );

  //layoutSizeChanged(...): adjusts display binning if necessary
  virtual void layoutSizeChanged ( int width, int height );
  
  SpectrumChart *m_chart;
  SpectrumDataModel *m_model;

  int m_layoutWidth;
  int m_layoutHeight;
  bool m_autoAdjustDisplayBinnning;
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  friend class SpectrumViewerTester;
#endif
};//class SpectrumDisplayDiv


#endif
