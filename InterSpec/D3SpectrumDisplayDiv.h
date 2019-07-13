#ifndef D3SpectrumDisplayDiv_h
#define D3SpectrumDisplayDiv_h

#include "InterSpec_config.h"

#include <map>
#include <memory>
#include <vector>
#include <utility>

#include <Wt/WColor>
#include <Wt/WEvent>
#include <Wt/WSignal>
#include <Wt/WCssStyleSheet>
#include <Wt/WContainerWidget>

#include <boost/thread.hpp>

#include "InterSpec/SpectrumChart.h"
#include "SpecUtils/SpectrumDataStructs.h"

static_assert( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE, "RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE must be enabled when USE_SPECTRUM_CHART_D3 is enabled" );

#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
#include "InterSpec/ReferenceLineInfo.h"
#endif

//Forward declarations
class SpecMeas;
class PeakModel;
class SpectrumChart;
class InterSpec;
class SpectrumDataModel;
class CanvasForDragging;
namespace Wt
{
  class WGridLayout;
  class WCssTextRule;
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


class D3SpectrumDisplayDiv : public Wt::WContainerWidget
{
public:
  D3SpectrumDisplayDiv( Wt::WContainerWidget *parent = 0 );
  virtual ~D3SpectrumDisplayDiv();
  
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
  
  Wt::Signal<double/*keV*/,double/*counts*/,int/*pageX*/,int/*pageY*/> &chartClicked();
  Wt::Signal<double/*kev*/,double/*counts*/,int/*pageX*/,int/*pageY*/> &rightClicked();
  Wt::Signal<double/*keV*/,double/*counts*/> &doubleLeftClick();
  Wt::Signal<double/*keV start*/,double/*keV end*/,int/*last pageX*/,int/*last pageY*/> &controlKeyDragged();
  Wt::Signal<double/*keV start*/,double/*keV end*/> &shiftKeyDragged();
  
  Wt::Signal<double /*new roi lower energy*/,
             double /*new roi upper energy*/,
             double /*new roi lower px*/,
             double /*new roi upper px*/,
             double /*original roi lower energy*/,
             bool /*isFinalRange*/> &roiDragUpdate();
  
  void setPeakModel( PeakModel *model );
  
  void setData( std::shared_ptr<Measurement> data_hist,
               float liveTime,
               float realTime,
               float neutronCounts,
               bool keep_curent_xrange );
  void setSecondData( std::shared_ptr<Measurement> hist,
                     float liveTime,
                     float realTime,
                     float neutronCounts,
                     bool ownAxis );
  void setBackground( std::shared_ptr<Measurement> background,
                     float liveTime,
                     float realTime,
                     float neutronCounts );
  
  //updateData(): updates the data JSON for the D3 spectrum on the JS side
  void updateData();
  void updateBackground();
  void updateSecondData();
  
  
  void setForegroundSpectrumColor( const Wt::WColor &color );
  void setBackgroundSpectrumColor( const Wt::WColor &color );
  void setSecondarySpectrumColor( const Wt::WColor &color );
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  void setChartMarginColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  void setDefaultPeakColor( const Wt::WColor &color );
  
  
  // These 8 functions retrieve the corresponding info from the model.
  std::shared_ptr<Measurement> data();
  std::shared_ptr<const Measurement> data()       const;
  std::shared_ptr<Measurement> secondData();
  std::shared_ptr<const Measurement> secondData() const;
  std::shared_ptr<Measurement> background();
  std::shared_ptr<const Measurement> background() const;
  
  float foregroundLiveTime() const;
  float foregroundRealTime() const;
  
  float backgroundLiveTime() const;
  float backgroundRealTime() const;
  
  float secondForegroundLiveTime() const;
  float secondForegroundRealTime() const;
  
  std::shared_ptr<const Measurement> histUsedForXAxis() const;
  
  //displayScaleFactor():  This is the multiple
  float displayScaleFactor( const SpectrumType spectrum_type ) const;
  
  //setDisplayScaleFactor(): set the effective live time of 'spectrum_type'
  //  to be 'sf' timess the live time of 'spectrum_type'.
  void setDisplayScaleFactor( const float sf,
                             const SpectrumType spectrum_type );
  
  
  void visibleRange( double &xmin, double &xmax,
                    double &ymin, double &ymax ) const;
  
  virtual void setPlotAreaPadding( int left, int top, int right, int bottom );
  virtual int plotAreaPadding( const Wt::Side side ) const;
  virtual void setXAxisTitle( const std::string &title );
  virtual void setYAxisTitle( const std::string &title );
  
  const std::string xAxisTitle() const;
  const std::string yAxisTitle() const;

  
  float xUnitsPerPixel() const;
  
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
  
  void showHistogramIntegralsInLegend( const bool show );

  
  virtual void setHidden( bool hidden,
                         const Wt::WAnimation &animation = Wt::WAnimation() );
  

  Wt::Signal<double,double> &xRangeChanged();
  Wt::Signal<double,double> &rightMouseDragg();
  
  Wt::Signal<double,double> &shiftAltKeyDragged();

  /** Set search energies.  First number is energy, second number is window;
     you can have how ever many pairs you want.  Calling with zero pairs removes
     all things on the chart.
   */
  void setSearchEnergies( const std::vector<std::pair<double,double>> &energy_windows );
  
  bool removeDecorativeHighlightRegion( size_t regionid );
  size_t addDecorativeHighlightRegion( const float lowerx,
                                      const float upperx,
                                      const Wt::WColor &color );
  
  
  //By default SpectrumDisplayDiv has setLayoutSizeAware(true) set, so if the
  //  widget is being sized by a Wt layout manager, layoutWidth() and
  //  layoutHeight() will return this widget width and height respectively
  int layoutWidth() const;
  int layoutHeight() const;
  
  //For the case of auto-ranging x-axis, the below _may_ return 0 when auto
  //  range is set, but chart hasnt been rendered  (although maybe +-DBL_MAX)
  double xAxisMinimum() const;
  double xAxisMaximum() const;
  
  double chartWidthInPixels() const;
  double chartHeightInPixels() const;
  
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
  
  void setFeatureMarkerOption( InterSpec::FeatureMarkerType option, bool show );
  void setComptonPeakAngle( int angle );
  
  
  void setXAxisMinimum( const double minimum );
  void setXAxisMaximum( const double maximum );
  void setXAxisRange( const double minimum, const double maximum );
  
  void setYAxisMinimum( const double minimum );
  void setYAxisMaximum( const double maximum );
  void setYAxisRange( const double minimum, const double maximum );
    
  //peakLabel should be of type SpectrumChart::PeakLabels, but I didnt want
  //  to include the SpectrumChart header for just this
  void setShowPeakLabel( int peakLabel, bool show );
  bool showingPeakLabel( int peakLabel ) const;
  
#if( BUILD_AS_UNIT_TEST_SUITE || BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE || BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
  SpectrumDataModel *model(){ return m_model; }
#endif
  
    
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc );
  void persistCurrentReferncePhotoPeakLines();
  void clearAllReferncePhotoPeakLines();
  void updateReferncePhotoPeakLines();
#endif
  
  //setShowRefLineInfoForMouseOver(): set wether or not the text information
  //  should be shown for the line that the mouse is currently over.  Default is
  //  to show the information.
  void setShowRefLineInfoForMouseOver( const bool show );
  
  void removeAllPeaks();
  
protected:
  void defineJavaScript();
  
  void initUserTools();
  
  /** In order to change some chart colors after intial load, we have to use
      WCssTextRule's, which dont seem to overide the text css files loaded, and
      using our own WCssStyleSheet (instead of WApplications) didnt seem to work
      out super easily on first try...
   */
  void initChangeableCssRules();
  
  /** Sets the highlight regions to client - currently unimplemented. */
  void setHighlightRegionsToClient();
  
  //layoutSizeChanged(...): adjusts display binning if necessary
  virtual void layoutSizeChanged ( int width, int height );
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  //ToDo: should eliminate use of SpectrumDataModel in this class
  SpectrumDataModel *m_model;
  PeakModel *m_peakModel;
  
  int m_layoutWidth;
  int m_layoutHeight;
  bool m_autoAdjustDisplayBinnning;
  
  bool m_compactAxis;
  bool m_legendEnabled;
  bool m_yAxisIsLog;
  bool m_backgroundSubtract;
  
  bool m_showVerticalLines;
  bool m_showHorizontalLines;
  bool m_showHistogramIntegralsInLegend;
  
  std::vector<std::pair<double,double> > m_searchEnergies;
  std::vector<SpectrumChart::HighlightRegion> m_highlights;
  
  std::map<SpectrumChart::PeakLabels,bool> m_peakLabelsToShow;
  
  std::string m_xAxisTitle;
  std::string m_yAxisTitle;
  
  // JSignals
  //for all the bellow, the doubles are all the <x,y> coordinated of the action
  //  where x is in energy, and y is in counts.
  boost::scoped_ptr<Wt::JSignal<double,double,int/*pageX*/,int/*pageY*/> > m_controlKeyDraggJS;
  boost::scoped_ptr<Wt::JSignal<double, double> > m_shiftKeyDraggJS;
  boost::scoped_ptr<Wt::JSignal<double, double> > m_shiftAltKeyDraggJS;
  boost::scoped_ptr<Wt::JSignal<double, double> > m_rightMouseDraggJS;
  boost::scoped_ptr<Wt::JSignal<double, double> > m_doubleLeftClickJS;
  boost::scoped_ptr<Wt::JSignal<double,double,int/*pageX*/,int/*pageY*/> > m_leftClickJS;
  boost::scoped_ptr<Wt::JSignal<double,double,int/*pageX*/,int/*pageY*/> > m_rightClickJS;
  /** Currently including chart area in pixels in xRange changed from JS; this
      size in pixels is only approximate, since chart may not have been totally layed out
      and rendered when this signal was emmitted.
   ToDo: Should create dedicated signals for chart size in pixel, and also Y-range.
   */
  boost::scoped_ptr<Wt::JSignal<double,double,double,double> > m_xRangeChangedJS;
  boost::scoped_ptr<Wt::JSignal<double,double,double,double,double,bool> > m_roiDraggedJS;
  
  boost::scoped_ptr<Wt::JSignal<> > m_legendClosedJS;
  
  // Wt Signals
  //for all the bellow, the doubles are all the <x,y> coordinated of the action
  //  where x is in energy, and y is in counts.
  Wt::Signal<> m_legendEnabledSignal;
  Wt::Signal<> m_legendDisabledSignal;
  Wt::Signal<double/*xlow*/,double/*xhigh*/> m_xRangeChanged;
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> m_controlKeyDragg;
  Wt::Signal<double,double> m_shiftKeyDragg;
  Wt::Signal<double,double> m_shiftAltKeyDragg;
  Wt::Signal<double,double> m_rightMouseDragg;
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> m_leftClick;
  Wt::Signal<double,double> m_doubleLeftClick;
  Wt::Signal<double,double,int/*pageX*/,int/*pageY*/> m_rightClick;
  
  Wt::Signal<double /*new roi lower energy*/,
             double /*new roi upper energy*/,
             double /*new roi lower px*/,
             double /*new roi upper px*/,
             double /*original roi lower energy*/,
             bool /*isFinalRange*/> m_roiDrag;
  
  
  // Signal Callbacks
  void chartControlKeyDragCallback( double x0, double x1, int pageX, int pageY );
  void chartShiftKeyDragCallback( double x0, double x1 );
  void chartShiftAltKeyDragCallback( double x0, double x1 );
  void chartRightMouseDragCallback( double x0, double x1 );
  void chartLeftClickCallback( double x, double y, int pageX, int pageY );
  void chartDoubleLeftClickCallback( double x, double y );
  void chartRightClickCallback( double x, double y, int pageX, int pageY );
  void chartRoiDragedCallback( double new_lower_energy, double new_upper_energy,
                               double new_lower_px, double new_upper_px,
                               double original_lower_energy,
                               bool isfinal );
  
  //chartXRangeChangedCallback(...): rebins the displayed data, and sets the
  //  y-axis to be auto-range
  void chartXRangeChangedCallback( double x, double y, double chart_width_px, double chart_height_px );
  
  /** The javascript variable name used to refer to the SpecrtumChartD3 object.
      Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;
  
  // X-axis and Y-axis values
  double m_xAxisMinimum;
  double m_xAxisMaximum;
  double m_yAxisMinimum;
  double m_yAxisMaximum;
  
  /** The width of the plotting area. */
  double m_chartWidthPx;
  double m_chartHeightPx;
  
  // Plot Area Padding values
  int m_plotAreaPaddingLeft;
  int m_plotAreaPaddingRight;
  int m_plotAreaPaddingTop;
  int m_plotAreaPaddingBottom;
  
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  ReferenceLineInfo m_referencePhotoPeakLines;
  std::vector<ReferenceLineInfo> m_persistedPhotoPeakLines;
#endif
  bool m_showRefLineInfoForMouseOver;
  
  bool m_showFeatureMarker[InterSpec::NumFeatureMarkers];
  
  /** Current compton angle used - note there is a bug in the WSpinBox inside
     WMenu (at least in Wt 3.3.4) so this value is likely not valid, and instead
     a JS only mehtod is used for setting this value client-side.
   */
  int m_comptonPeakAngle;
  
  Wt::WColor m_foregroundLineColor;
  Wt::WColor m_backgroundLineColor;
  Wt::WColor m_secondaryLineColor;
  Wt::WColor m_textColor;
  Wt::WColor m_axisColor;
  Wt::WColor m_chartMarginColor;
  Wt::WColor m_chartBackgroundColor;
  Wt::WColor m_defaultPeakColor;
  
  std::map<std::string,Wt::WCssTextRule *> m_cssRules;
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
     Note that not all calls to the D3 chart before Wt's rendering need to go
     here as they will be options set to the D3 chart during first rendering.
   */
  std::vector<std::string> m_pendingJs;
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  friend class SpectrumViewerTester;
#endif
};//class SpectrumDisplayDiv


#endif
