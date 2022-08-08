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
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/InterSpec.h"  //Only for FeatureMarkerType::NumFeatureMarkers
#include "InterSpec/ReferenceLineInfo.h"


//Forward declarations
class SpecMeas;
class PeakModel;
class InterSpec;
struct ColorTheme;
class SpectrumDataModel;
namespace Wt
{
  class WCssTextRule;
}//namespace Wt
enum class FeatureMarkerType : int;
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }

/**
 ToDo:
 - When new foreground/background/secondary histogram is set, should set
   internal flag, call WWidget::scheduleRender(0), and not load the JS/JSON to
   client until D3SpectrumDisplayDiv::render() is called.  Same thing with
   colors and scale factors.
 - Get rid of SpectrumDataModel.  Will also require modifying PeakModel.
 - The y-axis range is not propogated from the client to server after many
   operations; this should be fixed.  Also, not sure if x-axis range is always
   propogated.
 - showHistogramIntegralsInLegend() is not implemented client side
 - setTextInMiddleOfChart() is not implemented client side
 - assign unique spectrum ID's in the JSON for each spectrum; convert JS to
   use this info.
 - Upgrade to use v5 of d3.js.
 */


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
  Wt::Signal<double/*keV start*/,double/*keV end*/> &shiftKeyDragged();
  
  /** When a previously existing ROI gets dragged by its edge, this signal will be emitted as it
   is being dragged, as well as when the mouse is finally let up.
   
   Note that by default #performExistingRoiEdgeDragWork is not hooked up to this signal, which you
   want to do for the primary spectrum display in InterSpec
   
   \sa performExistingRoiEdgeDragWork
   */
  Wt::Signal<double /*new roi lower energy*/,
             double /*new roi upper energy*/,
             double /*new roi lower px*/,
             double /*new roi upper px*/,
             double /*original roi lower energy*/,
             bool /*isFinalRange*/> &existingRoiEdgeDragUpdate();
  
  /** When a ROI is being created by holding the ctrl-key and dragging, this signal is emitted as
   as the user drags, and when the user lets up.
   
   Note that by default #performDragCreateRoiWork is not hooked up to the signal, which you want
   to do for the primary spectrum display in InterSpec
   
   \sa performDragCreateRoiWork
   */
  Wt::Signal<double /*lower energy*/,
             double /*upper energy*/,
             int    /*num peaks to force*/,
             bool /*isFinalRange*/,
             double /*window_xpx*/,
             double /*window_ypx*/> &dragCreateRoiUpdate();
  
  Wt::Signal<double,SpecUtils::SpectrumType> &yAxisScaled();
  
  
  /** Performs the work for the primary spectrum display in InterSpec that causes the peaks
   in an existing ROI to get re-fit as the user drags the edge.  To get this behavior, you
   need to hook this function up to the #existingRoiEdgeDragUpdate() signal.
   */
  void performExistingRoiEdgeDragWork( double new_lower_energy, double new_upper_energy,
                                      double new_lower_px, double new_upper_px,
                                      double original_lower_energy,
                                      bool isfinal );
  
  /** Does the the work for the primary spectrum display in InterSpec that lets you ctrl-drag
   to define a ROI with peaks in it.  To enable this functionality, you must hook this function
   up to the #dragCreateRoiUpdate() signal.
   
   \sa dragCreateRoiUpdate
   */
  void performDragCreateRoiWork( double lower_energy, double upper_energy,
                                int npeaks, bool isfinal,
                                double window_xpx, double window_ypx );
  
  
  
  
  
  void setPeakModel( PeakModel *model );
  
  void setData( std::shared_ptr<const SpecUtils::Measurement> data_hist, const bool keep_curent_xrange );
  void setSecondData( std::shared_ptr<const SpecUtils::Measurement> hist );
  void setBackground( std::shared_ptr<const SpecUtils::Measurement> background );
  
  void scheduleUpdateForeground();
  void scheduleUpdateBackground();
  void scheduleUpdateSecondData();
  
  
  /** Schedules the foreground peaks to be re-loaded to the client during the
   next call to #render (which Wt takes care of calling).
   */
  void scheduleForegroundPeakRedraw();
  
  /** Applies the current color theme.
   if nullptr, then sets to default colors.
   */
  void applyColorTheme( std::shared_ptr<const ColorTheme> theme );
  
  void setForegroundSpectrumColor( const Wt::WColor &color );
  void setBackgroundSpectrumColor( const Wt::WColor &color );
  void setSecondarySpectrumColor( const Wt::WColor &color );
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  void setChartMarginColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  void setDefaultPeakColor( const Wt::WColor &color );
  
  
  // These 3 functions retrieve the corresponding info from the model.
  std::shared_ptr<const SpecUtils::Measurement> data()       const;
  std::shared_ptr<const SpecUtils::Measurement> secondData() const;
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
  
  virtual void setXAxisTitle( const std::string &title );
  virtual void setYAxisTitle( const std::string &title );
  
  const std::string xAxisTitle() const;
  const std::string yAxisTitle() const;

  
  void enableLegend();
  void disableLegend();
  bool legendIsEnabled() const;
  
  Wt::Signal<> &legendEnabled();
  Wt::Signal<> &legendDisabled();
  
  void showHistogramIntegralsInLegend( const bool show );

  
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
  void removeAllDecorativeHighlightRegions();
  
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
  
  void setFeatureMarkerOption( FeatureMarkerType option, bool show );
  void setComptonPeakAngle( int angle );
  
  void showXAxisSliderChart( const bool show );
  bool xAxisSliderChartIsVisible() const;
  
  void showYAxisScalers( const bool show );
  bool yAxisScalersIsVisible() const;
  
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
  
    
  void setReferncePhotoPeakLines( const ReferenceLineInfo &nuc );
  void persistCurrentReferncePhotoPeakLines();
  void clearAllReferncePhotoPeakLines();
  void updateReferncePhotoPeakLines();
  
  //setShowRefLineInfoForMouseOver(): set wether or not the text information
  //  should be shown for the line that the mouse is currently over.  Default is
  //  to show the information.
  void setShowRefLineInfoForMouseOver( const bool show );
  
  /** Highlights a peak, at the specified energy, as if you had moused over it.
   
   Energy must match the peak mean to to places past the decimal; otherwise no action is taken.
   */
  void highlightPeakAtEnergy( const double energy );
  
  void updateRoiBeingDragged( const std::vector<std::shared_ptr<const PeakDef> > &roiBeingDragged );
  
  void removeAllPeaks();
  
  
  /** Executes appropriate javascript to generate and download a PNG or SVG based on
   the currently showing spectrum.  PNG or SVG generation is done client side.
   */
  void saveChartToImg( const std::string &name, const bool asPng );
  
  
protected:

  //updates the data JSON for the D3 spectrum on the JS side
  void renderForegroundToClient();
  void renderBackgroundToClient();
  void renderSecondDataToClient();
  
  
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
  
  void setForegroundPeaksToClient();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  /** Flags */
  enum D3RenderActions
  {
    UpdateForegroundPeaks = 0x01,
    
    UpdateForegroundSpectrum = 0x02,
    UpdateBackgroundSpectrum = 0x04,
    UpdateSecondarySpectrum = 0x08,
    
    ResetXDomain = 0x10,
    
    UpdateHighlightRegions = 0x20
    
    //ToDo: maybe add a few other things to this mechanism.
  };//enum D3RenderActions
  
  Wt::WFlags<D3RenderActions> m_renderFlags;
  
  //ToDo: should eliminate use of SpectrumDataModel in this class
  SpectrumDataModel *m_model;
  PeakModel *m_peakModel;
  
  bool m_compactAxis;
  bool m_legendEnabled;
  bool m_yAxisIsLog;
  bool m_backgroundSubtract;
  
  bool m_showVerticalLines;
  bool m_showHorizontalLines;
  bool m_showHistogramIntegralsInLegend;  //Not currently used/implemented
  
  bool m_showXAxisSliderChart;
  bool m_showYAxisScalers;
  
  std::vector<std::pair<double,double> > m_searchEnergies;
  std::vector<SpectrumChart::HighlightRegion> m_highlights;
  
  std::map<SpectrumChart::PeakLabels,bool> m_peakLabelsToShow;
  
  std::string m_xAxisTitle;
  std::string m_yAxisTitle;
  
  // JSignals
  //for all the bellow, the doubles are all the <x,y> coordinated of the action
  //  where x is in energy, and y is in counts.
  std::unique_ptr<Wt::JSignal<double, double> > m_shiftKeyDraggJS;
  std::unique_ptr<Wt::JSignal<double, double> > m_shiftAltKeyDraggJS;
  std::unique_ptr<Wt::JSignal<double, double> > m_rightMouseDraggJS;
  std::unique_ptr<Wt::JSignal<double, double> > m_doubleLeftClickJS;
  std::unique_ptr<Wt::JSignal<double,double,double/*pageX*/,double/*pageY*/> > m_leftClickJS;
  std::unique_ptr<Wt::JSignal<double,double,double/*pageX*/,double/*pageY*/> > m_rightClickJS;
  /** Currently including chart area in pixels in xRange changed from JS; this
      size in pixels is only approximate, since chart may not have been totally layed out
      and rendered when this signal was emmitted.
   ToDo: Should create dedicated signals for chart size in pixel, and also Y-range.
   */
  std::unique_ptr<Wt::JSignal<double,double,double,double> > m_xRangeChangedJS;
  std::unique_ptr<Wt::JSignal<double,double,double,double,double,bool> > m_existingRoiEdgeDragJS;
  std::unique_ptr<Wt::JSignal<double,double,int,bool,double,double> > m_dragCreateRoiJS;
  std::unique_ptr<Wt::JSignal<double,std::string> > m_yAxisDraggedJS;
  
  std::unique_ptr<Wt::JSignal<> > m_legendClosedJS;
  
  // Wt Signals
  //for all the bellow, the doubles are all the <x,y> coordinated of the action
  //  where x is in energy, and y is in counts.
  Wt::Signal<> m_legendEnabledSignal;
  Wt::Signal<> m_legendDisabledSignal;
  Wt::Signal<double/*xlow*/,double/*xhigh*/> m_xRangeChanged;
  Wt::Signal<double,double> m_controlKeyDragg;
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
             bool /*isFinalRange*/> m_existingRoiEdgeDrag;
  
  Wt::Signal<double /*lower energy*/,
             double /*upper energy*/,
             int /*force n peaks*/,
             bool /*isFinalRange*/,
             double /*window_xpx*/,
             double /*window_ypx*/> m_dragCreateRoi;
  
  Wt::Signal<double,SpecUtils::SpectrumType> m_yAxisScaled;
  
  // Signal Callbacks
  void chartShiftKeyDragCallback( double x0, double x1 );
  void chartShiftAltKeyDragCallback( double x0, double x1 );
  void chartRightMouseDragCallback( double x0, double x1 );
  void chartLeftClickCallback( double x, double y, double pageX, double pageY );
  void chartDoubleLeftClickCallback( double x, double y );
  void chartRightClickCallback( double x, double y, double pageX,
                                double pageY );
  
  void existingRoiEdgeDragCallback( double new_lower_energy, double new_upper_energy,
                                           double new_lower_px, double new_upper_px,
                                           double original_lower_energy,
                                           bool isfinal );
  
  void dragCreateRoiCallback( double lower_energy, double upper_energy,
                                int npeaks, bool isfinal,
                                double window_xpx, double window_ypx );
  
  void yAxisScaled( const double scale, const std::string &spectrum );
  
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
  
  ReferenceLineInfo m_referencePhotoPeakLines;
  std::vector<ReferenceLineInfo> m_persistedPhotoPeakLines;

  bool m_showRefLineInfoForMouseOver;
  
  bool m_showFeatureMarker[static_cast<int>(FeatureMarkerType::NumFeatureMarkers)];
  
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
