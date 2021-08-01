#ifndef D3TimeChart_h
#define D3TimeChart_h

#include "InterSpec_config.h"

#include <map>
#include <memory>
#include <vector>
#include <utility>

//#include <boost/optional.hpp>

#include <Wt/WColor>
#include <Wt/WEvent>
#include <Wt/WSignal>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpectrumChart.h"


//Forward declarations
struct ColorTheme;

namespace Wt
{
  class WCssTextRule;
}//namespace Wt

namespace SpecUtils
{
  class SpecFile;
  class Measurement;
  enum class SpectrumType : int;
}//namespace SpecUtils


class D3TimeChartFilters; //defined in D3TimeChart.cpp

/**
 Things to handle:
 - Display time series of gamma + neutron gross-count data.
 - Neutron and gamma use diff y-axis (neutron axis should disapear if no neutron info avalable).
 - Handle case where background is like ~5 minutes long, but item of interest is only a few seconds.
 - Option to have x-axis either be real time (fixing up issue with previous point), or have each
   sample take up the same number of pixels.
 - Handle case where there are more time-segments than pixels.
 - Indicate time segments being used for foreground, background, secondary spectrum (e.g., yellow,
   blue, etc fill).  Segments may not be continuos.
 - Indicate where vehicle occupancy begins and ends, and area designated as background or intrinsic.
 - Handle allowing user to select time regions to show spectrum for foreground, background,
   secondary (currently in InterSpec if user holds option when sleecting region, it will be used for
   background).  Regions may be discontinuous, and user may want to add/remove from currently used
   time regions (currently shift key while selecting will add to region, or remove if already
   selected).
 - Use d3.v3.min.js as used by SpecUtils.
 - Touch device compatible (details to be tested/worked out later).
 - Chart may be dynamically resized (e.g., user changes screen size).
 - X and Y axis units should be human-friendly numbers (e.g., time is at {1, 1.5, 2.0,...} and not
   {1.05, 1.55, 2.05, ..}, and similarly for y-axis counts.
 - X-axis should indicate start/end of time interval, not middle of time interval.
 - Should be histogram, not smooth graph.
 - Support displaying counts/time based on mouse position.
 - Make x-axis label compact and/or disapear.
 - Display multiple gamma lines corresponding to detectors in the system (stacked or overlaid); this
   can be ignored at first maybe.
 - Line, axis, background, time highlight region colors settable.
 - Adjust x-axis and y-axis ranges automatatically.  Default to have y-axis always go down to zero,
   and maybe have an option to adjust to minimum data height.
 - Show horizantal and vertical grid lines.
 - Adjustable padding around chart
 - Optional: be able to save as a PNG/JPEG/static-SVG.
 */


class D3TimeChart : public Wt::WContainerWidget
{
public:
  D3TimeChart( Wt::WContainerWidget *parent = nullptr );
  virtual ~D3TimeChart();
  
  /** Set the spectrum file to display the time history for.
   
   Will remove any existing highlighted intervals.
   */
  void setData( std::shared_ptr<const SpecUtils::SpecFile> data );
  
  
  void setHighlightedIntervals( const std::set<int> &sample_numbers,
                                const SpecUtils::SpectrumType type );
  
  void saveChartToPng( const std::string &filename );

  /** Signal when the user clicks on the chart.
   Gives the sample numebr user clicked on and a bitwise or of Wt::KeyboardModifiers.
   */
  Wt::Signal<int/*sample number*/,Wt::WFlags<Wt::KeyboardModifier>> &chartClicked();
  
  /** When the user drags on the chart to change the time range the spectrum is displayed for. */
  Wt::Signal<int/*start sample number*/,int/*end sample number*/,Wt::WFlags<Wt::KeyboardModifier>> &chartDragged();
  
  /**  Signal emitted when the displayed x-axis range changes via a user action; e.g., when zooming
   into or out of a region of interest.
   */
  Wt::Signal<int/*start sample number*/,int/*end sample number*/,int/*samples per channel*/> &displayedXRangeChange();
  
  /** Signal emitted when the user changes the gamma energy range that should be summed for to create this gross count chart. */
  //Wt::Signal<boost::optional<float>,boost::optional<float>> &energyRangeFilterChanged();
  
  
  
  static std::vector<std::pair<int,int>> sampleNumberRangesWithOccupancyStatus(
                                                const SpecUtils::OccupancyStatus status,
                                                std::shared_ptr<const SpecUtils::SpecFile> spec );
  
  /** Schedules (re)-rendering data + highlight regions. */
  void scheduleRenderAll();
  
  /** Schedules rendering the highlight regions. */
  void scheduleHighlightRegionRender();
  
  void applyColorTheme( std::shared_ptr<const ColorTheme> theme );
  
  void setGammaLineColor( const Wt::WColor &color );
  void setNeutronLineColor( const Wt::WColor &color );
  
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  void setChartMarginColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  
  void setXAxisTitle( const std::string &normalTitle, const std::string &compactTitle );
  void setY1AxisTitle( const std::string &title );
  void setY2AxisTitle( const std::string &title );
  
  void setCompactAxis( const bool compact );
  bool isAxisCompacted() const;
  
  void showGridLines( const bool draw );
  void showVerticalLines( const bool draw );
  void showHorizontalLines( const bool draw );
  bool verticalLinesShowing() const;
  bool horizontalLinesShowing() const;
  
  void setDontRebin( const bool dontRebin );
  bool dontRebin() const;
  
  void setXAxisRangeSamples( const int min_sample_num, const int max_sample_num );
  
  /** Returns the current user-entered gamma energy range that should be summed to create this gross count chart. */
  //std::pair<boost::optional<float>,boost::optional<float>> energyRangeFilters();
  
  
  /** Override WWebWidget::doJavaScript() to wait until this widget has been rendered before
   executing javascript so we can be sure all the JS objects we need are created.
   */
  virtual void doJavaScript( const std::string &js );
  
protected:
  
  void defineJavaScript();
  
  /** In order to change some chart colors after intial load, we have to use
      WCssTextRule's, which dont seem to overide the text css files loaded, and
      using our own WCssStyleSheet (instead of WApplications) didnt seem to work
      out super easily on first try...
   */
  void initChangeableCssRules();
  
  void setDataToClient();
  void setHighlightRegionsToClient();
  
  /** Shows or hides the user-selectable filters to control what the mouse/touch selects and energy range. */
  void showFilters( const bool show );
  
  /** The different modes users can interact with the chart.
   
   Note:
   - in all modes the mouse-wheel continues to zoom/pan the chart.
   -
   */
  enum class UserInteractionMode
  {
    /** The "normal" user mode where different modifier keys, or right mouse buttons do different things. */
    Default,
    
    /** Clicking and dragging with the left mouse button, or dragging with single finger will cause chart to zoom in/out (like the
     right-mouse button, or ctrl+left mouse normally does).
     */
    Zoom,
    
    /** Clicking and dragging, or tap-and-dragging causes the chart to pan left/right.*/
    Pan,
    
    /** Clicking and dragging with left-mouse button, or a single finger will cause the foreground/background/secondary spectrum
     to be defined by the drug region.  In this mode, holding the shift key will cause samples to be added to region.
     */
    SelectForeground, SelectBackground, SelectSecondary,
    
    /** Clicking and dragging with left-mouse button, or a single finger will cause the drug time region to be added to the
     foreground/background/secondary spectrum (wether or not shift key is held).
     */
    AddForeground, AddBackground, AddSecondary,
    
    /**
     
     */
    RemoveForeground, RemoveBackground, RemoveSecondary
  };//enum class UserSelectionMode
  
  /** Sets the current user interaction mode.*/
  void setUserInteractionMode( const UserInteractionMode mode );
  
  /** Called when the user changes energy range to display gross counts for. */
  //void userChangedEnergyRangeFilterCallback( const boost::optional<float> lowerEnergy,
  //                                           const boost::optional<float> upperEnergy );
  void userChangedEnergyRangeFilterCallback();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  /** Flags */
  enum TimeRenderActions
  {
    UpdateData = 0x01,
    
    UpdateHighlightRegions = 0x02,
    
    //ResetXDomain = 0x10
    
    //ToDo: maybe add a few other things to this mechanism.
  };//enum D3RenderActions
  
  Wt::WFlags<TimeRenderActions> m_renderFlags;
  
  /** The width, in pixels, of this entire widget. */
  int m_layoutWidth;
  /** The height, in pixels, of this entire widget. */
  int m_layoutHeight;
  
  /** The width of the plotting area in pixels. */
  double m_chartWidthPx;
  
  
  bool m_compactXAxis;
  bool m_showVerticalLines;
  bool m_showHorizontalLines;
  bool m_dontRebin;
  
  std::shared_ptr<const SpecUtils::SpecFile> m_spec;
  
  struct HighlightRegion
  {
    int start_sample_number;
    int end_sample_number;
    SpecUtils::SpectrumType type;
    Wt::WColor color;
  };//HighlightRegion
  
  
  std::vector<D3TimeChart::HighlightRegion> m_highlights;
  
  std::string m_xAxisTitle, m_compactXAxisTitle;
  std::string m_y1AxisTitle;
  std::string m_y2AxisTitle;
  
  // Signals to hook C++ code to, to be notified when a user action happens
  Wt::Signal<int/*sample number*/,Wt::WFlags<Wt::KeyboardModifier>> m_chartClicked;
  Wt::Signal<int/*start sample number*/,int/*end sample number*/,Wt::WFlags<Wt::KeyboardModifier>> m_chartDragged;
  Wt::Signal<int/*start sample number*/,int/*end sample number*/,int/*samples per channel*/> m_displayedXRangeChange;
  //Wt::Signal<boost::optional<float> /*lower keV*/,boost::optional<float>/*upper keV*/> m_energyRangeFilterChanged;
  
  // Signals called from JS to propogate infromation to the C++
  std::unique_ptr<Wt::JSignal<int,int>>       m_chartClickedJS;
  std::unique_ptr<Wt::JSignal<int,int,int>>   m_chartDraggedJS;
  std::unique_ptr<Wt::JSignal<int,int,int>>   m_displayedXRangeChangeJS;
  
  // Functions connected to the JSignal's
  void chartClickedCallback( int sample_number, int modifier_keys );
  void chartDraggedCallback( int first_sample_number, int last_sample_number, int modifier_keys );
  void displayedXRangeChangeCallback( int first_sample_number, int last_sample_number, int samples_per_channel );
  
  /** The javascript variable name used to refer to the SpecrtumChartD3 object.
      Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;

  Wt::WContainerWidget *m_chart;
  D3TimeChartFilters *m_options;
  Wt::WContainerWidget *m_showOptionsIcon;
  
  Wt::WColor m_gammaLineColor;
  Wt::WColor m_neutronLineColor;
  Wt::WColor m_foregroundHighlightColor;
  Wt::WColor m_backgroundHighlightColor;
  Wt::WColor m_secondaryHighlightColor;
  Wt::WColor m_occLineColor;
  Wt::WColor m_textColor;
  Wt::WColor m_axisColor;
  Wt::WColor m_chartMarginColor;
  Wt::WColor m_chartBackgroundColor;
  
  std::map<std::string,Wt::WCssTextRule *> m_cssRules;
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
     Note that not all calls to the D3 chart before Wt's rendering need to go
     here as they will be options set to the D3 chart during first rendering.
   */
  std::vector<std::string> m_pendingJs;
  
  friend class D3TimeChartFilters;
};//class D3TimeChart_h


#endif //D3TimeChart_h
