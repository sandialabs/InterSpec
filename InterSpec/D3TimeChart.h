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
 - Handle case where background is like ~5 minutes long, but item of interest is only a few seconds.
 - Display multiple gamma lines corresponding to detectors in the system (stacked or overlaid); this
   can be ignored at first maybe.
 - Optional: be able to save as a PNG/JPEG/static-SVG.
 */


/** Wait until the chart is visible before loading JavaScript or CSS files, or defining the JS for this class.
 
 This relies on overloading `virtual void setHidden(bool,WAnimation&)` so that it loads
 and defines things only the first time it is set to not hidden.

 This saves about 200 KB (e.g., requires 1.5 MB, instead of 1.7 MB for full first-load, with a HPGe foreground+background state)
 loading because most common case is the time chart is initially hidden.

 As of 20240316, this has only been barely tested, so we wont enable until after releasing v1.0.12, so
 this way it can be well-tested for v1.0.13.
*/
#define OPTIMIZE_D3TimeChart_HIDDEN_LOAD 0


class D3TimeChart : public Wt::WContainerWidget
{
public:
  D3TimeChart( Wt::WContainerWidget *parent = nullptr );
  virtual ~D3TimeChart();
  
#if( OPTIMIZE_D3TimeChart_HIDDEN_LOAD )
  virtual void setHidden( bool hidden, const Wt::WAnimation& animation = Wt::WAnimation() );
#endif

  /** Set the spectrum file to display the time history for.
   @param data The data to be displayed.  May be nullptr to not display any data.
   @param det_to_display Detectors to use for displaying.  If empty, will display all detectors.  If any detectors specified are not
          in \p data, then an exception will be thrown.
   
   Will remove any existing highlighted intervals.
   */
  void setData( std::shared_ptr<const SpecUtils::SpecFile> data,
                std::vector<std::string> det_to_display );
  
  
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
  
  void setNeutronsHidden( const bool hide );
  bool neutronsHidden() const;
  
  void setGammaLogY( const bool logy );
  bool gammaLogY() const;
  
  
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
  
  /** Constructs and initializes everything. */
  void init();
  
#if( OPTIMIZE_D3TimeChart_HIDDEN_LOAD )
  bool m_inited;
#endif
  
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
  bool m_hideNeutrons;
  bool m_gammaLogY;
  
  std::shared_ptr<const SpecUtils::SpecFile> m_spec;
  std::vector<std::string> m_detectors_to_display;
  
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
  std::unique_ptr<Wt::JSignal<int,int,int,double,double,bool>>   m_displayedXRangeChangeJS;
  
  // Functions connected to the JSignal's
  void chartClickedCallback( int sample_number, int modifier_keys );
  void chartDraggedCallback( int first_sample_number, int last_sample_number, int modifier_keys );
  
  /** Callback to handle signal from JS whenever the displayed time-range is changed.
   @param first_sample_number The first sample number that is displayed (i.e., on far left of chart).
   @param last_sample_number The last sample number that is displayed (i.e., on far right of chart).
   @param compression_index Givens how many samples per display bin are used, via  `pow(2,compression_index)`
   @param start_time The starting time of the charted data (i.e., the left-most value of the chart x-axis); I think this is sum,med real-time.
   @param end_time The ending time of the charted data (i.e., the right-most value of the chart x-axis); I think this is sum,med real-time.
   @param is_user_action If this signal was emitted due to the user explicitly changing the time range, then will be true.  If the change is from setting the data, or being resized, or whatever, this will be false.
   */
  void displayedXRangeChangeCallback( const int first_sample_number,
                                     const int last_sample_number,
                                     const int compression_index,
                                     const double start_time,
                                     const double end_time,
                                     const bool is_user_action );
  
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
  
  /** We will track the displayed sample number range, as well as real-time range, for undo/redo purposes. */
  std::array<int,2> m_displayedSampleNumbers;
  std::array<double,2> m_displayedTimes;
  
  friend class D3TimeChartFilters;
};//class D3TimeChart_h


#endif //D3TimeChart_h
