#include "InterSpec_config.h"

#include <memory>
#include <vector>
#include <utility>

#include <Wt/WServer>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/D3TimeChart.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace Wt;
using namespace std;

using SpecUtils::SpecFile;
using SpecUtils::Measurement;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
  const std::string &jsbool( bool val )
  {
    static const std::string t = "true";
    static const std::string f = "false";
    
    return val ? t : f;
  };
}//namespace



D3TimeChart::D3TimeChart( Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
  m_renderFlags( 0 ),
  m_layoutWidth( 0 ),
  m_layoutHeight( 0 ),
  m_chartWidthPx( 0.0 ),
  m_chartHeightPx( 0.0 ),
  m_compactXAxis( false ),
  m_showVerticalLines( false ),
  m_showHorizontalLines( false ),
  m_spec( nullptr ),
  m_highlights(),
  m_xAxisTitle( "Real Time (s)"),
  m_y1AxisTitle( "&gamma; counts"),
  m_y2AxisTitle( "n counts"),
  m_chartClicked( this ),
  m_chartDragged( this ),
  m_chartResized( this ),
  m_displayedXRangeChange( this ),
  m_chartClickedJS( nullptr ),
  m_chartDraggedJS( nullptr ),
  m_chartResizedJS( nullptr ),
  m_displayedXRangeChangeJS( nullptr ),
  m_jsgraph( jsRef() + ".chart" ),
  m_gammaLineColor( 0x00, 0x00, 0x00 ),
  m_neutronLineColor( 0x00, 0x00, 0x00 ),
  m_foregroundHighlightColor( 0x00, 0x00, 0x00 ),
  m_backgroundHighlightColor( 0x00, 0x00, 0x00 ),
  m_secondaryHighlightColor( 0x00, 0x00, 0x00 ),
  m_occupancyLineColor( 0x00, 0x00, 0x00 ),
  m_textColor( 0x00, 0x00, 0x00 ),
  m_axisColor( 0x00, 0x00, 0x00 ),
  m_chartMarginColor( 0x00, 0x00, 0x00 ),
  m_chartBackgroundColor( 0x00, 0x00, 0x00 )
{
  setLayoutSizeAware( true );
  addStyleClass( "D3TimeChart" );
    
  // Cancel right-click events for the div, we handle it all in JS
  setAttributeValue( "oncontextmenu",
                     "event.cancelBubble = true; event.returnValue = false; return false;" );
  
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/D3TimeChart.js" );
  wApp->useStyleSheet( "InterSpec_resources/D3TimeChart.css" );
  initChangeableCssRules();
}//D3TimeChart(...)


D3TimeChart::~D3TimeChart()
{
}


void D3TimeChart::defineJavaScript()
{
  string options = "{";
  options += "xtitle: '" + m_xAxisTitle + "'";
  options += ", y1title: '" + m_y1AxisTitle + "'";
  options += ", y2title: '" + m_y2AxisTitle + "'";
  options += ", compactXAxis: " + jsbool(m_compactXAxis);
  options += ", gridx: " + jsbool(m_showHorizontalLines);
  options += ", gridy: " + jsbool(m_showVerticalLines);
  options += ", chartLineWidth: 1.0";  //ToDo: Let this be specified in C+
  options += "}";
  
  setJavaScriptMember( "chart", "new D3TimeChart(" + jsRef() + "," + options + ");");
  setJavaScriptMember( "wtResize", "function(self, w, h, layout){" + m_jsgraph + ".handleResize();}" );
  
  setHighlightRegionsToClient();
  
  if( !m_chartClickedJS )
  {
    m_chartClickedJS.reset( new Wt::JSignal<int,int>(this, "timeclicked", false) );
    m_chartDraggedJS.reset( new Wt::JSignal<int,int,int>(this, "timedragged", false) );
    m_chartResizedJS.reset( new Wt::JSignal<double,double>(this, "timeresized", false ) );
    m_displayedXRangeChangeJS.reset( new Wt::JSignal<int,int,int>(this,"timerangechange",false) );
    
    m_chartClickedJS->connect( this, &D3TimeChart::chartClickedCallback );
    m_chartDraggedJS->connect( this, &D3TimeChart::chartDraggedCallback );
    m_chartResizedJS->connect( this, &D3TimeChart::chartResizedCallback );
    m_displayedXRangeChangeJS->connect( this, &D3TimeChart::displayedXRangeChangeCallback );
  }//if( !m_xRangeChangedJS )
  
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void D3TimeChart::doJavaScript( const std::string& js )
{
  if( isRendered() )
    WContainerWidget::doJavaScript( js );
  else
    m_pendingJs.push_back( std::move(js) );
}//doJavaScript(...)


void D3TimeChart::setData( std::shared_ptr<const SpecUtils::SpecFile> data )
{
  if( !data || data->sample_numbers().empty() )
  {
    if( m_spec )
      doJavaScript( m_jsgraph + ".setData(null);" );
    m_spec = data;
    return;
  }//if( no data to load );
  
  m_spec = data;
  
  /** Description of JSON format sent to client JS charting
   {
     // All arrays of numbers (realTimes, sampleNumbers, and various counts) will be the same length
     // Time in seconds of each time interval being plotted
     realTimes: [302.4,0.1,0.1,0.9,1,...],
   
     // The sample number for each time interval.  These usually be monotonically increasing, but
     //  this isnt garunteed, and the starting value isnt garunteed either.  These values link-up
     //  to the SampleNumber of the SpecFile the data is being loaded from.
     sampleNumbers: [-1,2,3,4,5,...],
   
     // The gamma counts to plot.  Usually we will only plot one line, but if it makes sense to
     //  plan ahead now, being able to plot multiple gamma lines would be useful (either stacked, or
     //  all lines just drawn on top of eachother); if there are multiple lines than the mouse going
     //  over a line could maybe pop-up the counts value for that time interval, and detector name
     //  (unless only a single detector, then name isnt necassary).
     gammaCounts: [{detName: 'Aa1', color: 'rgb(13,55,19)', counts: [1022,974,333,827,921,...]},
                   {detName: 'Ba1', color: 'green', counts: [55,60,18,99,1000,...]}],
   
     // The neutron counts to plot, similar to gamma, but the field may not exist, or be null, in
     //  which case there will be no neutron counts to plot, and the Y2-axis should disapear.
     neutronCounts: [{detName: 'Aa1n', color: '#145366', counts: [1,0,4,10,12,...]},
                     {detName: 'VD2n', color: '#ff0000', counts: [8,12,13,19,0,...]}]
   }
   */
  
  //Variable to control if we will only plot a single gamma and neutron line, or if we will plot
  //  each detector seperately.  This may become a user option at some point, or the idea of more
  //  than one line for gamma/neutron may get scrapped if it is to confusing or unhelpful.
  const bool plotDetectorsSeperate = false;
  
  vector<double> realTimes;
  vector<int> sampleNumbers;
  map<string,bool> hasGamma, hasNuetron;
  map<string,vector<double>> gammaCounts, neutronCounts;  //maps from detector name to counts
  
  const set<int> &sample_numbers = m_spec->sample_numbers();
  const vector<string> &detNames = m_spec->detector_names();

  
  for( const int sample_num : sample_numbers )
  {
    /// \TODO: check that all Measurements have same real times; right now we'll just use the max
    ///        value for each sample
    float realTime = 0.0f;
    sampleNumbers.push_back( sample_num );
    
    for( const string &detName : detNames )
    {
      const auto m = m_spec->measurement( sample_num, detName );
      
      if( !m )
      {
        gammaCounts[detName].push_back( 0.0 );
        neutronCounts[detName].push_back( 0.0 );
        
        if( !hasGamma.count(detName) )
          hasGamma[detName] = false;
        
        if( !hasNuetron.count(detName) )
          hasNuetron[detName] = false;
        
        continue;
      }//if( !m )
      
      realTime = std::max( realTime, m->real_time() );
      
      auto gammas = m->gamma_counts();
      if( gammas && gammas->size() )
        hasGamma[detName] = true;
      else if( !hasGamma.count(detName) )
        hasGamma[detName] = false;
               
      if( m->contained_neutron() )
        hasNuetron[detName] = true;
      else if( !hasNuetron.count(detName) )
        hasNuetron[detName] = false;
               
      gammaCounts[detName].push_back( m->gamma_count_sum() );
      neutronCounts[detName].push_back( m->neutron_counts_sum() );
    }//for( const string &name : m_spec->detector_names() )
    
    realTimes.push_back( realTime );
  }//for( loop over sample numbers )
  
  const size_t numSamples = sampleNumbers.size();
  //A quick sanity check all the arrays will be the same length.
  assert( numSamples == realTimes.size() );
  for( const auto &p : gammaCounts )
  {
    assert( p.second.size() == numSamples );
  }
  
  for( const auto &p : neutronCounts )
  {
    assert( p.second.size() == numSamples );
  }
  
  WStringStream js;
  js << m_jsgraph <<  ".setData( {\n";
  
  
  js << "\t\"realTimes\": [";
  for( size_t i = 0; i < numSamples; ++i )
    js << string(i ? "," : "") << realTimes[i];
  js << "],\n\t\"sampleNumbers\": [";
  for( size_t i = 0; i < numSamples; ++i )
    js << string(i ? "," : "") << sampleNumbers[i];
  js << "]";
  

  if( plotDetectorsSeperate )
  {
    int nwrote = 0;
    for( const auto &detName : detNames )
    {
      if( !hasGamma[detName] )
        continue;
    
      //For now we'll just make all detector lines the same color (if we are even doing multiple lines)
      js << string(nwrote++ ? "," : ",\n\t\"gammaCounts\": [" ) << "\n\t\t{\"detName\": \"" << detName << "\", \"color\": \""
         << (m_gammaLineColor.isDefault() ? string("#cfced2") :  m_gammaLineColor.cssText())
         << "\", \"counts\": [";
      
      const auto &counts = gammaCounts[detName];
      for( size_t i = 0; i < numSamples; ++i )
        js << string(i ? "," : "") << counts[i];
      js << "]}";
    }//for( const auto &detName : detNames )
    
    if( nwrote )
      js << "]";
    
    nwrote = 0;
    for( const auto &detName : detNames )
    {
      if( !hasNuetron[detName] )
        continue;
    
      //For now we'll just make all detector lines the same color (if we are even doing multiple lines)
      js << string(nwrote++ ? "," : ",\n\t\"neutronCounts\": [" ) << "\n\t\t{\"detName\": \"" << detName << "\", \"color\": \""
         << (m_neutronLineColor.isDefault() ? string("rgb(0,128,0)") :  m_neutronLineColor.cssText())
         << "\", \"counts\": [";

      const auto &counts = neutronCounts[detName];
      for( size_t i = 0; i < numSamples; ++i )
        js << string(i ? "," : "") << counts[i];
      js << "]}";
    }//for( const auto &detName : detNames )
    
    if( nwrote )
    js << "]";
  }else
  {
    bool haveAnyGamma = false, haveAnyNeutron = false;
    for( const auto &p : hasGamma )
      haveAnyGamma |= p.second;
    for( const auto &p : hasNuetron )
      haveAnyNeutron |= p.second;
      
    if( haveAnyGamma )
    {
      js << ",\n\t\"gammaCounts\": [{\"detName\": \"\", \"color\": \""
         << (m_gammaLineColor.isDefault() ? string("#cfced2") :  m_gammaLineColor.cssText())
         << "\", \"counts\": [";
      
      for( size_t i = 0; i < numSamples; ++i )
      {
        double sum = 0.0;
        for( const auto &p : gammaCounts )
          sum += p.second[i];
        js << string(i ? "," : "") << sum;
      }
      js << "]}]";
    }//if( haveAnyGamma )
    
    if( haveAnyNeutron )
    {
      js << ",\n\t\"neutronCounts\": [{\"detName\": \"\", \"color\": \""
         << (m_neutronLineColor.isDefault() ? string("#cfced2") :  m_neutronLineColor.cssText())
         << "\", \"counts\": [";
      
      for( size_t i = 0; i < numSamples; ++i )
      {
        double sum = 0.0;
        for( const auto &p : neutronCounts )
          sum += p.second[i];
        js << string(i ? "," : "") << sum;
      }
      js << "]}]";
    }//if( haveAnyNeutron )
  }//if( plotDetectorsSeperate ) / else
  
  js << "\n\t} );";
  
  cout << "\n\nWill set time data with JSON=" + js.str() + "\n\n" << endl;

  doJavaScript( js.str() );
}//void D3TimeChart::setData( std::vector<std::shared_ptr<const SpecUtils::SpecFile>> data )


void setHighlightedIntervals( const std::set<int> &sample_numbers, SpecUtils::SpectrumType type )
{
  //ToDo: implement
}


Wt::Signal<int,int> &D3TimeChart::chartClicked()
{
  return m_chartClicked;
}


Wt::Signal<int,int,int> &D3TimeChart::chartDragged()
{
  return m_chartDragged;
}


Wt::Signal<double,double> &D3TimeChart::chartResized()
{
  return m_chartResized;
}


Wt::Signal<int,int,int> &D3TimeChart::displayedXRangeChange()
{
  return m_displayedXRangeChange;
}


void D3TimeChart::setOccupancyStartAndStopSampleNumbers( const int first, const int last )
{
  // \TODO: implement
}


void D3TimeChart::setHighlightRegionsToClient()
{
  // \TODO: implement
}//setHighlightRegionsToClient()


void D3TimeChart::saveChartToPng( const std::string &filename )
{
  // \TODO: implement
}//void saveChartToPng( const std::string &filename )


void D3TimeChart::scheduleRenderAll()
{
  m_renderFlags |= TimeRenderActions::UpdateData;
  m_renderFlags |= TimeRenderActions::UpdateHighlightRegions;
  
  scheduleRender();
}//scheduleRenderAll()


void D3TimeChart::scheduleHighlightRegionRender()
{
  m_renderFlags |= TimeRenderActions::UpdateHighlightRegions;
  
  scheduleRender();
}//void scheduleHighlightRegionRender();


void D3TimeChart::applyColorTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return;
  
  setGammaLineColor( theme->timeChartGammaLine );
  setNeutronLineColor( theme->timeChartNeutronLine );
  setAxisLineColor( theme->timeAxisLines );
  setChartBackgroundColor( theme->timeChartBackground );
  setChartMarginColor( theme->timeChartMargins );
  setTextColor( theme->timeChartText );
  m_occupancyLineColor = theme->occupancyIndicatorLines;
  m_foregroundHighlightColor = theme->timeHistoryForegroundHighlight;
  m_backgroundHighlightColor = theme->timeHistoryBackgroundHighlight;
  m_secondaryHighlightColor = theme->timeHistorySecondaryHighlight;
  
  scheduleRenderAll();
}//void applyColorTheme( std::shared_ptr<const ColorTheme> theme );


void D3TimeChart::setGammaLineColor( const Wt::WColor &color )
{
  m_gammaLineColor = color.isDefault() ? WColor( 0x00, 0x00, 0x00 ) : color;
  scheduleRenderAll();
}

void D3TimeChart::setNeutronLineColor( const Wt::WColor &color )
{
  m_neutronLineColor = color.isDefault() ? WColor(0x00,0xff,0xff) : color;
  scheduleRenderAll();
}

void D3TimeChart::setTextColor( const Wt::WColor &color )
{
  m_textColor = color.isDefault() ? WColor(0,0,0) : color;
  const string c = m_textColor.cssText();
  
  const string rulename = "TimeTextColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ".D3TimeChart .xaxistitle,"
                                       " .D3TimeChart .yaxistitle,"
                                       " .D3TimeChart .yaxis,"
                                       " .D3TimeChart .yaxislabel,"
                                       " .D3TimeChart .xaxis",
                                       "stroke: " + c );
}//setTextColor(...)


void D3TimeChart::setAxisLineColor( const Wt::WColor &color )
{
  m_axisColor = color.isDefault() ? WColor(0,0,0) : color;
  
  string rulename = "TimeAxisColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ".D3TimeChart .xaxis > .domain,"
                                        " .D3TimeChart .yaxis > .domain,"
                                        " .D3TimeChart .xaxis > .tick > line,"
                                        " .D3TimeChart .yaxis > .tick,"
                                        " .D3TimeChart .yaxistick",
                                       "stroke: " + m_axisColor.cssText() );
}//setAxisLineColor(...)


void D3TimeChart::setChartMarginColor( const Wt::WColor &color )
{
  m_chartMarginColor = color;
  
  const string rulename = "TimeMarginColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  
  if( color.isDefault() )
  {
    if( m_cssRules.count(rulename) )
    {
      style.removeRule( m_cssRules[rulename] );
      m_cssRules.erase( rulename );
    }
    
    return;
  }//if( color.isDefault() )
  
  //Actualkly this will set the background for the entire chart...
  m_cssRules[rulename] = style.addRule( "#" + id() + " > svg", "background: " + color.cssText() );
}//setChartMarginColor(...)


void D3TimeChart::setChartBackgroundColor( const Wt::WColor &color )
{
  m_chartBackgroundColor = color;
  const string c = color.isDefault() ? "rgba(0,0,0,0)" : color.cssText();
  
  const string rulename = "TimeBackgroundColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  
  if( color.isDefault() )
  {
    if( m_cssRules.count(rulename) )
      style.removeRule( m_cssRules[rulename] );
  }//if( color.isDefault() )
  
  m_cssRules[rulename] = style.addRule( "#chartarea" + id(), "fill: " + c );
}//setChartBackgroundColor(...)


void D3TimeChart::setXAxisTitle( const std::string &title )
{
  m_xAxisTitle = title;
  SpecUtils::ireplace_all( m_xAxisTitle, "'", "&#39;" );
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle('" + m_xAxisTitle + "');" );
}//void setXAxisTitle( const std::string &title )


void D3TimeChart::setY1AxisTitle( const std::string &title )
{
  m_y1AxisTitle = title;
  SpecUtils::ireplace_all( m_y1AxisTitle, "'", "&#39;" );
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setY1AxisTitle('" + m_y1AxisTitle + "');" );
}//void setY1AxisTitle( const std::string &title )


void D3TimeChart::setY2AxisTitle( const std::string &title )
{
  m_y2AxisTitle = title;
  SpecUtils::ireplace_all( m_y2AxisTitle, "'", "&#39;" );
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setY2AxisTitle('" + m_y2AxisTitle + "');" );
}//void setY2AxisTitle( const std::string &title )


void D3TimeChart::layoutSizeChanged ( int width, int height )
{
  m_layoutWidth = width;
  m_layoutHeight = height;
}//void layoutSizeChanged ( int width, int height )


int D3TimeChart::layoutWidth() const
{
  return m_layoutWidth;
}


int D3TimeChart::layoutHeight() const
{
  return m_layoutHeight;
}


double D3TimeChart::chartWidthInPixels() const
{
  return m_chartWidthPx;
}


double D3TimeChart::chartHeightInPixels() const
{
  return m_chartHeightPx;
}


void D3TimeChart::setCompactAxis( const bool compact )
{
  m_compactXAxis = compact;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setCompactXAxis(" + jsbool(compact) + ");" );
}


bool D3TimeChart::isAxisCompacted() const
{
  return m_compactXAxis;
}

void D3TimeChart::showGridLines( bool show )
{
  m_showVerticalLines = show;
  m_showHorizontalLines = show;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridX(" + jsbool(show) + ");"
                  + m_jsgraph + ".setGridY(" + jsbool(show) + ");" );
}

void D3TimeChart::showVerticalLines( const bool draw )
{
  m_showVerticalLines = draw;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridX(" + jsbool(draw) + ");" );
}

void D3TimeChart::showHorizontalLines( const bool draw )
{
  m_showHorizontalLines = draw;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridY(" + jsbool(draw) + ");" );
}

bool D3TimeChart::verticalLinesShowing() const
{
  return m_showVerticalLines;
}

bool D3TimeChart::horizontalLinesShowing() const
{
  return m_showHorizontalLines;
}


void D3TimeChart::setXAxisRangeSamples( const int min_sample_num, const int max_sample_num )
{
  doJavaScript( m_jsgraph + ".setXAxisZoomSamples("
                + std::to_string(min_sample_num) + "," + std::to_string(max_sample_num) + ");" );
}


void D3TimeChart::initChangeableCssRules()
{
  WCssStyleSheet &style = wApp->styleSheet();

  m_cssRules["TimeGridColor"] = style.addRule( " .D3TimeChart .xgrid > .tick,  .D3TimeChart .ygrid > .tick", "stroke: #b3b3b3" );
  m_cssRules["TimeMinorGridColor"] = style.addRule( " .D3TimeChart .minorgrid", "stroke: #e6e6e6" );
}//void initChangeableCssRules()


void D3TimeChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
  
  if( m_renderFlags.testFlag(TimeRenderActions::UpdateData) )
    setData( m_spec );
  
  if( m_renderFlags.testFlag(TimeRenderActions::UpdateHighlightRegions) )
    setHighlightRegionsToClient();
  
  m_renderFlags = 0;
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void D3TimeChart::chartClickedCallback( int sample_number, int modifier_keys )
{
  m_chartClicked.emit(sample_number, modifier_keys);
}//chartClickedCallback(...)


void D3TimeChart::chartDraggedCallback( int first_sample_number, int last_sample_number, int modifier_keys )
{
  m_chartDragged.emit( first_sample_number, last_sample_number, modifier_keys );
}//chartDraggedCallback(...)


void D3TimeChart::chartResizedCallback( double chart_width_px, double chart_height_px )
{
  m_chartResized.emit( chart_width_px, chart_height_px );
}//chartResizedCallback(...)


void D3TimeChart::displayedXRangeChangeCallback( int first_sample_number, int last_sample_number, int samples_per_channel )
{
  m_displayedXRangeChange.emit( first_sample_number, last_sample_number, samples_per_channel );
}//displayedXRangeChangeCallback(...)
