#include "InterSpec_config.h"

#include <tuple>
#include <limits>
#include <memory>
#include <vector>
#include <utility>

#include <Wt/WServer>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/D3TimeChart.h"

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
  m_occLineColor(),
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
  // Do we need to cleanup any JS objects here?
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
    if( m_highlights.size() )
    {
      m_highlights.clear();
      scheduleHighlightRegionRender();
    }
    
    if( m_spec )
      doJavaScript( m_jsgraph + ".setData(null);" );
    
    m_spec = data;
    return;
  }//if( no data to load );
  
  if( !m_highlights.empty() )
  {
    m_highlights.clear();
    scheduleHighlightRegionRender();
  }//if( !m_highlights.empty() )
  
  if( m_spec != data )
    scheduleRenderAll();
  
  m_spec = data;
}//void setData(...)
  


void D3TimeChart::setDataToClient()
{
  if( !m_spec )
  {
    doJavaScript( m_jsgraph +  ".setData( null );" );
    return;
  }//if( !m_spec )
  
  /** Description of JSON format sent to client JS charting
   {
     // All arrays of numbers (realTimes, sampleNumbers, and various counts) will be the same length
     // Time in seconds of each time interval being plotted.  These numbers specify the length of
     //  each time interval on the x-axis.
     realTimes: [302.4,0.1,0.1,0.9,1,...],
   
     
     // The start times of each sample is represented as the number of milliseconds since the
     //  UNIX epoch, like in JavaScript, however, to save space we will define a 'startTimeOffset'
     //  that will need to be added to each of the actual 'startTimes' to get the final time.  Also,
     //  note that JS Date assumes UTC time, so make sure when you display the time, its in UTC so
     //  no timezone offsets get added to it
     //  'startTimeOffset' and 'startTimes' may not be present if start times are not known.  Also,
     //  entries in 'startTimes' may be null if start time wasnt available for that sample.
     startTimeOffset: 1264291704078,
     startTimes: [0,604800,100,100,900, ...],
   
     // The sample number for each time interval.  These usually be monotonically increasing, but
     //  this isnt garunteed, and the starting value isnt garunteed either.  These values link-up
     //  to the SampleNumber of the SpecFile the data is being loaded from.
     sampleNumbers: [-1,2,3,4,5,...],
   
     // Array, that if present, will be same length as reaTimes and sampleNumbers, and gives the
     //  source-type of each sample.  The values in this array coorespond to the numberical values
     //  of the SpecUtils::SourceType enum, specifically:
     //    IntrinsicActivity : 0, Calibration : 1, Background : 2, Foreground : 3, Unknown : 4.
     sourceTypes: [2,3,3,3,2,...],
   
     // The GPS coordinates of the timeslice.
     //  This array will be missing if no GPS goordinates are available, and if there isnt
     //  coordinates for a given sample, null will be provided.
     //  In the future may add timestamp as third element of array.
     gpsCoordinates: [[37.675821,-121.709863],[37.675821,-121.709863],null,[37.675821,-121.709863],
                      [37.675821,-121.709863],...],
   
     // The gamma counts to plot.  Usually we will only plot one line, but if it makes sense to
     //  plan ahead now, being able to plot multiple gamma lines would be useful (either stacked, or
     //  all lines just drawn on top of eachother); if there are multiple lines than the mouse going
     //  over a line could maybe pop-up the counts value for that time interval, and detector name
     //  (unless only a single detector, then name isnt necassary).
     //  The live-time array may not be present; if it isnt, then you can assume live-times are
     //  equal to real times, and use the realTimes array.  If live-time array is present, it will
     //  be the same length as the counts array; use the live-time to divide the counts by to get
     //  counts per second for each time interval (live-time is typically a little less than
     //  real-time, but may be up to like ~90% less in extreme situations).
     //  Entries in 'counts' and/or 'liveTimes' arrays may be null if that information is not
     //  available for the corresponding sample number.
     gammaCounts: [{detName: 'Aa1', color: 'rgb(13,55,19)', counts: [1022,974,333,827,921,...],
                    liveTimes: [300,0.1,0.1,0.99,0.9,...]},
                   {detName: 'Ba1', color: 'green', counts: [55,60,18,99,1000,...],
                    liveTimes: [299,0.092,0.1,1,1,...]}
                  ],
   
     // The neutron counts to plot, similar to gamma, but the field may not exist, or be null, in
     //  which case there will be no neutron counts to plot, and the Y2-axis should disapear.
     //  Notes under gammaCounts above for the live-times array also apply here
     neutronCounts: [{detName: 'Aa1n', color: '#145366', counts: [1,0,4,10,12,...],
                     liveTimes: [300,0.1,0.1,0.1,0.1,...]},
                     {detName: 'VD2n', color: '#ff0000', counts: [8,12,13,19,0,...]}],
   
     // The sample number ranges for "occupied" and non-occupied times of the data.
     //  This field may not be present, may be null, may be empty array, or have entries as
     //  described.  This field is typically used for describing radiation portal monitor data, and
     //  will be absent for portable detection systems.
     occupancies: [{
                    // Says whether this sample range cooresponds to being "occupied" (cooresponding
                    //  to (SpecUtils::OccupancyStatus::Occupied) or if false, "not occupied"
                    //  cooresponding to SpecUtils::OccupancyStatus::NotOccupied).  If occupancy
                    //  status for a sample is unknown or undefined, then that sample will not be
                    //  included in any of these ranges.
                    //  Typically non-occupied intervals will also be marked as Background in the
                    //  'sourceTypes' array, and occupied intervals will be marked as Foreground.
                    status: true,
   
                    // The sample number cooresponding to the first occupied sample
                    // In the previous implementation a vertical line was drawn, on left side of
                    // sample with the text "occ. start"
                    startSample: 4,
        
                    // The last occupied sample number.  Was previously drawn as a vertical line
                    //  with text "occ. end"
                    endSample: 10,
                   
                    // Color to draw vertical line and text
                    color: 'grey'
                   },
                   {...}
                  ]
   }
   */
  
  //Variable to control if we will only plot a single gamma and neutron line, or if we will plot
  //  each detector seperately.  This may become a user option at some point, or the idea of more
  //  than one line for gamma/neutron may get scrapped if it is to confusing or unhelpful.
  const bool plotDetectorsSeperate = false;
  
  vector<double> realTimes;
  vector<int> sampleNumbers;
  vector<SpecUtils::SourceType> sourceTypes;
  vector<std::tuple<double,double,boost::posix_time::ptime>> gpsCoordinates;
  int64_t startTimesOffset = 0;
  vector<int64_t> startTimes;
  bool allUnknownSourceType = true, anyStartTimeKnown = false, haveAnyGps = false;
  
  map<string,bool> hasGamma, hasNuetron;
  map<string,vector<double>> gammaCounts, neutronCounts, liveTimes;  //maps from detector name to counts
  
  const set<int> &sample_numbers = m_spec->sample_numbers();
  const vector<string> &detNames = m_spec->detector_names();

#define Q_DBL_NaN std::numeric_limits<double>::quiet_NaN()
  
  for( const int sample_num : sample_numbers )
  {
    /// \TODO: check that all Measurements have same real times; right now we'll just use the max
    ///        value for each sample
    float realTime = 0.0f;
    boost::posix_time::ptime startTime;
    SpecUtils::SourceType sourcetype = SpecUtils::SourceType::Unknown;
    tuple<double,double,boost::posix_time::ptime> coords{ Q_DBL_NaN, Q_DBL_NaN, {} };
    
    sampleNumbers.push_back( sample_num );
    
    for( const string &detName : detNames )
    {
      const auto m = m_spec->measurement( sample_num, detName );
      
      if( !m )
      {
        liveTimes[detName].push_back( Q_DBL_NaN );
        gammaCounts[detName].push_back( Q_DBL_NaN );
        neutronCounts[detName].push_back( Q_DBL_NaN );
        
        if( !hasGamma.count(detName) )
          hasGamma[detName] = false;
        
        if( !hasNuetron.count(detName) )
          hasNuetron[detName] = false;
        
        continue;
      }//if( !m )
      
      realTime = std::max( realTime, m->real_time() );
      if( startTime.is_special() )
        startTime = m->start_time();
      
      auto gammas = m->gamma_counts();
      if( gammas && gammas->size() )
        hasGamma[detName] = true;
      else if( !hasGamma.count(detName) )
        hasGamma[detName] = false;
               
      if( m->contained_neutron() )
        hasNuetron[detName] = true;
      else if( !hasNuetron.count(detName) )
        hasNuetron[detName] = false;
               
      if( m->source_type() != SpecUtils::SourceType::Unknown )
        sourcetype = m->source_type();
      
      gammaCounts[detName].push_back( m->gamma_count_sum() );
      neutronCounts[detName].push_back( m->neutron_counts_sum() );
      liveTimes[detName].push_back( m->live_time() );
      
      if( m->has_gps_info() )
      {
        haveAnyGps = true;
        std::get<0>(coords) = m->longitude();
        std::get<1>(coords) = m->latitude();
        std::get<2>(coords) = m->position_time();
      }
    }//for( const string &name : m_spec->detector_names() )
    
    realTimes.push_back( realTime );
    sourceTypes.push_back( sourcetype );
    allUnknownSourceType = (allUnknownSourceType && (sourcetype == SpecUtils::SourceType::Unknown));
    
    if( startTime.is_special() )
    {
      startTimes.push_back( std::numeric_limits<int64_t>::min() );
    }else
    {
      anyStartTimeKnown = true;
       
      static const boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
      
      // std::time_t (what boost uses for seconds type) may be 32-bit, so convert to int64_t to
      //  avoid overflow.
      int64_t nmilli = (startTime - epoch).total_seconds();
      nmilli *= 1000;
      
      const double frac = startTime.time_of_day().fractional_seconds()
                        / static_cast<double>(boost::posix_time::time_duration::ticks_per_second());
      nmilli += static_cast<int64_t>( std::round(1000*frac) );
      
      if( startTimesOffset <= 0 && nmilli > 0 )
        startTimesOffset = nmilli;
      
      startTimes.push_back( nmilli - startTimesOffset );
    }
    
    gpsCoordinates.push_back( coords );
  }//for( loop over sample numbers )
  
  const size_t numSamples = sampleNumbers.size();
  //A quick sanity check all the arrays will be the same length.
  assert( numSamples == realTimes.size() );
  assert( numSamples == sourceTypes.size() );
  assert( numSamples == startTimes.size() );
  assert( numSamples == gpsCoordinates.size() );
  
  for( const auto &p : gammaCounts )
  {
    assert( p.second.size() == numSamples );
  }
  
  for( const auto &p : neutronCounts )
  {
    assert( p.second.size() == numSamples );
  }
  
  for( const auto &p : liveTimes )
  {
    assert( p.second.size() == numSamples );
  }
  
  if( !numSamples )
  {
    doJavaScript( m_jsgraph +  ".setData( null );" );
    return;
  }//if( !numSamples )
  
  const char jssep[] = { '[', ',' };
  
  auto printNumberArray = [&jssep]( WStringStream &js, const vector<double> &arr ){
    
    if( arr.empty() )
    {
      js << "[]";
      return;
    }//if( !arr.empty() )
    
    for( size_t i = 0; i < arr.size(); ++i )
    {
      if( IsNan(arr[i]) )
        js << jssep[i ? 1 : 0] << "null";
      else
        js << jssep[i ? 1 : 0] << arr[i];
    }//for( loop over realTimes )
    
    js << "]";
  };//printNumberArray
  
  
  WStringStream js;
  js << m_jsgraph <<  ".setData( {\n";
  
  js << "\t\"realTimes\": ";
  printNumberArray( js, realTimes );
  
  if( anyStartTimeKnown )
  {
    //doubles have 53 bits of integer precision - should be fine
    js << ",\n\t\"startTimeOffset\": " << startTimesOffset;
    js << ",\n\t\"startTimes\": ";
    for( size_t i = 0; i < startTimes.size(); ++i )
    {
      js << jssep[i ? 1 : 0];
      if( startTimes[i] == std::numeric_limits<int64_t>::min() )
        js << "null";
      else
        js << startTimes[i];
    }
    js << "]";
  }//if( anyStartTimeKnown )
  
  js << ",\n\t\"sampleNumbers\": ";
  for( size_t i = 0; i < sampleNumbers.size(); ++i )
    js << jssep[i ? 1 : 0] << sampleNumbers[i];
  js << "]";
  
  if( !allUnknownSourceType )
  {
    js << ",\n\t\"sourceTypes\": ";
    for( size_t i = 0; i < sourceTypes.size(); ++i )
      js << jssep[i ? 1 : 0] << static_cast<int>(sourceTypes[i]);
    js << "]";
  }//if( !allUnknownSourceType )
  
  
  
  if( haveAnyGps )
  {
    js << ",\n\t\"gpsCoordinates\": ";
    for( size_t i = 0; i < gpsCoordinates.size(); ++i )
    {
      js << jssep[i ? 1 : 0];
      
      const auto &coords = gpsCoordinates[i];
      if( IsNan(std::get<0>(coords)) )
        js << "null";
      else
        js << "[" << std::get<0>(coords) << "," << std::get<1>(coords) << "]";
    }//for( size_t i = 0; i < gpsCoordinates.size(); ++i )
    
    js << "]";
  }//if( haveAnyGps )
  
  
  if( plotDetectorsSeperate )
  {
    int nwrote = 0;
    for( const auto &detName : detNames )
    {
      if( !hasGamma[detName] )
        continue;
    
      //For now we'll just make all detector lines the same color (if we are even doing multiple lines)
      js << string(nwrote++ ? "," : ",\n\t\"gammaCounts\": [" )
         << "\n\t\t{\"detName\": \"" << detName << "\", \"color\": \""
         << (m_gammaLineColor.isDefault() ? string("#cfced2") :  m_gammaLineColor.cssText())
         << "\", \"counts\": ";
      
      const auto &counts = gammaCounts[detName];
      printNumberArray( js, counts );
      
      
      const auto &lts = liveTimes[detName];
      js << ",\n\t\"liveTimes\": ";
      printNumberArray( js, lts );
      
      js << "}";
    }//for( const auto &detName : detNames )
    
    if( nwrote )
      js << "]";
    
    nwrote = 0;
    for( const auto &detName : detNames )
    {
      if( !hasNuetron[detName] )
        continue;
    
      //For now we'll just make all detector lines the same color (if we are even doing multiple lines)
      js << string(nwrote++ ? "," : ",\n\t\"neutronCounts\": [" ) << "\n\t\t{\"detName\": \""
         << detName << "\", \"color\": \""
         << (m_neutronLineColor.isDefault() ? string("rgb(0,128,0)") :  m_neutronLineColor.cssText())
         << "\", \"counts\": ";

      const auto &counts = neutronCounts[detName];
      printNumberArray( js, counts );
      
      js << "}";
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
         << "\",\n\t\t\"counts\": [";
      
      for( size_t i = 0; i < numSamples; ++i )
      {
        double sum = 0.0;
        bool anyNonNan = false;
        for( const auto &p : gammaCounts )
        {
          assert( p.second.size() == numSamples );
          
          if( !IsNan(p.second[i]) )
          {
            anyNonNan = true;
            sum += p.second[i];
          }
        }//for( const auto &p : gammaCounts )
        
        double liveTime = 0.0;
        bool haveLiveTime = false;
        for( const auto &p : liveTimes )
        {
          if( !IsNan(p.second[i]) )
          {
            haveLiveTime = true;
            liveTime += p.second[i];
          }
        }
        
        
        js << string(i ? "," : "");
        if( anyNonNan )
          js << sum;
        else
          js << "null";
      }//for( size_t i = 0; i < numSamples; ++i )
      js << "]";
      
      js << ",\n\t\t\"liveTimes\": [";
      for( size_t i = 0; i < numSamples; ++i )
      {
        double liveTime = 0.0;
        bool haveLiveTime = false;
        for( const auto &p : liveTimes )
        {
          if( !IsNan(p.second[i]) )
          {
            haveLiveTime = true;
            liveTime += p.second[i];
          }
        }
        
        js << string(i ? "," : "");
        if( haveLiveTime )
          js << liveTime;
        else
          js << "null";
      }//for( size_t i = 0; i < numSamples; ++i )
      js << "]";
      
      js << "\n\t}]";
    }//if( haveAnyGamma )
    
    if( haveAnyNeutron )
    {
      js << ",\n\t\"neutronCounts\": [{\"detName\": \"\", \"color\": \""
         << (m_neutronLineColor.isDefault() ? string("#cfced2") :  m_neutronLineColor.cssText())
         << "\", \"counts\": [";
      
      vector<double> neutronLiveTimes( numSamples, std::numeric_limits<double>::quiet_NaN() );
      for( size_t i = 0; i < numSamples; ++i )
      {
        bool haveNonNan = false;
        double sum = 0.0;
        for( const auto &p : neutronCounts )
        {
          if( !IsNan(p.second[i]) )
          {
            haveNonNan = true;
            sum += p.second[i];
            
            if( IsNan(neutronLiveTimes[i]) )
              neutronLiveTimes[i] = 0;
            if( !IsNan(realTimes[i]) )
              neutronLiveTimes[i] += realTimes[i];
          }//if( we have neutron counts )
        }//for( loop over neutron detectors to sum their counts )
        js << string(i ? "," : "");
        if( haveNonNan )
          js << sum;
        else
          js << "null";
      }//for( loop over samples )
      
      js << "],\n\t\t\"liveTimes\": ";
      printNumberArray( js, neutronLiveTimes );
      
      js << "\n\t}]";
    }//if( haveAnyNeutron )
  }//if( plotDetectorsSeperate ) / else
  
  auto occRanges = sampleNumberRangesWithOccupancyStatus( SpecUtils::OccupancyStatus::Occupied, m_spec );
  if( occRanges.size() )
  {
    js << ",\n\t\"occupancies\": [";
    for( size_t i = 0; i < occRanges.size(); ++i )
    {
      js << string(i ? ",\n\t\t" : "\n\t\t")
         << "{ startSample: " << occRanges[i].first
         << ", endSample: " << occRanges[i].second
         << ", color: \""
         << (m_occLineColor.isDefault() ? string("rgb(128,128,128)") :  m_occLineColor.cssText())
         << "\" }";
    }//for( size_t i = 0; i < occRanges.size(); ++i )
    js << "\n\t]";
  }//if( occRanges.size() )
  
  js << "\n\t} );";
  
  cout << "\n\nWill set time data with JSON=" + js.str() + "\n\n" << endl;

  doJavaScript( js.str() );
}//void setDataToClient()


void D3TimeChart::setHighlightedIntervals( const std::set<int> &sample_numbers,
                                           const SpecUtils::SpectrumType type )
{
  m_highlights.erase( std::remove_if( begin(m_highlights), end(m_highlights),
    [type](const D3TimeChart::HighlightRegion &region) -> bool {
      return (region.type == type);
  }), end(m_highlights) );
 
  if( sample_numbers.empty() )
  {
    scheduleHighlightRegionRender();
    return;
  }
  
  auto interspec = InterSpec::instance();
  if( !m_spec )
  {
    const char *msg = "Time chart passed sample numbers when no data is set - should not happen";
    if( interspec )
      interspec->logMessage( msg, "", 3 );
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, msg );
#endif
    return;
  }//if( !m_spec )

#if( PERFORM_DEVELOPER_CHECKS )
  for( const int sample : sample_numbers )
  {
    if( !m_spec->sample_numbers().count(sample) )
    {
      log_developer_error( __func__, ("Was passed an invalid samplenumber: "
                                       + std::to_string(sample)).c_str() );
    }//if( spec file doesnt have the input sample number )
  }//for( loop over all incomming sample numbers )
#endif
  
  shared_ptr<const ColorTheme> theme = (interspec ? interspec->getColorTheme() : nullptr);
  
  D3TimeChart::HighlightRegion region;
  region.type = type;
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      region.color = theme ? theme->timeHistoryForegroundHighlight : Wt::WColor(255,255,0,155);
      break;
    case SpecUtils::SpectrumType::SecondForeground:
      region.color = theme ? theme->timeHistoryBackgroundHighlight : Wt::WColor(0,255,255,75);
      break;
    case SpecUtils::SpectrumType::Background:
      region.color = theme ? theme->timeHistorySecondaryHighlight : Wt::WColor(0,128,0,75);
      break;
  }//switch( type )
  
  int firstInRange = *(sample_numbers.begin());
  int previous = firstInRange;
  
  for( auto iter = begin(sample_numbers); iter != end(sample_numbers); ++iter )
  {
    const int thisval = *iter;
      
    if( (thisval > (previous+1)) )
    {
      region.start_sample_number = firstInRange;
      region.end_sample_number = previous;
      m_highlights.push_back( region );
      
      firstInRange = thisval;
    }//if( thisval > (previous+1) )
      
    previous = thisval;
  }//for( loop over smaple_numbers )
    
  region.start_sample_number = firstInRange;
  region.end_sample_number = previous;
  m_highlights.push_back( region );
  
  scheduleHighlightRegionRender();
}//setHighlightedIntervals(...)


Wt::Signal<int,Wt::WFlags<Wt::KeyboardModifier>> &D3TimeChart::chartClicked()
{
  return m_chartClicked;
}


Wt::Signal<int,int,Wt::WFlags<Wt::KeyboardModifier>> &D3TimeChart::chartDragged()
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


vector<pair<int,int>> D3TimeChart::sampleNumberRangesWithOccupancyStatus(
                                                const SpecUtils::OccupancyStatus wanted_occ_status,
                                                std::shared_ptr<const SpecUtils::SpecFile> spec )
{
  vector< pair<int,int> > occ_ranges;
  
  if( !spec )
    return occ_ranges;
  
  std::set<int> wantedSamples;
  const set<int> &allsamples = spec->sample_numbers();
  const vector<string> &detnames = spec->detector_names();
  
  bool inWantedRegion = false;
  int startSample = numeric_limits<int>::min(), prevSample = numeric_limits<int>::min();
    
  for( const int sample : allsamples )
  {
    bool wantSample = false;
    for( size_t i = 0; !wantSample && i < detnames.size(); ++i )
    {
      auto m = spec->measurement( sample, detnames[i] );
      wantSample = (m && (m->occupied() == wanted_occ_status));
    }
    
    if( wantSample )
      wantedSamples.insert( sample );
    
    if( wantSample && !inWantedRegion )
      startSample = sample;
    else if( inWantedRegion && !wantSample )
      occ_ranges.push_back( {startSample,prevSample} );
    
    prevSample = sample;
    inWantedRegion = wantSample;
  }//for( const int sample : allsamples )
  
  if( inWantedRegion )
    occ_ranges.push_back( {startSample,prevSample} );

  if( wantedSamples == spec->sample_numbers() )
    return {};
  
  return occ_ranges;
}//sampleNumberRangesWithOccupancyStatus(...)


void D3TimeChart::setHighlightRegionsToClient()
{
  auto type_to_str = []( const SpecUtils::SpectrumType type ) -> std::string {
    switch (type )
    {
      case SpecUtils::SpectrumType::Foreground:       return "FOREGROUND";
      case SpecUtils::SpectrumType::SecondForeground: return "SECONDARY";
      case SpecUtils::SpectrumType::Background:       return "BACKGROUND";
    }
    return "";
  };//type_to_str
  
  
  WStringStream js;
  js << m_jsgraph <<  ".setHighlightRegions( [";
  
  int nadded = 0;
  for( const auto &region : m_highlights )
  {
    js << string(nadded ? "," : "")
       << "{startSample: " << region.start_sample_number
       << ", endSample: " << region.end_sample_number
       << ", fillColor: \"" << region.color.cssText() << "\""
       << ", type: \"" << type_to_str(region.type) << "\"}";
  }//for( loop over highlight regions )
  js << "] );";
  
  cout << "\n\nWill set highlight regions with JSON=" + js.str() + "\n\n" << endl;
  doJavaScript( js.str() );
}//setHighlightRegionsToClient()


void D3TimeChart::saveChartToPng( const std::string &filename )
{
  // \TODO: implement - see D3SpectrumDisplayDiv for a starting point, although it isnt finished
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
  m_occLineColor = theme->occupancyIndicatorLines;
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
    setDataToClient();
  
  if( m_renderFlags.testFlag(TimeRenderActions::UpdateHighlightRegions) )
    setHighlightRegionsToClient();
  
  m_renderFlags = 0;
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void D3TimeChart::chartClickedCallback( int sample_number, int modifier_keys )
{
  m_chartClicked.emit(sample_number, Wt::WFlags<Wt::KeyboardModifier>(Wt::KeyboardModifier(modifier_keys)) );
}//chartClickedCallback(...)


void D3TimeChart::chartDraggedCallback( int first_sample_number, int last_sample_number, int modifier_keys )
{
  m_chartDragged.emit( first_sample_number, last_sample_number, Wt::WFlags<Wt::KeyboardModifier>(Wt::KeyboardModifier(modifier_keys)) );
}//chartDraggedCallback(...)


void D3TimeChart::chartResizedCallback( double chart_width_px, double chart_height_px )
{
  m_chartResized.emit( chart_width_px, chart_height_px );
}//chartResizedCallback(...)


void D3TimeChart::displayedXRangeChangeCallback( int first_sample_number, int last_sample_number, int samples_per_channel )
{
  m_displayedXRangeChange.emit( first_sample_number, last_sample_number, samples_per_channel );
}//displayedXRangeChangeCallback(...)
