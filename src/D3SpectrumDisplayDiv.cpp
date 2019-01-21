#include "InterSpec_config.h"

#include <memory>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/math/constants/constants.hpp>

#include <boost/version.hpp>
#if(BOOST_VERSION >= 104800)
#include <boost/timer/timer.hpp>
#endif

#include <Wt/WPen>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WBreak>
#include <Wt/WServer>
#include <Wt/WLength>
#include <Wt/WSpinBox>
#include <Wt/WIconPair>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WMessageBox>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WStringStream>
#include <Wt/WBorderLayout>
#include <Wt/WDoubleSpinBox>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WDataSeries>
#include <Wt/Chart/WAbstractChart>


#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "SpecUtils/SpectrumDataStructs.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "SpecUtils/D3SpectrumExport.h"


using namespace Wt;
using namespace std;

#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

D3SpectrumDisplayDiv::D3SpectrumDisplayDiv( WContainerWidget *parent )
: WContainerWidget( parent ),
m_model( new SpectrumDataModel( this ) ),
m_peakModel( 0 ),
m_layoutWidth( 0 ),
m_layoutHeight( 0 ),
m_autoAdjustDisplayBinnning( false ),
m_compactAxis( false ),
m_legendEnabled( true ),
m_yAxisIsLog( true ),
m_xAxisUnits( SpectrumChart::kUndefinedUnits ),
m_showRefLineInfoForMouseOver( true ),
m_foregroundLineColor( 0x00, 0x00, 0x00 ),  //black
m_backgroundLineColor( 0x00, 0xff, 0xff ),  //cyan
m_secondaryLineColor( 0x00, 0x80, 0x80 ),   //dark green
m_textColor( 0x00, 0x00, 0x00 ),
m_axisColor( 0x00, 0x00, 0x00 ),
m_chartMarginColor(),
m_chartBackgroundColor(),
m_defaultPeakColor( 0, 51, 255, 155 )
{
  
  
  setLayoutSizeAware( true );
  addStyleClass( "SpectrumDisplayDiv" );
  
  // Cancel right-click events for the div, we handle it all in JS
  setAttributeValue( "oncontextmenu",
                     "event.cancelBubble = true; event.returnValue = false; return false;"
                    );
  
  wApp->useStyleSheet("external_libs/SpecUtils/d3_resources/SpectrumChartD3.css");
  wApp->require( "external_libs/SpecUtils/d3_resources/d3.v3.min.js" );
  wApp->require( "external_libs/SpecUtils/d3_resources/c.min.js" );
  
  // Turn peak fitting off temporarily (Christian: 05282018)
  // wApp->require( "external_libs/SpecUtils/d3_resources/numeric-1.2.6.min.js" );
  // wApp->require( "external_libs/SpecUtils/d3_resources/PeakFit.js" );
  
  wApp->require( "external_libs/SpecUtils/d3_resources/SpectrumChartD3.js" );
  
  for( SpectrumChart::PeakLabels label = SpectrumChart::PeakLabels(0);
      label < SpectrumChart::PeakLabels::kNumPeakLabels; label = SpectrumChart::PeakLabels(label+1) )
  {
    m_peakLabelsToShow[label] = false;
  }//for( loop over all labels )
  
  const char *js = INLINE_JAVASCRIPT(
                                     function(self, w, h, layout) {
                                       window.graph.handleResize();
                                     });
  
  setJavaScriptMember( "wtResize", js );
    
  m_xRangeChangedJS.reset( new JSignal<double,double>( this, "xrangechanged", true ) );
  m_xRangeChangedJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartXRangeChangedCallback, this, _1, _2 ) );
  
  m_controlKeyDraggJS.reset( new JSignal<double,double,int,int>( this, "controlkeydragged", true ) );
  m_controlKeyDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartControlKeyDragCallback, this, _1, _2, _3, _4 ) );
  
  m_shiftKeyDraggJS.reset( new JSignal<double,double>( this, "shiftkeydragged", true ) );
  m_shiftKeyDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartShiftKeyDragCallback, this, _1, _2 ) );
  
  m_shiftAltKeyDraggJS.reset( new JSignal<double,double>( this, "shiftaltkeydragged", true ) );
  m_shiftAltKeyDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartShiftAltKeyDragCallback, this, _1, _2 ) );
  
  m_rightMouseDraggJS.reset( new JSignal<double,double>( this, "rightmousedragged", true ) );
  m_rightMouseDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartRightMouseDragCallback, this, _1, _2 ) );
  
  m_leftClickJS.reset( new JSignal<double,double,int,int>( this, "leftclicked", true ) );
  m_leftClickJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartLeftClickCallback, this, _1, _2, _3, _4 ) );
  
  m_doubleLeftClickJS.reset( new JSignal<double,double>( this, "doubleclicked", true ) );
  m_doubleLeftClickJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartDoubleLeftClickCallback, this, _1, _2 ) );
  
  m_rightClickJS.reset( new JSignal<double,double,int,int>( this, "rightclicked", true ) );
  m_rightClickJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartRightClickCallback, this, _1, _2, _3, _4 ) );
  
  //need legend closed signal.
  
  // Create the Spectrum Chart D3 object to display the chart
  doJavaScript(
               "window.graph = new SpectrumChartD3(" + jsRef() + ",{title:'',xlabel:'',ylabel:''});"
               + "window.graph.setShowUserLabels(false);"
               + "window.graph.setShowPeakLabels(false);"
               + "window.graph.setShowNuclideNames(false);"
               + "window.graph.setShowNuclideEnergies(false);"
               );
}//D3SpectrumDisplayDiv constructor

void D3SpectrumDisplayDiv::setTextInMiddleOfChart( const Wt::WString &s )
{
  // TODO: D3 chart currently does not have functionality to set text in the middle of chart...
}

void D3SpectrumDisplayDiv::setCompactAxis( const bool compact )
{
  m_compactAxis = compact;
  const string isCompact = compact ? "true" : "false";
  doJavaScript(
               "window.graph.setCompactXAxis(" + isCompact + ");"
               );
}

bool D3SpectrumDisplayDiv::isAxisCompacted() const
{
  return m_compactAxis;
}

void D3SpectrumDisplayDiv::setPeakModel( PeakModel *model )
{
  if( !model )
    throw runtime_error( "setPeakModel(...): invalid input model" );
  
  model->setDataModel( m_model );
  m_peakModel = model;
}//void setPeakModel( PeakModel *model );


void D3SpectrumDisplayDiv::prefferPngRenderForChart()
{
}


bool D3SpectrumDisplayDiv::legendIsEnabled() const
{
  return m_legendEnabled;
}


void D3SpectrumDisplayDiv::enableLegend( const bool forceMobileStyle )
{
  m_legendEnabled = true;
  m_legendEnabledSignal.emit();
  doJavaScript(
               "window.graph.setShowLegend(true);"
               );
  
}//void D3SpectrumDisplayDiv::enableLegend()


void D3SpectrumDisplayDiv::disableLegend()
{
  m_legendEnabled = false;
  m_legendDisabledSignal.emit();
  doJavaScript(
               "window.graph.setShowLegend(false);"
               );
}//void disableLegend()


void D3SpectrumDisplayDiv::setIsEnergyDisplay()
{
  m_xAxisUnits = SpectrumChart::kkeV;
}

void D3SpectrumDisplayDiv::setIsTimeDisplay()
{
  m_xAxisUnits = SpectrumChart::kSeconds;
}

bool D3SpectrumDisplayDiv::isEnergyDisplay() const
{
  return (m_xAxisUnits == SpectrumChart::kkeV);
}

bool D3SpectrumDisplayDiv::isTimeDisplay() const
{
  return (m_xAxisUnits == SpectrumChart::kSeconds);
}


void D3SpectrumDisplayDiv::setControlDragDebouncePeriod( int milliseconds )
{
}


void D3SpectrumDisplayDiv::setControlDragContinuumPreview( double x0,
                                                        double x1 )
{
}//void handleControlMouseMove( int x_start )



void D3SpectrumDisplayDiv::showHistogramIntegralsInLegend( const bool show )
{
  m_showHistogramIntegralsInLegend = show;
  // TODO: No option in D3 to show/hide histogram integrals in legend
}

void D3SpectrumDisplayDiv::enableOverlayCanvas( bool outline,
                                             bool highlight,
                                             bool enalbeAltShiftHighlight  )
{
}//void enableOverlayCanvas()


CanvasForDragging *D3SpectrumDisplayDiv::overlayCanvas()
{
  return NULL;
}//CanvasForDragging *overlayCanvas()

bool D3SpectrumDisplayDiv::overlayCanvasEnabled() const
{
  return false;
}//bool overlayCanvasEnabled() const

void D3SpectrumDisplayDiv::setOverlayCanvasVisible( bool visible )
{
}//void setOverlayCanvasVisible( bool visible )

void D3SpectrumDisplayDiv::disableOverlayCanvas()
{
}//void disableOverlayCanvas()

void D3SpectrumDisplayDiv::setScrollingParent( Wt::WContainerWidget *parent )
{
}//void setScrollingParent( Wt::WContainerWidget *parent );

void D3SpectrumDisplayDiv::setScrollY( int scrollY )
{
}//void setScrollY( int scrollY )

Signal<> &D3SpectrumDisplayDiv::legendEnabled()
{
  return m_legendEnabledSignal;
}


Signal<> &D3SpectrumDisplayDiv::legendDisabled()
{
  return m_legendDisabledSignal;
}


void D3SpectrumDisplayDiv::setHidden( bool hidden, const Wt::WAnimation &anim )
{
  WContainerWidget::setHidden( hidden, anim );
}//void setHidden( bool hidden, const Wt::WAnimation &animation )


void D3SpectrumDisplayDiv::setShowPeakLabel( int label, bool show )
{
  SpectrumChart::PeakLabels peakLabel = SpectrumChart::PeakLabels(label);
  m_peakLabelsToShow[ peakLabel ] = show;
  
  const string shouldShow = show ? "true" : "false";
  
  switch ( peakLabel )
  {
    case SpectrumChart::PeakLabels::kShowPeakUserLabel:
      doJavaScript("window.graph.setShowUserLabels(" + shouldShow + ");");
      break;
    case SpectrumChart::PeakLabels::kShowPeakEnergyLabel:
      doJavaScript("window.graph.setShowPeakLabels(" + shouldShow + ");");
      break;
    case SpectrumChart::PeakLabels::kShowPeakNuclideLabel:
      doJavaScript("window.graph.setShowNuclideNames(" + shouldShow + ");");
      break;
    case SpectrumChart::PeakLabels::kShowPeakNuclideEnergies:
      doJavaScript("window.graph.setShowNuclideEnergies(" + shouldShow + ");");
      break;
    case SpectrumChart::PeakLabels::kNumPeakLabels:
    default: break;
  }
}//void setShowPeakLabel( int label, bool show )


bool D3SpectrumDisplayDiv::showingPeakLabel( int peakLabel ) const
{
  return m_peakLabelsToShow.at( SpectrumChart::PeakLabels(peakLabel) );
}//bool showingPeakLabel( int peakLabel ) const


Signal<double,double,int,int> &D3SpectrumDisplayDiv::chartClicked()
{
  return m_leftClick;
}//Signal<double,double,int,int> &chartClicked()


Wt::Signal<double,double,int,int> &D3SpectrumDisplayDiv::rightClicked()
{
  return m_rightClick;
}


Wt::Signal<double,double> &D3SpectrumDisplayDiv::doubleLeftClick()
{
  return m_doubleLeftClick;
}

Wt::Signal<double,double,int,int> &D3SpectrumDisplayDiv::controlKeyDragged()
{
  return m_controlKeyDragg;
}

Wt::Signal<double,double> &D3SpectrumDisplayDiv::shiftKeyDragged()
{
  return m_shiftKeyDragg;
}


Wt::Signal<double,double> &D3SpectrumDisplayDiv::rightMouseDragg()
{
  return m_rightMouseDragg;
}//Signal<double,double> &rightMouseDragg()


Wt::JSlot *D3SpectrumDisplayDiv::alignOverlayCanvas()
{
  return NULL;
}


Wt::JSignal<std::string> *D3SpectrumDisplayDiv::overlayCanvasJsException()
{
  return NULL;
}


void D3SpectrumDisplayDiv::connectWtMouseConnections()
{
}


void D3SpectrumDisplayDiv::disconnectWtMouseConnections()
{
}


#if( IOS || ANDROID )
void D3SpectrumDisplayDiv::forceOverlayAlign()
{
  Wt::JSlot *align = alignOverlayCanvas();
  if( align )
    align->exec();
}//void forceOverlayAlign()
#endif


#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
void D3SpectrumDisplayDiv::setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )
{
  m_referencePhotoPeakLines = nuc;
  updateReferncePhotoPeakLines();
}

void D3SpectrumDisplayDiv::persistCurrentReferncePhotoPeakLines()
{
  if( m_referencePhotoPeakLines.energies.empty() )
    return;
  
  vector<const ReferenceLineInfo>::iterator pos;
  pos = std::find( m_persistedPhotoPeakLines.begin(),
                  m_persistedPhotoPeakLines.end(),
                  m_referencePhotoPeakLines );
  
  if( pos == m_persistedPhotoPeakLines.end()
     && m_referencePhotoPeakLines.displayLines )
    m_persistedPhotoPeakLines.push_back( m_referencePhotoPeakLines );
  
  m_referencePhotoPeakLines.reset();
  
  updateReferncePhotoPeakLines();
}

void D3SpectrumDisplayDiv::clearAllReferncePhotoPeakLines()
{
  m_referencePhotoPeakLines.reset();
  m_persistedPhotoPeakLines.clear();
  doJavaScript("window.graph.clearReferenceLines();");
}

void D3SpectrumDisplayDiv::updateReferncePhotoPeakLines()
{
  string result = "[";
  const ReferenceLineInfo &showingNuclide = m_referencePhotoPeakLines;
  bool addComma = !showingNuclide.energies.empty();
  
  if (addComma)
    showingNuclide.toJson(result);
  
  for (const ReferenceLineInfo ref : m_persistedPhotoPeakLines) {
    if (showingNuclide.energies.empty() || ref.parentLabel() != showingNuclide.parentLabel()) {
      if ( addComma )
        result += ",";
      addComma = true;
      ref.toJson(result);
    }
  }
  result += "]";
  doJavaScript("window.graph.setReferenceLines(" + result + ")");
}
#endif //#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )

void D3SpectrumDisplayDiv::setShowRefLineInfoForMouseOver( const bool show )
{
  m_showRefLineInfoForMouseOver = show;
  const string showRefLineInfo = show ? "true" : "false";
  doJavaScript(
               "window.graph.setShowRefLineInfoForMouseOver(" + showRefLineInfo + ")"
               );
}//void setShowRefLineInfoForMouseOver( const bool show )


void D3SpectrumDisplayDiv::layoutSizeChanged ( int width, int height )
{
  m_layoutWidth = width;
  m_layoutHeight = height;
  
  if( m_autoAdjustDisplayBinnning )
  {
    guessAndUpdateDisplayRebinFactor();
  }//if( m_autoAdjustDisplayBinnning )
  
  
#if( IOS || ANDROID )
  //When the soft-keyboard disapears (on Android at a minimum), the overlays
  //  dont resize properly (until you change tab bellow the chart, or something)
  //  so we will force it.
  forceOverlayAlign();
#endif
}//void layoutSizeChanged ( int width, int height )


int D3SpectrumDisplayDiv::layoutWidth() const
{
  return m_layoutWidth;
}//int layoutWidth() const


int D3SpectrumDisplayDiv::layoutHeight() const
{
  return m_layoutHeight;
}//int layoutHeight() const


double D3SpectrumDisplayDiv::xAxisMinimum() const
{
  return m_xAxisMinimum;
}//double xAxisMinimum() const


double D3SpectrumDisplayDiv::xAxisMaximum() const
{
  return m_xAxisMaximum;
}//double xAxisMaximum() const


double D3SpectrumDisplayDiv::yAxisMinimum() const
{
  return m_yAxisMinimum;
}//double yAxisMinimum() const


double D3SpectrumDisplayDiv::yAxisMaximum() const
{
  return m_yAxisMaximum;
}//double yAxisMaximum() const


bool D3SpectrumDisplayDiv::yAxisIsLog() const
{
  return m_yAxisIsLog;
}//bool yAxisIsLog() const;


void D3SpectrumDisplayDiv::setYAxisLog( bool log )
{
  m_yAxisIsLog = log;
  doJavaScript(
               log ? "window.graph.setLogY();" : "window.graph.setLinearY();"
               );
}//void setYAxisLog( bool log )

void D3SpectrumDisplayDiv::showGridLines( bool show )
{
  const string shouldDraw = show ? "true" : "false";
  m_showVerticalLines = show;
  m_showHorizontalLines = show;
  doJavaScript(
               "window.graph.setGridX(" + shouldDraw + ");" +
               "window.graph.setGridY(" + shouldDraw + ");"
               );
}

void D3SpectrumDisplayDiv::showVerticalLines( const bool draw )
{
  const string shouldDraw = draw ? "true" : "false";
  m_showVerticalLines = draw;
  doJavaScript(
               "window.graph.setGridX(" + shouldDraw + ");"
               );
}

void D3SpectrumDisplayDiv::showHorizontalLines( const bool draw )
{
  const string shouldDraw = draw ? "true" : "false";
  m_showHorizontalLines = draw;
  doJavaScript(
               "window.graph.setGridY(" + shouldDraw + ");"
               );
}

bool D3SpectrumDisplayDiv::verticalLinesShowing() const
{
  return m_showVerticalLines;
}

bool D3SpectrumDisplayDiv::horizontalLinesShowing() const
{
  return m_showHorizontalLines;
}


bool D3SpectrumDisplayDiv::backgroundSubtract() const
{
  return m_model->backgroundSubtract();
}//bool backgroundSubtract() const


void D3SpectrumDisplayDiv::setBackgroundSubtract( bool subtract )
{
  if( subtract == m_model->backgroundSubtract() )
    return;
  
  m_model->setBackgroundSubtract( subtract );
  
  const string shouldSubtract = subtract ? "true" : "false";
  doJavaScript(
               "window.graph.setBackgroundSubtract(" + shouldSubtract + ");"
               );
}//void setBackgroundSubtract( bool subtract )

void D3SpectrumDisplayDiv::setXAxisMinimum( const double minimum )
{
  const string minimumStr = to_string( minimum );
  m_xAxisMinimum = minimum;
  doJavaScript(
               "window.graph.setXAxisMinimum(" + minimumStr + ");"
               );
}//void setXAxisMinimum( const double minimum )


void D3SpectrumDisplayDiv::setXAxisMaximum( const double maximum )
{
  const string maximumStr = to_string( maximum );
  m_xAxisMaximum = maximum;
  doJavaScript(
               "window.graph.setXAxisMaximum(" + maximumStr + ");"
               );
}//void setXAxisMaximum( const double maximum )


void D3SpectrumDisplayDiv::setYAxisMinimum( const double minimum )
{
  const string minimumStr = to_string( minimum );
  m_yAxisMinimum = minimum;
  doJavaScript(
               "window.graph.setYAxisMinimum(" + minimumStr + ");"
               );
}//void setYAxisMinimum( const double minimum )


void D3SpectrumDisplayDiv::setYAxisMaximum( const double maximum )
{
  const string maximumStr = to_string( maximum );
  m_yAxisMaximum = maximum;
  doJavaScript(
               "window.graph.setYAxisMaximum(" + maximumStr + ");"
               );
}//void setYAxisMaximum( const double maximum )


void D3SpectrumDisplayDiv::setXAxisRange( const double minimum, const double maximum )
{
  const string minimumStr = to_string( minimum );
  const string maximumStr = to_string( maximum );
  m_xAxisMinimum = minimum;
  m_xAxisMaximum = maximum;
  doJavaScript(
               "window.graph.setXAxisRange(" + minimumStr + "," + maximumStr + ");"
               );
  if( m_autoAdjustDisplayBinnning )
    m_xRangeChanged.emit( minimum, maximum );
}//void setXAxisRange( const double minimum, const double maximum );


void D3SpectrumDisplayDiv::setYAxisRange( const double minimum,
                                       const double maximum )
{
  const string minimumStr = to_string( minimum );
  const string maximumStr = to_string( maximum );
  m_yAxisMinimum = minimum;
  m_yAxisMaximum = maximum;
  doJavaScript(
               "window.graph.setYAxisRange(" + minimumStr + "," + maximumStr + ");"
               );
}//void setYAxisRange( const double minimum, const double maximum );


void D3SpectrumDisplayDiv::setDisplayRebinFactor( const int factor )
{
  m_model->setRebinFactor( factor );
}//void setDisplayRebinFactor( const int factor )


void D3SpectrumDisplayDiv::setAutoAdjustDisplayRebinFactor( bool auto_rebin )
{
  m_autoAdjustDisplayBinnning = auto_rebin;
}//void setAutoAdjustDisplayRebinFactor( bool auto_rebin )


int D3SpectrumDisplayDiv::displayRebinFactor() const
{
  return m_model->rebinFactor();
}//int displayRebinFactor() const


void D3SpectrumDisplayDiv::guessAndUpdateDisplayRebinFactor()
{
  //my guess is this function houldnt have ant real effect (in terms of
  //  server/client interactions) if it is being called due to itself (e.g.
  //  calling setDisplayRebinFactor(...) with a NEW rebin factor causes
  //  the model to change dimensions, which causes the chart to resize).
  std::shared_ptr<Measurement> axisH = m_model->histUsedForXAxis();
  
  if( !axisH || !layoutSizeAware() )
  {
    setDisplayRebinFactor( 1 );
    return;
  }//if( !axisH )
  
  const float displayedxmin = static_cast<float>( m_xAxisMinimum );
  const float displayedxmax = static_cast<float>( m_xAxisMaximum );
  const size_t displayednbin = axisH->find_gamma_channel( displayedxmax )
  - axisH->find_gamma_channel( displayedxmin );
  const int width = layoutWidth()
  - m_plotAreaPaddingLeft
  - m_plotAreaPaddingRight;
  const float bins_per_pixel = float(displayednbin) / float(width);
  const int rebin_factor = max( static_cast<int>(ceil(bins_per_pixel)), 1 );
  
  setDisplayRebinFactor( rebin_factor );
}//void guessAndUpdateDisplayRebinFactor()


void D3SpectrumDisplayDiv::setData( std::shared_ptr<Measurement> data_hist,
                                 float liveTime,
                                 float realTime,
                                 float neutronCounts,
                                 bool keep_curent_xrange )
{
  m_model->setDataHistogram( data_hist, liveTime, realTime, neutronCounts );
  
  const string resetDomain = keep_curent_xrange ? "false" : "true";
  
  // Set the data for the chart
  if ( data_hist ) {
    // Create the measurement array (should only have one measurement)
    std::ostringstream ostr;
    std::vector< std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> foregroundData;
    D3SpectrumExport::D3SpectrumOptions foregroundOptions;
    
    // Set options for the spectrum
    foregroundOptions.line_color = m_foregroundLineColor.isDefault() ? string("black") : m_foregroundLineColor.cssText();
    foregroundOptions.peak_color = m_defaultPeakColor.isDefault() ? string("blue") : m_defaultPeakColor.cssText();
    foregroundOptions.spectrum_type = kForeground;
    foregroundOptions.display_scale_factor = displayScaleFactor( kForeground );
    
    // Set the peak data for the spectrum
    if ( m_peakModel ) {
      std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
      vector< std::shared_ptr<const PeakDef> > inpeaks( peaks->begin(), peaks->end() );
      foregroundOptions.peaks_json = PeakDef::peak_json( inpeaks );
    }
    
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(data_hist.get(),foregroundOptions) );
    
    // Set the data on the JS side
    if ( D3SpectrumExport::write_and_set_data_for_chart(ostr, id(), measurements) ) {
      string data = ostr.str();
      size_t index = data.find( "spec_chart_" );
      data = data.substr( 0, index );
      doJavaScript( data + "window.graph.setSpectrumData(data_" + id() + ", " + resetDomain + ", 'FOREGROUND', 0, 1 );" );
    }
  } else {
    doJavaScript( "window.graph.removeSpectrumData(" + resetDomain + ", 'FOREGROUND' );" );
  }//if ( data_hist )
}//void setData( std::shared_ptr<Measurement> data_hist )


std::shared_ptr<Measurement> D3SpectrumDisplayDiv::data()
{
  return m_model->getData();
}//std::shared_ptr<Measurement> data()


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::data() const
{
  return m_model->getData();
}//std::shared_ptr<const Measurement> data() const


std::shared_ptr<Measurement> D3SpectrumDisplayDiv::secondData()
{
  return m_model->getSecondData();
}//std::shared_ptr<Measurement> secondData()


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::secondData() const
{
  return m_model->getSecondData();
}//std::shared_ptr<const Measurement> secondData() const


std::shared_ptr<Measurement> D3SpectrumDisplayDiv::background()
{
  return m_model->getBackground();
}//std::shared_ptr<Measurement> background()


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::background() const
{
  return m_model->getBackground();
}//std::shared_ptr<const Measurement> background() const


float D3SpectrumDisplayDiv::foregroundLiveTime() const
{
  return m_model->dataLiveTime();
}

float D3SpectrumDisplayDiv::foregroundRealTime() const
{
  return m_model->dataRealTime();
}

float D3SpectrumDisplayDiv::backgroundLiveTime() const
{
  return m_model->backgroundLiveTime();
}

float D3SpectrumDisplayDiv::backgroundRealTime() const
{
  return m_model->backgroundRealTime();
}

float D3SpectrumDisplayDiv::secondForegroundLiveTime() const
{
  return m_model->secondDataLiveTime();
}

float D3SpectrumDisplayDiv::secondForegroundRealTime() const
{
  return m_model->secondDataRealTime();
}


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::histUsedForXAxis() const
{
  return m_model->histUsedForXAxis();
}


void D3SpectrumDisplayDiv::setDisplayScaleFactor( const float sf,
                                               const SpectrumType spectrum_type )
{
  switch( spectrum_type )
  {
    case kForeground:
      throw runtime_error( "setDisplayScaleFactor can not be called for foreground" );
      
    case kSecondForeground:
      m_model->setSecondDataScaleFactor( sf );
      updateSecondData();
      break;
      
    case kBackground:
      m_model->setBackgroundDataScaleFactor( sf );
      updateBackground();
      break;
  }//switch( spectrum_type )
  
}//void setDisplayScaleFactor(...)


float D3SpectrumDisplayDiv::displayScaleFactor( const SpectrumType spectrum_type ) const
{
  switch( spectrum_type )
  {
    case kForeground:
      return 1.0f;
    case kSecondForeground:
      return m_model->secondDataScaledBy();
    case kBackground:
      return m_model->backgroundScaledBy();
      //  m_spectrumDiv->continuum();
  }//switch( spectrum_type )
  
  throw runtime_error( "D3SpectrumDisplayDiv::displayScaleFactor(...): invalid input arg" );
  
  return 1.0;
}//double displayScaleFactor( SpectrumType spectrum_type ) const;


void D3SpectrumDisplayDiv::setBackground( std::shared_ptr<Measurement> background,
                                       float liveTime,
                                       float realTime,
                                       float neutronCounts )
{
  m_model->setBackgroundHistogram( background, liveTime, realTime, neutronCounts );
  
  if( !background && m_model->backgroundSubtract() )
    m_model->setBackgroundSubtract( false );
  
  updateBackground();
}//void D3SpectrumDisplayDiv::setBackground(...);


void D3SpectrumDisplayDiv::setSecondData( std::shared_ptr<Measurement> hist,
                                       float liveTime,
                                       float realTime,
                                       float neutronCounts,
                                       bool ownAxis )
{
  m_model->setSecondDataHistogram( hist, liveTime, realTime, neutronCounts, ownAxis );
  
  updateSecondData();
}//void D3SpectrumDisplayDiv::setSecondData( std::shared_ptr<Measurement> background );




void D3SpectrumDisplayDiv::clearTimeHighlightRegions( const SpectrumType type )
{
}//void clearTimeHighlightRegions();


void D3SpectrumDisplayDiv::setTimeHighLightRegions( const vector< pair<double,double> > &p,
                                                 const SpectrumType type )
{
}//void setHighLightRegions(...)


bool D3SpectrumDisplayDiv::removeDecorativeHighlightRegion( size_t uniqueid )
{
  return true;
}//void removeDecorativeHighlightRegions()


size_t D3SpectrumDisplayDiv::addDecorativeHighlightRegion( const float lowerx,
                                                        const float upperx,
                                                        const Wt::WColor &color )
{
  return 0;
}//void addDecorativeHighlightRegion(...)


void D3SpectrumDisplayDiv::setAutoAxisRange()
{
}//void setAutoAxisRange()



void D3SpectrumDisplayDiv::setPlotAreaPadding( int left, int top, int right, int bottom )
{
  m_plotAreaPaddingLeft = left;
  m_plotAreaPaddingTop = top;
  m_plotAreaPaddingRight = right;
  m_plotAreaPaddingBottom = bottom;
}//void setPlotAreaPadding( int left, int top, int right, int bottom )


void D3SpectrumDisplayDiv::visibleRange( double &xmin, double &xmax,
                                      double &ymin, double &ymax ) const
{
  xmin = m_xAxisMinimum;
  xmax = m_xAxisMaximum;
  ymin = m_yAxisMinimum;
  ymax = m_yAxisMaximum;
}

int D3SpectrumDisplayDiv::plotAreaPadding( const Wt::Side side ) const
{
  switch ( side ) {
    case Left: return m_plotAreaPaddingLeft;
    case Right: return m_plotAreaPaddingRight;
    case Top: return m_plotAreaPaddingTop;
    case Bottom: return m_plotAreaPaddingBottom;
    default: return 0;
  }
}//int plotAreaPadding( const Wt::Side side ) const


const string D3SpectrumDisplayDiv::xAxisTitle() const
{
  return m_xAxisTitle;
}//const Wt::WString &xAxisTitle() const;


const string D3SpectrumDisplayDiv::yAxisTitle() const
{
  return m_yAxisTitle;
}//const Wt::WString &yAxisTitle() const;


const string D3SpectrumDisplayDiv::y2AxisTitle() const
{
  return m_y2AxisTitle;
}//const Wt::WString &y2AxisTitle() const;


void D3SpectrumDisplayDiv::setXAxisTitle( const std::string &title )
{
  m_xAxisTitle = title;
  doJavaScript(
               "window.graph.setXAxisTitle('" + title + "');"
               );
}//void setXAxisTitle( const std::string &title )


void D3SpectrumDisplayDiv::setYAxisTitle( const std::string &title )
{
  m_yAxisTitle = title;
  doJavaScript(
               "window.graph.setYAxisTitle('" + title + "');"
               );
}//void setYAxisTitle( const std::string &title )


void D3SpectrumDisplayDiv::setY2AxisTitle( const std::string &title )
{
  // Since no Y2 Axis title exists for current set of D3 charts, we ignore this (for now)
  m_y2AxisTitle = title;
}//void setY2AxisTitle( const std::string &title )


float D3SpectrumDisplayDiv::xUnitsPerPixel() const
{
  //Christian: We use 0.001 as a placeholder since we don't have reference to SpectrumChart anymore
  return 0.001;
}

void D3SpectrumDisplayDiv::setMouseDragZooms()
{
}//void D3SpectrumDisplayDiv::setMouseDragZooms()


void D3SpectrumDisplayDiv::setMouseDragHighlights( const bool allowMultiple,
                                                const bool allowSingleClick )
{
}//void setMouseDragHighlights( const bool allowMultiple )


void D3SpectrumDisplayDiv::allowArrowToMoveSingleClickRegion( bool allow )
{
}//void allowArrowToMoveSingleClickRegion( allow )


void D3SpectrumDisplayDiv::disableMouseDragActions()
{
}//void disableMouseDragActions()

//Christian: I'm not exactly if this signal is emitted with the D3 chart, since I saw
//  that it's main use was with time series. If D3 should support time series, then
//  I believe this should have use. But for now, this is just a placeholder.
Wt::Signal<double,double> &D3SpectrumDisplayDiv::xRangeChanged()
{
  return m_xRangeChanged;
}//xRangeChanged()


Wt::Signal<double,double> &D3SpectrumDisplayDiv::altKeyDragged()
{
  return m_altKeyDragg;
}

Wt::Signal<double,double> &D3SpectrumDisplayDiv::shiftAltKeyDragged()
{
  return m_shiftAltKeyDragg;
}

Wt::Signal<double,double> &D3SpectrumDisplayDiv::controlMouseMoved()
{
  return m_controlMouseMoved;
}

void D3SpectrumDisplayDiv::updateData()
{
  const std::shared_ptr<Measurement> data = m_model->getData();
  if ( data )
    setData( data,
             m_model->dataLiveTime(),
             m_model->dataRealTime(),
             m_model->dataNeutronCounts(),
             true );
}//void D3SpectrumDisplayDiv::updateData()

void D3SpectrumDisplayDiv::updateBackground()
{
  const std::shared_ptr<Measurement> background = m_model->getBackground();
  
  // Set the data for the chart
  if ( background ) {
    // Create the measurement array (should only have one measurement)
    std::ostringstream ostr;
    std::vector< std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> backgroundData;
    D3SpectrumExport::D3SpectrumOptions backgroundOptions;
    
    // Set options for the spectrum
    backgroundOptions.line_color = m_backgroundLineColor.isDefault() ? string("green") : m_backgroundLineColor.cssText();
    backgroundOptions.spectrum_type = kBackground;
    backgroundOptions.display_scale_factor = displayScaleFactor( kBackground );
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(background.get(), backgroundOptions) );
    
    // Set the data on the JS side
    if ( D3SpectrumExport::write_and_set_data_for_chart(ostr, id(), measurements) ) {
      string data = ostr.str();
      size_t index = data.find( "spec_chart_" );
      data = data.substr( 0, index );
      doJavaScript( data + "window.graph.setSpectrumData(data_" + id() + ", false, 'BACKGROUND', 1, -1);" );
    }
  } else {
    doJavaScript( "window.graph.removeSpectrumData(false, 'BACKGROUND' );" );
  }//if ( background )
}//void D3SpectrumDisplayDiv::updateBackground()

void D3SpectrumDisplayDiv::updateSecondData()
{
  const std::shared_ptr<Measurement> hist = m_model->getSecondData();
  
  // Set the data for the chart
  if ( hist ) {
    // Create the measurement array (should only have one measurement)
    std::ostringstream ostr;
    std::vector< std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> secondaryData;
    D3SpectrumExport::D3SpectrumOptions secondaryOptions;
    
    // Set options for the spectrum
    secondaryOptions.line_color = m_secondaryLineColor.isDefault() ? string("steelblue") : m_backgroundLineColor.cssText();
    secondaryOptions.spectrum_type = kSecondForeground;
    secondaryOptions.display_scale_factor = displayScaleFactor( kSecondForeground );
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(hist.get(), secondaryOptions) );
    
    // Set the data on the JS side
    if ( D3SpectrumExport::write_and_set_data_for_chart(ostr, id(), measurements) ) {
      string data = ostr.str();
      size_t index = data.find( "spec_chart_" );
      data = data.substr( 0, index );
      doJavaScript( data + "window.graph.setSpectrumData(data_" + id() + ", false, 'SECONDARY', 2, 1);" );
    }
  } else {
    doJavaScript( "window.graph.removeSpectrumData(false, 'SECONDARY' );" );
  }//if ( hist )
}//void D3SpectrumDisplayDiv::updateSecondData()


void D3SpectrumDisplayDiv::setForegroundSpectrumColor( const Wt::WColor &color )
{
  m_foregroundLineColor = color.isDefault() ? WColor( 0x00, 0x00, 0x00 ) : color;
  updateData();
}

void D3SpectrumDisplayDiv::setBackgroundSpectrumColor( const Wt::WColor &color )
{
  m_backgroundLineColor = color.isDefault() ? WColor(0x00,0xff,0xff) : color;
  updateBackground();
}

void D3SpectrumDisplayDiv::setSecondarySpectrumColor( const Wt::WColor &color )
{
  m_secondaryLineColor = color.isDefault() ? WColor(0x00,0x80,0x80) : color;
  updateSecondData();
}

void D3SpectrumDisplayDiv::setTextColor( const Wt::WColor &color )
{
  m_textColor = color.isDefault() ? WColor(0,0,0) : color;
  #warning "D3SpectrumDisplayDiv::setTextColor() Not Implemented"
}


void D3SpectrumDisplayDiv::setAxisLineColor( const Wt::WColor &color )
{
  const string rulename = "SpectrumChartAxisLineColor";
  m_axisColor = color.isDefault() ? WColor(0,0,0) : color;
  
  //".tick > line"
  //.attr( 'stroke', 'green' );
  
  //$('.xaxis').css( 'stroke', 'red' );
  //$('.yaxis').css( 'stroke', 'red' );
  
  
  //WCssStyleSheet &style = wApp->styleSheet();
  //WCssTextRule *rule = style.addRule(<#const std::string &selector#>, <#const Wt::WString &declarations#>)
  
  //style.ruleModified(...)
#warning "D3SpectrumDisplayDiv::setAxisLineColor() Not Implemented"
}

void D3SpectrumDisplayDiv::setChartMarginColor( const Wt::WColor &color )
{
  m_chartMarginColor = color;
  const string c = color.isDefault() ? string("rgba(0,0,0,0)") : color.cssText();
  doJavaScript( "$('#" + id() + " > svg').css('background','" + c + "');" );
}

void D3SpectrumDisplayDiv::setChartBackgroundColor( const Wt::WColor &color )
{
  m_chartBackgroundColor = color;
  const string c = color.isDefault() ? string("rgba(0,0,0,0)") : color.cssText();
  doJavaScript( "$('#" + id() + " > svg > g > rect').css('background','" + c + "');" );
}

void D3SpectrumDisplayDiv::setDefaultPeakColor( const Wt::WColor &color )
{
  m_defaultPeakColor = color.isDefault() ? WColor(0,51,255,155) : color;
  updateData();
}



void D3SpectrumDisplayDiv::removeAllPeaks()
{
  if ( m_peakModel ) {
    m_peakModel->removeAllPeaks();
    updateData();
  }
}

void D3SpectrumDisplayDiv::setFeatureMarkerOption( InterSpec::FeatureMarkerType option, bool show )
{
  const string shouldShow = show ? "true" : "false";
  switch ( option )
  {
    case InterSpec::FeatureMarkerType::EscapePeakMarker:
      doJavaScript(
                   "window.graph.setEscapePeaks(" + shouldShow + ");"
                   );
      break;
    case InterSpec::FeatureMarkerType::ComptonPeakMarker:
      doJavaScript(
                   "window.graph.setComptonPeaks(" + shouldShow + ");"
                   );
      break;
    case InterSpec::FeatureMarkerType::ComptonEdgeMarker:
      doJavaScript(
                   "window.graph.setComptonEdge(" + shouldShow + ");"
                   );
      break;
    case InterSpec::FeatureMarkerType::SumPeakMarker:
      doJavaScript(
                   "window.graph.setSumPeaks(" + shouldShow + ");"
                   );
      break;
    case InterSpec::FeatureMarkerType::NumFeatureMarkers:
    default: break;
  }
}//void D3SpectrumDisplayDiv::setFeatureMarkerOption(...)

void D3SpectrumDisplayDiv::chartControlKeyDragCallback( double x0, double x1, int pageX, int pageY )
{
  m_controlKeyDragg.emit( x0, x1, pageX, pageY );
}//void D3SpectrumDisplayDiv::chartControlKeyDragCallback(...)

void D3SpectrumDisplayDiv::chartShiftKeyDragCallback( double x0, double x1 )
{
  m_shiftKeyDragg.emit( x0, x1 );
}//void D3SpectrumDisplayDiv::chartShiftKeyDragCallback(...)

void D3SpectrumDisplayDiv::chartShiftAltKeyDragCallback( double x0, double x1 )
{
  m_shiftAltKeyDragg.emit( x0, x1 );
}//void D3SpectrumDisplayDiv::chartShiftAltKeyDragCallback(...)

void D3SpectrumDisplayDiv::chartRightMouseDragCallback( double x0, double x1 )
{
  m_rightMouseDragg.emit( x0, x1 );
}//void D3SpectrumDisplayDiv::chartRightMouseDragCallback(...)

void D3SpectrumDisplayDiv::chartLeftClickCallback( double x, double y, int pageX, int pageY )
{
  m_leftClick.emit( x, y, pageX, pageY );
}//void D3SpectrumDisplayDiv::chartDoubleLeftClickCallback(...)

void D3SpectrumDisplayDiv::chartDoubleLeftClickCallback( double x, double y )
{
  m_doubleLeftClick.emit( x, y );
}//void D3SpectrumDisplayDiv::chartDoubleLeftClickCallback(...)

void D3SpectrumDisplayDiv::chartRightClickCallback( double x, double y, int pageX, int pageY )
{
  m_rightClick.emit( x, y, pageX, pageY );
}//void D3SpectrumDisplayDiv::chartRightClickCallback(...)

void D3SpectrumDisplayDiv::chartXRangeChangedCallback( double x, double y )
{
  m_xRangeChanged.emit( x, y );
}//void D3SpectrumDisplayDiv::chartXRangeChangedCallback(...)

D3SpectrumDisplayDiv::~D3SpectrumDisplayDiv()
{
}//~D3SpectrumDisplayDiv()


