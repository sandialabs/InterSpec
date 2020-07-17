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

#include <boost/math/constants/constants.hpp>

#include <boost/version.hpp>

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
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/SpectrumDisplayDiv.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace Wt;
using namespace std;

using SpecUtils::Measurement;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

SpectrumDisplayDiv::SpectrumDisplayDiv( WContainerWidget *parent )
 : WContainerWidget( parent ),
#if(!SpectrumDisplayDiv_USE_CUSTOM_LAYOUT)
   m_chart( new SpectrumChart( NULL ) ),
#else
   m_chart( new SpectrumChart( this ) ),
#endif
   m_model( new SpectrumDataModel( this ) ),
   m_layoutWidth( 0 ),
   m_layoutHeight( 0 ),
   m_autoAdjustDisplayBinnning( false )
{
  setLayoutSizeAware( true );
  addStyleClass( "SpectrumDisplayDiv" );
  
  m_chart->setModel( m_model );
  m_chart->setXSeriesColumn( 0 );
  m_chart->xRangeChanged().connect( boost::bind( &SpectrumDisplayDiv::chartXRangeChangedCallback, this, _1, _2 ) );
  
#if(!SpectrumDisplayDiv_USE_CUSTOM_LAYOUT)
  WGridLayout *layout = new WGridLayout();
  setLayout( layout );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->addWidget( m_chart, 0, 0 );
  
//#if( ANDROIS || IOS )
//#else
  const char *js = INLINE_JAVASCRIPT(
    function(self, w, h, layout) {
      if( self.children.length === 1 && self.children[0].children.length===1 ){
        var c = self.children[0].children[0];
        try{ Wt.WT.AlignOverlay(c.id+'Cover',c.id); }catch(e){}
        self.style.overflow = 'hidden';
        try{ Wt.WT.AlignLegend(c.id); }catch(e){}
        self.style.overflow = 'hidden';
        var problemNode = document.getElementById('p'+c.id);
        if(problemNode){try{problemNode.style.overflowY='hidden';}catch(e){};}
      }else{
        console.log( "SpectrumDisplayDiv and SpectrumChart should have exactly one child (have " + self.children.length + ")" );
      }
  });
  
  setJavaScriptMember( "wtResize", js );
//#endif
  
#else
  
  const string js =
  """function(self, w, h,layout) {"
  ""  "if (!self.wtWidth || self.wtWidth!=w "
  ""      "|| !self.wtHeight || self.wtHeight!=h) {"
  ""    "self.wtWidth=w; self.wtHeight=h;"
  ""    "self.style.height=h + 'px';"
  ""    "Wt.emit(self,'resized',Math.round(w),Math.round(h));"
  ""  "}"
  ""  "var child = $('#" + m_chart->id() + "').get(0);"
  ""  "if(child && child.wtResize) {"
  ""     "child.wtResize(child,w,h,layout);"
  ""     "try{ Wt.WT.AlignOverlay(child.id + 'Cover', child.id ); }catch(e){}"
  ""     "self.style.overflow = 'hidden';"
  ""     "try{ Wt.WT.AlignLegend('" + m_chart->id() + "'); }catch(e){}"
  ""     "self.style.overflow = 'hidden';"
  ""     "var probNode = document.getElementById('p'+child.id);"
  ""     "if(probNode){try{probNode.style.overflowY='hidden';}catch(e){};}"
  ""   "}"
  "" "}";
  setJavaScriptMember( "wtResize", js );
#endif

//  setPlotAreaPadding( 80, 2, 10, 42 );

  const vector<Chart::WDataSeries> series = m_model->suggestDataSeries();
  m_chart->setSeries( series );
  m_chart->axis(Chart::Y2Axis).setVisible( m_model->secondDataOwnAxis() );
}//SpectrumDisplayDiv constructor

void SpectrumDisplayDiv::setTextInMiddleOfChart( const Wt::WString &s )
{
  m_chart->setTextInMiddleOfChart(s);
}

void SpectrumDisplayDiv::setCompactAxis( const bool compact )
{
  m_chart->setCompactAxis( compact );
}

bool SpectrumDisplayDiv::isAxisCompacted() const
{
  return m_chart->isAxisCompacted();
}

void SpectrumDisplayDiv::setAvoidMobileMenu( const bool avoid )
{
  m_chart->setAvoidMobileMenu( avoid );
}

void SpectrumDisplayDiv::setPeakModel( PeakModel *model )
{
  if( !model )
    throw runtime_error( "setPeakModel(...): invalid input model" );
  
//  model->setDataModel( dynamic_cast<SpectrumDataModel *>(m_chart->model()) );
  model->setDataModel( m_model );
  m_chart->setPeakModel( model );
}//void setPeakModel( PeakModel *model );


void SpectrumDisplayDiv::prefferPngRenderForChart()
{
  m_chart->setPreferredMethod( WPaintedWidget::PngImage );
}


bool SpectrumDisplayDiv::legendIsEnabled() const
{
  return m_chart->legendIsEnabled();
}


void SpectrumDisplayDiv::enableLegend( const bool forceMobileStyle )
{
  m_chart->enableLegend( forceMobileStyle );
}//void SpectrumDisplayDiv::enableLegend()


void SpectrumDisplayDiv::disableLegend()
{
  m_chart->disableLegend();
}//void disableLegend()


void SpectrumDisplayDiv::setIsEnergyDisplay()
{
  m_chart->setXAxisUnits( SpectrumChart::kkeV );
}

void SpectrumDisplayDiv::setIsTimeDisplay()
{
  m_chart->setXAxisUnits( SpectrumChart::kSeconds );
}


void SpectrumDisplayDiv::saveChartToPng( const std::string &name )
{
  m_chart->saveChartToPng( name );
}


bool SpectrumDisplayDiv::isEnergyDisplay() const
{
  return (m_chart->xAxisUnits() == SpectrumChart::kkeV);
}

bool SpectrumDisplayDiv::isTimeDisplay() const
{
  return (m_chart->xAxisUnits() == SpectrumChart::kSeconds);
}


void SpectrumDisplayDiv::setControlDragDebouncePeriod( int milliseconds )
{
  m_chart->setControlDragDebouncePeriod( milliseconds );
}


void SpectrumDisplayDiv::setControlDragContinuumPreview( double x0,
                                                         double x1 )
{
  SpectrumDataModel *theModel
                      = dynamic_cast<SpectrumDataModel *>( m_chart->model() );
  if( !theModel || !m_chart->overlayCanvas() )
    return;
  
  try
  {
    const int row0 = theModel->findRow( x0 );
    const int row1 = theModel->findRow( x1 );
    const double y0 = theModel->data( row0, 1 );
    const double y1 = theModel->data( row1, 1 );
    
    //Currently
    WStringStream js;
    js << "$('#c" << m_chart->overlayCanvas()->id() << "').data('ContEst',"
    << "{x0:"<< x0 << ", x1:"<< x1 << ", y0:" << y0 << ", y1:" << y1 << "});";
    
    CanvasForDragging *canvas = m_chart->overlayCanvas();
    if( canvas && canvas->controlDragDebounceTimePeriod() )
      js << "Wt.WT.DrawContEst('" << canvas->id() << "');";
    
    //Note: you must use WApplication::doJavaScript(...) here, or else the
    //      <canvas> will get erased - I'm not exactly sure as to why this is
    wApp->doJavaScript( js.str() );
  }catch(...)
  {
    cerr << "Caught exception in setControlDragContinuumPreview(...)" << endl;
  }
}//void handleControlMouseMove( int x_start )



void SpectrumDisplayDiv::showHistogramIntegralsInLegend( const bool show )
{
  m_chart->showHistogramIntegralsInLegend( show );
}

void SpectrumDisplayDiv::enableOverlayCanvas( bool outline,
                                              bool highlight,
                                              bool enalbeAltShiftHighlight  )
{
  m_chart->enableOverlayCanvas( outline, highlight, enalbeAltShiftHighlight );
}//void enableOverlayCanvas()


CanvasForDragging *SpectrumDisplayDiv::overlayCanvas()
{
  return m_chart->overlayCanvas();
}//CanvasForDragging *overlayCanvas()

bool SpectrumDisplayDiv::overlayCanvasEnabled() const
{
  return m_chart->overlayCanvasEnabled();
}//bool overlayCanvasEnabled() const

void SpectrumDisplayDiv::setOverlayCanvasVisible( bool visible )
{
  m_chart->setOverlayCanvasVisible( visible );
}//void setOverlayCanvasVisible( bool visible )

void SpectrumDisplayDiv::disableOverlayCanvas()
{
  m_chart->disableOverlayCanvas();
}//void disableOverlayCanvas()

void SpectrumDisplayDiv::setScrollingParent( Wt::WContainerWidget *parent )
{
  m_chart->setScrollingParent( parent );
}//void setScrollingParent( Wt::WContainerWidget *parent );

void SpectrumDisplayDiv::setScrollY( int scrollY )
{
  m_chart->setScrollY( scrollY );
}//void setScrollY( int scrollY )

Signal<> &SpectrumDisplayDiv::legendEnabled()
{
  return m_chart->legendEnabled();
}


Signal<> &SpectrumDisplayDiv::legendDisabled()
{
  return m_chart->legendDisabled();
}


void SpectrumDisplayDiv::setHidden( bool hidden, const Wt::WAnimation &anim )
{
  WContainerWidget::setHidden( hidden, anim );
  m_chart->setHidden( hidden, anim );
}//void setHidden( bool hidden, const Wt::WAnimation &animation )


void SpectrumDisplayDiv::setShowPeakLabel( int label, bool show )
{
  m_chart->setShowPeakLabel( SpectrumChart::PeakLabels(label), show );
  m_chart->update();
}//void setShowPeakLabel( int label, bool show )


bool SpectrumDisplayDiv::showingPeakLabel( int peakLabel ) const
{
  return m_chart->showingPeakLabel( SpectrumChart::PeakLabels(peakLabel) );
}//bool showingPeakLabel( int peakLabel ) const


Signal<double,double,int,int> &SpectrumDisplayDiv::chartClicked()
{
  return m_chart->leftClick();
}//Signal<double,double,int,int> &chartClicked()


Wt::Signal<double,double,int,int> &SpectrumDisplayDiv::rightClicked()
{
  return m_chart->rightClick();
}


Wt::Signal<double,double> &SpectrumDisplayDiv::doubleLeftClick()
{
  return m_chart->doubleLeftClick();
}

Wt::Signal<double,double,int,int> &SpectrumDisplayDiv::controlKeyDragged()
{
  return m_chart->controlKeyDragged();
}

Wt::Signal<double,double> &SpectrumDisplayDiv::shiftKeyDragged()
{
  return m_chart->shiftKeyDragged();
}


Wt::Signal<double,double> &SpectrumDisplayDiv::rightMouseDragg()
{
  return m_chart->rightMouseDragg();
}//Signal<double,double> &rightMouseDragg()


void SpectrumDisplayDiv::chartXRangeChangedCallback( double x, double y )
{
  if( m_autoAdjustDisplayBinnning )
    guessAndUpdateDisplayRebinFactor();

  //XXX - does are all boost slots/signals executed after the emit() statment?
  //XXX - when there is a continuum defined, the setAutoYAxisRange() call
  //      doesnt update the screen properly on the first go around
  m_chart->setAutoYAxisRange();
}//void chartXRangeChangedCallback( double x, double y )


Wt::JSlot *SpectrumDisplayDiv::alignOverlayCanvas()
{
  return m_chart->alignOverlayCanvas();
}


Wt::JSignal<std::string> *SpectrumDisplayDiv::overlayCanvasJsException()
{
  return m_chart->overlayCanvasJsException();
}


void SpectrumDisplayDiv::connectWtMouseConnections()
{
  m_chart->connectWtMouseConnections();
}


void SpectrumDisplayDiv::disconnectWtMouseConnections()
{
  m_chart->disconnectWtMouseConnections();
}


#if( IOS || ANDROID )
void SpectrumDisplayDiv::forceOverlayAlign()
{
  Wt::JSlot *align = m_chart->alignOverlayCanvas();
  if( align )
    align->exec();
}//void forceOverlayAlign()
#endif


#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
void SpectrumDisplayDiv::setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )
{
  m_chart->setReferncePhotoPeakLines( nuc );
}

void SpectrumDisplayDiv::persistCurrentReferncePhotoPeakLines()
{
  m_chart->persistCurrentReferncePhotoPeakLines();
}

void SpectrumDisplayDiv::clearAllReferncePhotoPeakLines()
{
  m_chart->clearAllReferncePhotoPeakLines();
}
#endif //#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )

void SpectrumDisplayDiv::setShowRefLineInfoForMouseOver( const bool show )
{
  m_chart->setShowRefLineInfoForMouseOver( show );
}//void setShowRefLineInfoForMouseOver( const bool show )


void SpectrumDisplayDiv::layoutSizeChanged ( int width, int height )
{
  m_layoutWidth = width;
  m_layoutHeight = height;

  if( m_autoAdjustDisplayBinnning )
  {
    guessAndUpdateDisplayRebinFactor();
    m_chart->setAutoYAxisRange();
  }//if( m_autoAdjustDisplayBinnning )
  
  
#if( IOS || ANDROID )
  //When the soft-keyboard disapears (on Android at a minimum), the overlays
  //  dont resize properly (until you change tab below the chart, or something)
  //  so we will force it.
  forceOverlayAlign();
#endif
}//void layoutSizeChanged ( int width, int height )


int SpectrumDisplayDiv::layoutWidth() const
{
  return m_layoutWidth;
}//int layoutWidth() const


int SpectrumDisplayDiv::layoutHeight() const
{
  return m_layoutHeight;
}//int layoutHeight() const


double SpectrumDisplayDiv::xAxisMinimum() const
{
  return m_chart->axis(Chart::XAxis).minimum();
}//double xAxisMinimum() const


double SpectrumDisplayDiv::xAxisMaximum() const
{
  return m_chart->axis(Chart::XAxis).maximum();
}//double xAxisMaximum() const


double SpectrumDisplayDiv::yAxisMinimum() const
{
  return m_chart->axis(Chart::YAxis).minimum();
}//double yAxisMinimum() const


double SpectrumDisplayDiv::yAxisMaximum() const
{
  return m_chart->axis(Chart::YAxis).maximum();
}//double yAxisMaximum() const


bool SpectrumDisplayDiv::yAxisIsLog() const
{
  return (m_chart->axis(Chart::YAxis).scale() == Chart::LogScale);
}//bool yAxisIsLog() const;


void SpectrumDisplayDiv::setYAxisLog( bool log )
{
  const Chart::AxisScale scale = log ? Chart::LogScale : Chart::LinearScale;
  m_chart->setYAxisScale( scale );
  m_chart->setAutoYAxisRange();
}//void setYAxisLog( bool log )


WPointF SpectrumDisplayDiv::energyCountsToPixels( double energy, double counts ) const
{
  return m_chart->mapToDevice( energy, counts );
}//Wt::WPointF energyCountsToPixels( double energy, double counts ) const


double SpectrumDisplayDiv::mapEnergyToXPixel( double energy ) const
{
  return m_chart->mapToDevice( energy, yAxisMaximum() ).x();
}


double SpectrumDisplayDiv::mapCountsToYPixel( double counts ) const
{
  return m_chart->mapToDevice( xAxisMaximum(), counts ).y();
}


void SpectrumDisplayDiv::showGridLines( bool show )
{
  m_chart->showGridLines( show );
}

void SpectrumDisplayDiv::showVerticalLines( const bool draw )
{
  m_chart->showVerticalLines( draw );
}

void SpectrumDisplayDiv::showHorizontalLines( const bool draw )
{
  m_chart->showHorizontalLines( draw );
}

bool SpectrumDisplayDiv::verticalLinesShowing() const
{
  return m_chart->verticalLinesShowing();
}

bool SpectrumDisplayDiv::horizontalLinesShowing() const
{
  return m_chart->horizontalLinesShowing();
}


bool SpectrumDisplayDiv::backgroundSubtract() const
{
  return m_model->backgroundSubtract();
}//bool backgroundSubtract() const


void SpectrumDisplayDiv::setBackgroundSubtract( bool subtract )
{
  if( subtract == m_model->backgroundSubtract() )
    return;
  
  m_model->setBackgroundSubtract( subtract );
  
#if( WT_VERSION >= 0x3030800 )
  for( auto series : m_chart->series() )
  {
    if( series->modelColumn() == SpectrumDataModel::BACKGROUND_COLUMN )
      m_chart->series(SpectrumDataModel::BACKGROUND_COLUMN).setHidden( subtract );
  }//for( size_t i = 0; i < series.size(); ++i )
#else
  const vector<Chart::WDataSeries> &series = m_chart->series();
  for( size_t i = 0; i < series.size(); ++i )
  {
    if( series[i].modelColumn() == SpectrumDataModel::BACKGROUND_COLUMN )
      m_chart->series(SpectrumDataModel::BACKGROUND_COLUMN).setHidden( subtract );
  }//for( size_t i = 0; i < series.size(); ++i )
#endif
  
  m_chart->setAutoYAxisRange();
}//void setBackgroundSubtract( bool subtract )

void SpectrumDisplayDiv::setXAxisMinimum( const double minimum )
{
  m_chart->axis(Chart::XAxis).setMinimum( minimum );
}//void setXAxisMinimum( const double minimum )


void SpectrumDisplayDiv::setXAxisMaximum( const double maximum )
{
  m_chart->axis(Chart::XAxis).setMaximum( maximum );
}//void setXAxisMaximum( const double maximum )


void SpectrumDisplayDiv::setYAxisMinimum( const double minimum )
{
  m_chart->axis(Chart::YAxis).setMinimum( minimum );
}//void setYAxisMinimum( const double minimum )


void SpectrumDisplayDiv::setYAxisMaximum( const double maximum )
{
  m_chart->axis(Chart::YAxis).setMaximum( maximum );
}//void setYAxisMaximum( const double maximum )


void SpectrumDisplayDiv::setXAxisRange( const double minimum, const double maximum )
{
  m_chart->axis(Chart::XAxis).setMinimum( minimum );
  m_chart->axis(Chart::XAxis).setMaximum( maximum );
  if( m_autoAdjustDisplayBinnning )
    m_chart->xRangeChanged().emit( minimum, maximum );
}//void setXAxisRange( const double minimum, const double maximum );


void SpectrumDisplayDiv::setYAxisRange( const double minimum,
                                        const double maximum )
{
  m_chart->axis(Chart::YAxis).setMinimum( minimum );
  m_chart->axis(Chart::YAxis).setMaximum( maximum );
}//void setYAxisRange( const double minimum, const double maximum );


void SpectrumDisplayDiv::setDisplayRebinFactor( const int factor )
{
  if( m_model->rebinFactor() != factor )
  {
    WStringStream msg;

    if( factor == 1 )
      msg << "Counts/Channel";
    else
      msg << "Counts per " << factor << " Channels";

    setYAxisTitle( msg.str() );

//    if( factor == 1 )
//      msg << "Displaying all gamma channels individually";
//    else
//      msg << "Combining every " << factor
//          << " gamma channels into 1 for display";
//    passMessage( msg.str(), "", 0 );
  }//if( m_model->rebinFactor() != factor )

  m_model->setRebinFactor( factor );
}//void setDisplayRebinFactor( const int factor )


void SpectrumDisplayDiv::setAutoAdjustDisplayRebinFactor( bool auto_rebin )
{
  m_autoAdjustDisplayBinnning = auto_rebin;
}//void setAutoAdjustDisplayRebinFactor( bool auto_rebin )


int SpectrumDisplayDiv::displayRebinFactor() const
{
  return m_model->rebinFactor();
}//int displayRebinFactor() const


void SpectrumDisplayDiv::guessAndUpdateDisplayRebinFactor()
{
  //my guess is this function houldnt have ant real effect (in terms of
  //  server/client interactions) if it is being called due to itself (e.g.
  //  calling setDisplayRebinFactor(...) with a NEW rebin factor causes
  //  the model to change dimensions, which causes the chart to resize).
  auto axisH = m_model->histUsedForXAxis();

  if( !axisH || !layoutSizeAware() )
  {
    setDisplayRebinFactor( 1 );
    return;
  }//if( !axisH )

  const float displayedxmin = static_cast<float>( m_chart->axis(Chart::XAxis).minimum() );
  const float displayedxmax = static_cast<float>( m_chart->axis(Chart::XAxis).maximum() );
  const size_t displayednbin = axisH->find_gamma_channel( displayedxmax )
                            - axisH->find_gamma_channel( displayedxmin );
  const int width = layoutWidth()
                    - m_chart->plotAreaPadding(Left)
                    - m_chart->plotAreaPadding(Right);
  const float bins_per_pixel = float(displayednbin) / float(width);
  const int rebin_factor = max( static_cast<int>(ceil(bins_per_pixel)), 1 );
  
  setDisplayRebinFactor( rebin_factor );
}//void guessAndUpdateDisplayRebinFactor()


void SpectrumDisplayDiv::setData( std::shared_ptr<Measurement> data_hist,
                                  float liveTime,
                                  float realTime,
                                  float neutronCounts,
                                  bool keep_curent_xrange )
{
  if( m_autoAdjustDisplayBinnning )
  {
    //Lets set the rebin factor before the server gets a chance to render
    //  the chart (it renders it 2 or 3 times before the user gets the final
    //  chart), to save on, at a minimum server CPU, but maybe also bandwidth.
    if( data_hist )
    {
      int rebin_factor = 1;
      const int width = layoutWidth()
                        - m_chart->plotAreaPadding(Left)
                        - m_chart->plotAreaPadding(Right);
      
      float nbin = static_cast<float>( data_hist->num_gamma_channels() );
      if( keep_curent_xrange )
      {
        const float lowerx = static_cast<float>( m_chart->axis(Chart::XAxis).minimum() );
        const float upperx = static_cast<float>( m_chart->axis(Chart::XAxis).maximum() );
        const size_t lowerbin = data_hist->find_gamma_channel( lowerx );
		    const size_t upperbin = data_hist->find_gamma_channel( upperx );
        nbin = static_cast<float>(upperbin - lowerbin);
      }//if( keep_curent_xrange )
      
      const float bins_per_pixel = nbin / static_cast<float>(width);
      rebin_factor = max( static_cast<int>(ceil(bins_per_pixel)), 1 );
      setDisplayRebinFactor( rebin_factor );
    }else
    {
      m_chart->axis(Chart::XAxis).setMinimum( 0.0 );
      m_chart->axis(Chart::XAxis).setMaximum( 3000.0 );
    }//if( data_hist ) / else
  }//if( m_autoAdjustDisplayBinnning )

  m_model->setDataHistogram( data_hist, liveTime, realTime, neutronCounts );
//  enableDisableBackgroundSubtractCB();
  //  setAutoAxisRange();

  if( !keep_curent_xrange )
    m_chart->setAutoXAxisRange();
  m_chart->setAutoYAxisRange();
  m_chart->update();
}//void setData( std::shared_ptr<Measurement> data_hist )


std::shared_ptr<Measurement> SpectrumDisplayDiv::data()
{
  return m_model->getData();
}//std::shared_ptr<Measurement> data()


std::shared_ptr<const Measurement> SpectrumDisplayDiv::data() const
{
  return m_model->getData();
}//std::shared_ptr<const Measurement> data() const


std::shared_ptr<Measurement> SpectrumDisplayDiv::secondData()
{
  return m_model->getSecondData();
}//std::shared_ptr<Measurement> secondData()


std::shared_ptr<const Measurement> SpectrumDisplayDiv::secondData() const
{
  return m_model->getSecondData();
}//std::shared_ptr<const Measurement> secondData() const


std::shared_ptr<Measurement> SpectrumDisplayDiv::background()
{
  return m_model->getBackground();
}//std::shared_ptr<Measurement> background()


std::shared_ptr<const Measurement> SpectrumDisplayDiv::background() const
{
  return m_model->getBackground();
}//std::shared_ptr<const Measurement> background() const


float SpectrumDisplayDiv::foregroundLiveTime() const
{
  return m_model->dataLiveTime();
}

float SpectrumDisplayDiv::foregroundRealTime() const
{
  return m_model->dataRealTime();
}

float SpectrumDisplayDiv::backgroundLiveTime() const
{
  return m_model->backgroundLiveTime();
}

float SpectrumDisplayDiv::backgroundRealTime() const
{
  return m_model->backgroundRealTime();
}

float SpectrumDisplayDiv::secondForegroundLiveTime() const
{
  return m_model->secondDataLiveTime();
}

float SpectrumDisplayDiv::secondForegroundRealTime() const
{
  return m_model->secondDataRealTime();
}


std::shared_ptr<const Measurement> SpectrumDisplayDiv::histUsedForXAxis() const
{
  return m_model->histUsedForXAxis();
}


void SpectrumDisplayDiv::setDisplayScaleFactor( const float sf,
                                                const SpecUtils::SpectrumType spectrum_type )
{
  switch( spectrum_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      throw runtime_error( "setDisplayScaleFactor can not be called for foreground" );

    case SpecUtils::SpectrumType::SecondForeground:
      m_model->setSecondDataScaleFactor( sf );
      break;
      
    case SpecUtils::SpectrumType::Background:
      m_model->setBackgroundDataScaleFactor( sf );
      break;
  }//switch( spectrum_type )
}//void setDisplayScaleFactor(...)


float SpectrumDisplayDiv::displayScaleFactor( const SpecUtils::SpectrumType spectrum_type ) const
{
  switch( spectrum_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      return 1.0f;
    case SpecUtils::SpectrumType::SecondForeground:
      return m_model->secondDataScaledBy();
    case SpecUtils::SpectrumType::Background:
      return m_model->backgroundScaledBy();
//  m_spectrumDiv->continuum();
  }//switch( spectrum_type )

  throw runtime_error( "SpectrumDisplayDiv::displayScaleFactor(...): invalid input arg" );

  return 1.0;
}//double displayScaleFactor( SpecUtils::SpectrumType spectrum_type ) const;


void SpectrumDisplayDiv::setBackground( std::shared_ptr<Measurement> background,
                                        float liveTime,
                                        float realTime,
                                        float neutronCounts )
{
  const bool hadBackground = !!(m_model->getBackground());

  m_model->setBackgroundHistogram( background, liveTime, realTime, neutronCounts );
//  enableDisableBackgroundSubtractCB();

  //XXX 2011-10-23 - I dont understand what I was doing below, it should be checked
  const bool hasBackground = !!(m_model->getBackground());

  if( hadBackground != hasBackground )
  {
    const vector<Chart::WDataSeries> series = m_model->suggestDataSeries();
    m_chart->setSeries( series );
    m_chart->axis(Chart::Y2Axis).setVisible( m_model->secondDataOwnAxis() );
  }//if( !hadBackground )


  if( !background && m_model->backgroundSubtract() )
    m_model->setBackgroundSubtract( false );
  
  m_chart->update();
}//void SpectrumDisplayDiv::setBackground(...);


void SpectrumDisplayDiv::setSecondData( std::shared_ptr<Measurement> hist,
                                        float liveTime,
                                        float realTime,
                                        float neutronCounts,
                                        bool ownAxis )
{
  const bool alreadyHad = !!(m_model->getSecondData());

  m_model->setSecondDataHistogram( hist, liveTime, realTime, neutronCounts, ownAxis );

  const bool nowHas = !!(m_model->getSecondData());

  if( alreadyHad != nowHas )
  {
    const vector<Chart::WDataSeries> series = m_model->suggestDataSeries();
    m_chart->setSeries( series );
    m_chart->axis(Chart::Y2Axis).setVisible( m_model->secondDataOwnAxis() );
  }//if( !hadBackground )
  
  m_chart->update();
}//void SpectrumDisplayDiv::setSecondData( std::shared_ptr<Measurement> background );




void SpectrumDisplayDiv::clearTimeHighlightRegions( const SpecUtils::SpectrumType type )
{
  SpectrumChart::HighlightRegionType rtype;
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      rtype = SpectrumChart::ForegroundHighlight;
      break;
    case SpecUtils::SpectrumType::Background:
      rtype = SpectrumChart::BackgroundHighlight;
      break;
    case SpecUtils::SpectrumType::SecondForeground:
      rtype = SpectrumChart::SecondForegroundHighlight;
      break;
  }//switch( type )
  
  m_chart->clearTimeHighlightRegions(rtype);
  m_chart->update();
}//void clearTimeHighlightRegions();


void SpectrumDisplayDiv::setTimeHighLightRegions( const vector< pair<double,double> > &p,
                                             const SpecUtils::SpectrumType type )
{
  SpectrumChart::HighlightRegionType rtype;
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      rtype = SpectrumChart::ForegroundHighlight;
      break;
    case SpecUtils::SpectrumType::Background:
      rtype = SpectrumChart::BackgroundHighlight;
      break;
    case SpecUtils::SpectrumType::SecondForeground:
      rtype = SpectrumChart::SecondForegroundHighlight;
      break;
  }//switch( type )
  
  m_chart->setTimeHighLightRegions( p, rtype );
  m_chart->update();
}//void setHighLightRegions(...)


void SpectrumDisplayDiv::clearOccupancyRegions()
{
  m_chart->clearOccupancyRegions();
  m_chart->update();
}//void clearTimeUnderlineRegions( const SpecUtils::SpectrumType type )


void SpectrumDisplayDiv::setOccupancyRegions( const std::vector< std::pair<double,double> > &p )
{
  m_chart->setOccupancyRegions( p );
  m_chart->update();
}//void setTimeUnderlineRegions(...)



bool SpectrumDisplayDiv::removeDecorativeHighlightRegion( size_t uniqueid )
{
  const bool success = m_chart->removeDecorativeHighlightRegion( uniqueid );
  if( success )
    m_chart->update();
  return success;
}//void removeDecorativeHighlightRegions()


size_t SpectrumDisplayDiv::addDecorativeHighlightRegion( const float lowerx,
                                  const float upperx,
                                  const Wt::WColor &color )
{
  const size_t regionid = m_chart->addDecorativeHighlightRegion( lowerx, upperx, color );
  m_chart->update();
  return regionid;
}//void addDecorativeHighlightRegion(...)


void SpectrumDisplayDiv::setAutoAxisRange()
{
  m_chart->setAutoXAxisRange();
  m_chart->setAutoYAxisRange();
}//void setAutoAxisRange()



void SpectrumDisplayDiv::setPlotAreaPadding( int left, int top, int right, int bottom )
{
  m_chart->setPlotAreaPadding( left, Left );
  m_chart->setPlotAreaPadding( top, Top );
  m_chart->setPlotAreaPadding( right, Right );
  m_chart->setPlotAreaPadding( bottom, Bottom );
}//void setPlotAreaPadding( int left, int top, int right, int bottom )


void SpectrumDisplayDiv::visibleRange( double &xmin, double &xmax,
                                        double &ymin, double &ymax ) const
{
  m_chart->visibleRange( xmin, xmax, ymin, ymax );
}

int SpectrumDisplayDiv::plotAreaPadding( const Wt::Side side ) const
{
  return m_chart->plotAreaPadding( side );
}//int plotAreaPadding( const Wt::Side side ) const


const string SpectrumDisplayDiv::xAxisTitle() const
{
  const WString &title = m_chart->axis(Chart::XAxis).title();
  if( title.empty() )
    return "";
  return title.toUTF8();
}//const Wt::WString &xAxisTitle() const;


const string SpectrumDisplayDiv::yAxisTitle() const
{
  const WString &title = m_chart->axis(Chart::YAxis).title();
  if( title.empty() )
    return "";
  return title.toUTF8();
}//const Wt::WString &yAxisTitle() const;


const string SpectrumDisplayDiv::y2AxisTitle() const
{
  const WString &title = m_chart->axis(Chart::Y2Axis).title();
  if( title.empty() )
    return "";
  return title.toUTF8();
}//const Wt::WString &y2AxisTitle() const;


void SpectrumDisplayDiv::setXAxisTitle( const std::string &title )
{
  m_chart->axis(Chart::XAxis).setTitle( title );
}//void setXAxisTitle( const std::string &title )


void SpectrumDisplayDiv::setYAxisTitle( const std::string &title )
{
  m_chart->axis(Chart::YAxis).setTitle( title );
}//void setYAxisTitle( const std::string &title )


void SpectrumDisplayDiv::setY2AxisTitle( const std::string &title )
{
  m_chart->axis(Chart::Y2Axis).setTitle( title );
}//void setY2AxisTitle( const std::string &title )

float SpectrumDisplayDiv::xUnitsPerPixel() const
{
  return m_chart->xUnitsPerPixel();
}

void SpectrumDisplayDiv::setForegroundSpectrumColor( const Wt::WColor &color )
{
  m_model->setForegroundSpectrumColor( color );
  m_chart->setSeries( m_model->suggestDataSeries() );
}

void SpectrumDisplayDiv::setBackgroundSpectrumColor( const Wt::WColor &color )
{
  m_model->setBackgroundSpectrumColor( color );
  m_chart->setSeries( m_model->suggestDataSeries() );
}

void SpectrumDisplayDiv::setSecondarySpectrumColor( const Wt::WColor &color )
{
  m_model->setSecondarySpectrumColor( color );
  m_chart->setSeries( m_model->suggestDataSeries() );
}

void SpectrumDisplayDiv::setDefaultPeakColor( const Wt::WColor &color )
{
  m_chart->setDefaultPeakColor( color );
}

void SpectrumDisplayDiv::setAxisLineColor( const Wt::WColor &color )
{
  m_chart->setAxisLineColor(color);
}

void SpectrumDisplayDiv::setChartMarginColor( const Wt::WColor &color )
{
  m_chart->setChartMarginColor(color);
}

void SpectrumDisplayDiv::setChartBackgroundColor( const Wt::WColor &color )
{
  m_chart->setChartBackgroundColor(color);
}

void SpectrumDisplayDiv::setTextColor( const Wt::WColor &color )
{
  m_chart->setTextColor(color);
}

void SpectrumDisplayDiv::setOccupiedTimeSamplesColor( const Wt::WColor &color )
{
  m_chart->setOccupiedTimeSamplesColor( color );
}

void SpectrumDisplayDiv::setForegroundHighlightColor( const Wt::WColor &color )
{
  m_chart->setForegroundHighlightColor( color );
}


void SpectrumDisplayDiv::setBackgroundHighlightColor( const Wt::WColor &color )
{
  m_chart->setBackgroundHighlightColor( color );
}


void SpectrumDisplayDiv::setSecondaryHighlightColor( const Wt::WColor &color )
{
  m_chart->setSecondaryHighlightColor( color );
}


void SpectrumDisplayDiv::setMouseDragZooms()
{
  m_chart->setDragAction( SpectrumChart::ZoomIn );
}//void SpectrumDisplayDiv::setMouseDragZooms()


void SpectrumDisplayDiv::setMouseDragHighlights( const bool allowMultiple,
                                              const bool allowSingleClick )
{
  m_chart->setDragAction( SpectrumChart::HighLight );
  m_chart->allowMultipleRegionHighlight( allowMultiple );
  m_chart->allowSingleClickHighlight( allowSingleClick );
}//void setMouseDragHighlights( const bool allowMultiple )


void SpectrumDisplayDiv::allowArrowToMoveSingleClickRegion( bool allow )
{
  m_chart->allowArrowToMoveSingleClickRegion( allow );
}//void allowArrowToMoveSingleClickRegion( allow )


void SpectrumDisplayDiv::disableMouseDragActions()
{
  m_chart->setDragAction( SpectrumChart::NoAction );
}//void disableMouseDragActions()


Wt::Signal<double,double> &SpectrumDisplayDiv::xRangeChanged()
{
  return m_chart->xRangeChanged();
}//xRangeChanged()


Wt::Signal<double,double> &SpectrumDisplayDiv::altKeyDragged()
{
  return m_chart->altKeyDragged();
}

Wt::Signal<double,double> &SpectrumDisplayDiv::shiftAltKeyDragged()
{
  return m_chart->shiftAltKeyDragged();
}

Wt::Signal<double,double> &SpectrumDisplayDiv::controlMouseMoved()
{
  return m_chart->controlMouseMoved();
}

SpectrumDisplayDiv::~SpectrumDisplayDiv()
{
}//~SpectrumDisplayDiv()


