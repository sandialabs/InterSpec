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

#include <Wt/WPoint>
#include <Wt/WServer>
#include <Wt/WLength>
#include <Wt/WIOService>
#include <Wt/WJavaScript>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace Wt;
using namespace std;

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
  
  //DeleteOnClosePopupMenu is a PopupDivMenu that deletes itself on close.
  //  Necassary because (with Wt 3.3.4 at least) using the aboutToHide() signal
  //  to delete the menu causes a crash.
  class DeleteOnClosePopupMenu : public PopupDivMenu
  {
    bool m_deleteWhenHidden;
  public:
    DeleteOnClosePopupMenu( WPushButton *p, const PopupDivMenu::MenuType t )
    : PopupDivMenu( p, t ), m_deleteWhenHidden( false ) {}
    virtual ~DeleteOnClosePopupMenu(){}
    void markForDelete(){ m_deleteWhenHidden = true; }
    virtual void setHidden( bool hidden, const WAnimation &a = WAnimation() )
    {
      PopupDivMenu::setHidden( hidden, a );
      if( hidden && m_deleteWhenHidden )
        delete this;
    }
  };//class PeakRangePopupMenu
  
}//namespace


WT_DECLARE_WT_MEMBER
(SvgToImgDownload, Wt::JavaScriptFunction, "SvgToImgDownload",
 function(chart,filename,asPng)
{
  // To copy to the pasteboard, and use the Clipboard API, an https connection is required, so we would have to send the image blob back to C++-land, and then put it in the clipboard there; see git commit 0db92a30928a3d8 for the sketched out Clipboard API try.
  try
  {
    // 'chart' is the JS SpectrumChartD3 object
    // 'chart.chart' is the <div> that holds the spectrum
    //  and chart.svg is a D3 array with one entry to the actual <svg> element
    //console.log( 'In SvgToImgDownload: svgchart' );
    if( filename.length === 0 )
      filename = "spectrum." + (asPng ? "png" : "svg"); //ToDo: add in date/time or something.
    
    let svgMarkup = chart.getStaticSvg();
    
    // A heavy-duty version of element.click()
    function simulateClick( node ) {
      try {
        node.dispatchEvent(new MouseEvent('click'));
      } catch (e) {
        var evt = document.createEvent('MouseEvents');
        evt.initMouseEvent('click', true, true, window,
                            0, 0, 0, 100, 100,
                            false, false, false, false, 0, null);
        node.dispatchEvent(evt);
      }
    };// function simulateClick()
    
    if( !asPng )
    {
      var dl = document.createElement("a");
      document.body.appendChild(dl);
      dl.target = "_blank";
      dl.setAttribute("href", "data:image/svg+xml," + encodeURIComponent(svgMarkup) );
      dl.setAttribute("download", filename);
      simulateClick(dl);
      document.body.removeChild(dl);
    }else
    {
      // The getStaticSvg() function does not include the energy slider chart, so we should
      //  compensate for that height if its showing
      let height = chart.svg.attr("height");
      let width = chart.svg.attr("width");
      
      if( chart.sliderChart && chart.size && chart.size.sliderChartHeight )
        height -= chart.size.sliderChartHeight;
      
      if( chart.scalerWidget )
        width -= chart.scalerWidget.node().getBBox().width;
      
      var canvas = document.createElement("canvas");
      canvas.width = width;
      canvas.height = height;
      chart.chart.appendChild(canvas);
    
      var ctx = canvas.getContext("2d");
      
      // The default style in InterSpec doesnt have a background for the SVG (but the dark theme
      //  does), so start with a solid white background instead of no fill (the default PNG viewer
      //  in Windows fills in a solid black background in this case, making the PNG look incorrect)
      ctx.fillStyle = 'white';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      
      var img = new Image();
      var svg = new Blob([svgMarkup], {type: "image/svg+xml"});
      var url = window.URL.createObjectURL(svg);
      
      // JIC something is messed up we can see it in debug console for development
      img.onerror = function(e){ console.log( 'img.error: '); console.log(e); };
      
      img.onload = function() {
        ctx.drawImage(img, 0, 0);
        window.URL.revokeObjectURL(url);
      
        var dt = canvas.toDataURL('image/png');
        dt = dt.replace(/^data:image\\/[^;]*/, 'data:application/octet-stream');
        dt = dt.replace(/^data:application\\/octet-stream/, 'data:application/octet-stream;headers=Content-Disposition%3A%20attachment%3B%20filename='+filename);
        canvas.remove();
      
        var link = document.createElement("a");
        link.style.display = "none";
        document.body.appendChild(link);
        link.download = filename;
        link.href = dt;
        link.target="_blank";
        simulateClick(link);
      
        window.URL.revokeObjectURL(link.href);
        document.body.removeChild(link);
      };
      img.src = url;
    }//if( !asPng ) / else
    
  }catch(e){
    console.log( 'Error saving PNG or SVG of spectrum: ' );
    console.log( e );
  }
}
);




D3SpectrumDisplayDiv::D3SpectrumDisplayDiv( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_renderFlags( 0 ),
  m_peakModel( nullptr ),
  m_foreground( nullptr ),
  m_secondary( nullptr ),
  m_background( nullptr ),
  m_secondaryScale( 1.0f ),
  m_backgroundScale( 1.0f ),
  m_compactAxis( false ),
  m_legendEnabled( true ),
  m_yAxisIsLog( true ),
  m_backgroundSubtract( false ),
  m_showVerticalLines( false ),
  m_showHorizontalLines( false ),
  m_showHistogramIntegralsInLegend( true ),
  m_showXAxisSliderChart( false ),
  m_compactXAxisWithSliderChart( true ),
  m_showYAxisScalers( false ),
  m_searchEnergies{},
  m_highlights{},
  m_peakLabelsToShow{},
  m_xAxisTitle( WString::tr("Energy (keV)") ),
  m_yAxisTitle( WString::tr("d3sdd-yAxisTitle") ),
  m_yAxisTitleMulti( WString::tr("d3sdd-yAxisTitleMulti") ),
  m_chartTitle(),
  // These next member vairables are all `std::unique_ptr` to `Wt::JSignal`s
  m_shiftKeyDraggJS{ nullptr },
  m_shiftAltKeyDraggJS{ nullptr },
  m_rightMouseDraggJS{ nullptr },
  m_doubleLeftClickJS{ nullptr },
  m_leftClickJS{ nullptr },
  m_rightClickJS{ nullptr },
  m_xRangeChangedJS{ nullptr },
  m_existingRoiEdgeDragJS{ nullptr },
  m_dragCreateRoiJS{ nullptr },
  m_yAxisDraggedJS{ nullptr },
  m_legendClosedJS{ nullptr },
  m_sliderDisplayed{ nullptr },
  m_yAxisTypeChanged{ nullptr },
  //These next member variables are all the C++ signals for when events happen
  m_legendEnabledSignal( this ),
  m_legendDisabledSignal( this ),
  m_xRangeChanged( this ),
  m_controlKeyDragg( this ),
  m_shiftKeyDragg( this ),
  m_shiftAltKeyDragg( this ),
  m_rightMouseDragg( this ),
  m_leftClick( this ),
  m_doubleLeftClick( this ),
  m_rightClick( this ),
  m_existingRoiEdgeDrag( this ),
  m_dragCreateRoi( this ),
  m_yAxisScaled( this ),
  m_xAxisSliderShown( this ),
  m_yAxisLogLinChanged( this ),
  m_xAxisCompactnessChanged( this ),
  // These next member variables roughly track state in the JS
  m_jsgraph( jsRef() + ".chart" ),
  m_xAxisMinimum(0.0),
  m_xAxisMaximum(0.0),
  m_yAxisMinimum(0.0),
  m_yAxisMaximum(0.0),
  m_chartWidthPx(0.0),
  m_chartHeightPx(0.0),
  m_showRefLineInfoForMouseOver( true ),
  m_comptonPeakAngle( 180 ),
  m_foregroundLineColor( 0x00, 0x00, 0x00 ),  //black
  m_backgroundLineColor( 0x00, 0xff, 0xff ),  //cyan
  m_secondaryLineColor( 0x00, 0x80, 0x80 ),   //dark green
  m_textColor( 0x00, 0x00, 0x00 ),
  m_axisColor( 0x00, 0x00, 0x00 ),
  m_chartMarginColor(),
  m_chartBackgroundColor(),
  m_defaultPeakColor( 0, 51, 255, 155 ),
  m_peakLabelFontSize{},
  m_peakLabelRotationDegrees(0.0),
  m_logYAxisMin(0.1),
  m_cssRules{},
  m_pendingJs{},
  m_last_drag_time{},
  m_continuum_being_drug( nullptr ),
  m_last_being_drug_peaks()
{
  addStyleClass( "D3SpectrumDisplayDiv" );
  
  // Cancel right-click events for the div, we handle it all in JS
  setAttributeValue( "oncontextmenu",
                     "event.cancelBubble = true; event.returnValue = false; return false;" );
  
  InterSpec::instance()->useMessageResourceBundle( "D3SpectrumDisplayDiv" );
    
  //For development it may be useful to directly use the original JS/CSS files,
  //  but normally we should use the resources CMake will copy into
  //  InterSpec_resources (not checked that Andorid build systems will
  //  grab these files); when developing note that CMake will only update
  //  files to InterSpec_resources when you run the "make" command.
#if( defined(NDEBUG) || IOS || ANDROID || BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
  //THe NDEBUG should be enough, but just making sure
  const string resource_base = "InterSpec_resources/";
#else
  const string resource_base = "external_libs/SpecUtils/d3_resources/";
#endif
  
  wApp->useStyleSheet( resource_base + "SpectrumChartD3.css" );
  initChangeableCssRules();
  
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( resource_base + "SpectrumChartD3.js" );
  
  
  for( SpectrumChart::PeakLabels label = SpectrumChart::PeakLabels(0);
      label < SpectrumChart::PeakLabels::kNumPeakLabels; label = SpectrumChart::PeakLabels(label+1) )
  {
    m_peakLabelsToShow[label] = false;
  }//for( loop over all labels )
  
  for( FeatureMarkerType type = FeatureMarkerType(0);
      type < FeatureMarkerType::NumFeatureMarkers;
      type = FeatureMarkerType( static_cast<int>(type) + 1) )
  {
    m_showFeatureMarker[static_cast<int>(type)] = false;
  }
}//D3SpectrumDisplayDiv constructor


void D3SpectrumDisplayDiv::defineJavaScript()
{
  string options = "{title: '', showAnimation: true, animationDuration: 200";
  options += ", showSliderCloseBtn: true";
  options += ", compactXAxis: " + jsbool(m_compactAxis);
  options += ", allowDragRoiExtent: true";
  options += ", showRefLineInfoForMouseOver: " + jsbool(m_showRefLineInfoForMouseOver);
  options += ", yscale: " + string(m_yAxisIsLog ? "'log'" : "'lin'");
  options += ", backgroundSubtract: " + jsbool( m_backgroundSubtract );
  options += ", showLegend: " + jsbool(m_legendEnabled);
  options += ", gridx: " + jsbool(m_showVerticalLines);
  options += ", gridy: " + jsbool(m_showHorizontalLines);
  options += ", showXAxisSliderChart: " + jsbool(m_showXAxisSliderChart);
  options += ", compactXAxisWithSliderChart: " + jsbool(m_compactXAxisWithSliderChart);
  options += ", scaleBackgroundSecondary: " + jsbool(m_showYAxisScalers);
  options += ", wheelScrollYAxis: true";
  options += ", sliderChartHeightFraction: 0.1";  //ToDo: track this in C++
  options += ", spectrumLineWidth: 1.0";  //ToDo: Let this be specified in C++
  options += ", showUserLabels: " + jsbool(m_peakLabelsToShow[SpectrumChart::kShowPeakUserLabel]);
  options += ", showPeakLabels: " + jsbool(m_peakLabelsToShow[SpectrumChart::kShowPeakEnergyLabel]);
  options += ", showNuclideNames: " + jsbool(m_peakLabelsToShow[SpectrumChart::kShowPeakNuclideLabel]);
  options += ", showNuclideEnergies: " + jsbool(m_peakLabelsToShow[SpectrumChart::kShowPeakNuclideEnergies]);
  
  for( FeatureMarkerType type = FeatureMarkerType(0);
      type < FeatureMarkerType::NumFeatureMarkers;
      type = FeatureMarkerType(static_cast<int>(type) + 1) )
  {
    switch( type )
    {
      case FeatureMarkerType::EscapePeakMarker:
        options += ", showEscapePeaks:" + jsbool(m_showFeatureMarker[static_cast<int>(type)]);
        break;
      case FeatureMarkerType::ComptonPeakMarker:
        options += ", showComptonPeaks:" + jsbool(m_showFeatureMarker[static_cast<int>(type)]);
        break;
      case FeatureMarkerType::ComptonEdgeMarker:
        options += ", showComptonPeaks:" + jsbool(m_showFeatureMarker[static_cast<int>(type)]);
        break;
      case FeatureMarkerType::SumPeakMarker:
        options += ", showSumPeaks:" + jsbool(m_showFeatureMarker[static_cast<int>(type)]);
        break;
      case FeatureMarkerType::NumFeatureMarkers:
      default:
        return;
    }//switch( option )
  }//for( loop over FeatureMarkerType's )
  options += ", comptonPeakAngle:" + std::to_string(m_comptonPeakAngle);
  
  const string labelFontSize = m_peakLabelFontSize.empty() ? string("null") : Wt::WWebWidget::jsStringLiteral(m_peakLabelFontSize,'\'');
  options += ", peakLabelSize: " + labelFontSize;
  options += ", peakLabelRotation: " + std::to_string(m_peakLabelRotationDegrees);
  if( m_logYAxisMin > 0 )
    options += ", logYAxisMin: " + std::to_string(m_logYAxisMin);
  
//Reading the configuration property from wt_config.xml doesnt seem to work; didnt spend time to debug
//  string double_click_time;
//  if( wApp->readConfigurationProperty( "double-click-timeout", double_click_time ) )
//    options += ", doubleClickDelay:" + double_click_time;
  options += ", txt: " + localizedStringsJson();
  options += "}";
  
  setJavaScriptMember( "chart", "new SpectrumChartD3(" + jsRef() + "," + options + ");");
  
  setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + id() + "') ){"
          // When we "Clear Session", jsRef() will give a null result temporarily, so we'll protect against that
          "const c=" + jsRef() + ";"
          "if(c && c.chart)"
            "c.chart.handleResize();"
        "}"
      "}"
    "});"
  );
  
  callJavaScriptMember( "resizeObserver.observe", jsRef() );
  
  setSearchEnergies( m_searchEnergies );
  
  if( !m_xRangeChangedJS )
  {
    m_xRangeChangedJS.reset( new JSignal<double,double,double,double,bool>( this, "xrangechanged", true ) );
    m_xRangeChangedJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartXRangeChangedCallback, this,
                                            boost::placeholders::_1, boost::placeholders::_2,
                                            boost::placeholders::_3, boost::placeholders::_4,
                                            boost::placeholders::_5) );
    
    m_shiftKeyDraggJS.reset( new JSignal<double,double>( this, "shiftkeydragged", true ) );
    m_shiftKeyDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartShiftKeyDragCallback, this,
                                            boost::placeholders::_1, boost::placeholders::_2 ) );
    
    m_shiftAltKeyDraggJS.reset( new JSignal<double,double>( this, "shiftaltkeydragged", true ) );
    m_shiftAltKeyDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartShiftAltKeyDragCallback,
                                               this, boost::placeholders::_1, boost::placeholders::_2 ) );
    
    m_rightMouseDraggJS.reset( new JSignal<double,double>( this, "rightmousedragged", true ) );
    m_rightMouseDraggJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartRightMouseDragCallback,
                                              this, boost::placeholders::_1, boost::placeholders::_2 ) );
    
    m_leftClickJS.reset( new JSignal<double,double,double,double>( this, "leftclicked", true ) );
    m_leftClickJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartLeftClickCallback, this,
                                        boost::placeholders::_1, boost::placeholders::_2,
                                        boost::placeholders::_3, boost::placeholders::_4 ) );
    
    m_doubleLeftClickJS.reset( new JSignal<double,double>( this, "doubleclicked", true ) );
    m_doubleLeftClickJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartDoubleLeftClickCallback,
                                              this, boost::placeholders::_1, boost::placeholders::_2 ) );
    
    m_rightClickJS.reset( new JSignal<double,double,double,double>( this, "rightclicked", true ) );
    m_rightClickJS->connect( boost::bind( &D3SpectrumDisplayDiv::chartRightClickCallback, this,
                                         boost::placeholders::_1, boost::placeholders::_2,
                                         boost::placeholders::_3, boost::placeholders::_4 ) );
    
    m_yAxisDraggedJS.reset( new Wt::JSignal<double,std::string>( this, "yscaled", true ) );
    m_yAxisDraggedJS->connect( boost::bind( &D3SpectrumDisplayDiv::yAxisScaled, this,
                                           boost::placeholders::_1, boost::placeholders::_2 ) );
    
     
    m_existingRoiEdgeDragJS.reset( new JSignal<double,double,double,double,double,bool>( this, "roiDrag", true ) );
    m_existingRoiEdgeDragJS->connect( boost::bind( &D3SpectrumDisplayDiv::existingRoiEdgeDragCallback, this,
                                         boost::placeholders::_1, boost::placeholders::_2,
                                         boost::placeholders::_3, boost::placeholders::_4,
                                         boost::placeholders::_5, boost::placeholders::_6 ) );
    
    m_dragCreateRoiJS.reset( new JSignal<double,double,int,bool,double,double>( this, "fitRoiDrag", true ) );
    m_dragCreateRoiJS->connect( boost::bind( &D3SpectrumDisplayDiv::dragCreateRoiCallback, this,
                                         boost::placeholders::_1, boost::placeholders::_2,
                                         boost::placeholders::_3, boost::placeholders::_4,
                                         boost::placeholders::_5, boost::placeholders::_6 ) );
    
    //need legend closed signal.
    m_legendClosedJS.reset( new JSignal<>( this, "legendClosed", true ) );
    m_legendClosedJS->connect( std::bind( [this](){
      m_legendEnabled = false;
      m_legendDisabledSignal.emit();
    }) );
    
    
    m_sliderDisplayed.reset( new Wt::JSignal<bool>(this, "sliderChartDisplayed") );
    m_sliderDisplayed->connect( boost::bind( &D3SpectrumDisplayDiv::sliderChartDisplayedCallback, this, boost::placeholders::_1 ) );
    
    m_yAxisTypeChanged.reset( new Wt::JSignal<string>(this, "yAxisTypeChanged") );
    m_yAxisTypeChanged->connect( boost::bind( &D3SpectrumDisplayDiv::yAxisTypeChangedCallback, this, boost::placeholders::_1 ) );
  }//if( !m_xRangeChangedJS )
  
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
  
  //I think x and y ranges should be taken care of via m_pendingJs... untested
  //m_xAxisMinimum, m_xAxisMaximum, m_yAxisMinimum, m_yAxisMaximum;
}//void defineJavaScript()


void D3SpectrumDisplayDiv::initChangeableCssRules()
{
  WCssStyleSheet &style = wApp->styleSheet();
  m_cssRules["GridColor"] = style.addRule( ":root", "--d3spec-grid-major-color: #b3b3b3;" );      //Not actually needed, as CSS will default to this
  m_cssRules["MinorGridColor"] = style.addRule( ":root", "--d3spec-grid-minor-color: #e6e6e6;" ); //Not actually needed, as CSS will default to this
}//void initChangeableCssRules()


void D3SpectrumDisplayDiv::setTextInMiddleOfChart( const Wt::WString &s )
{
  // TODO: D3 chart currently does not have functionality to set text in the middle of chart...
}

void D3SpectrumDisplayDiv::setCompactAxis( const bool compact )
{
  if( compact == m_compactAxis )
    return;
  
  m_compactAxis = compact;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setCompactXAxis(" + jsbool(compact) + ");" );
  
  m_xAxisCompactnessChanged.emit( compact );
}

bool D3SpectrumDisplayDiv::isAxisCompacted() const
{
  return m_compactAxis;
}

void D3SpectrumDisplayDiv::setPeakModel( PeakModel *model )
{
  if( !model )
    throw runtime_error( "setPeakModel(...): invalid input model" );
  
  m_peakModel = model;

  m_peakModel->dataChanged().connect( this, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->rowsRemoved().connect( this, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->rowsInserted().connect( this, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->layoutChanged().connect( this, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->modelReset().connect( this, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
}//void setPeakModel( PeakModel *model );


bool D3SpectrumDisplayDiv::legendIsEnabled() const
{
  return m_legendEnabled;
}


void D3SpectrumDisplayDiv::enableLegend()
{
  m_legendEnabled = true;
  m_legendEnabledSignal.emit();
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setShowLegend(true);" );
}//void D3SpectrumDisplayDiv::enableLegend()


void D3SpectrumDisplayDiv::disableLegend()
{
  m_legendEnabled = false;
  m_legendDisabledSignal.emit();
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setShowLegend(false);" );
}//void disableLegend()


void D3SpectrumDisplayDiv::showHistogramIntegralsInLegend( const bool show )
{
  m_showHistogramIntegralsInLegend = show;
  // TODO: No option in D3 to show/hide histogram integrals in legend
}


Signal<> &D3SpectrumDisplayDiv::legendEnabled()
{
  return m_legendEnabledSignal;
}


Signal<> &D3SpectrumDisplayDiv::legendDisabled()
{
  return m_legendDisabledSignal;
}


void D3SpectrumDisplayDiv::setShowPeakLabel( int label, bool show )
{
  SpectrumChart::PeakLabels peakLabel = SpectrumChart::PeakLabels(label);
  
  string js = m_jsgraph;
  switch( peakLabel )
  {
    case SpectrumChart::PeakLabels::kShowPeakUserLabel:
      js += ".setShowUserLabels";
      break;
    case SpectrumChart::PeakLabels::kShowPeakEnergyLabel:
      js += ".setShowPeakLabels";
      break;
    case SpectrumChart::PeakLabels::kShowPeakNuclideLabel:
      js += ".setShowNuclideNames";
      break;
    case SpectrumChart::PeakLabels::kShowPeakNuclideEnergies:
      js += ".setShowNuclideEnergies";
      break;
    case SpectrumChart::PeakLabels::kNumPeakLabels:
    default:
      return;
      break;
  }
  js += "(" + jsbool(show) + ");";
  
  m_peakLabelsToShow[peakLabel] = show;
  
  if( isRendered() )
    doJavaScript( js );
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

Wt::Signal<double,double,double,double,double,bool> &D3SpectrumDisplayDiv::existingRoiEdgeDragUpdate()
{
  return m_existingRoiEdgeDrag;
}

Wt::Signal<double, double, int, bool, double, double> &D3SpectrumDisplayDiv::dragCreateRoiUpdate()
{
  return m_dragCreateRoi;
}


Wt::Signal<double,double,SpecUtils::SpectrumType> &D3SpectrumDisplayDiv::yAxisScaled()
{
  return m_yAxisScaled;
}

Wt::Signal<double,double> &D3SpectrumDisplayDiv::doubleLeftClick()
{
  return m_doubleLeftClick;
}

Wt::Signal<double,double> &D3SpectrumDisplayDiv::shiftKeyDragged()
{
  return m_shiftKeyDragg;
}


Wt::Signal<double,double> &D3SpectrumDisplayDiv::rightMouseDragg()
{
  return m_rightMouseDragg;
}//Signal<double,double> &rightMouseDragg()


Wt::Signal<bool> &D3SpectrumDisplayDiv::xAxisSliderShown()
{
  return m_xAxisSliderShown;
}


Wt::Signal<bool> &D3SpectrumDisplayDiv::yAxisLogLinChanged()
{
  return m_yAxisLogLinChanged;
}


Wt::Signal<bool> &D3SpectrumDisplayDiv::xAxisCompactnessChanged()
{
  return m_xAxisCompactnessChanged;
}


void D3SpectrumDisplayDiv::setReferncePhotoPeakLines( const ReferenceLineInfo &nuc )
{
  m_referencePhotoPeakLines = nuc;
  m_renderFlags |= D3RenderActions::UpdateRefLines;
  scheduleRender();
}

void D3SpectrumDisplayDiv::persistCurrentReferncePhotoPeakLines()
{
  const string input_txt = m_referencePhotoPeakLines.m_input.m_input_txt;
  if( input_txt.empty() )
    return;
  
  ReferenceLineInfo *prev_line = nullptr;
  for( size_t i = 0; !prev_line && (i < m_persistedPhotoPeakLines.size()); ++i )
  {
    if( m_persistedPhotoPeakLines[i].m_input.m_input_txt == input_txt )
      prev_line = &(m_persistedPhotoPeakLines[i]);
  }
  
  if( prev_line )
    *prev_line = m_referencePhotoPeakLines;
  else if( m_referencePhotoPeakLines.m_validity == ReferenceLineInfo::InputValidity::Valid )
    m_persistedPhotoPeakLines.push_back( m_referencePhotoPeakLines );
  
  m_referencePhotoPeakLines.reset();
  
  m_renderFlags |= D3RenderActions::UpdateRefLines;
  scheduleRender();
}

void D3SpectrumDisplayDiv::clearAllReferncePhotoPeakLines()
{
  m_referencePhotoPeakLines.reset();
  m_persistedPhotoPeakLines.clear();
  
  if( isRendered() )
    doJavaScript(m_jsgraph + ".clearReferenceLines();");
}


void D3SpectrumDisplayDiv::setReferenceLinesToClient()
{
  string result = "[";
  const ReferenceLineInfo &showingNuclide = m_referencePhotoPeakLines;
  
  bool addComma = false;
  if( !showingNuclide.m_ref_lines.empty() )
  {
    showingNuclide.toJson(result);
    addComma = true;
  }
  
  for (const ReferenceLineInfo &ref : m_persistedPhotoPeakLines)
  {
    if( showingNuclide.m_ref_lines.empty()
        || (ref.m_input.m_input_txt != showingNuclide.m_input.m_input_txt) )
    {
      if( addComma )
        result += ",";
      addComma = true;
      ref.toJson(result);
    }
  }
  result += "]";
  
  const string js =
  "try{"
    + m_jsgraph + ".setReferenceLines(" + result + ");"
  "}catch(e){ console.log('Exception setting ref lines: ' + e ); }";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setReferenceLinesToClient()


void D3SpectrumDisplayDiv::setShowRefLineInfoForMouseOver( const bool show )
{
  m_showRefLineInfoForMouseOver = show;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setShowRefLineInfoForMouseOver(" + jsbool(show) + ")" );
}//void setShowRefLineInfoForMouseOver( const bool show )


void D3SpectrumDisplayDiv::highlightPeakAtEnergy( const double energy )
{
  if( isRendered() )
    doJavaScript( m_jsgraph + ".highlightPeakAtEnergy(" + std::to_string(energy) + ")" );
}


void D3SpectrumDisplayDiv::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  WContainerWidget::render( flags );
  
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  if( renderFull )
  {
    // Note that we are doing any JS calls we've made (like set x-axis range),
    //  before setting the data, so they may get squashed.
    defineJavaScript();
  }
  
  // We will render background and secondary spectra first; this way
  //  if we are loading a new foreground, and unloading the background
  //  (e.g., its a different detector type), and the new energy range
  //  is less/different than previous, we will update the viewable
  //  range correctly.  This is because:
  //    - The energy range is only adjusted by `renderForegroundToClient`
  //    - The x-range bounds are set by largest loaded spectra in JS
  //    - If we load foreground first, then background will still be
  //      loaded in JS, so the adjusted energy range will reflect
  //      background, which may go to like 12 MeV, while new foreground
  //      only goes to like 3 MeV
  
  if( m_renderFlags.testFlag(UpdateBackgroundSpectrum) )
    renderBackgroundToClient();
  
  if( m_renderFlags.testFlag(UpdateSecondarySpectrum) )
    renderSecondDataToClient();
  
  if( m_renderFlags.testFlag(UpdateForegroundSpectrum) )
    renderForegroundToClient();
  
  if( m_renderFlags.testFlag(UpdateForegroundPeaks) )
    setForegroundPeaksToClient();
  
  if( m_renderFlags.testFlag(UpdateHighlightRegions) )
    setHighlightRegionsToClient();
  
  if( m_renderFlags.testFlag(D3RenderActions::UpdateRefLines) )
    setReferenceLinesToClient();
  
  m_renderFlags = 0;
}//void render( flags )


std::string D3SpectrumDisplayDiv::localizedStringsJson() const
{
  // Strings used in the chart that we will localize, most, but probably not quite all:
  // "Real Time", "Live Time", "Dead Time", "Scaled by", "counts", "chan", "Zoom In", "Neutrons", "neutrons", "cps"
  // "ROI counts", "Gamma Counts", "cont. area", "peak cps", "peak area", "FWHM", "mean", "Spectrum"
  // "Compton Edge", "Sum Peak", "Click to set sum peak first energy.", "Single Escape", "Double Escape"
  // "Will Erase Peaks In Range", "Zoom-In on Y-axis", "Zoom-out on Y-axis"
  
  return "{\n\t"
  "xAxisTitle: " + m_xAxisTitle.jsStringLiteral('\'') + ",\n\t"                                   // "Energy (keV)"
  "yAxisTitle: " + m_yAxisTitle.jsStringLiteral('\'') + ",\n\t"                                   // "Counts"
  "yAxisTitleMulti: " + m_yAxisTitleMulti.jsStringLiteral('\'') + ",\n\t"                         // "Counts per {1} Channels"
  "title: " + (m_chartTitle.empty() ? string("null") : m_chartTitle.jsStringLiteral('\'')) + ",\n\t" //""
  "realTime: " + WString::tr("d3sdd-realTime").jsStringLiteral() + ",\n\t"                        // "Real Time"
  "liveTime: " + WString::tr("d3sdd-liveTime").jsStringLiteral() + ",\n\t"                        // "Live Time"
  "deadTime: " + WString::tr("d3sdd-deadTime").jsStringLiteral() + ",\n\t"                        // "Dead Time"
  "scaledBy: " + WString::tr("d3sdd-scaledBy").jsStringLiteral() + ",\n\t"                        // "Scaled by"
  "zoomIn: " + WString::tr("d3sdd-zoomIn").jsStringLiteral() + ",\n\t"                            // "Zoom In"
  "neutrons: " + WString::tr("d3sdd-neutrons").jsStringLiteral() + ",\n\t"                        // "neutrons"
  "Neutrons: " + WString::tr("d3sdd-Neutrons").jsStringLiteral() + ",\n\t"                        // "Neutrons"
  "cps: " + WString::tr("d3sdd-cps").jsStringLiteral() + ",\n\t"                                  // "cps"
  "roiCounts: " + WString::tr("d3sdd-roiCounts").jsStringLiteral() + ",\n\t"                      // "ROI counts"
  "gammaCounts: " + WString::tr("d3sdd-gammaCounts").jsStringLiteral() + ",\n\t"                  // "Gamma Counts"
  "contArea: " + WString::tr("d3sdd-contArea").jsStringLiteral() + ",\n\t"                        // "cont. area"
  "peakCps: " + WString::tr("d3sdd-peakCps").jsStringLiteral() + ",\n\t"                          // "peak cps"
  "peakArea: " + WString::tr("d3sdd-peakArea").jsStringLiteral() + ",\n\t"                        // "peak area"
  "fwhm: " + WString::tr("FWHM").jsStringLiteral() + ",\n\t"                                      // "FWHM"
  "mean: " + WString::tr("d3sdd-mean").jsStringLiteral() + ",\n\t"                                // "mean"
  "spectrum: " + WString::tr("d3sdd-spectrum").jsStringLiteral() + ",\n\t"                        // "Spectrum"
  "comptonEdge: " + WString::tr("d3sdd-comptonEdge").jsStringLiteral() + ",\n\t"                  // "Compton Edge"
  "sumPeak: " + WString::tr("d3sdd-sumPeak").jsStringLiteral() + ",\n\t"                          // "Sum Peak"
  "firstEnergyClick: " + WString::tr("d3sdd-firstEnergyClick").jsStringLiteral() + ",\n\t"        // "Click to set sum peak first energy."
  "singleEscape: " + WString::tr("d3sdd-singleEscape").jsStringLiteral() + ",\n\t"                // "Single Escape"
  "doubleEscape: " + WString::tr("d3sdd-doubleEscape").jsStringLiteral() + ",\n\t"                // "Double Escape"
  "eraseInRange: " + WString::tr("d3sdd-eraseInRange").jsStringLiteral() + ",\n\t"                // "Will Erase Peaks In Range"
  "zoomInY: " + WString::tr("d3sdd-zoomInY").jsStringLiteral() + ",\n\t"                          // "Zoom-In on Y-axis"  (make sure to search for "Zoom-In")
  "zoomOutY: " + WString::tr("d3sdd-zoomOutY").jsStringLiteral() + ",\n\t"                        // "Zoom-out on Y-axis" (make sure to search for "Zoom-Out")
  "touchDefineRoi: " + WString::tr("d3sdd-touchDefineRoi").jsStringLiteral() + ",\n\t"            // "Move 2 fingers to right to define ROI"
  "foreNSigmaBelowBack: " + WString::tr("d3sdd-foreNSigmaBelowBack").jsStringLiteral() + ",\n\t"  // "Foreground is {1} σ below background."
  "foreNSigmaAboveBack: " + WString::tr("d3sdd-foreNSigmaAboveBack").jsStringLiteral() + ",\n\t"  // "Foreground is {1} σ above background."
  "backSubCounts: " + WString::tr("d3sdd-backSubCounts").jsStringLiteral() + ",\n\t"              // "counts (BG sub)"
  "afterScalingBy: " + WString::tr("d3sdd-afterScalingBy").jsStringLiteral() + ",\n\t"            // "after scaling by "
  "clickedPeak: " + WString::tr("d3sdd-clickedPeak").jsStringLiteral() + ",\n\t"                  // "Clicked Peak"
  "recalFromTo: " + WString::tr("d3sdd-recalFromTo").jsStringLiteral() + ",\n\t"                  // "Recalibrate data from {1} to {2} keV"
  "sumFromTo: " + WString::tr("d3sdd-sumFromTo").jsStringLiteral() + ",\n\t"                      // "{1} to {2} keV"
  "comptonPeakAngle: " + WString::tr("d3sdd-comptonPeakAngle").jsStringLiteral() + ",\n\t"        // "{1}° Compton Peak"
  "}";
  
  // In JS, probably need to call: `self.handleResize( false );`
}//std::string localizedStringsJson() const


void D3SpectrumDisplayDiv::refresh()
{
  // If we are here, the language has probably been changed - we'll update the localization
  //  phrases the JS uses.
  WContainerWidget::refresh();

  if( isRendered() )
    doJavaScript( m_jsgraph + ".setLocalizations( " + localizedStringsJson() + ",\n  false);" );
  
  m_renderFlags |= UpdateHighlightRegions;
}//void refresh();


double D3SpectrumDisplayDiv::xAxisMinimum() const
{
  return m_xAxisMinimum;
}//double xAxisMinimum() const


double D3SpectrumDisplayDiv::xAxisMaximum() const
{
  return m_xAxisMaximum;
}//double xAxisMaximum() const


double D3SpectrumDisplayDiv::chartWidthInPixels() const
{
  return m_chartWidthPx;
}//double chartWidthInPixels() const

double D3SpectrumDisplayDiv::chartHeightInPixels() const
{
  return m_chartHeightPx;
}//double chartHeightInPixels() const

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
  if( m_yAxisIsLog == log )
    return;
  
  m_yAxisIsLog = log;
  if( isRendered() )
    doJavaScript( m_jsgraph + (log ? ".setLogY();" : ".setLinearY();") );
  
  m_yAxisLogLinChanged.emit( log );
}//void setYAxisLog( bool log )

void D3SpectrumDisplayDiv::showGridLines( bool show )
{
  m_showVerticalLines = show;
  m_showHorizontalLines = show;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridX(" + jsbool(show) + ");"
                  + m_jsgraph + ".setGridY(" + jsbool(show) + ");" );
}

void D3SpectrumDisplayDiv::showVerticalLines( const bool draw )
{
  m_showVerticalLines = draw;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridX(" + jsbool(draw) + ");" );
}

void D3SpectrumDisplayDiv::showHorizontalLines( const bool draw )
{
  m_showHorizontalLines = draw;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridY(" + jsbool(draw) + ");" );
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
  return m_backgroundSubtract;
}//bool backgroundSubtract() const


void D3SpectrumDisplayDiv::setBackgroundSubtract( bool subtract )
{
  if( subtract == m_backgroundSubtract )
    return;
  
  m_backgroundSubtract = subtract;
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setBackgroundSubtract(" + jsbool(subtract) + ");" );
}//void setBackgroundSubtract( bool subtract )

void D3SpectrumDisplayDiv::setXAxisMinimum( const double minimum )
{
  const string minimumStr = to_string( minimum );
  m_xAxisMinimum = minimum;
  
  string js = m_jsgraph + ".setXAxisMinimum(" + minimumStr + ");";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setXAxisMinimum( const double minimum )


void D3SpectrumDisplayDiv::setXAxisMaximum( const double maximum )
{
  const string maximumStr = to_string( maximum );
  m_xAxisMaximum = maximum;
  
  string js = m_jsgraph + ".setXAxisMaximum(" + maximumStr + ");";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setXAxisMaximum( const double maximum )


void D3SpectrumDisplayDiv::setYAxisMinimum( const double minimum )
{
  const string minimumStr = to_string( minimum );
  m_yAxisMinimum = minimum;
  
  string js = m_jsgraph + ".setYAxisMinimum(" + minimumStr + ");";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setYAxisMinimum( const double minimum )


void D3SpectrumDisplayDiv::setYAxisMaximum( const double maximum )
{
  const string maximumStr = to_string( maximum );
  m_yAxisMaximum = maximum;
  
  string js = m_jsgraph + ".setYAxisMaximum(" + maximumStr + ");";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setYAxisMaximum( const double maximum )


void D3SpectrumDisplayDiv::setXAxisRange( const double minimum, const double maximum )
{
  const string minimumStr = to_string( minimum );
  const string maximumStr = to_string( maximum );
  m_xAxisMinimum = minimum;
  m_xAxisMaximum = maximum;
  
  string js = m_jsgraph + ".setXAxisRange(" + minimumStr + "," + maximumStr + ",false);"
              + m_jsgraph + ".redraw()();";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setXAxisRange( const double minimum, const double maximum );


std::tuple<double,double,Wt::WString> D3SpectrumDisplayDiv::setYAxisRange( double lower_counts,
                                       double upper_counts )
{
  if( upper_counts < lower_counts )
    std::swap( lower_counts, upper_counts );
  
  if( upper_counts == lower_counts )
    return {m_yAxisMinimum, m_yAxisMaximum, WString::tr("d3sdd-yaxis-lower-upper-equal")};
  
  WString errmsg;
  const bool isLogY = yAxisIsLog();
  if( isLogY )
  {
    const auto hist = m_foreground;
    if( !hist || !hist->gamma_counts() || (hist->gamma_counts()->size() < 2) )
      return {m_yAxisMinimum, m_yAxisMaximum, WString::tr("d3sdd-yaxis-no-spectrum")};
    
    const shared_ptr<const vector<float>> &channels = hist->gamma_counts();
    assert( channels );
    
    // Lets check the y-range - ignoring we may be showing multiple channels per display-bin,
    // to see if the requested range is reasonable
    double xmin, xmax, ymin, ymax;
    visibleRange( xmin, xmax, ymin, ymax );
    
    const size_t lower_channel = hist->find_gamma_channel(xmin);
    const size_t upper_channel = hist->find_gamma_channel(xmax);
    float min_nonzero_value = std::numeric_limits<float>::max();
    float max_value = -std::numeric_limits<float>::max();
    
    for( size_t channel = lower_channel; channel <= upper_channel; ++channel )
    {
      if( channel < channels->size() )
      {
        const float val = (*channels)[channel];
        max_value = std::max( max_value, val );
        if( val > 0.0f )
          min_nonzero_value = std::min( min_nonzero_value, val );
      }//if( channel < channels->size() )
    }//for( loop over visible channels )
    
    if( (max_value <= 0.0f) || IsInf(max_value) || IsNan(max_value) )
      max_value = 1.0f;
    
    if( (min_nonzero_value == std::numeric_limits<float>::max())
       || IsInf(min_nonzero_value)
       || IsNan(min_nonzero_value) )
    {
      min_nonzero_value = 0.1f*max_value;
    }
    
    if( upper_counts <= 0.0f )
    {
      errmsg = WString::tr("d3sdd-yaxis-upper-counts-below-zero");
      upper_counts = (max_value > 0.0f) ? 2.0f*max_value : 1.0f;
    }
    
    if( upper_counts <= min_nonzero_value )
    {
      errmsg = WString::tr("d3sdd-yaxis-below-min-non-zero");
      upper_counts = 2.0f*min_nonzero_value;
    }
    
    if( (lower_counts <= 0.0) || (lower_counts > max_value) )
    {
      errmsg = WString::tr("d3sdd-yaxis-lower-counts-below-zero");
      lower_counts = min_nonzero_value;
    }
  }//if( isLogY && (lower_counts <= 0.0f) )
  
  
  if( isLogY )
  {
    const double prevMinLogY = logYAxisMin();
    if( lower_counts < prevMinLogY )
      setLogYAxisMin( lower_counts );
  }//if( isLogY )
  
  
  const string minimumStr = SpecUtils::printCompact(lower_counts, 8);
  const string maximumStr = SpecUtils::printCompact(upper_counts, 8);
  m_yAxisMinimum = lower_counts;
  m_yAxisMaximum = upper_counts;
  
  string js = m_jsgraph + ".setYAxisRange(" + minimumStr + "," + maximumStr + ");";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
  
  return {m_yAxisMinimum, m_yAxisMaximum, errmsg};
}//void setYAxisRange( const double minimum, const double maximum );


void D3SpectrumDisplayDiv::doBackgroundLiveTimeNormalization()
{
  if( !m_background || !m_foreground
     || (m_background->live_time() <= 0.0f)
     || (m_foreground->live_time() <= 0.0f) )
  {
    m_backgroundScale = 1.0f;
  }else
  {
    m_backgroundScale = m_foreground->live_time() / m_background->live_time();
  }
}//void doBackgroundLiveTimeNormalization()


void D3SpectrumDisplayDiv::doSecondaryLiveTimeNormalization()
{
  if( !m_secondary || !m_foreground
     || (m_secondary->live_time() <= 0.0f)
     || (m_foreground->live_time() <= 0.0f) )
  {
    m_secondaryScale = 1.0f;
  }else
  {
    m_secondaryScale = m_foreground->live_time() / m_secondary->live_time();
  }
}//void doSecondaryLiveTimeNormalization()


void D3SpectrumDisplayDiv::setForegroundPeaksToClient()
{
  string js;
  
  if( m_peakModel )
  {
    std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
    if( peaks )
    {
      vector< std::shared_ptr<const PeakDef> > inpeaks( peaks->begin(), peaks->end() );
      const std::shared_ptr<const Measurement> &foreground = m_foreground;
      js = PeakDef::peak_json( inpeaks, foreground );
    }
  }
  
  if( js.empty() )
    js = "[]";
  
  js = m_jsgraph + ".setRoiData(" + js + ", 'FOREGROUND');";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setForegroundPeaksToClient();

void D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw()
{
  m_renderFlags |= UpdateForegroundPeaks;
  
  scheduleRender();
}//void scheduleForegroundPeakRedraw()


void D3SpectrumDisplayDiv::setData( std::shared_ptr<const Measurement> data_hist, const bool keep_curent_xrange )
{
  const float oldBackSF = m_backgroundScale;
  const float oldSecondSF = m_secondaryScale;
  
  m_foreground = data_hist;
  if( m_peakModel )
    m_peakModel->setForeground( data_hist );
  
  if( !keep_curent_xrange )
    m_renderFlags |= ResetXDomain;
  
  scheduleUpdateForeground();
  
  // Update background and second SF
  doBackgroundLiveTimeNormalization();
  doSecondaryLiveTimeNormalization();
  
  
  //If the background or secondary data spectrum scale factors changed, we need
  //  to update those to the client as well.
  if( m_background && (fabs(m_backgroundScale - oldBackSF) > 0.000001f*std::max(m_backgroundScale,oldBackSF)) )
    scheduleUpdateBackground();
  
  if( m_secondary && (fabs(m_secondaryScale - oldSecondSF) > 0.000001f*std::max(m_secondaryScale,oldSecondSF)) )
    scheduleUpdateSecondData();
}//void setData( std::shared_ptr<const Measurement> data_hist )


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::data() const
{
  return m_foreground;
}//std::shared_ptr<const Measurement> data() const


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::secondData() const
{
  return m_secondary;
}//std::shared_ptr<const Measurement> secondData() const


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::background() const
{
  return m_background;
}//std::shared_ptr<const Measurement> background() const


float D3SpectrumDisplayDiv::foregroundLiveTime() const
{
  return m_foreground ? m_foreground->live_time() : 0.0f;
}


float D3SpectrumDisplayDiv::foregroundRealTime() const
{
  return m_foreground ? m_foreground->real_time() : 0.0f;
}


float D3SpectrumDisplayDiv::backgroundLiveTime() const
{
  return m_background ? m_background->live_time() : 0.0f;
}


float D3SpectrumDisplayDiv::backgroundRealTime() const
{
  return m_background ? m_background->real_time() : 0.0f;
}


float D3SpectrumDisplayDiv::secondForegroundLiveTime() const
{
  return m_secondary ? m_secondary->live_time() : 0.0f;
}


float D3SpectrumDisplayDiv::secondForegroundRealTime() const
{
  return m_secondary ? m_secondary->real_time() : 0.0f;
}


std::shared_ptr<const Measurement> D3SpectrumDisplayDiv::histUsedForXAxis() const
{
  if( m_foreground )
    return m_foreground;
  return nullptr;
}


void D3SpectrumDisplayDiv::setDisplayScaleFactor( const float sf,
                                               const SpecUtils::SpectrumType spectrum_type )
{
  switch( spectrum_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      throw runtime_error( "setDisplayScaleFactor can not be called for foreground" );
      
    case SpecUtils::SpectrumType::SecondForeground:
      m_secondaryScale = sf;
      scheduleUpdateSecondData();
      break;
      
    case SpecUtils::SpectrumType::Background:
      m_backgroundScale = sf;
      scheduleUpdateBackground();
      break;
  }//switch( spectrum_type )
}//void setDisplayScaleFactor(...)


float D3SpectrumDisplayDiv::displayScaleFactor( const SpecUtils::SpectrumType spectrum_type ) const
{
  switch( spectrum_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      return 1.0f;
    case SpecUtils::SpectrumType::SecondForeground:
      return m_secondaryScale;
    case SpecUtils::SpectrumType::Background:
      return m_backgroundScale;
  }//switch( spectrum_type )
  
  throw runtime_error( "D3SpectrumDisplayDiv::displayScaleFactor(...): invalid input arg" );
  
  return 1.0;
}//double displayScaleFactor( SpecUtils::SpectrumType spectrum_type ) const;


void D3SpectrumDisplayDiv::setBackground( std::shared_ptr<const Measurement> background )
{
  m_background = background;
  doBackgroundLiveTimeNormalization();
  if( !background && m_backgroundSubtract )
    setBackgroundSubtract( false );
  scheduleUpdateBackground();
}//void D3SpectrumDisplayDiv::setBackground(...);


void D3SpectrumDisplayDiv::setSecondData( std::shared_ptr<const Measurement> hist )
{
  m_secondary = hist;
  doSecondaryLiveTimeNormalization();
  scheduleUpdateSecondData();
}//void D3SpectrumDisplayDiv::setSecondData( std::shared_ptr<const Measurement> background );


void D3SpectrumDisplayDiv::visibleRange( double &xmin, double &xmax,
                                      double &ymin, double &ymax ) const
{
  xmin = m_xAxisMinimum;
  xmax = m_xAxisMaximum;
  ymin = m_yAxisMinimum;
  ymax = m_yAxisMaximum;
}


const Wt::WString &D3SpectrumDisplayDiv::xAxisTitle() const
{
  return m_xAxisTitle;
}//const Wt::WString &xAxisTitle() const;


const Wt::WString &D3SpectrumDisplayDiv::yAxisTitle() const
{
  return m_yAxisTitle;
}//const Wt::WString &yAxisTitle() const;


void D3SpectrumDisplayDiv::setXAxisTitle( const Wt::WString &title )
{
  m_xAxisTitle = title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle(" + m_xAxisTitle.jsStringLiteral('\'') + "');" );
}//void setXAxisTitle( const std::string &title )


void D3SpectrumDisplayDiv::setYAxisTitle( const Wt::WString &title, const Wt::WString &titleMulti )
{
  m_yAxisTitle = title;
  m_yAxisTitleMulti = titleMulti;
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setYAxisTitle(" + m_yAxisTitle.jsStringLiteral('\'')
                 + ", " + m_yAxisTitleMulti.jsStringLiteral('\'') + ");" );
}//void setYAxisTitle( const std::string &title )



Wt::Signal<double,double,double,double,bool> &D3SpectrumDisplayDiv::xRangeChanged()
{
  return m_xRangeChanged;
}//xRangeChanged()


Wt::Signal<double,double> &D3SpectrumDisplayDiv::shiftAltKeyDragged()
{
  return m_shiftAltKeyDragg;
}


void D3SpectrumDisplayDiv::setChartTitle( const Wt::WString &title )
{
  m_chartTitle = title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setTitle(" + title.jsStringLiteral('\'') + ");" );
}//void setChartTitle( const Wt::WString &title )


void D3SpectrumDisplayDiv::setSearchEnergies( const vector<pair<double,double>> &energies )
{
  m_searchEnergies = energies;
  
  string js = "[";
  for( size_t i = 0; i < energies.size(); ++i )
  {
    js += string(i ? "," : "") + "{energy:" + std::to_string(energies[i].first)
          + ", window:" + std::to_string(energies[i].second) + "}";
  }
  js += "]";

  if( isRendered() )
    doJavaScript( m_jsgraph + ".setSearchWindows(" + js + ");" );
}//void setSearchEnergies( const std::vector<std::pair<double,double>> &energies )


bool D3SpectrumDisplayDiv::removeDecorativeHighlightRegion( size_t uniqueid )
{
  if( uniqueid < 3 )
    return false;
  
  const size_t nprev = m_highlights.size();
  for( size_t i = 0; i < nprev; ++i )
  {
    if( m_highlights[i].hash == uniqueid )
    {
      m_highlights.erase( m_highlights.begin() + i );
      
      m_renderFlags |= UpdateHighlightRegions;
      scheduleRender();
      return true;
    }
  }
  
#if( PERFORM_DEVELOPER_CHECKS )
  log_developer_error( __func__, ("D3SpectrumDisplayDiv::removeDecorativeHighlightRegion:"
                                  " Didnt find uniqueid=" + std::to_string(uniqueid)).c_str() );
#endif
  
  return false;
}//void removeDecorativeHighlightRegions()


void D3SpectrumDisplayDiv::removeAllDecorativeHighlightRegions()
{
  if( m_highlights.empty() )
    return;
  
  m_highlights.clear();
  m_renderFlags |= UpdateHighlightRegions;
  scheduleRender();
}//void removeAllDecorativeHighlightRegions()


size_t D3SpectrumDisplayDiv::addDecorativeHighlightRegion( const float lowerx,
                                                          const float upperx,
                                                          const Wt::WColor &color,
                                                          const HighlightRegionFill fillType,
                                                          const Wt::WString &txt  )
{
  D3SpectrumDisplayDiv::HighlightRegion region;
  region.lowerx = lowerx;
  region.upperx = upperx;
  region.color = color;
  region.fill_type = fillType;
  region.display_txt = txt;
  region.hash = 0;
  
  boost::hash_combine( region.hash, lowerx );
  boost::hash_combine( region.hash, upperx );
  boost::hash_combine( region.hash, color.red() );
  boost::hash_combine( region.hash, color.green() );
  boost::hash_combine( region.hash, color.blue() );
  boost::hash_combine( region.hash, color.alpha() );
  boost::hash_combine( region.hash, static_cast<int>(fillType) );
  boost::hash_combine( region.hash, hash<string>()(txt.toUTF8()) );
  
  if( region.hash <= 2 )
    region.hash += 3;
  
  //should in principle check for collision, but whatever
  
  m_highlights.push_back( region );
  
  m_renderFlags |= UpdateHighlightRegions;
  scheduleRender();
  
  return region.hash;
}//void addDecorativeHighlightRegion(...)


void D3SpectrumDisplayDiv::setHighlightRegionsToClient()
{
  WStringStream jsstrm;
  jsstrm << m_jsgraph << ".setHighlightRegions";
  
  if( m_highlights.empty() )
  {
    jsstrm << "(null);";
  }else
  {
    jsstrm << "([";
    for( size_t i = 0; i < m_highlights.size(); ++i )
    {
      const D3SpectrumDisplayDiv::HighlightRegion &region = m_highlights[i];
      jsstrm << std::string(i ? ",{" : "{")
             << "lowerEnergy:" << static_cast<double>(region.lowerx)
             << ", upperEnergy:" << static_cast<double>(region.upperx)
             << ", fill: 'rgba(" << region.color.red() << "," << region.color.green()  //bug in Wt 3.3.4 means need to manually write out rgba()
                                 << "," << region.color.blue()
             << "," << (region.color.alpha()/255.0) << ")'";
      
      switch( region.fill_type )
      {
        case HighlightRegionFill::Full:
          // We'll skip putting this in - previous to 20240527 this field didnt exist, so we'll keep compatiblity
          //jsstrm << ", drawRegion: 'All'";
          break;
          
        case HighlightRegionFill::BelowData:
          jsstrm << ", drawRegion: 'BelowData'";
          break;
      }//switch( region.fill_type )
      
      if( !region.display_txt.empty() )
        jsstrm << ", text: " << region.display_txt.jsStringLiteral('\'');
      
      jsstrm << ", hash:" << std::to_string(region.hash)
             << "}";
    }//for( loop over m_highlights )
  
    jsstrm << "]);";
  }//if( no highlights ) / else
  
  if( isRendered() )
    doJavaScript( jsstrm.str() );
}//void setHighlightRegionsToClient()


void D3SpectrumDisplayDiv::scheduleUpdateForeground()
{
  m_renderFlags |= UpdateForegroundSpectrum;
  scheduleRender();
}


void D3SpectrumDisplayDiv::scheduleUpdateBackground()
{
  m_renderFlags |= UpdateBackgroundSpectrum;
  scheduleRender();
}


void D3SpectrumDisplayDiv::scheduleUpdateSecondData()
{
  m_renderFlags |= UpdateSecondarySpectrum;
  scheduleRender();
}


void D3SpectrumDisplayDiv::renderForegroundToClient()
{
  const std::shared_ptr<const Measurement> &data_hist = m_foreground;
  
  string js;
  
  const string resetDomain = m_renderFlags.testFlag(ResetXDomain) ? "true" : "false";
  
  // Set the data for the chart
  if ( data_hist )
  {
    // Create the measurement array (should only have one measurement)
    std::ostringstream ostr;
    std::vector< std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> foregroundData;
    D3SpectrumExport::D3SpectrumOptions foregroundOptions;
    
    // Set options for the spectrum
    foregroundOptions.peak_color = m_defaultPeakColor.isDefault() ? string("blue") : m_defaultPeakColor.cssText();
    foregroundOptions.spectrum_type = SpecUtils::SpectrumType::Foreground;
    foregroundOptions.display_scale_factor = displayScaleFactor( SpecUtils::SpectrumType::Foreground );  //will always be 1.0f
    // Leave line color blank, as we set a CSS rule for `--d3spec-fore-line-color`
    //foregroundOptions.line_color = m_foregroundLineColor.isDefault() ? string("black") : m_foregroundLineColor.cssText();
    
    // Set the peak data for the spectrum
    if ( m_peakModel )
    {
      std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
      vector< std::shared_ptr<const PeakDef> > inpeaks;
      if( peaks )
        inpeaks.insert( end(inpeaks), peaks->begin(), peaks->end() );
      
      foregroundOptions.peaks_json = PeakDef::peak_json( inpeaks, data_hist );
    }
    
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(data_hist.get(),foregroundOptions) );
    
    // Set the data on the JS side
    if ( D3SpectrumExport::write_and_set_data_for_chart(ostr, id(), measurements) ) {
      string data = ostr.str();
      size_t index = data.find( "spec_chart_" );
      data = data.substr( 0, index );
      js = data + m_jsgraph + ".setSpectrumData(data_" + id() + ", " + resetDomain + ", 'FOREGROUND', 0, 1 );";
    }
  } else {
    //js = m_jsgraph + ".removeSpectrumDataByType(" + resetDomain + ", 'FOREGROUND' );";
    js = m_jsgraph + ".setData(null,true);";
  }//if ( data_hist ) / else
  
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void D3SpectrumDisplayDiv::updateData()


void D3SpectrumDisplayDiv::renderBackgroundToClient()
{
  string js;
  const std::shared_ptr<const Measurement> &background = m_background;
  
  // Set the data for the chart
  if ( background ) {
    // Create the measurement array (should only have one measurement)
    std::ostringstream ostr;
    std::vector< std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> backgroundData;
    D3SpectrumExport::D3SpectrumOptions backgroundOptions;
    
    // Set options for the spectrum
    backgroundOptions.spectrum_type = SpecUtils::SpectrumType::Background;
    backgroundOptions.display_scale_factor = displayScaleFactor( SpecUtils::SpectrumType::Background );
    // Leave line color blank, as we set a CSS rule for `--d3spec-back-line-color`
    //backgroundOptions.line_color = m_backgroundLineColor.isDefault() ? string("green") : m_backgroundLineColor.cssText();
    
    // We cant currently access the parent InterSpec class, but if we could, then
    //  we could draw the background peaks by doing something like:
    //const std::set<int> &backSample = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
    //std::shared_ptr<SpecMeas> backgroundMeas = m_interspec->measurment(SpecUtils::SpectrumType::Background);
    //std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > backpeaks = backgroundMeas->peaks( backSample );
    //vector< std::shared_ptr<const PeakDef> > inpeaks( backpeaks->begin(), backpeaks->end() );
    //backgroundOptions.peaks_json = PeakDef::peak_json( inpeaks );
    
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(background.get(), backgroundOptions) );
    
    // Set the data on the JS side
    if ( D3SpectrumExport::write_and_set_data_for_chart(ostr, id(), measurements) ) {
      string data = ostr.str();
      size_t index = data.find( "spec_chart_" );
      data = data.substr( 0, index );
      js = data + m_jsgraph + ".setSpectrumData(data_" + id() + ", false, 'BACKGROUND', 1, -1);";
    }
  } else {
    js = m_jsgraph + ".removeSpectrumDataByType(false, 'BACKGROUND' );";
  }//if ( background )
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void D3SpectrumDisplayDiv::updateBackground()


void D3SpectrumDisplayDiv::renderSecondDataToClient()
{
  string js;
  const std::shared_ptr<const Measurement> &hist = m_secondary;
  
  // Set the data for the chart
  if ( hist ) {
    // Create the measurement array (should only have one measurement)
    std::ostringstream ostr;
    std::vector< std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    std::pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> secondaryData;
    D3SpectrumExport::D3SpectrumOptions secondaryOptions;
    
    // Set options for the spectrum
    secondaryOptions.spectrum_type = SpecUtils::SpectrumType::SecondForeground;
    secondaryOptions.display_scale_factor = displayScaleFactor( SpecUtils::SpectrumType::SecondForeground );
    // Leave line color blank, as we set a CSS rule for `--d3spec-second-line-color`
    //secondaryOptions.line_color = m_secondaryLineColor.isDefault() ? string("steelblue") : m_secondaryLineColor.cssText();
    
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(hist.get(), secondaryOptions) );
    
    // Set the data on the JS side
    if ( D3SpectrumExport::write_and_set_data_for_chart(ostr, id(), measurements) ) {
      string data = ostr.str();
      size_t index = data.find( "spec_chart_" );
      data = data.substr( 0, index );
      js = data + m_jsgraph + ".setSpectrumData(data_" + id() + ", false, 'SECONDARY', 2, 1);";
    }
  } else {
    js = m_jsgraph + ".removeSpectrumDataByType(false, 'SECONDARY' );";
  }//if ( hist )
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void D3SpectrumDisplayDiv::updateSecondData()


void D3SpectrumDisplayDiv::applyColorTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( theme )
  {
    setForegroundSpectrumColor( theme->foregroundLine );
    setBackgroundSpectrumColor( theme->backgroundLine );
    setSecondarySpectrumColor( theme->secondaryLine );
    setDefaultPeakColor( theme->defaultPeakLine );
    
    setAxisLineColor( theme->spectrumAxisLines );
    setChartMarginColor( theme->spectrumChartMargins );
    setChartBackgroundColor( theme->spectrumChartBackground );
    setTextColor( theme->spectrumChartText );
    
    setPeakLabelSize( theme->spectrumPeakLabelSize.toUTF8() );
    setPeakLabelRotation( theme->spectrumPeakLabelRotation );
    setLogYAxisMin( theme->spectrumLogYAxisMin );
  }else
  {
    setForegroundSpectrumColor( {} );
    setBackgroundSpectrumColor( {} );
    setSecondarySpectrumColor( {} );
    setDefaultPeakColor( {} );
    
    setAxisLineColor( {} );
    setChartMarginColor( {} );
    setChartBackgroundColor( {} );
    setTextColor( {} );
    
    setPeakLabelSize( "" );
    setPeakLabelRotation( 0.0 );
    setLogYAxisMin( 0.1 );
  }
}//void applyColorTheme( std::shared_ptr<const ColorTheme> theme );


void D3SpectrumDisplayDiv::setForegroundSpectrumColor( const Wt::WColor &color )
{
  m_foregroundLineColor = color.isDefault() ? WColor( 0x00, 0x00, 0x00 ) : color;
  
  string rulename = "ForeSpecLineColor";
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-fore-line-color: " + m_foregroundLineColor.cssText() );
  
  //scheduleUpdateForeground();
}

void D3SpectrumDisplayDiv::setBackgroundSpectrumColor( const Wt::WColor &color )
{
  m_backgroundLineColor = color.isDefault() ? WColor(0x00,0xff,0xff) : color;
  
  string rulename = "BackSpecLineColor";
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-back-line-color: " + m_backgroundLineColor.cssText() );
  
  //scheduleUpdateBackground();
}

void D3SpectrumDisplayDiv::setSecondarySpectrumColor( const Wt::WColor &color )
{
  m_secondaryLineColor = color.isDefault() ? WColor(0x00,0x80,0x80) : color;
  
  string rulename = "SecondSpecLineColor";
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-second-line-color: " + m_secondaryLineColor.cssText() );
  
  //scheduleUpdateSecondData();
}

void D3SpectrumDisplayDiv::setTextColor( const Wt::WColor &color )
{
  m_textColor = color.isDefault() ? WColor(0,0,0) : color;
  const string c = m_textColor.cssText();
  const string rulename = "TextColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-text-color: " + c + ";" );
}


void D3SpectrumDisplayDiv::setAxisLineColor( const Wt::WColor &color )
{
  m_axisColor = color.isDefault() ? WColor(0,0,0) : color;
  
  string rulename = "AxisColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  // Sets axis colors, and ".peakLine, .escapeLineForward, .mouseLine, .secondaryMouseLine"
  
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-axis-color: " + m_axisColor.cssText() );
  
  
  //ToDo: figure out how to make grid lines look okay.
  //rulename = "GridColor";
  //if( m_cssRules.count(rulename) )
  //  style.removeRule( m_cssRules[rulename] );
  //m_cssRules[rulename] = style.addRule( ".xgrid > .tick, .ygrid > .tick", "stroke: #b3b3b3" );
  //
  //rulename = "MinorGridColor";
  //if( m_cssRules.count(rulename) )
  //  style.removeRule( m_cssRules[rulename] );
  //m_cssRules[rulename] = style.addRule( ".minorgrid", "stroke: #e6e6e6" );
}

void D3SpectrumDisplayDiv::setChartMarginColor( const Wt::WColor &color )
{
  m_chartMarginColor = color;
  
  const string rulename = "MarginColor";
  
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
  
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-background-color: " + color.cssText() + ";" );
  //Actually this will set the background for the entire chart...
  //m_cssRules[rulename] = style.addRule( "#" + id() + " > svg", "background: " + color.cssText() ); // to set background color for just this chart
}//setChartMarginColor(...)


void D3SpectrumDisplayDiv::setChartBackgroundColor( const Wt::WColor &color )
{
  m_chartBackgroundColor = color;
  const string c = color.isDefault() ? "rgba(0,0,0,0)" : color.cssText();
  
  const string rulename = "BackgroundColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
    
  m_cssRules[rulename] = style.addRule( ":root", "--d3spec-chart-area-color: " + c + ";" );
  //m_cssRules[rulename] = style.addRule( "#chartarea" + id(), "fill: " + c ); //If we wanted to apply this color to only this chart
}

void D3SpectrumDisplayDiv::setDefaultPeakColor( const Wt::WColor &color )
{
  m_defaultPeakColor = color.isDefault() ? WColor(0,51,255,155) : color;
  
  //The default peak color is specified as part of the foreground JSON, so we
  //  need to reload the foreground to client to update the color.
  scheduleUpdateForeground();
}


/** Sets the charts font size.
 
 Example values: "8px", "smaller", "12", "10px", "x-small", etc.
 An empty string sets to default chart font size.
 */
void D3SpectrumDisplayDiv::setPeakLabelSize( const std::string &fontSize )
{
  m_peakLabelFontSize = fontSize;
  
  if( isRendered() )
  {
    string str = fontSize.empty() ? string("null") : Wt::WWebWidget::jsStringLiteral(fontSize,'\'');
    doJavaScript( m_jsgraph + ".setPeakLabelSize(" + str + ");" );
  }
  
  scheduleUpdateForeground(); //JIC, the JS setPeakLabelSize(...) wont cause a re-draw
}//void setPeakLabelSize( const std::string &fontSize );


void D3SpectrumDisplayDiv::setPeakLabelRotation( const double rotation )
{
  m_peakLabelRotationDegrees = rotation;
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setPeakLabelRotation(" + std::to_string(rotation) + ");" );
  
  scheduleUpdateForeground(); //JIC, the JS setPeakLabelRotation(...) wont cause a re-draw
}//void setPeakLabelRotation( const double rotation )


void D3SpectrumDisplayDiv::setLogYAxisMin( const double ymin )
{
  assert( ymin >= 0 );
  
  if( ymin <= 0.0 )
    throw runtime_error( "Log Y-Axis min must be larger than zero." );
  
  m_logYAxisMin = ymin;
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setLogYAxisMin(" + std::to_string(ymin) + "); "
                 + m_jsgraph + ".redraw()();" );
}//void setLogYAxisMin( const double ymin )


double D3SpectrumDisplayDiv::logYAxisMin() const
{
  return m_logYAxisMin;
}

void D3SpectrumDisplayDiv::saveChartToImg( const std::string &filename, const bool asPng )
{
  // For now we grab CSS rules dynamically in the JS, which is tedious and error-prone, but we could
  //  use the following to get (some of) the rules from C++ land.
  //string cssrules;
  //for( const auto &p : m_cssRules )
  //  cssrules += p.second->selector() + "{ " + p.second->declarations() + " }  ";
  //cout << "CSS rules: " << cssrules << endl;
  
  LOAD_JAVASCRIPT(wApp, "src/D3SpectrumDisplayDiv.cpp", "D3SpectrumDisplayDiv", wtjsSvgToImgDownload);
  
  if( isRendered() )
    doJavaScript( "Wt.WT.SvgToImgDownload(" + m_jsgraph + ",'" + filename + "', "
                  + string(asPng ? "true" : "false") + ");" );
}//void saveChartToPng( const std::string &name )



void D3SpectrumDisplayDiv::setFeatureMarkerOption( FeatureMarkerType option, bool show )
{
  string js = m_jsgraph;
  
  switch( option )
  {
    case FeatureMarkerType::EscapePeakMarker:
      js += ".setEscapePeaks";
      break;
    case FeatureMarkerType::ComptonPeakMarker:
      js += ".setComptonPeaks";
      break;
    case FeatureMarkerType::ComptonEdgeMarker:
      js += ".setComptonEdge";
      break;
    case FeatureMarkerType::SumPeakMarker:
      js += ".setSumPeaks";
      break;
    case FeatureMarkerType::NumFeatureMarkers:
    default:
      return;
  }//switch( option )
  
  m_showFeatureMarker[static_cast<int>(option)] = show;
  js += "(" + jsbool(show) + ");";
  
  if( isRendered() )
    doJavaScript( js );
}//void D3SpectrumDisplayDiv::setFeatureMarkerOption(...)


void D3SpectrumDisplayDiv::setComptonPeakAngle( int angle )
{
  m_comptonPeakAngle = angle;
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setComptonPeakAngle(" + std::to_string(angle) + ");" );
}//void D3SpectrumDisplayDiv::setComptonPeakAngle( int angle )


int D3SpectrumDisplayDiv::comptonPeakAngle() const
{
  return m_comptonPeakAngle;
}//int comptonPeakAngle() const


void D3SpectrumDisplayDiv::showXAxisSliderChart( const bool show )
{
  if( m_showXAxisSliderChart == show )
    return;
  
  m_showXAxisSliderChart = show;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setShowXAxisSliderChart(" + jsbool(show) + ");" );
  
  m_xAxisSliderShown.emit( show );
}//void showXAxisSliderChart( const bool show )


bool D3SpectrumDisplayDiv::xAxisSliderChartIsVisible() const
{
  return m_showXAxisSliderChart;
}//void xAxisSliderChartIsVisible() const;



void D3SpectrumDisplayDiv::showYAxisScalers( const bool show )
{
  if( m_showYAxisScalers == show )
    return;
  
  m_showYAxisScalers = show;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setShowSpectrumScaleFactorWidget(" + jsbool(show) + ");" );
}//void showXAxisSliderChart( const bool show )


bool D3SpectrumDisplayDiv::yAxisScalersIsVisible() const
{
  return m_showYAxisScalers;
}//void xAxisSliderChartIsVisible() const;


void D3SpectrumDisplayDiv::chartShiftKeyDragCallback( double x0, double x1 )
{
  // cout << "chartShiftKeyDragCallback" << endl;
  m_shiftKeyDragg.emit( x0, x1 );
}//void D3SpectrumDisplayDiv::chartShiftKeyDragCallback(...)

void D3SpectrumDisplayDiv::chartShiftAltKeyDragCallback( double x0, double x1 )
{
  // cout << "chartShiftAltKeyDragCallback" << endl;
  m_shiftAltKeyDragg.emit( x0, x1 );
}//void D3SpectrumDisplayDiv::chartShiftAltKeyDragCallback(...)

void D3SpectrumDisplayDiv::chartRightMouseDragCallback( double x0, double x1 )
{
  // cout << "chartRightMouseDragCallback" << endl;
  m_rightMouseDragg.emit( x0, x1 );
}//void D3SpectrumDisplayDiv::chartRightMouseDragCallback(...)

void D3SpectrumDisplayDiv::chartLeftClickCallback( double x, double y, double pageX, double pageY )
{
  // cout << "chartLeftClickCallback" << endl;
  m_leftClick.emit( x, y, pageX, pageY );
}//void D3SpectrumDisplayDiv::chartDoubleLeftClickCallback(...)

void D3SpectrumDisplayDiv::chartDoubleLeftClickCallback( double x, double y )
{
  // cout << "chartDoubleLeftClickCallback" << endl;
  m_doubleLeftClick.emit( x, y );
}//void D3SpectrumDisplayDiv::chartDoubleLeftClickCallback(...)

void D3SpectrumDisplayDiv::chartRightClickCallback( double x, double y, double pageX, double pageY )
{
  // cout << "chartRightClickCallback" << endl;
  m_rightClick.emit( x, y, pageX, pageY );
}//void D3SpectrumDisplayDiv::chartRightClickCallback(...)


void D3SpectrumDisplayDiv::existingRoiEdgeDragCallback( double new_lower_energy, double new_upper_energy,
                                                  double new_lower_px, double new_upper_px,
                                                  double original_lower_energy, bool isfinal )
{
  m_existingRoiEdgeDrag.emit( new_lower_energy, new_upper_energy, new_lower_px, new_upper_px,
                                             original_lower_energy, isfinal );
}//void D3SpectrumDisplayDiv::chartRoiDragedCallback(...)


void D3SpectrumDisplayDiv::performExistingRoiEdgeDragWork(
                                                  double new_lower_energy, double new_upper_energy,
                                                  double new_lower_px, double new_upper_px,
                                                  double original_lower_energy, bool isfinal )
{
  std::unique_ptr<UndoRedoManager::PeakModelChange> peak_undo_creator;
  if( isfinal )
    peak_undo_creator.reset( new UndoRedoManager::PeakModelChange() );
  
//  cout << "chartRoiDragedCallback: energy={" << new_lower_energy << "," << new_upper_energy << "}, "
//       << "newPx={" << new_lower_px << "," << new_upper_px << "}, original_lower_energy=" << original_lower_energy
//       << ", isfinal=" << isfinal << endl;
  D3SpectrumDisplayDiv *spectrum = this;
  
  //if( !spectrum )
  //  return;
  
  PeakModel *peakModel = spectrum->m_peakModel;
  shared_ptr<const Measurement> foreground = spectrum->data();
  
  try
  {
    if( !peakModel || !foreground )  //Shouldnt ever happen
      return;
    
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> origpeaks = peakModel->peaks();
    if( !origpeaks )
      return;  //shouldnt ever happen
    
    double minDe = 999999.9;
    std::shared_ptr<const PeakContinuum> continuum;
    for( auto p : *origpeaks )
    {
      const double de = fabs( p->continuum()->lowerEnergy() - original_lower_energy );
      if( de < minDe )
      {
        minDe = de;
        continuum = p->continuum();
      }
    }//for( auto p : *origpeaks )
    
    if( !continuum || minDe > 1.0 )  //0.001 would probably be fine instead of 1.0
    {
      spectrum->updateRoiBeingDragged( {} );
      
      throw runtime_error( "Couldnt find a continuum with lower energy " + std::to_string(original_lower_energy)
                          + " (de=" + std::to_string(minDe) + ")" );
    }//if( failed to find continuum )
    
  
    auto new_continuum = std::make_shared<PeakContinuum>( *continuum );
    
    //Lets re-use the c++ ROI value that we arent dragging, to avoid some cycle
    //  of rounding that could lead us to the end not being changed ending up
    //  with a different value than initially.
    const bool dragginUpperEnergy = (fabs(new_lower_energy - continuum->lowerEnergy()) < fabs(new_upper_energy - continuum->upperEnergy()));
    
    if( dragginUpperEnergy )
      new_lower_energy = continuum->lowerEnergy();
    else
      new_upper_energy = continuum->upperEnergy();
    
    new_continuum->setRange( new_lower_energy, new_upper_energy );
    
    bool allPeaksInRoiAreGaus = true;
    vector< shared_ptr<const PeakDef>> new_roi_initial_peaks, orig_roi_peaks;
    for( auto p : *origpeaks )
    {
      if( p->continuum() == continuum )
      {
        orig_roi_peaks.push_back( p );
        auto newpeak = make_shared<PeakDef>(*p);
        newpeak->setContinuum( new_continuum );
        allPeaksInRoiAreGaus = (allPeaksInRoiAreGaus && newpeak->gausPeak());
        new_roi_initial_peaks.push_back( newpeak );
      }
    }
    
    //if( dragginUpperEnergy && newpeak->gau && new_upper_energy < (newpeak->mean() - newpeak->sigma()) )
    if( allPeaksInRoiAreGaus && (new_roi_initial_peaks.size() > 1) )
    {
      map<double,shared_ptr<const PeakDef>> distance_to_roi;
      
      //Add all peaks that are at least 1 sigma outside of the ROI bounds to distance_to_roi.
      for( auto p : new_roi_initial_peaks )
      {
        const double m = p->mean();
        if( (m+p->sigma()) < new_lower_energy || (m-p->sigma()) > new_upper_energy )
          distance_to_roi[-std::min(m - new_lower_energy, m-new_upper_energy)] = p;
      }
      
      //Erase peaks in distance_to_roi from new_roi_initial_peaks, starting with
      //  ones furthest outside of ROI bounds.  Make sure that there is at least
      //  one peak left though.
      for( const auto &dp : distance_to_roi )
      {
        if( new_roi_initial_peaks.size() > 1 )
          new_roi_initial_peaks.erase( std::find(begin(new_roi_initial_peaks),end(new_roi_initial_peaks),dp.second) );
      }
    }//if( new_roi_initial_peaks.size() > 1 )
    
    //If region to narrow, or fit fails, pass back null if !isFinal, or original ROI if isFinal
    //cout << "new_upper_px=" << new_upper_px << ", new_lower_px=" << new_lower_px << endl;
    if( !allPeaksInRoiAreGaus )
    {
      if( isfinal )
      {
        // TODO: add a `PeakModel::replacePeaks(orig,newer)` function to allow updating peak quanitites in one step, or allow setting kLowerX/kUpperX columns in setData(...) function
        peakModel->removePeaks( orig_roi_peaks );
        
        std::vector<PeakDef> peaks_to_add;
        for( auto p : new_roi_initial_peaks )
          peaks_to_add.push_back( *p );
        
        peakModel->addPeaks( peaks_to_add );
      }else
      {
        spectrum->updateRoiBeingDragged( new_roi_initial_peaks );
      }
    }else if( new_upper_px >= (new_lower_px+10) )  //perhaps this should be by percentage of ROI?
    {
      //Need to check that all peaks are Gaussian.
      //  Actually should make sure a ROI can only have data defined or Gausian peaks only (what happens now)
      std::shared_ptr<const DetectorPeakResponse> detector;
      
      // As user drag a ROI, the parameters may change substantially, but we are always comparing against the
      //  original peaks chi2/dof; instead we should compare to the last peaks fit that were used, and the original ones
      const auto now = chrono::steady_clock::now();
      const auto dt = now - m_last_drag_time;
      if( (m_continuum_being_drug == continuum)
         && (new_roi_initial_peaks.size() == m_last_being_drug_peaks.size())
         && !m_last_being_drug_peaks.empty()
         && (dt < std::chrono::seconds(15)) )
      {
        //cout << "Using previously drug-n-fit peaks" << endl;
        
        const auto old_continuum = m_last_being_drug_peaks.front()->continuum();
        assert( old_continuum );
        auto new_continuum = std::make_shared<PeakContinuum>( *old_continuum );
        new_continuum->setRange( new_lower_energy, new_upper_energy );
        
        new_roi_initial_peaks.clear();
        for( const auto &p : m_last_being_drug_peaks )
        {
          auto newpeak = make_shared<PeakDef>( *p );
          newpeak->setContinuum( new_continuum );
          new_roi_initial_peaks.push_back( newpeak );
        }
      }//if( we should start from the peaks we last fit while dragging )
      
      vector<shared_ptr<const PeakDef>> refitpeaks
                  = refitPeaksThatShareROI( foreground, detector, new_roi_initial_peaks, 3.0 );
      
      m_continuum_being_drug = continuum;
      m_last_being_drug_peaks = refitpeaks;
      m_last_drag_time = chrono::steady_clock::now();
      
      
      //If the fit failed, use the old peaks, but with the ROI changed to what the user has
      const auto &newpeaks = refitpeaks.empty() ? new_roi_initial_peaks : refitpeaks;
      
      if( isfinal )
      {
        // TODO: add a `PeakModel::replacePeaks(orig,newer)` function to allow updating peak quanitites in one step, or allow setting kLowerX/kUpperX columns in setData(...) function
        
        // TODO: or rather, to consider: we could `m_continuum_being_drug` use, so that we get exactly what is shown
        
        std::vector<PeakDef> peaks_to_add;
        for( auto p : newpeaks )
          peaks_to_add.push_back( *p );
        
        peakModel->removePeaks( orig_roi_peaks );
        peakModel->addPeaks( peaks_to_add );
        
        m_continuum_being_drug.reset();
        m_last_being_drug_peaks.clear();
        m_last_drag_time = chrono::steady_clock::time_point{};
      }else
      {
        spectrum->updateRoiBeingDragged( newpeaks );
      }
      
    }else
    {
      cout << "User wants to erase peaks" << endl;
      spectrum->updateRoiBeingDragged( {} );
      
      if( isfinal )
        peakModel->removePeaks( orig_roi_peaks );
    }//if( not narrow region ) / else
  }catch( std::exception &e )
  {
    cerr << "Caught exception: " << e.what() << endl;
  }//try / catch
  
}//performExistingRoiEdgeDragWork(...)


void D3SpectrumDisplayDiv::dragCreateRoiCallback( double lower_energy, double upper_energy,
                                                    int nForcedPeaks, bool isfinal,
                                                    double window_xpx, double window_ypx )
{  
  m_dragCreateRoi.emit( lower_energy, upper_energy, nForcedPeaks, isfinal, window_xpx, window_ypx );
}//dragCreateRoiCallback(...)


void D3SpectrumDisplayDiv::performDragCreateRoiWork( double lower_energy, double upper_energy,
                                                          int nForcedPeaks, bool isfinal,
                                                          double window_xpx, double window_ypx )
{
  /* ToDo:
     - try to use RSP if available (need to figure out how to get it here).
     - Put all of this in a separate function to post to the Wt Worker pool, and then push results
     - maybe keep state between calls to speed up subsequent calls
     - Really punish peaks being close together with no dip in-between to avoid the tendency to fit lots of peaks.
   */
  
  D3SpectrumDisplayDiv *spectrum = this;
  const bool allowAsync = true;
  
  assert( spectrum && spectrum->m_peakModel );
  if( !spectrum || !spectrum->m_peakModel )
  {
    updateRoiBeingDragged( {} );
    return;
  }
  
  if( upper_energy < lower_energy )
    std::swap( lower_energy, upper_energy );
  

  InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
  InterSpec *viewer = app ? app->viewer() : nullptr;
  PeakModel *peakModel = spectrum->m_peakModel;
  std::shared_ptr<const SpecMeas> meas = viewer ? viewer->measurment(SpecUtils::SpectrumType::Foreground) : nullptr;
  std::shared_ptr<const DetectorPeakResponse> detector = meas ? meas->detector() : nullptr;

  std::shared_ptr<const Measurement> foreground = spectrum->data();
  
  if( !foreground )
  {
    updateRoiBeingDragged( {} );
    return;
  }
  
  auto fcnworker = [foreground,detector,lower_energy,upper_energy,nForcedPeaks,isfinal,window_xpx,window_ypx,app,spectrum,peakModel](){
  
    const bool isHpge = PeakFitUtils::is_high_res( foreground );
    const float erange = upper_energy - lower_energy;
    const float midenergy = 0.5f*(lower_energy + upper_energy);
    
    try
    {
      //const auto start_cpu_time = SpecUtils::get_cpu_time();
      //const auto start_wall_time = SpecUtils::get_wall_time();
      
      float min_sigma_width_kev, max_sigma_width_kev;
      expected_peak_width_limits( midenergy, isHpge, foreground, min_sigma_width_kev, max_sigma_width_kev );
      
      if( erange < min_sigma_width_kev )
        throw runtime_error( "to small range" );
      
      int nPeaks = 1;
      const size_t start_channel = foreground->find_gamma_channel(lower_energy);
      const size_t end_channel = foreground->find_gamma_channel(upper_energy);
      
      std::vector< std::tuple<float,float,float> > candidate_peaks;
      secondDerivativePeakCanidates( foreground, start_channel, end_channel, candidate_peaks );
      const size_t ncandidates = candidate_peaks.size();
      
      //const auto derivative_peaks = secondDerivativePeakCanidatesWithROI( foreground, start_channel, end_channel );
      //const size_t ncandidates = derivative_peaks.size();
      
      nPeaks = static_cast<int>( ncandidates );
      nPeaks = std::max( nPeaks, 1 );
      nPeaks = std::min( nPeaks, static_cast<int>(2*erange/min_sigma_width_kev) );
      //Always allow for fitting at least three different number of peaks, even if only single core.
      // But note! Results are currently dependant on how many cores the users computer has!
      //  (I dont like this, in principle - need to do some benchmarking to see how long things
      //   take (on low end devices) to see if we can make the same for all devices)
      const size_t ncores = static_cast<size_t>( std::max(3, SpecUtilsAsync::num_physical_cpu_cores()) );
      
      vector<int> npeakstry;
      
      if( nForcedPeaks > 0 && nForcedPeaks < 10 )
      {
        npeakstry.push_back( nForcedPeaks );
      }else
      {
        if( nPeaks > 1 && ncores > 1 )
          npeakstry.push_back( nPeaks - 1 );
        
        npeakstry.push_back( nPeaks );
        
        if( npeakstry.size() < ncores && ((erange/(nPeaks + 1)) > min_sigma_width_kev) )
          npeakstry.push_back( nPeaks + 1 );
        
        if( (npeakstry.size() < ncores) && nPeaks > 2 )
          npeakstry.push_back( nPeaks - 2 );
        
        //if( npeakstry.size() < ncores && ((erange/(nPeaks + 2)) > min_sigma_width_kev) )
        //  npeakstry.push_back( nPeaks + 2 );
        
        //if( (npeakstry.size() < ncores) && nPeaks > 3 )
        //  npeakstry.push_back( nPeaks - 3 );
        
        //if( npeakstry.size() < ncores  && ((erange/(nPeaks + 3)) > min_sigma_width_kev)  )
        //  npeakstry.push_back( nPeaks + 3 );
      }//if( npeaks > 0 && npeaks < 10 )
      
      std::sort( begin(npeakstry), end(npeakstry) );
      
      vector<double> chi2s( npeakstry.size() );
      vector<vector<std::shared_ptr<const PeakDef> > > results( npeakstry.size() );
      
      //cout << "ncandidates=" << ncandidates << endl;
      
      auto worker = [&chi2s, &results, &npeakstry, ncandidates,
                     nForcedPeaks, lower_energy, upper_energy, foreground, detector]( const size_t index ){
        std::vector<std::shared_ptr<PeakDef> > newpeaks;
        const auto method = (static_cast<size_t>(npeakstry[index])==ncandidates)
        ? MultiPeakInitialGuessMethod::FromDataInitialGuess
        :  MultiPeakInitialGuessMethod::UniformInitialGuess;
        
        
        //Note: InterSpec::findPeakFromControlDrag(...) tries all methods, e.g.,
        //for( MultiPeakInitialGuessMethod method = MultiPeakInitialGuessMethod(0);
        //    method < FromInputPeaks;
        //    method = MultiPeakInitialGuessMethod(method+1) )
        //{
        //  ...
        //}
        //(we should do this, at least when the number of peaks is forced, however
        // then we could run into an issue of for the final fit, we could get a
        // slightly different answer since this (currently) will always force a
        // number of peaks)
        
        
        //This next function call duplicates a lot of work between threads - may
        //  be a place that can be optimized if it turns out to be needed.
        findPeaksInUserRange( lower_energy, upper_energy, npeakstry[index], method,
                             foreground, detector, newpeaks, chi2s[index] );
        
        //findPeaksInUserRange_linsubsolve( lower_energy, upper_energy, npeakstry[index], method,
        //                     foreground, detector, newpeaks, chi2s[index] );
        
        //Note: InterSpec::findPeakFromControlDrag(...) adjusts Chi2 as:
        //const size_t nbin = end_channel - start_channel;
        //const double dof = (nbin + 3*npeakstry[index]);
        //const double chi2Dof = chi2s[index] / dof;
        //for( size_t i = 0; i < newpeaks.size(); ++i )
        //  newpeaks[i]->set_coefficient( chi2Dof, PeakDef::Chi2DOF );
        //(havent checked if this is necassary)
        
        for( const auto &p : newpeaks )
          results[index].push_back( p );
      };//worker
      
      SpecUtilsAsync::ThreadPool pool;
      for( size_t index = 0; index < npeakstry.size(); ++index )
        pool.post( [worker,index](){ worker(index); } );
      
      pool.join();
      
      //const auto stop_cpu_time = SpecUtils::get_cpu_time();
      //const auto stop_wall_time = SpecUtils::get_wall_time();
      
      //Peak fitting can take up to a third of a second or so
      //cout << "Peak fitting took " << (stop_cpu_time-start_cpu_time)
      //     << " cpu seconds, and " << (stop_wall_time-start_wall_time)
      //     << " wall seconds" << endl;
      
      //pick the best chi2, but have some weighting to prefer less numbers of peaks
      int best_choice = -1;
      for( size_t i = 0; i < chi2s.size(); ++i )
      {
        if( results[i].empty() )
          continue;
        
        const int npeaksdiff = npeakstry[i] - ((best_choice >= 0) ? npeakstry[best_choice] : 0);
        
        //Require at least 10% better chi2 for each peak you want to add
        //  - this was completely drawn from the air - needs actual testing/optimizing.
        const double weight = std::max(0.25,(1.0 - 0.10*npeaksdiff));
        //cout << "npeakstry=" << npeakstry[i] << " had chi2=" << chi2s[i]
        //     << ", npeaksdiff=" << npeaksdiff << ", weight=" << weight
        //     << ", chi2s[" << best_choice << "]=" << chi2s[best_choice]
        //     << " npeakstry[" << best_choice << "]=" << npeakstry[best_choice]
        //     << endl;
        
        
        if( !IsNan(chi2s[i]) && !IsInf(chi2s[i]) 
           && ((best_choice < 0) || (chi2s[i] < weight*chi2s[best_choice])) )
        {
          best_choice = static_cast<int>(i);
        }
      }//for( loop over results to find best number of peaks )
      
      //cout << "Chose " << npeakstry[best_choice] << " peaks" << endl;
      
      if( (best_choice < 0) || results[best_choice].empty() )
        throw runtime_error( "Failed to fit for peaks." );
      
      WApplication::UpdateLock lock( app );
      
      if( !lock )
      {
        cerr << "Failed to get WApplication::UpdateLock in D3SpectrumDisplayDiv::dragCreateRoiCallback(...)" << endl;
        return;
      }
      
      InterSpec *viewer = app ? app->viewer() : nullptr;
      
      if( isfinal )
      {
        UndoRedoManager::PeakModelChange peak_undo_creator;
        
        deque< PeakModel::PeakShrdPtr > preaddpeaks, postaddpeaks;
        if( peakModel->peaks() ) //should always be true, but JIC
          preaddpeaks = *peakModel->peaks();
        
        std::vector<PeakDef> peaks_to_add;
        for( auto p : results[best_choice] )
        {
          PeakDef peak = *p;
          
          //Assign nuclide from reference lines
          if( viewer )
          {
            const bool showingEscape = viewer->showingFeatureMarker(FeatureMarkerType::EscapePeakMarker);
            const auto refwidget = viewer->referenceLinesWidget();
            const bool colorFromRefLines = viewer->colorPeaksBasedOnReferenceLines();
            PeakSearchGuiUtils::assign_nuclide_from_reference_lines( peak, peakModel,
                            foreground, refwidget, colorFromRefLines, showingEscape );
          }//if( viewer )
          
          peaks_to_add.emplace_back( std::move(peak) );
        }//for( loop over fit peaks and add them to peaks_to_add and assign nuclides )
        
        peakModel->addPeaks( peaks_to_add );
          
        if( peakModel->peaks() ) //should always be true, but JIC
          postaddpeaks = *peakModel->peaks();
        
        vector< PeakModel::PeakShrdPtr > added_peaks;
        for( const auto &p : postaddpeaks )
        {
          if( !std::count(begin(preaddpeaks), end(preaddpeaks), p) )
            added_peaks.push_back( p );
        }
        
        // If the coordinates of the event that triggered this fitting are on the page
        // (i.e. at least one not negative), then pop-up a menu to allow the user to
        //  select number of peaks.
        if( window_xpx >= 0 || window_ypx >= 0 )
        {
          DeleteOnClosePopupMenu *menu = new DeleteOnClosePopupMenu( nullptr, PopupDivMenu::TransientMenu );
          menu->aboutToHide().connect( menu, &DeleteOnClosePopupMenu::markForDelete );

          PopupDivMenuItem *item = nullptr, *selecteditem = nullptr;
          
          menu->setPositionScheme( Wt::Absolute );
          if( menu->isMobile() )
            item = menu->addPhoneBackItem( NULL );
          
          item = menu->addMenuItem( "Peaks To Keep In ROI:" );
          item->disable();
          item->setSelectable( false );
          menu->addSeparator();
          
          const vector<const char *> numnames{ "Cancel", "Single Peak", "Two Peaks",
            "Three Peaks", "Four Peaks", "Five Peaks", "Six Peaks", "Seven Peaks",
            "Eight Peaks", "Nine Peaks", "Ten Peaks", "Eleven Peaks"
          };
          
          for( size_t i = 0; i < (peaks_to_add.size() + 3); ++i )
          {
            string name = ((i < numnames.size()) ? std::string(numnames[i]) : (std::to_string(i) + " Peaks"));
            if( i == npeakstry[best_choice] )
              name = "Keep " + name;
            item = menu->addMenuItem( name );
            if( i == npeakstry[best_choice] )
            {
              item->setFocus();
              item->decorationStyle().font().setWeight( Wt::WFont::Weight::Bold );
              //item->decorationStyle().setBackgroundColor(<#WColor color#>)
              selecteditem = item;
            }
            
            
            item->triggered().connect( std::bind( [=](){
              menu->setHidden( true );
              
              if( i == npeakstry[best_choice] )
              {
                //m_peakModel->addPeaks( peaks_to_add );
                return;
              }
              
              try
              {
                peakModel->removePeaks( added_peaks );
              }catch(std::exception &e)
              {
                cerr << "Unexpected error removing peaks - must not be a valid peak any more...: "
                     << e.what() << endl;
              }
                
              if( i > 0 )
                spectrum->dragCreateRoiCallback( lower_energy, upper_energy, static_cast<int>(i), true, window_xpx, window_ypx );
            }) );
          }//for( size_t i = 0; i < (peaks_to_add.size() + 3); ++i )
          
          if( selecteditem )
            menu->select( selecteditem );
          
          if( menu->isMobile() )
          {
            menu->addStyleClass( " Wt-popupmenu Wt-outset" );
            menu->showMobile();
          }else
          {
            menu->addStyleClass( " Wt-popupmenu Wt-outset NumPeakSelect" );
            menu->popup( WPoint(window_xpx-30,window_ypx-30) );
            //menu->popup( WPoint(200,200) );
          }
        }
      }else
      {
        vector<shared_ptr<const PeakDef> > peaks( begin(results[best_choice]), end(results[best_choice]) );
        spectrum->updateRoiBeingDragged( peaks );
      }
      
      app->triggerUpdate();
    }catch( std::exception &e )
    {
      WApplication::UpdateLock lock( app );
      
      if( !lock )
      {
        cerr << "Failed to get WApplication::UpdateLock in (2) D3SpectrumDisplayDiv::dragCreateRoiCallback(...)" << endl;
        return;
      }
      
      cerr << "D3SpectrumDisplayDiv::dragCreateRoiCallback: caught exception  " << e.what() << endl;
      
      // If the user hasnt dragged at least 1 keV, dont try to fit things
      if( (lower_energy + 1.0) < upper_energy )
      {
        try
        {
          PeakDef tmppeak(midenergy, 0.5*erange, 0);
          std::shared_ptr<PeakContinuum> cont = tmppeak.continuum();
          cont->calc_linear_continuum_eqn( foreground, midenergy, lower_energy, upper_energy, 2, 2 );
          
          std::vector<std::shared_ptr<const PeakDef> > peaks{ make_shared<const PeakDef>(tmppeak) };
          
          spectrum->updateRoiBeingDragged( peaks );
        }catch( std::exception &e )
        {
          cerr << "D3SpectrumDisplayDiv::dragCreateRoiCallback: couldnt ." << endl;
        }//try / catch
      }else
      {
        spectrum->updateRoiBeingDragged( {} );
      }
      
      app->triggerUpdate();
    }//try / catch
  };//fcnworker
  
  if( allowAsync )
  {
    //Wt::WServer::instance()->post( wApp->sessionId(), fcnworker );
    WServer::instance()->ioService().boost::asio::io_service::post( fcnworker );
  }else
  {
    fcnworker();
  }
}//void performDragCreateRoiWork(...)


void D3SpectrumDisplayDiv::updateRoiBeingDragged( const vector<shared_ptr<const PeakDef> > &peaks )
{
  const shared_ptr<const Measurement> &fore = m_foreground;
  const string json = peaks.empty() ? string("null") : PeakDef::gaus_peaks_to_json( peaks, fore );
  
  doJavaScript( "try{" + m_jsgraph + ".updateRoiBeingDragged(" + json + ");}catch(error){}" );
}//void updateRoiBeingDragged( vector<shared_ptr<const PeakDef> > &roiBeingDragged )


void D3SpectrumDisplayDiv::yAxisScaled( const double scale, const std::string &spectrum )
{
  SpecUtils::SpectrumType type;

  //Dont call D3SpectrumDisplayDiv::setDisplayScaleFactor(...) since we dont
  //  have to re-load data to client, but we should keep all the c++ up to date.
  double previous_scale = 1.0;
  if( spectrum == "FOREGROUND" )
  {
    type = SpecUtils::SpectrumType::Foreground;
    previous_scale = 1.0;
  }else if( spectrum == "BACKGROUND" )
  {
    type = SpecUtils::SpectrumType::Background;
    previous_scale = m_backgroundScale;
    m_backgroundScale = scale;
  }else if( spectrum == "SECONDARY" )
  {
    type = SpecUtils::SpectrumType::SecondForeground;
    previous_scale = m_secondaryScale;
    m_secondaryScale = scale;
  }else
  {
    cerr << "Received yscaled signal with scale " << scale << " and spectrum = '"
    << spectrum << "', which is invalid" << endl;
    return;
  }
  
  // Add Undo/Redo step here.
  
  m_yAxisScaled.emit( scale, previous_scale, type );
}//void yAxisScaled( const double scale, const std::string &spectrum )



void D3SpectrumDisplayDiv::chartXRangeChangedCallback( double x0, double x1,
                                                      double chart_width_px, double chart_height_px,
                                                      bool user_action )
{
  if( (fabs(m_xAxisMinimum - x0) < 0.001)
     && (fabs(m_xAxisMaximum - x1) < 0.001)
     && (fabs(m_chartWidthPx - chart_width_px) < 0.001)
     && (fabs(m_chartHeightPx - chart_height_px) < 0.001) )
  {
    cout << "No appreciable change in x-range or chart pixel, not emitting" << endl;
    return;
  }
  
  const double oldXmin = m_xAxisMinimum, oldXmax = m_xAxisMaximum;
  
  // cout << "chartXRangeChangedCallback{" << x0 << "," << x1 << "," << chart_width_px << "," << chart_height_px << "}" << endl;
  m_xAxisMinimum = x0;
  m_xAxisMaximum = x1;
  m_chartWidthPx = chart_width_px;
  m_chartHeightPx = chart_height_px;
  
  m_xRangeChanged.emit( x0, x1, oldXmin, oldXmax, user_action );
}//void D3SpectrumDisplayDiv::chartXRangeChangedCallback(...)


void D3SpectrumDisplayDiv::sliderChartDisplayedCallback( const bool madeVisisble )
{
  m_showXAxisSliderChart = madeVisisble;
  m_xAxisSliderShown.emit( madeVisisble );
}//void sliderChartDisplayedCallback( const bool madeVisisble );


void D3SpectrumDisplayDiv::yAxisTypeChangedCallback( const std::string &type )
{
  const bool isLogY = (type == "log");
  if( isLogY == m_yAxisIsLog )
    return;
  
  m_yAxisIsLog = isLogY;
  m_yAxisLogLinChanged.emit( isLogY );
}//void yAxisTypeChanged( const std::string &type )


D3SpectrumDisplayDiv::~D3SpectrumDisplayDiv()
{
  //doJavaScript( "try{" + m_jsgraph + "=null;}catch(){}" );
}//~D3SpectrumDisplayDiv()


