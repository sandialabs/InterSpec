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

// Disable some warnings in xutility
#pragma warning(disable:4996)

#include <fstream>
#include <iostream>
#include <iterator>

#include <Wt/WColor>
#include <Wt/WBorder>
#include <Wt/WResource>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WPaintDevice>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WPaintedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WAbstractChart>

#include "SpecUtils/StringAlgo.h"
#include "InterSpec/CanvasForDragging.h"
#include "js/CanvasForDragging.js"


using namespace Wt;
using namespace std;

/*
 * See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
 */
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


#if( USE_OverlayDragEvent )
void OverlayDragEvent::clear()
{
  button = keyCode = charCode = -1;
  x0 = x1 = y0 = y1 = offsetLeft = offsetTop = wheelDelta = keyModifiers = 0;
}//void OverlayDragEvent::clear()


//USE_FAST_EVENT_DECODE: faster decoding time by using c funtions rather than
//  boost.  Also results in about a 3.5 kb binary size reduction.  Not
//  incredably well tested as of 20140210.
//  The slower function is  abit safer in terms of recognizing invalid input.
#define USE_FAST_EVENT_DECODE 1

namespace
{
#if( USE_FAST_EVENT_DECODE )
bool isnullstr( const char *str )
{
  if( !str )
    return true;
  
  size_t len = 0;
  while( str[len] && str[len]!='&' )
    len++;
  
  if( len != 4 )
    return false;
  return (str[0]=='n' && str[1]=='u' && str[2]=='l' && str[3]=='l');
}//bool isnullstr( const char *str )
#else
bool str_to_bool( const string &str )
{
  return (str=="true" || str=="1");
}
#endif //#if( USE_FAST_EVENT_DECODE )
}//namespace


bool operator>>( std::istream &is, OverlayDragEvent& t )
{
#if( USE_FAST_EVENT_DECODE )
  t.clear();
  OverlayDragEvent &a = t;
  string input;
  is >> input;
  
  if( input.empty() )
    return false;
  
  vector<const char *> fieldstrs;
  fieldstrs.reserve( 14 );
  fieldstrs.push_back( input.c_str() );
  for( size_t i = input.find('&',1); i < string::npos; i = input.find('&',i+1) )
    fieldstrs.push_back( input.c_str() + i  + 1 );
  
  if( fieldstrs.size() != 15 )
  {
    cerr << "CanvasForDragging operator>>:"
         << "\n\tbool operator>>( std::istream &is, OverlayDragEvent& t ):\n\t"
         << "Recieved an input with " << fieldstrs.size() << " fields; I expected 15"
         << endl;
    return false;
  }//if( fieldstrs.size() != 15 )
  
  //'null&230&null&585&1&0&0&0&0&0&0&-1'
  //Note that we cast first to a float, then to an int because we nominally
  //  expect a int, but some browsers give pixel positions in floats, and I
  //  think I round to ints everywhere in JavaScript, but just to be surea
  int i = 0;
  const char *arg = fieldstrs[i++];
  a.dt = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
  arg = fieldstrs[i++];
  a.x0 = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
  arg = fieldstrs[i++];
  a.x1 = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
  arg = fieldstrs[i++];
  a.y0 = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
  arg = fieldstrs[i++];
  a.y1 = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
  arg = fieldstrs[i++];
  a.offsetLeft = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
  arg = fieldstrs[i++];
  a.offsetTop = (isnullstr(arg) ? -1 : static_cast<int>(atof(arg)) );
    
  a.button   = atoi( fieldstrs[i++] );
  a.keyCode  = atoi( fieldstrs[i++] );
  a.charCode = atoi( fieldstrs[i++] );
  const bool altKey   = (atoi(fieldstrs[i++]) > 0);
  const bool ctrlKey  = (atoi(fieldstrs[i++]) > 0);
  const bool metaKey  = (atoi(fieldstrs[i++]) > 0);
  const bool shiftKey = (atoi(fieldstrs[i++]) > 0);
    
  a.keyModifiers = 0;
  if( altKey )
    a.keyModifiers |= Wt::AltModifier;
  if( ctrlKey )
    a.keyModifiers |= Wt::ControlModifier;
  if( metaKey )
    a.keyModifiers |= Wt::MetaModifier;
  if( shiftKey )
    a.keyModifiers |= Wt::ShiftModifier;
  a.wheelDelta = atoi( fieldstrs[i++] );
  
  return !!is;
#else //#if( USE_FAST_EVENT_DECODE )
  
  t.clear();
  
  string arg;
  is >> arg;
  vector<string> fields;
  SpecUtils::split( fields, arg, "&" );
  
  if( fields.size() != 15 )
  {
    cerr << "\nbool CanvasForDragging::operator>>( std::istream &is, OverlayDragEvent& t ):\n\t"
    << "Recieved an input with " << fields.size() << " fields; I expected 15"
    << endl;
    return !!is;
  }//if( fields.size() != 15 )
  
  OverlayDragEvent a;
  
  try
  {
    //'null&230&null&585&1&0&0&0&0&0&0&-1'
    //Note that we cast first to a float, then to an int because we nominally
    //  expect a int, but some browsers give pixel positions in floats, and I
    //  think I round to ints everywhere in JavaScript, but just to be surea
    int i = 0;
    string arg = fields[i++];
    a.dt = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    arg = fields[i++];
    a.x0 = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    arg = fields[i++];
    a.x1 = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    arg = fields[i++];
    a.y0 = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    arg = fields[i++];
    a.y1 = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    arg = fields[i++];
    a.offsetLeft = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    arg = fields[i++];
    a.offsetTop = (arg=="null" ? -1 : static_cast<int>(std::stod( arg )) );
    
    a.button   = static_cast<int>( std::stod(fields[i++]) );
    a.keyCode  = static_cast<int>( std::stod(fields[i++]) );
    a.charCode = static_cast<int>( std::stod(fields[i++]) );
    const bool altKey   = str_to_bool( fields[i++] );
    const bool ctrlKey  = str_to_bool( fields[i++] );
    const bool metaKey  = str_to_bool( fields[i++] );
    const bool shiftKey = str_to_bool( fields[i++] );
    
    a.keyModifiers = 0;
    if( altKey )
      a.keyModifiers |= Wt::AltModifier;
    if( ctrlKey )
      a.keyModifiers |= Wt::ControlModifier;
    if( metaKey )
      a.keyModifiers |= Wt::MetaModifier;
    if( shiftKey )
      a.keyModifiers |= Wt::ShiftModifier;
    
    a.wheelDelta = std::stoi( fields[i++] );
  }catch(...)
  {
    cerr << endl
    << "CanvasForDragging.cpp"
    << "\n\tbool operator>>( std::istream &is, OverlayDragEvent& t ):\n\t"
    << "There was an error converting the stream to a OverlayDragEvent.\n"
    << "String from client: '" << arg << "'\n"
    << endl;
    return !!is;
  }//try / catch
  
  t = a;
  
  return !!is;
#endif
}//std::istream& operator<<( std::istream&, OverlayDragEvent& t)
#endif  //#if( USE_OverlayDragEvent )


#if( USE_HIGH_BANDWIDTH_INTERACTIONS )

void dummy_touchmoved( const WTouchEvent & )
{
}
#endif  //#if( USE_HIGH_BANDWIDTH_INTERACTIONS )

CanvasForDragging::CanvasForDragging( Wt::Chart::WAbstractChart *parent,
                                      bool outline, bool highlight,
                                      bool altShiftHighlight )
  : WPaintedWidget(),
    m_userDraggedSignal( NULL ),
    m_userSingleClickedSignal( NULL ),
    m_rightClickSignal( NULL ),
    m_controlMouseDown( NULL ),
    m_controlMouseMove( NULL ),
    m_dblTap( NULL ),
    m_jsException( NULL ),
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
    m_mousedown( NULL ),
    m_mouseup( NULL ),
    m_leftDownMouseMove( NULL ),
    m_altLeftDownMouseMove( NULL ),
    m_pinchZoomChange( NULL ),
#endif
    m_parent( parent )
{
  addStyleClass( "CanvasForDragging" );
  setPositionScheme( Absolute );
  setLoadLaterWhenInvisible( false );
  setPreferredMethod( WPaintedWidget::HtmlCanvas );
  setId( parent->id() + "Cover" );  //Doesnt work with Wt<3.2.0
  wApp->domRoot()->addWidget( this );
  show();
  
  //make it so we can use control-drag and right click on the cavas without
  //  the browsers context menu poping up
  setAttributeValue( "oncontextmenu", "return false;" );
  
#if( USE_OverlayDragEvent )
  m_userDraggedSignal             = new JSignal<OverlayDragEvent>( this, "userDragged", true );
#else
  m_userDraggedSignal             = new JSignal<int,int,WMouseEvent>( this, "userDragged", true );
#endif
  
  m_userSingleClickedSignal       = new JSignal<int,int,int,int,int>( this, "userSingleClicked", true );
  m_rightClickSignal              = new JSignal<int,int,int,int,int>( this, "rightClick", true );
//  m_keyPressWhileMousedOverSignal.reset( new JSignal<WKeyEvent>( this, "keyPressWhileMousedOver", true ) );
  m_controlMouseDown              = new JSignal<int>( this, "cntrlMouseDown", true );
  m_controlMouseMove              = new JSignal<int,int>( this, "cntrlMouseMove", true );
  m_dblTap                        = new JSignal<int,int>( this, "dblTap", true );
  m_jsException                   = new JSignal<std::string>( this, "jsException", true );

  
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  m_mousedown = new JSignal<int,int,int>( this, "UserMouseDown", true );
  m_mouseup = new JSignal<int,int,int,int>( this, "UserMouseUp", true );
  m_altLeftDownMouseMove = new JSignal<int,int,int>( this, "UserMouseAltMove", true );
  m_leftDownMouseMove = new JSignal<int,int,int>( this, "UserMouseLeftMove", true );
  m_pinchZoomChange = new JSignal<int,int,int,int,int>( this, "PinchZoom", true );
  
  //For some reason if we dont connect dummy_touchmoved here, then the
  //  'UserMouseAltMove' signal wont get emmitted from OverlayTouchChange
  //  reliably on iOS (seems fine on android, but well enable it JIC)
#if( ANDROID || IOS )
  touchMoved().connect( boost::bind(&dummy_touchmoved, _1) );
#endif
#endif

  clicked().preventPropagation();
  clicked().preventDefaultAction();
  doubleClicked().preventPropagation();
  doubleClicked().preventDefaultAction();
  
  clicked().connect(        boost::bind( &EventSignal< WMouseEvent >::emit,
                                         &(parent->clicked ()), _1 ) );
  doubleClicked().connect(  boost::bind( &EventSignal< WMouseEvent >::emit,
                                         &(parent->doubleClicked ()), _1 ) );
  
  //hide the peak info tip when the mouse leaves the canvas
  mouseWentOut().connect( "function(s,e){$('#'+s.id+'_tip').hide();}" );
  
  setChartPadding();
  setInteractionMode( outline, highlight, altShiftHighlight );
  loadJs();
  
  
  //Note that if a WPaintedWidget is in a layout, then some additional
  //  instructions will be needed, namely:
//  """var u = $(self).find('canvas, img');"
//  """if (w >= 0) "
//  ""  "u.width(w);"
//  """if (h >= 0) "
//  ""  "u.height(h);"
          
  //Right now I'm only redrawing the gamma lines, but in principle there could
  //  be other things that should be redrawn (and that will happen now only once
  //  the mouse is moved)
  string js = "function(self, w, h, layout){Wt.WT.DrawGammaLines('c" + id() + "', true);}";
  setJavaScriptMember( WT_RESIZE_JS, js );
}//CanvasForDragging constructos


WPaintedWidget::Method CanvasForDragging::getMethod() const
{
  return Wt::WPaintedWidget::HtmlCanvas; 
}


void CanvasForDragging::setInteractionMode( bool outline, bool highlight,
                                            bool altShiftHighlight )
{

  if( !m_setInteractionModeSlot )
  {
    const char *js = INLINE_JAVASCRIPT
    (
     function(o,mode)
     {
       try{$('#c'+o.id).data('drawMode',mode);}
       catch(e){console.log('Failed to set chart drawing mode');}
     }
    );
    
    m_setInteractionModeSlot.reset( new JSlot( js, this ) );
  }//if( !m_setInteractionModeSlot )

  
  const string drawMode
        = "{ highlight: " + string(highlight?"true":"false")
           + ", outline: " + string(outline?"true":"false")
           + ", altShiftHighlight: " + string(altShiftHighlight?"true":"false")
       + " }";
  m_setInteractionModeSlot->exec( id(), drawMode );
//  doJavaScript( "$('#c" + id() + "').data('drawMode'," + drawMode + ");" );
}//void setInteractionMode(...)


void CanvasForDragging::setShowRefLineInfoForMouseOver( const bool show )
{
  const string jsref = "$('#c" + id() + "')";
  const string jsshow = show ? "true" : "null";
  doJavaScript( jsref + ".data('NoRefLineInfo'," + jsshow + ");" );
}//setShowRefLineInfoForMouseOver( const bool show )


void CanvasForDragging::setNoSpectrumManipulationMode( const bool nomanip )
{
  const string jsref = "$('#c" + id() + "')";
  const string jsshow = nomanip ? "true" : "null";
  doJavaScript( jsref + ".data('NoSpecManip'," + jsshow + ");" );
}//void setNoSpectrumManipulationMode( const bool nomanip )


void CanvasForDragging::setChartPadding()
{
  if( !m_parent )
    throw runtime_error( "CanvasForDragging::setChartPadding(): expected "
                         "m_parent to be valid" );

  if( !m_setChartPaddingSlot )
  {
    const char *js = INLINE_JAVASCRIPT
    (
      function(o,padding)
      {
        try{$('#c'+o.id).data('chartPadding',padding);}
        catch(e){ console.log('Failed to set chart padding');}
      }
    );
    m_setChartPaddingSlot.reset( new JSlot( js, this ) );
  }//if( !m_setChartPaddingSlot )
  
  // Left off here: 3/26/12 - change the absolute css position (bottom, top) if
  // the user scrolls - need to do this to not occlude the tabs from being
  // selected.
  const string cbottom = std::to_string( m_parent->plotAreaPadding(Bottom) );
  const string ctop    = std::to_string( m_parent->plotAreaPadding(Top) );
  const string cleft   = std::to_string( m_parent->plotAreaPadding(Left) );
  const string cright  = std::to_string( m_parent->plotAreaPadding(Right) );
  const string padding = "{ top: " + ctop + ", bottom: " + cbottom
  + ", left: " + cleft + ", right: " + cright + " }";
  
  m_setChartPaddingSlot->exec( id(), padding );
//  doJavaScript( "$('#c" + id() + "').data('chartPadding'," + padding + ");" );
}//void CanvasForDragging::setChartPadding()
  

void CanvasForDragging::loadJs()
{
  WApplication *app = WApplication::instance();

  //LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsAlignOverlay);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayUpdateMouseCoords);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsDrawExpectedPeakConsequences);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsDrawSearchEnergies);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsDrawGammaLines);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayShowPeakTip);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsDrawContEst);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsIsDeletePeakSwipe);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsIsControlDragSwipe);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsIsAltShiftSwipe);
//  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsDrawZoomingOut);
  
  
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnClick);
  const string onclickSlotJS = "function(s,e){this.WT.OverlayOnClick(s,e);}";
  std::shared_ptr<JSlot> onclickSlot( new JSlot( onclickSlotJS, this ) );
  clicked().connect( *onclickSlot );
  m_jslots.push_back( onclickSlot );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnMouseDown);
  const string onMouseDownSlotJS = "function(s,e){this.WT.OverlayOnMouseDown(s,e);}";
  std::shared_ptr<JSlot> onmousedown( new JSlot( onMouseDownSlotJS, this ) );
  mouseWentDown().connect( *onmousedown );
  m_jslots.push_back( onmousedown );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnMouseUp);
  const string onMouseUpSlotJS = "function(s,e){this.WT.OverlayOnMouseUp(s,e);}";
  std::shared_ptr<JSlot> onmouseup( new JSlot( onMouseUpSlotJS, this ) );
  mouseWentUp().connect( *onmouseup );
  m_jslots.push_back( onmouseup );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnMouseOut);
  const string onMouseOutSlotJS = "function(s,e){this.WT.OverlayOnMouseOut(s.id);}";
  std::shared_ptr<JSlot> onmouseout( new JSlot( onMouseOutSlotJS, this ) );
  mouseWentOut().connect( *onmouseout );
  m_jslots.push_back( onmouseout );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnMouseOver);
  const string onMouseOverSlotJS = "function(s,e){this.WT.OverlayOnMouseOver(s,e);}";
  std::shared_ptr<JSlot> onmouseoverSlot( new JSlot( onMouseOverSlotJS, this ) );
  mouseWentOver().connect( *onmouseoverSlot );
  m_jslots.push_back( onmouseoverSlot );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnMouseMove);
  const string onMouseMoveSlotJS = "function(s,e){this.WT.OverlayOnMouseMove(s,e);}";
  std::shared_ptr<JSlot> onmousemoveSlot( new JSlot( onMouseMoveSlotJS, this ) );
  mouseMoved().connect( *onmousemoveSlot );
  m_jslots.push_back( onmousemoveSlot );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnKeyPress);
  const string onEscKeyPressSlotJS = "function(s,e){this.WT.OverlayOnKeyPress('" + id() + "',e.target,e,false);}";
  std::shared_ptr<JSlot> onEscPressSlot( new JSlot( onEscKeyPressSlotJS, this ) );
  app->globalEscapePressed().connect( *onEscPressSlot );
  m_jslots.push_back( onEscPressSlot );

  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayTouchBegin);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayTouchChange);
  LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayTouchEnd);
    
  std::shared_ptr<JSlot> onTouchBeginSlot( new JSlot( "function(s,e){this.WT.OverlayTouchBegin(s,e);}", this ) );
  std::shared_ptr<JSlot> onTouchChangeSlot( new JSlot( "function(s,e){this.WT.OverlayTouchChange(s,e);}", this ) );
  std::shared_ptr<JSlot> onTouchEndSlot( new JSlot( "function(s,e){this.WT.OverlayTouchEnd(s,e);}", this ) );
  touchStarted().connect( *onTouchBeginSlot );
  touchMoved().connect( *onTouchChangeSlot );
  touchEnded().connect( *onTouchEndSlot );
  m_jslots.push_back( onTouchBeginSlot );
  m_jslots.push_back( onTouchChangeSlot );
  m_jslots.push_back( onTouchEndSlot );
  
  
  doJavaScript("$('#c" + id() + "').data('startDragX',null);"
               "$('#c" + id() + "').data('startDragY',null);");
  
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  doJavaScript( "$('#c" + id() + "').data('HighBandwidth',true);" );
#else
  doJavaScript( "$('#c" + id() + "').data('HighBandwidth',false);" );
#endif
}//void CanvasForDragging::loadJs()

void CanvasForDragging::setControlDragDebounceTimePeriod( int milliseconds )
{
  if( milliseconds <= 0 )
  {
    m_controlDragDebounceTimePeriod = 0;
    doJavaScript("$('#c" + id() + "').data('ctrlKeyUpdateDt',null);" );
  }else
  {   
    m_controlDragDebounceTimePeriod = milliseconds;
    doJavaScript("$('#c" + id() + "').data('ctrlKeyUpdateDt',"
                 + std::to_string(milliseconds) + ");" );
  }
}//void setControlDragDebounceTimePeriod( int milliseconds )

int CanvasForDragging::controlDragDebounceTimePeriod() const
{
  return m_controlDragDebounceTimePeriod;
}

void CanvasForDragging::setScrollingParent( Wt::WContainerWidget *parent )
{
  string js;
  if( parent )
    js = "$('#" + id() + "').data('scrollParent','" + parent->id() + "');";
  else
    js = "$('#" + id() + "').data('scrollParent',null);";

  doJavaScript( js );
}//void CanvasForDragging::setScrollingParent( Wt::WContainerWidget *parent )

#if( USE_OverlayDragEvent )
Wt::JSignal<OverlayDragEvent> &CanvasForDragging::userDragged()
{
  return *m_userDraggedSignal;
}
#else
Wt::JSignal<int,int,WMouseEvent> &CanvasForDragging::userDragged()
{
  return *m_userDraggedSignal;
}
#endif

Wt::JSignal<int,int,int,int,int> &CanvasForDragging::userSingleClicked()
{
  return *m_userSingleClickedSignal;
}

Wt::JSignal<int,int,int,int,int> &CanvasForDragging::rightClickSignal()
{
  return *m_rightClickSignal;
}

Wt::JSignal<Wt::WKeyEvent> &CanvasForDragging::keyPressWhileMousedOver()
{
  if( !m_keyPressWhileMousedOverSignal )
  {
    m_keyPressWhileMousedOverSignal.reset( new JSignal<WKeyEvent>( this, "keyPressWhileMousedOver", true ) );
  
    //A canvas element cant have focus, so we actually need to use the document
    //  for this.  The javascript has a check to make sure the mouse is in this
    //  canvas before it emits a signal
    //LOAD_JAVASCRIPT(app, "js/CanvasForDragging.js", "CanvasForDragging", wtjsOverlayOnKeyPress);  //already taken care of in the constructor
    const string onKeyPressSlotJQuery = "$(document).keyup( function(e){ " + wApp->javaScriptClass() + ".WT.OverlayOnKeyPress('" + id() + "',e.target,e,true); } );";
    doJavaScript( onKeyPressSlotJQuery );
  }//if( !m_keyPressWhileMousedOverSignal )
  
  return *m_keyPressWhileMousedOverSignal;
}

Wt::JSignal<int> &CanvasForDragging::controlMouseDown()
{
  return *m_controlMouseDown;
}

Wt::JSignal<int,int> &CanvasForDragging::controlMouseMove()
{
  return *m_controlMouseMove;
}


Wt::JSignal<int,int> *CanvasForDragging::dblTap()
{
  return m_dblTap;
}

Wt::JSignal<std::string> *CanvasForDragging::jsException()
{
  return m_jsException;
}

#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
Wt::JSignal<int/*x*/,int/*y*/,int/*mods*/>  &CanvasForDragging::mousedown()
{
  return *m_mousedown;
}//mousedown()


Wt::JSignal<int/*x*/,int/*y*/,int/*mods*/,int/*dt in ms*/>  &CanvasForDragging::mouseup()
{
  return *m_mouseup;
}//mouseup()


Wt::JSignal<int/*x*/,int/*y*/,int/*dt in ms*/>  &CanvasForDragging::leftDownMouseMove()
{
  return *m_leftDownMouseMove;
}//leftDownMouseMove()


Wt::JSignal<int/*x*/,int/*y*/,int/*dt in ms*/>  &CanvasForDragging::altLeftDownMouseMove()
{
  return *m_altLeftDownMouseMove;
}//altLeftDownMouseMove()


Wt::JSignal<int,int/*x_t0*/,int/*x0_t1*/,int/*x_t1*/,int/*dir*/>  &CanvasForDragging::pinchZoomChange()
{
  return *m_pinchZoomChange;
}
#endif


CanvasForDragging::~CanvasForDragging()
{
  if( m_userDraggedSignal )
    delete m_userDraggedSignal;
  m_userDraggedSignal = NULL;

  if( m_userSingleClickedSignal )
    delete m_userSingleClickedSignal;
  m_userSingleClickedSignal = NULL;
  
  if( m_rightClickSignal )
    delete m_rightClickSignal;
  m_rightClickSignal = NULL;

  if( m_controlMouseDown )
    delete m_controlMouseDown;
  m_controlMouseDown = NULL;

  if( m_controlMouseMove )
    delete m_controlMouseMove;
  m_controlMouseMove = NULL;

  if( m_dblTap )
    delete m_dblTap;
  m_dblTap = NULL;
  
  if( m_jsException )
    delete m_jsException;
  m_jsException = NULL;

#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  if( m_mousedown )
    delete m_mousedown;
  m_mousedown = NULL;
  
  if( m_mouseup )
    delete m_mouseup;
  m_mouseup = NULL;
  
  if( m_leftDownMouseMove )
    delete m_leftDownMouseMove;
  m_leftDownMouseMove = NULL;
  
  if( m_altLeftDownMouseMove )
    delete m_altLeftDownMouseMove;
  m_altLeftDownMouseMove = NULL;

  if( m_pinchZoomChange )
    delete m_pinchZoomChange;
  m_pinchZoomChange = NULL;
#endif

  
}//~WPaintedWidget()


void CanvasForDragging::paintEvent( Wt::WPaintDevice * /*paintDevice*/ )
{
}//void paintEvent( Wt::WPaintDevice *paintDevice )

