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

#include "InterSpec/AuxWindow.h"
#include <Wt/WText.h>
#include <Wt/WTheme.h>
#include <Wt/WPanel.h>
#include <Wt/WImage.h>
#include <Wt/WIconPair.h>
#include <Wt/WTemplate.h>
#include <Wt/WBoxLayout.h>
#include <Wt/WBorderLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WJavaScript.h>
#include <Wt/WEnvironment.h>
#include <Wt/WApplication.h>
#include <Wt/WStringStream.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/WarningWidget.h"

#if( BUILD_AS_WX_WIDGETS_APP || BUILD_AS_ELECTRON_APP )
#include "InterSpec/InterSpecApp.h"
#endif


using namespace std;
using namespace Wt;


// An experiment to handle adjusting window position/size when browser window changes
//  So far seems to work well, except it creates problems on Android when the keyboard
//  shows, which resizes the whole screen, and causes problems, atm.
#if( ANDROID )
  #define AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE 0
#else
  #define AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE 1
#endif


#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


WT_DECLARE_WT_MEMBER
(AuxWindowResizeToFitOnScreen, Wt::JavaScriptFunction, "AuxWindowResizeToFitOnScreen",
function(id)
{
  try
  {
    const el = document.getElementById(id);
    const cel = document.querySelector('#' + id + ' .body.AuxWindow-content');
    if( !el || !cel )
    {
      console.log( 'Error in resizeToFitOnScreen' );
      return;
    }

    const ws = Wt.WT.windowSize();

    let h = el.offsetHeight;
    let w = el.offsetWidth;
    const do_resize = ((h > ws.y) || (w > ws.x));
    if(h > ws.y)
    {
      h = ws.y - 10;
      el.style.top = '5px';
      el.style.bottom = "";
      cel.style.overflowY = "auto";
    }

    if(w > ws.x)
    {
      w = ws.x - 10;
      el.style.left = '5px';
      el.style.right = "";
      cel.style.overflowX = "auto";
    }

    if( do_resize )
    {
      el._isData = el._isData || {};
      const was_centered = el._isData.centered;
      el.wtObj.onresize(w, h, true);

      if( was_centered )
        Wt.WT.AuxWindowCenter( el );

      el._isData.centered = was_centered;
    }
  }catch(err){
    console.error( 'AuxWindowResizeToFitOnScreen error:', err );
  }
}
 );//WT_DECLARE_WT_MEMBER( AuxWindowResizeToFitOnScreen )


WT_DECLARE_WT_MEMBER
(AuxWindowCenter, Wt::JavaScriptFunction, "AuxWindowCenter",
function( el )
{
  try
  {
    const olddispl = el.style.display;
    el.style.display = "";
    const ws = Wt.WT.windowSize();
    const hh = Math.max(0,Math.round( ( ws.y - el.offsetHeight ) / 2 ));
    const hw = Math.max(0,Math.round( ( ws.x - el.offsetWidth ) / 2 ));

    el.style.top = hh+'px';
    el.style.left = hw+'px';
    el.style.display = olddispl;
    el._isData = el._isData || {};
    el._isData.centered = true;
  }catch(err)
  {
    console.log( "Failed in AuxWindowCenter: " + err );
  }
}
);




WT_DECLARE_WT_MEMBER
(AuxWindowBringToFront, Wt::JavaScriptFunction, "AuxWindowBringToFront",
 function( id )
{
  var maxz = 0;
  document.querySelectorAll('.Wt-dialog,.MobileMenuButton').forEach( function(value){
    let z = parseInt(getComputedStyle(value).zIndex) || 0;
    //if( value.offsetParent !== null ) //seems to not be reliable, so leaving commented out.
    maxz = Math.max(maxz,z);
  } );
  var idEl = document.getElementById(id);
  var idZ = idEl ? (parseInt(getComputedStyle(idEl).zIndex) || 0) : 0;
  if( maxz > idZ )
  { if( idEl ) idEl.style.zIndex = maxz+1; }
  document.querySelectorAll('.window-controls-container').forEach(function(e){ e.style.zIndex = maxz+1; }); //for wxWidgets and Electron builds
  document.querySelectorAll('.suggestion').forEach(function(e){ e.style.zIndex = maxz+2; }); //This is for InterSpec::m_shieldingSuggestion
}
);


WT_DECLARE_WT_MEMBER
(AuxWindowTitlebarTouchStart, Wt::JavaScriptFunction, "AuxWindowTitlebarTouchStart",
function( sender, event, id )
{
  var e = event||window.event;
  var el = this.getElement(id);
  var titlebar = el.querySelector(".titlebar");

  el._isData = el._isData || {};
  if(el._isData.touchmoving == null)
  {
    Wt.WT.capture(titlebar);
    var touch = event.targetTouches[0];
    el._isData.dsx = touch.pageX;
    el._isData.dsy = touch.pageY;
    el._isData.touchmoving = true;
  }
}
);

WT_DECLARE_WT_MEMBER
(AuxWindowTitlebarTouchEnd, Wt::JavaScriptFunction, "AuxWindowTitlebarTouchEnd",
function( sender, e, id )
{
  if( e.targetTouches.length > 0 )
    return;

  var el = this.getElement(id);
  el._isData = el._isData || {};
  el._isData.dsx = null;
  el._isData.dsy = null;
  Wt.WT.capture(null);
  el._isData.touchmoving = null;
}
);



WT_DECLARE_WT_MEMBER
(AuxWindowTitlebarHandleMove, Wt::JavaScriptFunction, "AuxWindowTitlebarHandleMove",
function( sender, event, id )
{
  var e = event||window.event;
  var el = this.getElement(id);

  if( e.targetTouches.length < 1 )
    return;

  el._isData = el._isData || {};
  if( el._isData.touchmoving == null )
    return;

  var touch = e.targetTouches[0];
  var nowxy = {}, wxy = {};
  nowxy.x = touch.pageX;
  nowxy.y = touch.pageY;
  wxy.x = touch.pageX;
  wxy.y = touch.pageY;
  var wsize = Wt.WT.windowSize();

  var dsx = el._isData.dsx;
  var dsy = el._isData.dsy;

  if(wxy.x > 0 && wxy.x < wsize.x && wxy.y > 0 && wxy.y < wsize.y)
  {
    el.style.left = (Wt.WT.px(el, 'left') + nowxy.x - dsx) + 'px';
    el.style.top = (Wt.WT.px(el, 'top') + nowxy.y - dsy) + 'px';
    el.style.right = "";
    el.style.bottom = "";
    el._isData.dsx = nowxy.x;
    el._isData.dsy = nowxy.y;
    el._isData.centered = false;
  }
}
);

WT_DECLARE_WT_MEMBER
(AuxWindowCollapse, Wt::JavaScriptFunction, "AuxWindowCollapse",
function( sender, event, id, contentid, minId, maxId, fctncall )
{
  var jsthis = document.getElementById(id);
  var jscontent = document.getElementById(contentid);
  if( !jsthis || !jscontent ) return;
  jsthis._isData = jsthis._isData || {};
  var jsfooter = jsthis.querySelector('.modal-footer');
  var dl = jsthis.querySelector(':scope > .dialog-layout');

  var w = jsthis.offsetWidth;
  var h = jsthis.offsetHeight;
  var ch = jscontent.offsetHeight;
  var dlw = dl ? dl.offsetWidth : 0;

  if( jsfooter )
    ch += jsfooter.offsetHeight;

  if( jscontent.offsetParent !== null && jscontent.style.display !== 'none' )
  {
    var oldminheight = jsthis.style.minHeight;
    if( oldminheight )
      jsthis._isData.oldminheight = oldminheight;

    if(jsthis.getAttribute('style').indexOf(' width') !== -1)
      jsthis._isData.oldwidth = jsthis.offsetWidth;
    else
      jsthis._isData.oldwidth = null;

    if(jsthis.getAttribute('style').indexOf(' height') !== -1)
      jsthis._isData.oldheight = jsthis.offsetHeight;
    else
      jsthis._isData.oldheight = null;

    if(dl && dl.getAttribute('style') && dl.getAttribute('style').indexOf(' height') !== -1)
      jsthis._isData.oldlayoutheight = dl.offsetHeight;
    else
      jsthis._isData.oldlayoutheight = null;

    if(dl && dl.getAttribute('style') && dl.getAttribute('style').indexOf(' width') !== -1)
      jsthis._isData.oldlayoutwidth = dl.offsetWidth;
    else
      jsthis._isData.oldlayoutwidth = null;

    jscontent.style.display = 'none';
    if( jsfooter ) jsfooter.style.display = 'none';

    jsthis.style.width = w + 'px';
    if( dl ) dl.style.width = dlw + 'px';

    //jsthis.style.height = Math.max(25, h-ch) + 'px';
    jsthis.style.height = "24px";
    jsthis.style.minHeight = "24px";
    if( dl ){ dl.style.height = "24px"; dl.style.minHeight = "24px"; }
//    if( dl.wtResize )
//      dl.wtResize(dl, dlw, 25);

    if( jsthis.classList.contains('Wt-resizable') )
    {
      jsthis.classList.add('Wt-resizableHolder');
      jsthis.classList.remove('Wt-resizable');
    }
    try{ fctncall(); }catch(e){}
    var minEl = document.getElementById(minId); if( minEl ) minEl.style.display = "none";
    var maxEl = document.getElementById(maxId); if( maxEl ) maxEl.style.display = "";
  }
}
);


WT_DECLARE_WT_MEMBER
(AuxWindowExpand, Wt::JavaScriptFunction, "AuxWindowExpand",
 function( sender, event, id, contentid, minId, maxId, fctncall )
{
  var jsthis = document.getElementById(id);
  if( !jsthis )
    return;

  jsthis._isData = jsthis._isData || {};
  var jscontent = document.getElementById(contentid);
  var dl = jsthis.querySelector(':scope > .dialog-layout');

  if( jscontent && (jscontent.style.display === 'none' || jscontent.offsetParent === null) )
  {
    if( jsthis._isData.oldminheight )
      jsthis.style.minHeight = jsthis._isData.oldminheight;
    else
      jsthis.style.minHeight = "";

    if(jsthis._isData.oldwidth)
      jsthis.style.width = jsthis._isData.oldwidth + 'px';
    else
      jsthis.style.width = "";

    if(jsthis._isData.oldheight)
      jsthis.style.height = jsthis._isData.oldheight + 'px';
    else
      jsthis.style.height = "";

    if( jsthis._isData.oldlayoutheight )
    { if( dl ) dl.style.height = jsthis._isData.oldlayoutheight + 'px'; }
    else
    { if( dl ) dl.style.height = ""; }
    if( jsthis._isData.oldlayoutwidth )
    { if( dl ) dl.style.width = jsthis._isData.oldlayoutwidth + 'px'; }
    else
    { if( dl ) dl.style.width = ""; }

    jsthis._isData.oldwidth = null;
    jsthis._isData.oldheight = null;
    jsthis._isData.oldminheight = null;
    jsthis._isData.oldlayoutwidth = null;
    jsthis._isData.oldlayoutheight = null;

    if( jsthis.classList.contains('Wt-resizableHolder') )
    {
      jsthis.classList.remove('Wt-resizableHolder');
      jsthis.classList.add('Wt-resizable');
    }

    jscontent.style.display = "";
    var mf = jsthis.querySelector(".modal-footer");
    if( mf ) mf.style.display = "";

//    if( dl && dl.wtResize )
//      dl.wtResize(dl, dlw, 25);

    try{ fctncall(); }catch(e){}
    var minEl = document.getElementById(minId); if( minEl ) minEl.style.display = "";
    var maxEl = document.getElementById(maxId); if( maxEl ) maxEl.style.display = "none";
  }
}
);

#if( AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE )
/** This JavaScript function sets a listener (only the first time it is called), that
 when the user stops adjusting the the window size, it will make sure the AuxWindow
 is visible, and if necessary's centered.
 */
WT_DECLARE_WT_MEMBER
(AuxWindowOnDomResize, Wt::JavaScriptFunction, "AuxWindowOnDomResize",
function()
{
  if( window._IS.AWDomResize )
    return;
  window._IS.AWDomResize = true;

  window.addEventListener('resize', function(event) {

    // If we created this resize event, via like Wt.WT.TriggerResizeEvent(), then return,
    //  there is nothing we should do here.
    if( !event.isTrusted )
      return;

    // We will rate-limit this function to a max of every 100ms (a number pulled out of the air)
    let resizeTO = window._IS.ResizeTO;
    if(resizeTO) //If there is a pending event, dont disrupt it, let it keep ticking down
      return;

    const timenow = new Date();
    let lastTO = window._IS.LastResizeTO;
    if( typeof lastTO === 'undefined' ) // No previous resize event
      lastTO = 0; //Would also be reasonable to set lastTO = timenow;

    const dt = Math.max(0, 100 - (timenow - lastTO)); //

    resizeTO = setTimeout(function() {
      window._IS.ResizeTO = null;  //Mark that there is not a pending event
      window._IS.LastResizeTO = new Date(); // Note when the event was processed

      const matches = document.querySelectorAll(".Wt-dialog.AuxWindow");
      matches.forEach( function(dialog){

        // Skip dealing with hidden AuxWindows
        if( (dialog.style.display === 'none') || (dialog.style.visibility === 'hidden') )
          return;

        const ws = Wt.WT.windowSize();
        dialog._isData = dialog._isData || {};
        const was_centered = dialog._isData.centered;
        //console.log( 'Dialog id ', dialog.id, ' centered:', was_centered );

        // Make sure window isnt bigger than our screen.
        //  Note that we are calling the Wt WDialog resize utility so this way the change will
        //  be propagated and scroll bars or whatever is needed will show up
        const h = dialog.offsetHeight;
        const w = dialog.offsetWidth;
        if( (h > ws.y) || (w > ws.x) )
          dialog.wtObj.onresize( Math.min(w, ws.x), Math.min(h, ws.y), true);

        // And finally, center the window, if it was centered before
        if( was_centered )
        {
          Wt.WT.AuxWindowCenter(dialog);
        }else
        {
          // Make sure the window is somewhat visible - and if not, bring it back visible,
          //  before resizing or centering
          //
          // TODO: track if user has it part way off screen, and if not, reposition so whole window is visible
          const pos = {top: dialog.offsetTop, left: dialog.offsetLeft};

          if( pos.top > (ws.y - 25) ) //Title bar is approx 25 px
          {
            dialog.style.top = (ws.y - 25) + 'px';
            dialog.style.bottom = null;
          }

          if( pos.left > (ws.x - 35) ) //Use 35px since the collapse icon could be blocking the first ~20 px of titlebar
          {
            dialog.style.left = (ws.x - 35) + 'px';
            dialog.style.right = null;
          }
        }
      });
    }, dt);
    window._IS.ResizeTO = resizeTO;
  });
}
);//WT_DECLARE_WT_MEMBER( AuxWindowOnDomResize )
#endif //#if( AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE )




AuxWindow::AuxWindow(const Wt::WString& windowTitle, Wt::WFlags<AuxWindowProperties> properties)
  : WDialog(),
  m_auxIsHidden(true),
  m_modalOrig(false),
  m_titleText(NULL),
  m_collapseIcon(NULL),
  m_expandIcon(NULL),
  m_closeIcon(NULL),
  m_collapseSlot(),
  m_expandSlot(),
  m_toggleExpandedStatusSlot(),
  m_closedSignal(),
  m_openedSignal(),
  m_collapsedSignal(),
  m_expandedSignal(),
  m_contentStretcher(NULL),
  m_destructing(false),
  m_escapeIsReject(false),
  m_isPhone(false),
  m_isTablet(false),
  m_isAndroid(false),
  m_footer(NULL)
{
  InterSpec* viewer = InterSpec::instance();

  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowCollapse);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowExpand);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowResizeToFitOnScreen);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowBringToFront);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowTitlebarHandleMove);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowTitlebarTouchEnd);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowTitlebarTouchStart);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowCenter);

#if( AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE )
  addStyleClass("AuxWindow");
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowOnDomResize);
#endif


  const bool isPhone = viewer ? viewer->isPhone() : false;
  const bool isTablet = viewer ? viewer->isTablet() : false;

  const bool isPhoneNotFullScreen = properties.test(AuxWindowProperties::PhoneNotFullScreen);
  const bool isTabletNotFullScreen = properties.test(AuxWindowProperties::TabletNotFullScreen);

  const bool phoneFullScreen = (isPhone && !isPhoneNotFullScreen)
    || (isTablet && !isTabletNotFullScreen && !isPhoneNotFullScreen);

  m_modalOrig = properties.test(AuxWindowProperties::IsModal);
  if (phoneFullScreen)
  {
    m_modalOrig = false;  //Sometimes the Welcome screen modal underlay seems a bit sticky.
    setResizable(false);
    resizeScaledWindow(1.0, 1.0);
    m_isPhone = true; //disables any of future AuxWindow calls to change behavior
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    addStyleClass( "PhoneFullScreenDialog" );
#endif
  }

  setModal(m_modalOrig);

  if (isTablet && !isTabletNotFullScreen && !isPhoneNotFullScreen)
    m_isTablet = true;

  m_isAndroid = (viewer && viewer->isAndroid());

  WContainerWidget* title = titleBar();
  WContainerWidget* content = contents();
  title->clear();
  title->setContentAlignment(AlignmentFlag::Middle);
  title->setPadding(WLength(0));

  content->clear();

  content->addStyleClass("AuxWindow-content");

  if (isPhone || isTablet)
    content->addStyleClass(phoneFullScreen ? "MobileFullScreen" : "MobileModal");

  string fcntcall;
  std::string resources = WApplication::resourcesUrl();

  m_closedSignal = std::make_unique<JSignal<>>(this, "closed", true);
  m_openedSignal = std::make_unique<JSignal<>>(this, "opened", true);

  m_collapsedSignal = std::make_unique<JSignal<>>(this, "collapsed", true);
  m_expandedSignal = std::make_unique<JSignal<>>(this, "expanded", true);

  if (m_isPhone
    || properties.test(AuxWindowProperties::DisableCollapse)
    || ((isPhone || isTablet) && (isTabletNotFullScreen || isPhoneNotFullScreen)))
  {
    m_collapseSlot = std::make_unique<JSlot>("function(){}", this);
    m_expandSlot = std::make_unique<JSlot>("function(){}", this);
    m_toggleExpandedStatusSlot = std::make_unique<JSlot>("function(){}", this);
  }
  else
  {
    auto collapseIconOwner = std::make_unique<WImage>( WLink(resources + "collapse.gif") );
    m_collapseIcon = collapseIconOwner.get();
    m_collapseIcon->setMinimumSize(WLength(16), WLength(16));
    m_collapseIcon->setStyleClass("HeaderIcon");
    m_collapseIcon->addStyleClass("titleVertMiddle");

    auto expandIconOwner = std::make_unique<WImage>( WLink(resources + "expand.gif") );
    m_expandIcon = expandIconOwner.get();
    m_expandIcon->setMinimumSize(WLength(16), WLength(16));
    m_expandIcon->setStyleClass("HeaderIcon");
    m_expandIcon->addStyleClass("titleVertMiddle");
    m_collapseIcon->decorationStyle().setCursor(Cursor::PointingHand);
    m_expandIcon->decorationStyle().setCursor(Cursor::PointingHand);

    title->insertWidget(0, std::move(collapseIconOwner));
    title->insertWidget(0, std::move(expandIconOwner));


    fcntcall = "function(){"
      + m_collapsedSignal->createEventCall(id(), "null", {}) + "}";
    m_collapseSlot = std::make_unique<JSlot>("function(s,e){ Wt.WT.AuxWindowCollapse(s,e,'"
      + id() + "','" + content->id() + "','"
      + m_collapseIcon->id() + "','"
      + m_expandIcon->id() + "',"
      + fcntcall
      + ");}", this);

    m_collapseIcon->clicked().connect(*m_collapseSlot);
    m_collapsedSignal->connect(*m_collapseSlot);

    fcntcall = "function(){"
      + m_expandedSignal->createEventCall(id(), "null", {}) + "}";
    m_expandSlot = std::make_unique<JSlot>("function(s,e){ Wt.WT.AuxWindowExpand(s,e,'"
      + id() + "','" + content->id() + "','"
      + m_collapseIcon->id() + "','"
      + m_expandIcon->id() + "',"
      + fcntcall
      + ");}", this);

    m_expandIcon->clicked().connect(*m_expandSlot);
    m_expandIcon->clicked().connect( this, [this](){ refresh(); } );
    m_expandedSignal->connect(*m_expandSlot);

    stringstream toggleExpandJs;
    toggleExpandJs << "function(s,e)"
      "{"
      "var cel = document.getElementById('" << content->id() << "');"
      "if( cel && (cel.style.display === 'none' || cel.offsetParent === null) )"
      "{"
      "" << m_expandSlot->execJs("s", "e") << ";"
      "}"
      "else"
      "{"
      "" << m_collapseSlot->execJs("s", "e") << ";"
      "}"
      "}";
    m_toggleExpandedStatusSlot = std::make_unique<JSlot>(toggleExpandJs.str(), this);
    title->doubleClicked().connect(*m_toggleExpandedStatusSlot);
    title->doubleClicked().connect( this, [this](){ refresh(); } );
  }//if( m_isPhone ) ' else

  m_titleText = title->addNew<WText>( windowTitle, Wt::TextFormat::XHTML );
  m_titleText->addStyleClass("titleVertMiddle");

  // Centering is done using either centerWindow() or repositionWindow()

  // We wont bring dialog to top, if the user clicked a button on the title (primarily
  //  for mobile, where we might have buttons in the title)
  const string bringToFront = "function(el){"
    "if( el && el.target && !el.target.classList.contains('Wt-btn') && !el.target.classList.contains('FooterHelpBtn'))"
      "Wt.WT.AuxWindowBringToFront('" + id() + "');"
  "}";
  // Use addEventListener instead of mouseWentDown().connect() to avoid overwriting
  // Wt's WDialog onmousedown drag handler on the titlebar element
  doJavaScript("document.getElementById('" + title->id() + "').addEventListener('mousedown'," + bringToFront + ");");

  if( !m_isPhone )
  {
    const bool setCloseable = properties.test(AuxWindowProperties::SetCloseable);
    AuxWindow::setClosable( setCloseable );
    title->touchStarted().connect("function(s,e){ Wt.WT.AuxWindowTitlebarTouchStart(s,e,'" + id() + "'); }");
    title->touchMoved().connect("function(s,e){ Wt.WT.AuxWindowTitlebarHandleMove(s,e,'" + id() + "'); }");
    title->touchEnded().connect("function(s,e){ Wt.WT.AuxWindowTitlebarTouchEnd(s,e,'" + id() + "'); }");
  }

  //so we can just call show from javascript
  WWidget* impw = implementation();
  WTemplate* impl = dynamic_cast<WTemplate*>(impw);
  if (impl)
  {
    impl->setLoadLaterWhenInvisible(false);
    content->setLoadLaterWhenInvisible(false);
  }

  setOffsets(WLength(0, WLength::Unit::Pixel), Wt::Side::Left | Wt::Side::Top);

  //so the signals of all the descendant widgets will be connected
  WDialog::setHidden(false, WAnimation());
  if (m_expandIcon)
    doJavaScript("document.getElementById('" + m_expandIcon->id() + "').style.display='none';");


  if (properties.test(AuxWindowProperties::IsHelpWindow))
  {
    if (!m_isPhone)
    {
      if( m_collapseIcon )
      {
        m_collapseIcon->setImageLink( Wt::WLink( "InterSpec_resources/images/help_minimal.svg" ) );
        m_collapseIcon->setStyleClass( "helpIconTitleBar Wt-icon" );
      }

      if( m_expandIcon )
      {
        m_expandIcon->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg"));
        m_expandIcon->setStyleClass("helpIconTitleBar Wt-icon");
      }
    }//if( !m_isPhone )

    title->setStyleClass("helptitlebar");
  }//helpWindow


  if (!m_isPhone && properties.test(AuxWindowProperties::EnableResize))
  {
    // We have to set minimum size before calling setResizable, or else Wt's Resizable.js functions
    //  will be called first, which will then default to using the initial size as minimum
    //  allowable.
    // So we will set a min size of 200x50 px, which is about the smallest that is useable with
    //  a footer.
    // Client code can set this min width/height to whatever they want later, and it should be fine.
    const WLength mw = minimumWidth();
    const WLength mh = minimumWidth();
    double mw_px = mw.isAuto() ? 0.0 : mw.toPixels();
    double mh_px = mh.isAuto() ? 0.0 : mh.toPixels();

    if ((mw_px < 1.0) || (mh_px < 1.0))
    {
      mw_px = (mw_px < 1.0) ? 200.0 : mw_px;
      mh_px = (mh_px < 1.0) ? 50.0 : mh_px;
      setMinimumSize(mw_px, mh_px);
    }//if( either min height or width is not defined )

    setResizable(true);
  }

  if ((isPhone && isPhoneNotFullScreen)
    || (isTablet && isTabletNotFullScreen))
  {
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }

  show();

#if( AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE )
  doJavaScript("Wt.WT.AuxWindowOnDomResize();");

  // Now make a little code to set the 'centered' variable to false, once user moves or resizes
  //  this AuxWindow.
  if (!m_isPhone)
  {
    // moved().canAutoLearn() and resized().canAutoLearn() are both false, meaning we cant
    //  connect them to JS functions or even JSlot, so we have to go through the C++.

    // Note that this JS will be called a number (maybe even 6) when window is first created,
    // but since the centering code is done via setTimeout(...), the 'centered' data will be
    // end up being true, if the window is centered.
    //
    // However, it is also called when the whole window is resized, so we cant actually rely on this
    //const string js = "$('#" + id() + "').data('centered',false); console.log('AuxWindow Resized or moved');"; //
    //moved().connect( std::bind([js](){ wApp->doJavaScript( js ); }) );
    //resized().connect( std::bind([js](){ wApp->doJavaScript( js ); }) );

    const string js = "document.getElementById('" + title->id() + "').addEventListener('mousedown', function(){"
      "var el=document.getElementById('" + id() + "'); el._isData = el._isData || {}; el._isData.centered = false;"
      "console.log('setting AuxWindow moved');"
      "});";
    doJavaScript(js);

    title->touchMoved().connect("function(){var el=document.getElementById('" + id() + "'); el._isData = el._isData || {}; el._isData.centered = false;}"); //untested
  }
#endif //AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE
}//AuxWindow()


void AuxWindow::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  WDialog::render( flags );
  
  // Sometimes the contents of the AuxWindow will not get the contents laid out correctly, so we'll
  //  just always trigger a rather heavy handed resize event.
  // Definitely a hack
  if (flags.test(Wt::RenderFlag::Full))
  {
    wApp->doJavaScript(wApp->javaScriptClass() + ".TriggerResizeEvent();");

#if(  BUILD_AS_WX_WIDGETS_APP || BUILD_AS_ELECTRON_APP )
    // For wxWidgets, we have a bit of a hack to allow dragging the application window around
    //  by the titlebar - but when the dialog cover is being shown (i.e., a modal WDialog),
    //  the the titlebar cant be clicked.  So here we'll add a listener to the dialog cover
    //  to do the same hack, if the user clicks the titlebar region.
    if (InterSpecApp::isPrimaryWindowInstance())
    {
#if(  BUILD_AS_WX_WIDGETS_APP )
      // TODO: a similar thing for BUILD_AS_ELECTRON_APP, but probably need to add something like ".app-titlebar" to the dialog-cover
      WWidget* coverw = wApp->findWidget("dialog-cover");
      WContainerWidget* dialog_cover = dynamic_cast<WContainerWidget*>(coverw);
      if (dialog_cover && !dialog_cover->mouseWentDown().isConnected())
        dialog_cover->mouseWentDown().connect(wApp->javaScriptClass() + ".MouseDownOnDialogCover");
#endif

      // Raise windows controls (minimize, maximize, close), to above the dialog-cover, so we 
      //  dont block people from being able to close the app.
      wApp->doJavaScript( wApp->javaScriptClass() + ".RaiseWinCntrlsAboveCover();" );
    }//if (InterSpecApp::isPrimaryWindowInstance())
#endif
  }//if (flags & Wt::RenderFlag::RenderFull)
}//render( flags )

void AuxWindow::setResizable(bool resizable)
{
    if (m_isPhone)
      return; //disable running calls after AuxWindow initialized to be mobile
    WDialog::setResizable(resizable);
} //void AuxWindow::setResizable(bool resizable)

void AuxWindow::setModal(bool modal)
{
    m_modalOrig = modal; //keeps the original value
    WDialog::setModal(modal);
} //void AuxWindow::setModal(bool modal)

WContainerWidget* AuxWindow::footer()
{
  if( m_footer )
    return m_footer;
  
  if( m_isPhone )
  {
    WContainerWidget *bar = titleBar();
    bar->setHeight( WLength(35,WLength::Unit::Pixel) );
    bar->setStyleClass( "PhoneAuxWindowHeader" );
    
    //ToDo: for Android should show window title, but not in iOS
    //      (currently hiding everywhere)
    for( size_t i = 0; i < bar->children().size(); ++i )
    {
      WText *caption = static_cast<WText *>( bar->children()[i] );
      if( caption )
        caption->setHidden(true);
    }
    
    auto footerOwner = std::make_unique<WContainerWidget>();
    m_footer = footerOwner.get();
    m_footer->addStyleClass( "PhoneAuxWindowFooter" );

    m_footer->setMargin( 12, Wt::Side::Left );
    m_footer->setMargin( 12, Wt::Side::Right );
    bar->insertWidget( 0, std::move(footerOwner) );
  }else
  {
    m_footer = WDialog::footer();
    m_footer->setHeight(WLength(45,WLength::Unit::Pixel));
    m_footer->setStyleClass("modal-footer");
  }
  
  return m_footer;
}

AuxWindow::~AuxWindow()
{
  WDialog::setModal( false );
  m_destructing = true;
}//~AuxWindow()


void AuxWindow::emitReject()
{
  finished().emit( Wt::DialogCode::Rejected );
}

bool AuxWindow::isPhone() const
{
  return m_isPhone;
}

Wt::WPushButton *AuxWindow::addCloseButtonToFooter( Wt::WString override_txt,
                                                   const bool float_right,
                                                   Wt::WContainerWidget *footerOverride )
{
  auto closeOwner = std::make_unique<WPushButton>();
  WPushButton *close = closeOwner.get();
  
  if( override_txt.key().empty() && (override_txt.toUTF8() == "Close") )
    override_txt = WString::tr("Close");
  
  if( m_isPhone )
  {
    if( override_txt.empty() )
      override_txt = WString::tr("Back");
    
    close->addStyleClass( "MobileBackBtn" );
    if( m_isAndroid )
      close->addStyleClass( "InvertInDark" );
    
    close->setText( override_txt );
  }else
  {
    if( override_txt.empty() )
      override_txt = WString::tr("Close");
    
    close->setText( override_txt );
    if( float_right )
      close->addStyleClass( "DialogClose" );
  }//if( phone ) / else
 
  //Sometimes, the footer may not be footer(), so we allow user to override
  if( !footerOverride )
    footerOverride = footer();
  footerOverride->insertWidget( footerOverride->count(), std::move(closeOwner) );

  return close;
}//Wt::WPushButton *addCloseButtonToFooter()


void AuxWindow::setClosable( bool closable )
{
  //prevent mobile from having icon?
  if( m_isPhone )
     return;
    
  if( closable )
  {
    if( m_closeIcon )
      return;
    m_closeIcon = titleBar()->addNew<WText>();
    WApplication::instance()->theme()->apply(this, m_closeIcon,
                                              WidgetThemeRole::DialogCloseIcon);
    // Expand first (in case collapsed), then hide
    m_closeIcon->clicked().connect( *m_expandSlot );
    m_closeIcon->clicked().connect( this, [this](){ hide(); } );
  }else
  {
    if( !m_closeIcon )
      return;
    auto owned = m_closeIcon->removeFromParent();
    m_closeIcon = nullptr;
  }//if( closable ) / else
}//void setClosable( bool closable );


bool AuxWindow::isVisible() const
{
  return !m_auxIsHidden;
}//bool isVisible() const


bool AuxWindow::isHidden() const
{
  return m_auxIsHidden;
}//bool isHidden() const



WContainerWidget *AuxWindow::contents() const
{
  assert( !m_contentStretcher );
  return WDialog::contents();
}

WGridLayout *AuxWindow::stretcher()
{
  if( !m_contentStretcher )
  {
    WContainerWidget *container = WDialog::contents();
    if( container->count() )
    {
      cerr << endl << endl << endl << "AuxWindow::stretcher():"
           << " Warning, the window has children, but I'm adding a layout anyway"
           << endl;
      container->clear();
    }//if( container->count() )

    m_contentStretcher = container->setLayout( std::make_unique<WGridLayout>() );
    container->setOverflow( Wt::Overflow::Hidden );
  }//if( !m_contentStretcher )

  return m_contentStretcher;
}//WGridLayout *AuxWindow::stretcher()


void AuxWindow::resizeToFitOnScreen()
{
  if( m_isPhone )
    return; //disable running calls after AuxWindow initialized to be mobile
  
  //We will put in a setTimeout to call the resize function to first give Wts
  //  layout a chance to shrink things into appropriate sizes.
  //TODO: A timeout of 0 ms seems to work, but using 50 JIC - revisit.
  const string js = "setTimeout( function(){ "
    "Wt.WT.AuxWindowResizeToFitOnScreen('" + id() + "');"
  "},50);";
  
  doJavaScript( js );
}//void AuxWindow::resizeToFitOnScreen()


void AuxWindow::centerWindow()
{
  if( m_isPhone )
      return; //disable running calls after AuxWindow initialized to be mobile

  doJavaScript( "Wt.WT.AuxWindowCenter(" + this->jsRef() + ");" );
}

void AuxWindow::centerWindowHeavyHanded()
{
  // Heavy-handed centering is no longer needed with Wt 4
  centerWindow();
}//void centerWindowHeavyHanded()


std::string AuxWindow::repositionWindowJs( int x, int y )
{
  if( m_isPhone )
    return ""; //disable running calls after AuxWindow initialized to be mobile
  
  return string("var el = ") + this->jsRef() + ";"
  "var olddispl = el.style.display; el.style.display='';"
  + ((y==-32768) ? string() : "el.style.top = '"  + std::to_string(y) + "px';")
  + ((x==-32768) ? string() :  "el.style.left = '" + std::to_string(x) + "px';")
  + "el.style.display = olddispl;"
  "el._isData = el._isData || {}; el._isData.centered = false;";
}//std::string repositionWindowJs( int x, int y )


void AuxWindow::repositionWindow( int x, int y )
{
  const string js = repositionWindowJs( x, y );
  if( !js.empty() )
    doJavaScript( js );
} // void AuxWindow::repositionWindow( int x, int y )


void AuxWindow::resizeWindow( int width, int height )
{
  if( m_isPhone )
    return; //disable running calls after AuxWindow initialized to be mobile
    
  WStringStream sizeJs;
  sizeJs << "var el = " << this->jsRef() << ";"
         << "var olddispl = el.style.display; el.style.display='';"
         << "var ws = " << wApp->javaScriptClass() << ".WT.windowSize();";
  if( width > 0 )
    sizeJs << "el.style.width = Math.round(Math.min(" << width << ",0.95*ws.x)) + 'px';";
  if( height > 0 )
    sizeJs << "el.style.height = Math.round(Math.min(" << height << ",0.95*ws.y)) + 'px';";
  sizeJs << "el.style.display = olddispl;";
  doJavaScript( sizeJs.str() );
} // void AuxWindow::resizeWindow( int width, int height )


std::string AuxWindow::resizeScaledWindowJs( double xRatio, double yRatio ) const
{
  if( m_isPhone )
    return ""; //disable running calls after AuxWindow initialized to be mobile
  
  if( xRatio <= 0.0 && yRatio <= 0.0 )
    return "";
  
  stringstream sizeJs;
  sizeJs << "var el = " << this->jsRef() << ";"
         << "var olddispl = el.style.display; el.style.display='';";
  sizeJs << "var ws = " << wApp->javaScriptClass() << ".WT.windowSize();";

  if( xRatio > 0.0 )
    sizeJs << "if(ws.x>2) el.style.width = Math.round(" << xRatio << "*ws.x - 2) + 'px';";
  if( yRatio > 0.0 )
    sizeJs << "if(ws.y>2) el.style.height = Math.round(" << yRatio << "*ws.y - 2) + 'px';";
  sizeJs << "el.style.display = olddispl;";

  return sizeJs.str();
}//std::string resizeScaledWindowJs( double xRatio, double yRatio ) const

void AuxWindow::resizeScaledWindow( double xRatio, double yRatio )
{
  const string js = resizeScaledWindowJs(xRatio,yRatio);
 
  if( !js.empty() )
    doJavaScript( js );
}//void resizeScaledWindow( double xRatio, double yRatio )


void AuxWindow::deleteAuxWindow( AuxWindow *window )
{
  if( !window )
    return;

  // In Wt 4, dialog lifecycle must be managed carefully:
  //  1. Dismiss modal state so the dialog cover is properly updated
  //  2. Hide the dialog through the proper WDialog path (which calls popDialog on the cover)
  //  3. Remove from parent to destroy
  //
  // Note: this function is often called during a `finished()` signal handler of the window
  //  being deleted (e.g., from the Close button). Destroying the window during its own signal
  //  is generally safe in Wt 4 as long as the dialog cover is properly updated first.

  if( window->m_modalOrig )
    window->WDialog::setModal( false );

  if( !window->isHidden() )
  {
    // Use the AuxWindow::hide() path which calls WDialog::setHidden(true) to properly
    //  trigger popDialog() on the DialogCover.
    window->m_destructing = true;  // prevent emitReject() from firing finished() again
    window->WDialog::setHidden( true, WAnimation() );
  }

  // Use wApp->removeChild (not removeFromParent) because the factory
  //  gives ownership to WApplication via addChild. removeFromParent goes
  //  through domRoot_ which does not hold ownership for global widgets,
  //  and would leak the object.
  Wt::WApplication::instance()->removeChild( window );
}//void deleteAuxWindow( AuxWindow *window )


void AuxWindow::deleteSelf()
{
  AuxWindow::deleteAuxWindow( this );
}


void AuxWindow::rejectWhenEscapePressed( bool enable )
{
  if( enable == m_escapeIsReject )
    return;
  
  m_escapeIsReject = enable;
  if( m_escapeIsReject )
  {
    WWidget *impw = implementation();
    WTemplate *impl = dynamic_cast<WTemplate *>( impw );
    if( impl )
      m_escapeConnection1 = impl->escapePressed().connect(this, &AuxWindow::hide);
    
    WApplication *app = WApplication::instance();
    if( app )
      m_escapeConnection2 = app->globalEscapePressed()
                                              .connect(this, &AuxWindow::hide);
  }else
  {
    m_escapeConnection1.disconnect();
    m_escapeConnection2.disconnect();
  }
}//void rejectWhenEscapePressed()


void AuxWindow::show()
{
    doJavaScript( "Wt.WT.AuxWindowBringToFront('" + id() + "');" );

    //Check original modal value, as modal doesn't work well with hide/show
    if( m_modalOrig )
      WDialog::setModal( true );
    setHidden( false, WAnimation() );
}//void show()


void AuxWindow::hide()
{
    //Check original modal value, as modal doesn't work well with hide/show
    if( m_modalOrig )
      WDialog::setModal( false );

    setHidden( true, WAnimation() );
}//void hide()


void AuxWindow::expand() //executes jsExpand()
{
  m_expandSlot->exec();
}//void expand()


void AuxWindow::setWindowTitle( const Wt::WString& windowTitle )
{
  if( m_titleText )
    m_titleText->setText( windowTitle );
}//void setWindowTitle( const Wt::WString& windowTitle )


void AuxWindow::setMaximumSize( const Wt::WLength &width, const Wt::WLength &height )
{
  if( m_isPhone )
    return;
  WDialog::setMaximumSize( width, height );
}

void AuxWindow::setHidden( bool hide, const WAnimation &/*animation*/ ) //TODO: incorportate WAnimation
{
  if( m_destructing )
  {
    WDialog::setHidden( hide );
    if( m_footer )
      m_footer->setHidden(hide);
    return;
  }//if( m_destructing )

  const bool changingState = (m_auxIsHidden != hide);

  m_auxIsHidden = hide;

  if( m_footer )
    m_footer->setHidden(hide); //Is this really necassary?

  if( hide )
  {
    if( changingState && m_escapeIsReject )
    {
      m_escapeConnection1.disconnect();
      m_escapeConnection2.disconnect();
    }//if( m_escapeIsReject )

    WDialog::setHidden( true, WAnimation() );

    emitReject();
  }else
  {
    WDialog::setHidden( false, WAnimation() );

    doJavaScript( "Wt.WT.AuxWindowBringToFront('" + id() + "');" );

    if( changingState && m_escapeIsReject )
    {
      WWidget *impw = implementation();
      WTemplate *impl = dynamic_cast<WTemplate *>( impw );
      if( impl )
        m_escapeConnection1 = impl->escapePressed().connect(this, &AuxWindow::hide);

      WApplication *app = WApplication::instance();
      if( app )
        m_escapeConnection2 = app->globalEscapePressed().connect(this, &AuxWindow::hide);
    }//if( m_escapeIsReject )
  }//if( hide ) / else
}//void setHidden( bool hide, const Wt::WAnimation &animation )


void AuxWindow::disableCollapse()
{
  if( m_collapseIcon )
    m_collapseIcon->hide();
  if( m_expandIcon )
    m_expandIcon->hide();
  if( !!m_collapseSlot )
    m_collapseSlot->setJavaScript("function(){}");
  if( !m_isPhone )
  {
    titleBar()->setPadding(WLength(5), Wt::Side::Left);
    titleBar()->setPadding(WLength(2), Wt::Side::Top);
  }
}



JSlot &AuxWindow::jsExpand()
{
  return *m_expandSlot;
}


JSlot &AuxWindow::jsCollapse()
{
  return *m_collapseSlot;
}


JSlot &AuxWindow::jsToggleExpandState()
{
  return *m_toggleExpandedStatusSlot;
}


JSignal<> &AuxWindow::opened()
{
  return *m_openedSignal;
}


JSignal<> &AuxWindow::collapsed()
{
  return *m_collapsedSignal;
}


JSignal<> &AuxWindow::expanded()
{
  return *m_expandedSignal;
}


void AuxWindow::addHelpInFooter( WContainerWidget *footer, std::string page )
{
  InterSpec *viewer = InterSpec::instance();
  
  Wt::WImage *image = footer->addNew<Wt::WImage>( Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
  image->setStyleClass("Wt-icon FooterHelpBtn");
  image->setFloatSide( (viewer && viewer->isMobile()) ? Wt::Side::Right : Wt::Side::Left );

  image->setAlternateText("Help");
  image->clicked().connect( image, [page](){ HelpSystem::createHelpWindow( page ); } );
}
