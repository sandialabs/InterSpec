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
#include <Wt/WText>
#include <Wt/WTheme>
#include <Wt/WPanel>
#include <Wt/WImage>
#include <Wt/WIconPair>
#include <Wt/WTemplate>
#include <Wt/WBoxLayout>
#include <Wt/WBorderLayout>
#include <Wt/WPushButton>
#include <Wt/WJavaScript>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>

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
    const el = $("#" + id);
    const cel = $("#" + id + " .body.AuxWindow-content");
    if( !el || !cel )
    {
      console.log( 'Error in resizeToFitOnScreen' );
      return;
    }
    
    const eel = el.get(0);
    const ecel = cel.get(0);
    const ws = Wt.WT.windowSize();
    
    let h = el.height();
    let w = el.width();
    const do_resize = ((h > ws.y) || (w > ws.x));
    if(h > ws.y)
    {
      h = ws.y - 10;
      eel.style.top = '5px';
      eel.style.bottom = "";
      ecel.style.overflowY = "auto";
    }
    
    if(w > ws.x)
    {
      w = ws.x - 10;
      eel.style.left = '5px';
      eel.style.right = "";
      ecel.style.overflowX = "auto";
    }
    
    if( do_resize )
    {
      const was_centered = el.data('centered');
      eel.wtObj.onresize(w, h, true);
      
      if( was_centered )
        Wt.WT.AuxWindowCenter( el.get(0) );
      
      el.data('centered', was_centered);
    }
  }catch(err){
    console.error( 'AuxWindowResizeToFitOnScreen error:', err );
  }
}
 );//WT_DECLARE_WT_MEMBER( AuxWindowResizeToFitOnScreen )


#if( !USE_NEW_AUXWINDOW_ISH )
WT_DECLARE_WT_MEMBER
(AuxWindowCenter, Wt::JavaScriptFunction, "AuxWindowCenter",
function( el )
{
  try
  {
    const jsthis = $(el);
    const olddispl = el.style.display;
    el.style.display = "";
    const ws = Wt.WT.windowSize();
    const hh = Math.max(0,Math.round( ( ws.y - jsthis.height() ) / 2 ));
    const hw = Math.max(0,Math.round( ( ws.x - jsthis.width() ) / 2 ));
    
    el.style.top = hh+'px';
    el.style.left = hw+'px';
    el.style.display = olddispl;
    $(el).data('centered', true);
  }catch(err)
  {
    console.log( "Failed in AuxWindowCenter: " + err );
  }
}
);
#endif //#if( !USE_NEW_AUXWINDOW_ISH )




WT_DECLARE_WT_MEMBER
(AuxWindowBringToFront, Wt::JavaScriptFunction, "AuxWindowBringToFront",
 function( id )
{
  var maxz = 0;
  $('.Wt-dialog,.MobileMenuButton').each( function(index, value){
    let z = $(value).css('z-index');
    //if( $(value).is(":visible") ) //seems to not be reliable, so leaving commented out.
    maxz = Math.max(maxz,z);
  } );
  if( maxz > $('#'+id).css('z-index') )
    $('#'+id).css('z-index',maxz+1);
  $('.window-controls-container').css('z-index', maxz+1); //for wxWidgets and Electron builds
  $('.suggestion').css('z-index',maxz+2); //This is for InterSpec::m_shieldingSuggestion
}
);


WT_DECLARE_WT_MEMBER
(AuxWindowTitlebarTouchStart, Wt::JavaScriptFunction, "AuxWindowTitlebarTouchStart",
function( sender, event, id )
{
  var e = event||window.event;
  var el = this.getElement(id);
  var titlebar = $(el).find(".titlebar").first().get(0);
  
  if($(el).data('touchmoving') == null)
  {
    Wt.WT.capture(titlebar);
    var touch = event.targetTouches[0];
    $(el).data('dsx', touch.pageX);
    $(el).data('dsy', touch.pageY);
    $(el).data('touchmoving',true);
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
  var titlebar = $(el).find(".titlebar").first().get(0);
  $(el).data('dsx', null);
  $(el).data('dsy', null);
  Wt.WT.capture(null);
  $(el).data('touchmoving',null);
}
);



WT_DECLARE_WT_MEMBER
(AuxWindowTitlebarHandleMove, Wt::JavaScriptFunction, "AuxWindowTitlebarHandleMove",
function( sender, event, id )
{
  var e = event||window.event;
  var el = this.getElement(id);
  var titlebar = $(el).find(".titlebar").first().get(0);
  
  if( e.targetTouches.length < 1 )
    return;
  
  if( $(el).data('touchmoving')==null )
    return;
  
  var touch = e.targetTouches[0];
  var nowxy = {}, wxy = {};
  nowxy.x = touch.pageX;
  nowxy.y = touch.pageY;
  wxy.x = touch.pageX;
  wxy.y = touch.pageY;
  var wsize = Wt.WT.windowSize();
  
  var dsx = $(el).data('dsx');
  var dsy = $(el).data('dsy');

  if(wxy.x > 0 && wxy.x < wsize.x && wxy.y > 0 && wxy.y < wsize.y)
  {
    el.style.left = (Wt.WT.px(el, 'left') + nowxy.x - dsx) + 'px';
    el.style.top = (Wt.WT.px(el, 'top') + nowxy.y - dsy) + 'px';
    el.style.right = "";
    el.style.bottom = "";
    $(el).data('dsx', nowxy.x);
    $(el).data('dsy', nowxy.y);
    $(el).data('centered', false);
  }
}
);

WT_DECLARE_WT_MEMBER
(AuxWindowCollapse, Wt::JavaScriptFunction, "AuxWindowCollapse",
function( sender, event, id, contentid, minId, maxId, fctncall )
{
  var jsthis = $('#'+id);
  var jscontent = $('#'+contentid);
  var jsfooter = jsthis.find('.modal-footer');
  var dl = jsthis.children('.dialog-layout');
  
  var w = jsthis.width();
  var h = jsthis.height();
  var ch = jscontent.height();
  var dlw = dl.width();
  
  if( jsfooter.length !== 0 )
    ch += jsfooter.height();
  
  if( jscontent.is(':visible') )
  {
    var oldminheight = jsthis.css('min-height');
    if( oldminheight )
      jsthis.data('oldminheight',oldminheight);
    
    if(jsthis.attr('style').indexOf(' width') !== -1)
      jsthis.data('oldwidth',jsthis.width());
    else
      jsthis.data('oldwidth',null);
    
    if(jsthis.attr('style').indexOf(' height') !== -1)
      jsthis.data('oldheight',jsthis.height());
    else
      jsthis.data('oldheight',null);
    
    if(dl.attr('style').indexOf(' height') !== -1)
      jsthis.data('oldlayoutheight',dl.height());
    else
      jsthis.data('oldlayoutheight',null);

    if(dl.attr('style').indexOf(' width') !== -1)
      jsthis.data('oldlayoutwidth',dl.width());
    else
      jsthis.data('oldlayoutwidth',null);
    
    jscontent.hide();
    jsfooter.hide();

    jsthis.width(w);
    dl.width(dlw);
    
    //jsthis.height(Math.max(25, h-ch));
    jsthis.css( 'height', "24px" );
    jsthis.css ('min-height', "24px" );
    dl.css( 'height', "24px" );
    dl.css( 'min-height', "24px" );
//    if( dl.get(0).wtResize )
//      dl.get(0).wtResize(dl.get(0), dlw, 25);
    
    if( jsthis.hasClass('Wt-resizable'))
    {
      jsthis.addClass('Wt-resizableHolder');
      jsthis.removeClass('Wt-resizable');
    }
    try{ fctncall(); }catch(e){}
    $('#'+minId).hide();
    $('#'+maxId).show();
  }
}
);


WT_DECLARE_WT_MEMBER
(AuxWindowExpand, Wt::JavaScriptFunction, "AuxWindowExpand",
 function( sender, event, id, contentid, minId, maxId, fctncall )
{
  var jsthis = $('#'+id);
  if( jsthis.length === 0 )
    return;
  
  var jscontent = $('#'+contentid);
  var dl = jsthis.children('.dialog-layout');
  
  if( jscontent.is(':hidden') )
  {
    if( jsthis.data('oldminheight') )
      jsthis.css('min-height', jsthis.data('oldminheight') );
    else
      jsthis.css('min-height', "" );
    
    if(jsthis.data('oldwidth'))
      jsthis.css( 'width', jsthis.data('oldwidth') );
    else
      jsthis.css( 'width', "" );
    
    if(jsthis.data('oldheight'))
      jsthis.css( 'height', jsthis.data('oldheight') );
    else
      jsthis.css( 'height', "" );
    
    if( jsthis.data('oldlayoutheight') )
      dl.css( 'height', jsthis.data('oldlayoutheight') );
    else
      dl.css( 'height', "" );
    if( jsthis.data('oldlayoutwidth') )
      dl.css( 'width', jsthis.data('oldlayoutwidth') );
    else
      dl.css( 'width', "" );
    
    jsthis.data('oldwidth',null);
    jsthis.data('oldheight',null);
    jsthis.data('oldminheight',null);
    jsthis.data('oldlayoutwidth',null);
    jsthis.data('oldlayoutheight',null);
    
    if( jsthis.hasClass('Wt-resizableHolder'))
    {
      jsthis.removeClass('Wt-resizableHolder');
      jsthis.addClass('Wt-resizable');
    }
    
    jscontent.show();
    jsthis.find('.modal-footer').show();
    
//    if( dl.get(0).wtResize )
//      dl.get(0).wtResize(dl.get(0), dlw, 25);
    
    try{ fctncall(); }catch(e){}
    $('#'+minId).show();
    $('#'+maxId).hide();
  }
}
);

WT_DECLARE_WT_MEMBER
(AuxWindowHide, Wt::JavaScriptFunction, "AuxWindowHide",
 function( sender, event, id, fcntcall )
{
  var jsthis = $('#'+id);
  if( jsthis.length === 0 )
    return;
  
  var currleft = jsthis.offset().left;
  var currtop = jsthis.offset().top;
  
  if( (typeof currleft !== null) && currleft > -999 )
  {
    try{if(!event || !event.quietly) fcntcall(); }catch(e){}
  }
  
  var prevleft = jsthis.data('prevleft');
  var prevtop = jsthis.data('prevtop');
  if( !prevleft && prevleft!==0 )
  {
    jsthis.data('prevleft', currleft );
    jsthis.data('prevtop', currtop );
  }
  
  jsthis.offset( { top: -1999, left: -1999 } );
  
  jsthis.data('notshown',true);
  
//  if(!jsthis.is(':hidden') ){ try{if(!event || !event.quietly) fcntcall(); }catch(e){} }
//  if( el )
//  {
//    var prevleft = jsthis.offset().left;
//    var prevtop = jsthis.offset().top;
//    jsthis.offset( { top: -999, left: -999 } );
//    if( event.delay )
//    {
//      var prevleft = jsthis.offset().left;
//      var prevtop = jsthis.offset().top;
//      jsthis.offset( { top: -999, left: -999 } );
      //XXX - hack, in order for things to render okay when we want to show the
      //      File Manager, we have to let everythign render okay, before hiding
      //      it the first time, but inorder for it to not pop-up and be visible
      //      we first have to move it off the screen, and then move it back.
//      setTimeout(
//        function()
//        {
//          jsthis.offset( { left: prevleft, top: prevtop } );
//          el.style.display = 'none';
//        }, 50 );
//    }else
//    {
//      el.style.display = "none";
//    }
//  }
}
 );


WT_DECLARE_WT_MEMBER
(AuxWindowShow, Wt::JavaScriptFunction, "AuxWindowShow",
 function( sender, event, id, fcntcall )
{
  var jsthis = $('#'+id);
  if( jsthis.length === 0 )
    return;
  
  var prevleft = jsthis.data('prevleft');
  var prevtop  = jsthis.data('prevtop');
  var currleft = jsthis.offset().left;
  var currtop  = jsthis.offset().top;
  
  if( prevleft || prevleft===0 )
  {
    try{if(!event || !event.quietly) fcntcall(); }catch(e){}
    jsthis.offset( { top: prevtop, left: prevleft } );
  }//else jsthis.offset( { top: null, left: null } );
  
  jsthis.data('prevleft', null );
  jsthis.data('prevtop', null );

//  var jsthis = $('#'+id);
//  if( jsthis.is(':hidden') )
//  {
//    try{if(!event || !event.quietly) fcntcall();}catch(e){}
//  }
  
  var el = jsthis.get(0);
  
  //Make sure the user can at least close the window
  if( jsthis.height() > $(window).height() )
  {
    el.style.top = '0px';
//    setTimeout( function(){
//      el.style.top = '5px';
//      jsthis.css("min-height",null);
//      jsthis.css("height",Math.round(0.95*$(window).height()));
//      jsthis.css("overflow-y", "auto");
//    }, 100 );
  }
  
  if( jsthis.width() > $(window).width() )
  {
    el.style.left = '0px';
//    setTimeout( function(){
//       el.style.left = '5px';
//       jsthis.css("overflow-x", "auto");
//       jsthis.css("min-width",null);
//       jsthis.css("width",Math.round(0.95*$(window).width()));
//      }, 100 );
  }
  
  
  if( el )
  {
    //see note below for the reason this js code block is necessary
    //This next line is a hack to force things to be re-flowed to keep things from being 'sticky'
    $('<style></style>').appendTo($(document.body)).remove();
    if( jsthis.data('notshown') )
    {
      jsthis.data('notshown',false);
      var ws = Wt.WT.windowSize();
      el.style.marginLeft = '0px';
      el.style.marginTop = '0px';
      el.style.bottom = null; el.style.right=null;
    }
  }
  
  Wt.WT.AuxWindowBringToFront(id);
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
  if( $('.Wt-domRoot').data('AWDomResize') )
    return;
  $('.Wt-domRoot').data('AWDomResize', true);
  
  window.addEventListener('resize', function(event) {
    
    // If we created this resize event, via like Wt.WT.TriggerResizeEvent(), then return,
    //  there is nothing we should do here.
    if( !event.isTrusted )
      return;
    
    // We will rate-limit this function to a max of every 100ms (a number pulled out of the air)
    let resizeTO = $('.Wt-domRoot').data('ResizeTO');
    if(resizeTO) //If there is a pending event, dont disrupt it, let it keep ticking down
      return;
    
    const timenow = new Date();
    let lastTO = $('.Wt-domRoot').data('LastResizeTO');
    if( typeof lastTO === 'undefined' ) // No previous resize event
      lastTO = 0; //Would also be reasonable to set lastTO = timenow;
    
    const dt = Math.max(0, 100 - (timenow - lastTO)); //
    
    resizeTO = setTimeout(function() {
      $('.Wt-domRoot').data('ResizeTO',null);  //Mark that there is not a pending event
      $('.Wt-domRoot').data('LastResizeTO', new Date()); // Note when the event was processed
      
      const matches = document.querySelectorAll(".Wt-dialog.AuxWindow");
      matches.forEach( function(dialog){
        
        // Skip dealing with hidden AuxWindows
        if( (dialog.style.display === 'none') || (dialog.style.visibility === 'hidden') )
          return;
        
        const ws = Wt.WT.windowSize();
        const was_centered = $(dialog).data('centered');
        //console.log( 'Dialog id ', dialog.id, ' centered:', was_centered );
        
        // Make sure window isnt bigger than our screen.
        //  Note that we are calling the Wt WDialog resize utility so this way the change will
        //  be propagated and scroll bars or whatever is needed will show up
        const h = $(dialog).height();
        const w = $(dialog).width();
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
          const pos = $(dialog).position();
          
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
    $('.Wt-domRoot').data('ResizeTO',resizeTO);
  });
}
);//WT_DECLARE_WT_MEMBER( AuxWindowOnDomResize )
#endif //#if( AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE )


#if( USE_NEW_AUXWINDOW_ISH )
WT_DECLARE_WT_MEMBER
(AuxWindowScaleResize, Wt::JavaScriptFunction, "AuxWindowScaleResize",
 function( el, e )
{
  try
  {
    var xRatio = e.x;
    var yRatio = e.y;
    var jsthis = $(el);
    var olddispl = el.style.display;
    el.style.display="";
    var ws = Wt.WT.windowSize();
    jsthis.width( xRatio*ws.x );
    jsthis.height( yRatio*ws.y );
    /*  jsthis.data('notshown',true);*/
    el.style.display = olddispl;
  }catch(err)
  {
    console.log( "Failed AuxWindowScaleResize: " + err );
  }
}
 );



WT_DECLARE_WT_MEMBER
(AuxWindowResize, Wt::JavaScriptFunction, "AuxWindowResize",
 function( el, e )
{
  try
  {
    var jsthis = $(el);
    var olddispl = el.style.display;
    el.style.display="";
    var ws = Wt.WT.windowSize();
    jsthis.width( Math.min(e.x,0.95*ws.x) );
    jsthis.height( Math.min(e.y,0.95*ws.y) );
    jsthis.data('notshown',true);
    el.style.display = olddispl;
  }catch(err)
  {
    console.log( "Failed AuxWindowResize: " + err );
  }
}
 );


WT_DECLARE_WT_MEMBER
(AuxWindowCenter, Wt::JavaScriptFunction, "AuxWindowCenter",
 function( el, e )
{
  try
  {
    var jsthis = $(el);
    var olddispl = el.style.display;
    el.style.display="";
    var ws = Wt.WT.windowSize();
    var hh = Math.round( ( ws.y - jsthis.height() ) / 2 );
    var hw = Math.round( ( ws.x - jsthis.width() ) / 2 );
    if( hh < 0 )
      hh = 0;
    if( hw < 0 )
      hw = 0;
    if( jsthis.data('prevtop') )
    {
      jsthis.data( 'prevtop', hh );
      jsthis.data( 'prevleft', hw );
    }else
    {
      el.style.top = hh+'px';
      el.style.left = hw+'px';
    }
    el.style.display = olddispl;
  }catch(err)
  {
    console.log( "Failed in AuxWindowCenter: " + err );
  }
}
 );

WT_DECLARE_WT_MEMBER
(AuxWindowReposition, Wt::JavaScriptFunction, "AuxWindowReposition",
 function( el, e )
{
  try
  {
    var jsthis = $(el);
    var olddispl = el.style.display;
    el.style.display="";
    if( jsthis.data('prevtop') )
    {
      jsthis.data( 'prevtop', e.top );
      jsthis.data( 'prevleft', e.left );
    }else
    {
      el.style.top = e.top + "px";
      el.style.left = e.left + "px";
    }
    el.style.display = olddispl;
  }catch(err)
  {
    console.log( "Failed in AuxWindowReposition: " + err );
  }
}
 );
#endif


AuxWindow::AuxWindow(const Wt::WString& windowTitle, Wt::WFlags<AuxWindowProperties> properties)
  : WDialog(InterSpec::instance()),
  m_auxIsHidden(true),
  m_modalOrig(false),
  m_titleText(NULL),
  m_collapseIcon(NULL),
  m_expandIcon(NULL),
  m_closeIcon(NULL),
  m_collapseSlot(),
  m_expandSlot(),
  m_showSlot(),
  m_hideSlot(),
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
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowShow);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowHide);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowCenter);

#if( AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE )
  addStyleClass("AuxWindow");
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowOnDomResize);
#endif

#if( USE_NEW_AUXWINDOW_ISH )
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowScaleResize);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowResize);
  LOAD_JAVASCRIPT(wApp, "AuxWindow.cpp", "AuxWindow", wtjsAuxWindowReposition);
#endif


  const bool isPhone = viewer ? viewer->isPhone() : false;
  const bool isTablet = viewer ? viewer->isTablet() : false;

  const bool isPhoneNotFullScreen = properties.testFlag(AuxWindowProperties::PhoneNotFullScreen);
  const bool isTabletNotFullScreen = properties.testFlag(AuxWindowProperties::TabletNotFullScreen);

  const bool phoneFullScreen = (isPhone && !isPhoneNotFullScreen)
    || (isTablet && !isTabletNotFullScreen && !isPhoneNotFullScreen);

  m_modalOrig = properties.testFlag(AuxWindowProperties::IsModal);
  if (phoneFullScreen)
  {
    m_modalOrig = false;  //Sometimes the Welcome screen modal underlay seems a bit sticky.
    setResizable(false);
    resizeScaledWindow(1.0, 1.0);
    m_isPhone = true; //disables any of future AuxWindow calls to change behavior
  }

  setModal(m_modalOrig);

  if (isTablet && !isTabletNotFullScreen && !isPhoneNotFullScreen)
    m_isTablet = true;

  m_isAndroid = (viewer && viewer->isAndroid());


  //  setJavaScriptMember( "wtResize", "function(ignored,w,h){console.log('In My wtResize');}" );

  WContainerWidget* title = titleBar();
  WContainerWidget* content = contents();
  title->clear();
  title->setContentAlignment(AlignMiddle);
  title->setPadding(WLength(0));

  content->clear();

  content->addStyleClass("AuxWindow-content");

  if (isPhone || isTablet)
    content->addStyleClass(phoneFullScreen ? "MobileFullScreen" : "MobileModal");

  string fcntcall;
  std::string resources = WApplication::resourcesUrl();

  m_closedSignal.reset(new JSignal<void>(this, "closed", true));
  m_openedSignal.reset(new JSignal<void>(this, "opened", true));

  m_collapsedSignal.reset(new JSignal<void>(this, "collapsed", true));
  m_expandedSignal.reset(new JSignal<void>(this, "expanded", true));

  if (m_isPhone
    || properties.testFlag(AuxWindowProperties::DisableCollapse)
    || ((isPhone || isTablet) && (isTabletNotFullScreen || isPhoneNotFullScreen)))
  {
    m_collapseSlot.reset(new JSlot("function(){}", this));
    m_expandSlot.reset(new JSlot("function(){}", this));
    m_toggleExpandedStatusSlot.reset(new JSlot("function(){}", this));
  }
  else
  {
    m_collapseIcon = new WImage(resources + "collapse.gif");
    m_collapseIcon->setMinimumSize(WLength(16), WLength(16));
    m_collapseIcon->setStyleClass("HeaderIcon");
    m_collapseIcon->addStyleClass("titleVertMiddle");

    m_expandIcon = new WImage(resources + "expand.gif");
    m_expandIcon->setMinimumSize(WLength(16), WLength(16));
    m_expandIcon->setStyleClass("HeaderIcon");
    m_expandIcon->addStyleClass("titleVertMiddle");
    m_collapseIcon->decorationStyle().setCursor(PointingHandCursor);
    m_expandIcon->decorationStyle().setCursor(PointingHandCursor);

    title->insertWidget(0, m_collapseIcon);
    title->insertWidget(0, m_expandIcon);

    const string jsthis = "$('#" + id() + "')";
    const string jscontent = "$('#" + content->id() + "')";

    fcntcall = "function(){"
      + m_collapsedSignal->createEventCall(id(), "null") + "}";
    m_collapseSlot.reset(new JSlot("function(s,e){ Wt.WT.AuxWindowCollapse(s,e,'"
      + id() + "','" + content->id() + "','"
      + m_collapseIcon->id() + "','"
      + m_expandIcon->id() + "',"
      + fcntcall
      + ");}", this));

    m_collapseIcon->clicked().connect(*m_collapseSlot);
    m_collapsedSignal->connect(*m_collapseSlot);

    fcntcall = "function(){"
      + m_expandedSignal->createEventCall(id(), "null") + "}";
    m_expandSlot.reset(new JSlot("function(s,e){ Wt.WT.AuxWindowExpand(s,e,'"
      + id() + "','" + content->id() + "','"
      + m_collapseIcon->id() + "','"
      + m_expandIcon->id() + "',"
      + fcntcall
      + ");}", this));

    m_expandIcon->clicked().connect(*m_expandSlot);
    m_expandIcon->clicked().connect(boost::bind(&AuxWindow::refresh, this));
    m_expandedSignal->connect(*m_expandSlot);

    stringstream toggleExpandJs;
    toggleExpandJs << "function(s,e)"
      "{"
      ""  "if(" << jscontent << ".is(':hidden') )"
      ""   "{"
      "" << m_expandSlot->execJs("s", "e") << ";"
      ""   "}"
      ""  "else"
      ""  "{"
      ""  "" << m_collapseSlot->execJs("s", "e") << ";"
      ""  "}"
      "}";
    m_toggleExpandedStatusSlot.reset(new JSlot(toggleExpandJs.str(), this));
    title->doubleClicked().connect(*m_toggleExpandedStatusSlot);
    title->doubleClicked().connect(boost::bind(&AuxWindow::refresh, this));
  }//if( m_isPhone ) ' else

  m_titleText = new WText(windowTitle, XHTMLText, title);
  m_titleText->addStyleClass("titleVertMiddle");

  fcntcall = "function(){" + m_closedSignal->createEventCall(id(), "null") + "}";
  m_hideSlot.reset(new JSlot("function(s,e){ Wt.WT.AuxWindowHide(s,e,'" + id()
    + "'," + fcntcall + ");}", this));


  fcntcall = "function(){" + m_openedSignal->createEventCall(id(), "null") + "}";
  m_showSlot.reset(new JSlot("function(s,e){ Wt.WT.AuxWindowShow(s,e,'"
    + id() + "'," + fcntcall + ");}", this));

#if( USE_NEW_AUXWINDOW_ISH )
  m_repositionSlot.reset(new JSlot("function(s,e){Wt.WT.AuxWindowReposition(s,e);}", this));
  m_centerSlot.reset(new JSlot("function(s,e){Wt.WT.AuxWindowCenter(s,e);}", this));
  m_resizeSlot.reset(new JSlot("function(s,e){Wt.WT.AuxWindowResize(s,e);}", this));
  m_resizeScaledSlot.reset(new JSlot("function(s,e){Wt.WT.AuxWindowScaleResize(s,e);}", this));
#endif

  // Centering is done using either centerWindow() or repositionWindow()
  m_closedSignal->connect(boost::bind(&AuxWindow::setHidden, this, true, WAnimation()));
  m_openedSignal->connect(boost::bind(&AuxWindow::setHidden, this, false, WAnimation()));

  // We wont bring dialog to top, if the user clicked a button on the title (primarily
  //  for mobile, where we might have buttons in the title)
  const string bringToFront = "function(el){"
    "if( el && el.target && !$(el.target).hasClass('Wt-btn') && !$(el.target).hasClass('FooterHelpBtn'))"
      "Wt.WT.AuxWindowBringToFront('" + id() + "');"
  "}";
  title->clicked().connect(bringToFront);
  title->mouseWentDown().connect(bringToFront); //XXX - doesnt seem to work
  doJavaScript("$('#" + title->id() + "').mousedown(" + bringToFront + ");");

  if (!m_isPhone)
  {
    // TODO: actually respect/use the AuxWindowProperties::SetCloseable option!
    if (!properties.testFlag(AuxWindowProperties::SetCloseable))
      cout << "Warning: !AuxWindowProperties::SetCloseable not currently being respected for '"
      << windowTitle.toUTF8() << "' window." << endl;

    AuxWindow::setClosable(true);
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
    //    impl->setHiddenKeepsGeometry( true );
  }

  //footer()->setStyleClass("");
  //footer()->setHeight( 0 );

  setOffsets(WLength(0, WLength::Pixel), Wt::Left | Wt::Top);

  // By default it'll just spawn at 50% width and height of the host
//  resizeScaledWindow( 0.5, 0.5 );

  //so the signals of all the descendant widgets will be connected
  WDialog::setHidden(false, WAnimation());
  //  m_hideSlot->exec( "null", "{quietly:true, delay:5}" );
  if (m_expandIcon)
    doJavaScript("$('#" + m_expandIcon->id() + "').hide();");


  if (properties.testFlag(AuxWindowProperties::IsHelpWIndow))
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


  if (!m_isPhone && properties.testFlag(AuxWindowProperties::EnableResize))
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

    const string js = "$('#" + title->id() + "').mousedown( function(){"
      "$('#" + id() + "').data('centered',false);"
      "console.log('setting AuxWindow moved');"
      "} );";
    doJavaScript(js);

    title->touchMoved().connect("function(){$('#" + id() + "').data('centered',false);}"); //untested
  }
#endif //AUX_WINDOW_RE_CENTER_SIZE_ON_WINDOW_CHANGE
}//AuxWindow()


void AuxWindow::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  WDialog::render( flags );
  
  // Sometimes the contents of the AuxWindow will not get the contents laid out correctly, so we'll
  //  just always trigger a rather heavy handed resize event.
  // Definitely a hack
  if (flags & Wt::RenderFlag::RenderFull)
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

      /*
      // This next bit of JS accomplishes the same thing as the above, and being left commented
      //  out until
      const char* js = INLINE_JAVASCRIPT(
        setTimeout(function() {
        let c = $('.Wt-dialogcover');
        if (!c.length || c.data('ClickConnected')) return;
        c.data('ClickConnected', true);
        c.on('mousedown', function(evt){
          if ((evt.button != = 0) || (c[0].id != = evt.target.id) || (evt.pageY > 30)) return;
          window.wx.postMessage('MouseDownInTitleBar');
          evt.stopPropagation();
          evt.preventDefault();
        });
      }, 100);
      );
      doJavaScript(js);
      */
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
    bar->setHeight( WLength(35,WLength::Pixel) );
    bar->setStyleClass( "PhoneAuxWindowHeader" );
    
    //ToDo: for Android should show window title, but not in iOS
    //      (currently hiding everywhere)
    for( size_t i = 0; i < bar->children().size(); ++i )
    {
      WText *caption = static_cast<WText *>( bar->children()[i] );
      if( caption )
        caption->setHidden(true);
    }
    
    m_footer = new WContainerWidget;
    m_footer->addStyleClass( "PhoneAuxWindowFooter" );
    
    m_footer->setMargin( 12, Wt::Left );
    m_footer->setMargin( 12, Wt::Right );
    bar->insertWidget( 0, m_footer );
  }else
  {
    m_footer = WDialog::footer();
    m_footer->setHeight(WLength(45,WLength::Pixel));
    m_footer->setStyleClass("modal-footer");
  }
  
  return m_footer;
}

AuxWindow::~AuxWindow()
{
  WDialog::setModal( false ); //
//  if( m_titleText )
//    cerr << "AuxWindow destructing for: '"
//         << m_titleText->text().toUTF8() << "'" << endl;
  
  m_destructing = true;
}//~AuxWindow()


void AuxWindow::emitReject()
{
  finished().emit(WDialog::Rejected);
}

bool AuxWindow::isPhone() const
{
  return m_isPhone;
}

Wt::WPushButton *AuxWindow::addCloseButtonToFooter(string override_txt, bool floatright,  Wt::WContainerWidget * footerOverride)
{
  WPushButton *close = new WPushButton();
  
  if( m_isPhone )
  {
    if( override_txt.empty() )
      override_txt = "Back";
    
    close->addStyleClass( "MobileBackBtn" );
    if( m_isAndroid )
      close->addStyleClass( "InvertInDark" );
    
    close->setText( override_txt );
  }else
  {
    if( override_txt.empty() )
      override_txt = "Close";
    
    //If not fuule
    
    close->setText( override_txt );
    if( floatright )
      close->addStyleClass( "DialogClose" );
  }//if( phone ) / else
 
  //Sometimes, the footer may not be footer(), so we allow user to override
  if( !footerOverride )
    footerOverride = footer();
  footerOverride->insertWidget( footerOverride->count(), close );

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
    m_closeIcon = new WText( titleBar() );
    WApplication::instance()->theme()->apply(this, m_closeIcon,
                                              DialogCloseIconRole);
    //order of connections below matter
    m_closeIcon->clicked().connect( *m_expandSlot );
    m_closeIcon->clicked().connect( *m_hideSlot );
    
    //Trial Change 20181125: The m_hideSlot, when initiated from JS, will call
    //          back to c++ causing AuxWindow::setHidden to be called, which
    //          will then call emitReject(), so there *should* be no reason to
    //          explicitly call emitReject() when the closeIcon is clicked, and
    //          in fact the next line would cause the finished() signal to be
    //          emmitted twice.  (all this should be removed after more testing)
    //m_closeIcon->clicked().connect( this, &AuxWindow::emitReject );
  }else
  {
    if( !m_closeIcon )
      return;
    delete m_closeIcon;
    m_closeIcon = (WInteractWidget *)0;
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
    m_contentStretcher = new WGridLayout();
    WContainerWidget *container = WDialog::contents();
    if( container->count() )
    {
      cerr << endl << endl << endl << "AuxWindow::stretcher():"
           << " Warning, the window has children, but I'm adding a layout anyway"
           << endl;
      container->clear();
    }//if( container->count() )

    container->setLayout( m_contentStretcher );
    container->setOverflow( WContainerWidget::OverflowHidden );
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
  
#if( USE_NEW_AUXWINDOW_ISH )
  m_centerSlot->exec(id(),"null");
#else
  doJavaScript( "Wt.WT.AuxWindowCenter(" + this->jsRef() + ");" );
#endif
}

void AuxWindow::centerWindowHeavyHanded()
{
  if( m_isPhone )
    return;
  
  const string js = "var cntrfcn = function(){Wt.WT.AuxWindowCenter(" + this->jsRef() + ");};"
  "cntrfcn();"
  "setTimeout(cntrfcn,0);"
  "setTimeout(cntrfcn,100);"
  "setTimeout(cntrfcn,500);"
  "setTimeout(function(){Wt.WT.AuxWindowResizeToFitOnScreen('" + id() + "');},500);";
  
  doJavaScript( js );
}//void centerWindowHeavyHanded()


std::string AuxWindow::repositionWindowJs( int x, int y )
{
  if( m_isPhone )
    return ""; //disable running calls after AuxWindow initialized to be mobile
  
#if( USE_NEW_AUXWINDOW_ISH )
  const string pos = "{top:" + std::to_string(y)
  + ",left:" + std::to_string(x) + "}";
  return m_repositionSlot->execJs( id(), pos );
#else
  return string("var el = ") + this->jsRef() + ";"
  "var olddispl = el.style.display; el.style.display='';"
  + ((y==-32768) ? string() : "el.style.top = '"  + std::to_string(y) + "px';")
  + ((x==-32768) ? string() :  "el.style.left = '" + std::to_string(x) + "px';")
  + "el.style.display = olddispl;"
  "$(el).data('centered', false);";
#endif //#if( USE_NEW_AUXWINDOW_ISH )
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
    
#if( USE_NEW_AUXWINDOW_ISH )
  const string size = "{x:" + std::to_string(width)
  + ",y:" + std::to_string(height) + "}";
  m_resizeSlot->exec(id(),size);
#else
  const string jsthis = "$('#" + id() + "')";
  WStringStream sizeJs;
  sizeJs << "var el = " << this->jsRef() << ";"
         << "var olddispl = el.style.display; el.style.display='';"
         << "var ws = " << wApp->javaScriptClass() << ".WT.windowSize();";
  if( width > 0 )
    sizeJs << jsthis << ".width( Math.min("  << width << ",0.95*ws.x) );";
  if( height > 0 )
    sizeJs << jsthis << ".height( Math.min(" << height << ",0.95*ws.y) );";
  sizeJs << jsthis << ".data('notshown',true);"
         << "el.style.display = olddispl;";
  doJavaScript( sizeJs.str() );
#endif//#if( USE_NEW_AUXWINDOW_ISH )
} // void AuxWindow::resizeWindow( int width, int height )


std::string AuxWindow::resizeScaledWindowJs( double xRatio, double yRatio ) const
{
  if( m_isPhone )
    return ""; //disable running calls after AuxWindow initialized to be mobile
  
  if( xRatio <= 0.0 && yRatio <= 0.0 )
    return "";
  
#if( USE_NEW_AUXWINDOW_ISH )
  const string size = "{x:" + std::to_string(xRatio)
  + ",y:" + std::to_string(yRatio) + "}";
  return m_resizeScaledSlot->execJs(id(),size);
#else
  const string jsthis = "$('#" + id() + "')";

  stringstream sizeJs;
  sizeJs << "var el = " << this->jsRef() << ";"
         << "var olddispl = el.style.display; el.style.display='';";
//#if( IOS )
//  //At app startup, the first "Welcome..." screen wont be sized correctly if we
//  //  use the Wt.windowSize() method
//  sizeJs << "var ww = window.screen.width;"
//         << "var wh = window.screen.height;"
//         << "if(window.orientation==90 || window.orientation==-90) wh = [ww, ww = wh][0];";
  
//  if( xRatio > 0.0 )
//    sizeJs << "if(ww>2) " << jsthis << ".width( "  << xRatio << "*ww - 2 );";
//  if( yRatio > 0.0 )
//    sizeJs << "if(wh>2) " << jsthis << ".height( " << yRatio << "*wh - 2 );";
//#else
  sizeJs << "var ws = " << wApp->javaScriptClass() << ".WT.windowSize();";
  
  //TODO: Currently assuming a border with of 1px always, should evaluate this
  //TODO: would it be better to use $(window/document).width()/.height() below?
  if( xRatio > 0.0 )
    sizeJs << "if(ws.x>2) " << jsthis << ".width( "  << xRatio << "*ws.x - 2 );";
  if( yRatio > 0.0 )
    sizeJs << "if(ws.y>2) " << jsthis << ".height( " << yRatio << "*ws.y - 2 );";
//#endif
  sizeJs << jsthis << ".data('notshown',true);"
         << "el.style.display = olddispl;";
  
  return sizeJs.str();
#endif  //#if( USE_NEW_AUXWINDOW_ISH )
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
  
  if( window->m_modalOrig )
    window->WDialog::setModal(false);
  
  delete window;
}//void deleteAuxWindow( AuxWindow *window )


void AuxWindow::deleteSelf()
{
  AuxWindow::deleteAuxWindow( this );
}


void AuxWindow::rejectWhenEscapePressed( bool enable )
{
//  WDialog::rejectWhenEscapePressed( enable );
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
    /* see if this brings welcome screen above hamburger? */
    doJavaScript( "Wt.WT.AuxWindowBringToFront('" + id() + "');" );
    
    //Check original modal value, as modal doesn't work well with hide/show
    if (m_modalOrig)
        WDialog::setModal(true);
  setHidden( false, WAnimation() );
}//void show()


void AuxWindow::hide()
{
    //Check original modal value, as modal doesn't work well with hide/show
    if (m_modalOrig)
        WDialog::setModal(false);
    
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

//  WStringStream ss;
//	ss << WT_CLASS << ".animateDisplay('" << id()
//  << "'," << animation_.effects().value()
//  << "," << animation_.timingFunction()
//  << "," << animation_.duration()
//  << ",'" << element.getProperty(PropertyStyleDisplay)
//  << "');";
//	element.callJavaScript(ss.str());
  
  if( hide )
  {
    m_hideSlot->exec( "null", "{quietly: true}" );
    
    if( changingState && m_escapeIsReject )
    {
      m_escapeConnection1.disconnect();
      m_escapeConnection2.disconnect();
    }//if( m_escapeIsReject )
    
    emitReject();
  }else
  {
    m_showSlot->exec( "null", "{quietly: true}" );
    
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
    titleBar()->setPadding(WLength(5),Wt::Left);
    titleBar()->setPadding(WLength(2),Wt::Top);
  }
}



JSlot &AuxWindow::jsHide()
{
  return *m_hideSlot;
}


JSlot &AuxWindow::jsShow()
{
  return *m_showSlot;
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


//JSignal<> &AuxWindow::closed()
//{
//  return *m_closedSignal;
//}


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
  
  Wt::WImage *image = nullptr;
  
  image = new Wt::WImage(Wt::WLink("InterSpec_resources/images/help_minimal.svg"), footer);
  image->setStyleClass("Wt-icon FooterHelpBtn");
  image->setFloatSide( (viewer && viewer->isMobile()) ? Wt::Right : Wt::Left );
  
  image->setAlternateText("Help");
  image->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, page ) );
}
//
//void AuxWindow::openHelpWindow(std::string page,AuxWindow *parent)
//{
//  AuxWindow* m_helpWindow = new AuxWindow("Help", true);
//  m_helpWindow->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, m_helpWindow ) );
//  
//  WBorderLayout* layout = new Wt::WBorderLayout();
//  m_helpWindow->contents()->setLayout(layout);
//
//  WContainerWidget* m_helpWindowContent= new WContainerWidget();
//  m_helpWindowContent->setOverflow( WContainerWidget::OverflowAuto );
//  
//  layout->addWidget(m_helpWindowContent,Wt::WBorderLayout::Center);
//  layout->setContentsMargins(0, 0, 0, 0);
//  
//  m_helpWindow->footer()->setStyleClass("modal-footer");
//  WPushButton *closeButton = new WPushButton( "Close", m_helpWindow->footer() );
//  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, m_helpWindow ) );
//  closeButton->setFloatSide(Right);
//  
//  new WTemplate( WString::tr(page), m_helpWindowContent);
// 
//  m_helpWindow->footer()->setHeight(WLength(50,WLength::Pixel));
//  m_helpWindow->rejectWhenEscapePressed();
//  m_helpWindow->setMinimumSize( WLength(300, WLength::Pixel), WLength(300,WLength::Pixel));
//  m_helpWindow->resize( WLength(500, WLength::Pixel), parent->height());
//  m_helpWindow->setResizable( true );
//  m_helpWindow->show();
//  m_helpWindow->centerWindow();
//
//  //TODO: don't know why repositioning the parent to the top left
//  //then reposition the help window with positionAt() doesn't work
//  //
//  //  parent->repositionWindow(50,50);
//  //  m_helpWindow->positionAt(parent,Wt::Horizontal);
//}
