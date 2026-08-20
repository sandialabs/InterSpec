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

#include <Wt/WLength.h>
#include <Wt/WDialog.h>
#include <Wt/WServer.h>
#include <Wt/WTemplate.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include <string>
#include <cassert>
#include <algorithm>

#include "InterSpec/InterSpec.h"  //for InterSpec::instance()
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/SimpleDialog.h"

#if( BUILD_AS_WX_WIDGETS_APP || BUILD_AS_ELECTRON_APP )
#include "InterSpec/InterSpecApp.h"  //for InterSpecApp::isPrimaryWindowInstance()
#endif

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

WT_DECLARE_WT_MEMBER
 (SimpleDialogBringToFront, Wt::JavaScriptFunction, "SimpleDialogBringToFront",
  function( id )
  {
   const maxz = Array.from(document.querySelectorAll('.Wt-dialog,.MobileMenuButton')).reduce( function(result, item){
     if( item.id === id ) return result;
     const z = parseInt( getComputedStyle(item).zIndex );
     return isNaN(z) ? result : Math.max(result, z);
    }, 0);

   const el = document.getElementById(id);
   const z = el ? parseInt( getComputedStyle(el).zIndex ) : NaN;

   if( el && (isNaN(z) || maxz >= z) )
     el.style.zIndex = maxz + 1;
   var wcc = document.querySelector('.window-controls-container');
   if( wcc ) wcc.style.zIndex = maxz + 2; //for wxWidgets and Electron builds
   document.querySelectorAll('.suggestion').forEach( function(s){ s.style.zIndex = maxz + 2; } );
 }
);


namespace
{
  /** Default height of everything in the dialog that is not the scrollable body - the title bar plus
   the footer.  See SimpleDialog::setBodyChromeHeight().
   */
  const int sm_defaultBodyChromePx = 90;

  /** The most of the browser window a dialog may take; mirrors the `max-height` on `.simple-dialog`
   in InterSpec_resources/SimpleDialog.css, and must stay in sync with it.
   */
  const double sm_maxFractionOfWindow = 0.95;

  /** Resolves a WLength to pixels for the current browser window size, or returns -1 if it cannot
   be (an `auto` length, or a percentage of a parent we have no size for here).
   */
  double length_in_px( const Wt::WLength &length, const int windowWidth, const int windowHeight )
  {
    if( length.isAuto() )
      return -1.0;

    switch( length.unit() )
    {
      case WLength::Unit::ViewportWidth:
        return (windowWidth > 100) ? (0.01 * length.value() * windowWidth) : -1.0;

      case WLength::Unit::ViewportHeight:
        return (windowHeight > 100) ? (0.01 * length.value() * windowHeight) : -1.0;

      case WLength::Unit::ViewportMin:
        return (windowWidth > 100 && windowHeight > 100)
                 ? (0.01 * length.value() * std::min(windowWidth,windowHeight)) : -1.0;

      case WLength::Unit::ViewportMax:
        return (windowWidth > 100 && windowHeight > 100)
                 ? (0.01 * length.value() * std::max(windowWidth,windowHeight)) : -1.0;

      case WLength::Unit::Percentage:
        return -1.0;

      default:
        return length.toPixels();  //px, in, cm, mm, pt, pc, em, ex
    }//switch( length.unit() )
  }//double length_in_px(...)
}//namespace


SimpleDialog::SimpleDialog()
: Wt::WDialog(),
  m_title( nullptr ),
  m_msgContents( nullptr ),
  m_multipleBringToFront( true ),
  m_bodyChromeHeight( sm_defaultBodyChromePx ),
  m_bodyPreferredHeight( -1.0 )
{
  init( "", "" );
}


SimpleDialog::SimpleDialog( const Wt::WString &title )
 : Wt::WDialog(),
  m_title( nullptr ),
  m_msgContents( nullptr ),
  m_multipleBringToFront( true ),
  m_bodyChromeHeight( sm_defaultBodyChromePx ),
  m_bodyPreferredHeight( -1.0 )
{
  init( title, "" );
}


SimpleDialog::SimpleDialog( const Wt::WString &title, const Wt::WString &content )
 : Wt::WDialog(),
  m_title( nullptr ),
  m_msgContents( nullptr ),
  m_multipleBringToFront( true ),
  m_bodyChromeHeight( sm_defaultBodyChromePx ),
  m_bodyPreferredHeight( -1.0 )
{
  init( title, content );
}


void SimpleDialog::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  Wt::WDialog::render( flags );
  
  if( flags.test( Wt::RenderFlag::Full ) )
  {
    // WDialog::setMaximumSize will silently not use dimensions if WLength::Unit::Percentage, so we use CSS.
    //  Note that page dimensions wont be available during initial rendering of the webapp
    
    // The below seems to be necessary or else sometimes the window doesnt resize to fit its content
    wApp->doJavaScript( wApp->javaScriptClass() + ".TriggerResizeEvent();" );

#if( BUILD_AS_WX_WIDGETS_APP || BUILD_AS_ELECTRON_APP )
    // To allow moving window around when dialog showing; see note in AuxWindow::render 
    //  for the same code snippet
    if( InterSpecApp::isPrimaryWindowInstance() )
    {
#if( BUILD_AS_WX_WIDGETS_APP )
      WWidget* coverw = wApp->findWidget("dialog-cover");
      WContainerWidget* dialog_cover = dynamic_cast<WContainerWidget*>(coverw);
      if (dialog_cover && !dialog_cover->mouseWentDown().isConnected())
        dialog_cover->mouseWentDown().connect(wApp->javaScriptClass() + ".MouseDownOnDialogCover");
#endif

      // Raise windows controls (minimize, maximize, close), to above the dialog-cover.
      wApp->doJavaScript(wApp->javaScriptClass() + ".RaiseWinCntrlsAboveCover();");
    }
#endif
    
    // On mobile, it seems Wt.WT.AuxWindowBringToFront(...) may get called after this window is
    //  created (happens on the "QR code" link on Nuclide Decay Tool - since the user clicks
    //  a button in the titlebar), which will bring that dialog above this one - which isnt wanted,
    //  so we'll manually bring this dialog to the top on a delay.
    //  We'll add this JS, even on non-mobile, JIC
    LOAD_JAVASCRIPT(wApp, "SimpleDialog.cpp", "SimpleDialog", wtjsSimpleDialogBringToFront);
    
    const string time_delays_array = m_multipleBringToFront ? "[5,100,500]" : "[5]";
    
    doJavaScript( "for( const d of " + time_delays_array + "){"
                    "setTimeout( function(){ Wt.WT.SimpleDialogBringToFront('" + id() + "');}, d);"
                  "}");
  }//if( flags & RenderFull )
}//render( flags )


void SimpleDialog::init( const Wt::WString &title, const Wt::WString &content )
{
  wApp->useStyleSheet( "InterSpec_resources/SimpleDialog.css" );
  
  addStyleClass( "simple-dialog" );
  
  setModal( true );
  
  setMovable( false );
  
  if( title.empty() )
  {
    setTitleBarEnabled( false );
  }else
  {
    setTitleBarEnabled( true );
    // Wt4's WDialog creates a caption WTemplate with an unbound ${title} variable
    //  that renders as "??title??". Hide the caption element so it doesn't show;
    //  we use our own m_title WText instead.
    for( WWidget *child : titleBar()->children() )
    {
      if( dynamic_cast<WTemplate *>( child ) )
      {
        child->setHidden( true );
        break;
      }
    }
    titleBar()->removeStyleClass( "titlebar" );  //Avoid the Wt changing of text color and such
    titleBar()->addStyleClass( "title" );
    m_title = titleBar()->addNew<WText>( title );
    m_title->setInline( false );
    //m_title->addStyleClass( "title" );
  }
  
  if( !content.empty() )
  {
    m_msgContents = contents()->addNew<WText>( content );
    m_msgContents->addStyleClass( "content" );
    m_msgContents->setInline( false );
  }
  
  // We need to set the minimum size in C++; the dialogs maximum size is set in CSS, but the bodys
  //  has to come from here - see updateBodySizeForWindow().
  setMinimumSize( WLength(260,WLength::Unit::Pixel), WLength::Auto );
  updateBodySizeForWindow();

  show();
  finished().connect( this, &SimpleDialog::startDeleteSelf );
  
#if( WT_VERSION > 0x3040000 )
  // I havent checked version of Wt that does include `raiseToFront()`, but 3.3.4 doesnt.
  //20250228: raiseToFront() call removed do to encountering some JS exception, that _maybe_ have something to do with this (the object not found), on Windows - I dont actually know if this function call is the problem.
  //raiseToFront();
#endif
}//init(...)


SimpleDialog::~SimpleDialog()
{
}


void SimpleDialog::setMaximumSize( const Wt::WLength &width, const Wt::WLength &height )
{
  // Native sizing of the outer dialog and its inner layout.  Wt4's WLength supports viewport units
  //  (vw/vh) that Wt3 lacked, so callers can pass e.g. WLength(95, LengthUnit::ViewportWidth).
  //  (WDialog::setMaximumSize silently drops Percentage units for the inner layout, but not vw/vh.)
  WDialog::setMaximumSize( width, height );

  // The call above does not reach the scrollable `.body` element, whose height has to be limited
  //  separately.  Its width needs no help: the dialog layout gives it `max-width: 100%`, so it
  //  already follows whatever the dialog itself is limited to.
  updateBodySizeForWindow();
}//setMaximumSize(...)


void SimpleDialog::setMaxWidth( const Wt::WLength &width )
{
  // Width-only convenience; preserves any existing maximum height.
  setMaximumSize( width, maximumHeight() );
}//setMaxWidth(...)


void SimpleDialog::setBodyChromeHeight( const int pixels )
{
  m_bodyChromeHeight = std::max( 0, pixels );
  updateBodySizeForWindow();
}//setBodyChromeHeight(...)


void SimpleDialog::setBodyPreferredHeight( const double pixels )
{
  m_bodyPreferredHeight = pixels;
  updateBodySizeForWindow();
}//setBodyPreferredHeight(...)


void SimpleDialog::updateBodySizeForWindow()
{
  const InterSpec * const viewer = InterSpec::instance();
  const int windowWidth = viewer ? viewer->renderedWidth() : 0;
  const int windowHeight = viewer ? viewer->renderedHeight() : 0;

  // How much of the window is left for the body, once the dialogs own 95vh cap and its title bar
  //  and footer are taken off.  Negative means we dont know the window size yet (the very first
  //  render), in which case we leave the body to size to its content.
  double maxHeight = (windowHeight > 100)
                       ? (sm_maxFractionOfWindow*windowHeight - m_bodyChromeHeight) : -1.0;

  // A caller-supplied maximum for the whole dialog wins, if it is the more restrictive of the two.
  const double dialogMaxHeight = length_in_px( maximumHeight(), windowWidth, windowHeight );
  if( dialogMaxHeight > 0.0 )
    maxHeight = (maxHeight > 0.0) ? std::min( maxHeight, dialogMaxHeight - m_bodyChromeHeight )
                                  : (dialogMaxHeight - m_bodyChromeHeight);

  // Dont let a tiny window (or a large chrome allowance) collapse the body to nothing.
  if( maxHeight > 0.0 )
    maxHeight = std::max( maxHeight, 100.0 );

  WContainerWidget * const body = contents();
  assert( body );
  if( !body )
    return;

  body->setMaximumSize( body->maximumWidth(),
                        (maxHeight > 0.0) ? WLength(maxHeight,WLength::Unit::Pixel) : WLength::Auto );

  if( m_bodyPreferredHeight > 0.0 )
  {
    const double height = (maxHeight > 0.0) ? std::min( m_bodyPreferredHeight, maxHeight )
                                            : m_bodyPreferredHeight;
    body->setHeight( WLength(height,WLength::Unit::Pixel) );
  }
}//updateBodySizeForWindow()


void SimpleDialog::doNotUseMultpleBringstoFront()
{
  m_multipleBringToFront = false;
}


void SimpleDialog::rejectWhenEscapePressed( bool enable )
{
  WDialog::rejectWhenEscapePressed( enable );
  
  if( enable == m_escapeConnection1.isConnected() )
    return;
  
  if( enable )
  {
    WWidget * const implw = implementation();
    WTemplate * const impl = dynamic_cast<WTemplate *>( implw );
    if( impl )
      m_escapeConnection1 = impl->escapePressed().connect( this, &SimpleDialog::reject );
  }else
  {
    m_escapeConnection1.disconnect();
  }
}//void rejectWhenEscapePressed( bool enable )

Wt::WPushButton *SimpleDialog::addButton( const Wt::WString &txt )
{
  Wt::WPushButton *b = footer()->addNew<WPushButton>( txt );
  b->setStyleClass( "simple-dialog-btn" );

  // TODO: closing the dialog seems a little laggy; check if WDialog::hide is faster, or if we
  //       should stick to using the JS
  //b->clicked().connect( this, &WDialog::hide );
  b->clicked().connect( "function(){document.getElementById('" + id() + "').style.display='none'; document.querySelectorAll('.Wt-dialogcover').forEach(function(e){e.style.display='none';});}" );

  b->clicked().connect( this, [this](){ done( Wt::DialogCode::Accepted ); } );
  return b;
}//addButton(...)


void SimpleDialog::startDeleteSelf()
{
  if( isModal() )
    setModal(false);
  
  // We'll actually delete the windows later on in the event loop incase the order of connections
  //  to its signals is out of intended order, but also we will protect against being deleted in the
  //  current event loop as well.
  //  Only the dialog id crosses into the task: `post` copy-constructs the std::function and only
  //  guarantees the *session* is still alive - another path (an owner's destructor,
  //  deleteSimpleDialog) may have destroyed us in the meantime.  An observing_ptr must not be used
  //  here; it would mutate `observable::observers_`, which has no mutex, from the io-service thread.
  const string sessionId = wApp->sessionId();
  const WidgetUtils::WidgetHandle self( this );
  WServer::instance()->post( sessionId, [self, sessionId](){
    auto *app = WApplication::instance();
    if( !app || (app->sessionId() != sessionId) )
      return;

    SimpleDialog * const dialog = self.resolve_as<SimpleDialog>();
    if( !dialog )
      return;  //already destroyed

    dialog->deleteSelf();
    app->triggerUpdate();
  } );
}//startDeleteSelf()


void SimpleDialog::deleteSelf()
{
  // The SimpleDialog::make() factory gives ownership to wApp via addChild(),
  //  so we use removeChild() to properly release and destruct this dialog.
  Wt::WApplication::instance()->removeChild( this );
}


void SimpleDialog::deleteSimpleDialog( SimpleDialog *dialog )
{
  if( !dialog )
    return;

  // Clear modal first so the dialog cover is popped (same reason AuxWindow::deleteAuxWindow does),
  //  then hand the widget back from wApp and destroy it.  Deliberately *not* done()/reject(): those
  //  emit finished(), whose handlers usually call back into the owner that is tearing us down.
  if( dialog->isModal() )
    dialog->setModal( false );

  WApplication * const app = WApplication::instance();
  if( app )
    app->removeChild( dialog );
}//void deleteSimpleDialog( SimpleDialog *dialog )

