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

#include <Wt/WDialog>
#include <Wt/WServer>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include <string>

#include "InterSpec/InterSpec.h"  //to get rendered width/height
#include "InterSpec/SimpleDialog.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

SimpleDialog::SimpleDialog()
: Wt::WDialog( nullptr ), m_title( nullptr ), m_msgContents( nullptr )
{
  init( "", "" );
}

SimpleDialog::SimpleDialog( const Wt::WString &title )
 : Wt::WDialog( nullptr ), m_title( nullptr ), m_msgContents( nullptr )
{
  init( title, "" );
}


SimpleDialog::SimpleDialog( const Wt::WString &title, const Wt::WString &content )
 : Wt::WDialog( nullptr ), m_title( nullptr ), m_msgContents( nullptr )
{
  init( title, content );
}


void SimpleDialog::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  Wt::WDialog::render( flags );
  
  /*
  if( flags & RenderFull )
  {
    // WDialog::setMaximumSize will silently not use dimensions if WLength::Percentage
    //  Note that page dimensions wont be available during initial rendering of the webapp
    Wt::WDialog::render( flags );
    
    // Override some WDialog settings
    
    // The below seems to be necessary or else sometimes the window doesnt resize to fit its content
    doJavaScript( "let a = function(ms){"
                  + wApp->javaScriptClass() + ".layouts2.scheduleAdjust();"
                  " setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, ms );"
                  "};"
                  "a(0); a(50); a(250);"
    );
  }else
  {
    Wt::WDialog::render( flags );
  }//if( flags & RenderFull )
   */
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
    //setWindowTitle( title ); //Dont use the default Wt <h4>Some Text...</h4> stuff
    titleBar()->removeStyleClass( "titlebar" );  //Avoid the Wt changing of text color and such
    titleBar()->addStyleClass( "title" );
    m_title = new WText( title, titleBar() );
    m_title->setInline( false );
    //m_title->addStyleClass( "title" );
  }
  
  if( !content.empty() )
  {
    m_msgContents = new WText( content, contents() );
    m_msgContents->addStyleClass( "content" );
    m_msgContents->setInline( false );
  }
  
  setMinimumSize( WLength(260,WLength::Pixel), WLength::Auto );
  // We will set maximum size in a hacky way in SimpleDialog::render(...)
  
  show();
  finished().connect( this, &SimpleDialog::startDeleteSelf );
  
  // TODO: Wt 3.3.4 doesnt have WDialog::raiseToFront(), should simulate, or fix AuxWindow to
  //       better respect the DialogCover stuff defined in 
}//init(...)


SimpleDialog::~SimpleDialog()
{
  //cerr << "Deleting simpledialog" << endl;
}


Wt::WPushButton *SimpleDialog::addButton( const Wt::WString &txt )
{
  Wt::WPushButton *b = new WPushButton( txt, footer() );
  b->setStyleClass( "simple-dialog-btn" );
  
  // TODO: closing the dialog seems a little laggy; check if WDialog::hide is faster, or if we
  //       should stick to using the JS
  //b->clicked().connect( this, &WDialog::hide );
  b->clicked().connect( "function(){$('#" + id() + "').hide(); $('.Wt-dialogcover').hide();}" );
  
  b->clicked().connect( boost::bind( &WDialog::done, this, Wt::WDialog::DialogCode::Accepted ) );
  return b;
}


void SimpleDialog::startDeleteSelf()
{
  if( isModal() )
    setModal(false);
  
  // We'll actually delete the windows later on in the event loop incase the order of connections
  //  to its signals is out of intended order, but also we will protect against being deleted in the
  //  current event loop as well
  boost::function<void(void)> updater = wApp->bind( boost::bind( &SimpleDialog::deleteSelf, this ) );
  WServer::instance()->post( wApp->sessionId(), std::bind( [updater](){
    updater();
    wApp->triggerUpdate();
  } ) );
}//startDeleteSelf()


void SimpleDialog::deleteSelf()
{
  // calling 'delete this' makes me feel uneasy
  delete this;
}

