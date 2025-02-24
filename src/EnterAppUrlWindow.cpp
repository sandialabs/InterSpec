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

#include <cctype>
#include <string>
#include <iostream>

#include <Wt/WText>
#include <Wt/WTextArea>
#include <Wt/WPushButton>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UndoRedoManager.h"

using namespace std;
using namespace Wt;

namespace
{
  
  template<typename charT>
  struct char_iequal
  {
    char_iequal( const std::locale &loc ) : m_loc(loc) {}
    bool operator()(charT ch1, charT ch2) {
      return std::toupper(ch1, m_loc) == std::toupper(ch2, m_loc);
    }
  private:
    const std::locale &m_loc;
  };
  
}//namespace

namespace EnterAppUrlWindow
{
SimpleDialog *createEntryWindow( InterSpec *viewer )
{
  if( !viewer )
    return nullptr;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  
  SimpleDialog *window = new SimpleDialog( "Enter URL", "" );
  
  WTextArea *text = new WTextArea( window->contents() );
  text->setObjectName( "txtarea" );
  text->setInline( false );
  
  int w = viewer->renderedWidth();
  //int h = viewer->renderedHeight();
  //if( viewer->isMobile() && (w < 100) )
  //{
  //  w = wApp->environment().screenWidth();
  //  h = wApp->environment().screenHeight();
  //}
  
  text->setWidth( ((w < 100) || (w > 530)) ? 450 : w - 80 );
  
  const char *desctxt = "<div style=\"margin-top: 10px;\">Enter, usually through copy/paste, InterSpec URLs.</div>"
                        "<div>These URLs start with either <code>interspec://</code>, or <code>raddata://</code>.</div>";
  WText *desc = new WText( desctxt, window->contents() );
  desc->addStyleClass( "content" );
  desc->setInline( false );
  
  WPushButton *cancel = window->addButton( "Cancel" );
  WPushButton *okay = window->addButton( "Okay" );
  
  // To be generous to the user, we will allow them to have stuff before the relevant part of the
  //  URI, so we'll serach for URIs beginning with following paths:
  const string acceptable_paths[] = {
    "G0", "drf", "specsum", "flux", "specexport", "decay", "dose", "gammaxs", "1overr2", "unit", "about", "welcome"
#if( USE_REMOTE_RID )
    , "remoterid"
#endif
#if( USE_DETECTION_LIMIT_TOOL )
    , "detection-limit"
    , "simple-mda"
#endif
  };//acceptable_paths[]
    
  
  // We will validate the URL starts with 'interspec://' or 'raddata://', end enable/disable the
  //  "Okay" button in javascript
  const string jsokay = okay->jsRef();
  string enable_fcn = "function(value){"
  "const matches = /^((interspec:\\/\\/)|(raddata:\\/\\/))";
  for( const string &path : acceptable_paths )
    enable_fcn += "|(.*\\/" + path + "\\/)";
  enable_fcn += ".+/i.test(value);"
    "$(" + jsokay + ").prop('disabled', !matches);"
  "}";
  
  string validate_js = "function(textarea,event){ "
   "(" + enable_fcn + ")(textarea.value);"
  " }";
  text->keyWentUp().connect( validate_js );
  
  text->doJavaScript( "$(" + text->jsRef() + ").on('paste', function(){"
      "setTimeout(function(){"
          "(" + enable_fcn + ")(" + text->jsRef() + ".value);"
      "},1);"
  "} );"
  );
  
  // Dont disable okay button during undo/redo - would be better to validate after we might put text
  //  into the textarea - but not bothering to do it correctly, at the moment.
  if( !undoRedo || !undoRedo->isInUndoOrRedo() )
    okay->doJavaScript( "$(" + jsokay + ").prop('disabled', true);" );
  
  // TODO: All this undo/redo stuff is still probably not quite right - could use a little more work
  
  okay->clicked().connect( std::bind([viewer, text, acceptable_paths](){
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    
    string uri = text->text().toUTF8();
    
    // Remove control characters
    SpecUtils::erase_any_character( uri, "\n\r\t\b\f\a" );
    //
    // note: we are not currently using `std::iscntrl(...)`, like below, because it is locale
    //       dependent, and not tested
    //uri.erase( std::remove_if(begin(uri), end(uri),
    //             [](char ch){return std::iscntrl(static_cast<unsigned char>(ch));}),
    //            end(uri) );
    
    // Remove leading/trailing white spaces.
    //  Will leave spaces within uri for the moment though, as these could be valid if uri has
    //  already been uri-decoded
    SpecUtils::trim( uri );
    
    viewer->handleAppUrlClosed();
    
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      // Note: if you loaded a spectrum file, through a "RadData G0" URI, then undo will be
      //       disconnected anyway, so it wont work.
      //       Otherwise, the undo should close the decay, xs, or whatever window, then you
      //       will have to do another undo to bring up this dialog.
      auto undo = [viewer,uri](){
        SimpleDialog *window = viewer->makeEnterAppUrlWindow();
        if( !window )
          return;
        WTextArea *t = dynamic_cast<WTextArea *>( window->contents()->find( "txtarea" ) );
        if( t )
          t->setText( WString::fromUTF8(uri) );
      };
      
      auto redo = [viewer,uri](){
        SimpleDialog *window = viewer->makeEnterAppUrlWindow();
        if( window )
          window->accept();
        
        try
        {
          viewer->handleAppUrl( uri );
        }catch( std::exception &e )
        {
          passMessage( "Error handling URL: " + string(e.what()), WarningWidget::WarningMsgHigh );
        }
      };//redo lamda
      
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Enter URL Window" );
    }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
    
    try
    {
      if( SpecUtils::istarts_with(uri, "interspec://") )
      {
        viewer->handleAppUrl( uri );
      }else if( SpecUtils::istarts_with(uri, "raddata://") )
      {
        viewer->handleAppUrl( uri );
      }else
      {
        bool found_sub = false;
        for( const string &path : acceptable_paths )
        {
          const string key = "/" + path + "/";
          const auto it = std::search( begin(uri), end(uri), begin(key), end(key),
                                       char_iequal<char>(std::locale()) );
          if( it != end(uri) )
          {
            const string sub_uri = "interspec:/" + string(it, end(uri));
            viewer->handleAppUrl( sub_uri );
            found_sub = true;
            break;
          }//if( it != end(uri) )
        }//for( const string &path : acceptable_paths )
       
        if( !found_sub )
          throw runtime_error( "URL must start with 'interspec://', 'raddata://', or have a /G0/ component." );
      }//
    }catch( std::exception &e )
    {
      SimpleDialog *errdialog = new SimpleDialog( "Error with entered URL", e.what() );
      errdialog->addButton( "Okay" );
    }
  }) );
  
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    cancel->clicked().connect( std::bind([viewer,text](){
      viewer->handleAppUrlClosed();
      
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      
      if( !undoRedo || !undoRedo->canAddUndoRedoNow() )
        return;
      
      const string txtval = text->text().toUTF8();
      auto undo = [viewer,txtval](){
        SimpleDialog *window = viewer->makeEnterAppUrlWindow();
        if( !window )
          return;
        
        WTextArea *t = dynamic_cast<WTextArea *>( window->contents()->find( "txtarea" ) );
        if( t )
          t->setText( WString::fromUTF8(txtval) );
      };
      auto redo = [viewer](){
        SimpleDialog *window = viewer->makeEnterAppUrlWindow();
        viewer->handleAppUrlClosed();
        if( window )
          window->accept();
      };
      
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Enter URL" );
    }) );
    
    auto undo = [viewer](){
      SimpleDialog *window = viewer->makeEnterAppUrlWindow();
      viewer->handleAppUrlClosed();
      if( window )
        window->accept();
    };
    
    auto redo = [viewer](){ viewer->makeEnterAppUrlWindow(); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Show Enter URL Window" );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  
  return window;
}//SimpleDialog *createEntryWindow( InterSpec *viewer )
  
}//namespace EnterAppUrlWindow

