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

#include <iostream>
#include <algorithm>

#include <Wt/WLink.h>
#include <Wt/WText.h>
#include <Wt/Utils.h>
#include <Wt/WLabel.h>
#include <Wt/WImage.h>
#include <Wt/WCheckBox.h>
#include <Wt/WGridLayout.h>
#include <Wt/WJavaScript.h>
#include <Wt/WApplication.h>
#include <Wt/WStandardItem.h>
#include <Wt/WStringStream.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WStandardItemModel.h>

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RowStretchTreeView.h"

using namespace Wt;
using namespace std;

const char *iconUrl( const WarningWidget::WarningMsgLevel level )
{
  switch( level )
  {
    case WarningWidget::WarningMsgLevel::WarningMsgInfo:     return "InterSpec_resources/images/information.svg";
    case WarningWidget::WarningMsgLevel::WarningMsgLow:      return "InterSpec_resources/images/caution.svg";
    case WarningWidget::WarningMsgLevel::WarningMsgMedium:   return "InterSpec_resources/images/error.svg";
    case WarningWidget::WarningMsgLevel::WarningMsgHigh:     return "InterSpec_resources/images/exclamation.svg";
    case WarningWidget::WarningMsgLevel::WarningMsgSave:     return "InterSpec_resources/images/disk.svg";
    case WarningWidget::WarningMsgLevel::WarningMsgShowOnBoardRiid:
#if( USE_REMOTE_RID )
    case WarningWidget::WarningMsgLevel::WarningMsgExternalRiid:
#endif
      return "InterSpec_resources/images/search_results.svg";
      
    case WarningWidget::WarningMsgLevel::NumWarningMsgType:  return "";
  }//switch( WarningMsgLevel(level) )
  
  return "";
}//iconUrl(...)


const char *WarningWidget::tostr( const WarningMsgLevel level )
{
  switch( level )
  {
    case WarningWidget::WarningMsgInfo:
      return "BlockNotificationInfo";
    case WarningWidget::WarningMsgLow:
      return "BlockNotificationLow";
    case WarningWidget::WarningMsgMedium:
      return "BlockNotificationMedium";
    case WarningWidget::WarningMsgHigh:
      return "BlockNotificationHigh";
    case WarningWidget::WarningMsgSave:
      return "BlockNotificationSave";
    case WarningMsgLevel::WarningMsgShowOnBoardRiid:
      return "BlockNotificationRiid";
#if( USE_REMOTE_RID )
    case WarningWidget::WarningMsgLevel::WarningMsgExternalRiid:
      assert( 0 );
      return "BlockNotificationExternalRiid";
#endif
    case WarningMsgLevel::NumWarningMsgType:
      return "";
  }//switch( i )
  
  throw runtime_error( "WarningWidget::tostr(...): invalid level" );
  return "";
}//tostr(...)


const char *WarningWidget::popupToStr( const WarningMsgLevel level )
{
  switch( level )
  {
    case WarningWidget::WarningMsgInfo:
      return "PopupBlockNotificationInfo";
    case WarningWidget::WarningMsgLow:
      return "PopupBlockNotificationLow";
    case WarningWidget::WarningMsgMedium:
      return "PopupBlockNotificationMedium";
    case WarningWidget::WarningMsgHigh:
      return "PopupBlockNotificationHigh";
    case WarningWidget::WarningMsgSave:
      return "PopupBlockNotificationSave";
    case WarningWidget::WarningMsgShowOnBoardRiid:
      return "PopupBlockNotificationRiid";
#if( USE_REMOTE_RID )
    case WarningWidget::WarningMsgLevel::WarningMsgExternalRiid:
      assert( 0 );
      return "PopupBlockNotificationExternalRiid";
#endif
    case WarningWidget::NumWarningMsgType:
      return "PopupBlockNotificationNumTypes";
  }//switch( i )
  
  throw runtime_error( "WarningWidget::tostr(...): invalid level" );
  return "";
}//tostr(...)


const char *WarningWidget::description( const WarningMsgLevel level )
{
  switch( level )
  {
    case WarningWidget::WarningMsgInfo:
      return "info";
    case WarningWidget::WarningMsgLow:
      return "low";
    case WarningWidget::WarningMsgMedium:
      return "medium";
    case WarningWidget::WarningMsgHigh:
      return "high";
    case WarningWidget::WarningMsgSave:
      return "save";
    case WarningWidget::WarningMsgShowOnBoardRiid:
#if( USE_REMOTE_RID )
      return "On-board RIID";
#else
      return "RIID";
#endif
      
#if( USE_REMOTE_RID )
    case WarningWidget::WarningMsgLevel::WarningMsgExternalRiid:
      return "External RID";
#endif
      
    case WarningWidget::NumWarningMsgType:
      return "";
  }//switch( i )
  
  throw runtime_error( "WarningWidget::description(...): invalid level" );
  return "";
}//const char *description( const WarningMsgLevel level )


WarningWidget::WarningWidget( InterSpec *hostViewer )
: WContainerWidget(),
  m_hostViewer( hostViewer ),
  m_totalMessages(0),
  m_popupActive{ false },
  m_active{ false },
  m_layout(NULL),
  m_messageModel( nullptr ),
  m_tableView(NULL),
  m_description(NULL)
{
  setOffsets( WLength(0, WLength::Unit::Pixel), Wt::Side::Left | Wt::Side::Top );
  m_messageModel = std::make_shared<WStandardItemModel>( 1, 3 );
  
  // Find which messages should be active.
  for( WarningMsgLevel i = WarningMsgLevel(0); i <= WarningMsgHigh; i = WarningMsgLevel(i+1) )
    m_active[i] = !UserPreferences::preferenceValue<bool>(tostr(i), m_hostViewer);
  
  for( WarningMsgLevel i = WarningMsgLevel(0); i <= WarningMsgHigh; i = WarningMsgLevel(i+1) )
    m_popupActive[i] = !UserPreferences::preferenceValue<bool>(popupToStr(i), m_hostViewer);
  
  //Force WarningMsgSave to always be true
  m_active[WarningMsgSave] = true;
  m_popupActive[WarningMsgSave] = true;
  
  //Force ShowRiid to always be true
  m_active[WarningMsgShowOnBoardRiid] = true;
  m_popupActive[WarningMsgShowOnBoardRiid] = true;
  
#if( USE_REMOTE_RID )
  m_active[WarningMsgExternalRiid] = true;
  m_popupActive[WarningMsgExternalRiid] = true;
#endif
  
  // Hook it up to the message handler in InterSpec
  hostViewer->messageLogged().connect( this, &WarningWidget::addMessage );
}//WarningWidget::WarningWidget()


bool WarningWidget::active( WarningWidget::WarningMsgLevel level ) const
{
  return m_active[level];
} //bool WarningWidget::active( WarningWidget::WarningMsgLevel level ) const


WarningWidget::~WarningWidget()
{
  // no-op
} // WarningWidget::~WarningWidget()


void WarningWidget::createContent()
{
  if (!m_layout) //Create UI for first time
  {
    m_layout = setLayout( std::make_unique<WGridLayout>() );
    m_layout->setContentsMargins(0,0,0,0);
    
    setPadding(0);
    setMargin(0);
    m_messageModel->setItem( 0, 0, std::make_unique<WStandardItem>("0") );
    m_messageModel->setItem( 0, 1, std::make_unique<WStandardItem>("") );
    m_messageModel->setItem( 0, 2, std::make_unique<WStandardItem>("InterSpec started.") );
    m_messageModel->setItem( 0, 3, std::make_unique<WStandardItem>("") );
    m_messageModel->setHeaderData(  0, Wt::Orientation::Horizontal, WString("#"), Wt::ItemDataRole::Display );
    m_messageModel->setHeaderData(  1, Wt::Orientation::Horizontal, WString("Type"), Wt::ItemDataRole::Display );
    m_messageModel->setHeaderData(  2, Wt::Orientation::Horizontal, WString("Message"), Wt::ItemDataRole::Display );
    
    m_tableView = m_layout->addWidget( std::make_unique<RowStretchTreeView>(), 0, 0, 1, 6 );
    m_tableView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
    m_tableView->setModel(m_messageModel);
    m_tableView->setOffsets( WLength(0, WLength::Unit::Pixel), Wt::Side::Left | Wt::Side::Top );
    m_tableView->setColumnResizeEnabled(true);
    m_tableView->setColumnAlignment(0, Wt::AlignmentFlag::Center);
    m_tableView->setHeaderAlignment(0, Wt::AlignmentFlag::Center);
    m_tableView->setColumnWidth(0,50);
    m_tableView->setColumnWidth(1,100);
    m_tableView->setColumnWidth(2,400);

    m_tableView->setAlternatingRowColors(true);
    m_tableView->setSelectionMode(Wt::SelectionMode::None);
    m_tableView->setEditTriggers(Wt::WFlags<Wt::EditTrigger>{Wt::EditTrigger::None});

    m_tableView->setSelectionMode( Wt::SelectionMode::Single );
    m_tableView->setSelectionBehavior( Wt::SelectionBehavior::Rows );
    m_tableView->selectionChanged().connect( this, &WarningWidget::resultSelectionChanged );

    m_layout->setRowStretch(0,1);
    m_layout->setColumnStretch(0,1);

    {
      Wt::WImage *image = m_layout->addWidget( std::make_unique<Wt::WImage>(Wt::WLink( iconUrl(WarningMsgLevel::WarningMsgInfo) )), 1, 2, Wt::AlignmentFlag::Center );
      image->setMaximumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
      image->setMinimumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
    }

    {
      Wt::WImage *image = m_layout->addWidget( std::make_unique<Wt::WImage>(Wt::WLink( iconUrl(WarningMsgLevel::WarningMsgLow) )), 1, 3, Wt::AlignmentFlag::Center );
      image->setMaximumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
      image->setMinimumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
    }

    {
      Wt::WImage *image = m_layout->addWidget( std::make_unique<Wt::WImage>(Wt::WLink( iconUrl(WarningMsgLevel::WarningMsgMedium) )), 1, 4, Wt::AlignmentFlag::Center );
      image->setMaximumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
      image->setMinimumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
    }

    {
      Wt::WImage *image = m_layout->addWidget( std::make_unique<Wt::WImage>(Wt::WLink( iconUrl(WarningMsgLevel::WarningMsgHigh) )), 1, 5, Wt::AlignmentFlag::Center );
      image->setMaximumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
      image->setMinimumSize(WLength(16,WLength::Unit::Pixel), WLength(16,WLength::Unit::Pixel));
    }

    m_layout->addWidget( std::make_unique<Wt::WLabel>("info"), 2, 2, Wt::AlignmentFlag::Center );
    m_layout->addWidget( std::make_unique<Wt::WLabel>("low"), 2, 3, Wt::AlignmentFlag::Center );
    m_layout->addWidget( std::make_unique<Wt::WLabel>("medium"), 2, 4, Wt::AlignmentFlag::Center );
    m_layout->addWidget( std::make_unique<Wt::WLabel>("high"), 2, 5, Wt::AlignmentFlag::Center );

    {
      auto descOwner = std::make_unique<WContainerWidget>();
      WContainerWidget *desc = descOwner.get();
      desc->setOverflow(Wt::Overflow::Auto);
      m_description = desc->addNew<WText>();
      m_description->setWordWrap(true);
      m_layout->addWidget( std::move(descOwner), 1, 0, 4, 1 );
    }

    m_layout->addWidget( std::make_unique<WLabel>("Don't log:"), 3, 1, Wt::AlignmentFlag::Right );

    int i = 2;
    for( WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgLevel(0);
        level <= WarningWidget::WarningMsgHigh;
        level = WarningWidget::WarningMsgLevel(level+1) )
    {
      WCheckBox *warnToggle = m_layout->addWidget( std::make_unique<WCheckBox>(), 3, i++, Wt::AlignmentFlag::Center );

      const char *str = tostr( level );
      UserPreferences::associateWidget( str, warnToggle,  m_hostViewer );

      warnToggle->checked().connect( [this, level](){ setActivity( level, false ); } );
      warnToggle->unChecked().connect( [this, level](){ setActivity( level, true ); } );
    }//for( loop over

    m_layout->addWidget( std::make_unique<WLabel>("Hide notifications:"), 4, 1, Wt::AlignmentFlag::Right );

    i = 2;
    for( WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgLevel(0);
        level <= WarningWidget::WarningMsgHigh;
        level = WarningWidget::WarningMsgLevel(level+1) )
    {
      WCheckBox *warnToggle = m_layout->addWidget( std::make_unique<WCheckBox>(), 4, i++, Wt::AlignmentFlag::Center );

      const char *str  = popupToStr( level );
      UserPreferences::associateWidget( str, warnToggle,  m_hostViewer );

      warnToggle->checked().connect( [this, level](){ setPopupActivity( level, false ); } );
      warnToggle->unChecked().connect( [this, level](){ setPopupActivity( level, true ); } );
    }//for( loop over
  } //!m_layout //Create UI for first time
}//void WarningWidget::createContent()


//Displays the error message in m_description. This helps, as column may not be wide enough.
void WarningWidget::resultSelectionChanged()
{
  WModelIndexSet selected = m_tableView->selectedIndexes();
  if( selected.empty() )
  {
    m_description->setText( "" );
    return;
  }//if( selected.empty() )
  
  WStandardItem *item = m_messageModel->item( (*selected.begin()).row(), 2 );

  WString txt = item->text();
  m_description->setText( txt );
}//void resultSelectionChanged()


void WarningWidget::clearMessages()
{
  m_messageModel->removeRows(0, m_messageModel->rowCount());
}//void WarningWidget::clearMessages()


void WarningWidget::displayPopupMessageUnsafe( const Wt::WString &msg,
                                              WarningWidget::WarningMsgLevel level,
                                              int num_millies )
{
  if( num_millies <= 0 )
    num_millies = 5000;
  
  WStringStream strm;
  
  //We need the properly escaped string (e.g., not single quotes or unintended backslases), so we
  //  will use WString::jsStringLiteral(). Note that this will also place a single quote around the
  //  string.
  const string val = msg.jsStringLiteral( '\'' );
  string header, style;
  
  switch( level )
  {
    case WarningMsgInfo:
      header = "Info";
      style = "qtip-light";
    break;
      
    case WarningMsgLow:
      header = "Notice";
      style = "qtip-plain";
    break;
      
    case WarningMsgMedium:
      header = "Notice";
      style = "qtip-plain";
    break;
      
    case WarningMsgHigh:
      header = "Error";
      style = "qtip-red";
    break;
      
    case WarningMsgSave:
      header = "Save";
      style = "qtip-green";
    break;
      
    case WarningMsgShowOnBoardRiid:
#if( USE_REMOTE_RID )
      header = "On-board RIID Results";
#else
      header = "RIID Results";
#endif
      style = "qtip-blue qtip-riid";
      // remove all the other "show riid" results
      strm << "$('.qtip.qtip-riid').remove();";
      break;
    
#if( USE_REMOTE_RID )
    case WarningMsgExternalRiid:
      header = "External RID Results";
      style = "qtip-blue qtip-ext-riid";
      // remove all the other external RID results
      strm << "$('.qtip.qtip-ext-riid').remove();";
      break;
#endif
    
    case NumWarningMsgType:
      throw std::runtime_error( "Invalid warning message requested" );
  }//switch( level )
  
  //Create popup notifications
  strm << "var target = $('.qtip.jgrowl:visible:last'); $(document.body).qtip({ \
      content: { \
        title:  '<img style=\"vertical-align: middle; width: 16px;\" src=\\'" << string(iconUrl(level)) << "\\'/>&nbsp;&nbsp;" << header << "', \
        text:   "<< val <<", \
        button: true \
      }, \
      position: { \
        container: $('#qtip-growl-container') \
      }, \
      show: { \
        event: false, \
        ready: true, \
        delay: 500, \
        effect: function() { \
            $(this).stop(0, 1).animate({ height: 'toggle' }, 400, 'swing'); \
        }, \
      }, \
      hide: { \
        event:  false, \
        effect: function(api) { \
          $(this).stop(0, 1).animate({ height: 'toggle' }, 400, 'swing'); \
        } \
      }, \
      style: { \
        width: 250, \
        classes: 'jgrowl qtip-rounded qtip-shadow "<< style <<"', \
        tip: false \
      }, \
      events: { \
        render: function(event, api) { \
            if(!api.options.show.persistent) { \
                $(this).bind('mouseover mouseout', function(e) { \
                    var lifespan = " << num_millies << "; \
                    clearTimeout(api.timer); \
                    if (e.type !== 'mouseover') { \
                        api.timer = setTimeout(function() { api.hide(e) }, lifespan); \
                    } \
                }).triggerHandler('mouseout'); \
            } \
        } \
      } \
  }).removeData('qtip');";
  
  auto app = WApplication::instance();
  if( app )
    app->doJavaScript( strm.str() );
}//displayPopupMessageUnsafe(...)


void WarningWidget::addMessage( Wt::WString msg, int ilevel )
{
  if( ilevel < 0 || ilevel >= WarningMsgLevel::NumWarningMsgType )
    ilevel = WarningMsgHigh;

  const WarningMsgLevel level = WarningMsgLevel(ilevel);
  
  //Make sure the text is safe to display
  if( !Wt::Utils::removeScript(msg) )
    msg = Wt::Utils::htmlEncode( msg, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
  
  addMessageUnsafe( msg, level, 5000 );
} // void WarningWidget::addMessage(...)


void WarningWidget::addMessageUnsafe( const Wt::WString &msg,
                                     const WarningMsgLevel level,
                                     int num_millies )
{
  bool loggingActive = false, displayActive = false;
  
  switch( level )
  {
    case WarningMsgInfo:
    case WarningMsgLow:
    case WarningMsgMedium:
    case WarningMsgHigh:
    case WarningMsgSave:
    case WarningMsgShowOnBoardRiid:
#if( USE_REMOTE_RID )
    case WarningMsgExternalRiid:
#endif
      loggingActive =  m_active[level];
      displayActive = m_popupActive[level];
      break;
      
    case NumWarningMsgType:
      break;
  }//switch( level )
  
  if( loggingActive )
  {
    m_totalMessages++; //only count if logging
    
    vector<std::unique_ptr<WStandardItem>> message;
    message.push_back( std::make_unique<Wt::WStandardItem>( std::to_string(m_totalMessages) ) );
    message.push_back( std::make_unique<Wt::WStandardItem>( WarningWidget::description(level) ) );
    
    
    // Remove any script so we wont store that
    WString sanitized;
    
    // Wt::Utils::removeScript() will print out a message like
    //  "[secure] "XSS: discarding invalid attribute: onclick:..."
    // when it removes stuff - this is annoying, so lets manually remove
    // this attribute so I dont have to see this warning so much
    string msg_utf8 = msg.toUTF8();
    const auto onclick_pos = msg_utf8.find( "onclick=\"");
    if( onclick_pos != string::npos )
    {
      const auto end_pos = msg_utf8.find( "return false;\"", onclick_pos );
      if( end_pos != string::npos )
      {
        msg_utf8 = msg_utf8.substr(0,onclick_pos) + msg_utf8.substr(end_pos + 14);
        sanitized = WString::fromUTF8( msg_utf8 );
      }else
      {
        sanitized = msg;
      }
    }else
    {
      sanitized = msg;
    }//if( onclick_pos != string::npos ) / else
    
    if( !Wt::Utils::removeScript(sanitized) )
      sanitized = Wt::Utils::htmlEncode( sanitized, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
    
    if( (level == WarningMsgShowOnBoardRiid)
#if( USE_REMOTE_RID )
       || (level == WarningMsgExternalRiid)
#endif
       )
    {
      // Remove "buttons" we might be displaying
      string sanitized_utf8 = sanitized.toUTF8();
      const auto div_pos = sanitized_utf8.find("<div class=\"");
      if( div_pos != string::npos )
      {
        sanitized_utf8 = sanitized_utf8.substr( 0, div_pos );
        sanitized = WString::fromUTF8( sanitized_utf8 );
      }//if( we might have custom buttons on this message )
    }//if( a RIID message )
    
    message.push_back( std::make_unique<Wt::WStandardItem>( sanitized ) );

    m_messageModel->appendRow( std::move(message) );
  } // if( m_active[ level ] )
  
  if( displayActive )
    displayPopupMessageUnsafe( msg, level, 5000 );
}//addMessageUnsafe(...)


void WarningWidget::setActivity( WarningWidget::WarningMsgLevel priority, bool active )
{
  m_active[priority] = active;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo_redo = [priority,active]( const bool isUndo ){
      InterSpec *viewer = InterSpec::instance();
      WarningWidget *warn = viewer->warningWidget();
      if( !warn )
        return;
    
      try
      {
        const bool value = isUndo ? active : !active;
        UserPreferences::setBoolPreferenceValue( tostr(priority), value, viewer );
        warn->m_active[priority] = isUndo ? !active : active;
      }catch( std::exception &e )
      {
        Wt::log("error") << "Caught exception setting logging value for " << tostr(priority);
      }
    };//undo_redo
    
    auto undo = [undo_redo](){ undo_redo(true); };
    auto redo = [undo_redo](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Set logging preference" );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//WarningWidget::setActivity( WarningWidget::WarningMsgLevel priority, bool allowed )


void WarningWidget::setPopupActivity( WarningWidget::WarningMsgLevel priority, bool showPopup )
{
  m_popupActive[priority] = showPopup;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo_redo = [priority,showPopup]( const bool isUndo ){
      InterSpec *viewer = InterSpec::instance();
      WarningWidget *warn = viewer->warningWidget();
      if( !warn )
        return;
    
      try
      {
        const bool value = isUndo ? showPopup : !showPopup;
        UserPreferences::setBoolPreferenceValue( popupToStr(priority), value, viewer );
        warn->m_popupActive[priority] = isUndo ? !showPopup : showPopup;
      }catch( std::exception &e )
      {
        Wt::log("error") << "Caught exception setting logging value for " << popupToStr(priority);
      }
    };//undo_redo
    
    auto undo = [undo_redo](){ undo_redo(true); };
    auto redo = [undo_redo](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Set logging preference" );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//setPopupActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
