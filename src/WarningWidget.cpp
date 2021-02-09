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

#include <Wt/WLink>
#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStandardItem>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RowStretchTreeView.h"

using namespace Wt;
using namespace std;

const char *iconUrl( const WarningWidget::WarningMsgLevel level )
{
  switch( level )
  {
    case WarningWidget::WarningMsgLevel::WarningMsgInfo:     return "InterSpec_resources/images/information.png";
    case WarningWidget::WarningMsgLevel::WarningMsgLow:      return "InterSpec_resources/images/bullet_error.png";
    case WarningWidget::WarningMsgLevel::WarningMsgMedium:   return "InterSpec_resources/images/error.png";
    case WarningWidget::WarningMsgLevel::WarningMsgHigh:     return "InterSpec_resources/images/exclamation.png";
    case WarningWidget::WarningMsgLevel::WarningMsgSave:     return "InterSpec_resources/images/disk.png";
    case WarningWidget::WarningMsgLevel::WarningMsgShowRiid: return "InterSpec_resources/images/search_results.svg";
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
    case WarningMsgLevel::WarningMsgShowRiid:
      return "BlockNotificationRiid";
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
    case WarningWidget::WarningMsgShowRiid:
      return "PopupBlockNotificationRiid";
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
    case WarningWidget::WarningMsgShowRiid:
      return "RIID";
    case WarningWidget::NumWarningMsgType:
      return "";
  }//switch( i )
  
  throw runtime_error( "WarningWidget::description(...): invalid level" );
  return "";
}//const char *description( const WarningMsgLevel level )


WarningWidget::WarningWidget( InterSpec *hostViewer,
                              WContainerWidget *parent )
: WContainerWidget( parent ),
  m_hostViewer( hostViewer ),
  m_totalMessages(0),
  m_layout(NULL),
  m_messageModel(NULL),
  m_tableView(NULL),
  m_description(NULL)
{
  setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
  m_messageModel = new WStandardItemModel(1,4,this);
  
  // Find which messages should be active.
  for( WarningMsgLevel i = WarningMsgLevel(0); i <= WarningMsgHigh; i = WarningMsgLevel(i+1) )
    m_active[i] = !m_hostViewer->m_user->preferenceValue<bool>( tostr(i) );
  
  for( WarningMsgLevel i = WarningMsgLevel(0); i <= WarningMsgHigh; i = WarningMsgLevel(i+1) )
    m_popupActive[i] = !InterSpecUser::preferenceValue<bool>( WarningWidget::popupToStr(i), m_hostViewer );
  
  //Force WarningMsgSave to always be true
  m_active[WarningMsgSave] = true;
  m_popupActive[WarningMsgSave] = true;
  
  //Force ShowRiid to always be true
  m_active[WarningMsgShowRiid] = true;
  m_popupActive[WarningMsgShowRiid] = true;
  
  // Hook it up to the message handler in InterSpec
  hostViewer->messageLogged().connect( this, &WarningWidget::addMessage );
}//WarningWidget::WarningWidget( SpectrumDisplayDiv * )


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
    m_layout = new WGridLayout();
    m_layout->setContentsMargins(0,0,0,0);
    setLayout(m_layout);
    
    setPadding(0);
    setMargin(0);
    m_messageModel->setItem (0,0, new WStandardItem("0"));
    m_messageModel->setItem (0,1, new WStandardItem(""));
    m_messageModel->setItem (0,2, new WStandardItem("InterSpec started."));
    m_messageModel->setItem (0,3, new WStandardItem(""));
    m_messageModel->setHeaderData(  0, Horizontal, WString("#"), DisplayRole );
    m_messageModel->setHeaderData(  1, Horizontal, WString("Type"), DisplayRole );
    m_messageModel->setHeaderData(  2, Horizontal, WString("Message"), DisplayRole );
    m_messageModel->setHeaderData(  3, Horizontal, WString("Source"), DisplayRole );
    
    m_tableView = new RowStretchTreeView();
    m_tableView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
    m_tableView->setModel(m_messageModel);
    m_tableView->setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
    m_tableView->setColumnResizeEnabled(true);
    m_tableView->setColumnAlignment(0, Wt::AlignCenter);
    m_tableView->setHeaderAlignment(0, Wt::AlignCenter);
    m_tableView->setColumnWidth(0,50);
    m_tableView->setColumnWidth(1,100);
    m_tableView->setColumnWidth(2,400);
    m_tableView->setColumnWidth(3,200);
    
    m_tableView->setAlternatingRowColors(true);
    m_tableView->setSelectionMode(Wt::NoSelection);
    m_tableView->setEditTriggers(Wt::WAbstractItemView::NoEditTrigger);

    m_tableView->setSelectionMode( Wt::SingleSelection );
    m_tableView->setSelectionBehavior( Wt::SelectRows );
    m_tableView->selectionChanged().connect( this, &WarningWidget::resultSelectionChanged );
      
    m_layout->addWidget(m_tableView,0,0,1,6);
    m_layout->setRowStretch(0,1);
    m_layout->setColumnStretch(0,1);
    
    Wt::Dbo::ptr<InterSpecUser> m_user = m_hostViewer->m_user;
    
    
    WImage* image = new Wt::WImage(Wt::WLink("InterSpec_resources/images/information.png"));
    image->setMaximumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    image->setMinimumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    m_layout->addWidget(image,1,2, AlignCenter);

    image = new Wt::WImage(Wt::WLink("InterSpec_resources/images/bullet_error.png"));
    image->setMaximumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    image->setMinimumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    m_layout->addWidget(image,1,3, AlignCenter);
    
    image = new Wt::WImage(Wt::WLink("InterSpec_resources/images/error.png"));
    image->setMaximumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    image->setMinimumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    m_layout->addWidget(image,1,4, AlignCenter);
    
    image = new Wt::WImage(Wt::WLink("InterSpec_resources/images/exclamation.png"));
    image->setMaximumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    image->setMinimumSize(WLength(16,WLength::Pixel), WLength(16,WLength::Pixel));
    m_layout->addWidget(image,1,5, AlignCenter);
    
    m_layout->addWidget(new Wt::WLabel("info"),2,2, AlignCenter);
    m_layout->addWidget(new Wt::WLabel("low"),2,3, AlignCenter);
    m_layout->addWidget(new Wt::WLabel("medium"),2,4, AlignCenter);
    m_layout->addWidget(new Wt::WLabel("high"),2,5, AlignCenter);
    
    WContainerWidget* desc = new WContainerWidget();
    desc->setOverflow(Wt::WContainerWidget::OverflowAuto);
    m_description = new WText(desc);
    m_description->setWordWrap(true);
    m_layout->addWidget(desc,1,0, 4,1);
    WLabel* label = new WLabel("Record log types:");
    label->setImage(new WImage( "InterSpec_resources/images/database_edit.png", this ));
    m_layout->addWidget(label,3,1, AlignRight);
    
    int i = 2;
    for( WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgLevel(0);
        level <= WarningWidget::WarningMsgHigh;
        level = WarningWidget::WarningMsgLevel(level+1) )
    {
      WCheckBox *warnToggle = new WCheckBox();
      m_layout->addWidget(warnToggle,3,i++, AlignCenter);
      
      const char *str = tostr( level );
      InterSpecUser::associateWidget( m_user, str, warnToggle, m_hostViewer, true );
      
      warnToggle->checked().connect( boost::bind( &WarningWidget::setActivity, this,  level, true ) );
      warnToggle->unChecked().connect( boost::bind( &WarningWidget::setActivity, this,  level, false ) );
    }//for( loop over
    
    label = new WLabel("Show popup notification types:");
    label->setImage(new WImage( "InterSpec_resources/images/note.png", this ));
    
    m_layout->addWidget(label,4,1, AlignRight);
    
    i = 2;
    for( WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgLevel(0);
        level <= WarningWidget::WarningMsgHigh;
        level = WarningWidget::WarningMsgLevel(level+1) )
    {
      WCheckBox *warnToggle = new Wt::WCheckBox();
      m_layout->addWidget(warnToggle,4,i++, AlignCenter);
      
      const char *str  = popupToStr( level );
      InterSpecUser::associateWidget( m_user, str, warnToggle, m_hostViewer , true);
      
      warnToggle->checked().connect( boost::bind( &WarningWidget::setPopupActivity, this,  level, true ) );
      warnToggle->unChecked().connect( boost::bind( &WarningWidget::setPopupActivity, this,  level, false ) );
    }//for( loop over
  } //!m_layout //Create UI for first time
}//void WarningWidget::createContent()


//Displays the error message in m_description. This helps, as column may not be wide enough.
void WarningWidget::resultSelectionChanged()
{
    WModelIndexSet selected = m_tableView->selectedIndexes();
    if( selected.empty() )
    {
        return;
    }//if( selected.empty() )

    WStandardItem * item = m_messageModel->item(    (*selected.begin()).row(),2);

    m_description->setText(item->text());
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
      
    case WarningMsgShowRiid:
      header = "RIID Results";
      style = "qtip-blue qtip-riid";
      // remove all the other "show riid" results
      strm << "$('.qtip.qtip-riid').remove();";
      break;
    
    case NumWarningMsgType:
      throw std::runtime_error( "Invalid warning message requested" );
  }//switch( level )
  
  //Create popup notifications
  strm << "var target = $('.qtip.jgrowl:visible:last'); $(document.body).qtip({ \
      content: { \
        title:  '<img style=\"vertical-align: middle;\" src=\\'" << string(iconUrl(level)) << "\\'/>&nbsp;&nbsp;" << header << "', \
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




void WarningWidget::addMessage( Wt::WString msg, Wt::WString src, int ilevel )
{
  if( ilevel < 0 || ilevel > WarningMsgSave )
    ilevel = WarningMsgHigh;

  const WarningMsgLevel level = WarningMsgLevel(ilevel);
  
  //Make sure the text is safe to display
  if( !Wt::Utils::removeScript(msg) )
    msg = Wt::Utils::htmlEncode( msg, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
  
  if( m_active[ level ] )
  {
    m_totalMessages++; //only count if logging
    
    vector<WStandardItem*> message;
    Wt::WStandardItem * msgItem = new Wt::WStandardItem(std::to_string(m_totalMessages));
    message.push_back(msgItem);
    
    msgItem = new Wt::WStandardItem( string(iconUrl(level)), WarningWidget::description(level));
    message.push_back(msgItem);
    
    msgItem = new Wt::WStandardItem(msg);
    message.push_back(msgItem);
    msgItem = new Wt::WStandardItem(src);
    message.push_back(msgItem);
    
    m_messageModel->appendRow(message);
  } // if( m_active[ level ] )
    
  if( m_popupActive[level] )
    displayPopupMessageUnsafe( msg, level, 5000 );
} // void WarningWidget::addMessage(...)

void WarningWidget::setActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
{
  m_active[ priority ] = allowed;
} // WarningWidget::setActivity( WarningWidget::WarningMsgLevel priority, bool allowed )

void WarningWidget::setPopupActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
{
  m_popupActive[ priority ] = allowed;
} //WarningWidget::setPopupActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
