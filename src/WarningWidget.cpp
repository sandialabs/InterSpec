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

#include <sstream>
#include <iostream>
#include <algorithm>

#include <Wt/WLabel>
#include <Wt/WLink>
#include <Wt/WText>
#include <Wt/WTime>
#include <Wt/WImage>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WStandardItem>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>
#include <Wt/WCssDecorationStyle>

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/WarningWidget.h"
#include "SpecUtils/UtilityFunctions.h"
#if ( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#else
#include "InterSpec/SpectrumDisplayDiv.h"
#endif

using namespace Wt;
using namespace std;

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
  }//switch( i )
  
  throw runtime_error( "WarningWidget::description(...): invalid level" );
  return "";
}//const char *description( const WarningMsgLevel level )


WarningWidget::WarningWidget(
#if ( USE_SPECTRUM_CHART_D3 )
                             D3SpectrumDisplayDiv *spectrumDisplayDiv,
#else
                             SpectrumDisplayDiv *spectrumDisplayDiv,
#endif
                             InterSpec *hostViewer,
                             WContainerWidget *parent )
: WContainerWidget( parent ),
m_spectrumDisplayDiv( spectrumDisplayDiv ),
m_hostViewer( hostViewer ),
m_totalMessages(0),
m_app( wApp ),
m_layout(NULL),
m_messageModel(NULL),
m_tableView(NULL),
m_description(NULL)
{
  setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
  m_messageModel = new WStandardItemModel(1,4,this);
  
  // Find which messages should be active.
  for( WarningMsgLevel i = WarningMsgLevel(0);
       i <= WarningMsgHigh; i = WarningMsgLevel(i+1) )
    m_active[i] = !m_hostViewer->m_user->preferenceValue<bool>( tostr(i) );
  
  for( WarningMsgLevel i = WarningMsgLevel(0);
      i <= WarningMsgHigh; i = WarningMsgLevel(i+1) )
    m_popupActive[i] = !InterSpecUser::preferenceValue<bool>( WarningWidget::popupToStr(i), m_hostViewer );
  
  //Force WarningMsgSave to always be true
    m_active[WarningMsgSave]=true;
    m_popupActive[WarningMsgSave]=true;
    
  // Hook it up to the message handler in InterSpec
  hostViewer->messageLogged().connect( this, &WarningWidget::addMessage );
  
  //  addMessage( "Welcome to the InterSpec.",
  //           "", WarningWidget::WarningMsgInfo );
} // WarningWidget::WarningWidget( SpectrumDisplayDiv * )

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

void WarningWidget::addMessage( const Wt::WString &msg, const Wt::WString &src, int level )
{
  if( level < 0 || level > WarningMsgSave )
    level = WarningMsgHigh;

  string prefix;
  if( level <= WarningMsgInfo )
    prefix = "InterSpec_resources/images/information.png";
  else if( level <= WarningMsgLow )
    prefix = "InterSpec_resources/images/bullet_error.png";
  else if( level <= WarningMsgMedium )
    prefix = "InterSpec_resources/images/error.png";
  else if( level <= WarningMsgHigh )
    prefix = "InterSpec_resources/images/exclamation.png";
  else if( level <= WarningMsgSave )
    prefix = "InterSpec_resources/images/disk.png";
    
  
  if( m_active[ level ] )
  {
    m_totalMessages++; //only count if logging
    
    vector<WStandardItem*> message;
    Wt::WStandardItem * msgItem = new Wt::WStandardItem(std::to_string(m_totalMessages));
    message.push_back(msgItem);
    
    msgItem = new Wt::WStandardItem(prefix, WarningWidget::description(WarningWidget::WarningMsgLevel(level)));
    message.push_back(msgItem);
    
    msgItem = new Wt::WStandardItem(msg);
    message.push_back(msgItem);
    msgItem = new Wt::WStandardItem(src);
    message.push_back(msgItem);
    
    m_messageModel->appendRow(message);
  } // if( m_active[ level ] )
    
  if (m_popupActive[ level ] )
  {
//    //----------------------------------------------------------------------------
//    //need to escape the ' in the text message
//    string val = msg.toUTF8();
//    boost::replace_all(val, "'", "\\'");
//    
//    //Create popup notifications
//    std::stringstream strm;
//    if( level <= WarningMsgInfo )
//      strm <<    "new PNotify({text: '"<<val<<"', type:'info', delay:5000, opacity: 0.85});";
//    else if( level <= WarningMsgLow )
//      strm <<    "new PNotify({text: '"<<val<<"', type:'notice', delay:5000, opacity: 0.85});";
//    else if( level <= WarningMsgMedium )
//      strm <<    "new PNotify({text: '"<<val<<"', type:'notice', delay:5000, opacity: 0.85});";
//    else
//      strm <<    "new PNotify({text: '"<<val<<"', type:'error', delay:5000, opacity: 0.85});";
//    WApplication::instance()->doJavaScript(strm.str());
//    
    //----------------------------------------------------------------------------
    //need to escape the ' in the text message
    string val = msg.toUTF8();
    
    //Replace all single backslashes with a double backslash 
    //  -Assumes user has not already escaped the string
    UtilityFunctions::ireplace_all( val, "\\\\", "[%%%%%%%%]" );
    UtilityFunctions::ireplace_all( val, "\\", "\\\\" );
    UtilityFunctions::ireplace_all( val, "[%%%%%%%%]", "\\\\\\\\" );

    //Replace single quotes within the string - double should be fine to leave
    boost::replace_all(val, "'", "\\'");
    boost::replace_all(val, "\n", "<br />");

    string header;
    string style;
    
    //Create popup notifications
    std::stringstream strm;
    if( level <= WarningMsgInfo )
    {
      header = "Info";
      style = "qtip-light";
    }
    else if( level <= WarningMsgLow )
    {
      header = "Notice";
      style = "qtip-plain";
    }
    else if( level <= WarningMsgMedium )
    {
      header = "Notice";
      style = "qtip-plain";
    }
    else if( level <= WarningMsgHigh )
    {
      header = "Error";
      style = "qtip-red";
    }
    else if( level <= WarningMsgSave )
    {
        header = "Save";
        style = "qtip-green";
    }
    
    strm << "var target = $('.qtip.jgrowl:visible:last'); $(document.body).qtip({ \
      content: { \
        title:  '<img style=\"vertical-align: middle;\" src=\\' "<<prefix<<"\\'/>&nbsp;&nbsp;"<<header<<"', \
        text:   '"<<val<<"', \
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
        classes: 'jgrowl qtip-rounded qtip-shadow "<<style<<"', \
        tip: false \
      }, \
      events: { \
        render: function(event, api) { \
            if(!api.options.show.persistent) { \
                $(this).bind('mouseover mouseout', function(e) { \
                    var lifespan = 5000; \
                    clearTimeout(api.timer); \
                    if (e.type !== 'mouseover') { \
                        api.timer = setTimeout(function() { api.hide(e) }, lifespan); \
                    } \
                }).triggerHandler('mouseout'); \
            } \
        } \
      } \
  }).removeData('qtip');";
  
  WApplication::instance()->doJavaScript(strm.str());
  
  } //m_popupActive[ level ]
} // void WarningWidget::addMessage(...)

void WarningWidget::setActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
{
  m_active[ priority ] = allowed;
} // WarningWidget::setActivity( WarningWidget::WarningMsgLevel priority, bool allowed )

void WarningWidget::setPopupActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
{
  m_popupActive[ priority ] = allowed;
} //WarningWidget::setPopupActivity( WarningWidget::WarningMsgLevel priority, bool allowed )
