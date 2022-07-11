#ifndef WarningWidget_h
#define WarningWidget_h
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

#include <Wt/WContainerWidget>

// Some forward declarations
class InterSpec;
class RowStretchTreeView;
namespace Wt
{
  class WText;
  class WGridLayout;
  class WStandardItemModel;
}//namespace Wt


/// \TODO: separate this class into kinda the model that accepts error messages, flashes them to the
///        user, and keeps track of them, and the display (e.g., have #createContent() create a
///        view for this class), and then improve the memory management
class WarningWidget : public Wt::WContainerWidget
{
public:
  enum WarningMsgLevel
  {
    WarningMsgInfo = 0,
    WarningMsgLow = 1,
    WarningMsgMedium = 2,
    WarningMsgHigh = 3,
    WarningMsgSave = 4,
    WarningMsgShowOnBoardRiid = 5,
#if( USE_REMOTE_RID )
    WarningMsgExternalRiid = 6,
#endif
    
    NumWarningMsgType
  };//enum WarningMsgLevel
  
  static const char *tostr( const WarningMsgLevel level );
  static const char *popupToStr( const WarningMsgLevel level );
  static const char *description( const WarningMsgLevel level );
  
public:
  WarningWidget( InterSpec *hostViewer,
                 Wt::WContainerWidget *parent = 0 );
  virtual ~WarningWidget();

  void createContent();
  
  /** Displays message to the user in a q-tip popup if the user option is to display messages at
   that level.  Also stores the message in the model for later viewing.
   
   Input 'msg' will be sanitized before displaying (e.g., must be valid XHTML and not have any JS
   scripting) and storage.
   */
  void addMessage( Wt::WString msg, int level );

  /** Similar to #addMessage, but \c msg contents will not be sanitized before displaying (so
   see warnings from the #displayPopupMessageUnsafe function), and you can control how long the
   message will be displayed for.
   */
  void addMessageUnsafe( const Wt::WString &msg, const WarningMsgLevel level, int num_millies );
  
  /** Displays the message to the user in the q-tip popup style.
   
   Does not check the format or content of the msg, so you need to make sure any potentially bad
   script, or non-XHTML-ness is removed first.
   
   Also, does not include the message into the application log (e.g., m_messageModel)
   
   @param num_millies The number of milliseconds to have the popup stick around for.
   */
  static void displayPopupMessageUnsafe( const Wt::WString &msg, const WarningMsgLevel level,
                                         int num_millies );
  
  
  void setPopupActivity( WarningMsgLevel priority, bool allowed );
  void setActivity( WarningMsgLevel priority, bool allowed );

  bool active( WarningMsgLevel level ) const;

  void resultSelectionChanged();
    
  void clearMessages();
  
protected:  
  InterSpec *m_hostViewer;

  int m_totalMessages;

  bool m_popupActive[int(WarningMsgLevel::NumWarningMsgType)];
  bool m_active[int(WarningMsgLevel::NumWarningMsgType)];

  Wt::WGridLayout *m_layout;
  Wt::WStandardItemModel *m_messageModel;
    
  RowStretchTreeView *m_tableView;
  Wt::WText *m_description;
}; //class WarningWidget

#endif
