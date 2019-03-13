#ifndef HelpSystem_h
#define HelpSystem_h
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

#include <Wt/WObject>
#include <Wt/WMessageResourceBundle>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"


namespace Wt
{
  class WTree;
  class WTreeNode;
  class WLineEdit;
  class WPopupWidget;
  class WContainerWidget;
  class WMessageResourceBundle;
}//namespace Wt

namespace HelpSystem
{
  
  class HelpWindow: public AuxWindow
  {
  protected:
    /** Store a map from the ID of topic, to the tree node for that topic. */
    std::map<std::string, Wt::WTreeNode*> m_treeLookup;
    
    /** Store a map from the ID of topic, to the name of XML file that holds
        its info.
     */
    std::map<std::string, std::string> m_contentLookup;
    
    Wt::WContainerWidget* m_helpWindowContent;
    Wt::WTree* m_tree;
    Wt::WLineEdit *m_searchText;
    
    
    std::string m_helpLookupTable;
    
    void setPathVisible( Wt::WTreeNode *parent );
    
  public:
    HelpWindow(std::string preselect="");
    ~HelpWindow();
    void selectHelpToShow();
    void populateTree(Wt::Json::Array &res, Wt::WTreeNode* parent);
    void initialize();
    void handleArrowPress( const Wt::WKeyEvent e );
    void handleKeyPressInSearch( const Wt::WKeyEvent e );
  }; //class HelpWindow
  
  //createHelpWindow(...): just a convience function that creates a new
  //  HelpWindow instance.
  void createHelpWindow( const std::string &preselect = "" );
  
  //attachToolTipOn(...):  Creates a nice popup whenever the widget passed in
  //  gets focus or is moused over. The PopupToolBehavior argyment cooresponds
  //  to if this popup should obey the user preference to either show the popup
  //  instantly or delay showing it - or if the user prefference should be
  //  ignored and the popup aways instantly shown.
  //  The 'showInstantly' argument cooresponds to the users "ShowTooltips"
  //  preferance, and only matters if type==CanDelayShowing.
  enum PopupToolBehavior
  {
    AlwaysShow,
    CanDelayShowing
  };
  
  enum PopupPosition
  {
    Top,
    Bottom,
    Left,
    Right
  };
  void attachToolTipOn( Wt::WWebWidget* widget,
                       const std::string &text,
                       const bool showInstantly,
                       const PopupPosition pos = HelpSystem::Right,
                       const PopupToolBehavior override = HelpSystem::CanDelayShowing );
  
} //namespace HelpSystem
#endif //HelpSystem_h
