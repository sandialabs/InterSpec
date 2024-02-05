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

#include <initializer_list>

#include <Wt/WObject>
#include <Wt/WMessageResourceBundle>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"


namespace Wt
{
  class WTree;
  class WTemplate;
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
    //ToDo: a Boost.BiMap would be a lot better fit
    std::map<std::string, Wt::WTreeNode*> m_treeLookup;
    
    /** Store a map from the ID of topic, to the name of XML file that holds
        its info.
     */
    std::map<std::string, std::string> m_contentLookup;
    
    Wt::WContainerWidget* m_helpWindowContent;
    Wt::WTree* m_tree;
    Wt::WLineEdit *m_searchText;
    std::string m_displayedTopicName;
    
    std::string m_helpLookupTable;
    
    void setPathVisible( Wt::WTreeNode *parent );
    
  public:
    HelpWindow(std::string preselect="");
    ~HelpWindow();
    void selectHelpToShow();
    Wt::WTemplate *getContentToDisplay( const std::string &tag ) const;
    
    /** If there are not XML help files associated with a tag, but there are
        children to that tag, then the contents of all the (no-recursive)
        children will be returned.
     */
    std::string getHelpContents( const std::string &tag ) const;
    
    void populateTree(Wt::Json::Array &res, Wt::WTreeNode* parent);
    void initialize();
    void handleArrowPress( const Wt::WKeyEvent e );
    void handleKeyPressInSearch( const Wt::WKeyEvent e );
    
    void setTopic( const std::string &preselect );
    const std::string &currentTopic() const;
  }; //class HelpWindow
  
  //createHelpWindow(...): just a convience function that creates a new
  //  HelpWindow instance.
  void createHelpWindow( const std::string &preselect = "" );
  

  enum class ToolTipPrefOverride
  {
    AlwaysShow,
    RespectPreference
  };
  
  enum class ToolTipPosition
  {
    Top,
    Bottom,
    Left,
    Right
  };

  /** Creates a nice popup when the widget passed into this function gets focus, or is moused over.
   
   \param widget The widget to add the tooltip to
   \param text The xhtml text to place in the tooltip.
   \param enableShowing Should correspond to the users "ShowTooltips" preference.  If false, tooltip wont be trigger-able.
   \param pos The position, relative to the widget, to place the tooltip
   \param forceShowing If, for this particular widget, the tooltip should always be shown when moused-over. despite what
          enableShowing dictates
   
   Note: even if tooltip is not enabled, we still load it to the DOM, incase the user enables tooltips later - we could probably avoid
   cluttering up the DOM by caching server-side, but we'll leave that for the future.
   
   Note: when the widget is deleted from memory, server-side, the tooltip will be removed from the DOM client-side.
   */
  void attachToolTipOn( Wt::WWebWidget* widget,
                       const std::string &text,
                       const bool enableShowing,
                       const ToolTipPosition pos = HelpSystem::ToolTipPosition::Right,
                       const ToolTipPrefOverride forceShowing
                                            = HelpSystem::ToolTipPrefOverride::RespectPreference );
  
  /** Same as above, but attaches the same tooltip to multiple elements.
   
   The lifetime of the tooltip is governed by the first element in the list.
   */
  void attachToolTipOn( std::initializer_list<Wt::WWebWidget*> widgets,
                        const std::string &text,
                        const bool enableShowing,
                        const ToolTipPosition pos = HelpSystem::ToolTipPosition::Right,
                        const ToolTipPrefOverride forceShowing
                                  = HelpSystem::ToolTipPrefOverride::RespectPreference );
  
} //namespace HelpSystem
#endif //HelpSystem_h
