#ifndef UseInfoWindow_h
#define UseInfoWindow_h
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

#include <string>
#include <functional>

#include <Wt/WString>
#include <Wt/WMenuItem>
#include <Wt/WContainerWidget>
#include <Wt/WMessageResourceBundle>

#include "InterSpec/AuxWindow.h"

namespace Wt
{
  class WMenu;
  class WStandardItemModel;
}

class InterSpec;
class RowStretchTreeView;
#if( USE_DB_TO_STORE_SPECTRA )
class SnapshotBrowser;
#endif

namespace DataBaseUtils{ class DbSession; }

class SideMenuItem : public Wt::WMenuItem
{
public:
  SideMenuItem( const Wt::WString &text, Wt::WWidget *contents );
  
protected:
  virtual Wt::SignalBase& activateSignal();
};//class SideMenuItem


class UseInfoWindow : public AuxWindow
{
  //UseInfoWindow:
  //  When user clicks closed, or finished, the window will be hidden and the
  //  finsished signal emitted, but not deleted
public:
  UseInfoWindow(std::function<void(bool)> showAgainCallback, InterSpec* viewer );
  ~UseInfoWindow();
  
protected:
  struct VideoInformation
  {
    std::string key, fileMP4, fileOGV;
    Wt::WString title;
  };
  
protected:

  void initVideoInfoMap();
  void itemCreator( const std::string &resource, Wt::WContainerWidget *parent, Wt::WString title );
  void right_select_item(Wt::WMenuItem *item );
  void tab_select_item(Wt::WMenuItem *item);
  SideMenuItem * makeItem( const Wt::WString &title, const std::string &resource);
  
  /** Lazily makes a menu-item, whose contents are from a  #Wt::WMessageResourceBundle file.
   
   @param title The name to make on the side-menu button
   @param resource The filesystem path of the XML file that contains the contents for a #Wt::WMessageResourceBundle
   */
  SideMenuItem *makeTextItem( const Wt::WString &title, const std::string &resource );
  
  /** Creates a widget with the contents of the resource file.
   
   This function is called when the user clicks on the side-menu item created by the `makeTextItem` function.
   Each <message></message> in the resource file will be its own section.
   
   @param resource The filesystem path, relative to Wts docroot, to find the resource XML and CSS resources; this
          path should not have file extension, as the localized ending (e.x., "en-US") will attempt to be added, as well as
          ".xml", and if ".css" exists, it will be used/loaded as well.
   @param parent The parent of the WText widget that gets created.
   */
  void textItemCreator( const std::string &resource, Wt::WContainerWidget *parent );
  
public:
    void loadSample( const Wt::WModelIndex index );
    void loadSampleSelected( );
    void handleSampleSelectionChanged();
    void handleSampleDoubleClicked( Wt::WModelIndex index, Wt::WMouseEvent event );
    
  
  void showFaqTab();
protected:
  
  std::shared_ptr<DataBaseUtils::DbSession>           m_session;
  
#if( USE_DB_TO_STORE_SPECTRA )
  SnapshotBrowser                                     *m_snapshotModel;
#endif
  
  //Sample
  RowStretchTreeView                                  *m_tableSample;
  Wt::WStandardItemModel                              *m_messageModelSample;
  
  InterSpec                                           *m_viewer;

  std::map<std::string, Wt::WMediaPlayer *>            m_players;
  Wt::WMessageResourceBundle                           m_resourceBundle;
  Wt::WMenu                                           *m_menu;
  std::map<std::string,std::vector<VideoInformation> > m_videoInfos;
  
  //m_videoTab isnt currently being but its use is left commented out everywhere
  //  incase I ever get around to making videos that would be useful...
  //Wt::WTabWidget                                      *m_videoTab;
}; //class UseInfoWindow : public AuxWindow


//DeferredWidget: adapted from Wt Widget Gallery example
template <typename Function>
class DeferredWidget : public Wt::WContainerWidget
{
public:
  DeferredWidget( Function f ) : m_loaded( false ), m_f( f ) {}
  
  void docreate()
  {
    if( !m_loaded )
    {
      m_f( this );
      
      // On phone "Welcome" screen, at least, lazily loading content that needs to scroll, wont
      //  trigger the layout code give the loaded widget a height, so the widget wont be scrollable,
      //  so we'll do this
      wApp->doJavaScript( wApp->javaScriptClass() + ".TriggerResizeEvent();" );
    }//if( !m_loaded )
    m_loaded = true;
  }//void docreate()
  
private:
  void load()
  {
    docreate();
    WContainerWidget::load();
  }//void load()
  
  bool m_loaded;
  Function m_f;
};//class DeferredWidget


template <typename Function>
DeferredWidget<Function> *deferCreate( Function f )
{
  return new DeferredWidget<Function>( f );
}


#endif
