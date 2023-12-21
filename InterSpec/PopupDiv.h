#ifndef PopupDiv_h
#define PopupDiv_h
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

#include <Wt/WMenuItem>
#include <Wt/WPopupMenu>

class PopupDivMenuItem;

namespace Wt
{
  class WCheckBox;
  class WPushButton;
}

/* We want the application-level menus (i.e., "File", "View", "Tools", and "Help")
to appear instantly when you click on the buttons, but by default Wt needs a 
round-trip to C++ and back to JS, causing a slight (>100 ms), but perceptable delay.
So instead if `APP_MENU_STATELESS_FIX` is defined, we'll fix things up so the
menu will be shown using just JS, so the delay isnt perceptable.  However, because
of the design of the WMenu stuff, there are quite a bit of work-arounds needed so
the menus still behave well.
*/
#define APP_MENU_STATELESS_FIX 1

//Limitations of enabling USE_OSX_NATIVE_MENU:
//  -disabling the closing of the menu when an item is selected isn't supported.
//  -CheckBox items must be created through passing a WCheckBox to
//   PopupDivMenu::addWidget(...).
//  -Only WCheckBox is supported by PopupDivMenu::addWidget(...), adding any
//   other widget to, or as, a menu item is not supported.
//  -USE_OSX_NATIVE_MENU implementation casts the objective-c pointers to void*
//   pointers, I tried doing some forward declartions using things similar to
//   'typedef struct objc_object NSMenu', but then ran into linking errors.

class PopupDivMenu : public Wt::WPopupMenu
{
public:
  enum MenuType
  {
    //AppLevelMenu says the menu is for an applications level menu that should
    //  be customized for the device.  E.g. for OS X is USE_OSX_NATIVE_MENU is
    //  enabled, then this menu will get an entry on the OS X native toolbar.
    //  For phones it will slide in from the left.
    //For desktop style menus, you currently must still pass in a valid
    //  WPushButton parent in order for the menu to be placed in the applicaiton
    //  level menu role on OS X.  I'm not crazy about this somewhat ambigous
    //  initialization, so how these menus are created should be re-visited.
    AppLevelMenu,
    
    //TransientMenu: menu behaves as a popup, similar to WPopupMenu, but with
    //  a few enhancments like poping up on mouse over the menuParent (if
    //  specified).
    //  XXX - specialization for phone not yet handled
    TransientMenu
  };//enum MenuType
  
  //PopupDivMenu constructor: if a menuParent is non null, then clicking on, or
  //  mouseing over the button will cause activation of this PopupDivMenu.
  PopupDivMenu( Wt::WPushButton *menuParent, const MenuType menutype );
  
  //~PopupDivMenu(): currently a no-op function
  virtual ~PopupDivMenu();
  
    
#if( APP_MENU_STATELESS_FIX )
  static void pre_render(PopupDivMenu* menu);
#endif

  //Add separator; returns pointer to the sperator
  Wt::WMenuItem *addSeparator();
  
  /** Add a sepeartor at the specified index; if less than zero than adds to the
   end.
   */
  Wt::WMenuItem *addSeparatorAt( int index );
  
  /** For Electron and macOS native menues, we need to do some special handling
     to remove seperators.
   This function call does not delete the passed in seperator, just removes it
   from the menu.
   Returns true if the seperator was found and removed, false otherwise.
   */
  bool removeSeperator( Wt::WMenuItem *sepertor );
  
  
  virtual void setHidden( bool hidden,
                          const Wt::WAnimation &animation = Wt::WAnimation() );
  
  //addWidget(...): adds an arbitrary widget to the menu.  If the item is
  //  anything but a WCheckBox then the "PopupDivMenuWidget" class will be added
  //  to its styling so that mouse over events will not cause the whole widget
  //  to become highlighted.  Children of widget are also not currently
  //  highlighted when moused over
  PopupDivMenuItem *addWidget( Wt::WWidget *widget,
                               const bool closeOnClickInWidget = false );
  
  //addMenuItem(...): adds an additional menu item with desired text.
  PopupDivMenuItem *addMenuItem( const Wt::WString &text,
                                 const std::string &iconPath = "",
                                 const bool closeMenuOnTriggered = true);
  
  /** Inserts a menu item at the specified index.  An index less than zero will
      cause the item to be inserted as the last item.
   */
  virtual PopupDivMenuItem *insertMenuItem( const int index,
                                        const Wt::WString &text,
                                        const std::string &iconPath,
                                        const bool closeMenuOnActivation );
  
  //addPopupMenuItem(...): a better name would be addSubMenu(...)
  PopupDivMenu *addPopupMenuItem( const Wt::WString &text,
                                  const std::string &iconPath = "" );
  
  //addPhoneBackItem(...): if parent is specified, then only this menu will be
  //  be closed, but parent menu wont be.
  PopupDivMenuItem *addPhoneBackItem( PopupDivMenu *parent );
  
  //isHidden(): for m_mobile==true, this functions always returns true.  Else
  //  it returns WPopupMenu::isHidden().
  //  This is a hack to keep WPopupMenu from hiding this menu as soon as
  //  a submenu is opened (we want the sub-menu to open over the current one).
  //  This does create the issue where we have to manually close the menues
  //  on item activation though.
  virtual bool isHidden() const;
  
  
  //parentItem(): if this PopupDivMenu is a sub menu of another PopupDivMenu,
  //  and was created by calling addPopupMenuItem(...) on the parent, then
  //  parentItem() will return its cooresponding WMenuItem in its parent,
  //  otherwise NULL.  Note that there is no garuntee that pointer will be valid
  //  if the parent has been deleteed, or this widget reomved from it.
  Wt::WMenuItem *parentItem();

  //
  void showMobile();
  
  
  void setupDesktopMenuStuff();
  
  void parentClicked();

  void parentMouseWentOver();
  
  void parentTouchStarted();
  
  bool isMobile() const;
  
  Wt::WPushButton *parentButton();
  
protected:
  void desktopDoHide();

  void mobileDoHide();
  void mobileHideMenuAndParents();
  
  Wt::WMenuItem *m_parentItem;
  Wt::WPushButton *m_menuParent;
  std::string m_menuParentID;
  
#if( USE_OSX_NATIVE_MENU )
  void *m_nsmenu;
  friend class PopupDivMenuItem;
#endif

  bool m_mobile;
  const MenuType m_type;
};//class PopupDivMenu


//PopupDivMenuItem allows an abstraction over WMenuItem to allow using the
//  macOs or Electron native menu systems, as well as address a few minor
//  behavior preferences we have.
class PopupDivMenuItem : public Wt::WMenuItem
{
public:
  //Constructor: 'text' is what the menu item will say. 'iconPath' is the
  //  path to the icon file; the icon should have a width of 16 pixels.
  //  Note that text is PlainText by default.
  PopupDivMenuItem( const Wt::WString &text, const std::string &iconPath  );
  
  virtual ~PopupDivMenuItem();
  
  //checkBox(): returns null if checkable isnt set, or you didnt add a WCheckBox
  //  throough PopupDivMenu::addWidget(...), otherwise returns the WCheckBox
  //  pointer created by Wt if you called WMenuItem::setCheckable(), or the
  //  WCheckBox you passed into PopupDivMenu::addWidget(...).
  //  If you called setCheckable(), will throw exception if the check box cant
  //  be retrieved (this probably wont ever happen).
  Wt::WCheckBox *checkBox();

  //anchor(): returns the anchor of the menu item; will return NULL if this
  //  doesnt exist (it almost always should).
  Wt::WAnchor *anchor();
  
  //makeTextXHTML(): uses a hack to have Wt display the text as XHTML instead
  //  of the default PlainText.  There is no need to re-set the text.
  void makeTextXHTML();
  
  //nonAnchorClickHack(): the trigger() signal only gets emitted if the anchor
  //  in the menu item (of a normal item with text) gets clicked; clicking the
  //  area outside of anchor will not cause the trigger() signal to be emitted.
  //  To fix this, we will make it so the clicked() signal of the anchor will
  //  not be propogated, and if this->clicked() is emmited (which now wont
  //  happen if the anchor is clicked), we will call WMenuItem::select().
  void nonAnchorClickHack();
  
#if( USE_OSX_NATIVE_MENU )
  void *getNsMenuItem();
  virtual void setHidden( bool hidden,
                          const Wt::WAnimation &animation = Wt::WAnimation() );
#endif
  
protected:
#if( USE_OSX_NATIVE_MENU )
  void *m_nsmenu;
  void *m_nsmenuitem;
  friend class PopupDivMenu;
#endif
};//class PopupDivMenuItem

#endif //PopupDiv_h
