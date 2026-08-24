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

#include <string>

#include <Wt/WMenuItem.h>
#include <Wt/WPopupMenu.h>

class PopupDivMenu;
class PopupDivMenuItem;


/** Pushes a menu check box's current checked/enabled state to its native menu counterpart.

 `WCheckBox::setChecked(bool)`, `enable()` and `disable()` are not virtual and emit nothing, so
 changing a menu check box directly - rather than through its `checked()`/`unChecked()` signals -
 leaves a native menu item stale.  Call this after doing so.  Does nothing if `cb` is not inside a
 menu item, or for builds without native menus, so call sites need no `#if`.
 */
void syncNativeMenuCheckBox( Wt::WCheckBox *cb );

namespace Wt
{
  class WCheckBox;
  class WPushButton;
}


/** Creates an application-level PopupDivMenu (e.g., File, View, Tools, Help)
 and wires it to the given menu-bar button.

 On desktop: sets up stateless JS for instant display, hover-switching between
 menus, and keyboard navigation.  On macOS with USE_OSX_NATIVE_MENU: creates a
 native menu bar entry.  On phones: sets up slide-in animation with back/close.

 Lifetime is managed by the button (via button->addChild()); the menu is
 automatically destroyed when the button is destroyed.

 @param button  The menu-bar button that will trigger this menu. Must not be null.
 @return Non-owning pointer to the menu (lifetime managed by button).
 */
PopupDivMenu *makeAppLevelMenu( Wt::WPushButton *button );


/** Creates a transient/dropdown PopupDivMenu and wires it to the given button.

 The menu pops up on button click and auto-hides when the user clicks elsewhere
 or selects an item.

 Lifetime is managed by the button (via button->addChild()); the menu is
 automatically destroyed when the button is destroyed.

 @param button  The button that will trigger this popup menu. Must not be null.
 @return Non-owning pointer to the menu (lifetime managed by button).
 */
PopupDivMenu *makePopupMenu( Wt::WPushButton *button );


/** Like #makePopupMenu, but does NOT connect `button->clicked()` to `popup()`.

 For callers that already have their own `clicked()` handler on `button` (e.g. so the menu items
 can be built lazily on first click).  Those callers must call `menu->popup( button )` themselves;
 going through #makePopupMenu instead would pop the menu up twice per click.

 The menu is owned by `button` (via `WObject::addChild`), so it is destroyed with the button -
 never `delete` the returned pointer.
 */
PopupDivMenu *makeButtonOwnedPopupMenu( Wt::WPushButton *button );


/* We want the application-level menus (i.e., "File", "View", "Tools", and "Help")
to appear instantly when you click on the buttons, but by default Wt needs a
round-trip to C++ and back to JS, causing a slight (>100 ms), but perceptable delay.
So instead if `APP_MENU_STATELESS_FIX` is defined, we'll fix things up so the
menu will be shown using just JS, so the delay isnt perceptable.  However, because
of the design of the WMenu stuff, there are quite a bit of work-arounds needed so
the menus still behave well.
*/
#define APP_MENU_STATELESS_FIX 1

//Limitations shared by USE_OSX_NATIVE_MENU and USING_ELECTRON_NATIVE_MENU:
//  -disabling the closing of the menu when an item is selected isn't supported.
//  -CheckBox items must be created through passing a WCheckBox to
//   PopupDivMenu::addWidget(...); WMenuItem::setCheckable() is Wt-side only.
//  -Only WCheckBox is supported by PopupDivMenu::addWidget(...), adding any
//   other widget to, or as, a menu item is not supported.
//  -A WCheckBox changed from outside the menu (e.g. a preference toggled
//   elsewhere) does not push its new state to the native item; setChecked()
//   emits nothing we can hook.
//  -Items that rely on their anchor holding a URL (downloads, opening a new
//   window) do not work; see the note in target/macos/NativeMenu.h.
//  -XHTML labels (PopupDivMenuItem::makeTextXHTML()) render as plain text.
//Limitations specific to USE_OSX_NATIVE_MENU:
//  -sub-menus are only mirrored one level deep, below an app-level menu.
//  -changing an item's text after creation is not propagated; the Electron
//   path handles this through PopupDivMenuItem::setMenuText(...).
//  -implementation casts the objective-c pointers to void*
//   pointers, I tried doing some forward declartions using things similar to
//   'typedef struct objc_object NSMenu', but then ran into linking errors.

class PopupDivMenu : public Wt::WPopupMenu
{
public:
  enum MenuType
  {
    //AppLevelMenu: application-level menus (e.g., File, View, Tools, Help).
    //  On macOS with USE_OSX_NATIVE_MENU: gets native menu bar entry.
    //  On phones: slides in from the left.
    //  On desktop: uses stateless JS for instant display.
    AppLevelMenu,

    //TransientMenu: popup/dropdown menu, similar to WPopupMenu, with a few
    //  enhancements like auto-hide and z-index management.
    TransientMenu
  };//enum MenuType

  /** Construct a PopupDivMenu.  Does not wire up any button - call one of the
   setup methods after construction, or use makeAppLevelMenu()/makePopupMenu().
   */
  PopupDivMenu( const MenuType menutype = TransientMenu );

  virtual ~PopupDivMenu();

  /** Wire this menu as a transient popup on a button.
   Does NOT transfer ownership - caller manages lifetime (e.g., via addChild).
   Connects button click to popup() and sets up auto-hide and z-index fixing.
   */
  void setupAsTransientMenu( Wt::WPushButton *button );

  /** Wire this menu to a button using standard Wt4 setMenu() ownership.
   The button takes ownership of the menu.  Caller must NOT also addChild().
   Returns raw pointer to the menu (owned by button).
   */
  static PopupDivMenu *setupAsButtonOwnedMenu( std::unique_ptr<PopupDivMenu> menu,
                                                Wt::WPushButton *button );

  /** Wire this menu as a desktop app-level menu on a button.
   Sets up stateless JS click/hover/keyboard nav.
   Also handles USE_OSX_NATIVE_MENU if enabled.
   Caller manages lifetime (e.g., via addChild).
   */
  void setupAsAppMenu( Wt::WPushButton *button );

  /** Wire this menu as a mobile slide-in menu on a button.
   Sets up slide-in animation, back/close buttons, and overlay dismiss.
   Caller manages lifetime.
   */
  void setupAsMobileMenu( Wt::WPushButton *button );
  
    
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
  
  
#if( USING_ELECTRON_NATIVE_MENU )
  /** Tags this app-level menu with the part it plays in the native menu bar - "edit" or "help".
   \sa setElectronMenuRole
   */
  void setNativeMenuRole( const char *role );
  
#endif

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
  
  //Note: there used to be an `isHidden()` override here that returned true unconditionally when
  //  m_mobile, to keep WPopupMenu from hiding this menu as soon as a submenu was opened (we want the
  //  sub-menu to open over the current one).  `WPopupMenu::done()` starts with
  //  `if( isHidden() ) return;`, so that hack also silently disabled aboutToHide(), the menu-level
  //  triggered(), cancel(), and every isHidden() query on a phone menu - which leaked a popup menu
  //  per right-click, among other things.  `setHideOnSelect(false)` in the constructor now provides
  //  the sub-menu behaviour instead; we still close menus ourselves on item activation, via
  //  mobileHideMenuAndParents().


  //parentItem(): if this PopupDivMenu is a sub menu of another PopupDivMenu,
  //  and was created by calling addPopupMenuItem(...) on the parent, then
  //  parentItem() will return its cooresponding WMenuItem in its parent,
  //  otherwise NULL.  Note that there is no garuntee that pointer will be valid
  //  if the parent has been deleteed, or this widget reomved from it.
  Wt::WMenuItem *parentItem();

  void showMobile();

  void parentClicked();
  void parentMouseWentOver();
  void parentTouchStarted();

  bool isMobile() const;

  Wt::WPushButton *parentButton();
  
protected:
  void setupDesktopMenuStuff();
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

#if( USING_ELECTRON_NATIVE_MENU )
  /** Handle of this menu's counterpart in Electron's native Menu; empty if it has none.
   The handle is just this widget's id() - see target/electron/NativeMenu.h.
   */
  std::string m_electronMenu;
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
  
  /** Keyboard accelerator for this item, in Electron's format (e.g. "CmdOrCtrl+O").
   
   Only used when the item is mirrored into a native menu; ignored by the HTML menus, so call sites
   do not need to be guarded.  May be called before or after the item is added to a menu.
   */
  void setNativeAccelerator( const std::string &accelerator );
  const std::string &nativeAccelerator() const;
  
  /** Changes the item's text, keeping any native counterpart in sync.
   
   `WMenuItem::setText()` is not virtual, so this cannot be done by an override; use this instead of
   setText() for items whose label changes at runtime.  (The macOS native menu cannot follow label
   changes at all - see the limitations above.)
   */
  void setMenuText( const Wt::WString &text );
  
#if( USE_OSX_NATIVE_MENU )
  void *getNsMenuItem();
  virtual void setHidden( bool hidden,
                          const Wt::WAnimation &animation = Wt::WAnimation() );
  // Keep the native menu item's cached enabled state (read by validateMenuItem on the AppKit
  //  thread) in sync with this widget's enabled state.
  virtual void setDisabled( bool disabled ) override;
#endif
  
#if( USING_ELECTRON_NATIVE_MENU )
  /** Handle of this item's counterpart in Electron's native menu; empty if it has none. */
  const std::string &nativeMenuItemId() const;
  
  // Keep the native menu item in sync with this widget's hidden/enabled state.
  virtual void setHidden( bool hidden,
                          const Wt::WAnimation &animation = Wt::WAnimation() ) override;
  virtual void setDisabled( bool disabled ) override;
#endif
  
protected:
#if( USE_OSX_NATIVE_MENU )
  // Keep the native cache synchronized when an ancestor widget enables/disables this item.
  virtual void propagateSetEnabled( bool enabled ) override;

  void *m_nsmenu;
  void *m_nsmenuitem;
  // Opaque, thread-safe bridge retained by both this item and the native NSMenuItem. Keeping it
  //  separately lets the Wt session thread invalidate/update state without messaging AppKit.
  void *m_nsmenuitemtarget;
  friend class PopupDivMenu;
#endif
  
#if( USING_ELECTRON_NATIVE_MENU )
  virtual void propagateSetEnabled( bool enabled ) override;
  
  /** Handle of the menu this item was inserted into, and of the item itself (which is just this
   widget's id()); both empty if the item has no native counterpart.
   */
  std::string m_electronMenu;
  std::string m_electronItem;
  friend class PopupDivMenu;
#endif
  
  std::string m_accelerator;
};//class PopupDivMenuItem

#endif //PopupDiv_h
