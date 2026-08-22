#ifndef ElectronNativeMenu_h
#define ElectronNativeMenu_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#if( USING_ELECTRON_NATIVE_MENU )

#include <string>

/* Maps the PopupDivMenu / PopupDivMenuItem classes onto Electron's native Menu / MenuItem, the same
 way target/macos/NativeMenu.h maps them onto NSMenu / NSMenuItem.  The two headers are deliberately
 parallel, so the `#if` arms in src/PopupDiv.cpp read alike.
 
 Two things are simpler here than on macOS.  The "native" side lives in the Electron main process,
 reached only by message, so a handle is just the Wt widget id() of the menu or item - nothing
 native is ever dereferenced from C++, and there is no Target bridge / atomics / teardown ordering
 to get right.  And a click comes back through ElectronUtils::handle_message_from_nodejs(), which
 has already resolved the session and taken the WApplication::UpdateLock, so there is no
 WServer::post to do either.
 
 One thing is harder: Electron cannot insert into, or remove from, a Menu that is already installed
 (electron/electron#527) - only enabled/visible/checked/label/toolTip can be changed in place.  So
 main.js keeps the menu template as its source of truth: the state setters below take effect
 immediately, while the structural calls mutate the template and schedule a coalesced
 Menu.buildFromTemplate().  Callers do not need to do anything to batch them.
 
 Everything here is a no-op unless InterSpecApp::isPrimaryWindowInstance() - sessions opened in an
 external browser keep the HTML menus.
 
 Limitations shared with the macOS implementation are listed in InterSpec/PopupDiv.h.
 */

class PopupDivMenu;
class InterSpecApp;
class PopupDivMenuItem;

namespace Wt
{
  class WCheckBox;
}

/** Adds a top-level menu to the menu bar, in call order; returns its handle (empty if this is not
 the primary window).  \sa setElectronMenuRole
 */
std::string addElectronMenu( PopupDivMenu *menu, const char *name );

/** Tags an app-level menu with the part it plays in the native menu bar, so the main process can
 attach the pieces the OS supplies rather than InterSpec:
   - "app"  : this is the macOS application menu - InterSpec's File menu doubles as it, exactly as
              in the macOS app build - and gets the Services/Hide/Quit block appended.
   - "edit" : gets the cut/copy/paste/select-all roles appended; without them those shortcuts stop
              working in text fields once Electron's default menu is gone.
   - "help" : is what the Window menu is positioned in front of.
 Matching on the menu's label instead would break as soon as the user picks a different language.
 */
void setElectronMenuRole( const std::string &menu, const char *role );

/** Turns an already-inserted item into a sub-menu; returns the sub-menu's handle.
 'parentItem' is the item the sub-menu hangs off of - state changes (hide/disable) address that id,
 while items added to the sub-menu address the returned one.
 Unlike the macOS implementation, sub-menus may be nested to any depth.
 */
std::string addElectronSubMenu( const std::string &parentMenu, PopupDivMenu *menu,
                                const char *name, const std::string &parentItem );

/** Inserts an item into 'menu'; a negative index appends.  Returns the item's handle, and registers
 it with the session so a click can be routed back to it.
 */
std::string insertElectronMenuItem( const std::string &menu, PopupDivMenuItem *item, int index );

/** Inserts a separator into 'menu'; a negative index appends.  'item' is the Wt id of the
 WMenuItem Wt created for the separator, and becomes its handle.
 */
void addElectronSeparatorAt( int index, const std::string &menu, const std::string &item );

/** Converts an already-inserted item into a checkbox item driven by 'cb'.  The macOS analogue
 (addOsxCheckableMenuItem) destroys and re-creates the native item; here it is changed in place.
 */
void makeElectronMenuItemCheckable( const std::string &item, Wt::WCheckBox *cb );

/** Removes an item (or separator) and drops its click registration. */
void removeElectronMenuItem( const std::string &item );

void setElectronMenuItemHidden( const std::string &item, bool hidden );
void setElectronMenuItemEnabled( const std::string &item, bool enabled );
void setElectronMenuItemChecked( const std::string &item, bool checked );
void setElectronMenuItemLabel( const std::string &item, const std::string &label );

/** Unlike the other setters this one cannot be applied to a live MenuItem (Electron makes
 `accelerator` read-only), so main.js handles it by re-building the menu.
 */
void setElectronMenuItemAccelerator( const std::string &item, const std::string &accelerator );
void addElectronMenuItemToolTip( const std::string &item, const char *tooltip );

/** Dispatches a menu activation coming back from main.js; called from
 ElectronUtils::handle_message_from_nodejs(), on the session thread, under its UpdateLock.
 'payload' is "<itemid>", or "<itemid>|0" / "<itemid>|1" for a checkbox item.
 Returns whether the item was found.
 */
bool handleElectronMenuActivation( InterSpecApp *app, const std::string &payload );

#endif //USING_ELECTRON_NATIVE_MENU

#endif //ElectronNativeMenu_h
