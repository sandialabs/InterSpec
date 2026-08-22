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
#include <iostream>

#include <nlohmann/json.hpp>

#include <Wt/WString.h>
#include <Wt/WCheckBox.h>
#include <Wt/WMenuItem.h>
#include <Wt/WApplication.h>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpecApp.h"

#include "target/electron/NativeMenu.h"
#include "target/electron/ElectronUtils.h"

using namespace std;

namespace
{
/** Every function in this file is a no-op for sessions that are not the primary application
 window - i.e. sessions opened in an external browser, which keep the HTML menus.  This is the
 same guard addOsxMenu() applies.
 */
bool using_native_menu()
{
  return InterSpecApp::isPrimaryWindowInstance();
}//bool using_native_menu()


void send( const char *msg_name, const nlohmann::json &data )
{
  ElectronUtils::send_nodejs_message( msg_name, data.dump() );
}//void send(...)


/** Some menu labels are padded with a leading space to line the HTML check boxes up; a native menu
 does its own spacing, so that padding just looks like a stray indent.
 */
std::string native_label( const Wt::WString &text )
{
  return SpecUtils::trim_copy( text.toUTF8() );
}//std::string native_label(...)


/** Sends a single-property change for an item.  main.js applies these to the live MenuItem, so no
 menu rebuild happens.
 */
template<class T>
void set_item_state( const string &item, const char *property, const T &value )
{
  if( item.empty() || !using_native_menu() )
    return;
  
  nlohmann::json msg;
  msg["id"] = item;
  msg[property] = value;
  send( "MenuSetItemState", msg );
}//void set_item_state(...)
}//namespace


std::string addElectronMenu( PopupDivMenu *menu, const char *name )
{
  if( !menu || !name || !name[0] || !using_native_menu() )
    return string();
  
  nlohmann::json msg;
  msg["id"] = menu->id();
  msg["label"] = name;
  send( "MenuAddMenu", msg );
  
  return menu->id();
}//std::string addElectronMenu(...)


void setElectronMenuRole( const std::string &menu, const char *role )
{
  if( menu.empty() || !role || !using_native_menu() )
    return;
  
  nlohmann::json msg;
  msg["id"] = menu;
  msg["role"] = role;
  send( "MenuSetMenuRole", msg );
}//void setElectronMenuRole(...)


std::string addElectronSubMenu( const std::string &parentMenu, PopupDivMenu *menu,
                                const char *name, const std::string &parentItem )
{
  if( parentMenu.empty() || !menu || !using_native_menu() )
    return string();
  
  nlohmann::json msg;
  msg["parent"] = parentMenu;
  msg["id"] = menu->id();
  msg["item"] = parentItem;
  msg["label"] = (name ? name : "");
  send( "MenuAddSubMenu", msg );
  
  return menu->id();
}//std::string addElectronSubMenu(...)


std::string insertElectronMenuItem( const std::string &menu, PopupDivMenuItem *item, int index )
{
  if( menu.empty() || !item || !using_native_menu() )
    return string();
  
  InterSpecApp * const app = dynamic_cast<InterSpecApp *>( Wt::WApplication::instance() );
  if( !app )
    return string();
  
  const string id = item->id();
  
  nlohmann::json msg;
  msg["menu"] = menu;
  msg["id"] = id;
  msg["index"] = index;
  msg["label"] = native_label( item->text() );
  msg["enabled"] = item->isEnabled();
  msg["visible"] = !item->isHidden();
  if( !item->nativeAccelerator().empty() )
    msg["accelerator"] = item->nativeAccelerator();
  send( "MenuInsertItem", msg );
  
  app->registerElectronMenuItem( id, item );
  
  return id;
}//std::string insertElectronMenuItem(...)


void addElectronSeparatorAt( int index, const std::string &menu, const std::string &item )
{
  if( menu.empty() || item.empty() || !using_native_menu() )
    return;
  
  nlohmann::json msg;
  msg["menu"] = menu;
  msg["id"] = item;
  msg["index"] = index;
  send( "MenuInsertSeparator", msg );
}//void addElectronSeparatorAt(...)


void makeElectronMenuItemCheckable( const std::string &item, Wt::WCheckBox *cb )
{
  if( item.empty() || !cb || !using_native_menu() )
    return;
  
  nlohmann::json msg;
  msg["id"] = item;
  msg["label"] = native_label( cb->text() );
  msg["checked"] = cb->isChecked();
  msg["enabled"] = cb->isEnabled();
  send( "MenuSetCheckable", msg );
}//void makeElectronMenuItemCheckable(...)


void removeElectronMenuItem( const std::string &item )
{
  if( item.empty() || !using_native_menu() )
    return;
  
  InterSpecApp * const app = dynamic_cast<InterSpecApp *>( Wt::WApplication::instance() );
  if( app )
    app->unregisterElectronMenuItem( item );
  
  nlohmann::json msg;
  msg["id"] = item;
  send( "MenuRemoveItem", msg );
}//void removeElectronMenuItem(...)


void setElectronMenuItemHidden( const std::string &item, bool hidden )
{
  set_item_state( item, "visible", !hidden );
}

void setElectronMenuItemEnabled( const std::string &item, bool enabled )
{
  set_item_state( item, "enabled", enabled );
}

void setElectronMenuItemChecked( const std::string &item, bool checked )
{
  set_item_state( item, "checked", checked );
}

void setElectronMenuItemLabel( const std::string &item, const std::string &label )
{
  set_item_state( item, "label", label );
}

void setElectronMenuItemAccelerator( const std::string &item, const std::string &accelerator )
{
  set_item_state( item, "accelerator", accelerator );
}

void addElectronMenuItemToolTip( const std::string &item, const char *tooltip )
{
  set_item_state( item, "toolTip", string(tooltip ? tooltip : "") );
}


bool handleElectronMenuActivation( InterSpecApp *app, const std::string &payload )
{
  if( !app || payload.empty() )
    return false;
  
  // Payload is "<itemid>", or "<itemid>|0" / "<itemid>|1" for a checkbox item.
  const size_t sep_pos = payload.find( '|' );
  const string id = payload.substr( 0, sep_pos );
  
  PopupDivMenuItem * const item = app->electronMenuItem( id );
  if( !item )
  {
    cerr << "handleElectronMenuActivation: no menu item registered for '" << id << "'" << endl;
    return false;
  }
  
  // Mirrors doemit()/doemitcheck() in target/macos/NativeMenu.mm - the check box has to be brought
  //  in line with what the native menu is already showing before anything looks at it.
  if( sep_pos != string::npos )
  {
    const bool checked = (payload.substr(sep_pos + 1) == "1");
    Wt::WCheckBox * const cb = item->checkBox();
    if( cb )
    {
      cb->setChecked( checked );
      if( checked )
        cb->checked().emit();
      else
        cb->unChecked().emit();
    }//if( cb )
  }//if( a checkbox item )
  
  item->triggered().emit( static_cast<Wt::WMenuItem *>(item) );
  
  return true;
}//bool handleElectronMenuActivation(...)

#endif //USING_ELECTRON_NATIVE_MENU
