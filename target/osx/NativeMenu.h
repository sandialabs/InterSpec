#ifndef NativeMenu_h
#define NativeMenu_h

//The basic idea for the functions in this header/source files is to allow
//  maping the PopupDivMenu and PopupDivMenuItem classes into the NSMenu and
//  NSMenuItem objects of the native OS X menu system.  Right now it is only at
//  a proof of concept level, but kinda seems to be working okay.
//Also, I didnt know how to forward declare an objective-c class in c++, so
//  I took the (horrible) path of using void* for all objective-c objects.
//
//Two caveates, or special conditions are known so far:
//  1) There is no way to insert a SpinBox into a menu item, let alone
//     syncornize with Wt.
//  2) If you are relying on the PopupDivMenuItem containing a URL to redirect
//     the browser (or open a new window), then this wont work.  The workaround
//     is to introduce an extra PopupDivMenuItem item in the WPopoupMenu (but
//     not map to a native menu, use WPopupMenu::addItem( WMenuItem * ) ), which
//     instead maps its anchor to a URL/WLink, and have the item for the
//     original menu item, click the second link via javascript
//     (document.getElementById('second_id').click()) whenever the original item
//     is clicked.  InterSpec::addFileMenu(...) has an implentation of
//     this, as well as some notes why this awkwardness is necassary.

class PopupDivMenu;
class PopupDivMenuItem;

namespace Wt
{
  class WCheckBox;
}

//To avoid some compilation errors, we have to be a little unsafe and pass
//  around the NSMenu and NSMenuItem pointers as (void *)

//addOsxMenu(): returns a pointer to the OS X NSMenu added
void *addOsxMenu( PopupDivMenu *menu, const char *name );

//addOsxSubMenu(): returns a pointer to the OS X NSMenu added.
//  'parent' is the pointer to the NSMenu to add the sub menu to.
void *addOsxSubMenu( void *parent, PopupDivMenu *item, const char *name );

//addOsxMenuItem(): returns NSMenuItem pointer of item added.
//  'menu' is the pointer to the NSMenu to add to
void *addOsxMenuItem( void *menu, PopupDivMenuItem *item );

//addOsxCheckableMenuItem(): returns NSMenuItem pointer of item added.
//  'menu' is the pointer to the NSMenu to add to
//  The PopupDivMenuItem pointer is necassarry in order for the triggered()
//  signal to be sent.
void *addOsxCheckableMenuItem( void *menu, Wt::WCheckBox *cb,
                               PopupDivMenuItem *item );

//addOsxSeparator(): Add separator to the NSMenu 'voidmenu' passed in.
void *addOsxSeparator( void *voidmenu );

//removeOsxMenu(): /* Not Implemented */
//void removeOsxMenu( void *menu );

//removeOsxMenuItem(): removes the NSMenuItem passed in.
void removeOsxMenuItem( void *item, void *menu );

//setOsxMenuItemHidden(): sets the NSMenuItem passed in to hidden.
void setOsxMenuItemHidden( void *item, bool hidden );

//addOsxMenuItemToolTip(): adds tool tip to the NSMenuItem
void addOsxMenuItemToolTip( void *item, const char *tooltip );


#endif