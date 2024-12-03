#include <set>
#include <iostream>

#include <AppKit/NSCell.h>
#include <AppKit/NSImage.h>
#include <AppKit/NSMenu.h>
#include <AppKit/NSMenuItem.h>


//We gotta fix some wierd errors...
#ifdef check
#undef check
#endif
#ifdef require
#undef require
#endif

#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WApplication>

#include "InterSpec/PopupDiv.h"
#include "target/osx/NativeMenu.h"
#include "InterSpec/InterSpecApp.h"


void do_in_main_sync( std::function<void()> work )
{
  if( [NSThread isMainThread] )
  {
    work();
  }else
  {
    dispatch_sync(dispatch_get_main_queue(), ^{
      work();
    } );
  }
}//void do_in_main_sync(...)


void doemit( PopupDivMenuItem *item )
{
  item->triggered().emit( (Wt::WMenuItem *)item );
  wApp->triggerUpdate();
}//void doemit( PopupDivMenuItem *item )


void doemitcheck( Wt::WCheckBox *cb, PopupDivMenuItem *item, const bool checked )
{
  cb->setChecked( checked );
  if( checked )
    cb->checked().emit();
  else
    cb->unChecked().emit();
  
  if( item )
    item->triggered().emit( (Wt::WMenuItem *)item );
  
  wApp->triggerUpdate();
}//void doemitcheck( Wt::WCheckBox *cb, const bool checked )


@interface Target :NSObject {
  std::string m_appid;
  Wt::WCheckBox *m_cb;
  PopupDivMenuItem *m_item;
  NSMenuItem *m_nsitem;
}
- (id) initWithItem: (PopupDivMenuItem*)item;
- (id) initWithCb: (Wt::WCheckBox *)cb;
- (void) setNSItem: (NSMenuItem *)item;
- (void) setWtItem: (PopupDivMenuItem *)item;
- (void) clicked;
- (void) toggleChecked;
@end

@implementation Target  // <NSMenuValidation>
- (id) init {
  return [super init];
}

- (id) initWithItem: (PopupDivMenuItem*)wtItem {
  m_cb = 0;
  m_item = wtItem;
  m_nsitem = 0;
  m_appid = Wt::WApplication::instance()->sessionId();
  return self;
}

- (id) initWithCb: (Wt::WCheckBox*)cb {
  m_cb = cb;
  m_item = 0;
  m_nsitem = 0;
  m_appid = Wt::WApplication::instance()->sessionId();
  return self;
}

- (void) setNSItem: (NSMenuItem *)item {
  m_nsitem = item;
}

- (void) setWtItem: (PopupDivMenuItem *)item {
  m_item = item;
}

- (BOOL) validateMenuItem: (NSMenuItem*)menuItem {
  // This is called when the menu is shown.
  if( m_item )
    return m_item->isEnabled();
  if( m_cb )
    return m_cb->isEnabled();
  return NO;
}

- (void) clicked {
  if( m_item )
    Wt::WServer::instance()->post( m_appid, boost::bind( &doemit, m_item ) );
}


- (void) toggleChecked {
  if( !m_cb )
    return;
  if( [m_nsitem state] == NSOffState )
  {
    [m_nsitem setState:NSOnState];
    Wt::WServer::instance()->post( m_appid, boost::bind( &doemitcheck, m_cb, m_item, true ) );
  }else
  {
    [m_nsitem setState:NSOffState];
    Wt::WServer::instance()->post( m_appid, boost::bind( &doemitcheck, m_cb, m_item, false ) );
  }
}


@end


void *addOsxMenu( PopupDivMenu *menu, const char *name  )
{
  //Make it so the primary (first) InterSpecApp session can load
  //  menus to the OSX toolbar.  Note that menus are added during InterSpec
  //  initialization, which happens before the WApplication pointer is added
  //  to the running instances, so we require there to be no running instances
  //  in order to add the menus to the OSX menu bar.
  
  if( !InterSpecApp::isPrimaryWindowInstance() )
    return 0;
  
  NSMenu *newMenu = nil;
  
  //Cleans out the old menu if it exists before replacing it
  auto doWork = [name,&newMenu](){
    NSString *nsname = [NSString stringWithFormat:@"%s", name];
    NSMenu* rootMenu = [NSApp mainMenu];
    
    //check if menu already exists, if so, remove
    for( NSInteger i = 0; i < [[NSApp mainMenu] numberOfItems]; ++i )
    {
      NSMenu* menu = [[[NSApp mainMenu] itemAtIndex: i] submenu];
      //Checks the name, but also make sure it's not the main "InterSpec" menu
      if( [menu.title isEqualToString:nsname] )
      {
        if( [[menu title] isEqualToString:@"Edit"] )
        {
          newMenu = menu;
          return;
        }else if( ![[menu title] isEqualToString:@"InterSpec"] )
        {
          [rootMenu removeItemAtIndex:[rootMenu indexOfItemWithSubmenu:menu]];
        }
      }
    }
    
    
    if( std::strcmp(name,"InterSpec")==0 )
    {
      NSMenu* ret = [[[NSApp mainMenu] itemAtIndex: 0] submenu];
      //Go through all the menuitems
      for( NSMenuItem *item in [ret itemArray] )
      {
        NSString *menuString = item.title;
        if( ![menuString isEqualToString:@"Quit InterSpec"] )
        {
          //Remove everything if it is not Quit menuitem
          [ret removeItem:item];
        }
      }
      
      newMenu = ret;
      return;
    }//InterSpec menu
    
    NSMenuItem *newItem = [[NSMenuItem alloc] initWithTitle:@"" action:NULL keyEquivalent:@""];
    newMenu = [[NSMenu alloc] initWithTitle:nsname];
    
    [newItem setSubmenu:newMenu];
    [[NSApp mainMenu] addItem:newItem];
  };//doWork
  
  do_in_main_sync( doWork );
  
  return newMenu;
}//void *addOsxMenu( PopupDivMenu *menu, const char *name  )


void *addOsxSubMenu( void *parent, PopupDivMenu *item, const char *text )
{
  NSMenu *parentmenu = (NSMenu *)parent;
  NSString *nsname = [NSString stringWithFormat:@"%s", text];
  
  NSMenuItem *newItem = [[NSMenuItem alloc] initWithTitle:nsname action:NULL keyEquivalent:@""];
  NSMenu *newMenu = [[NSMenu alloc] initWithTitle:nsname];

  do_in_main_sync( [=](){
    NSInteger ind = [parentmenu indexOfItemWithTitle:@"Quit InterSpec"];
    if (ind!=-1)
    {
        //Adds the submenu before Quit menuitem
        [parentmenu insertItem:newItem atIndex:ind];
    }else
    {
        //Regularly add if there is no Quit menuitem in this menu
        //Should do this next part insdie the application thread?
        [parentmenu addItem:newItem];
    }

    [parentmenu setSubmenu:newMenu forItem:newItem];
  } );
  
  return newMenu;
}//void *addOsxSubMenu( void *parent, PopupDivMenu *item )


void *insertOsxMenuItem( void *voidmenu, PopupDivMenuItem *item, int position )
{
  NSMenu *menu;
    
  if( item->parentMenu()->parentItem() )
  {
    if( item->parentMenu()->parentItem()->text().toUTF8().compare("InterSpec")==0 )
    {
      //Special case: if the menu is 'InterSpec', instead add it to the default OSX 'InterSpec menu'
      menu = [[[NSApp mainMenu] itemAtIndex: 1] submenu];
    }
    
    menu = (NSMenu *)voidmenu;
  }else
  {
    //Regularly add if there is no Quit menuitem in this menu
    menu = (NSMenu *)voidmenu;
  }
  
  NSString *name = [[NSString alloc]initWithUTF8String: item->text().toUTF8().c_str()];
  
  NSMenuItem *itemnow = [[NSMenuItem alloc]
                         initWithTitle:name
                                action:@selector(clicked)
                         keyEquivalent:@""];

  const std::string iconpath = item->icon();
  NSString *nsiconpath = nil;
  if( iconpath.size() )
    nsiconpath = [NSString stringWithFormat:@"%s", iconpath.c_str()];
  Target* target = [[Target alloc] initWithItem:item];
  [target setNSItem:itemnow];
  [itemnow setTarget:target];
  [itemnow setEnabled:YES];
  
  item->setData( (void *)itemnow );

  do_in_main_sync( [=](){
    NSInteger ind = [menu indexOfItemWithTitle:@"Quit InterSpec"];
    if( position >= 0 )
    {
      [menu insertItem:itemnow atIndex:(position)];
    }else if( ind != -1 )
    {
      //Add menuitem before Quit InterSpec menuitem
      [menu insertItem:itemnow atIndex:(ind)];
    }else
    {
      [menu addItem:itemnow];
    }
    
    if( nsiconpath )
    {
      NSImage *image = [[NSImage alloc] initByReferencingFile:nsiconpath];
      // TODO: should resize the image to like 16x16px, if it isnt already
      [itemnow setImage:image];
    }
  } );
  
  return itemnow;
}//void *addOsxMenuItem( void *voidmenu, PopupDivMenuItem *item )


void removeOsxSeparator( void *voidmenu, void *voiditem )
{
  NSMenu *menu = (NSMenu *)voidmenu;
  NSMenuItem *item = (NSMenuItem *)voiditem;
  
  if( !menu || !item )
    return;

  // We will do this async, because if we are quitting the app, we actually dont care if it doesnt
  //  get removed, and in this case the async call to main thread wont ever happen (I dont think)
  dispatch_async(dispatch_get_main_queue(), ^{
    const NSInteger index = [menu indexOfItem: item];
    if( index >= 0 )
      [menu removeItem:item];
  } );
}//void removeOsxSeparator( ( void *voidmenu, void *voiditem )


void *addOsxCheckableMenuItem( void *voidmenu, Wt::WCheckBox *cb,
                               PopupDivMenuItem *wtItem )
{
  if( !voidmenu || !cb )
    return 0;

  NSMenu *menu = (NSMenu *)voidmenu;
  NSString *name = [NSString stringWithFormat:@"%s", cb->text().toUTF8().c_str()];
  Target* target = [[Target alloc] initWithCb:cb];
  
  NSMenuItem *itemnow = [[NSMenuItem alloc]
                         initWithTitle:name
                         action:@selector(toggleChecked)
                         keyEquivalent:@""];
  
  //Should do this next part insdie the application thread?
  [target setNSItem:itemnow];
  [target setWtItem:wtItem];
  
  [itemnow setTarget:target];
  [itemnow setEnabled:YES];
  [itemnow setState:cb->isChecked()];
  
  do_in_main_sync( [=](){
    [menu addItem:itemnow];
  } );
  
  return itemnow;
}//void *addOsxCheckableMenuItem( void *menu, Wt::WCheckBox *cb );


void *addOsxSeparatorAt( int index, void *voidmenu )
{
  if( !voidmenu  )
    return 0;
    
  NSMenu *menu = (NSMenu *)voidmenu;
  NSMenuItem *item = [NSMenuItem separatorItem];

  do_in_main_sync( [=](){
    if( index >= 0 )
      [menu insertItem:item atIndex:(index)];
    else
      [menu addItem:item];
  } );

  return item;
}//void *addOsxSeparatorAt( int index, void *voidmenu )


void *addOsxSeparator(void *voidmenu)
{
  if( !voidmenu  )
    return 0;
    
  NSMenu *menu = (NSMenu *)voidmenu;
  NSMenuItem *item = [NSMenuItem separatorItem];

  do_in_main_sync( [=](){
    NSInteger ind = [menu indexOfItemWithTitle:@"Quit InterSpec"];  
    if( ind != -1 )
    {
      //Add seperator before Quit InterSpec menuitem
      [menu insertItem:item atIndex:(ind)];
    }else
    {
      [menu addItem:item]; // Add seperator the normal way - a thin grey line
    }
  } );
  
  return item;
} //void *addOsxSeparator(void *voidmenu)


//void removeOsxMenu( void *menu )
//{
//  NSMenu *m = (NSMenu *)menu;
//  do_in_main_sync( [=](){
//    const NSInteger index = [[NSApp mainMenu] indexOfItem: m];
//    if( index >= 0 )
//      [[NSApp mainMenu] removeItem:m];
//  } );
//}//void removeOsxMenu( void *menu )


void removeOsxMenuItem( void *item, void *menu )
{
  if( !item || !menu )
    return;
    
  NSMenu *m = (NSMenu *)menu;
  NSMenuItem *i = (NSMenuItem *)item;
  
  // We will do this async, because if we are quitting the app, we actually dont care if it doesnt
  //  get removed, and in this case the async call to main thread wont ever happen (I dont think)
  dispatch_async(dispatch_get_main_queue(), ^{
    const NSInteger index = [m indexOfItem: i];
    if( index >= 0 )
      [m removeItem:i];
  } );
}//void removeOsxMenuItem( void *item )


void setOsxMenuItemHidden( void *item, bool hidden )
{
  if( !item )
    return;
  
  NSMenuItem *i = (NSMenuItem *)item;
  
  do_in_main_sync( [=](){
    [i setHidden:hidden];
  } );
}//void setOsxMenuItemHidden( void *item, bool hidden )



void addOsxMenuItemToolTip( void *item, const char *tooltip )
{
  if( !item || !tooltip )
    return;
  
  NSMenuItem *i = (NSMenuItem *)item;
  NSString *tip = [NSString stringWithFormat:@"%s", tooltip];
  
  do_in_main_sync( [=](){
    [i setToolTip:tip];
  } );
}//void addOsxMenuItemToolTip( void *item, const char *tooltip )

