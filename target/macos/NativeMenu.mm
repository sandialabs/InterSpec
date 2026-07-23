#include <set>
#include <atomic>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>

#include <AppKit/NSCell.h>
#include <AppKit/NSImage.h>
#include <AppKit/NSMenu.h>
#include <AppKit/NSMenuItem.h>
#include <objc/runtime.h>


//We gotta fix some wierd errors...
#ifdef check
#undef check
#endif
#ifdef require
#undef require
#endif

#include <Wt/WServer.h>
#include <Wt/WCheckBox.h>
#include <Wt/WApplication.h>

#include "InterSpec/PopupDiv.h"
#include "target/macos/NativeMenu.h"
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


namespace
{
  using CallbackValidity = std::shared_ptr<std::atomic<bool>>;
  using ClickedCallback = std::function<void()>;
  using CheckedCallback = std::function<void(bool)>;
  using ClickedCallbackHandle = std::shared_ptr<const ClickedCallback>;
  using CheckedCallbackHandle = std::shared_ptr<const CheckedCallback>;

  void doemit_if_valid( const CallbackValidity &valid, PopupDivMenuItem *item )
  {
    if( valid->load( std::memory_order_acquire ) )
      doemit( item );
  }

  void doemitcheck_if_valid( const CallbackValidity &valid, Wt::WCheckBox *cb,
                             PopupDivMenuItem *item, const bool checked )
  {
    if( valid->load( std::memory_order_acquire ) )
      doemitcheck( cb, item, checked );
  }

  char TargetAssociationKey;
}


@interface Target :NSObject {
  std::string m_appid;
  // bindSafe() owns Wt::Core::observing_ptr objects, whose copy/destruction mutates the observed
  //  WObject. Keep the actual bound functions behind immutable shared handles so AppKit and
  //  WServer::post only copy shared_ptrs, never observing_ptrs.
  CallbackValidity m_callbackValid;
  ClickedCallbackHandle m_clickedCallback;
  CheckedCallbackHandle m_checkedCallback;
  std::mutex m_callbackMutex;
  // validateMenuItem is called by AppKit and must never dereference a Wt object.
  std::atomic<bool> m_enabled;
  // Set and used only on the AppKit main thread.
  NSMenuItem *m_nsitem;
}
- (id) initWithItem: (PopupDivMenuItem*)item;
- (id) initWithCb: (Wt::WCheckBox *)cb wtItem: (PopupDivMenuItem *)item;
- (void) setNSItem: (NSMenuItem *)item;
- (void) setCachedEnabled: (bool)enabled;
- (void) invalidate;
- (void) invalidateAndClearCallbacks;
- (void) clicked;
- (void) toggleChecked;
@end


@implementation Target  // <NSMenuValidation>
- (id) init {
  return [super init];
}

- (id) initWithItem: (PopupDivMenuItem*)wtItem {
  self = [super init];
  if( !self )
    return nil;

  // This initializer runs on the Wt session thread. Build every Wt lifetime guard here.
  const bool valid = (wtItem != nullptr);
  m_callbackValid = std::make_shared<std::atomic<bool>>( valid );
  if( valid )
  {
    const CallbackValidity validity = m_callbackValid;
    ClickedCallback callback = [validity, wtItem](){
      doemit_if_valid( validity, wtItem );
    };
    ClickedCallback safeCallback = wtItem->bindSafe( callback );
    m_clickedCallback = std::make_shared<const ClickedCallback>( std::move(safeCallback) );
  }

  m_enabled.store( valid && wtItem->isEnabled(), std::memory_order_release );
  m_nsitem = nil;
  m_appid = wApp->sessionId();
  return self;
}

- (id) initWithCb: (Wt::WCheckBox*)cb wtItem: (PopupDivMenuItem *)wtItem {
  self = [super init];
  if( !self )
    return nil;

  const bool valid = (cb != nullptr) && (wtItem != nullptr);
  m_callbackValid = std::make_shared<std::atomic<bool>>( valid );
  if( valid )
  {
    const CallbackValidity validity = m_callbackValid;
    CheckedCallback callback = [validity, cb, wtItem]( const bool checked ){
      doemitcheck_if_valid( validity, cb, wtItem, checked );
    };
    CheckedCallback itemSafeCallback = wtItem->bindSafe( callback );
    CheckedCallback safeCallback = cb->bindSafe( itemSafeCallback );
    m_checkedCallback = std::make_shared<const CheckedCallback>( std::move(safeCallback) );
  }

  m_enabled.store( valid && wtItem->isEnabled(), std::memory_order_release );
  m_nsitem = nil;
  m_appid = wApp->sessionId();
  return self;
}

- (void) setNSItem: (NSMenuItem *)item {
  m_nsitem = item;
}

- (void) setCachedEnabled: (bool)enabled {
  m_enabled.store( enabled, std::memory_order_release );
}

- (void) invalidate {
  // Atomics-only and idempotent: safe from either the Wt session thread or AppKit main thread.
  if( m_callbackValid )
    m_callbackValid->store( false, std::memory_order_release );
  m_enabled.store( false, std::memory_order_release );
}

- (void) invalidateAndClearCallbacks {
  // Called synchronously on the Wt session thread before widget destruction or native replacement.
  // Mark invalid first so any callback already queued by AppKit becomes a no-op after it acquires
  // the session lock. The mutex ensures AppKit cannot still be copying a callback handle when the
  // Target's session-owned handles are reset.
  [self invalidate];
  std::lock_guard<std::mutex> lock( m_callbackMutex );
  m_clickedCallback.reset();
  m_checkedCallback.reset();
}

- (BOOL) validateMenuItem: (NSMenuItem*)menuItem {
  return m_enabled.load( std::memory_order_acquire ) ? YES : NO;
}

- (void) clicked {
  std::lock_guard<std::mutex> lock( m_callbackMutex );
  if( !m_callbackValid
      || !m_callbackValid->load( std::memory_order_acquire )
      || !m_clickedCallback )
    return;

  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
  {
    const ClickedCallbackHandle callback = m_clickedCallback;
    server->post( m_appid, [callback](){ (*callback)(); } );
  }
}


- (void) toggleChecked {
  std::lock_guard<std::mutex> lock( m_callbackMutex );
  if( !m_callbackValid
      || !m_callbackValid->load( std::memory_order_acquire )
      || !m_checkedCallback )
    return;

  const bool checked = ([m_nsitem state] == NSOffState);
  [m_nsitem setState:(checked ? NSOnState : NSOffState)];

  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
  {
    const CheckedCallbackHandle callback = m_checkedCallback;
    server->post( m_appid, [callback, checked](){ (*callback)( checked ); } );
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
  const std::string itemtext = text ? text : "";
  NSMenu *newMenu = nil;

  do_in_main_sync( [&](){
    NSString *nsname = [NSString stringWithUTF8String:itemtext.c_str()];
    NSMenuItem *newItem = [[NSMenuItem alloc] initWithTitle:nsname
                                                    action:NULL
                                             keyEquivalent:@""];
    newMenu = [[NSMenu alloc] initWithTitle:nsname];

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


void *insertOsxMenuItem( void *voidmenu, PopupDivMenuItem *item, int position,
                         void **targetOut )
{
  if( !targetOut )
    return nullptr;
  *targetOut = nullptr;

  NSMenu *menu = (NSMenu *)voidmenu;
  const std::string itemtext = item->text().toUTF8();
  const std::string iconpath = item->icon();
  Target* target = [[Target alloc] initWithItem:item];
  if( !target )
    return nullptr;

  NSMenuItem *itemnow = nil;
  do_in_main_sync( [&](){
    NSString *name = [NSString stringWithUTF8String:itemtext.c_str()];
    itemnow = [[NSMenuItem alloc]
               initWithTitle:name
                      action:@selector(clicked)
               keyEquivalent:@""];
    [target setNSItem:itemnow];
    [itemnow setTarget:target];
    [itemnow setEnabled:YES];
    objc_setAssociatedObject( itemnow, &TargetAssociationKey, target,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC );

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
    
    if( !iconpath.empty() )
    {
      NSString *nsiconpath = [NSString stringWithUTF8String:iconpath.c_str()];
      NSImage *image = [[NSImage alloc] initByReferencingFile:nsiconpath];
      // TODO: should resize the image to like 16x16px, if it isnt already
      [itemnow setImage:image];
      [image release];
    }

    // Keep the allocation retain as the Wt menu item's ownership. NSMenu may remove the item during
    // session recovery before PopupDivMenuItem is destroyed; this retain keeps the raw item pointer
    // valid until removeOsxMenuItem relinquishes ownership on the main queue.
  } );

  if( !itemnow )
  {
    [target invalidateAndClearCallbacks];
    [target release];
    return nullptr;
  }

  *targetOut = target;
  item->setData( (void *)itemnow );
  
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
    [item release];
  } );
}//void removeOsxSeparator( ( void *voidmenu, void *voiditem )


void *addOsxCheckableMenuItem( void *voidmenu, Wt::WCheckBox *cb,
                               PopupDivMenuItem *wtItem, void **targetOut )
{
  if( !targetOut )
    return nullptr;
  *targetOut = nullptr;
  if( !voidmenu || !cb || !wtItem )
    return nullptr;

  NSMenu *menu = (NSMenu *)voidmenu;
  const std::string itemtext = cb->text().toUTF8();
  const bool checked = cb->isChecked();
  Target* target = [[Target alloc] initWithCb:cb wtItem:wtItem];
  if( !target )
    return nullptr;

  NSMenuItem *itemnow = nil;
  do_in_main_sync( [&](){
    NSString *name = [NSString stringWithUTF8String:itemtext.c_str()];
    itemnow = [[NSMenuItem alloc]
               initWithTitle:name
                      action:@selector(toggleChecked)
               keyEquivalent:@""];
    [target setNSItem:itemnow];
    [itemnow setTarget:target];
    [itemnow setEnabled:YES];
    [itemnow setState:(checked ? NSOnState : NSOffState)];
    objc_setAssociatedObject( itemnow, &TargetAssociationKey, target,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC );
    [menu addItem:itemnow];

    // Keep the allocation retain as the Wt menu item's ownership; see insertOsxMenuItem().
  } );

  if( !itemnow )
  {
    [target invalidateAndClearCallbacks];
    [target release];
    return nullptr;
  }

  *targetOut = target;
  
  return itemnow;
}//void *addOsxCheckableMenuItem( void *menu, Wt::WCheckBox *cb );


void *addOsxSeparatorAt( int index, void *voidmenu )
{
  if( !voidmenu  )
    return 0;
    
  NSMenu *menu = (NSMenu *)voidmenu;
  NSMenuItem *item = nil;

  do_in_main_sync( [&](){
    // Keep an explicit retain until removeOsxSeparator's asynchronous block completes.
    item = [[NSMenuItem separatorItem] retain];
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
  NSMenuItem *item = nil;

  do_in_main_sync( [&](){
    item = [[NSMenuItem separatorItem] retain];
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
    Target *target = (Target *)objc_getAssociatedObject( i, &TargetAssociationKey );
    [i setTarget:nil];
    if( target )
      [target invalidate];

    const NSInteger index = [m indexOfItem: i];
    if( index >= 0 )
      [m removeItem:i];

    // Relinquish the NSMenuItem's ownership of Target. The session-side bridge is released
    // synchronously before this block is scheduled.
    objc_setAssociatedObject( i, &TargetAssociationKey, nil,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC );

    // Relinquish the explicit allocation retain. If this shutdown-time block never executes,
    // leaking until process exit is safer than synchronously crossing to AppKit and deadlocking.
    [i release];
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


void invalidateOsxMenuItemTarget( void *voidtarget )
{
  if( !voidtarget )
    return;

  // This bridge is separate from NSMenuItem, so the Wt session thread never messages an AppKit
  // object to obtain Target. Clear bindSafe callback ownership on the session thread, then
  // relinquish the Wt item's retain; the NSMenuItem association keeps Target alive until removal.
  Target *target = (Target *)voidtarget;
  [target invalidateAndClearCallbacks];
  [target release];
}//void invalidateOsxMenuItemTarget( void *target )


void setOsxMenuItemTargetEnabled( void *voidtarget, bool enabled )
{
  if( !voidtarget )
    return;

  [(Target *)voidtarget setCachedEnabled:enabled];
}//void setOsxMenuItemTargetEnabled( void *target, bool enabled )



void addOsxMenuItemToolTip( void *item, const char *tooltip )
{
  if( !item || !tooltip )
    return;
  
  NSMenuItem *i = (NSMenuItem *)item;
  const std::string tooltipText = tooltip;
  
  do_in_main_sync( [=](){
    NSString *tip = [NSString stringWithUTF8String:tooltipText.c_str()];
    [i setToolTip:tip];
  } );
}//void addOsxMenuItemToolTip( void *item, const char *tooltip )
