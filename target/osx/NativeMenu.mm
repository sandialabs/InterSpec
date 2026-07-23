#include <set>
#include <atomic>
#include <iostream>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

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


namespace
{
  typedef boost::shared_ptr<std::atomic<bool> > CallbackValidity;

  void doemit_if_valid( const CallbackValidity &valid, PopupDivMenuItem *item )
  {
    if( valid->load() )
      doemit( item );
  }

  void doemitcheck_if_valid( const CallbackValidity &valid, Wt::WCheckBox *cb,
                             PopupDivMenuItem *item, const bool checked )
  {
    if( valid->load() )
      doemitcheck( cb, item, checked );
  }
}


@interface Target :NSObject {
  std::string m_appid;
  // These immutable callbacks are created on the Wt session thread using WApplication::bind() and
  // boost::bind(). Wt therefore disconnects their slots automatically if a bound WObject is
  // destroyed. The shared atomic additionally suppresses callbacks when a native target is
  // invalidated without destroying its Wt item (for example, conversion to a checkable item).
  CallbackValidity m_callbackValid;
  boost::function<void()> m_clickedCallback;
  boost::function<void()> m_checkedCallback;
  boost::function<void()> m_uncheckedCallback;
  // Cached enabled state, pushed from the session thread (initWith.../setCachedEnabled) and read by
  // validateMenuItem on the AppKit main thread - so validateMenuItem never dereferences the Wt
  // widget (which would be a cross-thread data race / use-after-free).
  std::atomic<bool> m_enabled;
  NSMenuItem *m_nsitem;
}
- (id) initWithItem: (PopupDivMenuItem*)item;
- (id) initWithCb: (Wt::WCheckBox *)cb wtItem: (PopupDivMenuItem *)item;
- (void) setNSItem: (NSMenuItem *)item;
- (void) setCachedEnabled: (bool)enabled;
- (void) invalidate;
- (void) clicked;
- (void) toggleChecked;
@end


namespace
{
  char TargetAssociationKey;
}


@implementation Target  // <NSMenuValidation>
- (id) init {
  return [super init];
}

- (id) initWithItem: (PopupDivMenuItem*)wtItem {
  self = [super init];
  if( !self )
    return nil;

  // Runs on the Wt session thread during menu construction, so reading the widget is safe here.
  m_callbackValid.reset( new std::atomic<bool>( true ) );
  m_clickedCallback = wApp->bind( boost::bind( &doemit_if_valid, m_callbackValid, wtItem ) );
  m_enabled.store( wtItem && wtItem->isEnabled() );
  m_nsitem = 0;
  m_appid = wApp->sessionId();
  return self;
}

- (id) initWithCb: (Wt::WCheckBox*)cb wtItem: (PopupDivMenuItem *)wtItem {
  self = [super init];
  if( !self )
    return nil;

  m_callbackValid.reset( new std::atomic<bool>( true ) );
  m_checkedCallback = wApp->bind( boost::bind( &doemitcheck_if_valid, m_callbackValid,
                                              cb, wtItem, true ) );
  m_uncheckedCallback = wApp->bind( boost::bind( &doemitcheck_if_valid, m_callbackValid,
                                                cb, wtItem, false ) );
  m_enabled.store( wtItem && wtItem->isEnabled() );
  m_nsitem = 0;
  m_appid = wApp->sessionId();
  return self;
}

- (void) setNSItem: (NSMenuItem *)item {
  m_nsitem = item;
}

- (void) setCachedEnabled: (bool)enabled {
  m_enabled.store( enabled );
}

- (void) invalidate {
  // The shared flag is also captured by every posted callback, so a callback that was queued just
  // before invalidation will re-check it after WServer::post acquires the Wt session lock.
  if( m_callbackValid )
    m_callbackValid->store( false );
  m_enabled.store( false );
}

- (BOOL) validateMenuItem: (NSMenuItem*)menuItem {
  // Called by AppKit on the main thread whenever the menu is about to display. Read the cached
  // enabled flag instead of dereferencing the Wt widget (which lives on, and is mutated/destroyed
  // by, the session thread) - this avoids a cross-thread data race and a use-after-free.
  return m_enabled.load() ? YES : NO;
}

- (void) clicked {
  if( !m_callbackValid || !m_callbackValid->load() )
    return;

  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
    server->post( m_appid, m_clickedCallback );
}


- (void) toggleChecked {
  if( !m_callbackValid || !m_callbackValid->load() )
    return;

  bool checked = false;
  if( [m_nsitem state] == NSOffState )
  {
    [m_nsitem setState:NSOnState];
    checked = true;
  }else
  {
    [m_nsitem setState:NSOffState];
  }

  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
    server->post( m_appid, checked ? m_checkedCallback : m_uncheckedCallback );
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


void *insertOsxMenuItem( void *voidmenu, PopupDivMenuItem *item, int position,
                         void **targetOut )
{
  if( targetOut )
    *targetOut = nullptr;

  NSMenu *menu = (NSMenu *)voidmenu;
  const std::string itemtext = item->text().toUTF8();
  const std::string iconpath = item->icon();
  Target* target = [[Target alloc] initWithItem:item];
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
    // session recovery before PopupDivMenuItem is destroyed; this retain keeps its raw pointer valid
    // until removeOsxMenuItem relinquishes that ownership on the main queue.
  } );

  if( targetOut )
    *targetOut = target;
  else
    [target release];
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
  } );
}//void removeOsxSeparator( ( void *voidmenu, void *voiditem )


void *addOsxCheckableMenuItem( void *voidmenu, Wt::WCheckBox *cb,
                               PopupDivMenuItem *wtItem, void **targetOut )
{
  if( targetOut )
    *targetOut = nullptr;
  if( !voidmenu || !cb )
    return 0;

  NSMenu *menu = (NSMenu *)voidmenu;
  const std::string itemtext = cb->text().toUTF8();
  const bool checked = cb->isChecked();
  Target* target = [[Target alloc] initWithCb:cb wtItem:wtItem];
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
    [itemnow setState:checked];
    objc_setAssociatedObject( itemnow, &TargetAssociationKey, target,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC );
    [menu addItem:itemnow];

    // Keep the allocation retain as the Wt menu item's ownership; see insertOsxMenuItem().
  } );

  if( targetOut )
    *targetOut = target;
  else
    [target release];
  
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
    Target *target = (Target *)[i target];
    [i setTarget:nil];
    if( target )
      [target invalidate];

    const NSInteger index = [m indexOfItem: i];
    if( index >= 0 )
      [m removeItem:i];

    // Releasing the association drops the menu item's ownership of Target. Posted callbacks own
    // their Wt-bound function and shared validity token independently of Target.
    objc_setAssociatedObject( i, &TargetAssociationKey, nil,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC );

    // Relinquish the allocation retain held for the Wt menu item. If this shutdown-time block never
    // executes, leaking the native item until process exit is safer than synchronously crossing to
    // AppKit and risking the teardown deadlock this asynchronous path avoids.
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

  // This opaque bridge is separate from NSMenuItem. The method only writes std::atomics and is safe
  // to call synchronously from the Wt session thread. Relinquish the Wt item's ownership only after
  // invalidation; the NSMenuItem association keeps the Target alive until asynchronous removal.
  Target *target = (Target *)voidtarget;
  [target invalidate];
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
  NSString *tip = [NSString stringWithFormat:@"%s", tooltip];
  
  do_in_main_sync( [=](){
    [i setToolTip:tip];
  } );
}//void addOsxMenuItemToolTip( void *item, const char *tooltip )
