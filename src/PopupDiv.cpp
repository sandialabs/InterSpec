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

#include <Wt/WText>
#include <Wt/WPoint>
#include <Wt/WLabel>
#include <Wt/WString>
#include <Wt/WAnchor>
#include <Wt/WCheckBox>
#include <Wt/WJavaScript>
#include <Wt/WPushButton>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include "InterSpec/PopupDiv.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"

#if(USE_OSX_NATIVE_MENU)
#include "target/osx/NativeMenu.h"
#endif

using namespace Wt;
using namespace std;

#include <Wt/WPopupWidget>

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


WT_DECLARE_WT_MEMBER
(OnPopupDivShow, Wt::JavaScriptFunction, "OnPopupDivShow",
 function(id,p)
 {
   try{  //this try catch inserted for for debugging sake (exception on early menu showing when using websockets locally) at some point, and probably isnt necassary
   $('.PopupDivMenu').each(function(i,v){if(v.id!=id)$(v).hide();} );
   Wt.WT.BringAboveDialogs(id);
   $('.PopupMenuParentButton').removeClass('active');
   if(p) $('#'+p).addClass('active');
   }catch(e){
     console.log( 'OnPopupDivShow caught: ' + e );
   }
});

WT_DECLARE_WT_MEMBER
(BringAboveDialogs, Wt::JavaScriptFunction, "BringAboveDialogs",
 function( id )
 {
   var z = 0;
   $('.Wt-dialog').each(function(i,v){z=Math.max(z,$(v).css('z-index'));});
   if(z>$('#'+id).css('z-index'))$('#'+id).css('z-index',z+1);
});


WT_DECLARE_WT_MEMBER
(ShowPhone, Wt::JavaScriptFunction, "ShowPhone",
 function( id )
 {
   var w = Wt.WT.getElement(id);
   if(!w) return;
   
   var t0 = (new Date()).valueOf();
   var width = $(w).width();
   w.style['top'] = '0px';
   (function slidein(){
     var t = (new Date()).valueOf();
     var pos = Math.min(0,(((t-t0)/250) - 1)*width);
     w.style['left'] = pos + 'px';
     if(pos<0)
       setTimeout(slidein,1);
  })();
   
     
   var z = 0;
   $('.Wt-dialog').each(function(i,v){z=Math.max(z,$(v).css('z-index'));});
   $('.MobileMenuButton').each(function(i,v){z=Math.max(z,$(v).css('z-index'));});
   $('.PopupDivMenu').each(function(i,v){z=Math.max(z,$(v).css('z-index'));});
   if(z>$('#'+id).css('z-index'))
   {
       //if menu is below dialogs, set above, and also set overlay right below it
       $('#'+id).css('z-index',z+2);
       if ($('.mobilePopupMenuOverlay').is(':hidden'))
       {
           $('.mobilePopupMenuOverlay').css('z-index',z+1);
       } //check if the overlay is hidden already
   } //z>$('#'+id).css('z-index')
   else
   {
       if ($('.mobilePopupMenuOverlay').is(':hidden'))
       {
           $('.mobilePopupMenuOverlay').css('z-index',$('#'+id).css('z-index')-1);
       } //$('.mobilePopupMenuOverlay').is(':hidden')
   } //already above max z index, so not a problem.
     
   //Immediately make sure none of items for this menu have the active class
//   var k = w.childNodes;
//   for(var i = 0; i<k.length;++i)
//     if(k[i].className)
//       k[i].className = k[i].className.replace(/\\bactive\\b/, "");
   $(w).children('.active').removeClass("active");
   
   //make sure none of the items for any menu have the active class, but only
   //  after all the animations have occurred
   setTimeout( function(){
     $('.PopupDivMenuPhone').each( function(){
       $(this).children('.active').removeClass("active");
     } );
   }, 500 );
 });

WT_DECLARE_WT_MEMBER(HideOverlay, Wt::JavaScriptFunction, "HideOverlay", function()
{
  function hideverything(e)
  {
    var menus = $('.PopupDivMenuPhone');
    if( !menus.is(e.target) // if the target of the click isn't the menus...
        && menus.has(e.target).length === 0) // ... nor a descendant of the menus
    {
      //slide menus out to the left and then hide
      menus.each( function( index, el ){
        var jel = $(el);
        jel.children('.active').removeClass("active");
        if( !jel.is(":visible") ) return;
        
        //Scroll back to the top of the menu
        jel.scrollTop(0);
        
        //I cant actually tell the differnce in smoothness between my version
        //  and jQueries version.
        var a = -jel.width() - 10;
        jel.animate({left: a+'px'}, {queue: false, duration: 255}, "linear", function(){jel.hide();});
        
//        var w = jel.width(), t0 = (new Date()).valueOf();
//        ( function slideout()
//          {
//            var t = (new Date()).valueOf();
//            var pos = -w*((t-t0)/250);
//            el.style['left'] = pos + 'px';
//            if( pos>-w )
//              setTimeout(slideout,1);
//            else jel.hidde();
//          })()
      });
          
      $('.mobilePopupMenuOverlay').hide();
      $(document).off("touchstart.mobileOverlay");  //unbind event listener
    } //if(!menus.is(e.target)
  }//function hideverything(e)
    
  //define mobileOverlay namespace, so can easily remove event listener later
  $(document).on("touchstart.mobileOverlay",hideverything);
});

//XXX - this is a hack!  The Wt JavaScript function fitToWindow(...) doesnt
//  make it so the sub menu will always start with the top located at >=0px,
//  so we'll fix this when this happens.  Note we have to do the actual
//  repositioning in a timeout since fitToWindow(...) et al will do their
//  work after the mouseWentOver() signal is emitted.  Note I tried
//  placing this javascript into the overiden PopupDivMenu::setHidden(...)
//  function, but it didnt work (presumably for same reason we nee a timeout
//  function here).
WT_DECLARE_WT_MEMBER
(AdjustTopPos, Wt::JavaScriptFunction, "AdjustTopPos",
 function( id )
 {
   var w = Wt.WT.getElement(id);
   if(!w)
     return;
   
   var fcn = function()
   {
     if( Wt.WT.widgetPageCoordinates(w).y < 0 )
     {
       w.style['bottom'] = 'auto';
       w.style['top'] = '2px';
     }
     
     if( Wt.WT.widgetPageCoordinates(w).x < 0 )
     {
       w.style['right'] = 'auto';
       w.style['left'] = '2px';
     }
   };
   setTimeout( fcn, 50 );
 });

#if( USING_ELECTRON_NATIVE_MENU )
WT_DECLARE_WT_MEMBER(FindElectronMenu, Wt::JavaScriptFunction, "FindElectronMenu",
function( name, selfid )
{
  let appmenu = Menu.getApplicationMenu();
  //remote.getCurrentWindow()
  
  //console.log( 'appmenu:' );
  //console.log( appmenu );
  
  if(!appmenu){
    console.log('Couldnt get app menu');
    return;
  }
  
  for( var i of appmenu.items ){
    if( i.label == name )
    {
      //console.log( 'Found menu ' + name );
      //console.log( i );
      
      $(window).data('electronMenu'+selfid, i );
      //$('#'+selfid).data('electronMenu',i);
      //console.log('Set Electron Menu Item: ' + i.label );
      return;
    }
  }
  
  //Start of trying to dynamically add menu items...
  //let newmenu = new MenuItem( {type: 'submenu', label: name, id: selfid, submenu: appmenu } );
  //appmenu.items.push(newmenu);
  //$(window).data('electronMenu'+selfid, newmenu );
  //Menu.setApplicationMenu(appmenu);
  
  console.log('Failed to find Electron Menu Item data for : ' + name );
});

WT_DECLARE_WT_MEMBER(AddMenuItemToElectronMenu, Wt::JavaScriptFunction, "AddMenuItemToElectronMenu",
function( selfid, txt, iconstr, itemid, roleType )
{
  //let m = $('#'+selfid).data('electronMenu');
  let m = $(window).data('electronMenu'+selfid);
  
  if( !m )
  {
    console.log( 'No electron menu set when trying to add: ' + txt );
    return;
  }

  var newItem = new MenuItem({label: txt,
    role: roleType,
    icon: ((iconstr && (iconstr.length>0)) ? iconstr : null),
    click(menu,event,modifiers){ Wt.emit(itemid,'electron_clicked');}
  } );
  
  m.submenu.append( newItem );
  
  $(window).data('electronItem'+itemid, newItem);
  
  if( $(window).data('HaveTriggeredMenuUpdate') ){
    let appmenu = Menu.getApplicationMenu();
    if(appmenu)
      Menu.setApplicationMenu(appmenu);
    else
      console.log('Couldnt get app menu');  //dont think getting app menu has ever failed, that I've seen
  }
});

WT_DECLARE_WT_MEMBER(AddSeperatorToElectronMenu, Wt::JavaScriptFunction, "AddSeperatorToElectronMenu",
function( selfid )
{
  let m = $(window).data('electronMenu'+selfid);
  
  let appmenu = Menu.getApplicationMenu();
  if( !m || !appmenu ) return;

  m.submenu.append( new MenuItem({type: 'separator'}) );
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
});

WT_DECLARE_WT_MEMBER(InsertSeperatorInElectronMenu, Wt::JavaScriptFunction, "InsertSeperatorInElectronMenu",
                     function( selfid, pos )
{
  let m = $(window).data('electronMenu'+selfid);
  
  let appmenu = Menu.getApplicationMenu();
  if( !m || !appmenu ) return;
  
  if( pos >= 0 )
    m.submenu.insert( pos, new MenuItem({type: 'separator'}) );
  else
    m.submenu.append( new MenuItem({type: 'separator'}) );
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
});

WT_DECLARE_WT_MEMBER(HideElectronMenuItem, Wt::JavaScriptFunction, "HideElectronMenuItem",
  function( itemid, hidden )
{
  let appmenu = Menu.getApplicationMenu();
  if( !appmenu ) return;
  
  let eitem = $(window).data('electronItem'+itemid);
  if( eitem ){
    eitem.visible = !hidden;
  } else {
    console.log( 'Failed to get electronItem for ' + itemid );
  }
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
});


WT_DECLARE_WT_MEMBER(DisableElectronMenuItem, Wt::JavaScriptFunction, "DisableElectronMenuItem",
                     function( itemid, disabled )
{
  let appmenu = Menu.getApplicationMenu();
  if( !appmenu ) return;
  
  let eitem = $(window).data('electronItem'+itemid);
  if( eitem ){
    eitem.enabled = !disabled;
  } else {
    console.log( 'Failed to get electronItem for ' + itemid );
  }
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
});

/*
WT_DECLARE_WT_MEMBER(RemoveMenuItemFromElectronMenu, Wt::JavaScriptFunction, "RemoveMenuItemFromElectronMenu",
function( menuid, itemid )
{
  //https://github.com/electron/electron/issues/527
  
  let item = $(window).data('electronItem'+itemid);
  let parent = $(window).data('electronMenu'+menuid);
  if( !item ){
    console.log( 'RemoveMenuItemFromElectronMenu: Failed to get child, ' + itemid );
    return;
  }
  
  if( !parent ){
    console.log( 'RemoveMenuItemFromElectronMenu: Failed to get parent, ' + menuid );
    return;
  }
  
  var index = parent.submenu.items.indexOf(item);
  if (index > -1) {
    parent.submenu.items.splice(index, 1);
    $(window).data('electronMenu'+itemid,null);
  } else {
    console.log( 'DID NOt Got index to remove' );
  }
} );
*/

WT_DECLARE_WT_MEMBER(ClearElectronMenu, Wt::JavaScriptFunction, "ClearElectronMenu",
function( menuid )
{
  //https://github.com/electron/electron/issues/527
  let parent = $(window).data('electronMenu'+menuid);
  if( !parent ){
    console.log( 'ClearElectronMenu: Failed to get parent, ' + menuid );
    return;
  }
  parent.submenu.clear();
  
  let appmenu = Menu.getApplicationMenu();
  if( !appmenu ) return;
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
} );

WT_DECLARE_WT_MEMBER(AddCheckBoxItemToElectronMenu, Wt::JavaScriptFunction, "AddCheckBoxItemToElectronMenu",
function( menuid, itemid, txt, isChecked )
{
  let appmenu = Menu.getApplicationMenu();
  if(!appmenu) return;
  
  let m = $(window).data('electronMenu'+menuid);
  
  if( !m )
  {
    console.log( 'No electron menu set when trying to add checkbox item: ' + txt + ' menuid=' + menuid + ', itemid=' + itemid );
    return;
  }
  
  var newItem = new MenuItem({label: txt, type: 'checkbox', checked: isChecked, click: function(item, BrowserWindow){
    console.log( "Got click status: " + item.checked );
    Wt.emit(itemid,'electron_checked',item.checked);
  }} );
  
  m.submenu.append( newItem );
  
  $(window).data('electronItem'+itemid, newItem);
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
} );


WT_DECLARE_WT_MEMBER(AddElectronSubMenu, Wt::JavaScriptFunction, "AddElectronSubMenu",
function( menuid, submenuid, txt, iconPath )
{
  let appmenu = Menu.getApplicationMenu();
  if(!appmenu) return;
  
  //let m = $('#'+menuid).data('electronMenu');
  let m = $(window).data('electronMenu'+menuid);
  
  if( !m )
  {
    console.log( 'No electron menu for parent: ' + txt + ', menuid=' + menuid + ', submenuid=' + submenuid );
    return;
  }

  var newItem = new MenuItem({label: txt, icon: (iconPath.length>0 ? iconPath : null), submenu: []});
  m.submenu.append( newItem );
  
  $(window).data('electronMenu'+submenuid, newItem );
  //console.log( 'Added submenu for ' + submenuid + ', ' + newItem.label );
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
} );


WT_DECLARE_WT_MEMBER(AddRoleItemToElectronMenu, Wt::JavaScriptFunction, "AddRoleItemToElectronMenu",
function( menuid, roleName )
{
  let appmenu = Menu.getApplicationMenu();
  if( !appmenu ) return;
  
  let m = $(window).data('electronMenu'+menuid);
  if( m ){
    m.submenu.append( new MenuItem({role: roleName}) );
  } else {
    console.log( 'Failed to get electronMenu for ' + itemid + ' to implement ' + role );
  }
  
  if( $(window).data('HaveTriggeredMenuUpdate') )
    Menu.setApplicationMenu(appmenu);
});


//A hack to force the Electron Menu system to actually update its contents
WT_DECLARE_WT_MEMBER(TriggerElectronMenuUpdate, Wt::JavaScriptFunction, "TriggerElectronMenuUpdate",
function()
{
  //This function gets called once, after loading all the widgets into the
  //  GUI, and doing the initial defining of Electron menu items.
  //
  //Calling Menu.setApplicationMenu() for each item we add to the menu during
  //  program start up is pretty slow.  Takes initial load time from about
  //  1.5 seconds on macOS if we only do it once, to 3.4 seconds if we do it
  //  for every menu item addded, etc (first request Wt log prints out, to final
  //  request).  Since this function gets called once after initial load, we can
  //  defer all calls to Menu.setApplicationMenu() until now.  But after this
  //  we want changes to take effect immediately (ex enable/disable a menu item)
  //Note: requestNewCleanSession() will reset this variable to speedup resets.
  $(window).data('HaveTriggeredMenuUpdate',true);
  
  let appmenu = Menu.getApplicationMenu();
  if( appmenu )
    Menu.setApplicationMenu(appmenu);
});


//A hack since we apparetnyl need absolute path to icons, or else there is  a fata exception in the JS
string resolve_icon_path( string inputIconPath )
{
  if( inputIconPath.empty() )
    return "";
  
  const std::string docroot = (wApp ? wApp->docRoot() : string());
  string iconPath = SpecUtils::append_path( docroot, inputIconPath );
  
  if( !SpecUtils::is_file( iconPath ) )
    iconPath = SpecUtils::append_path( SpecUtils::get_working_path(), inputIconPath );
  
  if( !SpecUtils::is_file( iconPath ) )
  {
    cerr << "Couldnt find icon '" << inputIconPath
         << "' in docRoot ('" << docroot << "') or CWD ('"
         << SpecUtils::get_working_path() << "')" << endl;
    iconPath = "";
  }
  

#if( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  for (size_t pos = iconPath.find('\\'); pos != string::npos; pos = iconPath.find('\\'))
	  iconPath[pos] = '/';
#endif

  return iconPath;
}//string resolve_icon_path( string iconPath )
#endif



PopupDivMenu::PopupDivMenu( Wt::WPushButton *menuParent,
                            const PopupDivMenu::MenuType menutype )
: WPopupMenu(),
  m_parentItem( 0 ),
  m_menuParent( menuParent ),
#if(USE_OSX_NATIVE_MENU)
  m_nsmenu( 0 ),
#endif
#if( USING_ELECTRON_NATIVE_MENU )
  m_hasElectronCounterpart( false ),
#endif
  m_mobile( false ),
  m_type( menutype )
{
  addStyleClass( "wt-no-reparent" );
  setJavaScriptMember("wtNoReparent", "true");
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
  
  m_mobile = (app && app->isMobile());
  
  if( m_mobile )
  {
    addStyleClass( "PopupDivMenuPhone" );
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "ShowPhone", wtjsShowPhone);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "HideOverlay", wtjsHideOverlay);
  }else
  {
    addStyleClass( "PopupDivMenu" );
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsBringAboveDialogs);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsOnPopupDivShow);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsAdjustTopPos);
    
    setAutoHide( true, 500 );
    implementStateless( &PopupDivMenu::showFromClick, &PopupDivMenu::doHide );
    implementStateless( &PopupDivMenu::showFromMouseOver, &PopupDivMenu::doHide );
    implementStateless( &PopupDivMenu::doHide, &PopupDivMenu::showFromClick );
  }//if( mobile ) / else

#if( USING_ELECTRON_NATIVE_MENU || USE_OSX_NATIVE_MENU )
  const bool useNativeMenu = InterSpecApp::isPrimaryWindowInstance();
#endif
  
#if( USING_ELECTRON_NATIVE_MENU )
  if( useNativeMenu )
  {
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsFindElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsAddMenuItemToElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsAddSeperatorToElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsInsertSeperatorInElectronMenu);
    //LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsRemoveMenuItemFromElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsClearElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsAddCheckBoxItemToElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsAddElectronSubMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsAddRoleItemToElectronMenu);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsTriggerElectronMenuUpdate);
  }//if( useNativeMenu )
#endif
  
  if( menuParent )
  {
    if( m_mobile)
    {
      if( m_type != TransientMenu )
        menuParent->clicked().connect( this, &PopupDivMenu::showFromMouseOver );
    }else
    {
#if(USE_OSX_NATIVE_MENU)
      if( useNativeMenu && (menutype == AppLevelMenu) )
      {
        const string buttontxt = menuParent->text().toUTF8();
        if( buttontxt.length() )
          m_nsmenu = addOsxMenu( this, buttontxt.c_str() );
      } //AppLevelMenu
#endif
      
#if( USING_ELECTRON_NATIVE_MENU )
      if( useNativeMenu && (menutype == AppLevelMenu) )
      {
        //Here we are adding a menu at the topof the window.  "InterSpec", "View", "Tools", "Help"
        //In Js, we should find the menu, and get a reference to it.
        //app->doJavaScript( "console.log('Adding PopupDivMenu id=" + id() + " parentTxt=" + menuParent->text().toUTF8() + "');" );
        app->doJavaScript( "Wt.WT.FindElectronMenu('" + menuParent->text().toUTF8() + "', '" + id() + "');" );
        m_hasElectronCounterpart = true;
      }
#endif
      
      menuParent->setMenu( this );  //adds/removes active class to menuParent, etc
      menuParent->clicked().connect( this, &PopupDivMenu::showFromClick );
      menuParent->addStyleClass( "PopupMenuParentButton" );
      menuParent->mouseWentOver().connect( this, &PopupDivMenu::showFromMouseOver );
    } //else !m_mobile
  }//if( menuParent )
  else
  {
#if( USING_ELECTRON_NATIVE_MENU )
    if( menutype==AppLevelMenu )
      doJavaScript( "console.log('Not implemented: Adding PopupDivMenu (no parent) id=" + id() + (menutype==AppLevelMenu?" AppLevelMenu":" TransientMenu") + "');" );
#endif
  }
  // Note: need to call this to reset the menus when a submenu is
  //  selected/hidden.  No need to call triggered(), as it will be hidden too.
  aboutToHide().connect( this, &PopupDivMenu::doHide );
}//PopupDivMenu( constructor )



PopupDivMenu::~PopupDivMenu()
{
#if( USE_OSX_NATIVE_MENU )
//  if( m_nsmenu )
//    removeOsxMenu( m_nsmenu );
#endif

#if( USING_ELECTRON_NATIVE_MENU )
  if(  m_hasElectronCounterpart )
  {
    if( m_parentItem )
    {
  	//This is a sub-menu, so we should get rid of the electron sub-menu tooss
      //WApplication::instance()->doJavaScript( "Wt.WT.HideElectronMenuItem('" + id() + "',true);" );
      //WApplication::instance()->doJavaScript( "Wt.WT.ClearElectronMenu('" + id() + "');" );
    }else if( m_menuParent )
    {
    }
  }//if(  m_hasElectronCounterpart )
#endif
}

Wt::WMenuItem *PopupDivMenu::addSeparatorAt( int index )
{
  WMenuItem *item = WPopupMenu::addSeparator();
  
  if( (index >= 0) && (indexOf(item) != index) )
  {
    removeItem( item );
    WMenu::insertItem( index, item );
  }
  
#if( USING_ELECTRON_NATIVE_MENU )
  if( m_hasElectronCounterpart && (m_type == AppLevelMenu) )
  {
    //Need to modify Wt.WT.AddSeperatorToElectronMenu to also take the index
    WApplication::instance()->doJavaScript( "Wt.WT.InsertSeperatorInElectronMenu('" + id() + "'," + std::to_string(index) + ");" );
  }
#elif( USE_OSX_NATIVE_MENU )
  if( m_nsmenu )
  {
    void *macItem = addOsxSeparatorAt( index, m_nsmenu );
    item->setData( macItem );
  }
#else
#endif
  
  return item;
}


bool PopupDivMenu::removeSeperator( Wt::WMenuItem *sepertor )
{
  if( !sepertor )
    return false;
  
  int index = indexOf( sepertor );
  if( index < 0 )
  {
    cout << "PopupDivMenu::removeSeperator: !Couldnt find sepertor" << endl;
    return false;
  }
  removeItem( sepertor );
  
#if( USING_ELECTRON_NATIVE_MENU )
  //ToDo: Hack, just hide the electron menu since we cant easily remove items...
  if( m_hasElectronCounterpart )
    WApplication::instance()->doJavaScript( "Wt.WT.HideElectronMenuItem('" + sepertor->id() + "',true);" );
#elif( USE_OSX_NATIVE_MENU )
  if( m_nsmenu )
    removeOsxSeparator( m_nsmenu, sepertor->data() );
#endif

  return true;
}//bool removeSeperator( Wt::WMenuItem *sepertor )



Wt::WMenuItem *PopupDivMenu::addSeparator()
{
  return addSeparatorAt( -1 );
}//void PopupDivMenu::addSeparator(void)


#if( USING_ELECTRON_NATIVE_MENU )
void PopupDivMenu::addRoleMenuItem( MenuRole role )
{
  if( !m_hasElectronCounterpart )
    return;
  
  string rolestr;
  switch( role )
  {
    case MenuRole::Quit:             rolestr = "quit";             break;
    case MenuRole::ResetZoom:        rolestr = "resetzoom";        break;
    case MenuRole::ZoomIn:           rolestr = "zoomin";           break;
    case MenuRole::ZoomOut:          rolestr = "zoomout";          break;
    case MenuRole::ToggleFullscreen: rolestr = "togglefullscreen"; break;
    case MenuRole::Cut:              rolestr = "cut";              break;
    case MenuRole::Copy:             rolestr = "copy";             break;
    case MenuRole::Past:             rolestr = "paste";            break;
    case MenuRole::ToggleDevTools:   rolestr = "toggledevtools";   break;
#if defined(__APPLE__)
    case MenuRole::Hide:             rolestr = "hide";             break;
    case MenuRole::HideOthers:       rolestr = "hideothers";       break;
    case MenuRole::UnHide:           rolestr = "unhide";           break;
    case MenuRole::Front:            rolestr = "front";            break;
#endif
  }//switch( role )
  
  WApplication::instance()->doJavaScript( "Wt.WT.AddRoleItemToElectronMenu('" + id() + "','" + rolestr + "');" );
}//addRoleMenuItem( MenuRole role )


void PopupDivMenu::triggerElectronMenuUpdate()
{
  const bool useNativeMenu = InterSpecApp::isPrimaryWindowInstance();
  if( useNativeMenu )
    WApplication::instance()->doJavaScript( "Wt.WT.TriggerElectronMenuUpdate();" );
}//void PopupDivMenu::triggerElectronMenuUpdate()


void PopupDivMenu::clearElectronMenu()
{
  if( m_hasElectronCounterpart )
    WApplication::instance()->doJavaScript( "Wt.WT.ClearElectronMenu('" + id() + "');" );
}//void clearElectronMenu();


#if defined(__APPLE__)
PopupDivMenuItem *PopupDivMenu::createAboutThisAppItem()
{
  PopupDivMenuItem *item = new PopupDivMenuItem( "About", "" );
  addItem( item );
  
  if( !m_hasElectronCounterpart )
    return;
  
  item->m_electron_clicked.connect( item, &PopupDivMenuItem::emitClickFromElectronMenu );
  
  //I cant seem to figure out how to use use the native "role" property
  //  of 'about' and still costumize behaviour, so I'll hack it for the moment.
  //const string role = "'about'";
  //const string label = "null";
  const string role = "null";
  const string label = "'About InterSpec'";
  WApplication::instance()->doJavaScript( "Wt.WT.AddMenuItemToElectronMenu("
                                          "'" + id() + "'," + label + ",null,"
                                          + "'" + item->id() + "', " + role + ");" );
  return item;
}//createAboutThisAppItem()
#endif
#endif

void PopupDivMenu::setHidden( bool hidden, const Wt::WAnimation &animation )
{
  //WPopupMenu sets some javascript stuff I'de rather not meddle with for the
  //  phone version of the menus
  if( m_mobile )
  {
    if( hidden )
    {
      if( hasStyleClass("Root") )
      {
        doJavaScript("$('.mobilePopupMenuOverlay').hide();");
        doJavaScript("$(document).off(\"touchstart.mobileOverlay\");"); //unbind event listener
      }//hasStyleClass("Root") -- is the top parent
    }else
    {
      doJavaScript("$('.mobilePopupMenuOverlay').show();");
    }//if( hidden ) / else
    
    WCompositeWidget::setHidden(hidden, animation);
  }else
  {
    WPopupMenu::setHidden( hidden, animation );
  }//if( m_mobile ) / else
}//setHidden()


void PopupDivMenu::showFromClick()
{
  doShow( true );
}

void PopupDivMenu::showFromMouseOver()
{
  doShow( false );
}


void PopupDivMenu::doShow( bool clicked )
{
  if( !m_menuParent )
    return;

  if( !m_mobile )
    doJavaScript( "Wt.WT.OnPopupDivShow('" + id() + "','" + m_menuParent->id() + "');" );

  //Note/Hack: This is needed to check if clicked or mouseover.  If clicked,
  //  calling popup() will somehow cause the popupmenu to not show.
  //  I know, weird huh.
  if( !clicked )
    popup( m_menuParent, Wt::Vertical );

  if( m_mobile )
  {
    doJavaScript("$('.mobilePopupMenuOverlay').show();");
    doJavaScript( "Wt.WT.ShowPhone('" + id() + "');" );
    doJavaScript( "Wt.WT.HideOverlay();" );
  }//if( m_mobile )
}//void doShow()


void PopupDivMenu::doHide()
{
  if( m_menuParent )
    m_menuParent->removeStyleClass( "active" );
  
  if( m_mobile )
    setHidden(true, WAnimation(WAnimation::SlideInFromLeft,WAnimation::Linear,200));
  else
    hide();
}//void doHide()


PopupDivMenuItem *PopupDivMenu::addWidget( Wt::WWidget *widget,
                                           const bool closeMenuOnClick )
{
  PopupDivMenuItem *item = addMenuItem( "", "", false );

  if( !closeMenuOnClick )
  {
    WInteractWidget *w = dynamic_cast<WInteractWidget *>( widget );
    if( w )
      w->clicked().preventPropagation();
    item->clicked().preventPropagation();
    item->setSelectable( false );
  }
  
  if( !dynamic_cast<WCheckBox *>(widget) )
    item->addStyleClass( "PopupDivMenuWidget" );
  
  item->addWidget( widget );
  
#if(USE_OSX_NATIVE_MENU)
  if( m_nsmenu )
  {
    WCheckBox *cb = static_cast<WCheckBox *>( widget );
    if( cb )
    {
      removeOsxMenuItem( item->m_nsmenuitem, m_nsmenu );
      item->m_nsmenuitem = addOsxCheckableMenuItem( m_nsmenu, cb, item );
    }else
    {
      cerr << "PopupDivMenu::addWidget: Unsuppored Widget type on OS X" << endl;
    }
  }//if( m_nsmenu )
#endif
  
#if( USING_ELECTRON_NATIVE_MENU )
  if( m_hasElectronCounterpart && (m_type == PopupDivMenu::AppLevelMenu) )
  {
    WCheckBox *cb = static_cast<WCheckBox *>( widget );
    
    if( cb )
    {
      item->m_electron_checked.connect( item, &PopupDivMenuItem::toggleFromElectronMenu );
      
      //Theres a hack where we dont add electron menu items for blank labels, so we
      //  shouldnt need to remove item
      //WApplication::instance()->doJavaScript( "Wt.WT.RemoveMenuItemFromElectronMenu('" + id() + "',"
        //                                       + "'" + item->id() + "'"
          //                                     + ");" );
      //cerr << "ADDing check box for " << item->id() << ", " << cb->text().toUTF8() << endl;
      
      WApplication::instance()->doJavaScript( "Wt.WT.AddCheckBoxItemToElectronMenu('" + id() + "',"
                                             + "'" + item->id() + "',"
                                             + "'" + cb->text().toUTF8() + "',"
                                             + (cb->isChecked()?"true":"false")
                                             + ");" );
	  item->m_hasElectronItem = true;
    }
  }
#endif
  
  return item;
}//PopupDivMenuItem *addWidget( WWidget *widget, const bool closeMenuOnClick )


namespace
{
  void doTriggeredEmit( PopupDivMenuItem *item )
  {
    item->triggered().emit(item);
  }
}

PopupDivMenuItem *PopupDivMenu::addMenuItem( const Wt::WString &text,
                                             const std::string &iconPath,
                                            const bool closeMenuOnActivation )
{
  return insertMenuItem( -1, text, iconPath, closeMenuOnActivation );
}//addMenuItem(...)



PopupDivMenuItem *PopupDivMenu::insertMenuItem( const int index,
                                     const Wt::WString &text,
                                     const std::string &iconPath,
                                     const bool closeMenuOnActivation )
{
  PopupDivMenuItem *item = new PopupDivMenuItem( text, iconPath );
  
  if( !closeMenuOnActivation )
  {
    item->setSelectable( false );
    item->clicked().connect( boost::bind(&doTriggeredEmit, item) );
    item->clicked().preventPropagation();
    Wt::WAnchor *a = item->anchor();
    if( a )
    {
      item->anchor()->clicked().preventPropagation();
      item->anchor()->clicked().connect( boost::bind(&doTriggeredEmit, item) );
    }
  }
  

  WMenu::insertItem( ((index >= 0) ? index : count()), item );
  
  //addItem( item );
  
  if( m_mobile && closeMenuOnActivation )
    item->triggered().connect( this, &PopupDivMenu::hideMenuAndParents );
  
#if(USE_OSX_NATIVE_MENU)
  if( m_nsmenu )
  {
    item->m_nsmenuitem = insertOsxMenuItem( m_nsmenu, item, index );
    item->m_nsmenu = m_nsmenu;
    item->setData( item->m_nsmenuitem );
  }
#endif
  
#if( USING_ELECTRON_NATIVE_MENU )
  if( m_hasElectronCounterpart && (!text.empty() || !iconPath.empty())
     && (m_type == PopupDivMenu::AppLevelMenu) )
  {
    item->m_electron_clicked.connect( item, &PopupDivMenuItem::emitClickFromElectronMenu );
    //Need to edit Wt.WT.AddMenuItemToElectronMenu to take into accoutn the index were adding it at
     WApplication::instance()->doJavaScript( "Wt.WT.AddMenuItemToElectronMenu('" + id() + "',"
                                           + "'" + text.toUTF8() + "',"
                                           + "'" + resolve_icon_path(iconPath) + "',"
                                           + "'" + item->id() + "', null);" );
    item->m_hasElectronItem = true;
  }
#endif
  
  return item;
}//

Wt::WMenuItem *PopupDivMenu::parentItem()
{
  return m_parentItem;
}//Wt::WMenuItem *PopupDivMenu::parentItem()


void PopupDivMenu::hideMenuAndParents()
{
  doHide();
  if( m_parentItem )
  {
    PopupDivMenu *p = dynamic_cast<PopupDivMenu *>( m_parentItem->parentMenu() );
    if( p )
      p->hideMenuAndParents();
  }
}//void hideMenuAndParents()


bool PopupDivMenu::isHidden() const
{
  //See notes in header for how big of a hack this function is.
  return m_mobile || WPopupMenu::isHidden();
}


PopupDivMenuItem *PopupDivMenu::addPhoneBackItem( PopupDivMenu *parent )
{
  const char *txt = parent ? "Previous" : "Close";
  const char *icon = parent ? "InterSpec_resources/images/back-alt-16.png" : "InterSpec_resources/images/menuclose.png";
  PopupDivMenuItem *backitem = new PopupDivMenuItem( txt, icon );
  
  addItem( backitem );
  
  if( parent )
  {
    backitem->triggered().connect( this, &PopupDivMenu::doHide );
    backitem->addStyleClass( "PhoneMenuBack" );
  }else
  {
    backitem->triggered().connect( this, &PopupDivMenu::doHide );
    backitem->clicked().connect( this, &PopupDivMenu::doHide );
    
    backitem->addStyleClass( "PhoneMenuClose" );
    addStyleClass("Root");
  }
  
  return backitem;
}//PopupDivMenuItem *addPhoneBackItem()


PopupDivMenu *PopupDivMenu::addPopupMenuItem( const Wt::WString &text,
                                              const std::string &iconPath )
{
  PopupDivMenu *menu = new PopupDivMenu( 0, m_type );
  
  if( m_mobile )
  {
    menu->m_menuParent = m_menuParent;
    menu->m_parentItem = addItem( text );
    menu->m_parentItem->clicked().preventPropagation();
    menu->m_parentItem->clicked().preventDefaultAction();
    menu->m_parentItem->setStyleClass( "submenu" );
    menu->m_parentItem->setMargin(10,Wt::Right);
    menu->m_parentItem->triggered().connect( menu, &PopupDivMenu::showFromMouseOver );
    menu->addPhoneBackItem( this );
  }else
  {
    menu->m_parentItem = addMenu( iconPath, text, menu );
    menu->m_parentItem->clicked().preventPropagation();
    menu->m_parentItem->clicked().preventDefaultAction();
    
    string js;
    js = "function(){Wt.WT.BringAboveDialogs('" + menu->id() + "');"
         "Wt.WT.AdjustTopPos('" + menu->id() + "');}";
    menu->m_parentItem->mouseWentOver().connect( js );

#if(USE_OSX_NATIVE_MENU)
    PopupDivMenu *p = parentItem()
                      ? dynamic_cast<PopupDivMenu *>(parentItem()->menu())
                      : 0;
    
    if( m_nsmenu && (!p || p->m_type==PopupDivMenu::AppLevelMenu) )
      menu->m_nsmenu = addOsxSubMenu( m_nsmenu, menu, text.toUTF8().c_str() );
#endif
    
#if( USING_ELECTRON_NATIVE_MENU )
    if( m_hasElectronCounterpart )
    {
      PopupDivMenu *p = parentItem() ? dynamic_cast<PopupDivMenu *>(parentItem()->menu()) : 0;
      if( !p || p->m_type==PopupDivMenu::AppLevelMenu )
      {
        WApplication::instance()->doJavaScript( "Wt.WT.AddElectronSubMenu('" + id() + "',"
                                                + "'" + menu->id() + "',"
                                                + "'" + text.toUTF8() + "',"
                                                + "'" + resolve_icon_path(iconPath) + "'"
                                                + ");" );
      }
    }
#endif
  }
  
  return menu;
}//addPopupMenuItem( const Wt::WString &text, const std::string iconPath )



PopupDivMenuItem::PopupDivMenuItem( const Wt::WString &text,
                                    const std::string &iconPath )
  : WMenuItem( iconPath, text, 0, WMenuItem::PreLoading )
#if( USE_OSX_NATIVE_MENU )
   , m_nsmenu( 0 )
   , m_nsmenuitem( 0 )
#endif
#if( USING_ELECTRON_NATIVE_MENU )
   , m_hasElectronItem( false )
   , m_electron_clicked( this, "electron_clicked", false )
   , m_electron_checked( this, "electron_checked", false )
#endif

{
  WAnchor *a = anchor();
  if( a )
  {
    a->clicked().preventPropagation();
    clicked().connect( this, &PopupDivMenuItem::nonAnchorClickHack );
  }//if( a )
  
#if( USING_ELECTRON_NATIVE_MENU )
  const bool useNativeMenu = InterSpecApp::isPrimaryWindowInstance();
  if( useNativeMenu )
  {
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsHideElectronMenuItem);
    LOAD_JAVASCRIPT(wApp, "PopupDiv.cpp", "PopupDivMenu", wtjsDisableElectronMenuItem);
  }
#endif
}//PopupDivMenuItem constructor


PopupDivMenuItem::~PopupDivMenuItem()
{
#if( USE_OSX_NATIVE_MENU )
  if( m_nsmenu && m_nsmenuitem )
    removeOsxMenuItem( m_nsmenuitem, m_nsmenu );
#endif
  
#if( USING_ELECTRON_NATIVE_MENU )
//Right now, as a workaround for not being able to remove items from an Electron
//  menu (https://github.com/electron/electron/issues/527), we will hide the
//  item.  See clearElectronMenu().
  if( m_hasElectronItem )
    WApplication::instance()->doJavaScript( "Wt.WT.HideElectronMenuItem('" + id() + "',true);" );
#endif
}//~PopupDivMenuItem()


void PopupDivMenuItem::nonAnchorClickHack()
{
  if( !isCheckable() )
    select();
}//void PopupDivMenuItem::nonAnchorClickHack()



#if( USE_OSX_NATIVE_MENU )
void *PopupDivMenuItem::getNsMenuItem()
{
  return m_nsmenuitem;
}

void PopupDivMenuItem::setHidden( bool hidden, const Wt::WAnimation &animation )
{
  WMenuItem::setHidden( hidden, animation );
  if( m_nsmenuitem )
    setOsxMenuItemHidden( m_nsmenuitem, hidden );
}
#endif

#if( USING_ELECTRON_NATIVE_MENU )
void PopupDivMenuItem::emitClickFromElectronMenu()
{
  cout << "emitClickFromElectronMenu" << endl;
  this->triggered().emit( (Wt::WMenuItem *)this );
  wApp->triggerUpdate();
}//void emitClickFromElectronMenu( PopupDivMenuItem *item )

// #include <Wt/WServer>

void PopupDivMenuItem::toggleFromElectronMenu(bool checked)
{
  //Wt::WServer::instance()->post( WApplication::instance()->sessionId(), [this,checked](){
  
  cout << "toggleFromElectronMenu: " << checked << endl;
  WCheckBox *cb = checkBox();
  if( !cb )
  {
    cerr << "PopupDivMenuItem::toggleFromElectronMenu() called for non checkBox" << endl;
    return;
  }
  //const bool wasChecked = cb->isChecked();
  
  //For some reason the below doesnt seem to have an effect... not sure why yet
  cb->setChecked(checked);
  if( checked )
    cb->checked().emit();
  else
    cb->unChecked().emit();
  
  //if( checked != wasChecked )
    //cb->changed().emit();
  
  triggered().emit( (Wt::WMenuItem *)this );
  
  wApp->triggerUpdate();
//  } );
}//

void PopupDivMenuItem::setHidden( bool hidden, const Wt::WAnimation &animation )
{
  //cout << "Hiding " << text().toUTF8() << endl;
  WMenuItem::setHidden( hidden, animation );
  WApplication::instance()->doJavaScript( "Wt.WT.HideElectronMenuItem('" + id() + "',"
                                         + string(hidden ? "true" : "false") + ");" );
}

void PopupDivMenuItem::setDisabled(bool disabled)
{
  //cout << "DIsabling " << text().toUTF8() << endl;
  WMenuItem::setDisabled(disabled);
  WApplication::instance()->doJavaScript( "Wt.WT.DisableElectronMenuItem('" + id() + "',"
                                         + string(disabled ? "true" : "false") + ");" );
}

#endif





void PopupDivMenuItem::makeTextXHTML()
{
  WAnchor *a = anchor();
  if( a )
  {
    const int nkid = a->count();
    for( int i = 0; i < nkid; ++i )
    {
      WLabel *label = dynamic_cast<WLabel *>(a->widget(i));
      if( label )
      {
        label->setTextFormat( Wt::XHTMLText );
        return;
      }
    }
  }//if( a )
  
#if( PERFORM_DEVELOPER_CHECKS )
  //Making text XHTML is kinda a hack that could be Wt version specific/brittle.
  char buffer[512];
  snprintf( buffer, sizeof(buffer), "Unable to set menu item to XHTMLText."
           " Text='%s'.", text().toUTF8().c_str() );
  log_developer_error( __func__, buffer );
#endif
}//void makeTextXHTML()


Wt::WCheckBox *PopupDivMenuItem::checkBox()
{
  if( !isCheckable() )
  {
    //check if the user passed in a WCheckbox to PopupDivMenu::addWidget(...)
    const int nchild = count();
    for( int i = 0; i < nchild; ++i )
    {
      WCheckBox *cb = dynamic_cast<WCheckBox *>( widget(i) );
      if( cb )
        return cb;
    }//for( int i = 0; i < nchild; ++i )
    
    return 0;
  }//if( !isCheckable() )
  
  for (int i = 0; i < count(); ++i)
  {
    WAnchor *anchor = dynamic_cast<WAnchor *>(widget(i));
    if( anchor && anchor->children().size() )
      return dynamic_cast<WCheckBox *>(anchor->widget(0));
  }//for (int i = 0; i < count(); ++i)
  
  throw runtime_error( "Serious logic error in PopupDivMenuItem::checkBox()" );
  return 0;
}//checkBox()


Wt::WAnchor *PopupDivMenuItem::anchor()
{
  for (int i = 0; i < count(); ++i) {
    WAnchor *result = dynamic_cast<WAnchor *>(widget(i));
    if( result )
      return result;
  }

  return 0;
}//Wt::WAnchor *PopupDivMenuItem::anchor()



