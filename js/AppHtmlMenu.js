WT_DECLARE_WT_MEMBER
(TitleBarChangeMaximized, Wt::JavaScriptFunction, "TitleBarChangeMaximized",
function(maximized){
  if( maximized ){
    $('.resizer.top,.resizer.left').hide();
    $('.window-max-restore').removeClass('window-maximize').addClass('window-unmaximize');
  }else{
    $('.resizer.top,.resizer.left').show();
    $('.window-max-restore').removeClass('window-unmaximize').addClass('window-maximize');
  }

  let f = function(){ window.dispatchEvent(new Event('resize')); };
  setTimeout( f, 100 );
  setTimeout( f, 500 );
}
);


#if( BUILD_AS_ELECTRON_APP )
WT_DECLARE_WT_MEMBER
(SetupAppTitleBar, Wt::JavaScriptFunction, "SetupAppTitleBar",
function(){
  if( Wt.WT.IsElectronInstance() ){
    let currentWindow = remote.getCurrentWindow();
    $(window).data('ElectronWindow',currentWindow);
                                      
    Wt.WT.TitleBarChangeMaximized( currentWindow.isMaximized() );
                                      
    currentWindow.on( 'blur', function(){ $('.app-titlebar').addClass('inactive'); } );
    currentWindow.on( 'focus', function(){ $('.app-titlebar').removeClass('inactive'); } );
    currentWindow.on( 'unmaximize', function(){ Wt.WT.TitleBarChangeMaximized(false); } );
    currentWindow.on( 'maximize', function(){ Wt.WT.TitleBarChangeMaximized(true); } );
    currentWindow.on( 'leave-full-screen', function(){ console.log('currentWindow.leave-full-screen'); } );
    currentWindow.on( 'enter-full-screen', function(){ console.log('currentWindow.enter-full-screen'); } );
  }else{
    $(window).blur(function(){
      $('.app-titlebar').addClass('inactive');
    });
  
    $(window).focus(function(){
      $('.app-titlebar').removeClass('inactive');
    });
                                  
    document.addEventListener('fullscreenchange', (event) => {
      if (document.fullscreenElement) {
        console.log(`Element: ${document.fullscreenElement.id} entered full-screen mode.`);
      } else {
        console.log('Leaving full-screen mode.');
      }
      Wt.WT.TitleBarChangeMaximized(document.fullscreenElement);
    });
  }
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(SetupAppTitleBar, Wt::JavaScriptFunction, "SetupAppTitleBar",
function(){
  $(window).blur(function(){
    $('.app-titlebar').addClass('inactive');
  });

  $(window).focus(function(){
    $('.app-titlebar').removeClass('inactive');
  });
                                
  document.addEventListener('fullscreenchange', (event) => {
    if (document.fullscreenElement) {
      console.log(`Element: ${document.fullscreenElement.id} entered full-screen mode.`);
    } else {
      console.log('Leaving full-screen mode.');
    }
    Wt.WT.TitleBarChangeMaximized(document.fullscreenElement);
  });
}
);
#endif //BUILD_AS_ELECTRON_APP



#if( BUILD_AS_ELECTRON_APP )

WT_DECLARE_WT_MEMBER
(TitleBarHandleMaximizeClick, Wt::JavaScriptFunction, "TitleBarHandleMaximizeClick",
function(){
  if( Wt.WT.IsElectronInstance() ){
  let win = $(window).data('ElectronWindow');
  if( win.isMaximized() ) {
    win.unmaximize();
    Wt.WT.TitleBarChangeMaximized(false);
  } else {
    win.maximize();
    Wt.WT.TitleBarChangeMaximized(true);
  }
  }else{
    let elem = document.querySelector(".Wt-domRoot");
    if (!document.fullscreenElement) {
      elem.requestFullscreen().catch(err => {
        console.log( 'Error attempting to enable full-screen mode' );
    });
    } else {
      document.exitFullscreen();
    }
  }
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(TitleBarHandleMaximizeClick, Wt::JavaScriptFunction, "TitleBarHandleMaximizeClick",
function(){
  let elem = document.querySelector(".Wt-domRoot");
  if (!document.fullscreenElement) {
    elem.requestFullscreen().catch(err => {
      console.log( 'Error attempting to enable full-screen mode' );
    });
  } else {
      document.exitFullscreen();
  }
}
);
#endif //BUILD_AS_ELECTRON_APP



#if( BUILD_AS_ELECTRON_APP )
WT_DECLARE_WT_MEMBER
(ResetPageZoom, Wt::JavaScriptFunction, "ResetPageZoom",
function(){
  if( Wt.WT.IsElectronInstance() ){
    //let win = $(window).data('ElectronWindow');
    //win.webContents.setZoomFactor(1);
    let webFrame = require('electron').webFrame;
    webFrame.setZoomFactor(1.0);
  }else{
    $('body').css('zoom', 1.0);
  }
  
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(ResetPageZoom, Wt::JavaScriptFunction, "ResetPageZoom",
function(){
  $('body').css('zoom', 1.0);
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);
#endif //BUILD_AS_ELECTRON_APP



#if( BUILD_AS_ELECTRON_APP )
WT_DECLARE_WT_MEMBER
(IncreasePageZoom, Wt::JavaScriptFunction, "IncreasePageZoom",
function(){
  if( Wt.WT.IsElectronInstance() ){
    //let win = $(window).data('ElectronWindow');
    let webFrame = require('electron').webFrame;
    let level = webFrame.getZoomLevel();
    webFrame.setZoomLevel( level + 1 );
  }else{
    let current = 1.0;
    if( $('body').css('zoom') ) 
      current = parseFloat($('body').css('zoom'));
    $('body').css('zoom', 1.05*current );
  }
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(IncreasePageZoom, Wt::JavaScriptFunction, "IncreasePageZoom",
function(){
  let current = 1.0;
  if( $('body').css('zoom') ) 
    current = parseFloat($('body').css('zoom'));
  $('body').css('zoom', 1.05*current );
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);
#endif //BUILD_AS_ELECTRON_APP


#if( BUILD_AS_ELECTRON_APP )
WT_DECLARE_WT_MEMBER
(DecreasePageZoom, Wt::JavaScriptFunction, "DecreasePageZoom",
function(){
  if( Wt.WT.IsElectronInstance() ){
    //let win = $(window).data('ElectronWindow');
    let webFrame = require('electron').webFrame;
    let level = webFrame.getZoomLevel();
    webFrame.setZoomLevel( level - 1 );
  }else{
    let current = 1.0;
    if( $('body').css('zoom') ) 
      current = parseFloat($('body').css('zoom'));
    $('body').css('zoom', 0.95*current);
  }
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(DecreasePageZoom, Wt::JavaScriptFunction, "DecreasePageZoom",
function(){
  let current = 1.0;
  if( $('body').css('zoom') ) 
    current = parseFloat($('body').css('zoom'));
  $('body').css('zoom', 0.95*current);
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);
#endif //BUILD_AS_ELECTRON_APP


#if( BUILD_AS_ELECTRON_APP )
WT_DECLARE_WT_MEMBER
(ToggleDevTools, Wt::JavaScriptFunction, "ToggleDevTools",
function(){
  if( Wt.WT.IsElectronInstance() ){
    //let win = $(window).data('ElectronWindow');
    $(window).data('ElectronWindow').toggleDevTools();
  }else{
    alert( 'Can only show developer tools in primary instance of application' );
  }
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);
#endif  //BUILD_AS_ELECTRON_APP