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
    // Nothing to do here - nodejs will call back later to set things up
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
(ResetPageZoom, Wt::JavaScriptFunction, "ResetPageZoom",
function(){
  if( Wt.WT.IsElectronInstance() ){
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
