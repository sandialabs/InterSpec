WT_DECLARE_WT_MEMBER
(TitleBarChangeMaximized, Wt::JavaScriptFunction, "TitleBarChangeMaximized",
function(maximized){
  if( maximized ){
    document.querySelectorAll( '.resizer.top,.resizer.left' ).forEach( function( el ){ el.style.display = 'none'; } );
    document.querySelectorAll( '.window-max-restore' ).forEach( function( el ){ el.classList.remove( 'window-maximize' ); el.classList.add( 'window-unmaximize' ); } );
  }else{
    document.querySelectorAll( '.resizer.top,.resizer.left' ).forEach( function( el ){ el.style.display = ''; } );
    document.querySelectorAll( '.window-max-restore' ).forEach( function( el ){ el.classList.remove( 'window-unmaximize' ); el.classList.add( 'window-maximize' ); } );
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
    window.addEventListener( 'blur', function(){
      var tb = document.querySelector( '.app-titlebar' );
      if( tb ) tb.classList.add( 'inactive' );
    } );

    window.addEventListener( 'focus', function(){
      var tb = document.querySelector( '.app-titlebar' );
      if( tb ) tb.classList.remove( 'inactive' );
    } );
                                  
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
  window.addEventListener( 'blur', function(){
    var tb = document.querySelector( '.app-titlebar' );
    if( tb ) tb.classList.add( 'inactive' );
  } );

  window.addEventListener( 'focus', function(){
    var tb = document.querySelector( '.app-titlebar' );
    if( tb ) tb.classList.remove( 'inactive' );
  } );
                                
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
    document.body.style.zoom = 1.0;
  }

  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(ResetPageZoom, Wt::JavaScriptFunction, "ResetPageZoom",
function(){
  document.body.style.zoom = 1.0;
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
    if( document.body.style.zoom )
      current = parseFloat( document.body.style.zoom );
    document.body.style.zoom = 1.05 * current;
  }
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(IncreasePageZoom, Wt::JavaScriptFunction, "IncreasePageZoom",
function(){
  let current = 1.0;
  if( document.body.style.zoom )
    current = parseFloat( document.body.style.zoom );
  document.body.style.zoom = 1.05 * current;
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
    if( document.body.style.zoom )
      current = parseFloat( document.body.style.zoom );
    document.body.style.zoom = 0.95 * current;
  }
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);

#else  //BUILD_AS_ELECTRON_APP

WT_DECLARE_WT_MEMBER
(DecreasePageZoom, Wt::JavaScriptFunction, "DecreasePageZoom",
function(){
  let current = 1.0;
  if( document.body.style.zoom )
    current = parseFloat( document.body.style.zoom );
  document.body.style.zoom = 0.95 * current;
  setTimeout( function(){ window.dispatchEvent(new Event('resize')); }, 100 );
}
);
#endif //BUILD_AS_ELECTRON_APP
