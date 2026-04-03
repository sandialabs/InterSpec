
WT_DECLARE_WT_MEMBER
(FileUploadFcn, Wt::JavaScriptFunction, "FileUploadFcn",
function()
{
  const target = document.querySelector('.Wt-domRoot');
  
  var urlIdFromName = function(name)
  {
    var str = name.toLowerCase();
    if( str.indexOf('i-')===0 || str.indexOf('fore')>-1 || str.indexOf('item')>-1
        || str.indexOf('ipc')>-1 || str.indexOf('unk')>-1 || str.indexOf('inter')>-1 || str.indexOf('ioi')>-1  )
      return 'ForegroundUpUrl';
    if( str.indexOf('b-') === 0 || str.indexOf('back')>-1 || str.indexOf('bkg')>-1 || str.indexOf('bgd')>-1 || str.indexOf('bg')>-1 )
      return 'BackgroundUpUrl';
    if( str.indexOf('k-') === 0  || str.indexOf('kwn') > -1  || str.indexOf('known')>-1 || str.indexOf('cal')>-1 || str.indexOf('check')>-1 )
      return 'SecondUpUrl';
    return null;
  };
  
  var uploadFcn = function(evt,urlid){
    try
    {
      evt.preventDefault();
      evt.stopPropagation();
      
      const files = evt.dataTransfer.files;
      if( !files || !files.length )
      {
        console.log( "uploadFcn: no files for some reason." );
        return;
      }
      
      if( window._IS.BatchUploadOnly )
        urlid = 'BatchUpUrl';

      let files_to_upload = [];
      if( (files.length === 1) || (urlid === 'BatchUpUrl') )
      {
        for (let file of files)
          files_to_upload.push( {url: urlid, file: file} );
        console.log( "Batch upload detected" );
      }else
      {
        // Map files to upload slots (max one per type), maintaining upload order
        const uploadOrder = ["ForegroundUpUrl", "BackgroundUpUrl", "SecondUpUrl"];
        const toUpload = Object.fromEntries(uploadOrder.map(id => [id, null]));
        
        try
        {
          for( let file of files )
          {
            const uid = urlIdFromName(file.name);
            if ( !uid || (toUpload[uid] !== null) )
              throw 'invalid file';
            toUpload[uid] = file;
          }

          files_to_upload = uploadOrder
            .filter(uid => toUpload[uid])
            .map(uid => ({url: uid, file: toUpload[uid]}));
        }catch(error)
        {
          if( window._IS.BatchUploadEnabled )
          {
            // If batch upload is enabled, we'll just trigger that.
            urlid = 'BatchUpUrl';
            for (let file of files)
              files_to_upload.push( {url: urlid, file: file} );
          }else
          {
            // TODO: specialize how this error is sent to the server so it can put in proper internationalized text
            const msg = 'showMsg-error-You can only upload a single file at a time unless the name is'
            + ' prefixed with "i-", "b-", "k-" or has back/fore/item/calib'
            + ' in the file name, and there is at most one of each spectrum type.';
            
            Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'}, msg );
            return;
          }
        }//try / catch
      }
      
      function removeUploading(){
        //console.log( 'removeUploading' );
        var upEl = document.getElementById('Uploading');
        if( upEl ) upEl.remove();

        // Clear the time-out timer, if we have one
        clearTimeout( window._IS.UploadCoverTimer );
        window._IS.UploadCoverTimer = null;
      };
      
      
      /* Lets set a 1-minute timeout so we will eventually remove the uploading cover if no traffic. */
      function setUploadTimer(){
        clearTimeout( window._IS.UploadCoverTimer );
        const coverTimer = setTimeout( function(){
          console.error("The request to upload file timed out.");
          removeUploading();
        }, 60000 );
        window._IS.UploadCoverTimer = coverTimer;
      };
      
      function errfcn(){
        removeUploading();
        alert( "Error uploading file:" + "\n\t" + xhr.statusText + "\n\tResponse code " + xhr.status + "\n\t" + xhr.responseText );
      }
      
    
      function uploadNextFile(){
        if( files_to_upload.length === 0 ){
          removeUploading();
        }else{
          const thisUpload = files_to_upload.shift();
          uploadFileToUrl( window._IS[thisUpload.url], thisUpload.file, true );
        }
      };
      
      function uploadFileToUrl( uploadURL, file, lookForPath ){
        if( !uploadURL || !file )
        {
          console.log( 'uploadFileToUrl: - invalid url or file' );
          return;
        }
          
        //console.log( 'uploadFileToUrl: will upload to ' + uploadURL );
        
        let xhr = new XMLHttpRequest();
        xhr.open("POST", uploadURL, true);
        xhr.setRequestHeader("Cache-Control", "no-cache");
        xhr.setRequestHeader("X-Requested-With", "XMLHttpRequest");
        
        
        if( lookForPath && file.name && (file.name.length > 3) ){
          let fspath = null;
          if( (typeof file.path === "string") && file.path.length > 3 ){
            fspath = file.path;
          }else{
            const fns = window._IS.dragOverFilePaths;
            if( fns && Array.isArray(fns.filenames) && fns.time && (Math.abs(fns.time - (new Date())) < 30000) ){
              const filename_uri = encodeURIComponent(file.name);
              
              for( let i = 0; i < fns.filenames.length; ++i ) {
                if( encodeURIComponent(fns.filenames[i]).endsWith(filename_uri) ){
                  fspath = fns.filenames[i];
                  //remove this filename from the array, incase there are multiple files with same leaf-name, but diff paths
                  fns.filenames.splice(i,1);
                  console.log( "Matched filenames: ", file.name, " vs ", fspath );
                  console.log( "Remaining files: ", fns.filenames );
                  if( fns.filenames.length === 0 )
                    window._IS.dragOverFilePaths = null;
                  break;
                }//if( native fn ends with browsers fn )
              }//for( loop over fns.filenames.length )
            }//if( we have 'dragOverFilePaths' object, and its valid )
          }//if( we have `file.path` available ) else if( we have 'dragOverFilePaths' avaiable )
          
          if( (typeof fspath === "string") && (fspath.length > 2) && !fspath.toLowerCase().includes('fake') ) {
            // For Electron and macOS builds, fspath gives the full path (including filename), so we
            //   will just upload the path so we can read it in the c++ to avoid lots of copying and
            //   stuff.  But we also need to fallback, JIC c++ fails for some unforeseen reason.
            //   See #FileDragUploadResource::handleRequest for the other part of this logic
            xhr.setRequestHeader("Is-File-Path", "1" );
            xhr.setRequestHeader("Content-type", "application/x-spectrum");

            // Filename may have non-ISO-8859-1 code points, so we need to URI encode it
            const filename_uri = encodeURIComponent(file.name);
            xhr.setRequestHeader("X-File-Name", filename_uri);
          
            console.log( 'Will send native file-system path, instead of actual file; path: ', fspath );
          
            xhr.onload = function(){
              removeUploading();
            
              if( (xhr.readyState === 4) && (xhr.status !== 200) ){
                // Fallback to standard http upload, since reading in c++ failed
                console.log( 'Failed to upload via native file-system path.' );
                uploadFileToUrl( uploadURL, file, false );
              }else if( xhr.readyState === 4 ) {
                // Reading in c++ succeeded
                console.log( 'Successfully uploaded native file-system path.' );
                uploadNextFile();
              }
            };
          
            xhr.onloadend = removeUploading;
        
            const body = JSON.stringify( { fullpath: fspath} );
            xhr.send(body);
            return;
          }//if( have full path on Electron or macOS build )
        }
        
        xhr.onerror = errfcn;
        xhr.onload = function(){
          //console.log( 'xhr.readyState=' + xhr.readyState + ', xhr.status=' + xhr.status );
          if( (xhr.readyState === 4) && (xhr.status !== 200) )
            errfcn();
            
          if( xhr.readyState === 4 )
            uploadNextFile();
        };//xhr.onload
        
        
        // I'm not totally sure when functions for XMLHttpRequest vs XMLHttpRequest.upload get called
        // https://developer.mozilla.org/en-US/docs/Web/API/XMLHttpRequest/Using_XMLHttpRequest
        
        let lastUpdateTime = 0;
        xhr.upload.addEventListener('progress', function(pe) {
          if(pe.lengthComputable) {
            //console.log( 'pe.total=' + pe.total + ', pe.loaded=' + pe.loaded );
            const currentTime = Date.now();
            if( (currentTime - lastUpdateTime) > 500 ){
              lastUpdateTime = currentTime;
              var prog = document.getElementById('UploadingProgress');
              if( prog ) prog.textContent = Math.trunc(100*pe.loaded/pe.total) + '%';
            }
          }
          setUploadTimer(); //reset the upload timer, since we have recieved some data
        });
        
        xhr.upload.addEventListener("load", removeUploading);
        xhr.upload.addEventListener("error", removeUploading);
        
        xhr.addEventListener('abort', removeUploading);
        xhr.upload.addEventListener("abort", removeUploading);
        
        // The loadend event is fired when a request has completed, whether successfully (after
        //  load) or unsuccessfully (after abort or error).
        //  If server does not send response, then this function will not be called until
        //   xhr.timeout/xhr.ontimeout get triggered.
        xhr.addEventListener('loadend', function() {
          //console.log( 'onloadend called: ' + ', readyState=' + xhr.readyState
          //+ ', status=' + status + ', responseTxt=' + xhr.responseText );
          removeUploading();
        } );//
        

        // This next function will be called once data is sent to server, but the server doesnt
        //  have to send the repsonce.
        xhr.upload.addEventListener('loadend', function() {
          //console.log( 'xhr.upload.loadend called' );
          removeUploading();
        });
        
        // Instead of having different urls for foreground/background/secondary/batch, we
        //  could have a single upload resource, and instead use "foreground", "background", 
        // "secondary", "multiple", "batch" as the parameter name, and then sort this out on 
        // the server.  
        // This would also let us send multiple files at the same time.
        // But we will leave this to the future.
        const formData = new FormData();
        formData.append('file', file, file.name); 

        xhr.send(formData);
      };//uploadFileToUrl
      
      uploadNextFile();
      
      // Lets indicate to the user that something is going on; even on localhost, uploads speeds are
      //  only a few MB per second, so a sizable file can take quite a while.
      target.insertAdjacentHTML('beforeend',
        '<div id="Uploading" class="Wt-dialogcover UploadProgressCover">'
        + '<div class="UpFileContents"><div class="UpCenter"><div class="UpDivTxt">'
        + '<div>Copying file into InterSpec</div>'
        + '<div id="UploadingProgress">--%</div>'
        + '</div></div></div>'
        + '</div>' );
      
      setUploadTimer();
    }catch(error)
    {
      console.log( "Error in HandleDrop: " + error );
    }
  };
  
  target.addEventListener("dragover", function(event){
    event.preventDefault();
    window._IS.IsDragging = true;
  }, false);
  
  target.addEventListener("drop", function(event){
    event.preventDefault();
    var u = document.getElementById('Uploader');
    if( u ) u.remove();
  }, false);
  
  
  target.addEventListener("dragleave", function(event){
    window._IS.IsDragging = false;
    var el = event.target;
    if( el.classList ) el.classList.remove("HasFile"); //dragenter is called before dragleave
    var updiv = (el.classList && el.classList.contains("UpDiv")) ? el : el.closest(".UpDiv");
    if( updiv && !updiv.classList.contains("HasFile") && updiv.querySelector(".HasFile") === null )

    if( window._IS.HasForeground )
      if( updiv ) updiv.classList.remove('active');

    //    if(event.relatedTarget===null)
    //      remove uploader

    //WebKit browsers dont seem to have a valid event.relatedTarget a lot of times
    clearTimeout( window._IS.DragTimer );
    var timeout = setTimeout( function(){
      if( !window._IS.IsDragging ){ var u = document.getElementById('Uploader'); if( u ) u.remove(); }
    }, 200 );
    window._IS.DragTimer = timeout;

  }, false);
  
  
  target.addEventListener("dragenter", function(event){

      // If the "batch analysis" GUI is showing, then any file we drop will be uploaded to the batch upload URL.
    if( window._IS.BlockFileDrops ){
      event.preventDefault();
      event.stopPropagation();
      return;
    }

    var hasfore = window._IS.HasForeground;

    var thisel = event.target;
    window._IS.IsDragging = true;
    var updiv = (thisel.classList && thisel.classList.contains("UpDiv")) ? thisel : thisel.closest(".UpDiv");

    if( thisel.classList ) thisel.classList.add("HasFile"); //dragenter is called before dragleave

    if( hasfore && updiv )
      updiv.classList.add("active");

    const is_batch_upload = ((event.dataTransfer  && event.dataTransfer.items
                              && (event.dataTransfer.items.length > 2)
                              && window._IS.BatchUploadEnabled)
                            || window._IS.BatchUploadOnly );

    //For windows Qt version of app, when someone drags from Outlook onto app, the file contents
    //  dropped and intercepted in the C++ WebViewWindow::dropEvent() function, instead
    //  of by the browser. The fact this is happening is indicated by the JS "DropFileContents"
    //  data; so we will treat this case a little special (no user choice of background, etc, sonce this would require work in Qt).
    var fcdivtxt = '<div class="UpFileContents"><div class="UpCenterDiv"><div class="UpDivTxt">Drop to Open</div></div></div>';

    var isDropContents = window._IS.DropFileContents;
    var el = document.getElementById('Uploader');

    if( is_batch_upload && el )
    {
      var infodiv = document.querySelector('.UpInfo');
      if( infodiv ) infodiv.innerHTML = "Will start batch tool.";

      return;
    }else if( !isDropContents && el )
    {
        var infodiv = document.querySelector('.UpInfo');
        if( hasfore || (updiv && updiv.classList.contains("UpForeground")) )
        {
            if( event.target.id==="Uploader" )
            { if( infodiv ) infodiv.innerHTML = "Defaulting to upload as Foreground"; }
            else
            {
              var updivTxt = updiv ? updiv.querySelector(".UpDivTxt") : null;
              if( infodiv ) infodiv.innerHTML = "Will open as " + (updivTxt ? updivTxt.innerHTML : "");
            }

            /*Could add some text such as "Multiple files may be uploaded if the names contain information such as i-/item/fore/b-/back/bkg/k-/known/cal/check" */
        }else
        {
            if( infodiv ) infodiv.innerHTML = "Will open as the foreground spectrum since there is currently not one loaded";
        }

        return;
    }else if( el )
    {
      if( !el.querySelector(':scope > .UpFileContents') )
      {
        while( el.firstChild ) el.firstChild.remove();
        el.insertAdjacentHTML('beforeend', fcdivtxt);
      }

      return;
    }

    //The list of files (or even how many of them) does not appear to be available here for Safari.
    //  But were it is available, we could customize it for CALp files, Det. Eff. files, etc .
    //const files = event.dataTransfer && event.dataTransfer.files ? event.dataTransfer.files : null;
    //const isCALp = (files && files.length === 1 && files[0].includes(".CALp"));


    el = document.createElement('div');
    el.id = 'Uploader';
    el.className = 'FileUploadCover';
    target.appendChild(el);


    if( is_batch_upload )
    {
      el.insertAdjacentHTML('beforeend', '<div class="UpDiv UpBatch active"><div class="UpCenterDiv"><div class="UpDivTxt">Batch Upload</div></div></div>');
      var batch = el.lastElementChild;

      batch.insertAdjacentHTML('beforeend', '<div class="UpInfoWrapper"><div class="UpInfo"></div></div>');

      el.addEventListener("drop", function(event){
        var u = document.getElementById('Uploader'); if( u ) u.remove();
        uploadFcn(event,'BatchUpUrl')
      }, false);

      return;
    }

    el.addEventListener("drop", function(event){
      var u = document.getElementById('Uploader'); if( u ) u.remove();
      uploadFcn(event,'ForegroundUpUrl')
    }, false);

    if( isDropContents )
    {
      el.insertAdjacentHTML('beforeend', fcdivtxt);
      return;
    }

    el.insertAdjacentHTML('beforeend', '<div class="UpDiv UpForeground"><div class="UpCenterDiv"><div class="UpDivTxt">Foreground</div></div></div>');
    var up = el.lastElementChild;

    el.insertAdjacentHTML('beforeend', '<div class="UpDiv UpBackground"><div class="UpCenterDiv"><div class="UpDivTxt">Background</div></div></div>');
    var back = el.lastElementChild;

    el.insertAdjacentHTML('beforeend', '<div class="UpDiv Up2ndForeground"><div class="UpCenterDiv"><div class="UpDivTxt">2<sup>nd</sup> Foreground</div></div></div>');
    var sec = el.lastElementChild;

    var titletxt = "Drop file on how you want to view the spectrum";
    if( !hasfore )
      titletxt = "Drop file to view the spectrum";

    back.insertAdjacentHTML('beforeend', '<div class="UpInstrInfoWrapper"><div class="UpInstr">' + titletxt + '</div></div>');

    back.insertAdjacentHTML('beforeend', '<div class="UpInfoWrapper"><div class="UpInfo"></div></div>');


    if( hasfore )
    {
      up.addEventListener("drop", function(event){
        var u = document.getElementById('Uploader'); if( u ) u.remove();
        uploadFcn(event,'ForegroundUpUrl')
      }, false);
      back.addEventListener("drop", function(event){
        var u = document.getElementById('Uploader'); if( u ) u.remove();
          uploadFcn(event,'BackgroundUpUrl')
      }, false);
      sec.addEventListener("drop", function(event){
        var u = document.getElementById('Uploader'); if( u ) u.remove();
        uploadFcn(event,'SecondUpUrl')
      }, false);
    }else
    {
      up.classList.add("active");
      back.querySelector(".UpDivTxt").classList.add("disabled");
      sec.querySelector(".UpDivTxt").classList.add("disabled");
      back.querySelector(".UpCenterDiv").classList.add("disabled");
      sec.querySelector(".UpCenterDiv").classList.add("disabled");
    }
  }, false );
}
);

WT_DECLARE_WT_MEMBER
(DoOrientationChange, Wt::JavaScriptFunction, "DoOrientationChange",
function()
{
  // TODO: `screen.orientation.type` can take on values of
  //       'portrait-primary', 'portrait-secondary',
  //       'landscape-primary' (LandscapeLeft), 'landscape-secondary' (LandscapeRight)
  //       Should maybe use this instead - and also adjust to have this work for Android as well
  console.log( "DoOrientationChange: ", screen.orientation );
  
  let angle = window.orientation; //Apple devices, not sure if Android has this, I dont think so
  if( (typeof angle !== "number") && screen && screen.orientation )
    angle = screen.orientation.angle;
  
  if( (typeof angle === "number")
    //&& CSS
    //&& (CSS.supports('padding-bottom: env(safe-area-inset-left)')
    //    || CSS.supports('padding-bottom: constant(safe-area-inset-left)'))
    ) {
    
    switch( angle ) {
      case 0:
      {
        const dr = document.querySelector(".Wt-domRoot");
        dr.classList.remove("LandscapeRight", "LandscapeLeft");
        dr.classList.add("Portrait");
      }
      break;

      case 90:
      {
        const dr = document.querySelector(".Wt-domRoot");
        dr.classList.remove("LandscapeRight", "Portrait");
        dr.classList.add("LandscapeLeft");
      }
      break;

      case -90:
      {
        const dr = document.querySelector(".Wt-domRoot");
        dr.classList.remove("LandscapeLeft", "Portrait");
        dr.classList.add("LandscapeRight");
      }
      break;
      
      default:
        console.error( "Unknown orientation angle: ", angle );
      break;
    }
    window.dispatchEvent(new Event('resize')); //see also Wt.TriggerResizeEvent()
  } else {
    console.warn( "Unknown way to read orientation; screen:", screen, ", window.orientation:", window.orientation );
  }
}
);


//Requires {macOS 10.14, iOS 13.0, iPadOS 13.0, Windows 10} and {Chrome 76, Safari 12.1, Firefox 67}

WT_DECLARE_WT_MEMBER
(DetectOsColorThemeJs, Wt::JavaScriptFunction, "DetectOsColorThemeJs",
function(id)
{
  try
  {
    const darkQ = window.matchMedia('(prefers-color-scheme: dark)');
    const lightQ = window.matchMedia('(prefers-color-scheme: light)');
    const noPrefQ = window.matchMedia('(prefers-color-scheme: no-preference)');
    
    if( darkQ.matches )
      Wt.emit( id, {name: 'OsColorThemeChange' }, 'dark' );
    else if( lightQ.matches )
      Wt.emit( id, {name: 'OsColorThemeChange' }, 'light' );
    else if( noPrefQ.matches )
      Wt.emit( id, {name: 'OsColorThemeChange' }, 'no-preference' );
    else
      throw 'no-support';
  }catch(error)
  {
    Wt.emit( id, {name: 'OsColorThemeChange' }, 'no-support' );
  }
});

WT_DECLARE_WT_MEMBER
(SetupOsColorThemeChangeJs, Wt::JavaScriptFunction, "SetupOsColorThemeChangeJs",
function(id)
{
  try
  {
    const darkQ = window.matchMedia('(prefers-color-scheme: dark)');
    const lightQ = window.matchMedia('(prefers-color-scheme: light)');
    const noPrefQ = window.matchMedia('(prefers-color-scheme: no-preference)');
    
    if( !darkQ.matches && !lightQ.matches && !noPrefQ.matches )
      console.warn( 'No matches for color scheme found.' );
    
    darkQ.addEventListener( 'change', function(e){
      if( e && e.matches )
        Wt.emit( id, {name: 'OsColorThemeChange' }, 'dark' );
    } );

    lightQ.addEventListener( 'change', function(e){
      if( e && e.matches )
        Wt.emit( id, {name: 'OsColorThemeChange' }, 'light' );
    } );

    noPrefQ.addEventListener( 'change', function(e){
      if( e && e.matches )
        Wt.emit( id, {name: 'OsColorThemeChange' }, 'no-preference' );
    } );
    
  }catch(error)
  {
    console.warn( 'Color theme not supported:', error );
  }
}
);


WT_DECLARE_WT_MEMBER
(InitFlexResizer, Wt::JavaScriptFunction, "InitFlexResizer",
function(resizerid,elid) {
  
  // TODO:
  //  - [ ] Check height to be at least 10px?  Maybe make it be an option
  //  - [ ] The mouse can easily go out of the resizer, and if it goes into the spectrum chart, we
  //        will stop getting here, so we either need to up the order of calling events through
  //        some jQuery version specific hack (see https://stackoverflow.com/questions/2360655/jquery-event-handlers-always-execute-in-order-they-were-bound-any-way-around-t),
  //        or get the spectrum chart to stop canceling the events, which may take some work.
  //  - [ ] Allow passing option if element being resized is above or below resizer.
  //  - [ ] Could probably be made much more flexible, but good enough for now.
  //  - [ ] listen for 'esc' key, and if so, set back to original height and stop drag
  //  - [ ] IF you drag high-enough, the bottom chart will go off the bottom - this needs to be prevented
  //  - [ ] Unrelated to this, but got "dimensions of D3TimeChart div element are not set." error after resizing chart then switching to non-time chart more.
  let startPosAndSize;
  const allowTouch = (window.Touch || navigator.maxTouchPoints);
      
  function mousePos(e) {
    if( typeof e.clientX === "number" )
      return { x: e.clientX, y: e.clientY };
    else if( e.touches )
      return { x: e.touches[0].clientX, y: e.touches[0].clientY };
    return null;
  }
      
  function cancelEvent(e) {
    e.stopPropagation();
    e.preventDefault();
  }
      
  var flexController = null;

  function handleMouseDown(e) {
    const el = document.getElementById(elid);

    startPosAndSize = mousePos(e);
    startPosAndSize.height = el.offsetHeight;
    el.style.flexBasis = startPosAndSize.height + 'px';
    el.style.flexShrink = el.style.flexGrow = 0;

    // Abort any previous drag-session listeners
    if( flexController ) flexController.abort();
    flexController = new AbortController();
    var sig = flexController.signal;

    document.addEventListener('mousemove', handleMouseMove, {signal: sig});
    document.addEventListener('mouseup', handleMouseUp, {signal: sig});
    if( allowTouch ) {
      document.addEventListener('touchmove', handleMouseMove, {signal: sig});
      document.addEventListener('touchend', handleMouseUp, {signal: sig});
    }
    document.addEventListener('selectstart', cancelEvent, {signal: sig}); // disable selection
  }

  function handleMouseMove(e) {
    const pos = mousePos(e);
    const el = document.getElementById(elid);

    // See if our hight has already gone too large for the parent, and if so, stop making it go so large.
    //const parent = el.parentNode;
    //const parentHeight = parent.height(); //This is useable inner height, after accounting for padding and such.
    //var children = parent.children;
    //for( let i = 0; i < children.length; i++ ) {
    //  const child = children[i];
    //  const minHeight = window.getComputedStyle(child,null).getPropertyValue("min-height");
    //}

    const height = (startPosAndSize.height + startPosAndSize.y - pos.y);

    el.style.flexBasis = height + 'px'
  }

  function handleMouseUp(e) {
    // We could look at element above/bellow resizer, and got back to flex basline of zero but give the grows as the current ratio

    if( flexController ){ flexController.abort(); flexController = null; }

    return false;
  }

  // Remove any previous bindings, JIC this is a duplicate call to initialize the events
  //  (it doesnt look like this happens, but current code to show/hide tool-tabs is pretty
  //   fragile, so the off is really just to make sure).
  var resizerEl = document.getElementById(resizerid);
  if( resizerEl )
  {
    // Store AbortController on the element for cleanup
    if( resizerEl._flexInitCtrl ) resizerEl._flexInitCtrl.abort();
    resizerEl._flexInitCtrl = new AbortController();
    resizerEl.addEventListener('mousedown', handleMouseDown, {signal: resizerEl._flexInitCtrl.signal});
    resizerEl.addEventListener('touchstart', handleMouseDown, {signal: resizerEl._flexInitCtrl.signal});
  }
}
);
