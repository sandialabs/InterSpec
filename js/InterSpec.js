
WT_DECLARE_WT_MEMBER
(FileUploadFcn, Wt::JavaScriptFunction, "FileUploadFcn",
function()
{
  var target = $('.Wt-domRoot').get(0);
  
  var urlIdFromName = function(name)
  {
    var str = name.toLowerCase();
    if( str.indexOf('i-')===0 || str.indexOf('fore')>-1 || str.indexOf('item')>-1
        || str.indexOf('ipc')>-1 || str.indexOf('unk')>-1 || str.indexOf('inter')>-1  )
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
      var files = evt.dataTransfer.files;

      if( !files || !files.length )
      {
        console.log( evt );
        console.log( "uploadFcn: no files for some reason." );
        return;
      }
      
      
      //Go through and make sure we have either one file total, or a max of one
      //  for each foreground/background/second.
      var validin = true, uid, filelen = files.length, x, i;
      var uploadOrder = ["ForegroundUpUrl", "BackgroundUpUrl", "SecondUpUrl"];
      var toUpload = {"ForegroundUpUrl":null, "BackgroundUpUrl":null, "SecondUpUrl":null};
      
      for(i = 0; i < filelen; ++i)
      {
        uid = filelen>1 ? urlIdFromName(files[i].name) : urlid;
        if( filelen>1 && (uid===null || (toUpload[uid]!==null)) )
          validin = false;
        else
          toUpload[uid] = files[i];
      }
      
      if( !validin )
      {
        alert( 'You can only upload a single file at a time unless the name is'
               + ' prefixed with "i-", "b-", "k-" or has back/fore/item/calib'
               + ' in the file name, and there is at most one of each spectrum type.' );
        return;
      }
      
      // We want to upload the foreground, then background, then secondary, waiting for each upload
      //  to finish before starting the next.
      var uploadIndex = -1;
      
      function uploadNextFile(){
        for( ++uploadIndex; (uploadIndex < uploadOrder.length) && (toUpload[uploadOrder[uploadIndex]] === null); ++uploadIndex )
        {
        }
        
        if( uploadIndex >= uploadOrder.length )
          return;
        
        var thisUid = uploadOrder[uploadIndex];
        uploadFileToUrl( $(target).data(thisUid), toUpload[thisUid] );
      };
      
      function uploadFileToUrl( uploadURL, file ){
        if( !uploadURL || !file )
        {
          console.log( 'uploadFileToUrl: - invalid url or file' );
          return;
        }
          
        console.log( 'uploadFileToUrl: will upload to ' + uploadURL );
        
        var xhr;
        if (window.XMLHttpRequest)
          xhr = new XMLHttpRequest();
        else if (window.ActiveXObject)
          xhr = new ActiveXObject('MSXML2.XMLHTTP.3.0');
        else
          return;
      
        xhr.open("POST", uploadURL, true);
        xhr.setRequestHeader("Content-type", "application/x-spectrum");
        xhr.setRequestHeader("Cache-Control", "no-cache");
        xhr.setRequestHeader("X-Requested-With", "XMLHttpRequest");
        xhr.setRequestHeader("X-File-Name", file.name);
        
        var errfcn = function () {
          alert( "Error uploading file:" + "\n\t" + xhr.statusText + "\n\tResponse code " + xhr.status + "\n\t" + xhr.responseText );
        };
        
        xhr.onerror = errfcn;
        xhr.onload = function(){
          //console.log( 'xhr.readyState=' + xhr.readyState + ', xhr.status=' + xhr.status );
          if( (xhr.readyState === 4) && (xhr.status !== 200) )
            errfcn();
            
          if( xhr.readyState === 4 )
            uploadNextFile();
        };//xhr.onload
        
        //xhr.timeout = 5000;
        //xhr.ontimeout = function () {
        //  console.error("The request to upload file timed out.");
        //};
          
        xhr.send(file);
      };//uploadFileToUrl
      
      uploadNextFile();
    }catch(error)
    {
      console.log( "Error in HandleDrop: " + error );
    }
  };
  
  target.addEventListener("dragover", function(event){
    event.preventDefault();
    $(target).data('IsDragging',true);
  }, false);
  
  target.addEventListener("drop", function(event){
    event.preventDefault();
    $('#Uploader').remove();
  }, false);
  
  
  target.addEventListener("dragleave", function(event){
    $(target).data('IsDragging',false);
    var el=$(event.target);
    el.removeClass("HasFile"); //dragenter is called before dragleave
    var updiv = el.hasClass("UpDiv") ? el : el.parents(".UpDiv");
    if( !updiv.hasClass("HasFile") && updiv.find(".HasFile").length===0 )
    
    if( $('.Wt-domRoot').data("HasForeground") )
      updiv.removeClass('active');
    
    //    if(event.relatedTarget===null)
    //      $('#Uploader').remove();
    
    //WebKit browsers dont seem to have a valid event.relatedTarget a lot of times
    clearTimeout( $(target).data('DragTimer') );
    var timeout = setTimeout( function(){
      if( !$(target).data('IsDragging') ) $('#Uploader').remove();
    }, 200 );
    $(target).data('DragTimer',timeout);
    
  }, false);
  
  
  target.addEventListener("dragenter", function(event){
    var hasfore = $('.Wt-domRoot').data("HasForeground");
    
    var thisel=$(event.target);
    $(target).data('IsDragging',true);
    var updiv = thisel.hasClass("UpDiv") ? thisel : thisel.parents(".UpDiv");
    
    thisel.addClass("HasFile"); //dragenter is called before dragleave
    
    if( hasfore )
      updiv.addClass("active");
    
    //For windows Qt version of app, when someone drags from Outlook onto app, the file contents
    //  dropped and intercepted in the C++ WebViewWindow::dropEvent() function, instead
    //  of by the browser. The fact this is happening is indicated by the JS "DropFileContents"
    //  data; so we will treat this case a little special (no user choice of background, etc, sonce this would require work in Qt).
    var fcdivtxt = '<div class=\"UpFileContents\"><div class=\"UpCenterDiv\"><div class=\"UpDivTxt\">Drop to Open</div></div></div>';

    var isDropContents = $('.Wt-domRoot').data("DropFileContents");
    var el = $('#Uploader');
    
    if( !isDropContents && (el.length > 0) )
    {
        var infodiv = $('.UpInfo');
        if( hasfore || updiv.hasClass("UpForeground") )
        {
            if( event.target.id==="Uploader" )
                infodiv.html("Defaulting to upload as Foreground");
            else
                infodiv.html("Will open as " + updiv.find(".UpDivTxt").html() );
          
            /*Could add some text such as "Multiple files may be uploaded if the names contina information such as i-/item/fore/b-/back/bkg/k-/known/cal/check" */
        }else
        {
            infodiv.html("Will open as the foreground spectrum since there is currenly not one loaded");
        }
      
        return;
    }else if( el.length > 0 )
    {
      if( el.children('.UpFileContents').length < 1 )
      {
        el.children().remove();
        $(fcdivtxt).appendTo(el);
      }

      return;
    }

    //The list of files (or even how many of them) does not appear to be available here... :(
    
    el = $('<div id=\"Uploader\" class=\"FileUploadCover\"></div>');
    el.appendTo(target);
    el.get(0).addEventListener("drop", function(event){
      $('#Uploader').remove();
      uploadFcn(event,'ForegroundUpUrl')
    }, false);
    

    if( isDropContents )
    {
      $(fcdivtxt).appendTo(el);
      return;
    }
    
    var up = $('<div class=\"UpDiv UpForeground\"><div class=\"UpCenterDiv\"><div class=\"UpDivTxt\">Foreground</div></div></div>');
    up.appendTo(el);
      
    var back = $('<div class=\"UpDiv UpBackground\"><div class=\"UpCenterDiv\"><div class=\"UpDivTxt\">Background</div></div></div>');
    back.appendTo(el);
    
    var sec = $('<div class=\"UpDiv Up2ndForeground\"><div class=\"UpCenterDiv\"><div class=\"UpDivTxt\">2<sup>nd</sup> Foreground</div></div></div>');
    sec.appendTo(el);
    
    var titletxt = "Drop file on how you want to view the spectrum";
    if( !hasfore )
      titletxt = "Drop file to view the spectrum";
    
    var title = $('<div class=\"UpInstrInfoWrapper\"><div class=\"UpInstr\">' + titletxt + '</div></div>');
    title.appendTo(back);
    
    var info = $('<div class=\"UpInfoWrapper\"><div class=\"UpInfo\"></div></div>');
    info.appendTo(back);
      
      
    if( hasfore )
    {
      up.get(0).addEventListener("drop", function(event){
        $('#Uploader').remove();
        uploadFcn(event,'ForegroundUpUrl')
      }, false);
      back.get(0).addEventListener("drop", function(event){
        $('#Uploader').remove();
          uploadFcn(event,'BackgroundUpUrl')
      }, false);
      sec.get(0).addEventListener("drop", function(event){
        $('#Uploader').remove();
        uploadFcn(event,'SecondUpUrl')
      }, false);
    }else
    {
      up.addClass("active");
      back.find(".UpDivTxt").addClass("disabled");
      sec.find(".UpDivTxt").addClass("disabled");
      back.find(".UpCenterDiv").addClass("disabled");
      sec.find(".UpCenterDiv").addClass("disabled");
    }
  }, false );
}
);

WT_DECLARE_WT_MEMBER
(DoOrientationChange, Wt::JavaScriptFunction, "DoOrientationChange",
function()
{
  if( CSS && (CSS.supports('padding-bottom: env(safe-area-inset-left)')
     || CSS.supports('padding-bottom: constant(safe-area-inset-left)')) ) {
    
    switch( window.orientation ) {
      case 0:
        $(".Wt-domRoot").removeClass("LandscapeRight LandscapeLeft");
      break;
      
      case 90:
        $(".Wt-domRoot").removeClass("LandscapeRight").addClass("LandscapeLeft");
      break;
      
      case -90:
        $(".Wt-domRoot").removeClass("LandscapeLeft").addClass("LandscapeRight");
      break;
    }
    
    window.dispatchEvent(new Event('resize'));
  }
}
);


//Requires {macOS 10.14, iOS 13.0, iPadOS 13.0, Windows 10} and {Chrome 76, Safari 12.1, Firefox 67}
WT_DECLARE_WT_MEMBER
(SetupOsColorThemeChangeJs, Wt::JavaScriptFunction, "SetupOsColorThemeChangeJs",
function(id)
{
  try
  {
    var darkQ = window.matchMedia('(prefers-color-scheme: dark)');
    var lightQ = window.matchMedia('(prefers-color-scheme: light)');
    var noPrefQ = window.matchMedia('(prefers-color-scheme: no-preference)');
    var isSupported = darkQ.matches || lightQ.matches || noPrefQ.matches;
    
    if( !isSupported )
      throw 'no-support';
    
    if( darkQ.matches )
      Wt.emit( id, {name: 'OsColorThemeChange' }, 'dark' );
    
    if( lightQ.matches )
      Wt.emit( id, {name: 'OsColorThemeChange' }, 'light' );
    
    if( noPrefQ.matches )
      Wt.emit( id, {name: 'OsColorThemeChange' }, 'no-preference' );
    
    darkQ.addListener( function(e){
      if( e && e.matches )
        Wt.emit( id, {name: 'OsColorThemeChange' }, 'dark' );
    } );
    
    lightQ.addListener( function(e){
      if( e && e.matches )
        Wt.emit( id, {name: 'OsColorThemeChange' }, 'light' );
    } );
    
    noPrefQ.addListener( function(e){
      if( e && e.matches )
        Wt.emit( id, {name: 'OsColorThemeChange' }, 'no-preference' );
    } );
    
  }catch(error)
  {
    Wt.emit( id, {name: 'OsColorThemeChange' }, 'no-support' );
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
    else if( e.originalEvent.touches )
      return { x: e.originalEvent.touches[0].clientX, y: e.originalEvent.touches[0].clientY };
    return null;
  }
      
  function cancelEvent(e) {
    e.stopPropagation();
    e.preventDefault();
  }
      
  function handleMouseDown(e) {
    const el = document.getElementById(elid);
    const resizer = document.getElementById(resizerid);
    
    startPosAndSize = mousePos(e);
    startPosAndSize.height = el.offsetHeight;
    el.style.flexBasis = startPosAndSize.height + 'px';
    el.style.flexShrink = el.style.flexGrow = 0;
        
    $(document).bind('mousemove.FlexResize', handleMouseMove);
    $(document).bind('mouseup.FlexResize', handleMouseUp);
    if( allowTouch ) {
      $(document).bind('touchmove.FlexResize', handleMouseMove);
      $(document).bind('touchend.FlexResize', handleMouseUp);
    }
    $(document).bind('selectstart.FlexResize', cancelEvent); // disable selection
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
    
    $(document).unbind('mousemove.FlexResize', handleMouseMove);
    $(document).unbind('mouseup.FlexResize', handleMouseUp);
        
    if( allowTouch ) {
      $(document).unbind('touchmove.FlexResize', handleMouseMove);
      $(document).unbind('touchend.FlexResize', handleMouseUp);
    }
    $(document).unbind('selectstart.FlexResize', cancelEvent);
                      
    return false;
  }

  $('#' + resizerid).bind('mousedown.FlexResize touchstart.FlexResize', handleMouseDown);
}
);
