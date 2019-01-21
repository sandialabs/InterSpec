
WT_DECLARE_WT_MEMBER
(FileUploadFcn, Wt::JavaScriptFunction, "FileUploadFcn",
function()
{
  var target = $('.Wt-domRoot').get(0);
  
  var urlIdFromName = function(name)
  {
    var str = name.toLowerCase();
    if( str.indexOf('i-')===0 || str.indexOf('fore')>-1 || str.indexOf('item')>-1 || str.indexOf('ipc')>-1 || str.indexOf('unk')>-1 || str.indexOf('kwn')>-1 )
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
      //  for each foreground/background/second.  Also, make sure we upload
      //  foreground first
      var validin = true, uid, filelen = files.length, x, i;
      var dict = {"ForegroundUpUrl":null, "BackgroundUpUrl":null, "SecondUpUrl":null};
      
      for(i = 0; i < filelen; ++i)
      {
        uid = filelen>1 ? urlIdFromName(files[i].name) : urlid;
        if( filelen>1 && (uid===null || (dict[uid]!==null)) )
          validin = false;
        else
          dict[uid] = i;
      }
      
      if( !validin )
      {
        alert( 'You can only upload a single file at a time unless the name is'
               + ' prefixed with "i-", "b-", "k-" or has back/fore/item/calib'
               + ' in the file name, and there is at most one of each spectrum type.' );
        return;
      }

      for( uid in dict )
      {
        if( dict[uid]===null )
          continue;
        var xhr, uploadURL = $(target).data(uid), file = files[dict[uid]];
        if( !uploadURL )
          return;
        if (window.XMLHttpRequest)
          xhr = new XMLHttpRequest();
        else if (window.ActiveXObject)
          xhr = new ActiveXObject('MSXML2.XMLHTTP.3.0');
        else
          return;
      
        xhr.open("POST", uploadURL, false);
        xhr.setRequestHeader("Content-type", "application/x-spectrum");
        xhr.setRequestHeader("Cache-Control", "no-cache");
        xhr.setRequestHeader("X-Requested-With", "XMLHttpRequest");
        xhr.setRequestHeader("X-File-Name", file.name);
        xhr.send(file);
      
        if( xhr.status !== 200 && xhr.responseText )
          alert( xhr.responseText + "\nResponse code " + xhr.status );
      }
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
