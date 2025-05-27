/** Called when files  */
function onDragOverStart(event){
  event.preventDefault(); // to allow dropping
  event.stopPropagation(); // stop event from bubbling up
  console.log( "onDragOverStart" );
  
  try {
    const domRoot = document.querySelector('.Wt-domRoot');
    console.log( "onDragOverStart domRoot.dataset.batchUploadIds:", domRoot.dataset.batchUploadIds );
    const idsList = JSON.parse(domRoot.dataset.batchUploadIds);
    console.log( "onDragOverStart idsList:", idsList );
    for( const widget_id of idsList ){
      const el = document.getElementById(widget_id);
      if( el )
        el.classList.add('DomIsDrugOver');
     else
        console.log( "onDragOverStart: no element for id", widget_id );
    }
  }catch(err){
    console.log( "onDragOverStart error:", err );
  }
}//function onDragOverStart(event)


function onDragOverEnd(event){
  try {
    const domRoot = document.querySelector('.Wt-domRoot');
    const idsList = JSON.parse(domRoot.dataset.batchUploadIds);
    for( const widget_id of idsList ){
      const el = document.getElementById(widget_id);
      if( el )
        el.classList.remove('DomIsDrugOver');
      else
       console.log( "onDragOverEnd: no element for id", widget_id );
    }
  }catch(err){
    console.log( "onDragOverStart error:", err );
  }
}//function onDragOverEnd(event)


/** Setups DOM to look for files being started to be drug over.
 * @param target_id_list A list/array of ids (strings) of element that should
 *        have the `DomIsDrugOver` class added when this happens.
 */
function setupOnDragEnterDom( target_id_list ){
  console.log( "setupOnDragEnterDom:", target_id_list );
  const domRoot = document.querySelector('.Wt-domRoot');
  domRoot.addEventListener("dragenter", onDragOverStart, false );
  domRoot.addEventListener("dragover", onDragOverStart, false );
  domRoot.addEventListener("dragleave", onDragOverEnd, false);
  
  let oldIdsList = [];
  try{
    oldIdsList = JSON.parse(domRoot.dataset.batchUploadIds);
  }catch(err){
  }

  const idsList = [...oldIdsList, ...target_id_list];
  domRoot.dataset.batchUploadIds = JSON.stringify(idsList);
  console.log( "setupOnDragEnterDom: new idsList:", domRoot.dataset.batchUploadIds );
}//function setupOnDragEnterDom( target_id_list )


function removeOnDragEnterDom( target_id_list ){
  try {
    const domRoot = document.querySelector('.Wt-domRoot');
    domRoot.removeEventListener("dragenter", onDragOverStart, false );

    const batchUploadIdStrs = domRoot.dataset.batchUploadIds;
    let idsList = JSON.parse(batchUploadIdStrs);
    for( const widget_id of target_id_list ){
      const index = idsList.indexOf(widget_id);
      if (index > -1) { // Check if the string was found
        arr.splice(index, 1); // Remove 1 element at the found index
      }
    }
    console.log( "After removing ", target_id_list, " idsList=", idsList );
    domRoot.dataset.batchUploadIds = JSON.stringify(idsList);
  }catch(err){
    console.log( "removeOnDragEnterDom error:", err );
  }
}//function removeOnDragEnterDom( target_id_list )



function BatchInputDropUploadSetup( target )
{
  const uploadFcn = function(evt){
    try
    {
      const uploadURL = $('.Wt-domRoot').data('BatchUpUrl');
      let files_to_upload = [];
      
      for( const file of evt.dataTransfer.files )
          files_to_upload.push( file );

      function removeUploading(){
        target.classList.remove( "Uploading" );
        onDragOverEnd(null);
        // Clear the time-out timer, if we have one
        clearTimeout( $('.Wt-domRoot').data('BatchUploadTimer') );
        $('.Wt-domRoot').data( 'BatchUploadTimer', null );
      };


      /* Lets set a 1-minute timeout so we will eventually remove the uploading cover if no traffic. */
      function setUploadTimer(){
        clearTimeout( $('.Wt-domRoot').data('BatchUploadTimer') );
        const coverTimer = setTimeout( function(){
          console.error("The request to upload file timed out.");
          removeUploading();
        }, 60000 );
        $('.Wt-domRoot').data( 'BatchUploadTimer', coverTimer );
      };

      function errfcn(){
        removeUploading();
        alert( "Error uploading file:" + "\n\t" + xhr.statusText + "\n\tResponse code " + xhr.status + "\n\t" + xhr.responseText );
      }


      function uploadNextFile(){
        if( files_to_upload.length === 0 ){
          removeUploading();
        }else{
          uploadFileToUrl( files_to_upload.shift() );
        }
      };

      function uploadFileToUrl( file ){
        if( !uploadURL || !file )
        {
          console.log( 'uploadFileToUrl: - invalid url or file' );
          return;
        }

        let xhr = new XMLHttpRequest();
        xhr.open("POST", uploadURL, true);
        xhr.setRequestHeader("Cache-Control", "no-cache");
        xhr.setRequestHeader("X-Requested-With", "XMLHttpRequest");


        if( file.name && (file.name.length > 3) ){
          let fspath = null;
          if( (typeof file.path === "string") && file.path.length > 3 ){
            fspath = file.path;
          }else{
            const fns = $(document).data('dragOverFilePaths');
            if( fns && Array.isArray(fns.filenames) && fns.time && (Math.abs(fns.time - (new Date())) < 30000) ){
              for( let i = 0; i < fns.filenames.length; ++i ) {
                if( encodeURIComponent(fns.filenames[i]).endsWith(filename_uri) ){
                  fspath = fns.filenames[i];
                  //remove this filename from the array, incase there are multiple files with same leaf-name, but diff paths
                  fns.filenames = fns.filenames.splice(i,i);
                  console.log( "Matched filenames: ", file.name, " vs ", fspath );
                  console.log( "Remaining files: ", fns.filenames );
                  if( fns.filenames.length === 0 )
                    $(document).data('dragOverFilePaths', null);
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
        }//if( file.name && (file.name.length > 3) )

        xhr.onerror = errfcn;
        xhr.onload = function(){
          //console.log( 'xhr.readyState=' + xhr.readyState + ', xhr.status=' + xhr.status );
          if( (xhr.readyState === 4) && (xhr.status !== 200) )
            errfcn();

          if( xhr.readyState === 4 )
            uploadNextFile();
        };//xhr.onload

        let lastUpdateTime = 0;
        xhr.upload.addEventListener('progress', function(pe) {
          if(pe.lengthComputable) {
            //console.log( 'pe.total=' + pe.total + ', pe.loaded=' + pe.loaded );
            const currentTime = Date.now();
            if( (currentTime - lastUpdateTime) > 500 ){
              lastUpdateTime = currentTime;
              $('#UploadingProgress').text( Math.trunc(100*pe.loaded/pe.total) + '%' );
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

        const formData = new FormData();
        formData.append('file', file, file.name);

        xhr.send(formData);
      };//uploadFileToUrl

      uploadNextFile();

      target.classList.add("Uploading");

      setUploadTimer();
    }catch(error)
    {
      console.log( "Error in HandleDrop: " + error );
    }
  };

  target.addEventListener("drop", function(event){
    target.classList.remove("DragedOver");
    event.preventDefault(); // 
    event.stopPropagation();
    uploadFcn(event)
  }, false);

  target.addEventListener("dragleave", function(event){
   target.classList.remove("DragedOver");
  }, false);

  target.addEventListener("dragenter", function(event){
    event.preventDefault();
    target.classList.add("DragedOver");
  }, false);

  target.addEventListener("dragover", function(event){
    event.preventDefault();
    target.classList.add("DragedOver");
  }, false);
}//function setupOnDragEnterDom( target_id_list )