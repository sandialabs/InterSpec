/* Drop (or browse for) detector-efficiency files and folders, for the Detector Efficiency dialog's
   Import tab (DrfImportWidget).

   Every file is uploaded to `url` with its path within the drop ("Lab/HPGe-A/Detector.dat") in the
   X-Drf-Path header, and all the uploads of one drop share an X-Drf-Batch id; a final request with
   X-Drf-Done tells the server the batch is complete.  Folders are walked with webkitGetAsEntry (drag
   and drop) or <input webkitdirectory> ("Browse folder...").

   In the macOS app, the native code tells the page the file-system paths of whatever is dragged
   over it (window._IS.dragOverFilePaths, see target/macos/AppDelegate.mm); for a single dropped
   folder, its path is sent too (X-Drf-Native-Folder), so the folder can be kept where it is.
*/
function DrfImportDropSetup( target, url, filesBtn, folderBtn )
{
  if( !target )
    return;

  const MAX_FILES = 5000;                //DrfBatchUploadResource::sm_maxBatchFiles
  const MAX_FILE_BYTES = 32*1024*1024;   //DrfImport::sm_maxFileBytes
  const MAX_TOTAL_BYTES = 256*1024*1024; //DrfBatchUploadResource::sm_maxBatchBytes
  const NUM_PARALLEL = 4;
  let num_active = 0;  //batches uploading, for the 'Uploading' style
  // Of a drop of many files, only those that can describe a detector are sent.
  const DRF_FILE = /\.(dat|csv|tsv|txt|par|ecc|outx|detx|xml)$/i;

  const make_input = function( folder ){
    const input = document.createElement( 'input' );
    input.type = 'file';
    input.multiple = true;
    if( folder )
      input.webkitdirectory = true;
    input.style.display = 'none';
    target.appendChild( input );
    input.addEventListener( 'click', function(e){ e.stopPropagation(); } );
    input.addEventListener( 'change', function(){
      browsed( Array.from( input.files ).map( function(f){
        return { file: f, path: (f.webkitRelativePath || f.name) };
      } ) );
      input.value = '';
    } );
    return input;
  };

  const files_input = make_input( false );
  const folder_input = make_input( true );

  if( filesBtn )
    filesBtn.addEventListener( 'click', function(e){ e.stopPropagation(); e.preventDefault(); files_input.click(); } );
  if( folderBtn )
    folderBtn.addEventListener( 'click', function(e){ e.stopPropagation(); e.preventDefault(); folder_input.click(); } );

  // Clicking the area itself (not one of its buttons) browses for files.
  target.addEventListener( 'click', function(e){
    if( e.target && e.target.closest && e.target.closest( 'button' ) )
      return;
    files_input.click();
  } );

  target.addEventListener( 'dragenter', function(e){ e.preventDefault(); e.stopPropagation(); target.classList.add( 'DragOver' ); } );
  target.addEventListener( 'dragover', function(e){ e.preventDefault(); e.stopPropagation(); target.classList.add( 'DragOver' ); } );
  target.addEventListener( 'dragleave', function(e){
    if( !e.relatedTarget || !target.contains( e.relatedTarget ) )
      target.classList.remove( 'DragOver' );
  } );

  // The detector files in a dropped folder, recursively; `prefix` is the path of the folder within
  //  the drop.  Hidden files and folders, and other types of files, are counted in out.skipped.
  const walk = async function( entry, prefix, out ){
    if( entry.name.charAt(0) === '.' )
    {
      out.skipped += 1;
      return;
    }
    if( out.files.length >= MAX_FILES )
    {
      out.refused += 1;
      return;
    }
    if( entry.isFile )
    {
      if( !DRF_FILE.test( entry.name ) )
      {
        out.skipped += 1;
        return;
      }
      const file = await new Promise( function(resolve){ entry.file( resolve, function(){ resolve(null); } ); } );
      if( file )
        out.files.push( { file: file, path: prefix + entry.name } );
      else
        out.refused += 1;
    }else if( entry.isDirectory )
    {
      const reader = entry.createReader();
      let entries = [];
      do
      {
        // readEntries gives the entries a few at a time; an empty array means no more
        entries = await new Promise( function(resolve){ reader.readEntries( resolve, function(){ resolve([]); } ); } );
        for( const child of entries )
          await walk( child, prefix + entry.name + '/', out );
      }while( entries.length );
    }
  };

  target.addEventListener( 'drop', async function(e){
    e.preventDefault();
    e.stopPropagation();
    target.classList.remove( 'DragOver' );

    const dt = e.dataTransfer;
    const items = dt ? dt.items : null;
    if( items && items.length && items[0].webkitGetAsEntry )
    {
      // The entries must be taken before the first await: the DataTransfer is emptied after.
      const entries = [];
      for( let i = 0; i < items.length; ++i )
      {
        const entry = items[i].webkitGetAsEntry ? items[i].webkitGetAsEntry() : null;
        if( entry )
          entries.push( entry );
      }

      // A single folder, whose path the native app gave us (in the last 30 s)
      let native_folder = null;
      const fps = (window._IS && window._IS.dragOverFilePaths) ? window._IS.dragOverFilePaths : null;
      if( fps && fps.filenames && fps.time && ((new Date() - fps.time) < 30000)
          && (entries.length === 1) && entries[0].isDirectory )
      {
        const name = entries[0].name;
        for( const p of fps.filenames )
        {
          const trimmed = String(p).replace( /[\/\\]+$/, '' );
          const leaf = trimmed.split( /[\/\\]/ ).pop();
          if( leaf === name )
            native_folder = trimmed;
        }
      }

      // Dropped files (not folders) are all sent when there are only one or two - the server says
      //  what is wrong with a file of the wrong type - else only those that can be detector files.
      const out = { files: [], skipped: 0, refused: 0 };
      const loose = entries.filter( function(e){ return e.isFile; } ).length;
      for( const entry of entries )
      {
        if( entry.isFile && (loose <= 2) && (entries.length <= 2) && (entry.name.charAt(0) !== '.') )
        {
          const file = await new Promise( function(resolve){ entry.file( resolve, function(){ resolve(null); } ); } );
          if( file )
            out.files.push( { file: file, path: entry.name } );
          else
            out.refused += 1;
        }else
        {
          await walk( entry, '', out );
        }
      }
      upload( out, native_folder );
    }else if( dt && dt.files )
    {
      browsed( Array.from( dt.files ).map( function(f){ return { file: f, path: f.name }; } ) );
    }
  } );

  // Files from a browse (or a drop without folder support): filtered as for a drop.
  const browsed = function( list ){
    const many = (list.length > 2);
    const out = { files: [], skipped: 0, refused: 0 };
    for( const x of list )
    {
      if( /(^|\/)\./.test( x.path ) || (many && !DRF_FILE.test( x.path )) )
        out.skipped += 1;
      else if( out.files.length >= MAX_FILES )
        out.refused += 1;
      else
        out.files.push( x );
    }
    upload( out, null );
  };

  // Sends {files, skipped, refused} as one batch.  The final request (X-Drf-Done) is sent even
  //  when no file is, so the user is told that nothing in the drop could be used.
  const upload = function( out, native_folder ){
    let total = 0, refused = out.refused;
    const list = out.files.filter( function(x){
      if( !x.file || (x.file.size > MAX_FILE_BYTES) || ((total + x.file.size) > MAX_TOTAL_BYTES) )
      {
        refused += 1;
        return false;
      }
      total += x.file.size;
      return true;
    } );

    if( !list.length && !refused && !out.skipped )
      return;

    const batch = Date.now().toString(36) + Math.random().toString(36).slice(2);
    num_active += 1;
    target.classList.add( 'Uploading' );

    let next_index = 0, num_running = 0, num_done = 0, num_failed = 0;
    const finish = function(){
      const xhr = new XMLHttpRequest();
      xhr.open( 'POST', url, true );
      xhr.setRequestHeader( 'X-Drf-Batch', batch );
      xhr.setRequestHeader( 'X-Drf-Done', String( list.length ) );
      xhr.setRequestHeader( 'X-Drf-Failed', String( refused + num_failed ) );
      xhr.setRequestHeader( 'X-Drf-Skipped', String( out.skipped ) );
      if( native_folder )
        xhr.setRequestHeader( 'X-Drf-Native-Folder', encodeURIComponent( native_folder ) );
      xhr.onloadend = function(){
        num_active -= 1;
        if( num_active <= 0 )
          target.classList.remove( 'Uploading' );
      };
      xhr.send();
    };

    if( !list.length )
    {
      finish();
      return;
    }

    const start_next = function(){
      while( (num_running < NUM_PARALLEL) && (next_index < list.length) )
      {
        const item = list[next_index++];
        num_running += 1;
        const data = new FormData();
        data.append( 'file', item.file, item.file.name );
        const xhr = new XMLHttpRequest();
        xhr.open( 'POST', url, true );
        xhr.setRequestHeader( 'X-Drf-Batch', batch );
        xhr.setRequestHeader( 'X-Drf-Path', encodeURIComponent( item.path ) );
        xhr.onloadend = function(){
          if( (xhr.status < 200) || (xhr.status >= 300) )
            num_failed += 1;
          num_running -= 1;
          num_done += 1;
          if( num_done === list.length )
            finish();
          else
            start_next();
        };
        xhr.send( data );
      }
    };

    start_next();
  };
}
