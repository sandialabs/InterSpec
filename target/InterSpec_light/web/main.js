/* InterSpec Light - startup. */
'use strict';

/** Reads a data block embedded in the page by assemble_html.py.  Its `data-encoding` is "text"
 (returned as-is), "base64", or "gzip-base64" (decompressed with the browser's DecompressionStream);
 the latter two are returned as an ArrayBuffer. */
App.readAsset = async function( id ) {
  const el = document.getElementById( id );
  const encoding = el.getAttribute( 'data-encoding' );
  if( encoding === 'text' )
    return el.textContent;

  const bin = atob( el.textContent.trim() );
  const bytes = new Uint8Array( bin.length );
  for( let i = 0; i < bin.length; ++i )
    bytes[i] = bin.charCodeAt( i );
  if( encoding === 'base64' )
    return bytes.buffer;

  if( typeof DecompressionStream === 'undefined' )
    throw new Error( 'this browser can not decompress the page (it needs Safari 16.4+, Chrome 80+,'
                     + ' or Firefox 113+) - please use the uncompressed version of InterSpec Light' );
  const stream = new Blob( [ bytes ] ).stream().pipeThrough( new DecompressionStream( 'gzip' ) );
  return new Response( stream ).arrayBuffer();
};


(async function main() {
  let refLinesText;
  try
  {
    const [ wasmBinary, refLines ] = await Promise.all( [ App.readAsset( 'light-wasm' ), App.readAsset( 'light-ref-lines' ) ] );
    refLinesText = (typeof refLines === 'string') ? refLines : new TextDecoder().decode( refLines );

    // The fitting code prints diagnostics to stderr; keep them out of the error log.
    App.module = await InterSpecLightModule( { printErr: (txt) => console.debug( txt ), wasmBinary: wasmBinary } );
  }catch( e )
  {
    console.error( e );
    const msg = document.getElementById( 'empty-msg' );
    msg.textContent = 'Failed to load InterSpec Light: ' + (e.message || e);
    msg.style.color = 'var(--danger)';
    return;
  }

  App.module.FS.mkdirTree( '/in' );
  App.lightCall = App.module.cwrap( 'light_call', 'string', [ 'string' ] );

  App.initCharts();
  App.files.init();
  App.ref.init( refLinesText );
  App.search.init();
  App.peaks.init();
  App.ecal.init();

  const panel = document.getElementById( 'panel' );
  if( window.innerWidth < 800 )
    panel.classList.add( 'collapsed' );
  document.getElementById( 'panel-toggle' ).addEventListener( 'click', () => {
    panel.classList.toggle( 'collapsed' );
    requestAnimationFrame( () => App.spectrumChart.handleResize() );
  } );

  const initialState = App.call( 'getState' );
  App.applyState( initialState );
  App.features = initialState.features || {};
  // Text that depends on whether this build reads and writes peaks in N42/SPE files
  for( const el of document.querySelectorAll( '.peak-files' ) )
    el.hidden = !App.features.peakFileIo;
  for( const el of document.querySelectorAll( '.no-peak-files' ) )
    el.hidden = !!App.features.peakFileIo;
  await App.files.loadFromPageUrl();
})();
