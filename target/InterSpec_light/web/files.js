/* InterSpec Light - loading files (picker, drag-n-drop, URL), the spectrum slots, detectors, and export. */
'use strict';

App.files = {
  slots: {},
  loadCounter: 0,
  openType: 'FOREGROUND',  //spectrum type the file picker is loading
};

const SPEC_TYPES = [ 'FOREGROUND', 'BACKGROUND', 'SECONDARY' ];
const SPEC_LABELS = { FOREGROUND: 'Foreground', BACKGROUND: 'Background', SECONDARY: 'Secondary' };


/** "1,2,3,5,7,8" -> "1-3,5,7-8" */
App.files.formatSamples = function( samples ) {
  const s = (samples || []).slice().sort( (a, b) => a - b );
  const parts = [];
  for( let i = 0; i < s.length; )
  {
    let j = i;
    while( (j + 1 < s.length) && (s[j + 1] === s[j] + 1) )
      ++j;
    parts.push( (j > i) ? (s[i] + '-' + s[j]) : String( s[i] ) );
    i = j + 1;
  }
  return parts.join( ',' );
};


/** "1-3,5" -> [1,2,3,5]; throws on bad input. */
App.files.parseSamples = function( str ) {
  const answer = [];
  for( const part of str.split( ',' ) )
  {
    const t = part.trim();
    if( !t )
      continue;
    const m = t.match( /^(-?\d+)\s*(?:-\s*(-?\d+))?$/ );
    if( !m )
      throw new Error( 'Invalid sample range "' + t + '"' );
    const a = parseInt( m[1] ), b = (m[2] !== undefined) ? parseInt( m[2] ) : a;
    for( let i = Math.min( a, b ); i <= Math.max( a, b ); ++i )
      answer.push( i );
  }
  return answer;
};


/** Loads raw file bytes into the given display slot. */
App.files.loadBytes = async function( bytes, name, type ) {
  if( !App.hasForeground )
    type = 'FOREGROUND';

  const path = '/in/' + (App.files.loadCounter++) + '_' + name.replace( /[^\w.\-]/g, '_' );
  App.module.FS.writeFile( path, bytes );
  try
  {
    const result = await App.run( 'loadFile', { path: path, name: name, type: type }, { busy: true } );
    if( result && (type === 'FOREGROUND') )
      App.search.onForegroundChanged();
  }finally
  {
    try{ App.module.FS.unlink( path ); }catch( e ){ }
  }
};


App.files.loadFile = function( file, type ) {
  const reader = new FileReader();
  reader.onload = () => App.files.loadBytes( new Uint8Array( reader.result ), file.name, type );
  reader.onerror = () => App.toast( 'Could not read ' + file.name, true );
  reader.readAsArrayBuffer( file );
};


App.files.loadUrl = async function( url, type ) {
  try
  {
    const response = await fetch( url );
    if( !response.ok )
      throw new Error( 'HTTP ' + response.status );
    const bytes = new Uint8Array( await response.arrayBuffer() );
    let name = decodeURIComponent( new URL( url, window.location.href ).pathname.split( '/' ).pop() || 'spectrum' );
    await App.files.loadBytes( bytes, name, type );
  }catch( e )
  {
    App.toast( 'Failed to fetch ' + url + ': ' + e.message
               + ((window.location.protocol === 'file:') ? ' (the server must allow cross-origin requests)' : ''), true );
  }
};


App.files.renderSlots = function( files ) {
  const container = document.getElementById( 'slots' );
  container.innerHTML = '';

  for( const type of SPEC_TYPES )
  {
    const info = files[type];
    const head = App.el( 'div', { class: 'slot-head' }, [ App.el( 'span', { class: 'slot-type ' + type, text: SPEC_LABELS[type] } ) ] );
    const slot = App.el( 'div', { class: 'slot' }, [ head ] );

    if( !info )
    {
      head.appendChild( App.el( 'span', { class: 'slot-empty', text: 'none' } ) );
      if( App.hasForeground || (type === 'FOREGROUND') )
        head.appendChild( App.el( 'button', { class: 'btn-sm', text: 'Open…', onclick: () => App.files.openPicker( type ) } ) );
      container.appendChild( slot );
      continue;
    }

    head.appendChild( App.el( 'span', { class: 'slot-name', text: info.name, title: info.name } ) );
    head.appendChild( App.el( 'button', { class: 'icon-btn', text: '📂', title: 'Open a different ' + SPEC_LABELS[type].toLowerCase() + ' file',
                                          onclick: () => App.files.openPicker( type ) } ) );
    head.appendChild( App.el( 'button', { class: 'icon-btn', text: '✕', title: 'Remove',
                                          onclick: () => App.run( 'unload', { type: type } ) } ) );

    const details = [];
    if( typeof info.liveTime === 'number' )
      details.push( 'LT ' + App.fmt( info.liveTime, 1 ) + ' s' );
    if( typeof info.realTime === 'number' )
      details.push( 'RT ' + App.fmt( info.realTime, 1 ) + ' s' );
    if( info.numChannels )
      details.push( info.numChannels + ' ch' );
    if( info.instrument && (type === 'FOREGROUND') )
      details.push( info.instrument );
    slot.appendChild( App.el( 'div', { class: 'slot-detail', text: details.join( ' · ' ) } ) );

    if( info.numSamples > 1 )
    {
      const input = App.el( 'input', { type: 'text', value: App.files.formatSamples( info.samples ),
                                       title: 'Sample numbers to sum, e.g. 1-3,5' } );
      const apply = () => {
        try
        {
          App.run( 'setSamples', { type: type, samples: App.files.parseSamples( input.value ) }, { busy: true } );
        }catch( e )
        {
          App.toast( e.message, true );
        }
      };
      input.addEventListener( 'change', apply );
      input.addEventListener( 'keydown', (e) => { if( e.key === 'Enter' ) input.blur(); } );

      const row = App.el( 'div', { class: 'slot-samples' }, [
        App.el( 'span', { text: 'Samples' } ),
        input,
        App.el( 'button', { class: 'btn-sm', text: '◀', title: 'Previous sample',
                            onclick: () => App.run( 'stepSample', { type: type, delta: -1 }, { busy: true } ) } ),
        App.el( 'button', { class: 'btn-sm', text: '▶', title: 'Next sample',
                            onclick: () => App.run( 'stepSample', { type: type, delta: 1 }, { busy: true } ) } )
      ] );
      slot.appendChild( row );
      slot.appendChild( App.el( 'div', { class: 'slot-detail',
        text: info.numSamples + ' samples (' + info.firstSample + '–' + info.lastSample + ')'
              + ((info.passthrough && (type === 'FOREGROUND')) ? '; drag on the time chart to select' : '') } ) );
    }

    container.appendChild( slot );
  }//for( loop over types )

  const fore = files.FOREGROUND;
  document.getElementById( 'fg-summary' ).textContent = fore ? fore.name : '';
};


App.files.renderDetectors = function( dets ) {
  const section = document.getElementById( 'detectors-section' );
  const container = document.getElementById( 'detectors' );
  container.innerHTML = '';
  section.hidden = !dets.all || (dets.all.length < 2);

  for( const d of dets.all || [] )
  {
    const cb = App.el( 'input', { type: 'checkbox', value: d.name } );
    cb.checked = dets.shown.includes( d.name );
    cb.addEventListener( 'change', App.files.sendDetectors );
    container.appendChild( App.el( 'label', { class: d.neutronOnly ? 'neutron' : '',
                                   title: d.neutronOnly ? 'Neutron-only detector' : '' },
                                   [ cb, ' ' + d.name + (d.neutronOnly ? ' (n)' : '') ] ) );
  }
};


App.files.sendDetectors = function() {
  const shown = [];
  for( const cb of document.querySelectorAll( '#detectors input[type=checkbox]' ) )
  {
    if( cb.checked )
      shown.push( cb.value );
  }
  App.run( 'setDetectors', { shown: shown }, { busy: true } );
};


App.files.openPicker = function( type ) {
  App.files.openType = type || 'FOREGROUND';
  document.getElementById( 'file-input' ).click();
};


/** Asks the WASM to write a file (method returns {path, filename}), and downloads it. */
App.files.download = function( method, params ) {
  const path = '/tmp/interspec_light_download';
  try
  {
    const result = App.call( method, Object.assign( { path: path }, params || {} ) );
    const data = App.module.FS.readFile( result.path );
    App.module.FS.unlink( result.path );
    const url = URL.createObjectURL( new Blob( [ data ], { type: 'application/octet-stream' } ) );
    const a = App.el( 'a', { href: url, download: result.filename } );
    document.body.appendChild( a );
    a.click();
    a.remove();
    setTimeout( () => URL.revokeObjectURL( url ), 5000 );
    App.toast( 'Saved ' + result.filename );
  }catch( e )
  {
    App.toast( 'Export failed: ' + e.message, true );
  }
};


App.files.exportFile = function() {
  App.files.download( 'exportFile', { format: document.getElementById( 'export-format' ).value } );
};


App.files.init = function() {
  const input = document.getElementById( 'file-input' );
  input.addEventListener( 'change', () => {
    if( input.files.length )
      App.files.loadFile( input.files[0], App.files.openType );
    input.value = '';
  } );
  document.getElementById( 'empty-open' ).addEventListener( 'click', () => App.files.openPicker( 'FOREGROUND' ) );
  document.getElementById( 'export-btn' ).addEventListener( 'click', App.files.exportFile );
  document.getElementById( 'export-peaks-btn' ).addEventListener( 'click', () => App.files.download( 'exportPeakCsv' ) );
  document.getElementById( 'dets-all' ).addEventListener( 'click', () => {
    document.querySelectorAll( '#detectors input' ).forEach( (cb) => { cb.checked = true; } );
    App.files.sendDetectors();
  } );
  document.getElementById( 'dets-none' ).addEventListener( 'click', () => {
    document.querySelectorAll( '#detectors input' ).forEach( (cb) => { cb.checked = false; } );
    App.files.sendDetectors();
  } );

  // Drag-n-drop: the overlay offers a zone for each spectrum type
  const overlay = document.getElementById( 'drop-overlay' );
  let depth = 0;
  const hasFiles = (e) => e.dataTransfer && Array.from( e.dataTransfer.types || [] ).includes( 'Files' );
  document.addEventListener( 'dragenter', (e) => {
    if( !hasFiles( e ) )
      return;
    e.preventDefault();
    depth += 1;
    overlay.hidden = false;
    for( const zone of overlay.querySelectorAll( '.drop-zone' ) )
      zone.hidden = !App.hasForeground && (zone.dataset.type !== 'FOREGROUND');
  } );
  document.addEventListener( 'dragleave', (e) => {
    if( !hasFiles( e ) )
      return;
    depth = Math.max( 0, depth - 1 );
    if( !depth )
      overlay.hidden = true;
  } );
  document.addEventListener( 'dragover', (e) => {
    if( !hasFiles( e ) )
      return;
    e.preventDefault();
    for( const zone of overlay.querySelectorAll( '.drop-zone' ) )
      zone.classList.toggle( 'over', zone.contains( e.target ) );
  } );
  document.addEventListener( 'drop', (e) => {
    if( !hasFiles( e ) )
      return;
    e.preventDefault();
    depth = 0;
    overlay.hidden = true;
    const zone = e.target.closest ? e.target.closest( '.drop-zone' ) : null;
    const type = zone ? zone.dataset.type : 'FOREGROUND';
    if( e.dataTransfer.files.length )
      App.files.loadFile( e.dataTransfer.files[0], type );
  } );

  App.renderers.push( (state) => {
    if( state.files )
    {
      App.files.slots = state.files;
      App.hasForeground = !!state.files.FOREGROUND;
      App.files.renderSlots( state.files );
    }
    if( state.detectors )
      App.files.renderDetectors( state.detectors );
  } );
};


/** Loads any spectra given in the page URL: ?url=...&background=...&secondary=... */
App.files.loadFromPageUrl = async function() {
  const params = new URLSearchParams( window.location.search );
  const fore = params.get( 'url' ) || params.get( 'foreground' );
  if( fore )
    await App.files.loadUrl( fore, 'FOREGROUND' );
  if( fore && params.get( 'background' ) )
    await App.files.loadUrl( params.get( 'background' ), 'BACKGROUND' );
  if( fore && params.get( 'secondary' ) )
    await App.files.loadUrl( params.get( 'secondary' ), 'SECONDARY' );
};
