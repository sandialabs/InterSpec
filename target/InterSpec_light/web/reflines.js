/* InterSpec Light - reference photopeak lines: pick sources from the precomputed library. */
'use strict';

App.ref = {
  lib: [],
  byName: new Map(),
  shown: [],   // [{ parent, color }]
};

// InterSpec's default reference-line colors
const REF_COLORS = [ '#00b9ff', '#006600', '#cc3333', '#9933FF', '#FF66FF', '#830808', '#FF6633',
                     '#F1C232', '#0000CC', '#666666', '#003333', '#CCFFCC' ];
const REF_SOURCE_COLORS = { background: '#967f55' };
const REF_STORAGE_KEY = 'interspec-light-ref-lines';

const KIND_LABELS = { nuclide: 'nuclide', reaction: 'reaction', xray: 'x-rays', background: 'NORM' };


App.ref.normalize = function( name ) {
  return String( name ).toLowerCase().replace( /[\s\-]/g, '' );
};


App.ref.init = function( refLinesText ) {
  App.ref.lib = JSON.parse( refLinesText );
  for( const entry of App.ref.lib )
    App.ref.byName.set( App.ref.normalize( entry.parent ), entry );

  if( App.ref.lib.length )
    App.call( 'setRefLibrary', { lib: App.ref.lib } );

  const input = document.getElementById( 'ref-filter' );
  const sugg = document.getElementById( 'ref-suggestions' );
  let active = -1;

  const matches = () => {
    const q = App.ref.normalize( input.value );
    if( !q )
      return [];
    const shown = new Set( App.ref.shown.map( (s) => s.parent ) );
    const scored = [];
    for( const entry of App.ref.lib )
    {
      if( shown.has( entry.parent ) )
        continue;
      const n = App.ref.normalize( entry.parent );
      const pos = n.indexOf( q );
      if( pos >= 0 )
        scored.push( [ (n === q) ? 0 : ((pos === 0) ? 1 : 2), n.length, entry ] );
    }
    scored.sort( (a, b) => (a[0] - b[0]) || (a[1] - b[1]) );
    return scored.slice( 0, 40 ).map( (s) => s[2] );
  };

  const renderSuggestions = () => {
    const list = matches();
    sugg.innerHTML = '';
    sugg.hidden = !list.length;
    active = Math.min( active, list.length - 1 );
    list.forEach( (entry, i) => {
      const div = App.el( 'div', { class: (i === active) ? 'active' : '' }, [
        App.el( 'span', { text: entry.parent } ),
        App.el( 'span', { class: 'kind', text: (KIND_LABELS[entry.kind] || entry.kind) + (entry.age ? ', ' + entry.age : '') } )
      ] );
      div.addEventListener( 'mousedown', (e) => {
        e.preventDefault();
        App.ref.add( entry.parent );
        input.value = '';
        renderSuggestions();
      } );
      sugg.appendChild( div );
    } );
  };

  input.addEventListener( 'input', () => { active = 0; renderSuggestions(); } );
  input.addEventListener( 'blur', () => { sugg.hidden = true; } );
  input.addEventListener( 'focus', renderSuggestions );
  input.addEventListener( 'keydown', (e) => {
    const list = matches();
    if( e.key === 'ArrowDown' ) { active = Math.min( active + 1, list.length - 1 ); renderSuggestions(); e.preventDefault(); }
    else if( e.key === 'ArrowUp' ) { active = Math.max( active - 1, 0 ); renderSuggestions(); e.preventDefault(); }
    else if( (e.key === 'Enter') && list.length )
    {
      App.ref.add( list[Math.max( active, 0 )].parent );
      input.value = '';
      renderSuggestions();
    }else if( e.key === 'Escape' )
    {
      input.value = '';
      renderSuggestions();
    }
  } );

  // Restore the previous session's reference lines (a per-browser convenience)
  try
  {
    const saved = JSON.parse( localStorage.getItem( REF_STORAGE_KEY ) || '[]' );
    for( const s of saved )
    {
      const entry = App.ref.byName.get( App.ref.normalize( s.parent ) );
      if( entry )
        App.ref.shown.push( { parent: entry.parent, color: s.color } );
    }
  }catch( e )
  {
  }
  App.ref.update();
};


/** Names of the library's nuclides, elements, and reactions, for typing peak sources. */
App.ref.sourceNames = function() {
  if( !App.ref.names )
  {
    const names = new Set();
    for( const entry of App.ref.lib )
    {
      if( entry.kind === 'reaction' )
        (entry.desc_strs || []).forEach( (name) => names.add( name ) );
      else if( entry.kind !== 'background' )
        names.add( entry.parent );
    }
    App.ref.names = Array.from( names );
  }
  return App.ref.names;
};


App.ref.nextColor = function( parent ) {
  const fixed = REF_SOURCE_COLORS[ App.ref.normalize( parent ) ];
  if( fixed )
    return fixed;
  const used = new Set( App.ref.shown.map( (s) => s.color ) );
  return REF_COLORS.find( (c) => !used.has( c ) ) || REF_COLORS[ App.ref.shown.length % REF_COLORS.length ];
};


App.ref.isShown = function( parent ) {
  return App.ref.shown.some( (s) => s.parent === parent );
};


App.ref.add = function( parent ) {
  const entry = App.ref.byName.get( App.ref.normalize( parent ) );
  if( !entry || App.ref.isShown( entry.parent ) )
    return;
  App.ref.shown.push( { parent: entry.parent, color: App.ref.nextColor( entry.parent ) } );
  App.ref.update();
};


App.ref.remove = function( parent ) {
  App.ref.shown = App.ref.shown.filter( (s) => s.parent !== parent );
  App.ref.update();
};


App.ref.toggle = function( parent ) {
  if( App.ref.isShown( parent ) )
    App.ref.remove( parent );
  else
    App.ref.add( parent );
};


/** Pushes the shown sources to the chart, the WASM (for assigning peak sources), and the chips. */
App.ref.update = function() {
  const chartData = App.ref.shown.map( (s) => Object.assign( {}, App.ref.byName.get( App.ref.normalize( s.parent ) ), { color: s.color } ) );
  if( App.spectrumChart )
    App.spectrumChart.setReferenceLines( chartData.length ? chartData : null );

  try
  {
    App.call( 'setShownRefLines', { sources: App.ref.shown } );
  }catch( e )
  {
    console.error( e );
  }

  const chips = document.getElementById( 'ref-chips' );
  chips.innerHTML = '';
  for( const s of App.ref.shown )
  {
    chips.appendChild( App.el( 'span', { class: 'chip', style: 'border-color:' + s.color }, [
      s.parent,
      App.el( 'button', { text: '✕', title: 'Remove', onclick: () => App.ref.remove( s.parent ) } )
    ] ) );
  }

  try
  {
    localStorage.setItem( REF_STORAGE_KEY, JSON.stringify( App.ref.shown ) );
  }catch( e )
  {
  }

  if( App.search )
    App.search.renderResults();
};
