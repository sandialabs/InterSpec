/* InterSpec Light - search the reference-line library for sources with lines near one or two energies. */
'use strict';

App.search = {
  results: [],
};


/** Default half-width of the search window: ~FWHM for HPGe, or ~5% for lower resolution detectors. */
App.search.defaultWindow = function( energy ) {
  return App.isHighRes ? 1.5 : (1.0 + 0.05*energy);
};


App.search.inputs = function() {
  const val = (id) => parseFloat( document.getElementById( id ).value );
  const energies = [ val( 'search-e1' ), val( 'search-e2' ) ].filter( (e) => isFinite( e ) && (e > 0) );
  const userWindow = val( 'search-window' );
  return energies.map( (e) => ({ energy: e, window: (isFinite( userWindow ) && (userWindow > 0)) ? userWindow : App.search.defaultWindow( e ) }) );
};


/** Scores library sources: each wanted energy must have a line within its window.
 "closest" sums |dE|/window; "weighted" sums (0.25*window + |dE|)/intensity, the same metric used
 to assign sources to peaks (x-rays de-weighted by 0.2). */
App.search.compute = function( wanted, sort ) {
  const results = [];
  for( const entry of App.ref.lib )
  {
    let score = 0;
    const matched = [];
    let ok = true;
    for( const w of wanted )
    {
      let best = null, bestScore = Infinity;
      for( const line of entry.lines )
      {
        const de = Math.abs( line.e - w.energy );
        if( de > w.window )
          continue;
        const intensity = line.h * ((line.particle === 'xray') ? 0.2 : 1.0);
        const s = (sort === 'closest') ? (de / w.window) : ((0.25*w.window + de) / Math.max( intensity, 1e-12 ));
        if( s < bestScore )
        {
          bestScore = s;
          best = line;
        }
      }//for( loop over lines )

      if( !best )
      {
        ok = false;
        break;
      }
      score += bestScore;
      matched.push( { line: best, de: best.e - w.energy } );
    }//for( loop over wanted energies )

    if( ok )
      results.push( { entry: entry, score: score, matched: matched } );
  }//for( loop over library )

  results.sort( (a, b) => a.score - b.score );
  return results;
};


App.search.renderResults = function() {
  const container = document.getElementById( 'search-results' );
  const section = document.getElementById( 'search-section' );
  const wanted = App.search.inputs();
  container.innerHTML = '';

  if( App.spectrumChart )
    App.spectrumChart.setSearchWindows( (section.open && wanted.length) ? wanted : null );

  document.getElementById( 'search-window' ).placeholder =
    wanted.length ? ('auto: ' + App.fmt( wanted[0].window, 2 )) : 'auto';

  if( !wanted.length )
    return;

  const sort = document.getElementById( 'search-sort' ).value;
  const results = App.search.compute( wanted, sort ).slice( 0, 60 );
  if( !results.length )
  {
    container.appendChild( App.el( 'p', { class: 'hint', text: 'No sources with lines in the search window.' } ) );
    return;
  }

  // Two-energy searches list two values per cell, so let the cells wrap
  const table = App.el( 'table', { class: 'data wrap' } );
  table.appendChild( App.el( 'tr', {}, [
    App.el( 'th', { text: 'Source' } ), App.el( 'th', { text: 'Line (keV)' } ),
    App.el( 'th', { class: 'num', text: 'Rel. I' } ), App.el( 'th', { class: 'num', text: 'ΔE' } )
  ] ) );

  for( const r of results )
  {
    const shown = App.ref.isShown( r.entry.parent );
    const lines = r.matched.map( (m) => App.fmt( m.line.e, 2 ) + ((m.line.particle === 'xray') ? ' x' : '') ).join( ', ' );
    const intensities = r.matched.map( (m) => (m.line.h >= 0.01) ? m.line.h.toFixed( 2 ) : m.line.h.toExponential( 1 ) ).join( ', ' );
    const des = r.matched.map( (m) => ((m.de >= 0) ? '+' : '') + m.de.toFixed( 2 ) ).join( ', ' );
    const color = shown ? App.ref.shown.find( (s) => s.parent === r.entry.parent ).color : null;

    table.appendChild( App.el( 'tr', { class: 'clickable', title: shown ? 'Click to hide lines' : 'Click to show lines',
                                       onclick: () => App.ref.toggle( r.entry.parent ) }, [
      App.el( 'td', {}, [ color ? App.el( 'span', { class: 'swatch', style: 'background:' + color } ) : '',
                          App.el( shown ? 'b' : 'span', { text: r.entry.parent } ) ] ),
      App.el( 'td', { text: lines } ),
      App.el( 'td', { class: 'num', text: intensities } ),
      App.el( 'td', { class: 'num', text: des } )
    ] ) );
  }
  container.appendChild( table );
};


/** Opens the search section, searching for the given energy (e.g., from the right-click menu). */
App.search.setEnergy = function( energy ) {
  const section = document.getElementById( 'search-section' );
  section.open = true;
  document.getElementById( 'search-e1' ).value = energy.toFixed( 2 );
  document.getElementById( 'search-e2' ).value = '';
  App.search.renderResults();
  section.scrollIntoView( { block: 'nearest' } );
};


App.search.onForegroundChanged = function() {
  App.search.renderResults();
};


App.search.init = function() {
  let timer = null;
  const later = () => { clearTimeout( timer ); timer = setTimeout( App.search.renderResults, 150 ); };
  for( const id of [ 'search-e1', 'search-e2', 'search-window' ] )
    document.getElementById( id ).addEventListener( 'input', later );
  document.getElementById( 'search-sort' ).addEventListener( 'change', App.search.renderResults );
  document.getElementById( 'search-section' ).addEventListener( 'toggle', App.search.renderResults );
};
