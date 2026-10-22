/* InterSpec Light - peak table, fitting from the chart, and the right-click peak menu. */
'use strict';

App.peaks = {
  list: [],
};

// Values of PeakContinuum::OffsetType and PeakDef::SkewType
const CONTINUUM_TYPES = [
  [0, 'None'], [1, 'Constant'], [2, 'Linear'], [3, 'Quadratic'], [4, 'Cubic'],
  [5, 'Flat Step'], [6, 'Linear Step'], [7, 'Bi-Linear Step'],
  [8, 'Flat Step (CDF)'], [9, 'Linear Step (CDF)'], [10, 'Bi-Linear Step (CDF)'], [11, 'Global Continuum']
];
const SKEW_TYPES = [
  [0, 'None'], [1, 'Bortel'], [2, 'Gauss+Exp'], [3, 'Crystal Ball'], [4, 'Exp+Gauss+Exp'],
  [5, 'Double Crystal Ball'], [6, 'Voigt+Bortel'], [7, 'Gauss+Bortel'], [8, 'Double Bortel'],
  [9, 'GADRAS'], [10, 'GADRAS CZT']
];


App.peaks.typeLabel = function( types, value ) {
  const t = types.find( (tv) => tv[0] === value );
  return t ? t[1] : String( value );
};


App.peaks.fitAt = function( energy, refLineParent ) {
  if( !App.hasForeground )
    return;
  App.run( 'fitPeakAt', { energy: energy, pixPerKeV: App.pixelsPerKeV(), refLineParent: refLineParent || '' }, { busy: true } );
};


App.peaks.render = function( list ) {
  App.peaks.list = list;
  document.getElementById( 'peak-count' ).textContent = list.length ? '(' + list.length + ')' : '';
  for( const id of [ 'peaks-csv', 'peaks-clear', 'export-peaks-btn' ] )
    document.getElementById( id ).disabled = !list.length;
  document.getElementById( 'peaks-search' ).disabled = !App.hasForeground;

  const container = document.getElementById( 'peak-table' );
  container.innerHTML = '';
  if( !list.length )
  {
    container.appendChild( App.el( 'p', { class: 'hint', text: App.hasForeground
      ? 'Double-click on the spectrum to fit a peak.' : 'Load a spectrum to fit peaks.' } ) );
    return;
  }

  const table = App.el( 'table', { class: 'data' } );
  table.appendChild( App.el( 'tr', {}, [
    App.el( 'th', { text: 'Energy' } ), App.el( 'th', { class: 'num', text: 'FWHM' } ),
    App.el( 'th', { class: 'num', text: 'Area' } ), App.el( 'th', { text: 'Source' } ),
    App.el( 'th', { text: 'Cal', title: 'Use for energy calibration' } ), App.el( 'th' )
  ] ) );

  for( const p of list )
  {
    const useCb = App.el( 'input', { type: 'checkbox', title: 'Use for energy calibration (needs a source)' } );
    useCb.checked = !!p.source && p.wantUseForCal;
    useCb.disabled = !p.source;
    if( useCb.checked && !p.useForCal )  //e.g., a peak from a file InterSpec saved, with its mean fixed
    {
      useCb.indeterminate = true;
      useCb.title = 'Not used for energy calibration: this peak\'s mean was fixed (in InterSpec)';
    }
    useCb.addEventListener( 'click', (e) => e.stopPropagation() );
    useCb.addEventListener( 'change', () => App.run( 'setPeakProperty', { mean: p.mean, useForCal: useCb.checked } ) );

    const areaTxt = (p.area >= 1e5) ? p.area.toExponential( 3 ) : p.area.toPrecision( 4 );

    // The source column takes the remaining width, truncating long labels (full label on hover);
    //  clicking it lets the user type a source, as in InterSpec's peak table
    const sourceCell = App.el( 'td', { class: 'fill source', title: (p.source ? p.source + ' - click' : 'Click')
                                       + ' to type a source, e.g. Cs137, Pb 75, or Th232 S.E.' },
      p.source ? [ App.el( 'div', { class: 'fill-row' }, [
        p.color ? App.el( 'span', { class: 'swatch', style: 'background:' + p.color } ) : '',
        App.el( 'span', { class: 'ellipsis', text: App.peaks.shortSource( p ) } ),
        App.el( 'button', { class: 'icon-btn', text: '✕', title: 'Clear source',
          onclick: (e) => { e.stopPropagation(); App.run( 'setPeakProperty', { mean: p.mean, clearSource: true } ); } } )
      ] ) ] : [ App.el( 'span', { class: 'set-source', text: 'set…' } ) ] );
    sourceCell.addEventListener( 'click', (e) => { e.stopPropagation(); App.peaks.editSource( sourceCell, p ); } );

    const tip = 'Mean ' + App.fmt( p.mean, 3 ) + ' ± ' + App.fmt( p.meanUnc, 3 ) + ' keV; area ' + p.area.toPrecision( 6 )
                + ' ± ' + p.areaUnc.toPrecision( 3 ) + '; ' + App.peaks.typeLabel( CONTINUUM_TYPES, p.continuumType ) + ' continuum, '
                + App.peaks.typeLabel( SKEW_TYPES, p.skewType ) + ' skew'
                + ((p.chi2dof > 0) ? ('; χ²/dof ' + p.chi2dof.toFixed( 2 )) : '') + ' - click to zoom to peak';
    table.appendChild( App.el( 'tr', { class: 'clickable', title: tip,
                                       onclick: () => App.peaks.zoomTo( p ) }, [
      App.el( 'td', { class: 'num', text: App.fmt( p.mean, 2 ) } ),
      App.el( 'td', { class: 'num', text: App.fmt( p.fwhm, 2 ) } ),
      App.el( 'td', { class: 'num', text: areaTxt } ),
      sourceCell,
      App.el( 'td', {}, [ useCb ] ),
      App.el( 'td', {}, [ App.el( 'button', { class: 'icon-btn', text: '🗑', title: 'Delete peak',
        onclick: (e) => { e.stopPropagation(); App.run( 'deletePeakAt', { energy: p.mean } ); } } ) ] )
    ] ) );
  }//for( loop over peaks )

  container.appendChild( table );
};


/** The peak's source label, without the "keV" (still reads back as the same source). */
App.peaks.shortSource = function( p ) {
  return p.source.replace( / keV/, '' );
};


/** Replaces the peak's source cell with a text box to type a source into, or pick one of the
 suggestions (the shown reference lines' sources for the peak, then the library's names). */
App.peaks.editSource = function( cell, p ) {
  if( cell.querySelector( 'input' ) )
    return;

  const suggestions = document.getElementById( 'source-suggestions' );
  suggestions.innerHTML = '';
  let candidates = [];
  try
  {
    const info = App.call( 'peakInfoAt', { energy: p.mean } ).peak;
    candidates = (info && info.candidates) || [];
  }catch( e )
  {
    console.error( e );
  }
  const options = candidates.concat( App.ref.sourceNames() );
  for( const option of options )
    suggestions.appendChild( App.el( 'option', { value: option } ) );

  const original = App.peaks.shortSource( p );
  const input = App.el( 'input', { type: 'text', class: 'source-edit', list: 'source-suggestions',
                                   autocomplete: 'off', autocorrect: 'off', autocapitalize: 'off',
                                   spellcheck: 'false', placeholder: 'e.g. Cs137, Pb 75, Th232 S.E.' } );
  input.value = original;

  let done = false;
  const finish = (commit) => {
    if( done )
      return;
    done = true;
    const text = input.value.trim();
    if( commit && (text !== original) )
      App.run( 'setPeakProperty', { mean: p.mean, source: text } ).then( (r) => { if( !r ) App.peaks.render( App.peaks.list ); } );
    else
      App.peaks.render( App.peaks.list );
  };

  input.addEventListener( 'keydown', (e) => {
    e.stopPropagation();  //keep typing from reaching the charts' key handlers
    if( e.key === 'Enter' )
      finish( true );
    else if( e.key === 'Escape' )
      finish( false );
  } );
  // Picking a suggestion assigns it
  input.addEventListener( 'input', (e) => {
    if( (e.inputType === 'insertReplacementText') && options.includes( input.value ) )
      finish( true );
  } );
  input.addEventListener( 'blur', () => finish( true ) );
  input.addEventListener( 'click', (e) => e.stopPropagation() );

  cell.innerHTML = '';
  cell.appendChild( input );
  input.focus();
  input.select();
};


App.peaks.zoomTo = function( p ) {
  const chart = App.spectrumChart;
  const width = Math.max( p.upper - p.lower, 3*p.fwhm, 1 );
  chart.setXAxisRange( p.lower - 1.5*width, p.upper + 1.5*width, false, true );
  chart.redraw()();
};


/** Shows the right-click menu for the spectrum at `energy`. */
App.peaks.showContextMenu = function( energy, pageX, pageY ) {
  if( !App.hasForeground )
    return;

  let info = null;
  try
  {
    info = App.call( 'peakInfoAt', { energy: energy } ).peak;
  }catch( e )
  {
    App.toast( e.message, true );
    return;
  }

  const menu = document.getElementById( 'ctx-menu' );
  menu.innerHTML = '';
  const hide = () => { menu.hidden = true; };
  const item = (label, fn, cls) => App.el( 'div', { class: 'item' + (cls ? ' ' + cls : ''), text: label,
    onclick: (e) => { e.stopPropagation(); if( fn ){ hide(); fn(); } } } );

  if( info )
  {
    menu.appendChild( App.el( 'div', { class: 'title', text: 'Peak at ' + info.mean.toFixed( 2 ) + ' keV'
                                                     + (info.source ? (' — ' + info.source) : '') } ) );
    menu.appendChild( item( 'Delete peak', () => App.run( 'deletePeakAt', { energy: energy } ) ) );
    menu.appendChild( item( 'Refit ROI', () => App.run( 'refitRoi', { energy: energy }, { busy: true } ) ) );

    const contSub = App.el( 'div', { class: 'ctx-menu submenu' } );
    for( const [value, label] of CONTINUUM_TYPES )
      contSub.appendChild( item( label, () => App.run( 'setContinuumType', { energy: energy, type: value }, { busy: true } ),
                                 (value === info.continuumType) ? 'checked' : '' ) );
    const cont = item( 'Continuum', null, 'sub' );
    cont.appendChild( contSub );
    menu.appendChild( cont );

    if( info.gaussian )
    {
      const skewSub = App.el( 'div', { class: 'ctx-menu submenu' } );
      for( const [value, label] of SKEW_TYPES )
        skewSub.appendChild( item( label, () => App.run( 'setSkewType', { energy: energy, type: value }, { busy: true } ),
                                   (value === info.skewType) ? 'checked' : '' ) );
      const skew = item( 'Skew', null, 'sub' );
      skew.appendChild( skewSub );
      menu.appendChild( skew );
    }

    // Sources of the shown reference lines near the peak
    const candidates = info.candidates || [];
    if( candidates.length || info.source )
      menu.appendChild( App.el( 'div', { class: 'sep' } ) );
    for( const candidate of candidates )
      menu.appendChild( item( 'Assign as ' + candidate, () => App.run( 'setPeakProperty', { mean: info.mean, source: candidate } ) ) );
    if( info.source )
      menu.appendChild( item( 'Clear source', () => App.run( 'setPeakProperty', { mean: info.mean, clearSource: true } ) ) );
    menu.appendChild( App.el( 'div', { class: 'sep' } ) );
    menu.appendChild( item( 'Search energy ' + info.mean.toFixed( 1 ) + ' keV', () => App.search.setEnergy( info.mean ) ) );
  }else
  {
    menu.appendChild( App.el( 'div', { class: 'title', text: energy.toFixed( 1 ) + ' keV' } ) );
    menu.appendChild( item( 'Fit peak here', () => App.peaks.fitAt( energy, '' ) ) );
    menu.appendChild( item( 'Search energy ' + energy.toFixed( 1 ) + ' keV', () => App.search.setEnergy( energy ) ) );
  }

  menu.hidden = false;
  const rect = menu.getBoundingClientRect();
  const x = Math.min( pageX - window.scrollX, window.innerWidth - rect.width - 180 );
  const y = Math.min( pageY - window.scrollY, window.innerHeight - rect.height - 4 );
  menu.style.left = Math.max( 0, x ) + 'px';
  menu.style.top = Math.max( 0, y ) + 'px';
};


App.peaks.init = function() {
  document.getElementById( 'peaks-search' ).addEventListener( 'click', () => {
    if( App.hasForeground )
      App.run( 'searchPeaks', {}, { busy: true, busyText: 'Searching for peaks…' } );
  } );
  document.getElementById( 'peaks-csv' ).addEventListener( 'click', () => App.files.download( 'exportPeakCsv' ) );
  document.getElementById( 'peaks-clear' ).addEventListener( 'click', () => {
    if( App.peaks.list.length && window.confirm( 'Delete all ' + App.peaks.list.length + ' peaks?' ) )
      App.run( 'clearPeaks' );
  } );

  const menu = document.getElementById( 'ctx-menu' );
  document.addEventListener( 'mousedown', (e) => { if( !menu.contains( e.target ) ) menu.hidden = true; } );
  document.addEventListener( 'keydown', (e) => { if( e.key === 'Escape' ) menu.hidden = true; } );

  App.renderers.push( (state) => {
    if( state.peakList )
      App.peaks.render( state.peakList );
    else if( state.files && !state.files.FOREGROUND )
      App.peaks.render( [] );
  } );
};
