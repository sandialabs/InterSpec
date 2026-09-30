/* InterSpec Light - energy calibration: one editable coefficient set for the displayed foreground,
 * applied to all samples of the shown detectors; or fit from peaks with assigned sources; or
 * exported to, or applied from (drag-n-drop), a CALp file. */
'use strict';

App.ecal = {
  cal: null,
  fitFor: [ true, true ],
};

const MAX_CAL_COEFS = 5;
const POLY_LABELS = [ 'Offset (keV)', 'Gain (keV/ch)', 'Quadratic', 'Cubic', 'Quartic' ];
const FRF_LABELS = [ 'Offset (keV)', 'Gain (keV)', 'Quadratic', 'Cubic', 'Low-E term' ];


/** A link that downloads the applied calibration as a CALp file. */
const calpLink = () => App.el( 'a', { href: '#', class: 'calp-link', text: 'Export CALp',
  title: 'Download the applied energy calibration as a CALp file (InterSpec and PeakEasy read these);'
         + ' drop a CALp file on the page to apply one',
  onclick: (e) => { e.preventDefault(); App.files.download( 'exportCALp' ); } } );


/** Restores the calibrations the files were loaded with; enabled once any has changed. */
const revertButton = ( cal ) => App.el( 'button', { class: 'btn', text: 'Revert', disabled: !cal.changed,
  title: 'Restore the calibrations the files were loaded with',
  onclick: () => App.run( 'revertEnergyCal', {}, { busy: true } ) } );


App.ecal.render = function( cal ) {
  App.ecal.cal = cal;
  const container = document.getElementById( 'ecal' );
  container.innerHTML = '';

  if( !cal )
  {
    container.appendChild( App.el( 'p', { class: 'hint', text: 'No foreground spectrum.' } ) );
    return;
  }

  if( !cal.editable )
  {
    container.appendChild( App.el( 'p', { class: 'hint', text: 'The displayed spectrum uses a "' + cal.type
      + '" calibration, which can not be edited here.' } ) );
    container.appendChild( App.el( 'div', { class: 'row' }, [ revertButton( cal ), calpLink() ] ) );
    return;
  }

  const typeSel = App.el( 'select', {}, [
    App.el( 'option', { value: 'Polynomial', text: 'Polynomial' } ),
    App.el( 'option', { value: 'FullRangeFraction', text: 'Full Range Fraction' } )
  ] );
  typeSel.value = cal.type;
  container.appendChild( App.el( 'div', { class: 'row' }, [ App.el( 'span', { text: 'Type' } ), typeSel ] ) );

  const grid = App.el( 'div', { class: 'coefs' } );
  container.appendChild( grid );

  let coefs = cal.coefs.slice();
  while( coefs.length < 2 )
    coefs.push( 0 );

  const inputs = [];
  const fitCbs = [];
  const drawRows = () => {
    grid.innerHTML = '';
    inputs.length = 0;
    fitCbs.length = 0;
    grid.appendChild( App.el( 'span', { class: 'hint', text: 'Coefficient' } ) );
    grid.appendChild( App.el( 'span', { class: 'hint', text: 'Value' } ) );
    grid.appendChild( App.el( 'span', { class: 'hint', text: 'Fit', title: 'Fit this coefficient from peaks' } ) );
    const labels = (typeSel.value === 'FullRangeFraction') ? FRF_LABELS : POLY_LABELS;
    coefs.forEach( (c, i) => {
      // Coefficients are float32; show 7 significant digits, but keep the exact value if unedited
      const shown = String( Number( c.toPrecision( 7 ) ) );
      const input = App.el( 'input', { type: 'text', value: shown } );
      input.dataset.exact = String( c );
      input.dataset.shown = shown;
      input.addEventListener( 'keydown', (e) => { if( e.key === 'Enter' ) apply(); } );
      const cb = App.el( 'input', { type: 'checkbox' } );
      cb.checked = !!App.ecal.fitFor[i];
      cb.addEventListener( 'change', () => { App.ecal.fitFor[i] = cb.checked; } );
      inputs.push( input );
      fitCbs.push( cb );
      grid.appendChild( App.el( 'span', { text: labels[i] || ('Coef ' + i) } ) );
      grid.appendChild( input );
      grid.appendChild( cb );
    } );
  };

  const readCoefs = () => inputs.map( (inp) => {
    if( inp.value === inp.dataset.shown )
      return parseFloat( inp.dataset.exact );
    const v = parseFloat( inp.value );
    if( !isFinite( v ) )
      throw new Error( 'Invalid coefficient "' + inp.value + '"' );
    return v;
  } );

  const apply = () => {
    try
    {
      App.run( 'setEnergyCal', { type: typeSel.value, coefs: readCoefs() }, { busy: true } );
    }catch( e )
    {
      App.toast( e.message, true );
    }
  };

  typeSel.addEventListener( 'change', () => {
    // Re-express the same calibration in the other form; not applied until "Apply"
    try
    {
      const from = (typeSel.value === 'Polynomial') ? 'FullRangeFraction' : 'Polynomial';
      coefs = App.call( 'convertCalCoefs', { from: from, to: typeSel.value, coefs: readCoefs() } ).coefs;
      drawRows();
    }catch( e )
    {
      App.toast( e.message, true );
    }
  } );

  drawRows();

  const addBtn = App.el( 'button', { class: 'btn-sm', text: '+ coef', title: 'Add a higher order coefficient',
    onclick: () => { try{ coefs = readCoefs(); }catch( e ){} if( coefs.length < MAX_CAL_COEFS ){ coefs.push( 0 ); drawRows(); } } } );
  const rmBtn = App.el( 'button', { class: 'btn-sm', text: '− coef', title: 'Remove the highest order coefficient',
    onclick: () => { try{ coefs = readCoefs(); }catch( e ){} if( coefs.length > 2 ){ coefs.pop(); drawRows(); } } } );

  const applyBtn = App.el( 'button', { class: 'btn primary', text: 'Apply', onclick: apply } );
  const fitBtn = App.el( 'button', { class: 'btn', text: 'Fit from peaks',
    title: 'Fit the checked coefficients using peaks with an assigned source, and "Cal" checked',
    onclick: () => {
      const fitFor = fitCbs.map( (cb) => cb.checked );
      App.run( 'fitEnergyCal', { fitFor: fitFor }, { busy: true } );
    } } );
  const revertBtn = revertButton( cal );

  container.appendChild( App.el( 'div', { class: 'row' }, [ addBtn, rmBtn, calpLink() ] ) );
  container.appendChild( App.el( 'div', { class: 'row' }, [ applyBtn, fitBtn, revertBtn ] ) );

  const notes = [ cal.numChannels + ' channels, ' + App.fmt( cal.lowerEnergy, 1 ) + '–' + App.fmt( cal.upperEnergy, 1 ) + ' keV' ];
  if( cal.numDevPairs )
    notes.push( cal.numDevPairs + ' deviation pairs (kept unchanged)' );
  notes.push( 'Changes apply to all samples of the shown detectors of the foreground file.' );
  container.appendChild( App.el( 'p', { class: 'hint', text: notes.join( '. ' ) } ) );
};


/** If the bytes are a CALp file: its first non-blank line contains "CALp File" (as SpecUtils checks). */
App.ecal.isCALp = function( bytes ) {
  return /^\s*[^\r\n]*calp file/i.test( new TextDecoder().decode( bytes.subarray( 0, 256 ) ) );
};


/** Which spectra a CALp file should be applied to: just the foreground, unless a background or
 secondary from another file, with as many channels, is loaded, in which case a dialog asks
 (resolves null if cancelled). */
App.ecal.chooseCALpTargets = function( name ) {
  const slots = App.files.slots || {};
  const nchannel = slots.FOREGROUND ? slots.FOREGROUND.numChannels : undefined;
  const choices = [], fileIds = [];
  for( const type of SPEC_TYPES )
  {
    const info = slots[type];
    if( info && !fileIds.includes( info.fileId ) && (info.numChannels === nchannel) )
    {
      fileIds.push( info.fileId );
      choices.push( { type: type, cb: App.el( 'input', { type: 'checkbox', checked: 'checked' } ), name: info.name } );
    }
  }
  if( choices.length < 2 )
    return Promise.resolve( [ 'FOREGROUND' ] );

  return new Promise( (resolve) => {
    let answer = null;
    const applyBtn = App.el( 'button', { class: 'btn primary', text: 'Apply', autofocus: 'autofocus', onclick: () => {
      answer = choices.filter( (c) => c.cb.checked ).map( (c) => c.type );
      dlg.close();
    } } );
    for( const c of choices )
      c.cb.addEventListener( 'change', () => { applyBtn.disabled = !choices.some( (c) => c.cb.checked ); } );

    const dlg = App.el( 'dialog', { class: 'dlg' }, [
      App.el( 'p', {}, [ 'Apply the energy calibration in ', App.el( 'b', { text: name } ), ' to:' ] ),
      ...choices.map( (c) => App.el( 'label', { title: c.name }, [ c.cb, ' ' + SPEC_LABELS[c.type] + ' (' + c.name + ')' ] ) ),
      App.el( 'div', { class: 'row dlg-btns' }, [
        App.el( 'button', { class: 'btn', text: 'Cancel', onclick: () => dlg.close() } ),
        applyBtn
      ] )
    ] );
    // Escape is handled here, since SpectrumChartD3 cancels its default (closing the dialog) at the window
    dlg.addEventListener( 'keydown', (e) => {
      if( e.key === 'Escape' )
      {
        e.stopPropagation();
        dlg.close();
      }
    } );
    dlg.addEventListener( 'close', () => {
      dlg.remove();
      resolve( (answer && answer.length) ? answer : null );
    } );
    document.body.appendChild( dlg );
    dlg.showModal();
  } );
};


App.ecal.init = function() {
  App.renderers.push( (state) => {
    if( 'energyCal' in state )
      App.ecal.render( state.energyCal );
  } );
};
