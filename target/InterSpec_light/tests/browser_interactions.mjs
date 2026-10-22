// More browser checks of dist/InterSpecLight.html: drag-n-drop zones, NaI fitting, ctrl-drag ROI
//  fit, shift-drag erase, ROI-edge dragging, and assigning peak sources.
//  Usage: node tests/browser_interactions.mjs <out_dir> [page.html]   (page defaults to dist/InterSpecLight.html)

import fs from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

const here = path.dirname( fileURLToPath( import.meta.url ) );
const repo = path.resolve( here, '..', '..', '..' );
const require = createRequire( import.meta.url );
const { chromium } = require( path.join( repo, 'target/testing/PlaywrightPhoneEmulation/node_modules/playwright' ) );

const outDir = process.argv[2] || '/tmp';
const failures = [];
const check = ( cond, msg ) => { console.log( (cond ? 'PASS ' : 'FAIL ') + msg ); if( !cond ) failures.push( msg ); };

// Chrome: CHROME_PATH, else the standard macOS install, else Playwright's own browser
const macChrome = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const chromePath = process.env.CHROME_PATH || (fs.existsSync( macChrome ) ? macChrome : undefined);
const browser = await chromium.launch( { headless: true, executablePath: chromePath } );
const page = await browser.newPage( { viewport: { width: 1400, height: 850 } } );
const errors = [];
page.on( 'console', (m) => { if( m.type() === 'error' ) errors.push( m.text() ); } );
page.on( 'pageerror', (e) => errors.push( String( e ) ) );

await page.goto( 'file://' + (process.argv[3] ? path.resolve( process.argv[3] ) : path.join( here, '..', 'dist', 'InterSpecLight.html' )) );
await page.waitForFunction( () => window.App && App.lightCall && App.spectrumChart );

/** Drags a file over the page, and drops it on the zone for `type`, like a user would. */
const dropFile = async ( filePath, type ) => {
  const b64 = fs.readFileSync( filePath ).toString( 'base64' );
  await page.evaluate( ({ b64, name, type }) => {
    const bytes = Uint8Array.from( atob( b64 ), (c) => c.charCodeAt( 0 ) );
    const dt = new DataTransfer();
    dt.items.add( new File( [ bytes ], name ) );
    document.dispatchEvent( new DragEvent( 'dragenter', { dataTransfer: dt, bubbles: true } ) );
    const zone = document.querySelector( `.drop-zone[data-type="${type}"]` );
    zone.dispatchEvent( new DragEvent( 'dragover', { dataTransfer: dt, bubbles: true } ) );
    zone.dispatchEvent( new DragEvent( 'drop', { dataTransfer: dt, bubbles: true } ) );
  }, { b64, name: path.basename( filePath ), type } );
};

const pointFor = async ( energy, yfrac ) => page.evaluate( ({ e, yfrac }) => {
  const c = App.spectrumChart;
  const r = c.plot.node().getBoundingClientRect();
  return { x: r.left + c.xScale( e ), y: r.top + (yfrac || 0.5)*r.height };
}, { e: energy, yfrac } );
const zoom = async ( lo, hi ) => page.evaluate( ({ lo, hi }) => { App.spectrumChart.setXAxisRange( lo, hi, false, true ); App.spectrumChart.redraw()(); }, { lo, hi } );

// Drop a NaI Ba133 spectrum as the foreground
await dropFile( path.join( repo, 'target/testing/test_data/PeakFitLM/Ba133_Unshielded.n42' ), 'FOREGROUND' );
await page.waitForFunction( () => App.hasForeground, null, { timeout: 30000 } );
check( await page.evaluate( () => App.files.slots.FOREGROUND.name ) === 'Ba133_Unshielded.n42', 'dropped foreground' );
// Builds with the peak-file option load the 7 peaks InterSpec saved in this file; start without them
const nFilePeaks = await page.evaluate( () => App.peaks.list.length );
check( nFilePeaks === (await page.evaluate( () => App.features.peakFileIo ? 7 : 0 )), `peaks loaded from the file: ${nFilePeaks}` );
await page.evaluate( () => App.run( 'clearPeaks' ) );
check( await page.evaluate( () => App.isHighRes ) === false, 'NaI classified as not high-res' );

// Drop a background onto the background zone
await dropFile( path.join( repo, 'example_spectra/background_20100317.n42' ), 'BACKGROUND' );
await page.waitForFunction( () => !!App.files.slots.BACKGROUND, null, { timeout: 30000 } );
check( await page.evaluate( () => App.files.slots.BACKGROUND.name ) === 'background_20100317.n42', 'dropped background' );
check( await page.evaluate( () => document.querySelectorAll( '#spectrum-chart path.speclinepath' ).length ) >= 2, 'foreground and background drawn' );

// Drop CALp files.  The background has 16384 channels (the foreground 1024), so it is not offered,
//  and the calibration goes straight to the foreground.
const calpFile = path.join( here, 'data', 'peakeasy_dev_pairs.CALp' );
const specCal = ( type ) => page.evaluate( (type) => {
  const s = App.call( 'getState' ).spectra.find( (s) => s.type === type );
  return JSON.stringify( s.xeqn || s.x.slice( 0, 3 ) );
}, type );
const backBefore = await specCal( 'BACKGROUND' );
await dropFile( calpFile, 'BACKGROUND' );
await page.waitForFunction( () => App.ecal.cal && (App.ecal.cal.numDevPairs === 5), null, { timeout: 10000 } );
check( !(await page.evaluate( () => document.querySelector( 'dialog' ) )) && ((await specCal( 'BACKGROUND' )) === backBefore),
       'CALp applied to the foreground only, without asking, as the background has more channels' );
await page.evaluate( () => App.run( 'revertEnergyCal' ) );

// A secondary with as many channels (the same file, loaded again) is offered in a dialog
await dropFile( path.join( repo, 'target/testing/test_data/PeakFitLM/Ba133_Unshielded.n42' ), 'SECONDARY' );
await page.waitForFunction( () => !!App.files.slots.SECONDARY, null, { timeout: 30000 } );
const secondBefore = await specCal( 'SECONDARY' );
await dropFile( calpFile, 'FOREGROUND' );
await page.waitForSelector( 'dialog.dlg[open]' );
check( await page.evaluate( () => Array.from( document.querySelectorAll( 'dialog.dlg label' ) ).map( (l) => l.textContent.trim().split( ' ' )[0] + ':' + l.querySelector( 'input' ).checked ).join() )
       === 'Foreground:true,Secondary:true', 'CALp dialog offers the foreground and secondary, not the background' );
await page.screenshot( { path: path.join( outDir, 'inter_0_calp_dialog.png' ) } );
await page.locator( 'dialog.dlg label', { hasText: 'Foreground' } ).locator( 'input' ).uncheck();
await page.locator( 'dialog.dlg label', { hasText: 'Secondary' } ).locator( 'input' ).uncheck();
check( await page.evaluate( () => document.querySelector( 'dialog.dlg button.primary' ).disabled ), 'Apply is disabled with nothing checked' );
await page.locator( 'dialog.dlg label', { hasText: 'Foreground' } ).locator( 'input' ).check();
await page.click( 'dialog.dlg button.primary' );
await page.waitForFunction( () => App.ecal.cal && (App.ecal.cal.numDevPairs === 5), null, { timeout: 10000 } );
check( (await specCal( 'SECONDARY' )) === secondBefore, 'CALp applied to the foreground only, as chosen' );
await page.evaluate( () => App.run( 'revertEnergyCal' ) );
await dropFile( calpFile, 'FOREGROUND' );
await page.waitForSelector( 'dialog.dlg[open]' );
await page.keyboard.press( 'Escape' );
await page.waitForSelector( 'dialog.dlg', { state: 'detached' } );
await page.waitForTimeout( 200 );
check( await page.evaluate( () => !App.ecal.cal.changed ), 'cancelling the CALp dialog changes nothing' );
await page.evaluate( () => App.run( 'unload', { type: 'SECONDARY' } ) );

// A lower-channel-energy calibration can not be edited, but can be reverted
await dropFile( path.join( here, 'data', 'exact_energies_1024.CALp' ), 'FOREGROUND' );
await page.waitForFunction( () => App.ecal.cal && (App.ecal.cal.type === 'LowerChannelEdge'), null, { timeout: 10000 } );
await page.evaluate( () => { document.getElementById( 'ecal-section' ).open = true; } );
await page.locator( '#ecal button', { hasText: 'Revert' } ).click();
check( await page.waitForFunction( () => App.ecal.cal && (App.ecal.cal.type === 'Polynomial') && !App.ecal.cal.changed,
                                   null, { timeout: 10000 } ).then( () => true, () => false ),
       'Revert restores a calibration replaced by an "Exact Energies" CALp' );
await page.evaluate( () => { document.getElementById( 'ecal-section' ).open = false; } );

// Double-click the NaI 356 keV peak
await page.evaluate( () => App.ref.add( 'Ba133' ) );
await zoom( 150, 500 );
let pt = await pointFor( 356 );
await page.mouse.dblclick( pt.x, pt.y );
await page.waitForFunction( () => App.peaks.list.length >= 1, null, { timeout: 30000 } );
let peaks = await page.evaluate( () => App.peaks.list );
const p356 = peaks.find( (p) => Math.abs( p.mean - 356 ) < 15 );
check( !!p356 && (p356.fwhm > 15) && (p356.fwhm < 45), `NaI 356 peak: ${p356 ? p356.mean.toFixed(1) + ' keV, FWHM ' + p356.fwhm.toFixed(1) : 'none'}` );
check( !!p356 && /Ba133/.test( p356.source ), 'NaI peak source: ' + (p356 ? p356.source : '') );
await page.screenshot( { path: path.join( outDir, 'inter_1_nai.png' ) } );

// Ctrl-drag over the 60-110 keV region to fit peaks there
await zoom( 20, 160 );
const a = await pointFor( 65 ), b = await pointFor( 105 );
const nbefore = peaks.length;
await page.keyboard.down( 'Control' );
await page.mouse.move( a.x, a.y );
await page.mouse.down();
await page.mouse.move( b.x, b.y, { steps: 10 } );
await page.waitForTimeout( 300 );
await page.mouse.up();
await page.keyboard.up( 'Control' );
await page.waitForFunction( (n) => App.peaks.list.length > n, nbefore, { timeout: 30000 } ).catch( () => {} );
peaks = await page.evaluate( () => App.peaks.list );
check( peaks.length > nbefore, `ctrl-drag added ${peaks.length - nbefore} peak(s): ` + peaks.map( (p) => p.mean.toFixed( 1 ) ).join( ', ' ) );
await page.screenshot( { path: path.join( outDir, 'inter_2_ctrldrag.png' ) } );

// Drag the upper edge of the 356 keV ROI with the mouse, then shift-drag erase it
await zoom( 150, 500 );
const roi = await page.evaluate( () => App.peaks.list.find( (p) => Math.abs( p.mean - 356 ) < 15 ) );
const edge = await pointFor( roi.upper );
let handle = null;
for( let yfrac = 0.95; (yfrac > 0.05) && !handle; yfrac -= 0.05 )
{
  const pt2 = await pointFor( roi.upper, yfrac );
  await page.mouse.move( pt2.x - 1, pt2.y );
  handle = await page.evaluate( () => {
    const boxes = Array.from( document.querySelectorAll( '#spectrum-chart .roiDragBox' ) )
                       .filter( (b) => (getComputedStyle( b ).visibility !== 'hidden') && b.getBoundingClientRect().width );
    if( !boxes.length )
      return null;
    const r = boxes[boxes.length - 1].getBoundingClientRect();
    return { x: r.left + 0.5*r.width, y: r.top + 0.5*r.height };
  } );
}
check( !!handle, 'ROI drag handle shown when hovering the ROI edge' );
if( handle )
{
  const shift = (await pointFor( roi.upper + 20 )).x - edge.x;
  await page.mouse.move( handle.x, handle.y );
  await page.mouse.down();
  await page.mouse.move( handle.x + 0.5*shift, handle.y, { steps: 5 } );
  await page.waitForTimeout( 700 );  //chart throttles drag requests to 500 ms
  await page.mouse.move( handle.x + shift, handle.y, { steps: 5 } );
  await page.waitForTimeout( 700 );
  const preview = await page.evaluate( () => App.spectrumChart.roiBeingDrugUpdate );
  check( !!preview && (preview.upperEnergy > roi.upper + 5), 'ROI drag preview returned while dragging' );
  await page.mouse.up();
  await page.waitForTimeout( 500 );
}
const dragged = await page.evaluate( () => App.peaks.list.find( (p) => Math.abs( p.mean - 356 ) < 15 ) );
check( !!dragged && (dragged.upper > roi.upper + 5), `ROI upper edge moved ${roi.upper.toFixed(1)} -> ${dragged ? dragged.upper.toFixed(1) : '?'}` );

await zoom( 150, 500 );
const e0 = await pointFor( 340 ), e1 = await pointFor( 372 );
const nWithRoi = await page.evaluate( () => App.peaks.list.length );
await page.keyboard.down( 'Shift' );
await page.mouse.move( e0.x, e0.y );
await page.mouse.down();
await page.mouse.move( e1.x, e1.y, { steps: 8 } );
await page.mouse.up();
await page.keyboard.up( 'Shift' );
await page.waitForFunction( (n) => App.peaks.list.length < n, nWithRoi, { timeout: 30000 } ).catch( () => {} );
check( await page.evaluate( () => !App.peaks.list.some( (p) => Math.abs( p.mean - 356 ) < 10 ) ), 'shift-drag erased the 356 keV peak' );

// Peak sources: fit without reference lines shown, then type a source into the peak table
await page.evaluate( () => { App.run( 'clearPeaks' ); App.ref.remove( 'Ba133' ); } );
pt = await pointFor( 356 );
await page.mouse.dblclick( pt.x, pt.y );
await page.waitForFunction( () => App.peaks.list.length === 1, null, { timeout: 30000 } );
const sourceOf = () => page.evaluate( () => App.peaks.list[0].source );
check( await sourceOf() === '', 'peak fit without reference lines has no source' );
await page.click( '#peak-table td.source' );
await page.waitForSelector( '#peak-table input.source-edit' );
await page.keyboard.type( 'Ba133' );
await page.screenshot( { path: path.join( outDir, 'inter_3_type_source.png' ) } );
await page.keyboard.press( 'Enter' );
await page.waitForFunction( () => App.peaks.list[0].source !== '', null, { timeout: 10000 } ).catch( () => {} );
check( await sourceOf() === 'Ba133 356.02 keV', 'typed source: ' + await sourceOf() );
check( (await page.textContent( '#peak-table td.source' )).includes( 'Ba133 356.02' ), 'peak table shows the typed source' );

// Escape cancels an edit; an unknown source is an error, and leaves the source as it was
await page.click( '#peak-table td.source .ellipsis' );
await page.keyboard.type( 'Cs137' );
await page.keyboard.press( 'Escape' );
check( (await sourceOf() === 'Ba133 356.02 keV') && !(await page.$( '#peak-table input.source-edit' )), 'Escape cancels a source edit' );
await page.click( '#peak-table td.source .ellipsis' );
await page.keyboard.type( 'Xx123' );
await page.keyboard.press( 'Enter' );
await page.waitForSelector( '#toast.error:not([hidden])', { timeout: 5000 } ).catch( () => {} );
check( /not a source/.test( await page.textContent( '#toast' ) ), 'unknown source shows an error: ' + await page.textContent( '#toast' ) );
check( (await sourceOf() === 'Ba133 356.02 keV') && (await page.textContent( '#peak-table td.source' )).includes( 'Ba133 356.02' ),
       'unknown source leaves the source unchanged' );

// With reference lines shown, the right-click menu offers their sources near the peak
await page.evaluate( () => App.ref.add( 'Background' ) );
await page.mouse.click( pt.x, pt.y, { button: 'right' } );
await page.waitForSelector( '#ctx-menu:not([hidden])' );
const assignItems = await page.$$eval( '#ctx-menu .item', (els) => els.map( (e) => e.textContent ).filter( (t) => t.startsWith( 'Assign as' ) ) );
check( assignItems.includes( 'Assign as Ra226 351.93 keV' ), 'right-click menu: ' + assignItems.join( ', ' ) );
await page.screenshot( { path: path.join( outDir, 'inter_4_assign_menu.png' ) } );
await page.click( '#ctx-menu .item:text-is("Assign as Ra226 351.93 keV")' );
await page.waitForFunction( () => /Ra226/.test( App.peaks.list[0].source ), null, { timeout: 10000 } ).catch( () => {} );
const assigned = await page.evaluate( () => App.peaks.list[0] );
check( (assigned.source === 'Ra226 351.93 keV') && (assigned.color === '#967f55'), `assigned from the menu: ${assigned.source}, ${assigned.color}` );
await page.mouse.click( pt.x, pt.y, { button: 'right' } );
await page.waitForSelector( '#ctx-menu:not([hidden])' );
check( !(await page.textContent( '#ctx-menu' )).includes( 'Assign as Ra226' ), 'menu leaves out the nuclide the peak already has' );
await page.keyboard.press( 'Escape' );

// Skew change on a remaining peak via the context menu path
const any = await page.evaluate( () => App.peaks.list[0] );
if( any )
{
  const res = await page.evaluate( (m) => App.run( 'setSkewType', { energy: m, type: 2 } ), any.mean );
  check( !!res, 'skew type change ran' );
}

await page.screenshot( { path: path.join( outDir, 'inter_5_final.png' ) } );
const unexpected = errors.filter( (e) => !e.includes( 'Xx123' ) );  //App.run logs the error shown for the unknown source
check( unexpected.length === 0, 'no unexpected console errors' + (unexpected.length ? (': ' + unexpected.join( ' | ' )) : '') );
await browser.close();
console.log( failures.length ? `${failures.length} FAILURES` : 'ALL PASSED' );
process.exit( failures.length ? 1 : 0 );
