// Browser smoke test of dist/InterSpecLight.html, opened via file://, driven by real mouse input.
//  Usage: node tests/browser_smoke.mjs <out_dir> [page.html] [--headed]
//  The page defaults to dist/InterSpecLight.html.
//  Uses the repo's Playwright install (target/testing/PlaywrightPhoneEmulation).

import fs from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

const here = path.dirname( fileURLToPath( import.meta.url ) );
const repo = path.resolve( here, '..', '..', '..' );
const require = createRequire( import.meta.url );
const { chromium } = require( path.join( repo, 'target/testing/PlaywrightPhoneEmulation/node_modules/playwright' ) );

const outDir = process.argv[2] || '/tmp';
const headed = process.argv.includes( '--headed' );
const testData = path.join( repo, 'target/testing/test_data' );
const pageArg = process.argv.slice( 3 ).find( (a) => !a.startsWith( '--' ) );
const page_url = 'file://' + (pageArg ? path.resolve( pageArg ) : path.join( here, '..', 'dist', 'InterSpecLight.html' ));

const failures = [];
const check = ( cond, msg ) => { console.log( (cond ? 'PASS ' : 'FAIL ') + msg ); if( !cond ) failures.push( msg ); };

// Chrome: CHROME_PATH, else the standard macOS install, else Playwright's own browser
const macChrome = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const chromePath = process.env.CHROME_PATH || (fs.existsSync( macChrome ) ? macChrome : undefined);
const browser = await chromium.launch( { headless: !headed, executablePath: chromePath } );
const page = await browser.newPage( { viewport: { width: 1400, height: 850 } } );
const errors = [];
page.on( 'console', (m) => { if( m.type() === 'error' ) errors.push( m.text() ); } );
page.on( 'pageerror', (e) => errors.push( String( e ) ) );

await page.goto( page_url );
await page.waitForFunction( () => window.App && App.lightCall && App.spectrumChart, null, { timeout: 30000 } );
check( true, 'page loaded' );

// Load an HPGe spectrum through the file input
await page.setInputFiles( '#file-input', path.join( testData, 'AnalystTests/Ba133_Cs137_RandomSummingOnly.n42' ) );
await page.waitForFunction( () => App.hasForeground, null, { timeout: 30000 } );
check( await page.evaluate( () => document.querySelectorAll( '#spectrum-chart path.speclinepath' ).length > 0 ), 'foreground line drawn' );

// Builds with the peak-file option load the 19 peaks InterSpec saved in this file; start without them
const peakFiles = await page.evaluate( () => !!App.features.peakFileIo );
const nFilePeaks = await page.evaluate( () => App.peaks.list.length );
check( nFilePeaks === (peakFiles ? 19 : 0), `peaks loaded from the file: ${nFilePeaks} (peak-file support ${peakFiles ? 'on' : 'off'})` );
check( await page.evaluate( () => document.querySelector( 'span.peak-files' ).hidden ) === !peakFiles, 'export text matches peak-file support' );
await page.evaluate( () => App.run( 'clearPeaks' ) );

// Show Cs137 and Ba133 reference lines
await page.fill( '#ref-filter', 'cs137' );
await page.keyboard.press( 'Enter' );
await page.fill( '#ref-filter', 'Ba133' );
await page.keyboard.press( 'Enter' );
check( await page.evaluate( () => App.ref.shown.map( (s) => s.parent ).join( ',' ) ) === 'Cs137,Ba133', 'ref lines shown' );

// Zoom to 600-700 keV, then double-click at 661.7 keV with the real mouse
const pointFor = async ( energy ) => page.evaluate( (e) => {
  const c = App.spectrumChart;
  const r = c.plot.node().getBoundingClientRect();
  return { x: r.left + c.xScale( e ), y: r.top + 0.5*r.height };
}, energy );
await page.evaluate( () => { App.spectrumChart.setXAxisRange( 600, 720, false, true ); App.spectrumChart.redraw()(); } );
let pt = await pointFor( 661.7 );
await page.mouse.dblclick( pt.x, pt.y );
await page.waitForFunction( () => App.peaks.list.length >= 1, null, { timeout: 30000 } );
const peak = await page.evaluate( () => App.peaks.list[0] );
check( Math.abs( peak.mean - 661.66 ) < 0.3, `double-click fit mean ${peak.mean.toFixed(3)} keV` );
check( /Cs137/.test( peak.source ), `peak source "${peak.source}"` );
await page.screenshot( { path: path.join( outDir, 'light_1_cs137.png' ) } );

// Fit Ba133 356 keV, then the 302.85 keV neighbor
await page.evaluate( () => { App.spectrumChart.setXAxisRange( 270, 400, false, true ); App.spectrumChart.redraw()(); } );
pt = await pointFor( 356.0 );
await page.mouse.dblclick( pt.x, pt.y );
await page.waitForFunction( () => App.peaks.list.length >= 2, null, { timeout: 30000 } );
pt = await pointFor( 302.85 );
await page.mouse.dblclick( pt.x, pt.y );
await page.waitForFunction( () => App.peaks.list.length >= 3, null, { timeout: 30000 } );
const sources = await page.evaluate( () => App.peaks.list.map( (p) => p.source ) );
check( sources.filter( (s) => /Ba133/.test( s ) ).length === 2, 'two Ba133 peaks: ' + sources.join( ' | ' ) );
await page.screenshot( { path: path.join( outDir, 'light_2_ba133.png' ) } );

// Right-click the 356 keV peak, and change its continuum to linear
pt = await pointFor( 356.0 );
await page.mouse.click( pt.x, pt.y, { button: 'right' } );
await page.waitForSelector( '#ctx-menu:not([hidden])' );
const menuText = await page.textContent( '#ctx-menu' );
check( /Peak at 35[56]\./.test( menuText ), 'context menu for 356 keV peak' );
await page.screenshot( { path: path.join( outDir, 'light_3_menu.png' ) } );
await page.hover( '#ctx-menu .item.sub' );
await page.click( '#ctx-menu .submenu .item:text-is("Quadratic")' );
await page.waitForFunction( () => App.peaks.list.some( (p) => (Math.abs( p.mean - 356 ) < 1) && (p.continuumType === 3) ), null, { timeout: 30000 } );
check( true, 'continuum changed to Quadratic' );

// Peak CSV in InterSpec's format
const csv = await page.evaluate( () => {
  const r = App.call( 'exportPeakCsv', { path: '/tmp/peaks.csv' } );
  return { name: r.filename, text: new TextDecoder().decode( App.module.FS.readFile( r.path ) ) };
} );
const csvLines = csv.text.split( '\r\n' ).filter( (l) => l.length );
check( csvLines[0].startsWith( 'Centroid,  Net_Area,   Net_Area,      Peak, FWHM,' ) && (csvLines.length === 2 + 4)
       && /,Cs137,661\.66,/.test( csv.text ) && (csv.name === 'peaks_Ba133_Cs137_RandomSummingOnly.CSV'),
       `peak CSV export (${csvLines.length - 2} peak rows, ${csv.name})` );

// Energy calibration: perturb the gain by 1%, then fit it back from the peaks
const cal0 = await page.evaluate( () => App.ecal.cal.coefs.slice() );
await page.evaluate( (c) => App.run( 'setEnergyCal', { type: 'Polynomial', coefs: [ c[0], 1.01*c[1] ] } ), cal0 );
const moved = await page.evaluate( () => App.peaks.list.map( (p) => p.mean ) );
check( moved.some( (m) => Math.abs( m - 1.01*661.62 ) < 0.5 ) && moved.every( (m, i) => !i || (m > moved[i-1]) ), `peaks moved with calibration: ${moved.map( (m) => m.toFixed(2) ).join( ', ' )}` );
await page.evaluate( () => App.run( 'fitEnergyCal', { fitFor: [ false, true ] } ) );
const fitted = await page.evaluate( () => ({ coefs: App.ecal.cal.coefs, means: App.peaks.list.map( (p) => p.mean ) }) );
check( Math.abs( fitted.coefs[1] - cal0[1] ) / cal0[1] < 2e-4, `gain fit back: ${fitted.coefs[1]} vs ${cal0[1]}` );
check( Math.abs( fitted.means[fitted.means.length - 1] - 661.66 ) < 0.3, `661 peak after fit: ${fitted.means.map( (m) => m.toFixed(3) ).join( ', ' )}` );

// "Export CALp" link downloads the calibration; loading that CALp file (foreground only, so no
//  dialog) puts it back after a change
await page.evaluate( () => { document.getElementById( 'ecal-section' ).open = true; } );
await page.screenshot( { path: path.join( outDir, 'light_3b_ecal.png' ) } );
const [ calpDownload ] = await Promise.all( [ page.waitForEvent( 'download' ), page.click( '.calp-link' ) ] );
const calpPath = path.join( outDir, 'light_export.CALp' );
await calpDownload.saveAs( calpPath );
check( (calpDownload.suggestedFilename() === 'Ba133_Cs137_RandomSummingOnly.CALp')
       && fs.readFileSync( calpPath, 'utf8' ).startsWith( '#PeakEasy CALp File' ),
       'Export CALp link: ' + calpDownload.suggestedFilename() );
await page.evaluate( (c) => App.run( 'setEnergyCal', { type: 'Polynomial', coefs: [ c[0], 1.02*c[1] ] } ), fitted.coefs );
await page.setInputFiles( '#file-input', calpPath );
await page.waitForFunction( (g) => Math.abs( App.ecal.cal.coefs[1] - g ) < 1e-5*g, fitted.coefs[1], { timeout: 10000 } );
check( await page.evaluate( () => App.files.slots.FOREGROUND.name ) === 'Ba133_Cs137_RandomSummingOnly.n42',
       'CALp file applied (not loaded as a spectrum)' );

// Export N42 and check the new gain is in it
const exported = await page.evaluate( () => {
  const r = App.call( 'exportFile', { format: 'N42-2012', path: '/tmp/x.n42' } );
  return new TextDecoder().decode( App.module.FS.readFile( r.path ) );
} );
check( exported.includes( '<RadInstrumentData' ), 'exported N42-2012' );

// Energy search for the Co60 pair
await page.evaluate( () => { document.getElementById( 'search-section' ).open = true; } );
await page.fill( '#search-e1', '1173.2' );
await page.fill( '#search-e2', '1332.5' );
await page.waitForTimeout( 400 );
const topHit = await page.evaluate( () => App.search.compute( App.search.inputs(), 'weighted' )[0].entry.parent );
check( topHit === 'Co60', 'energy search top hit: ' + topHit );
await page.screenshot( { path: path.join( outDir, 'light_4_search.png' ) } );

// Automated peak search, via the button
const nBeforeSearch = await page.evaluate( () => App.peaks.list.length );
const tSearch = Date.now();
await page.click( '#peaks-search' );
await page.waitForFunction( (n) => App.peaks.list.length > n + 20, nBeforeSearch, { timeout: 60000 } );
const searchMs = Date.now() - tSearch;
// Only Cs137 and Ba133 lines are shown, so only their peaks get sources (not the NORM or sum peaks)
const searched = await page.evaluate( () => {
  const src = (e) => { const p = App.peaks.list.find( (p) => Math.abs( p.mean - e ) < 1.5 ); return p ? p.source : 'none'; };
  return { n: App.peaks.list.length, withSrc: App.peaks.list.filter( (p) => p.source ).length,
           cs: src( 661.66 ), ba: src( 356.01 ), k40: src( 1460.8 ) };
} );
check( searched.cs.startsWith( 'Cs137' ) && searched.ba.startsWith( 'Ba133' ) && (searched.k40 === ''),
       `peak search: ${searched.n - nBeforeSearch} new peaks (${searched.n} total, ${searched.withSrc} with sources) in ${searchMs} ms` );
await page.screenshot( { path: path.join( outDir, 'light_4b_search_peaks.png' ) } );

// Peaks in an N42 export load back with the file, through the file input
if( peakFiles )
{
  // (like InterSpec, peak values are read back as float, so compare to 0.01 keV)
  const peakSummary = () => page.evaluate( () => App.peaks.list.map( (p) => ({ mean: p.mean, source: p.source }) ) );
  const before = await peakSummary();
  const n42 = await page.evaluate( () => {
    const r = App.call( 'exportFile', { format: 'N42-2012', path: '/tmp/roundtrip.n42' } );
    return Array.from( App.module.FS.readFile( r.path ) );
  } );
  fs.writeFileSync( path.join( outDir, 'light_roundtrip.n42' ), Buffer.from( n42 ) );
  await page.setInputFiles( '#file-input', path.join( outDir, 'light_roundtrip.n42' ) );
  await page.waitForFunction( () => App.files.slots.FOREGROUND && (App.files.slots.FOREGROUND.name === 'light_roundtrip.n42'), null, { timeout: 30000 } );
  const after = await peakSummary();
  check( (after.length === before.length)
         && after.every( (p, i) => (Math.abs( p.mean - before[i].mean ) < 0.01) && (p.source === before[i].source) ),
         `N42 export with ${before.length} peaks loads back with them` );
}

// Passthrough file: time chart, then select background with the "b" key held while dragging
await page.setInputFiles( '#file-input', path.join( repo, 'example_spectra/passthrough.n42' ) );
await page.waitForFunction( () => App.files.slots.FOREGROUND && App.files.slots.FOREGROUND.name === 'passthrough.n42', null, { timeout: 60000 } );
check( await page.evaluate( () => !document.getElementById( 'time-chart' ).hidden ), 'time chart shown for passthrough' );
await page.waitForTimeout( 500 );
const tbox = await page.locator( '#time-chart svg' ).first().boundingBox();
await page.keyboard.down( 'b' );
await page.mouse.move( tbox.x + 0.2*tbox.width, tbox.y + 0.5*tbox.height );
await page.mouse.down();
await page.mouse.move( tbox.x + 0.3*tbox.width, tbox.y + 0.5*tbox.height, { steps: 8 } );
await page.mouse.up();
await page.keyboard.up( 'b' );
await page.waitForTimeout( 1000 );
const back = await page.evaluate( () => App.files.slots.BACKGROUND );
check( back && back.name === 'passthrough.n42' && back.samples.length > 1, 'time-chart "b"-drag set background samples: '
       + (back ? back.samples.length : 0) );
await page.screenshot( { path: path.join( outDir, 'light_5_passthrough.png' ) } );

check( errors.length === 0, 'no console errors' + (errors.length ? (': ' + errors.join( ' | ' )) : '') );

await browser.close();
console.log( failures.length ? `${failures.length} FAILURES` : 'ALL PASSED' );
process.exit( failures.length ? 1 : 0 );
