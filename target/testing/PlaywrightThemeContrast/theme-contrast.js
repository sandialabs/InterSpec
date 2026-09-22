// Measures WCAG text/border contrast of the running InterSpec app in the light and dark colour
// themes, and fails when anything is below AA (4.5:1 for text, 3:1 for borders and large text).
//
// The theme is switched through the app's own OS-dark-mode hook: a Playwright context created
// with `colorScheme: 'dark'` makes `prefers-color-scheme: dark` match, InterSpec's JavaScript
// reports that to the C++ side (InterSpec::osThemeChange), and the real predefined Dark theme is
// applied through the real code path.  This needs the "AutoDarkFromOs" preference on (its
// default).  The script checks the `--interspec-color-scheme` token before measuring, so a saved
// user theme that is not the default cannot silently make both runs measure the same palette.
//
// Backgrounds are composited before measuring: `getComputedStyle().backgroundColor` reports the
// declared value, so a translucent wash over a dark surface must be blended down to what the user
// actually sees.
//
// Usage (InterSpec must already be running, e.g. on http://127.0.0.1:8080/):
//   NODE_PATH=../PlaywrightPhoneEmulation/node_modules node theme-contrast.js [--url URL]
//         [--scheme light|dark|both] [--headed] [--screenshots DIR]
//
// One-time setup (shares the phone harness' Playwright install):
//   cd ../PlaywrightPhoneEmulation && npm install && npx playwright install chromium

const { chromium } = require( 'playwright' );

function parseArgs( argv )
{
  const opts = { url: 'http://127.0.0.1:8080/', scheme: 'both', headed: false, screenshots: '' };
  for( let i = 2; i < argv.length; ++i )
  {
    const a = argv[i];
    if( a === '--url' ) opts.url = argv[++i];
    else if( a === '--scheme' ) opts.scheme = argv[++i];
    else if( a === '--headed' ) opts.headed = true;
    else if( a === '--screenshots' ) opts.screenshots = argv[++i];
    else if( a === '--help' || a === '-h' )
    {
      console.log( 'See the header comment of this file for usage.' );
      process.exit( 0 );
    }else
    {
      console.error( 'Unknown arg: ' + a );
      process.exit( 1 );
    }
  }
  return opts;
}

// Runs inside the page: builds a probe of representative elements (so every semantic class and
// widget surface is measured even on an empty session), then measures each probe and a few real
// app elements.  Returns [{name, kind, fg, bg, ratio, need, size, weight}].
function measureInPage()
{
  const parse = ( str ) => {
    const m = str.match( /rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*(?:,\s*([\d.]+))?\s*\)/ );
    if( !m ) return null;
    return { r: +m[1], g: +m[2], b: +m[3], a: (m[4] === undefined) ? 1 : +m[4] };
  };
  const over = ( top, under ) => {
    const a = top.a + under.a * (1 - top.a);
    if( a <= 0 ) return { r: 0, g: 0, b: 0, a: 0 };
    const c = ( t, u ) => (t * top.a + u * under.a * (1 - top.a)) / a;
    return { r: c(top.r, under.r), g: c(top.g, under.g), b: c(top.b, under.b), a: a };
  };
  const lum = ( c ) => {
    const f = ( v ) => { v /= 255; return (v <= 0.03928) ? v / 12.92 : Math.pow( (v + 0.055) / 1.055, 2.4 ); };
    return 0.2126 * f(c.r) + 0.7152 * f(c.g) + 0.0722 * f(c.b);
  };
  const contrast = ( a, b ) => {
    const la = lum(a), lb = lum(b);
    return (Math.max(la, lb) + 0.05) / (Math.min(la, lb) + 0.05);
  };
  // Composite the effective background: walk up until the accumulated alpha is opaque.
  const effectiveBackground = ( el ) => {
    let acc = { r: 0, g: 0, b: 0, a: 0 };
    let layers = [];
    for( let e = el; e; e = e.parentElement )
    {
      const bg = parse( getComputedStyle( e ).backgroundColor );
      if( bg && bg.a > 0 ) layers.push( bg );
      if( layers.length && layers[layers.length - 1].a >= 0.999 ) break;
    }
    // layers are top-most first; composite from the bottom up
    for( let i = layers.length - 1; i >= 0; --i ) acc = over( layers[i], acc );
    if( acc.a < 0.999 )
    {
      // Nothing opaque behind it: assume the theme background (the html/body background)
      const themeBg = parse( getComputedStyle( document.body ).backgroundColor ) || { r: 255, g: 255, b: 255, a: 1 };
      acc = over( acc, themeBg.a > 0 ? themeBg : { r: 255, g: 255, b: 255, a: 1 } );
    }
    return acc;
  };

  const probe = document.createElement( 'div' );
  probe.id = 'ThemeContrastProbe';
  probe.style.cssText = 'position:fixed; left:8px; top:8px; z-index:999999; padding:8px;';
  probe.innerHTML = `
    <span data-probe="body text">Body text</span>
    <span class="Wt-label" data-probe=".Wt-label">Label</span>
    <a href="#" data-probe="link">Link</a>
    <button data-probe="button">Button</button>
    <button disabled data-probe="button[disabled]" data-need="none">Disabled</button>
    <input type="text" value="Input text" data-probe="input[type=text]">
    <span class="ErrorTxt" data-probe=".ErrorTxt">Error</span>
    <span class="WarnTxt" data-probe=".WarnTxt">Warning</span>
    <span class="OkTxt" data-probe=".OkTxt">OK</span>
    <span class="FainterTxt" data-probe=".FainterTxt">Fainter</span>
    <span class="SecondaryTxt" data-probe=".SecondaryTxt">Secondary</span>
    <span class="InfoTxt" data-probe=".InfoTxt">Info</span>
    <div class="Wt-tooltip" style="position:static" data-probe=".Wt-tooltip">Tooltip</div>
    <div class="titlebar" data-probe=".titlebar">Title bar</div>
    <ul class="Wt-popupmenu" style="position:static">
      <li><a><span><span data-probe=".Wt-popupmenu item">Menu item</span></span></a></li>
      <li class="Wt-disabled"><a><span><span data-probe=".Wt-popupmenu .Wt-disabled" data-need="none">Disabled item</span></span></a></li>
    </ul>
    <div class="Wt-itemview" style="position:static">
      <div class="Wt-headerdiv" data-probe=".Wt-itemview header">Header</div>
      <div class="Wt-selected" data-probe=".Wt-itemview .Wt-selected">Selected row</div>
    </div>
    <div class="DrfContentChips">
      <span class="DrfChip" data-probe=".DrfChip">Chip</span>
      <span class="DrfChip DrfChipAbsent" data-probe=".DrfChipAbsent">Absent chip</span>
    </div>`;
  // Inside the app's root container, so the probe composites onto the themed background
  ( document.querySelector( '.Wt-domRoot' ) || document.body ).appendChild( probe );

  // Real app elements worth measuring when present
  const real = [
    ['.Wt-tabs li a', 'tab strip'],
    ['.Wt-tabs .itemselected a', 'selected tab'],
    ['.MenuBar button', 'menu bar button'],
    ['.HeavyNavMenu li.item a', 'nav menu item'],
    ['.DetectorDisplay', 'detector display'],
  ];

  const results = [];
  const measure = ( el, name, opts ) => {
    const cs = getComputedStyle( el );
    const fg = parse( cs.color );
    const bg = effectiveBackground( el );
    if( !fg ) return;
    const size = parseFloat( cs.fontSize );
    const weight = parseInt( cs.fontWeight ) || (cs.fontWeight === 'bold' ? 700 : 400);
    const large = (size >= 24) || (size >= 18.66 && weight >= 700);
    const need = (opts.need === 'none') ? 0 : (large ? 3.0 : 4.5);
    results.push( { name: name, kind: 'text', fg: cs.color, bg: `rgb(${bg.r|0},${bg.g|0},${bg.b|0})`,
                    ratio: contrast( fg, bg ), need: need, size: size, weight: weight } );
    if( opts.border )
    {
      const bc = parse( cs.borderTopColor );
      if( bc && bc.a > 0 )
      {
        // the border sits on the parent's background
        const parentBg = effectiveBackground( el.parentElement );
        const eff = over( bc, parentBg );
        results.push( { name: name + ' border', kind: 'border', fg: cs.borderTopColor,
                        bg: `rgb(${parentBg.r|0},${parentBg.g|0},${parentBg.b|0})`,
                        ratio: contrast( eff, parentBg ), need: 3.0, size: 0, weight: 0 } );
      }
    }
  };

  probe.querySelectorAll( '[data-probe]' ).forEach( el => {
    measure( el, el.getAttribute('data-probe'), { need: el.getAttribute('data-need'), border: el.hasAttribute('data-border') } );
  } );
  for( const [sel, name] of real )
  {
    const el = document.querySelector( sel );
    if( el ) measure( el, name + ' (' + sel + ')', {} );
  }

  probe.remove();
  return results;
}

async function waitForTheme( page, scheme )
{
  const deadline = Date.now() + 20000;
  let seen = '';
  while( Date.now() < deadline )
  {
    seen = await page.evaluate( () => getComputedStyle( document.documentElement ).getPropertyValue( '--interspec-color-scheme' ).trim() );
    if( seen === scheme ) return;
    await page.waitForTimeout( 500 );
  }
  throw new Error( `Expected the ${scheme} theme but --interspec-color-scheme is '${seen}'.  If a saved` +
                   ` theme other than "Default" is selected (Help -> Color Themes), or the AutoDarkFromOs` +
                   ` preference is off, the OS-scheme hook cannot switch themes; reset those and retry.` );
}

( async () => {
  const opts = parseArgs( process.argv );
  const schemes = (opts.scheme === 'both') ? ['light', 'dark'] : [opts.scheme];
  const browser = await chromium.launch( { headless: !opts.headed } );
  let failures = 0;

  for( const scheme of schemes )
  {
    const context = await browser.newContext( { colorScheme: scheme, viewport: { width: 1400, height: 900 } } );
    const page = await context.newPage();
    await page.goto( opts.url );
    await page.waitForSelector( '.specviewer', { timeout: 60000 } );  // the app's root container
    await waitForTheme( page, scheme );
    await page.waitForTimeout( 1000 );

    const results = await page.evaluate( measureInPage );
    console.log( `\n=== ${scheme} theme ===` );
    console.log( 'result  ratio  need  fg                 bg                 probe' );
    for( const r of results )
    {
      const ok = (r.need === 0) || (r.ratio >= r.need);
      if( !ok ) failures += 1;
      const tag = (r.need === 0) ? 'info' : (ok ? 'PASS' : 'FAIL');
      console.log( `${tag.padEnd(6)}  ${r.ratio.toFixed(2).padStart(5)}  ${String(r.need).padStart(4)}  ` +
                   `${r.fg.padEnd(18)} ${r.bg.padEnd(18)} ${r.name}` );
    }
    if( opts.screenshots )
      await page.screenshot( { path: `${opts.screenshots}/theme-${scheme}.png`, fullPage: false } );
    await context.close();
  }

  await browser.close();
  console.log( failures ? `\n${failures} contrast failure(s)` : '\nAll contrast checks passed' );
  process.exit( failures ? 1 : 0 );
} )().catch( e => { console.error( e.message || e ); process.exit( 2 ); } );
