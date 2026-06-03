// Launches a Chromium window sized as a chosen device and points it at a
// running InterSpec server with the ?isphone=1 (or ?istablet=1) URL parameter.
//
// Why Playwright instead of Chrome DevTools "device mode"?
//   - The viewport is the actual rendered viewport (window.innerWidth matches).
//   - No browser chrome / DevTools panel eating into the test area.
//   - User-Agent, device-pixel-ratio, touch capability, and isMobile are all
//     set together via the device descriptor.
//   - Repeatable and scriptable.
//
// Usage:
//   node phone-test.js                              # iPhone 14 portrait, isphone=1
//   node phone-test.js --device "Pixel 7"
//   node phone-test.js --device "iPad Mini" --param istablet
//   node phone-test.js --url http://127.0.0.1:8080/ --device "iPhone SE"
//   node phone-test.js --headless                   # for CI / screenshot capture
//
// First-time setup:
//   npm install playwright            # local install in this dir, or
//   npm install -g playwright         # global install
//   npx playwright install chromium   # one-time browser download
//
// Full list of device presets:
//   node -e "console.log(Object.keys(require('playwright').devices).join('\n'))"

const { chromium, devices } = require( 'playwright' );

function parseArgs( argv )
{
  const opts = {
    url: 'http://127.0.0.1:8080/',
    device: 'iPhone 14',
    param: 'isphone',           // 'isphone' | 'istablet' | '' to send nothing
    headless: false,
  };
  for( let i = 2; i < argv.length; ++i )
  {
    const a = argv[i];
    if( a === '--url' )        opts.url = argv[++i];
    else if( a === '--device' ) opts.device = argv[++i];
    else if( a === '--param' )  opts.param = argv[++i];
    else if( a === '--headless' ) opts.headless = true;
    else if( a === '--help' || a === '-h' )
    {
      console.log( 'See header comments in this file for usage.' );
      process.exit( 0 );
    }
    else
    {
      console.error( 'Unknown arg: ' + a );
      process.exit( 1 );
    }
  }
  return opts;
}

function buildUrl( base, param )
{
  if( !param ) return base;
  const sep = base.includes( '?' ) ? '&' : '?';
  return base + sep + param + '=1';
}

( async () =>
{
  const opts = parseArgs( process.argv );
  const desc = devices[opts.device];
  if( !desc )
  {
    console.error( 'Unknown device: ' + opts.device );
    console.error( 'List options: node -e "console.log(Object.keys(require(\'playwright\').devices).join(\'\\n\'))"' );
    process.exit( 1 );
  }

  const target = buildUrl( opts.url, opts.param );
  console.log( 'Launching ' + opts.device + ' at ' + target );

  const browser = await chromium.launch( { headless: opts.headless } );
  const context = await browser.newContext( { ...desc } );
  const page = await context.newPage();

  await page.goto( target ).catch( e => {
    console.error( 'Navigation failed: ' + e.message );
  } );

  if( opts.headless )
  {
    // Headless mode: give the page a moment to render and exit so screenshots
    // / one-off checks can be done by piping or extending this script.
    await page.waitForTimeout( 3000 );
    await browser.close();
    return;
  }

  console.log( 'Browser is open. Press Ctrl-C in this terminal (or close the' );
  console.log( 'browser window) to exit. Hard-reload with Cmd-Shift-R inside' );
  console.log( 'the browser to pick up CSS changes.' );

  // Exit cleanly when the user closes the browser window.
  browser.on( 'disconnected', () => process.exit( 0 ) );

  // Keep the script alive.
  await new Promise( () => {} );
} )();
