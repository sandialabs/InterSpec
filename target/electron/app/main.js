/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* ToDo list (partial):
   - Handle fatal errors with dialog.showErrorBox(...)
   - Catch 'IntializeError' in stderr during startup, and handle
   - Look at creating a backup preferences file, and if the C++ fails to start
     2 or 3 times, go back to the previous preferences file (should be done for
     all targets maybe).
   - move checkWindowPosition() into its own file.
   - Test the window positon stuff with multiple displays.
 */

const electron = require('electron')

const interspec = require('./InterSpecAddOn.node');

const {dialog} = electron;
const {Menu, MenuItem} = electron;
//const {systemPreferences} = electron;

const http = require('http');
const path = require('path')
var fs = require("fs")
const url = require('url')

// Module to control application life.
const app = electron.app
// Module to create native browser window.
const BrowserWindow = electron.BrowserWindow

//Create the 'apptoken' we will designate for the main window
const crypto = require('crypto');


let initial_file_to_open = null;
let interspec_url = null;

global.__basedir = __dirname;


// Keep a global reference of the window objects, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
// We will keep windows stored with the most recently used window as the last 
//  element in the array
let openWindows = [];

// Just to debug, define a windowNumber that will increment with each newly created window
//  (This can be removed)
let windowNumber = 1;


if (process.defaultApp) {
  if (process.argv.length >= 2) {
    app.setAsDefaultProtocolClient('interspec', process.execPath, [path.resolve(process.argv[1])]);
    app.setAsDefaultProtocolClient('raddata', process.execPath, [path.resolve(process.argv[1])]);
  }
} else {
    app.setAsDefaultProtocolClient('interspec');
    app.setAsDefaultProtocolClient('raddata');
}

const gotTheLock = app.requestSingleInstanceLock();

if( !gotTheLock ) 
{
  // TODO: need to make sure the 'will-quit' and 'before-quit', and beforeunload and unload handlers are fine being called
  app.quit(); //Try to close all windows; however the quitting could be cancelled.
  app.exit(-1); //Exits immediately.
}else 
{
  app.on('second-instance', (event, commandLine, workingDirectory) => {  
    // Someone tried to run a second instance, we should focus our window.
    console.log( "Second instance: workingDirectory:", workingDirectory, ", commandLine:", commandLine, ", event:", event );

    const infiles = argvToPaths(commandLine,workingDirectory);


    if( infiles.length > 0 )   
    {
      let window = (openWindows.length ? openWindows[openWindows.length-1] : null);

      if( window ) 
      {
        if( window.isMinimized() ) 
          window.restore();
        window.focus();

        interspec.openFile( window.appSessionToken, JSON.stringify(infiles) );
      }else
      {
        // I dont think we will get here; so lets show an error, for development purposes, incase it happens, so I'll see it.
        dialog.showErrorBox('second instance with no previous window', `infiles.length: ${infiles.length}, workingDirectory: ${workingDirectory}, commandLine: ${argvstr}`)
        initial_file_to_open = infiles;
      }//if( window ) / else
    }else
    {
      // If the user double-clicked the executable again, then the only command line arguments will be executable 
      //  path and the "--allow-file-access-from-files" flag
      // Lets open a new window
      createWindow();
    }

    //let argvstr = commandLine.join(', ');
    //dialog.showErrorBox('second instance', `infiles.length: ${infiles.length}, workingDirectory: ${workingDirectory}, commandLine: ${argvstr}`)
  });

  
  app.on('open-url', (event, url) => {
    event.preventDefault();
    console.log( `'open-url': url: ${url}` );

    const lower = url.toLowerCase();
    if( !lower.startsWith("interspec://") && !lower.startsWith("raddata://g0/") )
      return;

    let window = (openWindows.length ? openWindows[openWindows.length-1] : null);

    if( window && window.pageHasLoaded && window.appHasLoadConfirmed )
      interspec.openAppUrl( window.appSessionToken, url );
    else
      appendInitialFileToLoad(url);
  })
}//if (!gotTheLock) /


/** Loops over argv_array arguments to look for actual files or "deeplinks" (urls starting with "interspec://") to open.
 * 
 * @param {*} argv_array Array of strings passes to the executable; the first entry is assumed to be executable path, so will be discarded.
 * @param {*} workingDir Optional argument of path to look for files relative to; if not provided, the CWD is used.
 * @returns Returns an array of strings of actual files (with absolute paths), and app-urls, that you can have InterSpec try to open.
 */
function argvToPaths( argv_array, workingDir ) {
  let infiles = [];
  if( !argv_array || (argv_array.length < 2) )
    return infiles;
    
  for( let path_string of argv_array.slice(1) )
  {
    if( path_string.toLowerCase().startsWith("interspec://") 
        || path_string.toLowerCase().startsWith("raddata://g0/") )
    {
      infiles.push( path_string );
    }else
    {
      let found = false;
      if( workingDir )
      {
        try {
          const abspath = path.join(workingDir, path_string);
          if( fs.lstatSync( abspath ).isFile() ){
            infiles.push( path.resolve( abspath ) );
            found = true;
          }
        }catch(e){

        }
      }

      if( !found )
      {
        try 
        {
          if( fs.lstatSync(path_string).isFile() ){
            infiles.push(path.resolve(path_string));
            found = true;
          }
        }catch( e ) 
        {
        }
      }//if( !found )
    }//if( deeplink ) / else (maybe a file )
  }//for( loop over command line ards )

  return infiles
}//argvToPaths


/** Appends either a single file-path to `initial_file_to_open` or an array of file-paths */
function appendInitialFileToLoad( fileToLoad ){
  if (typeof fileToLoad === 'string') {
    const singlefname = "" + fileToLoad;
    fileToLoad = [];
    fileToLoad.push( singlefname );
  }

  if( initial_file_to_open ){
    if( Array.isArray(initial_file_to_open) )
      fileToLoad = initial_file_to_open.concat(fileToLoad);
    else
      fileToLoad.push(initial_file_to_open);
  }
  
  initial_file_to_open = fileToLoad;
}//appendInitialFileToLoad(...)


function checkWindowPosition(state) {
  //Adapted from https://github.com/Sethorax/electron-window-state-manager/blob/master/src/lib/windowState.js (20171121)
  //  Not the best.  Should find which screen we're currently (mostly) on and shrink
  //  window to fit that display.

  const primaryDisplay = electron.screen.getPrimaryDisplay();

  if( (typeof state.width !== 'undefined') && state.width !== null && state.width > primaryDisplay.bounds.width )
    state.width = primaryDisplay.bounds.width - 10;

  if( (typeof state.height !== 'undefined') && state.height !== null && state.height > primaryDisplay.bounds.height )
    state.height = primaryDisplay.bounds.height - 20;

  //Check if window is in bounds of primary display
  if( (typeof state.x !== 'undefined')
      && (typeof state.y !== 'undefined')
      && (typeof state.width !== 'undefined')
      && (typeof state.height !== 'undefined')
      && state.x >= 0 && state.y >= 0
      && state.x < primaryDisplay.bounds.width
      && state.y < primaryDisplay.bounds.height)
    return;

  //Check if window is in bounds of any non-primary display
  const externalDisplays = electron.screen.getAllDisplays().filter((display) => {
    return display.bounds.x !== 0 || display.bounds.y !== 0;
  });

  for( const display of externalDisplays ){
    if( state.x >= display.bounds.x && state.y >= display.bounds.y
        && state.x < (display.bounds.x + display.bounds.width)
        && state.y < (display.bounds.y + display.bounds.height) )
      return;
  }

  //If we made it here, bounds are not valid
  const primaryDisplayBounds = primaryDisplay.bounds;
  state.width = 0.85 * primaryDisplayBounds.width;
  state.height = 0.85 * primaryDisplayBounds.height;
  state.x = 0.025 * state.width;
  state.y = 0.025 * state.height;
}

//Load single file or an array of files
function load_file(window, filename){
  if( !filename )
    return;

  if( !window.appHasLoadConfirmed || !window.pageHasLoaded ){

    appendInitialFileToLoad( filename );

    console.log( "Will open files " + JSON.stringify(filename) + " once page loads" );
    return;
  }

  // The native side expects a JSON array of file paths; normalize a single
  // string to a one-element array so JSON.stringify produces the right shape.
  const files_array = Array.isArray(filename) ? filename : [filename];
  const files_json = JSON.stringify(files_array);
  console.log( "To IPC Sending: openfile=" + files_json );

  interspec.openFile( window.appSessionToken, files_json );
}


/** True when started with `--test-load`; see the handling further down.  Read here because it
 * also has to suppress the sandbox-fallback relaunch, so CI sees a real failure rather than a
 * silently-retried one.
 */
const is_test_load = process.argv.some( arg => (arg === '--test-load') );

/** True when `--no-sandbox` is already on the real command line.
 *
 * Note we deliberately do *not* call `app.commandLine.appendSwitch('no-sandbox')` anywhere.
 * On Linux that has no effect: Chromium initializes the sandbox and forks the zygote during
 * its own startup, before this script is evaluated, so the switch only counts if it is on the
 * actual argv.  (It used to be appended unconditionally here as a workaround for a
 * "GPU process launch failed: error_code=18" error on one machine.  On Linux that accomplished
 * nothing - users had to pass `--no-sandbox` by hand anyway - while on Windows and macOS,
 * where the switch *is* honored from JS, it disabled a sandbox that works fine.  The narrower
 * `DisableGpuSandbox` setting below covers that original GPU problem.)
 */
const sandbox_disabled_on_cl = process.argv.some( arg => arg.startsWith('--no-sandbox') );


//File opening untested on OSX and Windows
if( process.platform == 'darwin' ) {
  app.on('open-file', function (event,path) {
     event.preventDefault();
   
     console.log( "Got request to open file: " + path );
     
     let window = openWindows.length ? openWindows[openWindows.length-1] : null;
     if( window && window.pageHasLoaded && window.appHasLoadConfirmed )
     {
       load_file(window, path);
     }else
     {
      appendInitialFileToLoad( path );
     }
  })
} else {
  
 // var myArgs = process.argv.slice(1);
 // console.log('myArgs: ', myArgs);
  
  const infiles = argvToPaths(process.argv,null);
  
  if( infiles.length )
  {
    let window = openWindows.length ? openWindows[openWindows.length-1] : null;
    if( window )
      load_file(window,infiles);
    else 
      appendInitialFileToLoad(infiles);
  }
}

/**
 * @returns Returns object with options from; default options overriden first
 * by what is specified in data/desktop_app_settings.json, then 
 * ~user_data/InterSpec_app_settings.json
 */
function get_interspec_options(){
  const settings = {
    proxy: "",
    httpPort: 0,
    restoreSession: true,
    requireToken: true,
    openDevTools: false,
    disableGpuSandbox: false,
    maxUndoSteps: 250
  };

  const getOptionsFromFile = function( filepath ){
    try{
      if( !fs.existsSync(filepath) )
        return;

      const config = JSON.parse(fs.readFileSync(filepath, 'utf8'));

      if( config.hasOwnProperty('ProxySetting') ){
        if( typeof config.ProxySetting !== 'string' )
          throw new Error("ProxySetting must be a string value");
        settings.proxy = config.ProxySetting;
      }

      if( config.hasOwnProperty('HttpPortToServeOn') ){
        if( !Number.isInteger(config.HttpPortToServeOn) || (config.HttpPortToServeOn < 0) )
          throw new Error("HttpPortToServeOn must be a non-negative integer");
        settings.httpPort = config.HttpPortToServeOn;
      }

      if( config.hasOwnProperty('RestorePreviousSession') ){
        if( typeof config.RestorePreviousSession !== "boolean" )
          throw new Error("RestorePreviousSession must be boolean");
        settings.restoreSession = config.RestorePreviousSession;
      }

      if( config.hasOwnProperty('AllowTokenFreeSessions') ){
        if( typeof config.AllowTokenFreeSessions !== "boolean" )
          throw new Error("AllowTokenFreeSessions must be boolean");
        settings.requireToken = !config.AllowTokenFreeSessions;
      }

      if( config.hasOwnProperty('OpenDevTools') ){
        if( typeof config.OpenDevTools !== "boolean" )
          throw new Error("OpenDevTools must be boolean");
        settings.openDevTools = config.OpenDevTools;
      }

      // Escape hatch for the "GPU process launch failed: error_code=18" failure that has been
      //  seen on the odd machine (possibly when running from a network drive).  This is much
      //  narrower than --no-sandbox: the renderer sandbox stays on, only the GPU process's does
      //  not.  Unlike --no-sandbox, this switch *is* honored when set from here.
      if( config.hasOwnProperty('DisableGpuSandbox') ){
        if( typeof config.DisableGpuSandbox !== "boolean" )
          throw new Error("DisableGpuSandbox must be boolean");
        settings.disableGpuSandbox = config.DisableGpuSandbox;
      }

      if( config.hasOwnProperty('MaxUndoSteps') ){
        if( !Number.isInteger(config.MaxUndoSteps) )
          throw new Error("MaxUndoSteps must be an integer");
        settings.maxUndoSteps = config.MaxUndoSteps;
      }
    }catch( error ){
      console.error( 'Error:', error );

      dialog.showErrorBox( "Error", "Error in " + path.basename(filepath) + ":", error.message );
    }
  };//const getOptionsFromFile

  const userdata = app.getPath('userData');
  const app_settings_file = path.join(userdata, "InterSpec_app_settings.json");

  const basedir = path.relative( process.cwd(), path.dirname(require.main.filename) );
  const user_settings_file = path.join(basedir, "data/config/InterSpec_app_settings.json");
  
  getOptionsFromFile( app_settings_file )
  getOptionsFromFile( user_settings_file );

  return settings;
};//function get_interspec_options()

const userdata = app.getPath('userData');
var guiOptionsPath = path.join(userdata, "init.json");
let allowRestorePath = path.join(userdata, "do_restore");
const app_options = get_interspec_options();

if( app_options.disableGpuSandbox )
  app.commandLine.appendSwitch('disable-gpu-sandbox');


/* Chromium's sandbox does not work everywhere, and when it fails it does so by killing the
   renderer, which leaves the user looking at a blank white window with nothing in the GUI to
   explain it.  Known cases:
     - Running from an extracted archive on Linux: `chrome-sandbox` needs to be owned by root
       with mode 4755, which zip/tar cannot convey.  (The .deb/.rpm set this in their postinst,
       and the tarball's launcher script probes for it.)
     - Ubuntu 23.10+ sets kernel.apparmor_restrict_unprivileged_userns=1, which blocks the
       namespace sandbox for binaries with no AppArmor profile.
     - Even with an AppArmor profile granting `userns`, Electron 41 on Ubuntu 24.04+ still
       fails, with the renderer dying on /dev/shm access - see
       https://github.com/sandialabs/InterSpec/issues/51.  That one looks like an upstream
       Chromium regression (Electron 13 was fine), so we cannot package around it.

   Rather than give up the sandbox everywhere for the sake of the systems where it is broken, we
   try it, and if the *first* window's renderer dies before that window ever loaded a session, we
   relaunch once with `--no-sandbox` on the real command line - which is the only place that
   switch has any effect on Linux.

   The three qualifying conditions below all exist to keep this from firing on an ordinary crash,
   because `app.relaunch()` + `app.exit()` takes the whole application down without running
   `before-quit` or the window close handlers, and so discards unsaved work:

     - First window only.  A later window (from "New app window") starts out with
       appHasLoadConfirmed false as well, so keying on that alone would let a second window's
       renderer crash kill a first window that has an hour of unsaved fitting in it.  If the
       first window has not loaded yet, there is by definition nothing open to lose.
     - Shortly after the window was created.  A broken sandbox kills the renderer immediately;
       a death minutes later is a real crash.
     - Only reasons a broken sandbox actually produces.  Notably NOT 'oom' or 'killed':
       InterSpec can legitimately exhaust memory on a large search-mode file, and that must not
       be misdiagnosed as a sandbox problem.

   And it takes two consecutive failed launches before we stop trying, so a single one-off
   renderer crash cannot silently disable the sandbox forever.  A successful sandboxed session
   clears the record, so the state self-heals once the underlying problem goes away.
 */
const sandboxFallback = require('./sandbox_fallback.js');

const sandboxStatePath = path.join(userdata, "sandbox_failures.json");

/** Have we already relaunched this run?  Guards against a relaunch loop. */
let hasRelaunchedWithoutSandbox = false;

/** True when the fallback should be suppressed, so a renderer crash stays a hard failure that CI
 * notices rather than being silently retried.  `--test-load-allow-sandbox-fallback` opts back in,
 * so CI can also test the fallback itself (see the "sandbox fallback" step in
 * .github/workflows/build_app.yml).
 */
const suppress_sandbox_fallback = is_test_load
      && !process.argv.some( arg => (arg === '--test-load-allow-sandbox-fallback') );

/** Relaunch with `--no-sandbox` once, recording the failure.  Returns true if relaunching. */
function relaunchWithoutSandbox( reason ){
  if( suppress_sandbox_fallback || sandbox_disabled_on_cl || hasRelaunchedWithoutSandbox )
    return false;

  hasRelaunchedWithoutSandbox = true;

  const failures = sandboxFallback.readSandboxFailures(sandboxStatePath) + 1;
  sandboxFallback.writeSandboxFailures( sandboxStatePath, failures, reason );

  console.error( "First window's renderer failed before its session loaded (" + reason + ")."
                 + " Relaunching with --no-sandbox.  Consecutive sandboxed-launch failures: "
                 + failures + " of " + sandboxFallback.SANDBOX_FAILURES_BEFORE_GIVING_UP
                 + " (recorded in " + sandboxStatePath + ")" );

  app.relaunch( { args: process.argv.slice(1).concat(['--no-sandbox']) } );
  app.exit(0);

  return true;
}//function relaunchWithoutSandbox

/** If previous runs established the sandbox does not work here, go straight to a relaunch rather
 * than making the user watch another blank window.  Returns true if relaunching.
 */
function relaunchIfSandboxKnownBroken(){
  if( sandbox_disabled_on_cl || suppress_sandbox_fallback
      || !sandboxFallback.sandboxKnownBroken(sandboxStatePath) )
    return false;

  console.log( "Sandboxed launch has failed "
               + sandboxFallback.readSandboxFailures(sandboxStatePath) + " times (see "
               + sandboxStatePath + ") - relaunching with --no-sandbox" );
  hasRelaunchedWithoutSandbox = true;
  app.relaunch( { args: process.argv.slice(1).concat(['--no-sandbox']) } );
  app.exit(0);

  return true;
}//function relaunchIfSandboxKnownBroken


let setAsMostRecentWindow = function(window){
  const index = openWindows.indexOf(window);
  if( index < 0 )
  {
    // We may have removed the window already - e.g., in 
    //dialog.showErrorBox('Error', 'Trying to set window as most recent that isnt in openWindows.');
    return;
  }

  console.log( 'Removing window ' + window.windowNumber + ' from active windows.' );

  openWindows.splice(index, 1);
  openWindows.push(window);
};


function createWindow() {
  
  app.setName( "InterSpec" );
  
  let guiConfig = {};
  try 
  {
    guiConfig = JSON.parse(fs.readFileSync(guiOptionsPath, 'utf8'));
  }catch(e) 
  {
  }
  
  if( (typeof guiConfig.bounds === 'undefined') || guiConfig.bounds == null )
   guiConfig.bounds = {};

  if( openWindows.length )
  {
    let lastWindow = openWindows[openWindows.length-1];
    guiConfig.bounds = lastWindow.getBounds();
    guiConfig.bounds.x += 20;
    guiConfig.bounds.y += 20;
  }

  checkWindowPosition(guiConfig.bounds);

  let windowPrefs = Object.assign({}, guiConfig.bounds);
  if( !windowPrefs.minWidth )
    windowPrefs.minWidth = 200;
  if( !windowPrefs.minHeight )  
    windowPrefs.minHeight = 200;

  //To get nodeIntegration to work, there is som JS hacks in
  //  InterSpecApp::setupDomEnvironment()
  windowPrefs.frame = (process.platform == 'darwin');
  windowPrefs.webPreferences = { nodeIntegration: false, contextIsolation: true, spellcheck: false };

  // Create the new window
  let newWindow = new BrowserWindow( windowPrefs );
  
  // Sequential window number, starting at 1.  Used for logging, and by isSandboxFailure() to
  //  tell the first window of the launch (which may still be starting the sandbox) from later
  //  ones (whose renderer dying is an ordinary crash).
  newWindow.windowNumber = windowNumber++;

  // When this window was created, so isSandboxFailure() can distinguish "died while starting"
  //  from "crashed later on".
  newWindow.createdAt = Date.now();

  // Set an indicator if the page has loaded, as messaged to us through Electron signals
  newWindow.pageHasLoaded = false;

  // Set an indicator if the InterSpec app has loaded, as messaged to us through out C++ code
  newWindow.appHasLoadConfirmed = false;
  
  // Only restore the previous session for the very first window of the app launch.
  //  A --test-load run never restores, so it tests loading rather than whatever state happened
  //  to be left behind.
  let allowRestore = false;
  if( (openWindows.length === 0) && !is_test_load )
  {
    try
    {
      if( fs.lstatSync(allowRestorePath).isFile() )
      {
        allowRestore = true;
        fs.unlinkSync( allowRestorePath );
      }
    }catch(e)
    {
      console.error( 'Exception checking on/deleting allow reload path ("' + allowRestorePath + '"): ' + e );
    }
  }//if( no other windows are open )

  // Add new window to the end (i.e., the most recently used) openWindows array
  openWindows.push( newWindow );

  console.log( 'Adding window ' + newWindow.windowNumber + ' to active windows.' );
  
  
  let hasSetInitialUrl = false;
  
  let setInitialUrl = function(){
    
    hasSetInitialUrl = true;
    
    if( interspec_url ) {
      const session_token_buf = crypto.randomBytes(16);
      newWindow.appSessionToken = session_token_buf.toString('hex');
      interspec.addPrimarySessionToken( newWindow.appSessionToken );
      
    
      let url_to_load = interspec_url  + "?apptoken=" + newWindow.appSessionToken;
      if( initial_file_to_open && ((typeof initial_file_to_open === 'string') || initial_file_to_open.length==1) ) {
        let filepath = (typeof initial_file_to_open === 'string') ? initial_file_to_open : initial_file_to_open[0];
        
        interspec.setInitialFileToLoad( newWindow.appSessionToken, filepath );
        
        initial_file_to_open = null;
      }

      //See https://github.com/electron/electron/blob/master/docs/tutorial/mojave-dark-mode-guide.md
      //  For implementing dark mode (leaving out for the moment since havent had time to test)
      //if( systemPreferences.isDarkMode() )
      //  url_to_load += "&colortheme=dark";
      //Actually should use nativeTheme.shouldUseDarkColors
      //  see https://github.com/electron/electron/blob/master/docs/api/native-theme.md#nativethemeshouldusedarkcolors-readonly
      
      if( !allowRestore || !app_options.restoreSession )
        url_to_load += "&restore=no";
      
      console.log('Will Load ' + url_to_load);
      
      newWindow.loadURL( url_to_load );
    } else {
      let workingdir = path.dirname(require.main.filename);
      newWindow.loadURL( "file://" + path.join(workingdir, "loading.html") );
    }
  };
  
  //If we are behind a proxy, we need to make sure we dont try to resolve the local address
  //  as proxies do all sorts of wierd things that will block us from loading our local 
  //  address (this for example happens to me at work)
  //ToDo: the csv list of local stuff is way overkill and can probably be reduced to one 
  //      value, but I was too lazy to test out what will work (requires transfering over 
  //      remote desktop currently), and also I guess should actually do the loadUrl() call
  //      from the setback, but the current way seems to work right now...
  
  const proxy_options = {
    proxyBypassRules: 'local,<local>,127.0.0.1,http://127.0.0.1,http://<local>,http://localhost'
  };

  if( !app_options.proxy || (app_options.proxy.length === 0) ){
    // Nothing to do here
  }else if( app_options.proxy.toLowerCase() == "direct" ){
    proxy_options.mode = "direct";
  }else if( app_options.proxy.toLowerCase() == "auto_detect" ){
    proxy_options.mode = "auto_detect";
  }else if( app_options.proxy.toLowerCase() == "system" ){
    proxy_options.mode = "system";
  }else{
    proxy_options.proxyRules = app_options.proxy;
  }

  
  let ses = newWindow.webContents.session;
  ses.setProxy( proxy_options ).then( function(){
    if( hasSetInitialUrl )
      return;
    console.log('Bypassing proxy for local');
    setInitialUrl();
  } );

  
  //JIC setProxy(...) doesnt call the callback or something... probably not actually needed.
  setTimeout( function() {
    if( !hasSetInitialUrl ){
      console.log('ses.setProxy Took longer than 500ms.');
      setInitialUrl();
    }
  }, 500 );
  
  
  // Open the developer tools.
  if( app_options.openDevTools )
    newWindow.webContents.openDevTools({mode: "bottom"});

  // A nice way to have the renderes console.log show up on the command line
  //  when running for development.
  //newWindow.webContents.on('console-message', (event, level, message, line, sourceId) => {
  //  //https://www.electronjs.org/docs/api/web-contents#event-console-message  
  //  //console.log( sourceId+ " ("+line+"): " + message );
  //  console.log( "From renderer: " + message );
  //});

  // Emitted when the window is closed.
  newWindow.on('closed', function () {
    //Could tell InterSpec to save user place...
    
    console.log( "Writing config: " + JSON.stringify(guiConfig) );
    fs.writeFileSync(guiOptionsPath, JSON.stringify(guiConfig));

    // Dereference the window object
    console.log( 'Removing window ' + newWindow.windowNumber + ' from active windows.' );
    const index = openWindows.indexOf(newWindow);
    if( index >= 0 )
      openWindows.splice(index, 1);
    else
      console.error( "on(closed): Trying to remove window not in openWindows" );
  });

  
  newWindow.on( 'blur', function(){ 
    //Emitted when the window loses focus.
    //Could save the work here.
    interspec.sendMessageToRenderer( newWindow.appSessionToken, "OnBlur"); 
  } );

  newWindow.on( 'focus', function(){ 
    //Emitted when the window gains focus.
    setAsMostRecentWindow(newWindow);
    interspec.sendMessageToRenderer( newWindow.appSessionToken, "OnFocus"); 
  } );

  newWindow.on( 'unmaximize', function(){ 
    //Emitted when the window exits from a maximized state.
    setAsMostRecentWindow(newWindow);
    interspec.sendMessageToRenderer( newWindow.appSessionToken, "OnUnMaximize"); 
  } );
  
  newWindow.on( 'maximize', function(){ 
    //Emitted when window is maximized.
    setAsMostRecentWindow(newWindow);
    interspec.sendMessageToRenderer( newWindow.appSessionToken, "OnMaximize"); 
  } );
  
  newWindow.on( 'leave-full-screen', function(){ 
    setAsMostRecentWindow(newWindow);
    interspec.sendMessageToRenderer( newWindow.appSessionToken, "OnLeaveFullScreen"); 
  } );
  
  newWindow.on( 'enter-full-screen', function(){ 
    setAsMostRecentWindow(newWindow);
    interspec.sendMessageToRenderer( newWindow.appSessionToken, "OnEnterFullScreen"); 
  } );

  newWindow.webContents.on('will-navigate', function(event, url){
    //Emitted when a user or the page wants to start navigation. It can happen
    //  when the window.location object is changed or a user clicks a link in
    //  the page.
    //This event will not emit when the navigation is started programmatically
    //  with APIs like webContents.loadURL and webContents.back.
    //It is also not emitted for in-page navigations, such as clicking anchor
    //  links or updating the window.location.hash.
    
    console.log('webContents: will-navigate');
    if( !url.startsWith(interspec_url) ) {
      console.log( "Will prevent Opening URL=" + url + ", newWindow.webContents.getURL()=" + newWindow.webContents.getURL() );
      event.preventDefault();
      electron.shell.openExternal(url)
    } else {
      //We seem to only get here if the JS application dies and the message saying
      // "The application has stopped running, would you like to restart?" and
      // the user clicks okay.
      // OR if we create a link in the app somewhere (like for a CSV file) and dont
      //   call WAnchor::setTarget(AnchorTarget::TargetNewWindow) for the link (and in 
      //   this case we can probably see if the URL looks somethign like:
      //   http://127.0.0.1:57851/?wtd=oiaGAdsaiwqAs&request=resource&resource=asSaEwq&rand=65
      //   but I didnt bother about this yet)
      //(as of 20191012 only tested by calling wApp->quit() from c++).
    }
  });


  newWindow.webContents.setWindowOpenHandler(({ url }) => {
    //console.log( 'url=' + url );
    //console.log( 'frameName=' + frameName );
    //console.log( 'additionalFeatures=' + additionalFeatures );

    if( url.startsWith(interspec_url) ) {
      //Lets prevent a weird popup window that the user has to close...
      //  I think this is because InterSpec targets a new window for downloads.
      newWindow.webContents.downloadURL(url);
    } else {
      //Keep the page from navigating away from the InterSpec, to say somewhere 
      //  like http://www.boost.org from on the "about" page.
      electron.shell.openExternal(url);
    }

    return { action: 'deny' };
  } );

  
  //For windows, if we want to customize the title of the download dialog
  //  could use the following (doesnt seem to be needed on macOS, but leaving in for consistency)
  newWindow.webContents.session.on('will-download', (event, item, webContents) => {
    //console.log( 'will-download: ' + item.getFilename() );

    // For some reason, I dont understand, clicking on a link to download a file in 
    //  one window, will then cause all the other open windows to get the 'will-download' 
    //  signal; so we'll detect this, and only do an action for the window that emits 
    //  the original  download.
    if( newWindow.webContents !== webContents )
      return;

    let fname = "Untitled";
    try 
    { 
      fname = item.getFilename(); 
    }catch( e ) 
    { 
    }

    item.once('done', (event, state) => {
      if (state === 'completed') {
        console.log('Download successfully')
      } else {
        console.log(`Download failed: ${state}`)
      }
    })
    
    let dialog_options = { title: "Save File (" + fname + ")" };

    try {
      fs.accessSync(guiConfig.defaultSavePath,fs.constants.F_OK|fs.constants.W_OK);
      dialog_options.defaultPath = guiConfig.defaultSavePath;
    } catch(e) {
      dialog_options.defaultPath = app.getPath('downloads');
    }

    dialog_options.defaultPath = path.join(dialog_options.defaultPath, fname);
    
    let filename = dialog.showSaveDialogSync( newWindow, dialog_options );
    
    if( typeof filename == "undefined" )
    {
      item.cancel();
      return;
    }
    
    console.log( "Saving to: " + filename );
    
    try 
    {
      guiConfig.defaultSavePath = path.dirname(filename);
      console.log( "Setting save path to: " + guiConfig.defaultSavePath );
      
      item.setSavePath(filename);
    }catch( e ) 
    {
      console.log( "Error saving file: " + e );
      item.cancel();
    }
  });

  
  newWindow.webContents.on('did-start-loading', (event) => {
    //Corresponds to the points in time when the spinner of the tab started spinning.
    console.log( "did-start-loading" );
    
  });


  newWindow.webContents.on('did-finish-load', (event) => {
    //Emitted when the navigation is done, i.e. the spinner of the tab has
    //  stopped spinning, and the onload event was dispatched.
    
    console.log( "In did-finish-load!" );
  
    let currentURL = newWindow.webContents.getURL();
    
    console.log("currentURL="+currentURL);
    if( currentURL.includes("loading.html") )
      return;
  
    if( currentURL.includes("apptoken="+newWindow.appSessionToken) ){
    }
  
    newWindow.pageHasLoaded = true;

    if( newWindow.appHasLoadConfirmed && initial_file_to_open )
    {
      load_file(newWindow,initial_file_to_open);
      initial_file_to_open = null;
    }
    
    //I think this is were we set any files to be opened
    
    // Lets make sure the titlebar is in the right state for the window size.
    //  Sending now because our C++ should have the newWindow.appSessionToken at this point.
    interspec.sendMessageToRenderer( newWindow.appSessionToken, newWindow.isMaximized() ? "OnMaximize" : "OnUnMaximize" );
  })

  newWindow.webContents.on('did-fail-load', (event, errorCode, errorDescription, validatedURL, isMainFrame) => {
    //This event is like did-finish-load but emitted when the load failed or was cancelled, e.g. window.stop() is invoked.
    console.log( "Unhandled did-fail-load event!" );
  })
 
  newWindow.webContents.on('render-process-gone', function(event, details){
    console.log('renderer process gone: reason=' + details.reason + ', exitCode=' + details.exitCode);

    // Only retry without the sandbox when this really looks like a sandbox failure - see
    //  isSandboxFailure() and the comment block above it.  Anything else is a genuine crash and
    //  must not take the application down with it.
    if( sandboxFallback.isSandboxFailure( newWindow, details ) )
      relaunchWithoutSandbox( 'reason=' + details.reason + ', exitCode=' + details.exitCode );
  });
  
  
  //newWindow.webContents.on('did-navigate-in-page',function(){
  //  console.log( 'did-navigate-in-page' );
  //});
  //newWindow.webContents.on('destroyed',function(){
  //  console.log( 'destroyed!' );
  //});
  //newWindow.webContents.on('unresponsive',function(){
  //  console.log( 'unresponsive!' );
  //});

  newWindow.on('session-end', function () {
    //Windows only here
    //Could save the work here.
  });
  
  
  newWindow.on('show', function () {
    //Emitted when the window is shown.
    setAsMostRecentWindow(newWindow);
  });
  
  newWindow.on('hide', function () {
    //Emitted when the window is hidden.
  });
  
  newWindow.on('ready-to-show', function () {
    //Emitted when the web page has been rendered (while not being shown) and window can be displayed without a visual flash.
  });
  
  
  newWindow.on('minimize', function () {
    //Emitted when the window is minimized.
  });
  
  newWindow.on('restore', function () {
    //Emitted when the window is restored from a minimized state.
    setAsMostRecentWindow(newWindow);
  });
  
  newWindow.on('resize', function () {
    //Emitted when the window is being resized.
    guiConfig.bounds = newWindow.getBounds();
    setAsMostRecentWindow(newWindow);
  });
  
  newWindow.on('move', function () {
    //Emitted when the window is being moved to a new position.
    guiConfig.bounds = newWindow.getBounds();
    setAsMostRecentWindow(newWindow);
  });
  
  newWindow.on('new-window-for-tab', function () { //macOS
    //Emitted when the native new tab button is clicked.
  });
  
  return newWindow;
}//createWindow


function messageToNodeJs( token, msg_name, msg_data ){
  console.log( 'In js messageToNodeJs, with msg_name="' + msg_name + '" for token="' + token + '"');
  
  let window = null;
  for( let i = 0; i < openWindows.length; ++i ){
    if( openWindows[i].appSessionToken === token ){
      window = openWindows[i];
      break;
    }
  }

  if( !window ){
    console.error( 'messageToNodeJs token: ' + token + "' did not coorespond to any open windows.");
    dialog.showErrorBox('Error', 'Recieved message from a window that couldnt be found - app state is suspect.\nThere are currently ' + openWindows.length + ' open windows.');
    
    return;
  }

  if( msg_name == 'NewCleanSession' ){
    //ToDo: Check that token is valid token; also identify window this token belongs to
    console.log( 'Got NewCleanSession from session with token: ' + token );
    const session_token_buf = crypto.randomBytes(16);
    window.appSessionToken = session_token_buf.toString('hex');
    window.appHasLoadConfirmed = false;
    console.log("New session token: " + window.appSessionToken );
    
    interspec.removeSessionToken( token );
    interspec.addPrimarySessionToken( window.appSessionToken );
    window.loadURL( interspec_url + "?apptoken=" + window.appSessionToken + "&restore=no");
  }else if( msg_name == 'SessionFinishedLoading' ){
    window.appHasLoadConfirmed = true;

    console.log( "Received SessionFinishedLoading for Token='" + token + "'" );

    // A session loaded while sandboxed, so whatever made a previous launch fail is gone - forget
    //  it, rather than letting one bad day disable the sandbox permanently.
    if( !sandbox_disabled_on_cl
        && sandboxFallback.clearSandboxFailures(sandboxStatePath) ){
      console.log( "Sandboxed session loaded - cleared " + sandboxStatePath );
    }

    if( is_test_load ){
      // This is the whole point of --test-load: the session came up, so we are done.
      console.log( "--test-load: session loaded successfully" );
      finishTestLoad( 0 );
      return;
    }

    try{
      fs.writeFileSync(allowRestorePath, ""+Date.now() );
      console.log( 'Wrote reload file: ' + allowRestorePath );
    }catch(e){
      console.log( "Error writing allow reload file " );
    }

    if( initial_file_to_open )
    {
      load_file(window,initial_file_to_open);
      initial_file_to_open = null;
    }
  }else if( msg_name == 'debug-msg' ){
    //Setup a debug message, mostly for timing things.
    var timestamp = '[' + Date.now() + '] ';
    console.log( timestamp + msg_data );
  }else if( msg_name == 'OpenInExternalBrowser' ){

    if( app_options.requireToken ){
      const ext_session_token_buf = crypto.randomBytes(16);
      const ext_session_token = ext_session_token_buf.toString('hex');
      interspec.addExternalSessionToken( ext_session_token );
    
      const ext_url = interspec_url + "?apptoken=" + ext_session_token + "&restore=no"
    
      console.log( `Will try to open ${ext_url} in external browser` );
    
      electron.shell.openExternal(ext_url);
    }else{
      electron.shell.openExternal( interspec_url + "?restore=no" );
    }

    console.log( `Should have opened external URL...` );
  }else if( msg_name == 'NewAppWindow' ){
    createWindow();
  }else if( msg_name == 'MinimizeWindow' ){
    window.minimize();
  }else if( msg_name == 'MaximizeWindow' ){
    window.maximize();
    interspec.sendMessageToRenderer(token, "OnMaximize");
  }else if( msg_name == 'CloseWindow' ){
    window.close();
  }else if( msg_name == "ToggleMaximizeWindow" ){
    if( window.isMaximized() ) {
      window.unmaximize();
      interspec.sendMessageToRenderer(token, "OnUnMaximize");
    } else {
      window.maximize();
      interspec.sendMessageToRenderer( token, "OnMaximize");
    }
  }else if( msg_name == 'ToggleDevTools' ){
    window.toggleDevTools();
  }else{
    console.log( "messageToNodeJs: unrecognized msg_name:", msg_name, ", token:", token, ", msg_data:", msg_data );
  }
  
  //"ServerKilled"

}//messageToNodeJs


function browseForDirectory( token, title, msg ){
  const { dialog } = require('electron');
  
  let window = (openWindows && openWindows.length) ? openWindows[openWindows.length - 1] : null;

  let dirs = dialog.showOpenDialogSync( window, {
    title: title,
    properties: ['openDirectory'],
    //defaultPath: '/my/previous/path',
    message: msg
  });
  
  if( typeof dirs==='undefined' )
    return null;
  
  return (dirs.length<1) ? '' : dirs[0];
};//function browseForDirectory

/* `--test-load` starts the app, waits for the session to report itself loaded, then exits - so
   CI can tell whether a packaged build actually runs.  The Windows wxWidgets target has had this
   for a while (see InterSpecWxApp.cpp); this is the Electron equivalent, and the exit codes are
   the same values, just positive, since POSIX exit codes are unsigned:
     0  - session loaded
     11 - timed out waiting for the session
     12 - the InterSpec server would not start
 */
const TEST_LOAD_TIMEOUT_MS = 120000;
let test_load_timer = null;
let test_load_finished = false;

function finishTestLoad( code ){
  // `SessionFinishedLoading` can arrive more than once for a session, so make this idempotent.
  if( test_load_finished )
    return;
  test_load_finished = true;

  if( test_load_timer ){
    clearTimeout( test_load_timer );
    test_load_timer = null;
  }
  console.log( "--test-load: exiting with code " + code );

  /* Record the outcome to a file as well as the exit code.  Needed because the sandbox fallback
     relaunches: app.relaunch() only starts the replacement once this process has exited, so the
     shell sees *this* process's exit code and the relaunched instance is detached, its own exit
     code lost.  CI therefore deletes this file, launches, and then polls for it - which is how
     the fallback path gets tested at all.  See the "sandbox fallback" step in
     .github/workflows/build_app.yml.
   */
  try{
    fs.writeFileSync( path.join(userdata, "test_load_result.json"), JSON.stringify({
      code: code,
      sandboxed: !sandbox_disabled_on_cl,
      argv: process.argv.slice(1),
      finished: new Date().toISOString()
    }, null, 1) );
  }catch(e){
    console.error( "--test-load: could not write result file: " + e );
  }

  // Arm the hard-exit backstop *first*.  Both killServer() and app.exit() can block - killServer
  //  waits on the Wt session, and app.exit() has been seen not to terminate when a child process
  //  is wedged (a zombie renderer, hit once while testing under CPU emulation; the process then
  //  sat idle indefinitely).  Registering this after those calls would mean it never gets armed
  //  in exactly the cases it exists for.  A test run that hangs is worse than one that fails.
  setTimeout( function(){
    console.error( "--test-load: did not terminate cleanly; forcing exit" );
    process.exit( code );
  }, 10000 ).unref();

  // Deferred, because we are usually called from inside `messageToNodeJs`, which the C++ addon
  //  invokes - exiting there would tear the process down with native frames still on the stack.
  setImmediate( function(){
    try{
      interspec.killServer();
    }catch(e){
      // Nothing useful to do; we are exiting anyway.
    }
    app.exit( code );
  } );
}//function finishTestLoad

if( is_test_load ){
  console.log( "--test-load: will wait up to " + (TEST_LOAD_TIMEOUT_MS/1000) + "s for the session to load" );
  test_load_timer = setTimeout( function(){
    console.error( "--test-load: timed out waiting for the session to load" );
    finishTestLoad( 11 );
  }, TEST_LOAD_TIMEOUT_MS );
}

// Check if we only want to run
for( let path_string of process.argv ) {
  if( path_string.startsWith("--batch") || path_string.startsWith("/batch") ) {
    console.log( "Will run batch");

    let has_docroot = false, has_userdata = false;
    for( let str of process.argv ) {
      has_docroot = (has_docroot || str.startsWith('--docroot'));
      has_userdata = (has_userdata || str.startsWith('--userdatadir'));
    }

    // These are passed straight through to the C++ option parser, not through a shell, so
    //  quoting the values would make the quotes part of the path.
    if( !has_docroot )
      process.argv.push( "--docroot=" + path.dirname(require.main.filename) );
    if( !has_userdata )
      process.argv.push( "--userdatadir=" + userdata );

    const rcode = interspec.runBatchAnalysis( process.argv );

    console.log( "Batch analysis returned code " + rcode );

    app.quit();
    return;
  }
}


// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', function(){
  // Do this before anything else, so we do not start the server just to throw it away.
  if( relaunchIfSandboxKnownBroken() )
    return;

  // On Windows, the default Electron application menu activates when the user presses Alt,
  // stealing focus from the web content.  The Alt key is used as a modifier for spectrum
  // interactions (e.g., Alt + double-click to fit a peak in the background spectrum),
  // so we remove the default menu on non-macOS platforms.
  if( process.platform !== 'darwin' )
    Menu.setApplicationMenu( null );

  const process_name = require.main.filename;
  //actually process.cwd()==path.dirname(require.main.filename) when running using node from command line

  //It looks like we dont need to change the CWD anymore (I think everywhere in
  // InterSpec no longer assumes a specific CWD), but lets do it anyway.
  process.chdir( path.dirname(require.main.filename) );

  const basedir = path.relative( process.cwd(), path.dirname(require.main.filename) );
  console.log( 'process.cwd()="' + process.cwd() + '"');
  console.log( 'path.dirname(require.main.filename)="' + path.dirname(require.main.filename) + '"');
  console.log( 'basedir="' + basedir + '"');
  
  interspec.setRequireSessionToken( app_options.requireToken );
  interspec.setMaxUndoRedoSteps( app_options.maxUndoSteps );
  interspec.setMessageToNodeJsCallback( messageToNodeJs );
  interspec.setBrowseForDirectoryCallback( browseForDirectory );
  
  const xml_config_path = path.join(basedir, "data/config/wt_config_electron.xml");
  let portnum = app_options.httpPort;
  
  try 
  {
    portnum = interspec.startServingInterSpec( process_name, userdata, basedir, xml_config_path, portnum );
  }catch( e )
  {
    if( is_test_load ){
      console.error( "--test-load: could not start the InterSpec server: " + e.message );
      finishTestLoad( 12 );
      return;
    }

    let window = createWindow();
    var html = [
      "<body>",
        "<h1>Error</h1>",
        e.message,
      "</body>",
    ].join("");

    window.loadURL( "data:text/html;charset=utf-8," + encodeURI(html) );

    return;
  }

  interspec_url = "http://127.0.0.1:" + portnum;
  
  createWindow();
});

// Quit when all windows are closed.
app.on('window-all-closed', function () {
  // On OS X it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  //if (process.platform !== 'darwin') {
    //app.quit()
  //}
  openWindows = [];
  
  interspec.killServer();

  app.quit();
})

app.on('before-quit', function() {
  //Emitted before the application starts closing its windows.
  console.log( "Sending Wt code command to exit" );
  interspec.killServer();
});

app.on('will-quit', function(){
  //Emitted when all windows have been closed and the application will quit.
  console.log( "will-quit" );
});

app.on('quit', function(){
  //Emitted when the application is quitting.
});


app.on('activate', function () {
  // On OS X it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if( openWindows.length )
    openWindows[openWindows.length-1].focus();
  else
    createWindow();
})

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.
