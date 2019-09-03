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
   -Finish setting up launch_options.json (see get_launch_options())
   -Get working, and package on Windows
   -Setup a decent way to develop and package the app, rather than abusing Cmake
     --See https://github.com/electron-userland/electron-packager
   -Test out opening files (macOS and Windows; somewhat works on mac)
   -Setup, or figure out, signing app on Windows and macOS
    --On mac you can 'codesign-electron.sh InterSpec.app' to sign
   -Handle errors in c++ by sending through IPC socket once its open
   -Handle fatal errors with dialog.showErrorBox(...)
   -Catch 'IntializeError' in stderr during startup, and handle
   -Catch death in C++ code and display an error
     --Could probably setup general error displaying mechanism of looking if
       there is a window showing loading.html, and if so, display there, and
       if not create a dialog.  Needs more thought.
   -Look at creating a backup preferences file, and if the C++ fails to start
    2 or 3 times, go back to the previous preferences file (should be done for
    all targets maybe).
   -Implement app.makeSingleInstance(...), see https://github.com/electron/electron/blob/master/docs/api/app.md
   -move checkWindowPosition() into its own file.
   -Test the window positon stuff with multiple displays.
   -Setup to allow multiple windows (but dont actually allow yet)
     --Make so a request for a new session is sent to C++, which then sends back
       a URL (which includes the externalid) to connect to
     --Change so externalid is is sent to c++ through IPC, who then sends a
       response when the session can then be loaded
   -Look at using node file menu integration
   -Look into making the c++ code a node addon, see http://www.benfarrell.com/2013/01/03/c-and-node-js-an-unholy-combination-but-oh-so-right/ and a nice guide at https://pspdfkit.com/blog/2018/running-native-code-in-electron-and-the-case-for-webassembly/ makes it look like it will be pretty easy.
 */

const electron = require('electron')

const interspec = require('./build/Release/InterSpecAddOn.node');


//To use IPC, just:
const { ipcMain } = electron;


const {dialog} = electron;
const {Menu, MenuItem} = electron;


const http = require('http');
const path = require('path')
var fs = require("fs")
const url = require('url')

// Module to control application life.
const app = electron.app
// Module to create native browser window.
const BrowserWindow = electron.BrowserWindow

var initial_file_to_open = null;
var page_loaded = false, wtapp_loaded = false;
let interspec_url = null;
let app_is_closing = false;

global.__basedir = __dirname;


exports.test = function(){
  console.log( "interspec_url=" + interspec_url );
};

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
//  (we could setup allowing multiple windows, but not currently properly
//   handled in most places)
let mainWindow

//Create the 'externalid' we will designate for the main window
const crypto = require('crypto');
var session_token = null;


function checkWindowPosition(state) {
  //Adapted from https://github.com/Sethorax/electron-window-state-manager/blob/master/src/lib/windowState.js (20171121)
  //  Not the best.  Should find which screen we're currently (mostly) on and shrink
  //  window to fit that display.
  
  var primaryDisplay = electron.screen.getPrimaryDisplay();
  
  if( (typeof state.width !== 'undefined') && state.width!==null && state.width > primaryDisplay.bounds.width )
    state.width = primaryDisplay.bounds.width - 10;
    
  if( (typeof state.height !== 'undefined') && state.height!==null && state.height > primaryDisplay.bounds.height )
    state.height = primaryDisplay.bounds.height - 20;
    
  //Check if window is in bounds of primary display
  if( (typeof state.x !== 'undefined')
      && (typeof state.y !== 'undefined')
      && (typeof state.x !== 'undefined')
      && (typeof state.width !== 'undefined')
      && (typeof state.height !== 'undefined')
      && state.x >= 0 && state.y >= 0
      && state.x < primaryDisplay.bounds.width
      && state.y < primaryDisplay.bounds.height)
    return
  
  //Find all external displays
  var externalDisplays = electron.screen.getAllDisplays().find((display) => {
    return display.bounds.x !== 0 || display.bounds.y !== 0;
  });
    
  //Check if there are external displays present
  if( externalDisplays ) {
    //Create an array if it is a single display
    if (typeof externalDisplays.length === 'undefined') {
      let singleExternal = externalDisplays;
      externalDisplays = [];
      externalDisplays.push(singleExternal);
    }
      
    //Iterate over each display
    for (let i = 0; i < externalDisplays.length; i++) {
      let display = externalDisplays[i];
      
    //Check if window is in bounds of this external display
    if( state.x >= display.bounds.x && state.y >= display.bounds.y
        && state.x < (display.bounds.x + display.bounds.width)
        && state.y < (display.bounds.y + display.bounds.height)){
          return
        }
    }
  }//if( externalDisplays )
  
  //If we made it here, bounds are not valid
  let primaryDisplayBounds = electron.screen.getPrimaryDisplay().bounds;
  state.width = 0.85*primaryDisplayBounds.width;
  state.height = 0.85*primaryDisplayBounds.height;
  state.x = 0.025*state.width;
  state.y = 0.025*state.height;
};

//Load single file or an array of files
function load_file(filename){
  if( !filename )
    return;
    
  //Make into array
  //if( !Array.isArray(filename) )
  if (typeof filename === 'string') {
    let singlefname = "" + filename;
    filename = [];
    filename.push( singlefname );
  }
  
  if( !wtapp_loaded || !page_loaded || !interspec_url ){
    if( initial_file_to_open ){
       if( Array.isArray(initial_file_to_open) )
         filename = initial_file_to_open.concat(filename);
       else
         filename.push(initial_file_to_open);
    }
    
    initial_file_to_open = filename;
    console.log( "Will open files " + JSON.stringify(filename) + " once page_loaded" );
    return;
  }

  var msg = "openfile=" + JSON.stringify(filename);
  console.log( "" + (typeof filename)+ "To IPC Sending: " +  msg );

  interspec.openFile( session_token, JSON.stringify(filename) );
}


//File opening untested on OSX and Windows
if( process.platform == 'darwin' ) {
  app.on('open-file', function (event,path) {
     event.preventDefault();
   
     console.log( "Got request to open file: " + path );
     
     if( page_loaded )
       load_file(path);
     else
      initial_file_to_open = path;
  })
} else {
  
  var myArgs = process.argv.slice(1);
  console.log('myArgs: ', myArgs);
  
  var infiles = [];
  for( var path_string of myArgs )
  {
      try {
          if (fs.lstatSync(path_string).isFile())
              infiles.push(path.resolve(path_string));
      } catch (e) {
      }
  }
  
  if( infiles.length )
    load_file(infiles);
}

app.on('open-url', function (event,url) {
   event.preventDefault()
   //...
})


const userdata = app.getPath('userData');
var guiOtionsPath = path.join(userdata, "init.json");

function get_launch_options(){
  //InterSpecResourceDir
  //ServerExe
  //UserDataDir
  
  //app.getPath('userData')
  //app.getPath('temp')
}


function doMenuStuff(currentwindow){
  
  currentwindow.setMenu(null);
  
  const template = [{label: 'Edit', submenu: [{role: 'cut'},{role: 'copy'},{role: 'paste'}]},
  { label: 'View', submenu: []},
  { label: 'Tools', submenu:[]},
  { role: 'window', submenu: [{role: 'minimize'},{role: 'close'}]},
  { role: 'help', submenu: [] }
  ];
  
  if (process.platform === 'darwin') {
    template.unshift({label: "InterSpec", submenu: [] });
    // Edit menu.
    //template[1].submenu.push({type: 'separator'},{label: 'Speech',submenu: [{role: 'startspeaking'},{role: 'stopspeaking'}]})
    
    // Window menu.
    template[4].submenu = [
    { label: 'Close', accelerator: 'CmdOrCtrl+W', role: 'close' },
    { label: 'Minimize', accelerator: 'CmdOrCtrl+M', role: 'minimize' },
    { label: 'Zoom', role: 'zoom' },
    { type: 'separator' },
    { label: 'Bring All to Front', role: 'front' }
    ]
  } else {
    template.unshift({ label: 'File', submenu: [] })
  }//if (process.platform === 'darwin') / else

  //Menu.getApplicationMenu()
  const menubar = Menu.buildFromTemplate(template);
  Menu.setApplicationMenu(menubar)
}



function createWindow () {
  
  app.setName( "InterSpec" );
  
  var guiConfig = {};
  try {
    guiConfig = JSON.parse(fs.readFileSync(guiOtionsPath, 'utf8'));
  }catch(e) {
  }
  if( (typeof guiConfig.bounds === 'undefined') || guiConfig.bounds == null )
   guiConfig.bounds = {};
  
  checkWindowPosition(guiConfig.bounds);
  
  var windowPrefs = Object.assign({}, guiConfig.bounds);
  
  //To get nodeIntegration to work, there is som JS hacks in
  //  InterSpecApp::setupDomEnvironment()
  windowPrefs.webPreferences = {nodeIntegration: true, nativeWindowOpen: true };

  mainWindow = new BrowserWindow( windowPrefs );
  

  //If we are behind a proxy, we need to make sure we dont try to resolve the local address
  //  as proxies do all sorts of wierd things that will block us from loading our local 
  //  address (this for example happens to me at work)
  //ToDo: the csv list of local stuff is way overkill and can probably be reduced to one 
  //      value, but I was too lazy to test out what will work (requires transfering over 
  //      remote desktop currently), and also I guess should actually do the loadUrl() call
  //      from the setback, but the current way seems to work right now...
  let ses = mainWindow.webContents.session;
  ses.setProxy( {proxyBypassRules: 'local,<local>,127.0.0.1,http://127.0.0.1,http://<local>,http://localhost'}, 
                function(){ console.log('Bypassing proxy for local'); } );

  doMenuStuff(mainWindow);
  
  
  if( interspec_url ) {
    let msg = interspec_url  + "?externalid=" + session_token;
    if( initial_file_to_open && ((typeof initial_file_to_open === 'string') || initial_file_to_open.length==1) ) {
      let filepath = (typeof initial_file_to_open === 'string') ? initial_file_to_open : initial_file_to_open[0];
      msg += "&specfilename=" + encodeURI(filepath);
      initial_file_to_open = null;
    }
  
    interspec.addSessionToken( session_token );

    console.log('Will Load ' + msg);

    mainWindow.loadURL( msg );
  } else {
    let workingdir = path.dirname(require.main.filename);
    mainWindow.loadURL( "file://" + path.join(workingdir, "loading.html") );
  }
  
  // Open the developer tools.
  //mainWindow.webContents.openDevTools({mode: "bottom"})

  // Emitted when the window is closed.
  mainWindow.on('closed', function () {
    //Could tell InterSpec to save user place...
    
    console.log( "Writing config: " + JSON.stringify(guiConfig) );
    fs.writeFileSync(guiOtionsPath, JSON.stringify(guiConfig));
    
    // Dereference the window object, usually you would store windows
    // in an array if your app supports multi windows, this is the time
    // when you should delete the corresponding element.
    mainWindow = null;
  });


  mainWindow.webContents.on('will-navigate', function(event, url){ 
    console.log('webContents: will-navigate');
    if( !url.startsWith(interspec_url) ) {
      console.log( "Will prevent Opening URL=" + url + ", mainWindow.webContents.getURL()=" + mainWindow.webContents.getURL() );
      event.preventDefault();
      electron.shell.openExternal(url)
    }
  });

  
  mainWindow.webContents.on('new-window', (event, url,frameName,disposition,options,additionalFeatures) => {
    event.preventDefault();
    event.defaultPrevented = true;

    //console.log( 'url=' + url );
    //console.log( 'frameName=' + frameName );
    //console.log( 'additionalFeatures=' + additionalFeatures );

    if( url.startsWith(interspec_url) ) {
      //Lets prevent a weird popup window that the user has to close...
      //  I think this is because InterSpec targets a new window for downloads.
      mainWindow.webContents.downloadURL(url)
    } else {
      //Keep the page from navigating away from the InterSpec, to say somewhere 
      //  like http://www.boost.org from on the "about" page.
      electron.shell.openExternal(url)
    }
  });


  //For windows, if we want to customize the title of the download dialog
  //  could use the following (doesnt seem to be needed on macOS, but leaving in for consistency)
  mainWindow.webContents.session.on('will-download', (event, item, webContents) => {
    //console.log( 'will-download: ' + item.getFilename() );
    var fname = "Untitled";
    try { fname = item.getFilename(); } catch(e) { }

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
    
    let filename = dialog.showSaveDialog( mainWindow, dialog_options );
    
    if (typeof filename == "undefined") {
      item.cancel();
      return;
    }
    console.log( "Saving to: " + filename );
    try {
      guiConfig.defaultSavePath = path.dirname(filename);
      console.log( "Setting save path to: " + guiConfig.defaultSavePath );
      
      item.setSavePath(filename);
    } catch( e ) {
      console.log( "Error saving file: " + e );
      item.cancel();
    }
  });


  mainWindow.webContents.on('did-finish-load', (event) => {
    //Emitted when the navigation is done, i.e. the spinner of the tab has
    //  stopped spinning, and the onload event was dispatched.
    
    console.log( "In did-finish-load!" );
  
    let currentURL = mainWindow.webContents.getURL();
    
    console.log("currentURL="+currentURL);
    if( currentURL.includes("loading.html") )
      return;
  
    if( currentURL.includes("externalid="+session_token) ){
      //We have found main window.  Lets be conservative though and not require
      //  the URL have this options.
    }
  
    page_loaded = true;
    if( wtapp_loaded && initial_file_to_open )
    {
      load_file(initial_file_to_open);
      initial_file_to_open = null;
    }
    //I think this is were we set any files to be opened
  })

  mainWindow.webContents.on('did-fail-load', (event, errorCode, errorDescription, validatedURL, isMainFrame) => {
    //This event is like did-finish-load but emitted when the load failed or was cancelled, e.g. window.stop() is invoked.
    console.log( "Unhandled did-fail-load event!" );
  })
 
  mainWindow.webContents.on('did-get-redirect-request', function(event,oldURL,newURL,isMainFrame,httpResponseCode,requestMethod,referrer,headers){
    console.log( 'did-get-redirect-request from ' + oldURL + ' to ' + newURL );
  });

  mainWindow.webContents.on('crashed', function(event,killed){
    console.log('renderer process ' + (killed ? 'killed' : 'crashed') );
  });
  


  //Could instead customize file download...
  //mainWindow.webContents.session.on('will-download', (event, item, webContents) => {
    // Set the save path, making Electron not to prompt a save dialog.
    //item.setSavePath('/tmp/save.pdf')
    //item.on('updated', (event, state) => { } )
   //item.once('done', (event, state) => { } )
  //})
  
  mainWindow.on('session-end', function () {
    //Windows only here
    //Could save the work here.
  });
  
  mainWindow.on('blur', function () {
    //Emitted when the window loses focus.
    //Could save the work here.
  });
  
  mainWindow.on('focus', function () {
    //Emitted when the window gains focus.
  });
  
  mainWindow.on('show', function () {
    //Emitted when the window is shown.
  });
  
  mainWindow.on('hide', function () {
    //Emitted when the window is hidden.
  });
  
  mainWindow.on('ready-to-show', function () {
    //Emitted when the web page has been rendered (while not being shown) and window can be displayed without a visual flash.
  });
  
  mainWindow.on('maximize', function () {
    //Emitted when window is maximized.
  });
  
  mainWindow.on('unmaximize', function () {
    //Emitted when the window exits from a maximized state.
  });
  
  mainWindow.on('minimize', function () {
    //Emitted when the window is minimized.
  });
  
  mainWindow.on('restore', function () {
    //Emitted when the window is restored from a minimized state.
  });
  
  mainWindow.on('resize', function () {
    //Emitted when the window is being resized.
    guiConfig.bounds = mainWindow.getBounds();
  });
  
  mainWindow.on('move', function () {
    //Emitted when the window is being moved to a new position.
    guiConfig.bounds = mainWindow.getBounds();
  });
  
  mainWindow.on('new-window-for-tab', function () { //macOS
    //Emitted when the native new tab button is clicked.
  });
  
}//createWindow



// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', function(){  

  const session_token_buf = crypto.randomBytes(16);
  session_token = session_token_buf.toString('hex');
  console.log("Session token: " + session_token );

  const process_name = require.main.filename;
  const basedir = path.dirname(require.main.filename);
  const xml_config_path = path.join(basedir, "/data/config/wt_config_electron.xml");

  const portnum = interspec.startServingInterSpec( process_name, userdata, basedir, xml_config_path );
  
  if( portnum <= 0 ){
    //Load
  }

  interspec_url = "http://127.0.0.1:" + portnum;
  
  createWindow();

  ipcMain.on('NewCleanSession', function(evt, oldtoken ){
    //ToDo: Check that oldtoken is valid token; also identify window this token belongs to
    console.log( 'Got NewCleanSession from session with token: ' + oldtoken );
    const session_token_buf = crypto.randomBytes(16);
    session_token = session_token_buf.toString('hex');
    console.log("New session token: " + session_token );
    
    interspec.removeSessionToken( oldtoken );
    interspec.addSessionToken( session_token );
    doMenuStuff(mainWindow);
    mainWindow.loadURL( interspec_url + "?externalid=" + session_token + "&restore=no");
  } );

  ipcMain.on('SessionFinishedLoading', function(evt, token ){
    wtapp_loaded = true;

    console.log( "Received SessionFinishedLoading for Token='" + token + "'" );
      
    //ToDo: should make sure token is what we want
    //ToDo: just aff &specfilename="UriEncodedFilename" to URL argument... (but make sure ALLOW_URL_TO_FILESYSTEM_MAP is enabled)
    
    if( initial_file_to_open )
    {
      load_file(initial_file_to_open);
      initial_file_to_open = null;
    }
  } );
  
  //Setup a debug message, mostly for timing things.
  ipcMain.on('debug-msg', function(evt,msg){
    var timestamp = '[' + Date.now() + '] ';
    console.log( timestamp + msg );
  } );
  
  //"ServerKilled"


////Then from the browsers JS, do:
//const { ipcRenderer } = require('electron');
//ipcRenderer.send('test-event', "Hello dude")

});

// Quit when all windows are closed.
app.on('window-all-closed', function () {
  // On OS X it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  //if (process.platform !== 'darwin') {
    //app.quit()
  //}
  
  interspec.killServer();

  app.quit();
})

app.on('before-quit', function() {
  //Emitted before the application starts closing its windows.
  app_is_closing = true;
  
  console.log( "Sending Wt code command to exit" );
  interspec.killServer();
});

app.on('will-quit', function(){
  //Emitted when all windows have been closed and the application will quit.
  app_is_closing = true;
  console.log( "will-quit" );
});

app.on('quit', function(){
  //Emitted when the application is quitting.
  app_is_closing = true;
});


app.on('activate', function () {
  // On OS X it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (mainWindow === null) {
    createWindow()
  }
})

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.
