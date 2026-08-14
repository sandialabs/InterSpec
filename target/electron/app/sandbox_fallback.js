/* Decides whether a renderer death is attributable to a broken Chromium sandbox, and tracks how
   many consecutive sandboxed launches have failed.

   This lives apart from main.js, and deliberately requires nothing but `fs`, so it can be tested
   with plain node - see ../linux/test_sandbox_fallback.js.  The logic is worth testing because
   acting on it is destructive: main.js responds by calling app.relaunch() + app.exit(), which
   takes the whole application down without running `before-quit` or the window close handlers,
   discarding any unsaved work.  A false positive therefore costs the user real data, and a
   false negative leaves them looking at a blank white window.

   See the comment block in main.js for why each condition below exists.
 */

'use strict';

const fs = require('fs');

/** How long after a window is created a renderer death can still be blamed on the sandbox.
 *  A broken sandbox kills the renderer immediately; a death minutes later is a real crash.
 */
const SANDBOX_FALLBACK_WINDOW_MS = 60000;

/** Consecutive failed sandboxed launches before we stop attempting one.  More than one, so a
 *  single freak crash cannot silently disable the sandbox for good.
 */
const SANDBOX_FAILURES_BEFORE_GIVING_UP = 2;

/** Renderer-exit reasons a broken sandbox actually produces.
 *
 *  The /dev/shm failure in https://github.com/sandialabs/InterSpec/issues/51 reports 'crashed';
 *  a sandbox that cannot start the process at all reports 'launch-failed'.
 *
 *  Electron can also report 'clean-exit', 'oom', 'killed' and 'integrity-failure'.  Do not add
 *  'oom' or 'killed': InterSpec can legitimately exhaust memory on a large search-mode file, and
 *  treating that as a sandbox problem would both lose the user's work and disable the sandbox for
 *  a reason that has nothing to do with it.
 */
const SANDBOX_FAILURE_REASONS = [ 'launch-failed', 'crashed', 'abnormal-exit' ];


/** Is this renderer death attributable to a broken sandbox?
 *
 * @param window   The BrowserWindow it happened to; needs `windowNumber`, `appHasLoadConfirmed`
 *                 and `createdAt` as set in main.js's createWindow().
 * @param details  Electron's `render-process-gone` details; needs `reason`.
 * @param now      Current epoch ms; injectable so tests need not sleep.
 */
function isSandboxFailure( window, details, now ){
  if( typeof now !== 'number' )
    now = Date.now();

  if( !window || !details )
    return false;

  // Only the first window of the launch.  A later window (from "New app window") also starts out
  //  with appHasLoadConfirmed false, so keying on that alone would let its renderer crash kill a
  //  first window holding hours of unsaved work.
  if( window.windowNumber !== 1 )
    return false;

  // If the session already loaded, the sandbox plainly works and this is a genuine crash.
  if( window.appHasLoadConfirmed )
    return false;

  if( (typeof window.createdAt !== 'number')
      || ((now - window.createdAt) > SANDBOX_FALLBACK_WINDOW_MS) )
    return false;

  return SANDBOX_FAILURE_REASONS.indexOf( details.reason ) >= 0;
}//function isSandboxFailure


/** Consecutive failed sandboxed launches recorded so far; 0 if unknown or unreadable.
 *
 * Any problem reading the file counts as 0, so a corrupt or truncated state file fails safe
 * toward *attempting* the sandbox rather than toward disabling it.
 */
function readSandboxFailures( statePath ){
  try{
    if( !fs.existsSync(statePath) )
      return 0;
    const state = JSON.parse( fs.readFileSync(statePath, 'utf8') );
    return Number.isInteger(state.failures) ? state.failures : 0;
  }catch(e){
    return 0;
  }
}//function readSandboxFailures


function writeSandboxFailures( statePath, failures, reason ){
  try{
    fs.writeFileSync( statePath, JSON.stringify({
      failures: failures,
      lastReason: reason,
      updated: new Date().toISOString(),
      note: "InterSpec's first window could not start with the Chromium sandbox enabled.  After "
            + SANDBOX_FAILURES_BEFORE_GIVING_UP + " consecutive failures it stops trying and"
            + " starts with --no-sandbox.  A successful sandboxed start deletes this file."
            + "  Delete it yourself to try the sandbox again now."
    }, null, 1) );
    return true;
  }catch(e){
    return false;
  }
}//function writeSandboxFailures


/** Called once a session has loaded while sandboxed - the sandbox works here after all, so the
 *  state self-heals instead of one bad day disabling it permanently.
 */
function clearSandboxFailures( statePath ){
  try{
    if( fs.existsSync(statePath) ){
      fs.unlinkSync( statePath );
      return true;
    }
  }catch(e){
    // Not important enough to bother the user about.
  }
  return false;
}//function clearSandboxFailures


function sandboxKnownBroken( statePath ){
  return readSandboxFailures(statePath) >= SANDBOX_FAILURES_BEFORE_GIVING_UP;
}


module.exports = {
  SANDBOX_FALLBACK_WINDOW_MS,
  SANDBOX_FAILURES_BEFORE_GIVING_UP,
  SANDBOX_FAILURE_REASONS,
  isSandboxFailure,
  readSandboxFailures,
  writeSandboxFailures,
  clearSandboxFailures,
  sandboxKnownBroken,
};
