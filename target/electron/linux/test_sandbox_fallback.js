#!/usr/bin/env node
/* Tests for app/sandbox_fallback.js - the logic deciding whether a renderer death is a broken
   Chromium sandbox, and how many failures it takes before InterSpec stops trying.

   Worth testing because main.js acts on it with app.relaunch() + app.exit(), which takes the
   whole application down without running `before-quit` or the window close handlers - so a false
   positive costs the user unsaved work, and over-eagerness silently disables a security boundary.
   The original version of this logic keyed only on "this window has not confirmed a session yet",
   which fired for a second window's renderer crash and for an out-of-memory kill.

   No dependencies, no Electron - run directly:
       node target/electron/linux/test_sandbox_fallback.js
 */

'use strict';

const assert = require('assert');
const fs = require('fs');
const os = require('os');
const path = require('path');

const sb = require('../app/sandbox_fallback.js');

let failures = 0;
let ran = 0;

function check( name, fn ){
  ran += 1;
  try{
    fn();
    console.log( "  PASS  " + name );
  }catch(e){
    failures += 1;
    console.log( "  FAIL  " + name + "\n          " + (e && e.message ? e.message : e) );
  }
}

const NOW = 1770000000000;   // fixed, so nothing depends on wall-clock time
const firstWindow = { windowNumber: 1, appHasLoadConfirmed: false, createdAt: NOW };


console.log( "isSandboxFailure - should fire:" );

check( "first window, 'crashed' immediately (the issue #51 /dev/shm signature)", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'crashed'}, NOW), true );
});

check( "first window, 'launch-failed' (sandbox could not start the process)", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'launch-failed'}, NOW), true );
});

check( "first window, 'abnormal-exit'", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'abnormal-exit'}, NOW), true );
});

check( "still inside the attribution window", () => {
  const w = { windowNumber:1, appHasLoadConfirmed:false,
              createdAt: NOW - (sb.SANDBOX_FALLBACK_WINDOW_MS - 1000) };
  assert.strictEqual( sb.isSandboxFailure(w, {reason:'crashed'}, NOW), true );
});


console.log( "\nisSandboxFailure - must NOT fire:" );

check( "second window's renderer dies while window 1 holds unsaved work", () => {
  const w = { windowNumber:2, appHasLoadConfirmed:false, createdAt: NOW };
  assert.strictEqual( sb.isSandboxFailure(w, {reason:'crashed'}, NOW), false );
});

check( "out of memory on a large file is not a sandbox problem", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'oom'}, NOW), false );
});

check( "killed by the OS is not a sandbox problem", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'killed'}, NOW), false );
});

check( "integrity-failure is not a sandbox problem", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'integrity-failure'}, NOW), false );
});

check( "clean exit", () => {
  assert.strictEqual( sb.isSandboxFailure(firstWindow, {reason:'clean-exit'}, NOW), false );
});

check( "crash well after the window was created is a real crash", () => {
  const w = { windowNumber:1, appHasLoadConfirmed:false,
              createdAt: NOW - (sb.SANDBOX_FALLBACK_WINDOW_MS + 1000) };
  assert.strictEqual( sb.isSandboxFailure(w, {reason:'crashed'}, NOW), false );
});

check( "crash after the session loaded (sandbox demonstrably works)", () => {
  const w = { windowNumber:1, appHasLoadConfirmed:true, createdAt: NOW };
  assert.strictEqual( sb.isSandboxFailure(w, {reason:'crashed'}, NOW), false );
});

check( "window with no createdAt stamp", () => {
  const w = { windowNumber:1, appHasLoadConfirmed:false };
  assert.strictEqual( sb.isSandboxFailure(w, {reason:'crashed'}, NOW), false );
});

check( "missing window or details does not throw", () => {
  assert.strictEqual( sb.isSandboxFailure(null, {reason:'crashed'}, NOW), false );
  assert.strictEqual( sb.isSandboxFailure(firstWindow, null, NOW), false );
});

check( "'oom' and 'killed' are not in the reason list", () => {
  assert.ok( sb.SANDBOX_FAILURE_REASONS.indexOf('oom') < 0,
             "'oom' must never be treated as a sandbox failure" );
  assert.ok( sb.SANDBOX_FAILURE_REASONS.indexOf('killed') < 0,
             "'killed' must never be treated as a sandbox failure" );
});


console.log( "\nfailure counting and self-healing:" );

const tmpdir = fs.mkdtempSync( path.join(os.tmpdir(), 'interspec-sandbox-test-') );
const statePath = path.join( tmpdir, 'sandbox_failures.json' );

try{
  check( "fresh install always attempts the sandbox", () => {
    assert.strictEqual( sb.readSandboxFailures(statePath), 0 );
    assert.strictEqual( sb.sandboxKnownBroken(statePath), false );
  });

  check( "one failure does NOT disable the sandbox", () => {
    sb.writeSandboxFailures( statePath, sb.readSandboxFailures(statePath) + 1, 'reason=crashed' );
    assert.strictEqual( sb.readSandboxFailures(statePath), 1 );
    assert.strictEqual( sb.sandboxKnownBroken(statePath), false );
  });

  check( "a successful sandboxed session clears the record", () => {
    assert.strictEqual( sb.clearSandboxFailures(statePath), true );
    assert.strictEqual( sb.readSandboxFailures(statePath), 0 );
    assert.strictEqual( sb.sandboxKnownBroken(statePath), false );
  });

  check( "two consecutive failures do disable it", () => {
    for( let i = 0; i < sb.SANDBOX_FAILURES_BEFORE_GIVING_UP; ++i )
      sb.writeSandboxFailures( statePath, sb.readSandboxFailures(statePath) + 1, 'reason=crashed' );
    assert.strictEqual( sb.readSandboxFailures(statePath),
                        sb.SANDBOX_FAILURES_BEFORE_GIVING_UP );
    assert.strictEqual( sb.sandboxKnownBroken(statePath), true );
  });

  check( "clearing after giving up re-enables attempts", () => {
    sb.clearSandboxFailures( statePath );
    assert.strictEqual( sb.sandboxKnownBroken(statePath), false );
  });

  check( "a corrupt state file fails safe toward attempting the sandbox", () => {
    fs.writeFileSync( statePath, "{ this is not json" );
    assert.strictEqual( sb.readSandboxFailures(statePath), 0 );
    assert.strictEqual( sb.sandboxKnownBroken(statePath), false );
  });

  check( "state file without a 'failures' integer reads as zero", () => {
    fs.writeFileSync( statePath, JSON.stringify({ failures: "two" }) );
    assert.strictEqual( sb.readSandboxFailures(statePath), 0 );
  });

  check( "clearing a file that is not there is not an error", () => {
    fs.unlinkSync( statePath );
    assert.strictEqual( sb.clearSandboxFailures(statePath), false );
  });

  check( "more than one failure is required, whatever the constant is set to", () => {
    assert.ok( sb.SANDBOX_FAILURES_BEFORE_GIVING_UP >= 2,
               "a single crash must not be able to disable the sandbox permanently" );
  });
}finally{
  fs.rmSync( tmpdir, { recursive: true, force: true } );
}

console.log( "\n" + (ran - failures) + "/" + ran + " passed" );
process.exit( failures ? 1 : 0 );
