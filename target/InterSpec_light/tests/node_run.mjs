// Replays a JSONL request script (same format as light_cli) through the WASM module in Node.
//  Usage: node tests/node_run.mjs <build_wasm/light_wasm.js> <script.jsonl> [ref_lines.json]
//  Input files named by "loadFile" and "importCALp" requests are copied into the in-memory filesystem
//  first (unless already there, e.g., written by an earlier "exportFile").
//  If a ref_lines.json is given, it is sent as the reference-line library before the script.

import fs from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';

const [ modulePath, scriptPath, refLinesPath ] = process.argv.slice( 2 );
if( !modulePath || !scriptPath )
{
  console.error( 'Usage: node node_run.mjs <light_wasm.js> <script.jsonl> [ref_lines.json]' );
  process.exit( 1 );
}

const require = createRequire( import.meta.url );
const factory = require( path.resolve( modulePath ) );
const Module = await factory();
const call = Module.cwrap( 'light_call', 'string', [ 'string' ] );

let nerrors = 0;
const send = ( req ) => {
  const response = call( JSON.stringify( req ) );
  if( response.startsWith( '{"error":' ) && !('expectError' in req) )  //see check_expectations.py
    nerrors += 1;
  return response;
};

if( refLinesPath )
{
  const lib = JSON.parse( fs.readFileSync( refLinesPath, 'utf8' ) );
  console.log( send( { method: 'setRefLibrary', params: { lib } } ) );
}

Module.FS.mkdirTree( '/in' );
for( const line of fs.readFileSync( scriptPath, 'utf8' ).split( '\n' ) )
{
  if( !line.trim() || line.startsWith( '#' ) )
    continue;

  const req = JSON.parse( line );
  // Files an earlier request exported are already in the in-memory filesystem
  if( ((req.method === 'loadFile') || (req.method === 'importCALp')) && req.params && req.params.path && !Module.FS.analyzePath( req.params.path ).exists )
  {
    const memPath = '/in/' + path.basename( req.params.path );
    Module.FS.writeFile( memPath, fs.readFileSync( req.params.path ) );
    req.params.path = memPath;
  }

  console.log( send( req ) );
}

// Not process.exit(): that can cut off large, piped, stdout writes
process.exitCode = nerrors ? 2 : 0;
