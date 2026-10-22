/* InterSpec Light - core: the App namespace, WASM bridge, Wt event shim, and state dispatch.
 *
 * All spectrum/peak/calibration state lives in the WASM module; every request returns JSON
 * holding just the parts of the display that changed, which applyState() hands to each
 * component's render function.
 */
'use strict';

const App = {
  module: null,
  lightCall: null,
  /** Components register { name: fn(state) } render callbacks here. */
  renderers: [],
  /** Chart element id -> { eventName: handler(...args) } */
  chartHandlers: {},
  isHighRes: true,
  hasForeground: false,
};
window.App = App;  //for debugging from the browser console


/** Synchronous call into the WASM module; throws on error. */
App.call = function( method, params ) {
  const response = App.lightCall( JSON.stringify( { method: method, params: params || {} } ) );
  const result = JSON.parse( response );
  if( result && (typeof result.error === 'string') )
    throw new Error( result.error );
  return result;
};


/** Calls into WASM after letting the browser paint a "working" indicator, then applies the
 result to the display.  Returns the result, or null on error (which is shown to the user). */
App.run = async function( method, params, opts ) {
  opts = opts || {};
  const busy = document.getElementById( 'busy' );
  if( opts.busy )
  {
    busy.textContent = opts.busyText || 'Working…';
    busy.hidden = false;
    await new Promise( (resolve) => requestAnimationFrame( () => setTimeout( resolve, 0 ) ) );
  }

  try
  {
    const result = App.call( method, params );
    App.applyState( result );
    return result;
  }catch( e )
  {
    console.error( method, e );
    App.toast( e.message || String(e), true );
    return null;
  }finally
  {
    busy.hidden = true;
  }
};


App.applyState = function( state ) {
  if( !state )
    return;

  if( typeof state.isHighRes === 'boolean' )
    App.isHighRes = state.isHighRes;

  for( const render of App.renderers )
  {
    try
    {
      render( state );
    }catch( e )
    {
      console.error( 'Render error', e );
    }
  }

  const msgs = [].concat( state.messages || [], state.message ? [state.message] : [] );
  if( msgs.length )
    App.toast( msgs.join( ' ' ) );
};


let toastTimer = null;
App.toast = function( msg, isError ) {
  const el = document.getElementById( 'toast' );
  el.textContent = msg;
  el.classList.toggle( 'error', !!isError );
  el.hidden = false;
  clearTimeout( toastTimer );
  toastTimer = setTimeout( () => { el.hidden = true; }, isError ? 7000 : 4000 );
};


App.el = function( tag, attrs, children ) {
  const e = document.createElement( tag );
  for( const [key, value] of Object.entries( attrs || {} ) )
  {
    if( key === 'class' ) e.className = value;
    else if( key === 'text' ) e.textContent = value;
    else if( key.startsWith( 'on' ) ) e.addEventListener( key.substring( 2 ), value );
    else if( value !== undefined && value !== null && value !== false ) e.setAttribute( key, value );
  }
  for( const c of [].concat( children || [] ) )
    e.appendChild( (typeof c === 'string') ? document.createTextNode( c ) : c );
  return e;
};


App.fmt = function( value, ndigits ) {
  if( (typeof value !== 'number') || !isFinite( value ) )
    return '';
  return value.toFixed( ndigits );
};


/** SpectrumChartD3/D3TimeChart send events with Wt.emit(elem, {name}, ...args); route them to
 the handlers registered for that chart element. */
window.Wt = {
  emit: function( elem, event ) {
    const args = Array.prototype.slice.call( arguments, 2 );
    const id = (typeof elem === 'string') ? elem : (elem && elem.id);
    const name = (event && event.name) ? event.name : String( event );
    const handler = App.chartHandlers[id] && App.chartHandlers[id][name];
    if( handler )
      handler.apply( null, args );
  }
};
