/* InterSpec: an application to analyze spectral gamma radiation data.

  Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
  Government retains certain rights in this software.

  Parent-side helper for the LLM `run_javascript` tool.

  Hosts a single off-screen, sandboxed iframe (origin = opaque, no
  same-origin) that runs untrusted JavaScript on behalf of the LLM, and
  routes results back to C++ via a Wt JSignal named "jsSandboxResult".

  Usage from C++:
    InterSpec.LlmJsSandbox.init('<wt-widget-id>');
    InterSpec.LlmJsSandbox.run('<requestId>', '<code>', <dataLiteral>, <timeoutMs>);
*/
(function(){
  if( typeof window.InterSpec === "undefined" )
    window.InterSpec = {};
  if( window.InterSpec.LlmJsSandbox )
    return;

  /* The HTML loaded into the sandboxed iframe via srcdoc.  The bootstrap
     installs a console shim, listens for postMessage, runs the supplied
     code inside a try/catch, and posts the result back to the parent. */
  const SANDBOX_HTML =
    "<!doctype html><html><head><meta charset=\"utf-8\"></head><body><script>\n" +
    "(function(){\n" +
    "  var __stdout = [], __stderr = [];\n" +
    "  function __fmt(args){\n" +
    "    try {\n" +
    "      return Array.prototype.map.call(args, function(a){\n" +
    "        if( typeof a === 'string' ) return a;\n" +
    "        try { return JSON.stringify(a); } catch(e){ return String(a); }\n" +
    "      }).join(' ');\n" +
    "    } catch(e){ return ''; }\n" +
    "  }\n" +
    "  window.console = {\n" +
    "    log:   function(){ __stdout.push(__fmt(arguments)); },\n" +
    "    info:  function(){ __stdout.push(__fmt(arguments)); },\n" +
    "    debug: function(){ __stdout.push(__fmt(arguments)); },\n" +
    "    warn:  function(){ __stderr.push(__fmt(arguments)); },\n" +
    "    error: function(){ __stderr.push(__fmt(arguments)); }\n" +
    "  };\n" +
    "  window.addEventListener('message', function(ev){\n" +
    " console.error(\"Error from insideiframe\");" +
    "    var req = ev.data || {};\n" +
    "    var requestId = req.requestId;\n" +
    "    var code = req.code;\n" +
    "    var data = req.data;\n" +
    "    __stdout = []; __stderr = [];\n" +
    "    var startMs = Date.now();\n" +
    "    var result; var resultType = 'undefined';\n" +
    "    var errorMessage = null, errorLine = null, errorName = null;\n" +
    "    try {\n" +
    //"      var fn = new Function('data', code);\n" +
    // If the JS code is a single line eval like `2 + 2;` we want to get the statments value - but if its multi-line, the `return (...);` will throw, which then hopefully its last line returns a value.
    "      var fn;\n" +
    "      try {\n" +
    "        fn = new Function('data', 'return (' + code + ');');\n" +
    "      } catch(exprErr) {\n" +
    "        fn = new Function('data', code);\n" +
    "      }\n" +
    "      result = fn(data);\n" +
    "      resultType = (result === null) ? 'null' : (typeof result);\n" +
    "    } catch(e){\n" +
    "      errorMessage = (e && e.stack) ? String(e.stack) : String(e);\n" +
    "      errorName = (e && e.name) ? String(e.name) : null;\n" +
    "      if( e && (typeof e.lineNumber === 'number') ) errorLine = e.lineNumber;\n" +
    "    }\n" +
    "    var resultJson = null, resultString = null;\n" +
    "    if( !errorMessage ){\n" +
    "      try {\n" +
    "        if( resultType === 'undefined' ){\n" +
    "          resultJson = null;\n" +
    "        } else if( resultType === 'function' || resultType === 'symbol' ){\n" +
    "          resultString = String(result);\n" +
    "        } else {\n" +
    "          resultJson = JSON.stringify(result, function(k, v){\n" +
    "            if( typeof v === 'number' && !isFinite(v) ) return String(v);\n" +
    "            if( typeof v === 'bigint' ) return v.toString();\n" +
    "            return v;\n" +
    "          });\n" +
    "          if( resultJson === undefined ){\n" +
    "            resultJson = null;\n" +
    "            resultString = String(result);\n" +
    "          }\n" +
    "        }\n" +
    "      } catch(e){\n" +
    "        resultJson = null;\n" +
    "        try { resultString = String(result); } catch(e2){ resultString = '[unserializable]'; }\n" +
    "      }\n" +
    "    }\n" +
    "    var elapsed = Date.now() - startMs;\n" +
    "    try {\n" +
    "      parent.postMessage({\n" +
    "        requestId: requestId,\n" +
    "        resultJson: resultJson,\n" +
    "        resultString: resultString,\n" +
    "        resultType: resultType,\n" +
    "        stdout: __stdout,\n" +
    "        stderr: __stderr,\n" +
    "        errorMessage: errorMessage,\n" +
    "        errorName: errorName,\n" +
    "        errorLine: errorLine,\n" +
    "        executionTimeMs: elapsed\n" +
    "      }, '*');\n" +
    "    } catch(e){}\n" +
    "  });\n" +
    "})();\n" +
    "<\/script></body></html>";

  let bridgeId = null;     // Wt id of the LlmJsSandboxBridge container
  let iframe = null;       // Current iframe element (rebuilt after kill)
  const pending = {};      // requestId -> { timeoutHandle, timeoutMs }

  function buildIframe()
  {
    const f = document.createElement( "iframe" );
    f.setAttribute( "sandbox", "allow-scripts" );
    f.style.cssText = "position:absolute;left:-9999px;top:-9999px;width:1px;height:1px;border:0;";
    f.__ready = false;
    f.__queue = [];
    f.addEventListener( "load", function(){
      console.log("Sandbox iframe loaded");
      f.__ready = true;
      while( f.__queue.length ){
        const m = f.__queue.shift();
        try { f.contentWindow.postMessage( m, "*" ); } catch(e){}
      }
    });
    f.setAttribute( "srcdoc", SANDBOX_HTML );
    const container = document.getElementById( bridgeId );
    if( container )
      container.appendChild( f );
    else
      document.body.appendChild( f );
    return f;
  }

  function killIframe()
  {
    if( iframe ){
      try { iframe.parentNode && iframe.parentNode.removeChild(iframe); } catch(e){}
    }
    iframe = null;
  }

  function postResult( payload )
  {
    try {
      Wt.emit( bridgeId, "jsSandboxResult", JSON.stringify(payload) );
    } catch(e){
      try { console.error( "LlmJsSandbox: failed to emit result:", e ); } catch(e2){}
    }
  }

  /* The parent's message listener routes responses from the iframe back to
     the C++ side via the JSignal.  Stale (post-timeout) responses are
     dropped silently. */
  window.addEventListener( "message", function(ev){
    const d = ev.data || {};
    const requestId = d.requestId;
    if( !requestId || !pending[requestId] ) return;
    if( iframe && ev.source !== iframe.contentWindow ) return;
    const entry = pending[requestId];
    clearTimeout( entry.timeoutHandle );
    delete pending[requestId];
    postResult({
      requestId: requestId,
      timedOut: false,
      resultJson: d.resultJson,
      resultString: d.resultString,
      resultType: d.resultType,
      stdout: d.stdout || [],
      stderr: d.stderr || [],
      errorMessage: d.errorMessage,
      errorName: d.errorName,
      errorLine: d.errorLine,
      executionTimeMs: d.executionTimeMs
    });
  });

  window.InterSpec.LlmJsSandbox = {
    init: function( id ){
      bridgeId = id;
    },

    run: function( requestId, code, data, timeoutMs ){
      console.log( "In run" );

      if( !bridgeId ){
        try { console.error( "LlmJsSandbox.run() called before init()" ); } catch(e){}
        return;
      }

      const safeTimeout = Math.min( Math.max( parseInt(timeoutMs, 10) || 5000, 1 ), 30000 );

      if( !iframe )
        iframe = buildIframe();

      console.log( "Built iframe run" );

      const handle = setTimeout( function(){
        console.log( "Timeout is firing" );
        if( !pending[requestId] ) return;
        console.log( "Timeout is doing somethign" );
        /* Snapshot all in-flight requests before tearing down the iframe -
           every one of them is now lost, so each gets a timedOut response. */
        const inflight = Object.keys( pending );
        killIframe();
        inflight.forEach( function( id ){
          const e = pending[id];
          if( !e ) return;
          clearTimeout( e.timeoutHandle );
          delete pending[id];
          const isMe = (id === requestId);
          postResult({
            requestId: id,
            timedOut: true,
            resultJson: null,
            resultString: null,
            resultType: "undefined",
            stdout: [],
            stderr: [],
            errorMessage: isMe
              ? ("Execution exceeded " + safeTimeout + "ms timeout; sandbox terminated")
              : "Aborted: a sibling request exceeded its timeout and the sandbox was terminated",
            errorName: "TimeoutError",
            errorLine: null,
            executionTimeMs: isMe ? safeTimeout : 0
          });
        });
      }, safeTimeout + 50 );

      pending[requestId] = { timeoutHandle: handle, timeoutMs: safeTimeout };

      const message = { requestId: requestId, code: code, data: data, timeoutMs: safeTimeout };
      if( iframe.__ready ){
        try {
          console.log( "Posting message to iframe.", message );
          iframe.contentWindow.postMessage( message, "*" );
        } catch(e){
          console.log( "Caught exception posting message to iframe." );
          clearTimeout( handle );
          delete pending[requestId];
          postResult({
            requestId: requestId,
            timedOut: false,
            resultJson: null,
            resultString: null,
            resultType: "undefined",
            stdout: [],
            stderr: [],
            errorMessage: "Failed to dispatch to sandbox: " + (e && e.message ? e.message : String(e)),
            errorName: "DispatchError",
            errorLine: null,
            executionTimeMs: 0
          });
        }
      } else {
        iframe.__queue.push( message );
      }
    }
  };
})();
