/**
 * InterSpecToast.js: lightweight in-app toast notification system.
 *
 * Replaces qTip2 jGrowl-style notifications with a zero-dependency,
 * cross-platform solution that works in all WebView environments.
 *
 * Usage from C++ (via doJavaScript):
 *   InterSpecToast.show(htmlContent, 'error', 5000);
 *   InterSpecToast.clearByClass('riid');  // remove previous RIID toasts
 */
window.InterSpecToast = (function() {
  "use strict";

  function getContainer() {
    var el = document.getElementById('toast-container');
    if( !el )
    {
      el = document.createElement('div');
      el.id = 'toast-container';
      document.body.appendChild(el);
    }
    return el;
  }

  /** Remove a single toast element with animation. */
  function dismiss( el ) {
    if( !el || el._dismissing )
      return;
    el._dismissing = true;
    el.classList.remove('show');
    // After the CSS transition finishes, remove from DOM
    setTimeout( function(){ if( el.parentNode ) el.parentNode.removeChild(el); }, 400 );
  }

  /**
   * Show a toast notification.
   * @param {string} html       - The message content (may contain HTML).
   * @param {string} type       - CSS type suffix: 'info', 'notice', 'error', 'save', 'riid'.
   * @param {number} durationMs - Auto-dismiss delay in ms (default 5000). 0 = no auto-dismiss.
   * @param {string} title      - Title bar text.
   * @param {string} iconUrl    - URL for the title bar icon (optional).
   */
  function show( html, type, durationMs, title, iconUrl ) {
    var container = getContainer();
    if( !container )
      return;

    if( typeof durationMs !== 'number' || durationMs <= 0 )
      durationMs = 5000;

    var toast = document.createElement('div');
    toast.className = 'ispec-toast ispec-toast-' + (type || 'info');

    // Build title bar
    var titleBar = document.createElement('div');
    titleBar.className = 'ispec-toast-title';

    if( iconUrl )
    {
      var img = document.createElement('img');
      img.src = iconUrl;
      titleBar.appendChild(img);
    }

    var titleText = document.createElement('span');
    titleText.textContent = title || '';
    titleBar.appendChild(titleText);

    // Close button
    var closeBtn = document.createElement('button');
    closeBtn.className = 'ispec-toast-close';
    closeBtn.innerHTML = '&times;';
    closeBtn.onclick = function() { dismiss(toast); };

    toast.appendChild(titleBar);
    toast.appendChild(closeBtn);

    // Body
    var body = document.createElement('div');
    body.className = 'ispec-toast-body';
    body.innerHTML = html;
    toast.appendChild(body);

    container.appendChild(toast);

    // Trigger show animation on next frame so the initial state (opacity:0, max-height:0) is painted
    requestAnimationFrame( function() {
      requestAnimationFrame( function() { toast.classList.add('show'); } );
    });

    // Auto-dismiss with pause-on-hover
    var timer = null;

    function startTimer() {
      if( durationMs > 0 )
        timer = setTimeout( function(){ dismiss(toast); }, durationMs );
    }

    toast.addEventListener( 'mouseenter', function() {
      if( timer ) { clearTimeout(timer); timer = null; }
    });

    toast.addEventListener( 'mouseleave', function() {
      startTimer();
    });

    startTimer();
  }

  /** Remove all visible toasts whose class list contains 'ispec-toast-' + type. */
  function clearByClass( type ) {
    var container = getContainer();
    if( !container )
      return;
    var selector = '.ispec-toast-' + type;
    var existing = container.querySelectorAll(selector);
    for( var i = 0; i < existing.length; i++ )
      dismiss(existing[i]);
  }

  /** Add a CSS class to the container (e.g., 'belowMenu' for electron/wxWidgets builds). */
  function addContainerClass( cls ) {
    getContainer().classList.add(cls);
  }

  return {
    show: show,
    clearByClass: clearByClass,
    addContainerClass: addContainerClass
  };
})();
