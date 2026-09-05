/* InterSpec: an open source tool for analyzing gamma spectrums.
 
 Renders LLM Assistant responses (Markdown + LaTeX) into HTML, client-side.
 
 Requires, loaded in this order, before this file is used:
   assets/js/katex-0.18.5/katex.min.js                    -> window.katex
   assets/js/marked-18.0.11/marked.umd.js                 -> window.marked
   assets/js/marked-katex-extension-5.1.12/index.umd.js   -> window.markedKatex
 and the stylesheet assets/js/katex-0.18.5/katex.min.css.
 
 This file is deliberately plain JavaScript (no Wt macros, no InterSpec globals) so that a future
 "download the conversation as a standalone HTML file" feature can inline it verbatim, along with the
 three libraries above, and get byte-identical rendering offline.
 
 SECURITY: LLM output is not trusted.  It is influenced by spectrum-file content, tool results, and
 whatever the model was trained on, and it ends up as innerHTML inside the InterSpec DOM.  Every byte
 of the emitted HTML comes from either marked's own escaping renderers, KaTeX's renderToString(), or
 filterRawHtml() here, and no <a>, <img>, <script>, <iframe> or <style> is ever produced - so rendered
 output cannot reference an off-machine resource or execute anything.

 The three renderer overrides in configure() below - link, image and html - are the complete set of
 places model text could otherwise reach an attribute value; keep them that way.  Some attributes are
 still emitted, but none of them carry unescaped model text: marked writes class="language-..." on a
 code fence (escaped by marked, so a fence language of `js" onload="x` yields class="language-js&quot;"),
 align on table cells (constrained to left|center|right by the tokenizer), start on ordered lists
 (coerced with unary +), and type/checked/disabled on a task-list checkbox (literal).  KaTeX emits
 only class, style, xmlns and aria-hidden, and is invoked with trust:false, which disables \href,
 \url and \includegraphics.
 */

var LlmMarkdown = (function(){
'use strict';

/** Raw HTML the model emitted is escaped to literal text, except for these tags, and then only when
    written with no attributes at all.  <sub>/<sup> are here because models write isotopes that way,
    e.g. "<sup>235</sup>U". */
const ALLOWED_TAGS = { b:1, i:1, em:1, strong:1, u:1, s:1, del:1,
                       sub:1, sup:1, code:1, kbd:1, mark:1, br:1 };

/** Matches only a bare open/close/self-closing tag - anything with an attribute fails to match, and
    so gets escaped. */
const BARE_TAG_RE = /^<\/?([a-zA-Z][a-zA-Z0-9]*)\s*\/?>$/;

const KATEX_OPTS = { throwOnError: false,  // bad math shows in red rather than aborting the parse
                     trust: false,         // disables \href, \url, \includegraphics
                     strict: "ignore",
                     output: "htmlAndMathml" };

let sm_configured = false;


function escapeHtml( text )
{
  return String(text).replace( /&/g, "&amp;" ).replace( /</g, "&lt;" )
                     .replace( />/g, "&gt;" ).replace( /"/g, "&quot;" )
                     .replace( /'/g, "&#39;" );
}//escapeHtml(...)


/** Like escapeHtml(), but leaves an already-valid entity ("&amp;", "&#65;") alone - the same rule
    marked's own text renderer uses, so ordinary prose renders identically. */
const VALID_ENTITY_AMP = /&(?!(#\d{1,7}|#[Xx][a-fA-F0-9]{1,6}|\w+);)/g;

function escapeText( text )
{
  return String(text).replace( VALID_ENTITY_AMP, "&amp;" ).replace( /</g, "&lt;" )
                     .replace( />/g, "&gt;" ).replace( /"/g, "&quot;" )
                     .replace( /'/g, "&#39;" );
}//escapeText(...)


/** Escapes a raw-HTML chunk to literal text, passing through only the attribute-free tags in
    ALLOWED_TAGS.  Used for both block-level and inline HTML tokens. */
function filterRawHtml( raw )
{
  let out = "", last = 0, m;
  const re = /<[^>]*>/g;

  while( (m = re.exec(raw)) !== null )
  {
    out += escapeHtml( raw.slice(last, m.index) );

    const tag = BARE_TAG_RE.exec( m[0] );
    if( tag && ALLOWED_TAGS[tag[1].toLowerCase()] )
      out += m[0].toLowerCase();
    else
      out += escapeHtml( m[0] );

    last = re.lastIndex;
  }//while( more tags )

  return out + escapeHtml( raw.slice(last) );
}//filterRawHtml(...)


function katexRenderer( displayMode )
{
  return function( token ){
    try
    {
      const opts = Object.assign( {}, KATEX_OPTS, {displayMode: displayMode} );
      return katex.renderToString( token.text, opts );
    }catch( e )
    {
      return escapeHtml( token.raw );  // show the LaTeX source rather than losing the content
    }
  };
}//katexRenderer(...)


/** Builds a marked inline extension that hands `regex`'s first capture group to KaTeX.

    Used to cover the math delimiters the vendored extension does not (see configure()).  Doing this
    as marked extensions, rather than rewriting the Markdown text first, matters: marked tries
    extensions before its built-in rules - including before the backslash-escape rule that would
    otherwise turn "\(" into "(" - but only after fenced blocks and code spans have been consumed, so
    math inside a code sample is left alone without us having to detect code fences ourselves. */
function latexExtension( name, regex, opener, displayMode )
{
  return {
    name: name,
    level: "inline",
    start( src ){
      // Report where the next real match could begin, so marked's text tokenizer stops there.
      let from = 0;
      for( ;; )
      {
        const i = src.indexOf( opener, from );
        if( i < 0 )
          return undefined;
        if( regex.test( src.slice(i) ) )
          return i;
        from = i + 1;
      }
    },
    tokenizer( src ){
      const m = regex.exec( src );
      if( m )
        return { type: name, raw: m[0], text: m[1].trim(), displayMode: displayMode };
    },
    renderer: katexRenderer( displayMode )
  };
}//latexExtension(...)


/** marked-katex-extension's standard $...$ rule requires the closing $ to be followed by whitespace
    or punctuation, so it does not fire on a letter-hugged span - and "$^{235}$U" (how everyone writes
    isotopes) is exactly that.  Its non-standard mode would fire, but then "$5 to $10 to ship" also
    becomes math, which is worse.

    So supplement it: match a single-dollar span the standard rule skipped, but only when the content
    holds a LaTeX control character (one of \ ^ _ { }).  Currency amounts never do, so the two cases
    stay cleanly separated.  Note marked unshifts tokenizers, so this rule is actually tried BEFORE
    the vendored one; that is harmless precisely because of the control-character requirement. */
const DOLLAR_LATEX_RE = /^\$(?!\$)([^$\n]*[\\^_{}][^$\n]*)\$/;


function isLoaded()
{
  return (typeof marked !== "undefined") && (typeof katex !== "undefined")
          && (typeof markedKatex !== "undefined");
}


function configure()
{
  if( sm_configured )
    return;

  marked.use( { gfm: true, breaks: true } );

  // $...$ and $$...$$ (vendored extension, standard delimiter rules - the non-standard mode would
  //  treat things like "$5 to $10" as math).
  marked.use( markedKatex(KATEX_OPTS) );

  // The delimiters the vendored extension leaves on the table: \(...\) and \[...\], which many models
  //  emit instead of dollars, and the letter-hugged $...$ described above.  KaTeX's displayMode output
  //  is a <span>, so the display form is safe to emit from an inline-level extension.
  marked.use( { extensions: [
    latexExtension( "inlineLatexParen",   /^\\\(([\s\S]+?)\\\)/, "\\(", false ),
    latexExtension( "inlineLatexBracket", /^\\\[([\s\S]+?)\\\]/, "\\[", true  ),
    latexExtension( "inlineLatexDollar",  DOLLAR_LATEX_RE,       "$",  false )
  ] } );

  // The allowlist - see the SECURITY note at the top of this file.
  marked.use( { renderer: {
    link: function( token ){
      // Never emit an <a>: show the link text, and append the target as inert text if it adds info.
      const text = token.tokens ? this.parser.parseInline( token.tokens )
                                : escapeHtml( token.text || "" );
      const href = token.href || "";
      if( !href || (text.indexOf(href) >= 0) )
        return text;
      return text + " <span class=\"LlmMarkdownUrl\">(" + escapeHtml(href) + ")</span>";
    },

    image: function( token ){
      // Never emit an <img>: show the alt text and the source as inert text.
      const alt = escapeHtml( token.text || "" );
      const href = token.href || "";
      if( !href )
        return alt;
      return "<span class=\"LlmMarkdownUrl\">[image: " + (alt ? alt + ", " : "") + escapeHtml(href) + "]</span>";
    },

    // Both block-level HTML and inline tags route here.
    html: function( token ){ return filterRawHtml( token.text ); },

    // Overriding `html` alone is NOT enough, and leaving this out is a hole straight through the
    //  allowlist.  marked's tokenizer flips lexer.state.inRawBlock on at <pre>, <code>, <kbd> or
    //  <script> - it does so from the raw source, so it happens even for a tag this renderer escapes,
    //  and it stays on for the rest of the message - and while it is on, inlineText stamps text
    //  tokens with escaped:true, which makes marked's default text renderer emit them verbatim.
    //  "<kbd>a <img src=x onerror=... </kbd>" then reaches the DOM as a live <img>: the unterminated
    //  tag never becomes an html token (marked's tag rule needs the closing >), so it rides through
    //  as raw text and the browser's parser closes it on the next > it finds.  Escape unconditionally.
    text: function( token ){
      return (token.tokens && token.tokens.length) ? this.parser.parseInline( token.tokens )
                                                   : escapeText( token.text );
    },

    // A GFM task list would otherwise emit a real <input type="checkbox">.
    checkbox: function( token ){ return token.checked ? "[x] " : "[ ] "; },

    // The fence info string becomes class="language-<model text>".  marked escapes the quotes, so it
    //  cannot break out of the attribute, but it still lets the model name any InterSpec CSS class.
    code: function( token ){
      const lang = String(token.lang || "").match( /^[A-Za-z0-9_+-]{0,32}/ )[0];
      const cls = lang ? " class=\"language-" + lang + "\"" : "";
      return "<pre><code" + cls + ">" + escapeText( token.text ) + "\n</code></pre>\n";
    }
  } } );

  sm_configured = true;
}//configure()


/** Renders Markdown to a sanitized HTML string.  Throws if the libraries are not loaded. */
function toHtml( markdown )
{
  configure();
  return marked.parse( markdown || "" );
}//toHtml(...)


// The libraries are pulled in with WApplication::require(), so on the very first render they may not
//  have finished loading yet; queue up until they have, then give up and show plain text.
const sm_pending = [];
let sm_pollTimer = null, sm_pollCount = 0;

function whenReady( fcn )
{
  if( isLoaded() )
  {
    fcn( false );
    return;
  }

  sm_pending.push( fcn );

  if( sm_pollTimer )
    return;

  sm_pollCount = 0;  // each wait gets its own five seconds
  sm_pollTimer = setInterval( function(){
    sm_pollCount += 1;
    const failed = (sm_pollCount > 100);  // ~5 seconds
    if( !isLoaded() && !failed )
      return;

    clearInterval( sm_pollTimer );
    sm_pollTimer = null;
    const waiting = sm_pending.splice( 0, sm_pending.length );
    for( const f of waiting )
      f( failed );
  }, 50 );
}//whenReady(...)


function showPlain( el )
{
  el.textContent = el.llmMarkdownSource || "";
  el.classList.add( "LlmMarkdownRaw" );
}


function showRendered( el )
{
  try
  {
    el.innerHTML = toHtml( el.llmMarkdownSource );
    el.classList.remove( "LlmMarkdownRaw" );
  }catch( e )
  {
    console.warn( "LlmMarkdown: failed to render, falling back to plain text:", e );
    showPlain( el );
  }
}//showRendered(...)


/** Renders `markdown` into the element with id `elemId`.  The Markdown source is kept on the element
    (llmMarkdownSource) for the raw/rendered toggle and for conversation export. */
function render( elemId, markdown )
{
  const el = document.getElementById( elemId );
  if( !el )
    return;

  el.llmMarkdownSource = markdown || "";

  whenReady( function( failed ){
    if( failed || el.llmMarkdownRaw )
      showPlain( el );
    else
      showRendered( el );
  } );
}//render(...)


/** Toggles the element between rendered HTML and the raw Markdown source. */
function setRawMode( elemId, raw )
{
  const el = document.getElementById( elemId );
  if( !el )
    return;

  el.llmMarkdownRaw = !!raw;

  if( raw )
    showPlain( el );
  else
    whenReady( function( failed ){ failed ? showPlain(el) : showRendered(el); } );
}//setRawMode(...)


return { render: render, setRawMode: setRawMode, toHtml: toHtml, isLoaded: isLoaded };
})();
