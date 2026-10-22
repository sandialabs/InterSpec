#!/usr/bin/env python3
"""review_page.py PART_GLOB IMG_DIR OUT_HTML TITLE : spectra grouped by verdict of the new fit, each with
its defect list, the before/after call, and the before/after image; filters by class and comparison."""
import sys, glob, os, html, collections
parts, img_dir, out, title = sys.argv[1:5]
verdict = {}; defects = collections.defaultdict( list ); cmp = {}
for f in sorted( glob.glob( parts ) ):
  for line in open( f ):
    c = line.rstrip( '\n' ).split( '\t' )
    if len( c ) < 4: continue
    if c[1] == 'CMP':
      cmp[c[0]] = ( c[2].upper(), c[3], '\t'.join( c[4:] ) )
    elif len( c ) >= 6:
      verdict[c[0]] = c[1].upper()
      if c[2].upper() != 'NONE': defects[c[0]].append( ( c[2].upper(), c[3], c[4], '\t'.join( c[5:] ) ) )
order = ['CATASTROPHIC', 'BAD', 'MINOR', 'GOOD']
classes = sorted( { d[0] for v in defects.values() for d in v } )
cmps = ['BETTER', 'WORSE', 'MIXED', 'SAME']
vc = collections.Counter( verdict.values() ); cc = collections.Counter( v[0] for v in cmp.values() )
h = ['<!doctype html><html><head><meta charset="utf-8"><title>%s</title><style>' % html.escape( title ),
     ':root{--bg:#fff;--fg:#222;--mut:#666;--card:#f6f6f6;--acc:#b30000;--good:#1a7f37;--bad:#b30000}',
     '@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#1b1b1b;--fg:#ddd;--mut:#999;--card:#262626;--acc:#ff7070;--good:#5fd07a;--bad:#ff7070}}',
     ':root[data-theme="dark"]{--bg:#1b1b1b;--fg:#ddd;--mut:#999;--card:#262626;--acc:#ff7070;--good:#5fd07a;--bad:#ff7070}',
     'body{background:var(--bg);color:var(--fg);font:14px system-ui,sans-serif;margin:0 16px}',
     '.card{background:var(--card);border-radius:6px;padding:8px 10px;margin:12px 0}',
     '.card img{width:100%;height:auto;display:block;margin-top:6px}',
     'h2{margin:24px 0 4px} .v{font-weight:600;color:var(--acc)} .c{font-weight:600} .BETTER{color:var(--good)} .WORSE{color:var(--bad)}',
     'table{border-collapse:collapse;font-size:13px} td{padding:1px 8px 1px 0;vertical-align:top}',
     '.bar{position:sticky;top:0;background:var(--bg);padding:8px 0;z-index:2} .bar label{margin-right:10px;white-space:nowrap}',
     '</style></head><body><h1>%s</h1>' % html.escape( title ),
     '<p>New fit: %s. &nbsp; New vs old: %s. &nbsp; Images: RED = new fit, BLUE = old fit.</p>' % (
        ', '.join( '%s %d' % ( v, vc.get( v, 0 ) ) for v in reversed( order ) ), ', '.join( '%s %d' % ( v, cc.get( v, 0 ) ) for v in cmps ) ),
     '<div class="bar">Class: ' + ' '.join( '<label><input type="checkbox" class="cf" value="%s">%s</label>' % ( c, c ) for c in classes ) +
     ' &nbsp;|&nbsp; New vs old: ' + ' '.join( '<label><input type="checkbox" class="mf" value="%s">%s</label>' % ( c, c ) for c in cmps ) + '</div>']
for v in order:
  ids = sorted( k for k, x in verdict.items() if x == v )
  h.append( '<h2>%s (%d)</h2>' % ( v, len( ids ) ) )
  for pid in ids:
    ds = defects.get( pid, [] ); cm = cmp.get( pid, ( '', '', '' ) )
    h.append( '<div class="card" data-cls="%s" data-cmp="%s"><span class="v">%s</span> &nbsp;<b>%s</b> &nbsp; <span class="c %s">%s</span> %s' % (
      ' '.join( sorted( { d[0] for d in ds } ) ), cm[0], v, html.escape( pid ), cm[0], cm[0], html.escape( cm[2] ) ) )
    if ds:
      h.append( '<table>' + ''.join( '<tr><td>%s</td><td>%s</td><td>sev %s</td><td>%s</td></tr>' % tuple( html.escape( x ) for x in d ) for d in ds ) + '</table>' )
    h.append( '<img loading="lazy" src="%s">' % html.escape( os.path.join( img_dir, pid + '.png' ) ) )
    h.append( '</div>' )
h.append( '''<script>
function upd(){const cs=[...document.querySelectorAll('.cf:checked')].map(e=>e.value);const ms=[...document.querySelectorAll('.mf:checked')].map(e=>e.value);
document.querySelectorAll('.card').forEach(c=>{const cl=c.dataset.cls.split(' ');const okc=!cs.length||cs.some(x=>cl.includes(x));const okm=!ms.length||ms.includes(c.dataset.cmp);c.style.display=(okc&&okm)?'':'none';});}
document.querySelectorAll('.cf,.mf').forEach(e=>e.addEventListener('change',upd));</script></body></html>''' )
open( out, 'w' ).write( '\n'.join( h ) )
print( 'wrote', out, len( verdict ), 'spectra' )
