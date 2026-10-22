#!/usr/bin/env python3
"""triage_page.py TRIAGE_GLOB IMG_DIR RUN_DIR OUT_HTML : one page, spectra grouped by verdict, each with
its defect list, spectrum path and whole-spectrum image; a class filter at the top."""
import sys, glob, json, os, html, collections
parts, img_dir, run_dir, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
verdict = {}; defects = collections.defaultdict( list )
for f in sorted( glob.glob( parts ) ):
  for line in open( f ):
    c = line.rstrip( '\n' ).split( '\t' )
    if len( c ) < 6: continue
    verdict[c[0]] = c[1].upper()
    if c[2].upper() != 'NONE': defects[c[0]].append( ( c[2].upper(), c[3], c[4], '\t'.join( c[5:] ) ) )
corpus = sys.argv[5] if len( sys.argv ) > 5 else ''
paths = { pid: os.path.join( corpus, pid + '.pcf' ) for pid in verdict } if corpus else {}
order = ['CATASTROPHIC', 'BAD', 'MINOR', 'GOOD']
classes = sorted( { d[0] for v in defects.values() for d in v } )
h = ['<!doctype html><html><head><meta charset="utf-8"><title>NaI Fit Triage</title><style>',
     ':root{--bg:#fff;--fg:#222;--mut:#666;--card:#f6f6f6;--acc:#b30000}',
     '@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#1b1b1b;--fg:#ddd;--mut:#999;--card:#262626;--acc:#ff7070}}',
     ':root[data-theme="dark"]{--bg:#1b1b1b;--fg:#ddd;--mut:#999;--card:#262626;--acc:#ff7070}',
     'body{background:var(--bg);color:var(--fg);font:14px system-ui,sans-serif;margin:0 16px}',
     '.card{background:var(--card);border-radius:6px;padding:8px 10px;margin:12px 0}',
     '.card img{width:100%;height:auto;display:block;margin-top:6px}',
     'h2{margin:24px 0 4px} .v{font-weight:600;color:var(--acc)} .p{color:var(--mut);font-family:monospace;font-size:12px;word-break:break-all}',
     'table{border-collapse:collapse;font-size:13px} td{padding:1px 8px 1px 0;vertical-align:top}',
     '.bar{position:sticky;top:0;background:var(--bg);padding:8px 0;z-index:2} .bar label{margin-right:10px;white-space:nowrap}',
     '</style></head><body><h1>NaI fit triage: IdentiFINDER-R500</h1>',
     '<div class="bar">Show spectra with class: ' + ' '.join( '<label><input type="checkbox" class="cf" value="%s">%s</label>' % ( c, c ) for c in classes ) +
     ' <label><input type="checkbox" id="sev3">severity 3 only</label></div>']
for v in order:
  ids = sorted( k for k, x in verdict.items() if x == v )
  h.append( '<h2>%s (%d)</h2>' % ( v, len( ids ) ) )
  for pid in ids:
    ds = defects.get( pid, [] )
    cls = ' '.join( sorted( { d[0] for d in ds } ) )
    s3 = ' '.join( sorted( { d[0] for d in ds if d[2].startswith( '3' ) } ) )
    h.append( '<div class="card" data-cls="%s" data-s3="%s"><span class="v">%s</span> &nbsp;<b>%s</b><div class="p">%s</div>' % (
      cls, s3, v, html.escape( pid ), html.escape( paths.get( pid, '' ) ) ) )
    if ds:
      h.append( '<table>' + ''.join( '<tr><td>%s</td><td>%s</td><td>sev %s</td><td>%s</td></tr>' % tuple( html.escape( x ) for x in d ) for d in ds ) + '</table>' )
    h.append( '<img loading="lazy" src="%s" alt="%s">' % ( html.escape( os.path.relpath( os.path.join( img_dir, pid + '.png' ), os.path.dirname( out ) ) ), html.escape( pid ) ) )
    h.append( '</div>' )
h.append( '''<script>
function upd(){const on=[...document.querySelectorAll('.cf:checked')].map(e=>e.value);const s3=document.getElementById('sev3').checked;
document.querySelectorAll('.card').forEach(c=>{const cl=(s3?c.dataset.s3:c.dataset.cls).split(' ');c.style.display=(on.length==0&&!s3)||(on.length==0?cl[0]!='':on.some(x=>cl.includes(x)))?'':'none';});}
document.querySelectorAll('input').forEach(e=>e.addEventListener('change',upd));
</script></body></html>''' )
open( out, 'w' ).write( '\n'.join( h ) )
print( 'wrote', out )
