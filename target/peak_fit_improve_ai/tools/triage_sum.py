#!/usr/bin/env python3
"""Merge the visual-triage TSVs: verdict counts, defects per class (by severity), worst spectra."""
import sys, glob, collections
files = sys.argv[1:]
verdict = {}
defects = []
for f in files:
  for line in open( f ):
    c = line.rstrip( '\n' ).split( '\t' )
    if len( c ) < 6 or c[0] in ( 'id', '' ): continue
    pid, ov, cls, roi, sev, note = c[0], c[1].upper(), c[2].upper(), c[3], c[4], '\t'.join( c[5:] )
    verdict[pid] = ov
    if cls != 'NONE':
      try: s = int( float( sev ) )
      except: s = 1
      defects.append( ( pid, ov, cls, roi, s, note ) )
print( 'spectra judged: %d' % len( verdict ) )
for v, n in collections.Counter( verdict.values() ).most_common(): print( '  %-13s %d' % ( v, n ) )
by = collections.defaultdict( lambda: [0, 0, 0, set()] )
for pid, ov, cls, roi, s, note in defects:
  by[cls][s - 1] += 1; by[cls][3].add( pid )
print( '\n%-10s %5s %5s %5s %8s' % ( 'class', 'sev3', 'sev2', 'sev1', 'spectra' ) )
for cls, v in sorted( by.items(), key = lambda kv: -( 3 * kv[1][2] + 2 * kv[1][1] + kv[1][0] ) ):
  print( '%-10s %5d %5d %5d %8d' % ( cls, v[2], v[1], v[0], len( v[3] ) ) )
