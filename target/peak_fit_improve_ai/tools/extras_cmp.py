#!/usr/bin/env python3
"""extras_cmp.py RUN_A RUN_B : significant extra peaks (verdict extra = matches nothing, extra_bkg =
a background line) in each run, and the ones B has that A lacks (no A extra within 10 keV), strongest
first - the list to look at when a change finds more truth lines but also more peaks."""
import sys, csv, collections, os

def extras( run ):
  out = collections.defaultdict( list )
  for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter='\t' ):
    if r['set'] == 'fitted' and r['verdict'] in ( 'extra', 'extra_bkg' ):
      out[r['id']].append( ( round( float( r['energy'] ), 1 ), r['verdict'], round( float( r['z_det'] ), 1 ) ) )
  return out

a, b = extras( sys.argv[1] ), extras( sys.argv[2] )
count = lambda ex, kind: sum( 1 for peaks in ex.values() for p in peaks if p[1] == kind )
print( 'extra: %d -> %d   extra_bkg: %d -> %d'
       % ( count( a, 'extra' ), count( b, 'extra' ), count( a, 'extra_bkg' ), count( b, 'extra_bkg' ) ) )
new = [ (pid,) + p for pid, peaks in b.items() for p in peaks
        if not any( abs( p[0] - q[0] ) < 10 for q in a.get( pid, [] ) ) ]
for pid, energy, verdict, z in sorted( new, key = lambda t: -t[3] ):
  print( '  new in B: %-24s %8.1f keV  %-9s z=%.1f' % ( pid, energy, verdict, z ) )
