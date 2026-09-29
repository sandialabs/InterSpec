#!/usr/bin/env python3
"""bitident.py RUN_A RUN_B : are the fitted peaks of two runs identical (energy, amplitude, FWHM, ROI
bounds for every peak of every problem)?  The guard for changes that must not move HPGe: equal summary
lines can hide compensating differences."""
import sys, csv, os
def load( run ):
  out = {}
  for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
    if r['set'] == 'fitted':
      out.setdefault( r['id'], [] ).append( ( r['energy'], r['amplitude'], r['fwhm'], r['roi_lower'], r['roi_upper'] ) )
  return out
a, b = load( sys.argv[1] ), load( sys.argv[2] )
diff = [k for k in sorted( set( a ) | set( b ) ) if a.get( k ) != b.get( k )]
print( '%s vs %s: %d problems differ%s' % ( os.path.basename( sys.argv[1] ), os.path.basename( sys.argv[2] ), len( diff ),
       ( ': ' + ' '.join( diff[:20] ) ) if diff else ' (bit-identical)' ) )
sys.exit( 1 if diff else 0 )
