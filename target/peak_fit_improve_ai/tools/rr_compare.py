#!/usr/bin/env python3
"""rr_compare.py 'OLD_GLOB' 'NEW_GLOB' : compare two visual reviews of the same spectra (TSVs in the
RUBRIC_BOTH.md format), e.g. a full review and a later re-review of the spectra it flagged.  Prints the
verdict and before/after (CMP) transitions over the spectra the NEW review covers, and the NEW review's
remaining WORSE notes.  Quote globs so the shell does not expand them."""
import glob, collections, sys
def load( pattern ):
  verdict, cmp_, note = {}, {}, {}
  for f in glob.glob( pattern ):
    for line in open( f ):
      c = line.rstrip( '\n' ).split( '\t' )
      if len( c ) < 3: continue
      if c[1].upper() == 'CMP':
        cmp_[c[0]] = c[2].upper(); note[c[0]] = c[4] if len( c ) > 4 else ''
      elif c[1].upper() in ( 'GOOD', 'MINOR', 'BAD', 'CATASTROPHIC' ):
        verdict[c[0]] = c[1].upper()
  return verdict, cmp_, note
vo, co, no = load( sys.argv[1] )
vn, cn, nn = load( sys.argv[2] )
ids = sorted( set( cn ) )
tv = collections.Counter( ( vo.get( i, '?' ), vn.get( i, '?' ) ) for i in ids )
tc = collections.Counter( ( co.get( i, '?' ), cn.get( i, '?' ) ) for i in ids )
print( 'spectra:', len( ids ) )
print( 'verdict old -> new:' )
for k, n in sorted( tv.items(), key = lambda t: -t[1] ): print( '  %-8s -> %-8s %d' % ( k[0], k[1], n ) )
print( 'CMP old -> new:' )
for k, n in sorted( tc.items(), key = lambda t: -t[1] ): print( '  %-7s -> %-7s %d' % ( k[0], k[1], n ) )
print( 'old totals:', dict( collections.Counter( vo.get( i, '?' ) for i in ids ) ), dict( collections.Counter( co.get( i, '?' ) for i in ids ) ) )
print( 'new totals:', dict( collections.Counter( vn.get( i, '?' ) for i in ids ) ), dict( collections.Counter( cn.get( i, '?' ) for i in ids ) ) )
print( '\nstill WORSE in the new review:' )
for i in ids:
  if cn.get( i ) == 'WORSE': print( '  %-22s %s' % ( i, nn[i][:170] ) )
