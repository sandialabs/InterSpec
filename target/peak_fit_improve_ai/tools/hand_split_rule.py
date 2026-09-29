#!/usr/bin/env python3
"""From the hand fits: for neighbouring significant truth lines, did the hand fit keep them in one ROI
or separate them?  Tabulate against separation in FWHM and the leak at the best split point."""
import sys, os, json, glob, math
run = sys.argv[1]; minz = float( sys.argv[2] ) if len( sys.argv ) > 2 else 5.0
def Phi( x ): return 0.5 * ( 1 + math.erf( x / math.sqrt( 2 ) ) )
rows = []
for path in sorted( glob.glob( os.path.join( run, 'plot_data', '*.json' ) ) ):
  d = json.load( open( path ) )
  hand = sorted( ( r['lower'], r['upper'] ) for r in d.get( 'ref', [] ) )
  tl = sorted( [t for t in d.get( 'truth', [] ) if t.get( 'signal' ) and t['z'] >= minz and t['e'] >= 30], key = lambda t: t['e'] )
  def roi_of( e ):
    for k, h in enumerate( hand ):
      if h[0] <= e <= h[1]: return k
    return None
  x = d['x']; y = d['y']
  for a, b in zip( tl, tl[1:] ):
    ra, rb = roi_of( a['e'] ), roi_of( b['e'] )
    if ra is None or rb is None: continue
    fw = 0.5 * ( a['fwhm'] + b['fwhm'] )
    sep = ( b['e'] - a['e'] ) / fw
    sa, sb = a['fwhm'] / 2.3548, b['fwhm'] / 2.3548
    # best split point: minimise total leak over a grid
    best = None
    for i in range( 1, 200 ):
      c = a['e'] + ( b['e'] - a['e'] ) * i / 200.0
      la = a['area'] * ( 1 - Phi( ( c - a['e'] ) / sa ) )
      lb = b['area'] * Phi( ( c - b['e'] ) / sb )
      if best is None or la + lb < best[0]: best = ( la + lb, c, la, lb )
    leak, c, la, lb = best
    # data counts within +-0.5 FWHM of the split point
    import bisect
    i0 = max( 0, bisect.bisect_left( x, c - 0.5 * fw ) - 1 ); i1 = bisect.bisect_left( x, c + 0.5 * fw )
    n = sum( y[i0:i1] )
    rows.append( ( d['id'], a['e'], b['e'], sep, leak / math.sqrt( max( n, 1 ) ), min( a['area'], b['area'] ) / max( a['area'], b['area'] ), ra == rb ) )
rows.sort( key = lambda r: r[3] )
print( '%-16s %7s %7s %6s %8s %6s %s' % ( 'id', 'E1', 'E2', 'sepF', 'leak/sd', 'ratio', 'hand' ) )
for r in rows:
  print( '%-16s %7.1f %7.1f %6.2f %8.2f %6.3f %s' % ( r[0], r[1], r[2], r[3], r[4], r[5], 'SAME' if r[6] else 'split' ) )
