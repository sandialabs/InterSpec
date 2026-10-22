#!/usr/bin/env python3
"""side_extent.py HAND_RUN FIT_RUN [N] : per-side ROI extent beyond the outermost peak, in FWHM, for the
hand ROIs (HAND_RUN plot_data 'ref', a run_hand.sh run) and for the fitted ROIs of FIT_RUN on the same
spectra.  The FWHM is the truth line's at that energy (nearest truth line), so both are measured with
the same ruler.  Sides that face a neighbouring ROI within 0.5 FWHM are 'shared' and left out of the
quantiles.  With N, also list the N fitted ROIs with the longest free sides."""
import sys, os, json, glob, statistics
hand_run, fit_run = sys.argv[1], sys.argv[2]
def truth_fwhm( d, e ):
  best = None
  for t in d.get( 'truth', [] ):
    if best is None or abs( t['e'] - e ) < abs( best['e'] - e ):
      best = t
  return best['fwhm'] if best else None
def sides( d, rois ):
  out = []
  rois = sorted( rois, key = lambda r: r['lower'] )
  for i, r in enumerate( rois ):
    pk = [ p for p in r.get( 'peaks', [] ) ]
    if not pk:
      continue
    lo_p = min( pk, key = lambda p: p['mean'] ); hi_p = max( pk, key = lambda p: p['mean'] )
    fl = truth_fwhm( d, lo_p['mean'] ) or lo_p['fwhm']; fh = truth_fwhm( d, hi_p['mean'] ) or hi_p['fwhm']
    lo = ( lo_p['mean'] - r['lower'] ) / fl; hi = ( r['upper'] - hi_p['mean'] ) / fh
    lo_shared = i > 0 and ( r['lower'] - rois[i-1]['upper'] ) < 0.5*fl
    hi_shared = i + 1 < len( rois ) and ( rois[i+1]['lower'] - r['upper'] ) < 0.5*fh
    out.append( ( d['id'], r['lower'], r['upper'], len( pk ), lo, hi, lo_shared, hi_shared, lo_p['mean'], hi_p['mean'] ) )
  return out
def q( v ):
  v = sorted( v )
  if not v: return 'n=0'
  f = lambda p: v[ min( len( v ) - 1, int( p*len( v ) ) ) ]
  return 'n=%3d q10 %.2f q25 %.2f q50 %.2f q75 %.2f q90 %.2f max %.2f' % ( len( v ), f(.1), f(.25), f(.5), f(.75), f(.9), v[-1] )
hand_sides, fit_sides = [], []
for path in sorted( glob.glob( os.path.join( hand_run, 'plot_data', '*.json' ) ) ):
  d = json.load( open( path ) )
  fp = os.path.join( fit_run, 'plot_data', os.path.basename( path ) )
  if not os.path.exists( fp ):
    continue
  hand_sides += sides( d, [ r for r in d['ref'] if r['upper'] > 20 ] )
  fd = json.load( open( fp ) )
  fd['truth'] = fd.get( 'truth' ) or d.get( 'truth' )
  fit_sides += sides( fd, fd['fit'] )
for name, s in ( ( 'hand', hand_sides ), ( 'fit', fit_sides ) ):
  free = [ x[4] for x in s if not x[6] ] + [ x[5] for x in s if not x[7] ]
  free1 = [ x[4] for x in s if not x[6] and x[3] == 1 ] + [ x[5] for x in s if not x[7] and x[3] == 1 ]
  lo_free = [ x[4] for x in s if not x[6] ]; hi_free = [ x[5] for x in s if not x[7] ]
  print( '%-5s free sides      %s' % ( name, q( free ) ) )
  print( '%-5s  1-peak ROIs    %s' % ( name, q( free1 ) ) )
  print( '%-5s  low side       %s' % ( name, q( lo_free ) ) )
  print( '%-5s  high side      %s' % ( name, q( hi_free ) ) )
if len( sys.argv ) > 3:
  for x in sorted( fit_sides, key = lambda x: -max( x[4] if not x[6] else 0, x[5] if not x[7] else 0 ) )[:int( sys.argv[3] )]:
    print( '%-22s %7.1f-%7.1f n=%d  low %.2f%s  high %.2f%s' % ( x[0], x[1], x[2], x[3], x[4], ' (shared)' if x[6] else '', x[5], ' (shared)' if x[7] else '' ) )
