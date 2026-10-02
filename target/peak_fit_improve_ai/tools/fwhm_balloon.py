#!/usr/bin/env python3
"""fwhm_balloon.py RUN [RUN_B] : how far each problem's accepted solve widths ran from the planner's width
model, from roi_plan_trace.txt (the first `fwhm model:` line, fwhm at 60/122/662/1332 keV, against the last
`solve fwhm:` line).  Prints, per energy band, how many solves exceed 1.3x / 1.6x the model or fall under
0.7x.  With RUN_B, also lists the problems whose worst low-energy ratio moved by more than 0.2.
SAM-Eagle's 12.5 keV channels hold every solve at a channel-width floor there, so its low bands read high
by design."""
import sys, re, math

BANDS = [ ( 55, 100 ), ( 100, 300 ), ( 300, 3000 ) ]


def load( run ):
  model, solve, prob = {}, {}, None
  for line in open( run + '/roi_plan_trace.txt' ):
    m = re.match( r'^==== (\S+)', line )
    if m:
      prob = m.group( 1 )
      continue
    m = re.search( r'fwhm model: .*?fwhm\(60\)=([0-9.]+) fwhm\(122\)=([0-9.]+) fwhm\(662\)=([0-9.]+) fwhm\(1332\)=([0-9.]+)', line )
    if m and prob not in model:
      model[prob] = [ float( x ) for x in m.groups() ]
    if line.startswith( '  solve fwhm:' ):
      solve[prob] = [ tuple( map( float, p.split( '=' ) ) ) for p in line.split( ':', 1 )[1].split() ]
  return model, solve


def model_at( fwhms, energy ):
  energies = [ 60.0, 122.0, 662.0, 1332.0 ]
  le = math.log( min( max( energy, 20.0 ), 3000.0 ) )
  xs = [ math.log( e ) for e in energies ]
  ys = [ math.log( f ) for f in fwhms ]
  for i in range( 3 ):
    if le <= xs[i + 1] or i == 2:
      return math.exp( ys[i] + ( ys[i + 1] - ys[i] ) / ( xs[i + 1] - xs[i] ) * ( le - xs[i] ) )


def worst_ratio( model, solve, prob, lo, hi ):
  if prob not in model or prob not in solve:
    return None
  ratios = [ w / model_at( model[prob], e ) for e, w in solve[prob] if lo <= e <= hi ]
  return max( ratios, key = lambda r: abs( math.log( r ) ) ) if ratios else None


def summary( run ):
  model, solve = load( run )
  print( run )
  for lo, hi in BANDS:
    rs = [ r for r in ( worst_ratio( model, solve, p, lo, hi ) for p in solve ) if r ]
    print( '  %4d-%-4d keV: %3d solves, >1.3x %3d, >1.6x %3d, <0.7x %3d' % ( lo, hi, len( rs ), sum( r > 1.3 for r in rs ),
           sum( r > 1.6 for r in rs ), sum( r < 0.7 for r in rs ) ) )
  return model, solve


a = summary( sys.argv[1] )
if len( sys.argv ) > 2:
  b = summary( sys.argv[2] )
  for p in sorted( set( a[1] ) & set( b[1] ) ):
    ra, rb = worst_ratio( *a, p, 55, 150 ), worst_ratio( *b, p, 55, 150 )
    if ra and rb and abs( ra - rb ) > 0.2:
      print( '  %-26s 55-150 keV worst ratio %.2f -> %.2f' % ( p, ra, rb ) )
