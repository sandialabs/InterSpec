#!/usr/bin/env python3
"""gate_rej.py RUN : for the LAST plan of each problem, the groups rejected by the keep gate
("z X < Y"), binned by predicted z, with whether a truth line (and its truth z) lies inside the group
(+- 0.5 FWHM).  Tells how many gate rejections are real lines vs phantoms."""
import sys, csv, os, re, collections
run = sys.argv[1]
plans = {}
prob = None
for line in open( os.path.join( run, 'roi_plan_trace.txt' ) ):
  if line.startswith( '==== ' ):
    prob = line.split()[1]
    plans[prob] = []
    continue
  if prob is None:
    continue
  if line.startswith( '  planning window' ):
    plans[prob] = []
    continue
  m = re.match( r'  rejected group ([0-9.]+)-([0-9.]+) keV \((\d+) lines, dominant ([0-9.]+) keV\): S=([0-9.eE+-]+) B=([0-9.eE+-]+) z=([0-9.eE+-]+) - z ', line )
  if m:
    plans[prob].append( dict( lo=float(m.group(1)), hi=float(m.group(2)), n=int(m.group(3)), dom=float(m.group(4)),
                              S=float(m.group(5)), B=float(m.group(6)), z=float(m.group(7)) ) )
truth = collections.defaultdict( list )
for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
  if r['set'] != 'truth':
    continue
  truth[r['id']].append( ( float( r['energy'] ), float( r['fwhm'] ), float( r['truth_z'] or 0 ), r['verdict'] ) )
bins = [ 0, 1, 2, 3, 4, 5, 6.39 ]
stats = collections.OrderedDict( (b, collections.Counter()) for b in bins[:-1] )
rows = []
for pid, rej in plans.items():
  for g in rej:
    best = None
    for e, f, tz, v in truth.get( pid, [] ):
      if g['lo'] - 0.5*f <= e <= g['hi'] + 0.5*f:
        if best is None or tz > best[2]:
          best = ( e, f, tz, v )
    b = max( [ x for x in bins[:-1] if g['z'] >= x ] )
    tz = best[2] if best else 0.0
    cls = 'truth>=3' if tz >= 3 else ( 'truth<3' if best else 'phantom' )
    stats[b][cls] += 1
    if best and tz >= 3 and best[3] == 'missed':
      stats[b]['missed>=3'] += 1
    rows.append( ( pid, g, best ) )
print( 'pred-z bin   n   phantom  truth<3  truth>=3  (missed>=3)' )
for b, c in stats.items():
  n = sum( c[k] for k in ( 'phantom', 'truth<3', 'truth>=3' ) )
  print( '%5.1f+   %4d   %5d    %5d    %5d     %5d' % ( b, n, c['phantom'], c['truth<3'], c['truth>=3'], c['missed>=3'] ) )
if len( sys.argv ) > 2:
  for pid, g, best in rows:
    if g['z'] >= float( sys.argv[2] ):
      print( '%-24s %7.1f-%7.1f dom %7.1f S=%7.0f B=%8.0f z=%5.2f  %s' % ( pid, g['lo'], g['hi'], g['dom'], g['S'], g['B'], g['z'],
             ( 'truth %.1f z=%.1f %s' % ( best[0], best[2], best[3] ) ) if best else 'PHANTOM' ) )
