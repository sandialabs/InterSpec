#!/usr/bin/env python3
"""turnon_rois.py RUN [v] : delivered ROIs that start at the detector turn-on - lower edge within
1.5 FWHM of the planning floor (the first data-alive energy) or below the spectroscopic extent - with
their peaks' verdicts and whether truth has a line (z>=3) inside.  Classes: REAL (a z>=3 truth line
matched in the ROI), NONE (no z>=3 truth line inside the ROI)."""
import sys, csv, os, re, collections
run = sys.argv[1]; verbose = len( sys.argv ) > 2
floor, extent = {}, {}
prob = None
for line in open( os.path.join( run, 'roi_plan_trace.txt' ) ):
  if line.startswith( '==== ' ):
    prob = line.split()[1]; continue
  m = re.search( r'planning window ([0-9.]+)-([0-9.]+) keV \(sub-extent below ([0-9.]+) keV\)', line )
  if m and prob and prob not in floor:
    floor[prob] = float( m.group(1) ); extent[prob] = float( m.group(3) )
peaks = collections.defaultdict( list ); truth = collections.defaultdict( list )
for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
  if r['set'] == 'fitted':
    peaks[r['id']].append( r )
  elif r['set'] == 'truth':
    truth[r['id']].append( ( float( r['energy'] ), float( r['fwhm'] ), float( r['truth_z'] or 0 ), r['verdict'] ) )
cls = collections.Counter(); rows = []
for pid, pk in peaks.items():
  if pid not in floor: continue
  rois = collections.defaultdict( list )
  for r in pk:
    rois[( float( r['roi_lower'] ), float( r['roi_upper'] ) )].append( r )
  lowest = min( rois ) if rois else None
  if not lowest: continue
  lo, hi = lowest
  fwhm_lo = min( float( r['fwhm'] ) for r in rois[lowest] )
  if not ( lo <= floor[pid] + 1.5*fwhm_lo or lo < extent[pid] ): continue
  real = [ t for t in truth[pid] if lo <= t[0] <= hi and t[2] >= 3 ]
  real_found = [ t for t in real if t[3].startswith( 'matched' ) ]
  c = 'REAL' if real_found else ( 'REAL-missed' if real else 'NONE' )
  cls[c] += 1
  rows.append( ( pid, floor[pid], extent[pid], lo, hi, c, ' '.join( '%.1f(%s,z%.0f)' % ( float( r['energy'] ), r['verdict'][:7], float( r['z_det'] or 0 ) ) for r in sorted( rois[lowest], key = lambda r: float( r['energy'] ) ) ),
                 ' '.join( '%.1f/z%.0f' % ( t[0], t[2] ) for t in real ) ) )
print( run, dict( cls ) )
if verbose:
  for r in sorted( rows, key = lambda r: ( r[5], r[0] ) ):
    print( '  %-22s floor %5.1f ext %5.1f  ROI %6.1f-%6.1f %-11s peaks: %s  truth: %s' % r )
