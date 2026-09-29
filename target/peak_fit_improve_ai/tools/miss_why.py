#!/usr/bin/env python3
"""miss_why.py RUN [minz] [maxz] : why each missed truth line (truth_z in [minz, maxz)) was missed,
from the run's roi_plan_trace.txt.  Categories: planned (an ROI of the LAST plan covered it, so the
line was lost in the solve/refit/filters), rejected:<reason> (a rejected group of the last plan
covered it), none (no group of the last plan came near it)."""
import sys, csv, os, re, collections
run = sys.argv[1]
minz = float( sys.argv[2] ) if len( sys.argv ) > 2 else 8.0
maxz = float( sys.argv[3] ) if len( sys.argv ) > 3 else 1e99

# Last plan per problem: the trace restarts a plan at each "planning window" line.
plans = {}
prob = None
for line in open( os.path.join( run, 'roi_plan_trace.txt' ) ):
  if line.startswith( '==== ' ):
    prob = line.split()[1]
    plans[prob] = { 'rois': [], 'rej': [] }
    continue
  if prob is None:
    continue
  if line.startswith( '  planning window' ):
    plans[prob] = { 'rois': [], 'rej': [] }
    continue
  m = re.match( r'  ROI ([0-9.]+)-([0-9.]+) keV', line )
  if m:
    plans[prob]['rois'].append( ( float( m.group(1) ), float( m.group(2) ) ) )
    continue
  m = re.match( r'  rejected group ([0-9.]+)-([0-9.]+) keV .*? - (.*)$', line )
  if m:
    reason = re.sub( r'[0-9.]+', '#', m.group(3) )[:70]
    plans[prob]['rej'].append( ( float( m.group(1) ), float( m.group(2) ), reason ) )

cats = collections.Counter()
rows = []
for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
  if r['set'] != 'truth' or r['verdict'] != 'missed':
    continue
  z = float( r['truth_z'] or 0 )
  if not ( minz <= z < maxz ):
    continue
  e, fwhm, pid = float( r['energy'] ), float( r['fwhm'] ), r['id']
  p = plans.get( pid, { 'rois': [], 'rej': [] } )
  cat = 'none'
  if any( lo <= e <= hi for lo, hi in p['rois'] ):
    cat = 'planned'
  else:
    best = None
    for lo, hi, reason in p['rej']:
      d = 0.0 if lo <= e <= hi else min( abs( e - lo ), abs( e - hi ) )
      if d <= 0.5*fwhm and ( best is None or d < best[0] ):
        best = ( d, reason )
    if best:
      cat = 'rejected: ' + best[1]
  cats[cat] += 1
  rows.append( ( pid, e, z, cat ) )

for pid, e, z, cat in sorted( rows, key = lambda t: t[3] ):
  print( '%-26s %8.1f keV z=%6.1f  %s' % ( pid, e, z, cat ) )
print()
for cat, n in cats.most_common():
  print( '%4d  %s' % ( n, cat ) )
