#!/usr/bin/env python3
"""invariants.py RUN [RUN...] : things that must never happen, whatever the scores say -
  fitted peaks whose mean lies outside their own ROI, ROIs that overlap another of the same spectrum,
  'ghost' peaks (fitted z < 1, a zero-amplitude peak the output should not carry), failures, timeouts and
  non-determinism flags.  Run it on every full run."""
import sys, csv, os, re
for run in sys.argv[1:]:
  outside, ghosts, overlaps, rois = [], [], [], {}
  for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
    if r['set'] != 'fitted':
      continue
    e, lo, hi = float( r['energy'] ), float( r['roi_lower'] ), float( r['roi_upper'] )
    rois.setdefault( r['id'], set() ).add( ( lo, hi ) )
    if e < lo - 0.01 or e > hi + 0.01:
      outside.append( '%s@%.1f' % ( r['id'], e ) )
    if r['verdict'] == 'ghost':
      ghosts.append( '%s@%.1f' % ( r['id'], e ) )
  for pid, ranges in rois.items():
    ranges = sorted( ranges )
    for ( lo0, hi0 ), ( lo1, hi1 ) in zip( ranges, ranges[1:] ):
      if lo1 < hi0 - 0.01:
        overlaps.append( '%s@%.0f-%.0f/%.0f-%.0f' % ( pid, lo0, hi0, lo1, hi1 ) )
  log = run + '.log'
  summary = timeouts = ''
  if os.path.exists( log ):
    text = open( log, errors = 'replace' ).read()
    timeouts = ' '.join( re.findall( r'\] (\S+) Timeout', text ) )
    m = re.findall( r'^problems=.*$', text, re.M )
    summary = m[-1][:90] if m else ''
  print( '%s\n  %s\n  outside-ROI peaks: %d %s\n  overlapping ROIs: %d %s\n  ghosts: %d %s\n  timeouts: %s' % ( run, summary,
         len( outside ), ' '.join( outside[:10] ), len( overlaps ), ' '.join( overlaps[:10] ), len( ghosts ),
         ' '.join( ghosts[:10] ), timeouts or 'none' ) )
