#!/usr/bin/env python3
"""turnon.py RUN : fitted peaks at or below the spectroscopic extent (planner's 'sub-extent below X'),
with their verdicts and the truth line they matched, so knee artefacts can be told from real lines."""
import sys, csv, os, re, collections
run = sys.argv[1]
extent = {}
prob = None
for line in open( os.path.join( run, 'roi_plan_trace.txt' ) ):
  if line.startswith( '==== ' ):
    prob = line.split()[1]; continue
  m = re.search( r'sub-extent below ([0-9.]+) keV', line )
  if m and prob and prob not in extent:
    extent[prob] = float( m.group(1) )
rows = collections.defaultdict( list )
for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
  rows[r['id']].append( r )
n = collections.Counter(); out = []
for pid, rs in rows.items():
  ext = extent.get( pid )
  if not ext: continue
  for r in rs:
    if r['set'] != 'fitted': continue
    e, f = float( r['energy'] ), float( r['fwhm'] )
    if e <= ext + 0.5*f:
      v = r['verdict']
      tz = float( r['truth_z'] or 0 ) if r.get( 'truth_z' ) else 0.0
      cls = 'real(tz>=3)' if ( v.startswith( 'matched' ) and tz >= 3 ) else ( 'weak-match' if v.startswith( 'matched' ) else v )
      n[cls] += 1
      out.append( ( pid, ext, e, f, float( r['amplitude'] ), float( r['z_det'] or 0 ), v, r.get( 'truth_energy', '' ), tz, float( r['roi_lower'] or 0 ), float( r['roi_upper'] or 0 ) ) )
print( n )
for o in sorted( out ):
  print( '%-24s ext %5.1f  peak %6.1f fwhm %5.1f area %8.0f z %6.1f  %-14s truth %s z=%.1f  roi %.1f-%.1f' % o )
