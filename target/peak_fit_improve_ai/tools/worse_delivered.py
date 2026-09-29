#!/usr/bin/env python3
"""worse_delivered.py RUN : ROIs of each problem's final solve that fit WORSE than no peaks, yet the
final filter neither dropped nor rescued - they were delivered on the solve's own strongest-peak
significance (compute_roi_chi2_significance's max_peak_significance), which a peak the data does not
support can still pass.  The final solve is taken from roi_plan_trace.txt: the last accepted
refinement challenger, else a kept zero-activity re-solve, else the first solve."""
import sys, re, collections, os

lines = collections.defaultdict( list )
problem = None
for l in open( os.path.join( sys.argv[1], 'roi_plan_trace.txt' ) ):
  m = re.match( r'^==== (\S+)', l )
  if m:
    problem = m.group( 1 )
  elif problem:
    lines[problem].append( l.rstrip( '\n' ) )

total, problems = 0, 0
for pid, trace in lines.items():
  final = None
  for l in trace:
    if ('first solve:' in l) and (final is None):
      final = l
    if ('zeroed' in l) and ('kept the re-solve' in l):
      final = l
    if 'delivering the best candidate' in l:
      final = l
    m = re.search( r'refinement \d+: challenger .*\| score .* -> (\S+)', l )
    if m and m.group( 1 ) in ( 'ACCEPTED', 'BEST' ):
      final = l
  if not final:
    continue
  filtered = [ tuple( map( float, re.search( r'ROI (\d+)-(\d+) keV', l ).groups() ) )
               for l in trace if l.strip().startswith( 'final filter: ROI' ) ]
  delivered = []
  for lo, hi, delta in re.findall( r'(\d+)-(\d+):WORSE\(([-\d.e+]+)\)', final ):
    lo, hi = float( lo ), float( hi )
    if not any( (abs( lo - a ) <= 3) and (abs( hi - b ) <= 3) for a, b in filtered ):
      delivered.append( '%d-%d(%s)' % ( lo, hi, delta ) )
  if delivered:
    problems += 1
    total += len( delivered )
    print( '%-24s %s' % ( pid, ' '.join( delivered ) ) )
print( 'TOTAL %d ROIs in %d problems' % ( total, problems ) )
