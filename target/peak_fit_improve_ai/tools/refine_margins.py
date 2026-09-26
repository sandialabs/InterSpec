#!/usr/bin/env python3
"""refine_margins.py RUN... : refinement decisions from roi_plan_trace.txt - how close the score
comparisons are, and whether rejected challengers were healthier (fewer ROIs worse than no peaks)."""
import sys, re, collections
for run in sys.argv[1:]:
  prob = None; inc = None
  rows = []
  for line in open( run + '/roi_plan_trace.txt' ):
    if line.startswith( '==== ' ):
      prob = line.split()[1]; continue
    m = re.match( r'  refinement (\d+): incumbent status \d+, (\d+) of (\d+) ROIs worse', line )
    if m:
      inc = ( int( m.group(2) ), int( m.group(3) ) ); continue
    m = re.match( r'  refinement (\d+): challenger status \d+, (\d+) of (\d+) ROIs worse.*score ([0-9.e+-]+) vs ([0-9.e+-]+) -> (\w+)', line )
    if m and inc:
      a, b = float( m.group(4) ), float( m.group(5) )
      rows.append( ( prob, int( m.group(1) ), inc, ( int( m.group(2) ), int( m.group(3) ) ), a, b, m.group(6) ) )
  acc = [ r for r in rows if r[6] == 'ACCEPTED' ]; rej = [ r for r in rows if r[6] != 'ACCEPTED' ]
  near = [ r for r in rej if r[5] <= 1.05*r[4] ]
  healthier = [ r for r in rej if r[3][0] < r[2][0] ]
  print( '%s: decisions %d accepted %d rejected %d (within 5%%: %d; challenger had fewer WORSE ROIs: %d, both: %d)' % (
    run, len( rows ), len( acc ), len( rej ), len( near ), len( healthier ), len( [ r for r in near if r[3][0] < r[2][0] ] ) ) )
  for r in rej:
    if r[3][0] < r[2][0]:
      print( '   %-24s iter %d  inc worse %d/%d  chal worse %d/%d  score %.4g vs %.4g (%.1f%%)' % ( r[0], r[1], r[2][0], r[2][1], r[3][0], r[3][1], r[4], r[5], 100*(r[5]/r[4]-1) ) )
