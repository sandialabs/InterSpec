#!/usr/bin/env python3
"""Solve health from the ROI-plan trace 'refinement N: incumbent|challenger ...' lines.
Per problem: the first solve's ROIs-worse-than-no-peaks count, zero activities, and how many
challengers were rejected while themselves broken."""
import sys, re, collections
txt = open( sys.argv[1] + '/roi_plan_trace.txt' ).read()
blocks = re.split( r'^==== ', txt, flags = re.M )[1:]
stats = collections.Counter()
rows = []
for b in blocks:
  pid = b.split()[0]
  inc = re.findall( r'refinement (\d+): incumbent status (\d+), (\d+) of (\d+) ROIs worse than no peaks;(.*?)\| act(.*)', b )
  cha = re.findall( r'refinement (\d+): challenger status (\d+), (\d+) of (\d+) ROIs worse than no peaks;(.*?)\| act(.*?)\| score (\S+) vs (\S+) -> (\w+)', b )
  if not inc:
    stats['no refinement trace'] += 1
    continue
  it0 = inc[0]
  worse, total = int( it0[2] ), int( it0[3] )
  acts = [float( a.split( '=' )[1] ) for a in it0[5].split() if '=' in a]
  zero_act = any( a == 0.0 for a in acts )
  broken = total > 0 and worse * 2 >= total
  stats['problems'] += 1
  if zero_act: stats['first solve: a zero activity'] += 1
  if worse > 0: stats['first solve: >=1 ROI worse than no peaks'] += 1
  if broken: stats['first solve: >=half ROIs worse than no peaks'] += 1
  rej_broken = sum( 1 for c in cha if c[8] == 'rejected' and int( c[2] ) > 0 )
  acc = sum( 1 for c in cha if c[8] == 'ACCEPTED' )
  if rej_broken: stats['a challenger rejected while it had ROIs worse than no peaks'] += 1
  rows.append( ( pid, worse, total, zero_act, acc, len( cha ), it0[4].strip()[:150], it0[5].strip()[:60] ) )
for k, v in stats.items(): print( '%-62s %d' % ( k, v ) )
print()
for r in sorted( rows, key = lambda r: -( r[1] / max( r[2], 1 ) ) ):
  if r[1] > 0 or r[3]:
    print( '%-24s worse %d/%d zero_act=%d accepted %d/%d  %s | %s' % r )
