#!/usr/bin/env python3
"""Per-problem width-model grade: fitted/truth FWHM over matched peaks, plus the planner's model at
the truth energies (parsed from the trace's 'fwhm model' line is not possible in general, so the
FITTED peaks are used).  Prints problems whose median ratio is outside [0.75, 1.33]."""
import sys, csv, collections, statistics
run = sys.argv[1]
rows = list( csv.DictReader( open( run + '/per_peak.tsv' ), delimiter='\t' ) )
truth = collections.defaultdict( list )
fit = collections.defaultdict( list )
for r in rows:
  if r['set'] == 'truth' and r['verdict'].startswith( 'matched' ) and float( r['energy'] ) >= 30 and r['match_energy']:
    truth[r['id']].append( r )
  elif r['set'] == 'fitted':
    fit[r['id']].append( r )
ratios = {}
for pid, ts in truth.items():
  fr = []
  for t in ts:
    me = float( t['match_energy'] )
    best = min( fit[pid], key = lambda f: abs( float( f['energy'] ) - me ), default = None )
    if best is None or abs( float( best['energy'] ) - me ) > 0.01: continue
    fr.append( float( best['fwhm'] ) / float( t['fwhm'] ) )
  if fr: ratios[pid] = ( statistics.median( fr ), len( fr ) )
vals = [v for v, n in ratios.values()]
print( '%s: %d problems graded, median ratio %.2f; outside [0.75,1.33]: %d' % ( run, len( vals ), statistics.median( vals ),
       sum( 1 for v in vals if v < 0.75 or v > 1.33 ) ) )
for pid, ( v, n ) in sorted( ratios.items(), key = lambda kv: -abs( __import__('math').log( kv[1][0] ) ) ):
  if v < 0.75 or v > 1.33: print( '  %-24s ratio %.2f  (%d peaks)' % ( pid, v, n ) )
