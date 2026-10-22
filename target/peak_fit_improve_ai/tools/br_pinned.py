#!/usr/bin/env python3
"""br_pinned.py RUN : for each branching-ratio nuisance (dBrNN) the solve pinned at its lower bound (zero
yield, or RelActCalcAuto::Options::additional_br_min_yield_fraction when set), the truth lines (z>=5)
within 0.6 FWHM of NN keV and their verdicts.  Pinned clusters on missed strong lines mean the solve is
switching lines off to buy chi2."""
import sys, csv, re, os, collections
run = sys.argv[1]
truth = collections.defaultdict(list)
for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter='\t' ):
  if r['set'] == 'truth':
    truth[r['id']].append( ( float(r['energy']), float(r['fwhm']), float(r['truth_z'] or 0), r['verdict'] ) )
tot = collections.Counter()
for r in csv.DictReader( open( os.path.join( run, 'per_problem.tsv' ) ), delimiter='\t' ):
  w = ' '.join( v for v in r.values() if v and 'pinned' in v )
  m = re.search( r'pinned at a bound[^:]*: ([^.]*)\.', w )
  if not m:
    continue
  for par in [x.strip() for x in m.group(1).split(',') if x.strip().startswith('dBr')]:
    e = float( par[3:] )
    hits = [t for t in truth[r['id']] if abs( t[0] - e ) <= 0.6*t[1] and t[2] >= 5]
    for t in hits:
      tot[t[3]] += 1
    print( '%-22s %-7s %s' % ( r['id'], par, ', '.join( '%.1f keV z=%.0f %s' % ( t[0], t[2], t[3] ) for t in hits ) or '-' ) )
print( dict( tot ) )
