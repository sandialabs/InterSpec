#!/usr/bin/env python3
"""run_table.py PREFIX [PREFIX...] : truth lines found (strong z>=8 / moderate z>=3 / weak) and
significant extra peaks (verdict extra or extra_bkg), per detector and in total, for the runs full.sh
names PREFIX_r500, PREFIX_ngh, PREFIX_sam, PREFIX_labr and PREFIX_czt (e.g. $FPR_WORK/runs/c77).
Missing runs print '-'.  Truth lines below 20 keV and dont-care lines are not counted."""
import sys, csv, collections, os

DETECTORS = [ 'r500', 'ngh', 'sam', 'labr', 'czt' ]

def stats( run ):
  path = os.path.join( run, 'per_peak.tsv' )
  if not os.path.exists( path ):
    return None
  s = collections.Counter()
  for r in csv.DictReader( open( path ), delimiter='\t' ):
    if r['set'] == 'truth' and r['dont_care'] != '1' and float( r['energy'] ) >= 20:
      z = float( r['z_det'] )
      cls = 'strong' if z >= 8 else ('moderate' if z >= 3 else 'weak')
      s[cls] += r['verdict'].startswith( 'matched' )
    elif r['set'] == 'fitted' and r['verdict'] in ( 'extra', 'extra_bkg' ):
      s['extras'] += 1
  return s

def cell( s ):
  return '%d/%d/%d x%d' % ( s['strong'], s['moderate'], s['weak'], s['extras'] )

prefixes = sys.argv[1:]
print( '%-6s' % 'det' + ''.join( '%24s' % os.path.basename( p ) for p in prefixes ) )
totals = { p: collections.Counter() for p in prefixes }
for det in DETECTORS:
  row = '%-6s' % det
  for p in prefixes:
    s = stats( p + '_' + det )
    if s is None:
      row += '%24s' % '-'
      continue
    totals[p] += s
    row += '%24s' % cell( s )
  print( row )
print( '%-6s' % 'total' + ''.join( '%24s' % cell( totals[p] ) for p in prefixes ) )
print( '(strong/moderate/weak truth lines found, x = significant extra peaks)' )
