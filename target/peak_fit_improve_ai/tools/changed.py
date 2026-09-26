#!/usr/bin/env python3
"""changed.py OLD_RUN NEW_RUN [--render OUT_DIR] : list the problems whose fitted ROIs or peaks changed,
with a one-line summary each (ROI extents, peak counts, strong truth found), most-changed first; with
--render, draw before/after images (NEW in red, OLD in blue) for exactly those problems."""
import sys, os, json, glob, csv, argparse, collections


def load( run ):
  out = {}
  for p in glob.glob( os.path.join( run, 'plot_data', '*.json' ) ):
    d = json.load( open( p ) )
    rois = sorted( ( round( r['lower'] ), round( r['upper'] ), len( r['peaks'] ) ) for r in d.get( 'fit', [] ) )
    peaks = sorted( round( pk['mean'] ) for r in d.get( 'fit', [] ) for pk in r['peaks'] )
    out[d['id']] = ( rois, peaks )
  return out


def strong_found( run ):
  found = collections.Counter(); total = collections.Counter()
  for r in csv.DictReader( open( os.path.join( run, 'per_peak.tsv' ) ), delimiter = '\t' ):
    if r['set'] != 'truth' or float( r['energy'] ) < 20 or r['verdict'].startswith( 'dontcare' ):
      continue
    if float( r['z_det'] ) >= 8:
      total[r['id']] += 1
      found[r['id']] += r['verdict'].startswith( 'matched' )
  return found, total


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument( 'old' ); ap.add_argument( 'new' )
  ap.add_argument( '--render', default = None )
  a = ap.parse_args()
  A, B = load( a.old ), load( a.new )
  fa, ta = strong_found( a.old )
  fb, tb = strong_found( a.new )
  rows = []
  for pid in sorted( set( A ) & set( B ) ):
    if A[pid] == B[pid]:
      continue
    ra, pa = A[pid]; rb, pb = B[pid]
    change = abs( len( ra ) - len( rb ) ) + len( set( pa ) ^ set( pb ) )
    rows.append( ( change, pid, ra, rb, fa[pid], fb[pid], ta[pid] ) )
  rows.sort( key = lambda r: -r[0] )
  print( '%d of %d problems changed' % ( len( rows ), len( set( A ) & set( B ) ) ) )
  for change, pid, ra, rb, sa, sb, t in rows:
    fmt = lambda rr: ' '.join( '%d-%d(%d)' % x for x in rr )
    print( '%-24s strong %d->%d of %d\n     old: %s\n     new: %s' % ( pid, sa, sb, t, fmt( ra ), fmt( rb ) ) )
  if a.render:
    sys.path.insert( 0, os.path.dirname( os.path.abspath( __file__ ) ) )
    import overview
    for _, pid, *_rest in rows:
      overview.render( a.new, pid, a.render, False, a.old )


if __name__ == '__main__':
  main()
