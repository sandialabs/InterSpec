#!/usr/bin/env python3
"""Render overview.py images for every problem of a run (parallel).  render_all.py RUN_DIR [OUT_DIR] [--ref] [--cmp RUN2]"""
import os, sys, glob, multiprocessing as mp
sys.path.insert( 0, os.path.dirname( os.path.abspath( __file__ ) ) )
import overview
def job( a ):
  run, pid, out, ref, cmp = a
  try:
    overview.render( run, pid, out, ref, cmp )
  except Exception as e:
    print( 'FAIL', pid, e )
if __name__ == '__main__':
  args = [a for a in sys.argv[1:]]
  ref = '--ref' in args
  cmp = None
  if '--cmp' in args:
    cmp = args[args.index( '--cmp' ) + 1]
    del args[args.index( '--cmp' ):args.index( '--cmp' ) + 2]
  args = [a for a in args if a != '--ref']
  run = args[0]
  out = args[1] if len( args ) > 1 else os.path.join( run, 'overview' )
  ids = [os.path.basename( p )[:-5] for p in sorted( glob.glob( os.path.join( run, 'plot_data', '*.json' ) ) )]
  with mp.Pool( 8 ) as pool:
    pool.map( job, [(run, pid, out, ref, cmp) for pid in ids] )
