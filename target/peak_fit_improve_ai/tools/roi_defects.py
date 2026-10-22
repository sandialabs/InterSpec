#!/usr/bin/env python3
"""Per-ROI defect indicators for a fit_peaks_corpus_eval run (reads RUN/plot_data/*.json).

  roi_defects.py RUN_DIR [--tsv OUT] [--top N CLASS]

For each fitted ROI, with T = continuum + peaks and sigma = sqrt(max(T,1)):
  width_f      ROI width in FWHM (FWHM from its peaks)
  valley_f     longest INTERIOR stretch (FWHM) where the ROI's significant peaks (z>=3) add
               less than 1 sigma, with significant signal on both sides  -> over-merged
  edge_lo_f/edge_hi_f  same, but from the ROI edge to the first significant signal -> overreach
  chi2         chi2 per channel of the data against T
  runs_z       Wald-Wolfowitz runs z of the residual signs (very negative = systematic misfit)
  cont_bias    mean (data - continuum)/sigma over the channels no significant peak covers:
               + = continuum too low (signal or curvature left over), - = continuum too high
  unres        peaks closer than 0.5 FWHM to a neighbour (unresolvable pairs -> comb)
  weak_frac    fraction of the ROI's fitted peak area carried by z<2 peaks
  max_z        largest peak z in the ROI (phantom ROI if < 3)
Flags (classes): WIDE (width_f>6), VALLEY (valley_f>=1.0), OVERREACH (edge>=2.0), MISFIT
(chi2>3 and runs_z<-4), CONT (|cont_bias|>2 over >=1 FWHM of uncovered channels), COMB (unres>=2),
PHANTOM (max_z<3), WEAKFILL (weak_frac>0.3).
"""
import os
import sys
import json
import math
import glob
import argparse
import collections

import numpy as np


def fwhm_of( r, d ):
  pk = r.get( "peaks", [] )
  if pk:
    return float( np.median( [p["fwhm"] for p in pk] ) )
  return 0.07 * math.sqrt( 662.0 * max( 0.5 * (r["lower"] + r["upper"]), 10.0 ) )


def runs_z( signs ):
  s = [v for v in signs if v != 0]
  n = len( s )
  if n < 10:
    return 0.0
  n1 = sum( 1 for v in s if v > 0 )
  n2 = n - n1
  if n1 == 0 or n2 == 0:
    return -10.0
  runs = 1 + sum( 1 for i in range( 1, n ) if s[i] != s[i - 1] )
  mu = 2.0 * n1 * n2 / n + 1.0
  var = 2.0 * n1 * n2 * (2.0 * n1 * n2 - n) / (n * n * (n - 1.0))
  return (runs - mu) / math.sqrt( var ) if var > 0 else 0.0


def analyse_roi( d, r ):
  x = np.asarray( d["x"], dtype = float )
  y = np.asarray( d["y"], dtype = float )
  first = r["first"]
  cont = np.asarray( r["cont"], dtype = float )
  m = len( cont )
  if m < 3:
    return None
  ys = y[first:first + m]
  xs = 0.5 * (x[first:first + m] + x[first + 1:first + m + 1])
  sig_sum = np.zeros( m )
  all_sum = np.zeros( m )
  weak_area = 0.0
  tot_area = 0.0
  for p in r["peaks"]:
    py = np.asarray( p["y"], dtype = float )
    all_sum += py
    a = float( py.sum() )
    tot_area += max( a, 0.0 )
    if (p.get( "z" ) or 0.0) >= 3.0:
      sig_sum += py
    elif (p.get( "z" ) or 0.0) < 2.0:
      weak_area += max( a, 0.0 )
  total = cont + all_sum
  sigma = np.sqrt( np.maximum( total, 1.0 ) )
  fw = fwhm_of( r, d )
  kev_per_ch = (r["upper"] - r["lower"]) / max( m, 1 )
  active = sig_sum >= sigma
  # interior valleys and edge stretches, in channels
  idx = np.where( active )[0]
  valley = 0
  edge_lo = edge_hi = m
  if len( idx ):
    edge_lo = idx[0]
    edge_hi = m - 1 - idx[-1]
    run = 0
    for i in range( idx[0], idx[-1] + 1 ):
      if not active[i]:
        run += 1
        valley = max( valley, run )
      else:
        run = 0
  resid = ys - total
  chi2 = float( np.sum( resid * resid / np.maximum( total, 1.0 ) ) / m )
  rz = runs_z( np.sign( resid ) )
  unc = ~active & (all_sum < sigma)
  cont_bias = float( np.mean( (ys[unc] - cont[unc]) / sigma[unc] ) ) if unc.sum() * kev_per_ch >= fw else 0.0
  means = sorted( (p["mean"], p["fwhm"]) for p in r["peaks"] )
  unres = sum( 1 for i in range( 1, len( means ) ) if (means[i][0] - means[i - 1][0]) < 0.5 * 0.5 * (means[i][1] + means[i - 1][1]) )
  zs = [(p.get( "z" ) or 0.0) for p in r["peaks"]]
  out = dict(
    lower = r["lower"], upper = r["upper"], type = r["type"], n_peaks = len( r["peaks"] ),
    n_sig = sum( 1 for z in zs if z >= 3 ), fwhm = fw,
    width_f = (r["upper"] - r["lower"]) / fw,
    valley_f = valley * kev_per_ch / fw,
    edge_lo_f = edge_lo * kev_per_ch / fw, edge_hi_f = edge_hi * kev_per_ch / fw,
    chi2 = chi2, runs_z = rz, cont_bias = cont_bias, unres = unres,
    weak_frac = (weak_area / tot_area) if tot_area > 0 else 0.0,
    max_z = max( zs ) if zs else 0.0 )
  flags = []
  if out["width_f"] > 6:
    flags.append( "WIDE" )
  if out["valley_f"] >= 1.0:
    flags.append( "VALLEY" )
  if max( out["edge_lo_f"], out["edge_hi_f"] ) >= 2.0 and out["max_z"] >= 3:
    flags.append( "OVERREACH" )
  if out["chi2"] > 3 and out["runs_z"] < -4:
    flags.append( "MISFIT" )
  if abs( out["cont_bias"] ) > 2:
    flags.append( "CONT" )
  if out["unres"] >= 2:
    flags.append( "COMB" )
  if out["max_z"] < 3:
    flags.append( "PHANTOM" )
  if out["weak_frac"] > 0.3:
    flags.append( "WEAKFILL" )
  out["flags"] = ",".join( flags )
  return out


def main():
  ap = argparse.ArgumentParser( description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter )
  ap.add_argument( "run_dir" )
  ap.add_argument( "--tsv", default = None )
  ap.add_argument( "--list", default = None, help = "print the ROIs carrying this flag" )
  a = ap.parse_args()
  rows = []
  for path in sorted( glob.glob( os.path.join( a.run_dir, "plot_data", "*.json" ) ) ):
    d = json.load( open( path ) )
    for k, r in enumerate( d.get( "fit", [] ) ):
      o = analyse_roi( d, r )
      if o is None:
        continue
      o["id"] = d["id"]
      o["roi"] = k
      rows.append( o )
  cols = ["id", "roi", "lower", "upper", "type", "n_peaks", "n_sig", "fwhm", "width_f", "valley_f", "edge_lo_f", "edge_hi_f",
          "chi2", "runs_z", "cont_bias", "unres", "weak_frac", "max_z", "flags"]
  out = a.tsv or os.path.join( a.run_dir, "roi_defects.tsv" )
  with open( out, "w" ) as f:
    f.write( "\t".join( cols ) + "\n" )
    for o in rows:
      f.write( "\t".join( ("%.4g" % o[c]) if isinstance( o[c], float ) else str( o[c] ) for c in cols ) + "\n" )
  nprob = len( { o["id"] for o in rows } )
  cnt = collections.Counter()
  probs = collections.defaultdict( set )
  for o in rows:
    for fl in filter( None, o["flags"].split( "," ) ):
      cnt[fl] += 1
      probs[fl].add( o["id"] )
  print( "%s: %d problems, %d ROIs, median width %.1f FWHM" % (a.run_dir, nprob, len( rows ),
         float( np.median( [o["width_f"] for o in rows] ) ) if rows else 0) )
  for fl in ["WIDE", "VALLEY", "OVERREACH", "MISFIT", "CONT", "COMB", "PHANTOM", "WEAKFILL"]:
    print( "  %-10s %4d ROIs in %3d problems" % (fl, cnt[fl], len( probs[fl] )) )
  if a.list:
    sel = [o for o in rows if a.list in o["flags"].split( "," )]
    key = {"WIDE": "width_f", "VALLEY": "valley_f", "MISFIT": "chi2", "COMB": "unres", "WEAKFILL": "weak_frac"}.get( a.list, "width_f" )
    for o in sorted( sel, key = lambda o: -o[key] ):
      print( "  %-26s R%-2d %7.1f-%-7.1f %-15s w=%5.1fF valley=%4.1fF edges=%4.1f/%4.1f chi2=%6.1f runs=%6.1f cb=%5.1f unres=%d weak=%.2f npk=%d flags=%s" % (
        o["id"], o["roi"], o["lower"], o["upper"], o["type"], o["width_f"], o["valley_f"], o["edge_lo_f"], o["edge_hi_f"],
        o["chi2"], o["runs_z"], o["cont_bias"], o["unres"], o["weak_frac"], o["n_peaks"], o["flags"] ) )


if __name__ == "__main__":
  main()
