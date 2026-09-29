#!/usr/bin/env python3
"""Whole-spectrum review images for fit_peaks_corpus_eval runs.

  overview.py RUN_DIR ID [ID ...] [--out DIR] [--ref] [--cmp RUN_DIR2]

One PNG per problem: a log-scale view of the whole spectrum, then linear panels over energy
segments.  Every fitted ROI is a shaded band labelled with its index, extent, width in FWHM,
peak count and continuum type, so over-wide or needless ROIs are obvious at a glance.  Fitted
continuum is dashed, fitted total solid.  Truth photopeaks (z >= 3) are faint dash-dot lines;
strong truth lines the fit missed get an orange triangle.  --ref also draws the reference ROIs
(the hand fit, for a manual-fit corpus) as green bars under the red fitted bars.  --cmp draws a
second run's ROIs/total in blue, for before/after review.
"""
import os
import sys
import json
import math
import argparse

import numpy as np
import matplotlib
matplotlib.use( "Agg" )
import matplotlib.pyplot as plt


def load( run_dir, pid ):
  with open( os.path.join( run_dir, "plot_data", pid + ".json" ) ) as f:
    return json.load( f )


def roi_fwhm( r ):
  pk = r.get( "peaks", [] )
  if pk:
    return float( np.median( [p["fwhm"] for p in pk] ) )
  return None


def fwhm_at( d, e ):
  """FWHM at energy e from all fitted peaks (nearest-neighbour interpolation), else NaI-ish."""
  pts = sorted( (p["mean"], p["fwhm"]) for r in d.get( "fit", [] ) for p in r.get( "peaks", [] ) )
  if not pts:
    return 0.07 * math.sqrt( 662.0 * max( e, 10.0 ) )
  es = np.array( [a for a, _ in pts] )
  ws = np.array( [b for _, b in pts] )
  # scale a sqrt-like curve through the nearest peak
  k = int( np.argmin( np.abs( es - e ) ) )
  return ws[k] * math.sqrt( max( e, 10.0 ) / max( es[k], 10.0 ) )


def segments( d, cmp = None ):
  sig = [t["e"] for t in d.get( "truth", [] ) if t.get( "signal" ) and t["z"] >= 3 and t["e"] >= 15]
  fits = [r["upper"] for r in d.get( "fit", [] )] + ([r["upper"] for r in cmp.get( "fit", [] )] if cmp else [])
  emax = max( sig + fits + [300.0] )
  segs = [(15.0, min( 330.0, emax + 40.0 ))]
  if emax > 300.0:
    segs.append( (280.0, min( 1100.0, emax + 60.0 )) )
  if emax > 1050.0:
    segs.append( (1000.0, emax + 120.0) )
  return segs


def draw_rois( ax, d, x, lo, hi, ymax, color, band_y, band_h, label_rois, alpha_fill ):
  rois = [r for r in d.get( "fit", [] ) if r["upper"] > lo and r["lower"] < hi]
  for k, r in enumerate( sorted( rois, key = lambda r: r["lower"] ) ):
    idx = d["fit"].index( r )
    first = r["first"]
    m = len( r["cont"] )
    cx = x[first:first + m + 1]
    cont = np.asarray( r["cont"], dtype = float )
    total = cont.copy()
    for p in r["peaks"]:
      total += np.asarray( p["y"], dtype = float )
    shade = "#fde0dd" if (idx % 2 == 0) else "#e7e1ef"
    if alpha_fill > 0:
      ax.axvspan( r["lower"], r["upper"], color = shade, alpha = alpha_fill, lw = 0, zorder = 0 )
    ax.step( cx, np.append( cont, cont[-1] ), where = "post", color = color, ls = "--", lw = 0.9, zorder = 3 )
    ax.step( cx, np.append( total, total[-1] ), where = "post", color = color, lw = 1.2, zorder = 3 )
    ax.plot( [r["lower"], r["upper"]], [band_y, band_y], color = color, lw = band_h, solid_capstyle = "butt",
             alpha = 0.55, zorder = 2 )
    if label_rois:
      w = roi_fwhm( r ) or fwhm_at( d, 0.5 * (r["lower"] + r["upper"]) )
      wf = (r["upper"] - r["lower"]) / w if w > 0 else 0
      ax.text( 0.5 * (max( lo, r["lower"] ) + min( hi, r["upper"] )), band_y,
               "R%d %.0f-%.0f  %.1fF  %dpk %s" % (idx, r["lower"], r["upper"], wf, len( r["peaks"] ), r["type"]),
               color = "black", fontsize = 6.5, ha = "center", va = "bottom", zorder = 5 )
    for p in r["peaks"]:
      if not (lo <= p["mean"] <= hi):
        continue
      strong = (p.get( "z" ) or 0) >= 2
      i = int( np.searchsorted( x, p["mean"] ) ) - 1 - first
      top = total[i] if 0 <= i < len( total ) else ymax * 0.3
      ax.plot( [p["mean"], p["mean"]], [0, top], color = color if strong else "#aaaaaa", lw = 0.6 if strong else 0.4,
               alpha = 0.8, zorder = 2 )
      if strong and label_rois:
        ax.text( p["mean"], min( top, ymax * 0.92 ), "%.0f z%.0f" % (p["mean"], p.get( "z", 0 )), color = color,
                 fontsize = 5.5, ha = "center", va = "bottom", rotation = 90, zorder = 4 )


def draw_ref_bars( ax, d, lo, hi, band_y, color = "#2ca02c" ):
  for r in d.get( "ref", [] ):
    if r["upper"] > lo and r["lower"] < hi:
      ax.plot( [r["lower"], r["upper"]], [band_y, band_y], color = color, lw = 5, solid_capstyle = "butt", alpha = 0.6 )
      first = r["first"]


def draw_ref_model( ax, d, x, lo, hi ):
  for r in d.get( "ref", [] ):
    if not (r["upper"] > lo and r["lower"] < hi):
      continue
    first = r["first"]
    m = len( r["cont"] )
    cx = x[first:first + m + 1]
    cont = np.asarray( r["cont"], dtype = float )
    total = cont.copy()
    for p in r["peaks"]:
      total += np.asarray( p["y"], dtype = float )
    ax.step( cx, np.append( cont, cont[-1] ), where = "post", color = "#2ca02c", ls = ":", lw = 0.9, zorder = 3 )
    ax.step( cx, np.append( total, total[-1] ), where = "post", color = "#2ca02c", lw = 0.9, alpha = 0.8, zorder = 3 )


def panel( ax, d, lo, hi, log_y, show_ref, cmp ):
  x = np.asarray( d["x"] )
  y = np.asarray( d["y"], dtype = float )
  i0 = max( 0, int( np.searchsorted( x, lo ) ) - 1 )
  i1 = min( len( y ), int( np.searchsorted( x, hi ) ) + 1 )
  xs = x[i0:i1 + 1]
  ys = y[i0:i1]
  if len( ys ) == 0:
    return
  ax.step( xs, np.append( ys, ys[-1] ), where = "post", color = "black", lw = 0.8, zorder = 1 )
  if "bg_y" in d:
    bg = np.asarray( d["bg_y"], dtype = float )[i0:i1] * d.get( "bg_scale", 1.0 )
    ax.step( xs, np.append( bg, bg[-1] ), where = "post", color = "#6baed6", lw = 0.6, alpha = 0.8, zorder = 1 )
  ymax = float( ys.max() ) if len( ys ) else 1.0
  top = ymax * (8.0 if log_y else 1.32)
  bands = (ymax * 2.2, ymax * 3.6) if log_y else (ymax * 1.10, ymax * 1.20)
  draw_rois( ax, d, x, lo, hi, ymax, "#d62728", bands[1] if cmp else bands[0], 5, not log_y, 0.35 if not log_y else 0.25 )
  if cmp is not None:
    draw_rois( ax, cmp, x, lo, hi, ymax, "#1f5fbf", bands[0], 5, True, 0.0 )
  if show_ref:
    draw_ref_model( ax, d, x, lo, hi )
    draw_ref_bars( ax, d, lo, hi, bands[0] * (0.93 if not log_y else 0.7) )
  for t in d.get( "truth", [] ):
    if t.get( "signal" ) and t["z"] >= 3 and lo <= t["e"] <= hi:
      ax.axvline( t["e"], color = "#2ca02c", lw = 0.5, alpha = 0.45, ls = "-.", zorder = 0 )
      if t["z"] >= 8 and not log_y:
        ax.text( t["e"], ymax * 0.02, "%.0f" % t["e"], color = "#2ca02c", fontsize = 5, ha = "center", va = "bottom",
                 rotation = 90, alpha = 0.8 )
  # strong truth peaks the fit missed (the inject reference carries the per-truth verdict)
  for r in d.get( "ref", [] ):
    for p in r.get( "peaks", [] ):
      if p.get( "verdict" ) == "missed" and (p.get( "z" ) or 0) >= 8 and lo <= p["mean"] <= hi:
        ax.plot( p["mean"], ymax * 0.05 if not log_y else max( 1, ys.min() ), marker = "^", color = "#ff7f0e", ms = 7,
                 zorder = 6 )
  ax.set_xlim( lo, hi )
  if log_y:
    pos = ys[ys > 0]
    ax.set_yscale( "log" )
    ax.set_ylim( max( 0.5, pos.min() * 0.5 ) if len( pos ) else 0.5, top )
  else:
    ax.set_ylim( 0, top )
  ax.tick_params( labelsize = 7 )
  ax.grid( axis = "x", alpha = 0.25, lw = 0.4 )


def render( run_dir, pid, out_dir, show_ref, cmp_dir ):
  d = load( run_dir, pid )
  cmp = load( cmp_dir, pid ) if cmp_dir else None
  segs = segments( d, cmp )
  x = d["x"]
  full = (15.0, max( s[1] for s in segs ))
  n = 1 + len( segs )
  fig, axes = plt.subplots( n, 1, figsize = (17, 3.3 * n) )
  panel( axes[0], d, full[0], full[1], True, show_ref, cmp )
  for k, (lo, hi) in enumerate( segs ):
    panel( axes[k + 1], d, lo, hi, False, show_ref, cmp )
  nroi = len( d.get( "fit", [] ) )
  npk = sum( len( r["peaks"] ) for r in d.get( "fit", [] ) )
  title = "%s  (%s)  %s  raw %.1f   %d ROIs, %d peaks" % (pid, d.get( "sources", "" ), d.get( "status", "" ),
                                                        d.get( "raw_cost", 0 ), nroi, npk)
  if cmp is not None:
    title += "    RED = %s, BLUE = %s" % (os.path.basename( run_dir.rstrip( "/" ) ), os.path.basename( cmp_dir.rstrip( "/" ) ))
  if show_ref:
    title += "    green = reference (hand fit)"
  fig.suptitle( title, fontsize = 10 )
  fig.tight_layout( rect = (0, 0, 1, 0.975) )
  os.makedirs( out_dir, exist_ok = True )
  out = os.path.join( out_dir, pid + ".png" )
  fig.savefig( out, dpi = 85 )
  plt.close( fig )
  print( "wrote", out )


def main():
  ap = argparse.ArgumentParser( description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter )
  ap.add_argument( "run_dir" )
  ap.add_argument( "ids", nargs = "+" )
  ap.add_argument( "--out", default = None )
  ap.add_argument( "--ref", action = "store_true" )
  ap.add_argument( "--cmp", default = None )
  a = ap.parse_args()
  out = a.out or os.path.join( a.run_dir, "overview" )
  for pid in a.ids:
    try:
      render( a.run_dir, pid, out, a.ref, a.cmp )
    except FileNotFoundError as e:
      print( "missing", pid, e )


if __name__ == "__main__":
  main()
