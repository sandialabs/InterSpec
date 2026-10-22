#!/usr/bin/env python3
"""Render the ROIs of fit_peaks_corpus_eval problems as PNG panels for visual review.

Reads `<run_dir>/plot_data/<id>.json` (written by fit_peaks_corpus_eval unless --no-plot-data) and
draws, per energy interval, the data, the live-time scaled background, the fitted continuum and
total (red), the reference continuum and total (green), the automated-search peaks, and the
GADRAS-inject truth photopeaks (vertical lines with expected areas).

  plot_fit_peaks_rois.py RUN_DIR ID [ID ...]        every ROI of each problem, worst panels first
  plot_fit_peaks_rois.py RUN_DIR ID --range 50,70   one panel over an energy range
  plot_fit_peaks_rois.py RUN_DIR --worst 12         the 12 worst problems of the run

Output: RUN_DIR/plots/<id>.png (or --out FILE for a single problem).
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


FIT_COLORS = { "matched": "#1f77b4", "extra": "#d62728", "extra_weak": "#d62728", "extra_bkg": "#d62728",
               "ghost": "#999999", "neutral": "#9467bd", "legit": "#17becf", "xray": "#17becf" }
# Excused reference peaks (brown) are real misses the score forgave - they must not read as matched.
REF_COLORS = { "matched": "#2ca02c", "missed": "#e6550d", "dontcare": "#bbbbbb", "matched_dontcare": "#bbbbbb",
               "dontcare_merged": "#8c564b", "dontcare_bkg": "#8c564b" }


def load_problem( run_dir, pid ):
  path = os.path.join( run_dir, "plot_data", pid + ".json" )
  with open( path ) as f:
    return json.load( f )


def merged_intervals( d, pad_kev ):
  """Union of the fitted and reference ROI extents, padded and merged where they overlap."""
  spans = []
  for key in ( "fit", "ref" ):
    for roi in d.get( key, [] ):
      spans.append( [roi["lower"] - pad_kev, roi["upper"] + pad_kev] )
  spans.sort()
  merged = []
  for lo, hi in spans:
    if merged and lo <= merged[-1][1]:
      merged[-1][1] = max( merged[-1][1], hi )
    else:
      merged.append( [lo, hi] )
  return merged


def rois_in( rois, lo, hi ):
  return [r for r in rois if (r["upper"] > lo) and (r["lower"] < hi)]


def interest( d, lo, hi ):
  """How much a panel deserves a look: presence errors, family disagreement, area pulls."""
  score = 0.0
  fit = rois_in( d.get( "fit", [] ), lo, hi )
  ref = rois_in( d.get( "ref", [] ), lo, hi )
  for r in fit:
    for p in r["peaks"]:
      if p["verdict"] in ( "extra", "extra_weak", "extra_bkg" ):
        score += 3.0
      elif p["verdict"] == "legit":
        score += 1.0
  for r in ref:
    for p in r["peaks"]:
      if p["verdict"] == "missed":
        score += 3.0
      elif p["verdict"] in ( "dontcare_merged", "dontcare_bkg" ):
        score += 1.5   # excused, but a human should still see it
  fams_fit = { r["family"] for r in fit }
  fams_ref = { r["family"] for r in ref }
  if fit and ref and (fams_fit != fams_ref):
    score += 1.5
  if len( fit ) != len( ref ):
    score += 1.0
  for t in d.get( "truth_area", [] ):
    if lo <= t["e"] <= hi:
      for key in ( "fit_pull", ):
        if t.get( key ) is not None:
          score += min( 4.0, abs( t[key] ) ) / 2.0
  return score


def draw_panel( ax, d, lo, hi, log_y = False, show_background_truth = False ):
  x = np.asarray( d["x"] )
  y = np.asarray( d["y"] )
  n = len( y )
  i0 = max( 0, int( np.searchsorted( x, lo ) ) - 1 )
  i1 = min( n, int( np.searchsorted( x, hi ) ) + 1 )
  if i1 <= i0:
    return
  xs = x[i0:i1 + 1]
  ys = y[i0:i1]
  ax.step( xs, np.append( ys, ys[-1] ), where = "post", color = "black", lw = 0.9, label = "data" )
  ymax = float( ys.max() ) if len( ys ) else 1.0

  if "bg_y" in d:
    bg = np.asarray( d["bg_y"] )[i0:i1] * d.get( "bg_scale", 1.0 )
    ax.step( xs, np.append( bg, bg[-1] ), where = "post", color = "steelblue", lw = 0.7, alpha = 0.7, label = "background (scaled)" )

  def draw_set( rois, cont_style, total_color, colors, y_frac ):
    for r in rois:
      first = r["first"]
      m = len( r["cont"] )
      cx = x[first:first + m + 1]
      cont = np.asarray( r["cont"], dtype = float )
      total = cont.copy()
      for p in r["peaks"]:
        total += np.asarray( p["y"], dtype = float )
      ax.step( cx, np.append( cont, cont[-1] ), where = "post", color = total_color, ls = cont_style, lw = 1.0 )
      ax.step( cx, np.append( total, total[-1] ), where = "post", color = total_color, lw = 1.3 )
      ax.plot( [r["lower"], r["upper"]], [ymax * y_frac] * 2, color = total_color, lw = 4, solid_capstyle = "butt", alpha = 0.6 )
      ax.text( 0.5 * (r["lower"] + r["upper"]), ymax * y_frac, r["type"], color = total_color, fontsize = 7,
               ha = "center", va = "bottom" )
      for p in r["peaks"]:
        c = colors.get( p["verdict"], total_color )
        ax.axvline( p["mean"], color = c, lw = 0.8, alpha = 0.8, ymin = 0.0, ymax = 0.06 if total_color == "#2ca02c" else 0.12 )
        label = "%.1f\n%s A=%.0f" % (p["mean"], p["verdict"], p["amp"])
        if p.get( "source" ):
          label += "\n" + p["source"]
        ax.text( p["mean"], ymax * (0.02 if total_color == "#2ca02c" else 0.13), label, color = c, fontsize = 6,
                 ha = "center", va = "bottom", rotation = 90 )

  draw_set( rois_in( d.get( "ref", [] ), lo, hi ), ":", "#2ca02c", REF_COLORS, 1.06 )
  draw_set( rois_in( d.get( "fit", [] ), lo, hi ), "--", "#d62728", FIT_COLORS, 1.12 )

  for t in d.get( "truth", [] ):
    if (lo <= t["e"] <= hi) and (t["signal"] or show_background_truth):
      c = "#2ca02c" if t["signal"] else "#7f7f7f"
      ax.axvline( t["e"], color = c, lw = 0.6, alpha = 0.5, ls = "-." )
      ax.text( t["e"], ymax * 0.98, "truth %.1f\nA=%.0f z=%.1f" % (t["e"], t["area"], t["z"]), color = c, fontsize = 6,
               ha = "center", va = "top", rotation = 90 )

  for a in d.get( "auto", [] ):
    if lo <= a["mean"] <= hi:
      ax.plot( a["mean"], ymax * 1.02, marker = "v", color = "orange", ms = 4 )

  ax.set_xlim( lo, hi )
  if log_y:
    ax.set_yscale( "log" )
    ax.set_ylim( max( 0.5, ys[ys > 0].min() * 0.5 ) if (ys > 0).any() else 0.5, ymax * 1.5 )
  else:
    ax.set_ylim( 0, ymax * 1.25 )
  ax.tick_params( labelsize = 7 )

  fams_fit = ", ".join( sorted( { r["type"] for r in rois_in( d.get( "fit", [] ), lo, hi ) } ) ) or "none"
  fams_ref = ", ".join( sorted( { r["type"] for r in rois_in( d.get( "ref", [] ), lo, hi ) } ) ) or "none"
  pulls = []
  for t in d.get( "truth_area", [] ):
    if lo <= t["e"] <= hi:
      fp = t.get( "fit_pull" )
      rp = t.get( "ref_pull" )
      pulls.append( "%.0f: fit %s / ref %s" % (t["e"], ("%.1f" % fp) if fp is not None else "-", ("%.1f" % rp) if rp is not None else "-") )
  title = "%.1f-%.1f keV   fit: %s   ref: %s" % (lo, hi, fams_fit, fams_ref)
  if pulls:
    title += "\npulls " + "; ".join( pulls[:4] ) + (" ..." if len( pulls ) > 4 else "")
  ax.set_title( title, fontsize = 7 )


def render_problem( run_dir, pid, out_path, energy_range = None, max_panels = 16, pad_kev = 4.0, log_y = False,
                    columns = 2 ):
  d = load_problem( run_dir, pid )
  if energy_range:
    panels = [list( energy_range )]
  else:
    panels = merged_intervals( d, pad_kev )
    panels.sort( key = lambda lh: -interest( d, lh[0], lh[1] ) )
    panels = panels[:max_panels]
    panels.sort()
  if not panels:
    print( "no ROIs for", pid )
    return
  rows = int( math.ceil( len( panels ) / float( columns ) ) )
  fig, axes = plt.subplots( rows, columns, figsize = (8.5 * columns, 3.6 * rows), squeeze = False )
  for k, (lo, hi) in enumerate( panels ):
    draw_panel( axes[k // columns][k % columns], d, lo, hi, log_y )
  for k in range( len( panels ), rows * columns ):
    axes[k // columns][k % columns].axis( "off" )
  fig.suptitle( "%s  (%s)  %s  raw %.2f   red = fitted (dashed continuum), green = reference (dotted), "
                "dash-dot = truth photopeaks, orange = auto-search" % (pid, d.get( "sources", "" ), d.get( "status", "" ),
                d.get( "raw_cost", 0.0 )), fontsize = 9 )
  fig.tight_layout( rect = (0, 0, 1, 0.97) )
  os.makedirs( os.path.dirname( out_path ) or ".", exist_ok = True )
  fig.savefig( out_path, dpi = 110 )
  plt.close( fig )
  print( "wrote", out_path )


def worst_problem_ids( run_dir, count ):
  rows = []
  with open( os.path.join( run_dir, "per_problem.tsv" ) ) as f:
    header = f.readline().rstrip( "\n" ).split( "\t" )
    id_col = header.index( "id" )
    raw_col = header.index( "raw_cost" )
    for line in f:
      cells = line.rstrip( "\n" ).split( "\t" )
      if len( cells ) > raw_col:
        rows.append( (float( cells[raw_col] ), cells[id_col]) )
  rows.sort( reverse = True )
  return [pid for _, pid in rows[:count]]


def main():
  ap = argparse.ArgumentParser( description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter )
  ap.add_argument( "run_dir" )
  ap.add_argument( "ids", nargs = "*" )
  ap.add_argument( "--range", help = "lo,hi keV: one panel over this range" )
  ap.add_argument( "--out", help = "output PNG (single problem only)" )
  ap.add_argument( "--worst", type = int, default = 0, help = "render the N worst problems by raw cost" )
  ap.add_argument( "--max-panels", type = int, default = 16 )
  ap.add_argument( "--pad", type = float, default = 4.0, help = "keV of padding around each ROI" )
  ap.add_argument( "--columns", type = int, default = 2 )
  ap.add_argument( "--log", action = "store_true" )
  args = ap.parse_args()

  ids = list( args.ids )
  if args.worst:
    ids += worst_problem_ids( args.run_dir, args.worst )
  if not ids:
    ap.error( "give problem ids or --worst N" )
  energy_range = None
  if args.range:
    lo, hi = args.range.split( "," )
    energy_range = (float( lo ), float( hi ))
  for pid in ids:
    pid = pid.replace( "|", "_" )
    out = args.out if (args.out and len( ids ) == 1) else os.path.join( args.run_dir, "plots", pid + (".range.png" if energy_range else ".png") )
    render_problem( args.run_dir, pid, out, energy_range, args.max_panels, args.pad, args.log, args.columns )


if __name__ == "__main__":
  main()
