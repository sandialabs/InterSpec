#!/usr/bin/env python3
"""search_review.py A B OUT_DIR [--verdict v1,v2] [--per-page N] : renders the search peaks on no truth
photopeak (`review_peaks.jsonl` of `fit_peaks_corpus_eval --search-only`) that run(s) B have and A do not.

A and B are run directories, or tags of tools/search_all.sh (every ${tag}_<set>_<dwell> pair).  Each panel
shows the measured counts (grey steps), the reference spectrum (the PCF's noise-free foreground, or its long
background scaled, for background/recovered peaks; blue), the search fit's model of the peak's ROI (red) and
the peak position, titled with the verdict, the reference excess z (>= 2: a real feature the truth list
does not hold) and the detection z.  Writes OUT_DIR/review_<k>.png and OUT_DIR/index.html.
Needs matplotlib (use $FPR_PY).
"""
import os, sys, glob, json, argparse
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt


def run_pairs( a, b, work ):
    if os.path.isdir( a ) and os.path.isdir( b ):
        return [(os.path.basename( b.rstrip( '/' ) ), a, b)]
    pairs = []
    for run_b in sorted( glob.glob( os.path.join( work, 'runs', b + '_*' ) ) ):
        suffix = os.path.basename( run_b )[len(b) + 1:]
        run_a = os.path.join( work, 'runs', a + '_' + suffix )
        if os.path.isdir( run_a ) and os.path.exists( os.path.join( run_b, 'review_peaks.jsonl' ) ):
            pairs.append( (suffix, run_a, run_b) )
    return pairs


def load( run ):
    path = os.path.join( run, 'review_peaks.jsonl' )
    if not os.path.exists( path ):
        return []
    with open( path ) as f:
        return [json.loads( line ) for line in f if line.strip()]


def search_energies( run ):
    """(problem, set) -> energies of every search peak in the run (search_peaks.tsv)."""
    import csv
    found = {}
    with open( os.path.join( run, 'search_peaks.tsv' ) ) as f:
        for r in csv.DictReader( f, delimiter = '\t' ):
            found.setdefault( (r['problem'], r['set']), [] ).append( float( r['energy'] ) )
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument( 'a' ); ap.add_argument( 'b' ); ap.add_argument( 'out' )
    ap.add_argument( '--verdict', default = 'unexplained,real_untabulated' )
    ap.add_argument( '--per-page', type = int, default = 24 )
    args = ap.parse_args()
    work = os.environ.get( 'FPR_WORK', os.path.expanduser( '~/fit_peaks_work' ) )
    verdicts = set( args.verdict.split( ',' ) )

    panels = []
    for name, run_a, run_b in run_pairs( args.a, args.b, work ):
        in_a = search_energies( run_a )
        for p in load( run_b ):
            if p['verdict'] not in verdicts:
                continue
            window = max( 1.0, 0.5*p['fwhm'] )
            if any( abs( e - p['energy'] ) < window for e in in_a.get( (p['problem'], p['set']), [] ) ):
                continue
            p['run'] = name
            panels.append( p )

    os.makedirs( args.out, exist_ok = True )
    ncol = 4
    pages = [panels[i:i + args.per_page] for i in range( 0, len( panels ), args.per_page )]
    html = ['<html><body><h3>%d peaks in %s and not in %s (%s)</h3>' % (len( panels ), args.b, args.a, args.verdict)]
    for k, page in enumerate( pages ):
        nrow = (len( page ) + ncol - 1) // ncol
        fig, axes = plt.subplots( nrow, ncol, figsize = (4.2*ncol, 3.0*nrow), squeeze = False )
        for ax in axes.flat:
            ax.set_visible( False )
        for ax, p in zip( axes.flat, page ):
            ax.set_visible( True )
            x, y = p['x'], p['y']
            ax.step( x, y, where = 'post', color = '0.55', lw = 0.8, label = 'data' )
            if any( v > 0 for v in p['ref'] ):
                ax.plot( x, p['ref'], color = 'tab:blue', lw = 1.0, label = 'reference' )
            model = [(xv, mv) for xv, mv in zip( x, p['model'] ) if mv > 0]
            if model:
                ax.plot( [m[0] for m in model], [m[1] for m in model], color = 'tab:red', lw = 1.0, label = 'search fit' )
            ax.axvline( p['energy'], color = 'tab:green', lw = 0.8, ls = '--' )
            ax.set_title( '%s %s %s\n%.1f keV %s ref_z=%.1f det_z=%.1f' % (p['run'], p['problem'], p['set'], p['energy'],
                          p['verdict'], p['reference_z'], p['det_z']), fontsize = 7 )
            ax.tick_params( labelsize = 6 )
        axes.flat[0].legend( fontsize = 6 )
        fig.tight_layout()
        png = 'review_%d.png' % k
        fig.savefig( os.path.join( args.out, png ), dpi = 90 )
        plt.close( fig )
        html.append( '<img src="%s" style="max-width:100%%"><br>' % png )
    html.append( '</body></html>' )
    with open( os.path.join( args.out, 'index.html' ), 'w' ) as f:
        f.write( '\n'.join( html ) )
    print( '%d panels on %d pages in %s' % (len( panels ), len( pages ), args.out) )


if __name__ == '__main__':
    main()
