#!/usr/bin/env python3
"""ROI geometry against the hand fits (manual-fit corpus run: ref = hand ROIs).
For each fitted ROI: how many hand ROIs it swallows (>=50 % of the hand ROI inside it); for each hand
ROI: how many fitted ROIs cut it.  Also fitted ROIs with no hand ROI under them (needless) and hand
ROIs with nothing fitted (missed)."""
import sys, os, json, glob
run = sys.argv[1]
tot = dict( fit = 0, hand = 0, merged_fit = 0, merged_hand = 0, split_hand = 0, needless = 0, missed = 0 )
lines = []
for path in sorted( glob.glob( os.path.join( run, 'plot_data', '*.json' ) ) ):
  d = json.load( open( path ) )
  fit = sorted( ( r['lower'], r['upper'] ) for r in d.get( 'fit', [] ) )
  hand = sorted( ( r['lower'], r['upper'] ) for r in d.get( 'ref', [] ) if r['upper'] > 20 )
  def ov( a, b ): return max( 0.0, min( a[1], b[1] ) - max( a[0], b[0] ) )
  msgs = []
  for f in fit:
    sw = [h for h in hand if ov( f, h ) >= 0.5 * ( h[1] - h[0] )]
    tot['fit'] += 1
    if len( sw ) >= 2:
      tot['merged_fit'] += 1; tot['merged_hand'] += len( sw )
      msgs.append( 'fit %.0f-%.0f swallows %d hand ROIs: %s' % ( f[0], f[1], len( sw ), ' '.join( '%.0f-%.0f' % h for h in sw ) ) )
    if not any( ov( f, h ) > 0 for h in hand ):
      tot['needless'] += 1; msgs.append( 'fit %.0f-%.0f has no hand ROI' % f )
  for h in hand:
    tot['hand'] += 1
    cut = [f for f in fit if ov( f, h ) >= 0.25 * ( h[1] - h[0] )]
    if len( cut ) >= 2:
      tot['split_hand'] += 1; msgs.append( 'hand %.0f-%.0f cut into %d fitted ROIs' % ( h[0], h[1], len( cut ) ) )
    if not any( ov( f, h ) > 0.25 * ( h[1] - h[0] ) for f in fit ):
      tot['missed'] += 1; msgs.append( 'hand %.0f-%.0f has no fitted ROI' % h )
  lines.append( '%-16s fit %2d  hand %2d   %s' % ( d['id'], len( fit ), len( hand ), ( '\n' + ' ' * 37 ).join( msgs ) ) )
print( '\n'.join( lines ) )
print( 'TOTAL fitted %(fit)d, hand %(hand)d; fitted ROIs swallowing >=2 hand ROIs: %(merged_fit)d (covering %(merged_hand)d hand ROIs); '
       'hand ROIs cut in two: %(split_hand)d; needless fitted ROIs: %(needless)d; hand ROIs with no fitted ROI: %(missed)d' % tot )
