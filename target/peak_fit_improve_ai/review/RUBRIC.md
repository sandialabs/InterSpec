# Visual triage of automated NaI peak fits

You are reviewing automated gamma-spectrum peak fits made by software for low-resolution NaI
detectors, the way an experienced spectroscopist reviewing the fits would.  A human expert who
hand-fits these spectra draws ROIs (regions of interest) around each VISIBLE peak feature or tight
cluster of overlapping peaks, leaves roughly 0.5-1.5 FWHM of continuum on each side, lets adjacent
ROIs butt against each other at the valley between features, never puts an ROI where the data show
no peak structure, and expects the continuum (the smooth baseline under the peaks) to follow the
data wherever no peak is present.

## How to read an image

Each PNG is one spectrum.  Title: spectrum id, (source nuclides), status, number of ROIs/peaks.
- Top panel: whole spectrum, LOG y.  Lower panels: LINEAR y over energy segments (15-330 keV,
  280-1100 keV, 1000+ keV); each panel is auto-scaled to the data in it.
- Black step line = data (counts per channel).  Light blue = background spectrum (live-time scaled).
- Each fitted ROI = a shaded band plus a thick red bar at the top, labelled
  `R<index> <lo>-<hi>  <width>F  <n>pk <continuum type>`, where `<width>F` is the ROI width in
  units of the peak FWHM.  Hand-fit NaI ROIs are typically 2-5 FWHM wide.
- Red dashed = fitted continuum; red solid = fitted total (continuum + peaks).  Thin red vertical
  lines = fitted peak centroids, labelled `<energy> z<significance>`.
- Faint green dash-dot vertical lines = true photopeak energies (from the simulation), with the energy
  printed at the bottom for strong ones.  An ORANGE TRIANGLE marks a strong true line the fit
  missed - but judge visually whether there is actually a visible peak there; many true lines are
  invisible in the data (buried), and missing an invisible line is fine.
- Energies below ~30-40 keV are usually the detector turn-on (the spectrum rising from zero); there
  is normally nothing to fit there.

## Physics you need (reviewers without it mis-called real peaks)

- NaI/CsI show iodine K x-ray ESCAPE peaks 28-33 keV below strong lines under ~250 keV (a small peak
  near 40 keV below a 69 keV x-ray peak, near 26 keV below a 55 keV one).  They are real; a fitted peak
  there is legitimate if the data show any bump, even a modest one on the rising turn-on slope.  The
  green true-line markers include them.
- True-line markers can be intensity-weighted centroids of blended lines (Ba133's 276+303 keV shows as
  one marker near 295), and some true lines are simply invisible in the data - judge what you see.
- High-energy lines (> 1.5 MeV) also have pair-production escape peaks 511 and 1022 keV below them.

## Defect classes (use these codes)

- `MERGE`    one ROI spans two or more separate features with a valley between them where the data
             come back down to (or near) the continuum; should be split.  Give the split energy(ies).
- `OVERREACH` an ROI edge extends >= ~2 FWHM past the last visible feature into plain continuum
             (often because a weak/invisible line sits there).
- `EMPTY`    an ROI (or most of it) covers a region with no visible peak structure at all.
- `PHANTOM`  a fitted peak where the data show no peak (e.g. on a Compton edge, on a smooth slope,
             or the red total rises visibly above the data).
- `COMB`     several Gaussians modelling a broad smooth hump (the low-energy scatter hump, typically
             50-300 keV, or a Compton continuum) instead of distinct peaks.
- `CONT`     the continuum does not follow the data in peak-free parts of the ROI: too high, too low,
             diving toward zero at an ROI edge, wrong curvature; often an edge peak compensates.
- `WIDTH`    a fitted peak clearly wider/narrower than the data feature, or two features that are
             visibly separate fit as one broad peak.
- `MISSED`   a clearly visible peak (a distinct bump above the local continuum) with no ROI / no
             fitted peak on it.
- `TURNON`   an ROI inside the detector turn-on region fitting the rising edge, not a peak.
- `MISFIT`   the red total is visibly off the data inside an ROI for another reason (peak centroid
             off the feature, amplitude clearly wrong, ...).
- `OTHER`    anything else wrong - describe it.

Overall verdict per spectrum:
- `GOOD`  an expert would accept it with at most cosmetic changes.
- `MINOR` usable, but some ROI boundaries/continua should be adjusted.
- `BAD`   at least one ROI is clearly wrong (nonsense geometry, phantom peaks, continuum far off).
- `CATASTROPHIC` the fit is mostly nonsense.

## Calibration examples (look at these three images first)

In the `calibration/` directory next to this file (R500 NaI fits from an early, poor version):
- `Cs137_Sh.png` -> GOOD.  One ROI 571-727 around the 662 keV peak, continuum follows the data.
  (The broad 50-250 keV scatter hump correctly has no ROI.)
- `Th228_Unsh.png` -> BAD.  R3 425-952 is `MERGE`: the 511/583 pair, the 727 feature and the
  860 feature are separated by valleys; an expert fits 436-642, 660-797, 813-922.  R1 53-150 is
  `OVERREACH` to 150 keV, driven by a weak 123 keV peak on flat data (expert stops ~110).
  R0 29-50 is a small ROI at the top of the turn-on (`TURNON`, low severity).
- `Br76_Phantom.png` -> BAD.  R0 340-797 is `MERGE` (expert: 433-618 and 624-715) and has a
  `PHANTOM` at 374 keV (it fits the Compton-edge shoulder, not a peak) and another at 778 keV where
  the continuum dives to zero near 780-797 (`CONT`).  R1 942-1456: `CONT`, continuum falls far
  below the data at the right edge and an edge peak at 1447 compensates.  R2 1603-2350: `CONT`,
  continuum starts at zero at its left edge with a compensating edge peak at 1629.  Peaks above
  1 MeV look ~1.5x too wide (`WIDTH`).

## What to produce

For EVERY spectrum in your list, view its PNG (`<id>.png` in the image directory you were given) and append lines
to your output TSV (tab-separated, no header):

    <id>  <overall>  <class or NONE>  <roi lo-hi keV, or energy>  <severity 1-3>  <one-line note>

One line per defect (so a spectrum with 3 defects gets 3 lines, each repeating the overall
verdict); a GOOD spectrum with no defects gets one line with class `NONE`.  Severity: 3 = an
expert would call it plainly wrong / nonsense, 2 = clearly should be changed, 1 = cosmetic/debatable.
For MERGE give the split energies in the note.  Be concrete and terse.  Judge what you SEE; do not
guess at the software's internals.  When unsure whether a bump is a real peak, say so in the note
and use severity 1.

## Other detector classes (read the section for your detector; NaI reviewers can skip this)

The rubric above was written for NaI.  The same defect classes and verdicts apply to the classes
below, with these differences:

- **LaBr3 (Radseeker-LaBr3)** - about 2.3x better resolution than NaI (FWHM ~20 keV at 662 keV), so
  hand-fit ROIs are narrower in keV but still 2-5 FWHM wide.  The crystal is itself radioactive:
  La138 gives a 1436 keV line that sums with the 32-37 keV Ba K x-rays into a broad, slightly
  high feature near 1440-1470 keV, and a ~790 keV line riding on a beta hump that runs to ~1000 keV;
  La/Ba K x-rays appear at 32-38 keV in every spectrum; alpha contamination can show bumps above
  1.6 MeV.  These are real but belong to the detector, not the source - a fit that ignores them is
  fine, and one that fits them as source peaks is not a phantom.  No iodine escape peaks
  (LaBr3 has no iodine); small La K escapes ~33-38 keV below strong lines under ~200 keV can exist.
- **CZT (Kromek GR1, 1 cm3)** - FWHM ~14-16 keV at 662 keV.  The detector is DEAD below ~39 keV
  (the discriminator): the spectrum is empty or ramps up there, and a line near the threshold is
  sliced in half.  Nothing below the first live channel can be fitted, so an ROI there is TURNON,
  and a truth marker there that is "missed" is not a miss.  Peaks have a low-energy tail (hole
  trapping) - a peak whose low side falls more slowly than its high side is normal, not a CONT or
  WIDTH defect.  The crystal is tiny: above ~1.5 MeV the double-escape peak (1022 keV below the
  line) can be bigger than the full-energy peak itself (Tl208 2614 keV -> 1592 keV).  Cd/Te K x-ray
  escapes appear ~23-27 keV below strong lines under ~150 keV.  No iodine escapes.
