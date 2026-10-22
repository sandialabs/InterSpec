# Review of automated NaI peak fits: absolute verdict + before/after

Read `RUBRIC.md` in this directory first.  It explains how to read the images, what a good fit
looks like, the defect classes (MERGE, OVERREACH, EMPTY, PHANTOM, COMB, CONT, WIDTH, MISSED,
TURNON, MISFIT, OTHER), the overall verdicts (GOOD, MINOR, BAD, CATASTROPHIC) and severities.
Look at its three calibration images too.

The images for this review are different in one way: each shows TWO fits of the same spectrum.
- RED = the NEW fit (thick red bars at the top labelled with the new ROIs, red dashed continuum, red
  solid total, red peak labels).  The ROI labels and widths printed are the RED fit's.
- BLUE = the OLD fit (blue bars under the red ones, blue dashed continuum, blue solid total, blue
  peak labels).
The title says `RED = <new run>, BLUE = <old run>`.

For EVERY spectrum in your list, view `<id>.png` in the image directory you were given and do two
things:

1. Judge the RED (new) fit on its own, exactly as RUBRIC.md describes: an overall verdict and one
   line per defect, ignoring the blue fit.
2. Compare RED with BLUE, the way an expert reviewing both would:
   `BETTER` (clearly better overall), `WORSE` (clearly worse overall), `MIXED` (clearly better in
   some places, clearly worse in others) or `SAME` (no meaningful difference).

Append to your output TSV (tab-separated, no header) one line per RED defect:

    <id>  <overall>  <class or NONE>  <roi lo-hi keV, or energy>  <severity 1-3>  <one-line note>

(a GOOD spectrum with no defects gets one line with class `NONE`), followed by exactly one
comparison line:

    <id>  CMP  <BETTER|WORSE|MIXED|SAME>  <severity 1-3 of the change>  <one-line note: what changed>

For WORSE and MIXED be specific (energy of the lost peak, the new phantom, the broken continuum).
Judge what you SEE; do not guess at the software's internals.
