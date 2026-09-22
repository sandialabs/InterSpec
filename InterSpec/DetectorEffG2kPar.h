#ifndef DetectorEffG2kPar_h
#define DetectorEffG2kPar_h
/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <array>
#include <string>
#include <memory>
#include <vector>
#include <cstdint>
#include <istream>

/*
 =============================================================================
  Reading a detector-characterization efficiency parameter file
 =============================================================================

 A number of commercial gamma-spectroscopy suites characterize a detector by
 shipping the result of a one-time Monte-Carlo run as a pair of files:

   * a binary ".PAR" file  - a spatial full-energy-peak (FEP) efficiency grid,
     tabulated over source position (radial range from the face x polar angle)
     at a set of characterization energies; and
   * an ASCII "DETECTOR.txt" record - the detector geometry (crystal size,
     endcap, gaps) and a layer stack, linked to the ".PAR" by filename.

 Together they let the suite report absolute FEP efficiency at an arbitrary
 point in space.  This reader converts such a pair into an InterSpec
 `DetectorPeakResponse` that is backed by a `ceelo::DetectorResponse`, so the
 imported detector supports the same arbitrary-position (distance + angle)
 queries, and can be saved as a standalone `.drf.xml`.

 -----------------------------------------------------------------------------
  PROVENANCE OF THE FORMAT KNOWLEDGE
 -----------------------------------------------------------------------------
 The format of the DETECTOR.txt file is fairly obvious and needed little decoding,
 and the .PAR file is also fairly obviously formatted with standard binary 
 representations; every field meaning, byte offset, quantization constant, and 
 grid axis semantic encoded below was determined only from inspecting the 
 binary/ASCII files themselves, only starting from the general idea the .PAR
 file encoded a grid of efficiencies and using simple 1/r^2 relationship to mostly 
 figure everything out, and then doing cross-checking against a efficiencies for 
 a handful of .ecc files to validate things.  No vendor documentation, specification 
 or other proprietary material was consulted for this work.  
 However, this most likely means the decoding is not 100% correct, or may not 
 generalize beyond the inital .txt/.PAR files (there are likely other variants,
 and I likely messed something up).

 -----------------------------------------------------------------------------
  WHAT THE GRID ENCODES, AND HOW IT IS EVALUATED
 -----------------------------------------------------------------------------
   * Cell value V (uint16) -> efficiency:  eff = 10^(-V/1000).
   * Grid columns = polar angle theta from the detector axis (0..180 deg).
   * Grid rows    = ln( RADIAL distance in mm ) - the straight-line range from the
     endcap-face centre to the source, NOT the perpendicular (axial) height above
     the face plane.  The grid is SPHERICAL in (radial range, polar angle), and
     the example ".gis" inputs say so outright: every off-axis run of theirs is
     tagged `~Geometry=SPHERE`, with the run's `~sd1` the radial range and `~sd4`
     only fixing the direction (theta = atan2(sd4, sd1)).  On the detector axis
     radial and axial coincide, which is why on-axis points match either way.

     I initially thought the axis was axial, but:
       1. Physics, needing no file and no MC: hold the grid row fixed and sweep
          theta 0->85 deg.  The stored efficiency stays the same ORDER: over rows
          of 100/250/1000 mm and 45-2000 keV it moves within -59%..+33% on one
          corpus detector and -4%..+68% on the other (at 250 mm, 300 keV: -2.7%
          and +47% respectively).  Were the row an axial height, that same sweep
          would carry the source from 25 cm to 287 cm, where 1/d^2 alone forces a
          ~132x DROP - two orders of magnitude, not tens of percent.  Staying
          within a factor of ~2 at a fixed row is only physical if the row is the
          radial range.
       2. Against each detector's own independent Monte-Carlo model, over a 2-D
          sweep (5-100 cm x 0-85 deg): the radial reading agrees to 3.1% / 9.5%
          worst-case, while the axial reading is off by up to 119x / 62x, with the
          discrepancy growing as 1/cos^2(theta).
       3. The reference outputs themselves: a run at sd1=250 mm, sd4=2500 mm reads
          2.925e-03 at 100 keV, essentially the on-axis 25 cm value (3.154e-03).
          A source truly 251 cm away would read ~94x smaller.

   * The grid is the intrinsic-crystal FEP efficiency (dead layers already
     folded in); it is interpolated bilinearly in the log-efficiency (raw-V)
     domain over (ln d_radial, theta).  Across energy it is interpolated with a
     shape-preserving monotone cubic (PCHIP) over ln(E) - a plain cubic spline
     overshoots the efficiency knee between nodes, while plain linear
     interpolation undershoots it badly (the knee is log-log concave, so the
     node-to-node chord runs below the true curve - up to ~9% low at 186 keV).
   * ONLY full-energy-peak efficiency is characterized.  Total efficiency does
     not look to be derivable from these files, so the assembled response is marked
     `ceelo::TotEffTier::NotCharacterized`: its eps_total queries refuse rather
     than return a number, `DetectorPeakResponse::hasAnyTotalEfficiencyInfo()`
     is false, and cascade-summing corrections are therefore unavailable on an
     imported detector.

 -----------------------------------------------------------------------------
  ENERGY INTERPOLATION: WE USE PCHIP AND DIFFER FROM ISOCS BY ~1-2% BETWEEN NODES
 -----------------------------------------------------------------------------
 A reference efficiency-output (".ecc") file reports efficiency at a FIXED,
 GENERIC energy list (observed: 45, 60, 80, 100, 150, 200, 300, 500, 700, 1000,
 1400, 2000 keV) that is the same for every detector and unrelated to that
 detector's ".PAR" energy nodes.  Only some of those energies are grid nodes:

   * AT a node (45/60/80/100/300/500 here) we agree with ".ecc" to between ~6
     digits (the best cases) and 0.0785% (the worst, which is what `EccMatch`
     gates, at 0.1%), so both tools demonstrably hold the same stored values and
     the residue is the vendor's own 4-significant-figure ASCII rounding plus its
     geometry handling, not a difference in the numbers.
   * OFF a node (150/200/700/1000/1400/2000) ".ecc" is the reference tool's own
     interpolation of those same values - not ground truth.  So the worst gap of 
     ~1-2% (at 150/200 keV), is two interpolators drawing different curves
     through identical samples.  Neither the 122 nor the 186 keV node is ever in
     the ".ecc" list, so both tools reconstruct that region from equally sparse
     data.

 What the off-node ".ecc" values agree best with: interpolating in the log-log
 domain (raw V versus ln E, the domain we also use - it beats linear-in-energy by
 ~5x), using a 3-point Lagrange polynomial through the FORWARD triple
 {E_i, E_i+1, E_i+2} where (E_i, E_i+1) brackets the query, plus a nearly
 constant ~+1 raw-V offset.  Fitted over 16 vacuum point-like positions (on- and
 off-axis) on two detectors, that rule holds to a standard deviation of ~0.4 raw-V
 per detector (worst ~2), against ~7.4 for PCHIP, ~9.0 for a plain bracketing
 chord and ~5 for global Akima / natural / not-a-knot splines.  Values are
 deterministic, not Monte-Carlo noise: at node energies ".ecc" matches our grid to
 ~0.001 raw-V, the file's own print precision.

 There is a hard agreement floor of about +/-1 raw-V (+/-0.23% in efficiency;
 eff = 10^(-V/1000), so one raw-V step is 0.230%).  The reference appears to
 interpolate each grid CELL in energy and re-quantize to uint16 BEFORE the
 spatial interpolation, making one step of its intermediate grid 0.230%.  No
 smooth formula can reverse-engineer ".ecc" more tightly than that, so claims of
 agreement at the 0.02% level are not supportable.

 WE USE PCHIP instead - same log-log domain, but shape-preserving.  It disagrees
 with ".ecc" worst at 200 keV (~6 raw-V, ~1.4% higher efficiency) and typically
 to ~0.2% elsewhere.  An independent Monte-Carlo transport (CeeLo, 0.1%
 precision, node-normalized) is the only unbiased arbiter of which is right, and
 it puts PCHIP closest to physical truth (median ~0.13% versus ~0.5% for
 ".ecc").  Measured as percent above the straight log-log chord of the bracketing
 nodes at 200 keV: truth +1.15%, PCHIP +1.6%, the reference +0.1% - we overshoot
 slightly, the reference undershoots by about 1%.

 So the tradeoff is: adopting the forward-triple rule above would match
 the reference ~3x more tightly while being LESS physically accurate.  We choose
 physical accuracy and overshoot-safety on real spectral knees.  Only more
 characterization energies through 150-300 keV, which the files do not carry,
 would satisfy both goals.

 A note that applies to the section below as well: when the assembled response's
 low-energy accuracy was fixed, node DENSITY was the lever, not the interpolator.
 PCHIP is unchanged, and the ".ecc" between-node comparison stayed bit-identical
 (worst 0.017222) across that work - which is the correct outcome, since the
 Monte-Carlo arbiter above says closing that gap would make us less accurate.

 -----------------------------------------------------------------------------
  HOW ACCURATELY THE ASSEMBLED RESPONSE REPRODUCES THE GRID
  (all figures below are printed by `GridReproductionByBand` in
   "target/testing/test_DetectorEffG2kPar.cpp"; regenerate, do not trust prose)
 -----------------------------------------------------------------------------
 The grid ITSELF is exact at its own nodes at every energy - the round-trip test
 reproduces stored values to ~6 digits down to the lowest node.  The ASSEMBLED
 CeeLo response does not store the grid: it factorizes efficiency into an
 analytic kernel times PCHIP-interpolated corrections (see step 7 of `makeDrf`),
 so all of its error is interpolation error in that factorization.  This section
 was previously a list of double-digit percentages attributed to the files; most
 of that was ours, and two defects and three under-resolved axes have since been
 fixed.  What remains, and why, is below.

 THREE METRICS, NEVER TO BE CONFLATED.  Every number here is one of:
   1. ON-LATTICE - a file energy node, an exact grid row, an exact grid column.
      Truth is the stored uint16 with NO interpolation anywhere.  This isolates
      our representation error and is the only thing gated.
   2. OFF-NODE IN ENERGY - geometric midpoints of adjacent response nodes, the
      worst case for cubic-Hermite error.  Truth still needs `ParEfficiency`'s
      own energy PCHIP, so this is us-versus-us.
   3. OFF-NODE IN d/theta - truth includes `interp_grid_V`'s bilinear-in-raw-V
      kinks, which are an artifact of the reader, not of the detector.  Reported
      for context; never gated.
 Relative error is also meaningless where the efficiency is eight decades below
 the crystal's peak, so the sub-22 keV bands are gated on ABSOLUTE error and the
 report prints the efficiency each worst-case relative error sits on.

 DEFECT 1, a real bug: a fabricated ~95x discontinuity at 11.107 keV.
 `EtaTable` interpolates in ln E segmented at the crystal K-edges, and its
 contract is that both flanks of every edge are nodes.  `makeDrf` asked for the
 edges and stored them but never added the flanks - the one producer in the tree
 that skipped it.  On 18211381 the Ge K-edge at 11.107 keV falls between the
 file's 10 and 12 keV nodes, leaving 10 keV ALONE in its segment, which
 `EtaTable::finalize` served as a constant stub: ln eta froze at the 10 keV value
 below the crossover and at the 12 keV value above it, 4.55 apart in ln.  The
 grid's own 10 and 12 keV pages are bit-identical, so the truth across that
 interval is exactly flat in eps - the entire jump was invented.  `makeDrf` now
 builds its energy axis (flanks included) before filling either table, and
 `EtaTable::finalize` folds a lone-node run into a neighbour when that node
 cannot stand in for its segment - see the comment there for why a lone FLANK
 legitimately keeps its stub.  The test measures the step across the edge:
 response 0.999508x against the grid's 1.000000x, where pre-fix it was ~3.8e8x.
 LAB06 never had this: its grid starts at 45 keV, so no edge is retained.

 DEFECT 2, under-resolution, on THREE axes rather than the one first suspected.
 Each was found the same way - hold the other axes node-exact and the error goes
 to exactly 0.00000, so the residual always needed two axes off-node together:
   * ENERGY.  Twenty log-spaced nodes cannot resolve an absorption knee that
     falls ~5 decades across a handful of them while the kernel's analytic cutoff
     falls on a DIFFERENT curve, leaving the whole mismatch in the ratio.  Fixed
     by density, not by changing the interpolator: 20 file energies -> 363 eta
     nodes and 63 near-field energies (LAB06: 15 -> 141 and 31).  The near-field
     axis in particular must contain the FILE's own nodes, not just a resampling.
   * DISTANCE.  The grid stores 440 radial rows at 3.03% each and truth is
     bilinear in raw V, so there is a kink on EVERY row; the original 18-rung
     ladder spanned ~10 of them per rung.  Now 60 rungs (58 on LAB06), uniform
     with the [5, 30] cm window subdivided - `distance_ladder` in
     "src/DetectorEffG2kPar.cpp" documents why the rungs are ADDED rather than
     moved, and why the window ends at 30 cm.
   * ANGLE.  `makeDrf` fills all 37 of the grid's columns on [0, 90] deg, and an
     earlier note concluded from that the angular axis was "already node-exact".
     IT IS NOT, and that claim is struck.  The response's angular nodes are
     CRYSTAL-frame cosines while the grid is indexed in the ENDCAP-FACE frame,
     and the two origins differ by `endcap_front_offset_cm`.  Measured AT 45 deg
     on 18211381 (offset 0.5700 cm), the skew against 2.5 deg columns is +0.077
     deg at 300 cm, +0.280 at 83 cm and +3.04 at 8 cm - it grows as the offset
     becomes a bigger fraction of the standoff, and by 1 cm (+34 deg at 45 deg,
     and larger at steeper angles) the frame correspondence is gone rather than
     merely perturbed.  State the angle when quoting these: the skew depends on
     it, rising toward theta = 90 deg (at 83 cm it is +0.28 at 45 deg but +0.39
     approaching 90).  Refining the eta angular axis does not help - it is cheap
     (1.13x bytes) and made a fixed 22-32 keV on-axis locus WORSE (0.00893 ->
     0.01033 at 289 nodes; that is one locus, NOT the band maximum in the table
     below), because the axis is not what is misaligned.

 ACHIEVED ON-LATTICE WORST |error|, and the efficiency it sits on.  MAIN is
 d >= 5 cm and theta <= 75 deg; CORNER is everything else and is gated loosely
 (see below).  18211381 (grid starts at 10 keV), then LAB06 (starts at 45 keV):

     band        MAIN      eps there   band peak  | CORNER
     10- 16 keV  0.49031   2.18e-09    2.18e-09   | 0.94210   <- abs 1.1e-09
     16- 22 keV  0.04210   1.45e-10    5.71e-08   | 0.60479   <- abs 2.7e-10
     22- 32 keV  0.02710   7.93e-07    1.77e-04   | 0.25616
     32- 45 keV  0.01734   3.42e-04    7.08e-03   | 0.11716
     45- 60 keV  0.00469   4.52e-03    2.38e-02   | 0.06269
     60-200 keV  0.00254   2.80e-02    4.67e-02   | 0.04767
     200-1k keV  0.00325   1.27e-03    3.03e-02   | 0.02642
     1k -7k keV  0.00434   4.16e-03    1.08e-02   | 0.02723

     LAB06:      MAIN                             | CORNER
     45- 60 keV  0.00344   4.33e-02    4.50e-02   | 0.08930
     60-200 keV  0.00365   2.47e-02    5.01e-02   | 0.06240
     200-1k keV  0.00594   1.37e-02    2.67e-02   | 0.03301
     1k -7k keV  0.00955   4.20e-04    9.23e-03   | 0.07592

 Read the 10-16 keV row with its magnitude: the band's PEAK efficiency is
 2.18e-09, just over eight decades below the crystal's ~3e-1, and the worst ABSOLUTE error
 is 1.07e-09.  A 49% relative error on a quantity that small cannot affect a
 spectrum, which is why that band and 16-22 keV are gated absolutely.

 WHAT STILL LIMITS EACH REGIME.  Three different things, and they are not
 interchangeable:
   * >= 45 keV, MAIN - the ladder was the cause and refinement fixes it fast: a
     fixed 16-22 keV locus went 0.171 -> 0.00071 over 36 -> 144 rungs, i.e. 241x
     for a 4x refinement, which is ~4th order (4^4 = 256) and is what cubic
     Hermite gives on a locally smooth curve.  The |D2|*(H/h)^2/8 amplifier that
     motivated densifying in the first place is a SECOND-order bound (it would
     predict only 16x), so it is the conservative side of this and not an
     explanation of the observed rate - do not quote the two as one result.
     Every MAIN band on both detectors is now under
     1%.  It is NOT quantization-limited: our error against a radially smoothed
     grid equals our error against the raw bilinear, while the grid's own
     self-inconsistency at those loci is 0.00007-0.0068.  Do not write the
     "quantization-limited" claim for these bands.
   * Sub-45 keV, MAIN - a sharp off-axis feature that the frame skew shears, and
     no ladder fixes it.  At 16 keV and 9 cm the grid falls by a factor of 44 out
     to 57.5 deg, with a second difference along theta of about
     -32 raw-V LSB per 2.5 deg column at that row; because the skew moves where a
     crystal-frame column lands in the face frame, and moves it BY A DIFFERENT
     AMOUNT at each distance, no separable crystal-frame lattice tracks that
     feature at all distances at once.  Contrast 122 keV at 5 cm, which departs
     from its on-axis value by only ~30 LSB over the whole 0-77.5 deg sweep with
     the largest single column-to-column step under 5 LSB - i.e. smooth on the
     scale of a column.  That contrast is why >= 60 keV was ladder-fixable and
     this is not.
   * CORNER - ladder- and frame-limited, and partly the deferred past-90-deg
     work.  Gated separately so a known-deferred regime cannot mask an on-axis
     regression.  The V == 0 no-data sentinel matters here and ONLY here: it is
     the endcap volume itself, which lies entirely behind the face plane (zero
     such cells at theta <= 90 deg, onset at 92.5 deg), so it costs near-field
     reach only in the grazing columns where a crystal-frame rung maps to a
     face-frame point inside the can.  Those rungs are TRUNCATED - the column
     holds its innermost valid value, which is what Pchip's flat clamp below its
     first node does anyway - rather than filled with the sentinel, which would
     anchor the ladder on an impossible number.  See `makeDrf`'s step for the
     measured cost (48 of 666 rungs on one detector, 39 on the other).

 MONTE-CARLO SPOT-CHECKS AT THE KNEE, and what they can and cannot settle.
 The band table above measures FIDELITY TO THE FILE.  That cannot tell us whether
 the file is right, so CeeLo transport through the SAME descriptor was run at the
 knee, on-axis at 25 cm in vacuum, seeded, asking 1% relative FEP precision but
 capped at 4e6 events and 180 s CPU - so the achieved precision is whatever the
 cap allowed, which is the "MC sigma" note on the 16 keV row below.  Percentages
 are against `ParEfficiency` at the same point; the energies are FILE NODES and
 on-axis is column 0 exactly, but 25 cm is NOT a grid row (the rows are log-spaced
 at 3.034%, and ln(250 mm)/r_step = 184.72), so the truth side still goes through
 the reader's own interpolation in d.  By the taxonomy above this table is
 therefore metric 3, not metric 1 (18211381, then LAB06):

  Only the RATIOS are recorded - the absolute stored efficiencies they were formed
  from are the vendor's characterization output, and this repository does not carry
  those files or transcribe their content.  Re-measure with a throwaway against a
  local corpus if the absolute numbers are ever needed.

     E keV    | resp/grid   MC/grid
      16.0    | +0.0064%    (unmeasurable, see below)
      22.0    | -0.0356%    -3.68%
      32.0    | -0.0175%    +1.03%
      45.0    | -0.0091%    -1.30%
     122.0    | +0.0171%    -1.36%
     LAB06  45.0  | +0.0088%   -18.17%
     LAB06 122.0  | +0.0323%    -6.56%

 Two separate readings, and conflating them is the trap:
   * `resp/grid` is our error against the reader at a point that is node-exact in
     ENERGY and in ANGLE but off-node in DISTANCE by 0.72 of a row: 0.006-0.036%,
     at or below the +-1 LSB storage floor.  It is a useful corroboration that
     "we add essentially nothing" survives a query a user would actually make -
     but it is NOT the on-lattice claim, which is the gated band table above, and
     it is not a tighter figure than that table because the two measure different
     things.  What it does show is that the distance ladder is dense enough here
     that being 0.72 of a row off-node costs less than storage rounding.
   * `MC/grid` is NOT our error.  `eps_fep_at` is grounded to the grid, so the
     response tracks the file by construction; MC is ungrounded absolute
     transport, and this column compares the vendor's characterization against
     first-principles physics through the descriptor the files describe.  A few
     percent at 22-122 keV is consistent with the 200 keV datum above and shows
     no gross error at the knee.  LAB06's -18% at 45 keV is a real
     descriptor-versus-vendor discrepancy, on the detector whose files pin the
     geometry least well; it is an absolute-grounding question, NOT something
     node density affects, and it is out of scope here.
   * 16 keV is genuinely unmeasurable this way, and it is the run that hit the cap
     rather than the 1% request: at eps ~6e-9 the whole 4e6-event budget yields a
     handful of FEP counts, 57.7% statistical, so no comparison is reportable.
     The other rows did reach ~1%, which is why they are quoted and this one is
     not.  Adjudicating that band against truth needs position-biased MC, not more
     nodes.
 So the MC confirms we converged toward the file AND that the file is not grossly
 wrong at the knee; it is not tight enough to adjudicate at the 0.1% level, and
 the gated numbers remain fidelity-to-file by design.

 THE FLOOR, AND WHY IT IS NOT WHAT LIMITS US.  The file stores
 eff = 10^(-V/1000) as a uint16, so one raw-V step is 0.002305 in ln (0.2305%)
 and +-1 LSB is 0.001152 (0.1152%).  That is the tightest claim definable against
 the file's content, and every MAIN band above sits ABOVE it - 2.2x (18211381
 60-200 keV) to 8.3x (LAB06 1k-7k) over the bands >= 45 keV, and 15x and 23x for
 32-45 and 22-32 keV where the sheared angular feature dominates.  So no MAIN
 band is storage-limited; each is limited by the method named for its regime
 above.

 ONE RESIDUE THAT IS NOT FIXABLE BY NODES, and is not a bug.  Within +-2% of a
 retained K-edge the kernel is MID-CLIFF: `K` falls ~9 decades (2.9e-09 ->
 4.5e-18 over 11.11-11.20 keV) as the dead layer's tau jumps through the edge,
 and `eta` has to cancel that to reproduce a grid which is flat there.  The
 cancellation is exact only where both sides share an interpolation basis, and
 they do not - eta is a cubic in ln E, while `K`'s attenuation is exp(-tau) off
 the mu table's own log-log grid.  What is left is a ~10% bump over ~0.09 keV
 that recovers to 1.0002 as soon as `K` clears the transition.  It was measured on
 the synthetic `KEdgeSegmentsHaveBothFlanks` fixture, so take the SHAPE from there
 and the MAGNITUDE from the band table above: on a real detector that window sits
 inside 10-16 keV, whose peak efficiency is 2.18e-09, so the residue is ~2e-10
 absolute at worst - which is why that band is gated absolutely.  Adding nodes
 does not remove it (both sides refine, the bases still differ); only putting eta
 on the kernel's basis would, which means changing the interpolator for that.
 Outside those windows the
 response is monotone across 10-45 keV to 0.1%, and the edge RATIO is pinned at
 0.999508x against the grid's 1.000000x.  `KEdgeSegmentsHaveBothFlanks` measures
 both figures, and excludes the window explicitly rather than by a loose gate.

 ANOTHER RESULT WORTH RECORDING, because it is counter-intuitive and will
 otherwise be re-discovered: the achieved worst error is NOT monotone in ladder
 density.
 LAB06's 1k-7k MAIN band measures 0.03111 / 0.00955 / 0.02110 at 2 / 3 / 4
 subdivisions, at three different loci.  Rung PHASE against the file's fixed
 3.03%/row grid matters as much as rung spacing, and row alignment is
 unachievable for the same reason the angular axis is skewed: a crystal-frame
 rung sitting on a grid row does not map to a grid row in the face frame.  Gates
 therefore carry phase headroom instead of being set tight against one run.

 COSTS.  The serialized response grows from 398,932 bytes to 3,554,015 for
 18211381 (8.91x) and 1,739,328 for LAB06 (4.36x) - far more than the ~2.5x first
 predicted, because `ln_n` is the dominant block and scales with the ladder.
 `makeDrf` fill cost is ne*nc*nd, and the gated corpus test runs ~7 min against
 its 7200 s timeout.

 KNOCK-ON, deliberate: `CeeLoUtils::setLegacyEfficiencyFromResponse` samples 48
 points from `valid_e_min_keV` (10 keV) through `eps_fep_at`, so its lowest
 points stop being garbage and the exported legacy curve - hence
 `intrinsicEfficiency` below ~30 keV - changes.  That is this fix working on a
 second user-visible path.

 Two candidate causes were checked and RULED OUT; do not re-test them:
   * 18211381's two lowest nodes (10 and 12 keV) are BIT-IDENTICAL across all
     32120 grid cells - the file duplicates its first node.  Dropping the
     duplicate changes the errors by exactly nothing.
   * It was never a near-field or grazing-angle effect at the original diagnosis
     level: the blowup was just as large on-axis at 25 and 100 cm as at 2 cm.

 Callers can now use the assembled response over the files' full range, reading
 the table above for the bound that applies.  `ParEfficiency::efficiency` remains
 available and is exact at grid nodes, so it is still the right call for anything
 that must reproduce the vendor's stored values bit-for-bit.
 =============================================================================
*/

class DetectorPeakResponse;

/** Reader for a detector-characterization efficiency parameter file (binary
 spatial FEP-efficiency grid) plus its ASCII geometry record.  See the file
 header comment for the (files-only) provenance of the format knowledge.
 */
namespace DetEffG2kPar
{
  //===========================================================================
  //  DETECTOR.txt  (ASCII geometry record)
  //===========================================================================

  /** One entry of the 9-slot layer stack (a 3x3 front/side/back triplet).
   An empty slot (a bare ",,," line in the file) has `material.empty()`.
   */
  struct LayerEntry
  {
    std::string material;    ///< elemental symbol as written, lower-case ("ge","al","mg",...); "" = empty slot
    double thickness_mm = 0.0;
    double density = 0.0;    ///< g/cm3 as written; 0 if absent
  };//struct LayerEntry

  /** A parsed DETECTOR.txt record: the friendly name, the linked `.par`
   filename, the crystal/endcap dimensions, and the 9-slot layer stack.

   Field roles (all derived from the files, see header): D1 crystal diameter,
   D2 crystal length, D3 thin-entrance-window diameter (0 for standard n-type
   coax), D4 endcap outer diameter, D5 endcap length, D6 front crystal-to-endcap
   gap, D7 side crystal-to-endcap gap - all in millimetres.
   */
  struct DetectorDef
  {
    std::string comment;     ///< raw "#..." header line preceding the record (if any)
    std::string serial;      ///< serial number parsed from the comment (best-effort; may be empty)
    std::string name;        ///< friendly name = definition-line field 0
    std::string parFile;     ///< linked .par filename from the definition line

    double d1_crystal_diam_mm = 0.0;
    double d2_crystal_len_mm = 0.0;
    double d3_window_diam_mm = 0.0;
    double d4_endcap_od_mm = 0.0;
    double d5_endcap_len_mm = 0.0;
    double d6_front_gap_mm = 0.0;
    double d7_side_gap_mm = 0.0;

    int typeCode = 0;        ///< the constant that precedes the .par field (role unknown; kept verbatim)
    int kCode = 0;           ///< the constant that follows the .par field (role unknown; NOT a layer count)

    std::array<LayerEntry,9> layers;   ///< slots 0/3/6 = front/side/back Ge dead layers; rest = housing

    /** Coarse crystal family, from the definition line + layer stack. */
    enum class Kind
    {
      NCoax,     ///< n-type coaxial (D3==0, empty window slot, thick Ge dead layer)
      PType,     ///< p-type / BEGe / extended-range (D3>0, thin window)
      Falcon,    ///< Falcon-type (D3>0, thick back Al)
      Generic    ///< does not fit the standard patterns
    };//enum class Kind

    Kind kind = Kind::Generic;

    double d( int i ) const;   ///< D1..D7 by 0-based index (i in [0,6]); NaN out of range
  };//struct DetectorDef

  /** Parses a DETECTOR.txt holding one or many records into a list, in file
   order.  Records are delimited as the format dictates (a definition line
   located by its `.par` token, followed by its layer lines).  Never throws for
   an empty/garbage stream - returns an empty vector.
   */
  std::vector<DetectorDef> parseDetectorTxt( std::istream &input );

  /** Picks the record whose definition-line `.par` field matches `parFileName`
   (case-insensitive, basename only).  If none match and there is exactly one
   record, that record is returned.  Throws std::runtime_error if no record can
   be chosen.
   */
  DetectorDef selectDetectorDef( const std::vector<DetectorDef> &defs,
                                 const std::string &parFileName );


  //===========================================================================
  //  .PAR  (binary spatial-efficiency grid)
  //===========================================================================

  /** One energy's spatial grid: `nrows * ncols` uint16 cells, row-major.
   Columns index polar angle theta (step `theta_step_rad`, from the axis);
   rows index ln(radial distance in mm) - the straight-line range from the
   face centre (step `r_step`, from ln(1)=0 at row 0); see the header on why the
   row is radial, not axial.
   */
  struct ParGrid
  {
    uint16_t ncols = 0;          ///< number of polar-angle columns
    uint16_t nrows = 0;          ///< number of radial-distance (log) rows
    double theta_step_rad = 0.0; ///< polar-angle step between columns (rad)
    double r_step = 0.0;         ///< ln(radial_distance_mm) step between rows
    std::vector<uint16_t> V;     ///< nrows*ncols cells, row-major; eff = 10^(-V/1000)
  };//struct ParGrid

  /** A decoded .par file: an ascending energy grid and one `ParGrid` per
   energy.  `energies_keV[i]` corresponds to `grids[i]`.
   */
  struct ParFile
  {
    double emin_keV = 0.0;
    double emax_keV = 0.0;
    std::vector<double> energies_keV;   ///< ascending
    std::vector<ParGrid> grids;         ///< parallel to `energies_keV`
  };//struct ParFile

  /** Decodes a .par file from its raw bytes.  The header/record framing is
   derived from the file (from the per-record marker), not hard-coded, so both
   known size variants decode.  Validates the size identity and per-record grid
   dimensions.  Throws std::runtime_error on any inconsistency.
   */
  ParFile parseParFile( const std::vector<uint8_t> &bytes );

  /** Reads and decodes a .par file from disk. */
  ParFile parseParFile( const std::string &path );


  //===========================================================================
  //  Efficiency evaluator (a direct evaluation of the grid, no CeeLo)
  //===========================================================================

  /** Evaluates FEP efficiency from a decoded `.par`, matching the reference
   suite's SPATIAL interpolation: bilinear in the log-efficiency domain over
   (ln distance, theta).  Across energy we use a monotone cubic (PCHIP) over
   ln(energy), which deliberately differs from the reference between nodes (by
   ~1.4% at 200 keV) - see the file header comment.

   This is the ground-truth evaluator the assembled CeeLo response is built to
   reproduce, and what the validation test compares against `.ecc` output.
   */
  class ParEfficiency
  {
  public:
    explicit ParEfficiency( ParFile par );

    /** Vacuum (in-medium-free) FEP efficiency at a point.
     @param energy_keV        query energy; clamped to the characterized range
     @param dist_from_face_mm RADIAL range from the endcap-face centre to the
                              source, mm (the grid row; the straight-line range,
                              NOT the perpendicular/axial height - see the
                              header).  On-axis the two are equal; off-axis pass
                              the full range, not its axial component.
     @param theta_rad         polar angle from the detector axis (0 = on-axis front)
     @param no_data           optional; set true when the position falls entirely
                              inside the grid's V == 0 sentinel zone - i.e. INSIDE
                              THE ENDCAP, a physically unreachable source position
                              the file marks rather than characterizes.  The
                              returned value is then meaningless (1.0) and must
                              not be used; callers should refuse or flag rather
                              than serve it.  The zone is energy-independent and
                              exists only behind the endcap face plane
                              (theta > 90 deg, at close range).
     Out-of-grid positions saturate at the grid edge (never extrapolated).
     */
    double efficiency( double energy_keV, double dist_from_face_mm,
                       double theta_rad, bool *no_data = nullptr ) const;

    /** Air-path transmission `exp(-mu_air(E) * d)`, scaled to `pressure_atm`
     (1.0 = the reference dry-air density).  The reference suite applies this
     only for non-vacuum runs; the grid itself is vacuum efficiency.  Returns
     1.0 for `pressure_atm <= 0`.

     `mu_air` EXCLUDES Rayleigh (coherent) scatter, which is elastic and so
     cannot remove a photon from the full-energy peak.  Verified against paired
     air/vacuum reference runs over 45-2000 keV and 15-1000 mm: rms residual
     0.0007-0.008%, worst single point ~0.03%.  See the implementation comment
     for the measurement and the physical argument.

     @param dist_mm The path length through air from the source point to the
            detector.  This is the same radial range `efficiency()` takes as its
            grid row, so for a point source both take the same value.  ONLY VALID
            FOR A POINT-LIKE SOURCE: an extended source has a
            different path per source element and the reference tool integrates
            transmission over the source, so no single scalar `dist_mm` can
            reproduce it (for a 25 m in-situ disk at 1 m standoff the implied
            effective path is 14-23 m, energy-dependent, so passing the 1 m
            standoff overestimates transmission by 13-42%).  Extended-source air
            attenuation needs a per-element integral, which this reader does not
            attempt.
     */
    double airTransmission( double energy_keV, double dist_mm,
                            double pressure_atm ) const;

    const ParFile &parFile() const { return m_par; }

  private:
    ParFile m_par;
  };//class ParEfficiency


  //===========================================================================
  //  Top-level: build a DetectorPeakResponse
  //===========================================================================

  /** Builds a DRF from an already-parsed `.par` + geometry record.

   Constructs a `ceelo::GeometryDescriptor` from `def`, fills a
   `ceelo::DetectorResponse` from the grid (full angular FEP table + near-field
   model, so off-axis and close-in queries reproduce the grid), attaches a
   sampled legacy efficiency curve so the DRF is valid/serializable, and marks
   the DRF source as `DetectorPeakResponse::DrfSource::CharacterizationParFile`.

   The response is FEP-only: its total-efficiency tier is
   `ceelo::TotEffTier::NotCharacterized`, so cascade summing declines to run.
   A layer whose material token is not a recognizable element cannot be modeled;
   it is omitted from the geometry and noted in the DRF description rather than
   dropped silently (the grid still holds the true efficiency at its own nodes,
   but the kernel that shapes the between-node and near-field model degrades).

   Throws std::runtime_error on failure.
   */
  std::shared_ptr<DetectorPeakResponse> makeDrf( const ParFile &par,
                                                 const DetectorDef &def );

  /** Convenience: parse both files from disk, select the matching geometry
   record, and build the DRF.  Throws std::runtime_error on failure.
   */
  std::shared_ptr<DetectorPeakResponse> makeDrfFromFiles(
                                          const std::string &parPath,
                                          const std::string &detectorTxtPath );
}//namespace DetEffG2kPar

#endif //DetectorEffG2kPar_h
