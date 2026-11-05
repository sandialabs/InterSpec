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
          within ~7% of the on-axis 25 cm value at 100 keV.  A source truly
          251 cm away would read ~94x smaller.

   * The grid is the ABSOLUTE point-source FEP efficiency of the whole detector
     as characterized (counts per emitted photon, source in vacuum): the Ge dead
     layers AND the housing - every layer of the DETECTOR.txt stack - are
     already in it.  Callers must not apply the layer stack on top; for a point
     source, air (`ParEfficiency::airTransmission`) is the only extra factor.
     The reference outputs say so: every vacuum point-source ".ecc" (on axis at
     1.5-100 cm, SPHERE runs out to 84 deg, both corpus detectors) equals the
     raw grid at the grid's energy nodes - to <=0.003%, and 0.08% on one
     off-axis run - where applying just one detector's 2.0 mm Al front would
     lower it by 22% at 45 keV and 5% at 500 keV; and Monte Carlo through the
     full descriptor, housing included, reproduces that detector's grid at 25 cm
     to within ~2% over 22-662 keV.
   * The grid is interpolated bilinearly in the log-efficiency (raw-V)
     domain over (ln d_radial, theta).  Across energy it is interpolated with a
     shape-preserving monotone cubic (PCHIP) over ln(E) - a plain cubic spline
     overshoots the efficiency knee between nodes, while plain linear
     interpolation undershoots it badly (the knee is log-log concave, so the
     node-to-node chord runs below the true curve - up to ~9% low at 186 keV).
   * ONLY full-energy-peak efficiency is characterized.  Total efficiency does
     not look to be derivable from these files, so the assembled response is marked
     `ceelo::TotEffTier::NotCharacterized`: its eps_total queries refuse rather
     than return a number, `DetectorPeakResponse::hasAnyTotalEfficiencyInfo()`
     is false, and cascade-summing corrections are unavailable - until the user
     generates a Monte-Carlo total, which `attachTotalEfficiency` grafts onto the
     grid without touching its FEP.

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
 with ".ecc" worst at 200 keV (~7 raw-V, 1.72% higher efficiency) and typically
 by ~0.3% (median) between nodes.  An independent Monte-Carlo transport (CeeLo, 0.1%
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
 that skipped it.  On DET1381 the Ge K-edge at 11.107 keV falls between the
 file's 10 and 12 keV nodes, leaving 10 keV ALONE in its segment, which
 `EtaTable::finalize` served as a constant stub: ln eta froze at the 10 keV value
 below the crossover and at the 12 keV value above it, 4.55 apart in ln.  The
 grid's own 10 and 12 keV pages are bit-identical, so the truth across that
 interval is exactly flat in eps - the entire jump was invented.  `makeDrf` now
 builds its energy axis (flanks included) before filling either table, and
 `EtaTable::finalize` folds a lone-node run into a neighbour when that node
 cannot stand in for its segment - see the comment there for why a lone FLANK
 legitimately keeps its stub.  The test measures the step across the edge:
 response 0.999499x against the grid's 1.000000x, where pre-fix it was ~3.8e8x.
 DET06 never had this: its grid starts at 45 keV, so no edge is retained.

 DEFECT 2, under-resolution, on THREE axes rather than the one first suspected.
 Each was found the same way - hold the other axes node-exact and the error goes
 to exactly 0.00000, so the residual always needed two axes off-node together:
   * ENERGY.  Twenty log-spaced nodes cannot resolve an absorption knee that
     falls ~5 decades across a handful of them while the kernel's analytic cutoff
     falls on a DIFFERENT curve, leaving the whole mismatch in the ratio.  Fixed
     by density, not by changing the interpolator: 20 file energies -> 363 eta
     nodes and 63 near-field energies (DET06: 15 -> 141 and 31).  The near-field
     axis in particular must contain the FILE's own nodes, not just a resampling.
   * DISTANCE.  The grid stores 440 radial rows at 3.03% each and truth is
     bilinear in raw V, so there is a kink on EVERY row; the original 18-rung
     ladder spanned ~10 of them per rung.  Now 68 rungs (70 on DET06), uniform
     with the [2, 30] cm window subdivided - `distance_ladder` in
     "src/DetectorEffG2kPar.cpp" documents why the rungs are ADDED rather than
     moved, why the window ends at 30 cm, and why it starts at 2 cm.
   * ANGLE.  `makeDrf` fills all 37 of the grid's columns on [0, 90] deg, and an
     earlier note concluded from that the angular axis was "already node-exact".
     IT IS NOT, and that claim is struck.  The response's angular nodes are
     CRYSTAL-frame cosines while the grid is indexed in the ENDCAP-FACE frame,
     and the two origins differ by `endcap_front_offset_cm`.  The skew is pure
     geometry, atan2(d sin(theta), d cos(theta) - offset) - theta: AT 45 deg on
     DET1381 (offset 0.37 cm), against 2.5 deg columns it is +0.050 deg at
     300 cm, +0.181 at 83 cm and +1.94 at 8 cm - it grows as the offset becomes a
     bigger fraction of the standoff, and by 1 cm (+19.5 deg at 45 deg, and
     larger at steeper angles) the frame correspondence is gone rather than
     merely perturbed.  State the angle when quoting these: the skew depends on
     it, rising toward theta = 90 deg (at 83 cm it is +0.18 at 45 deg but +0.26
     approaching 90).  Refining the eta angular axis does not help - it is cheap
     (1.13x bytes) and made a fixed 22-32 keV on-axis locus WORSE (0.00893 ->
     0.01033 at 289 nodes, measured when the offset was still 0.57 cm; that is
     one locus, NOT the band maximum in the table below), because the axis is
     not what is misaligned.

 ACHIEVED ON-LATTICE WORST relative |error|.  MAIN is d >= 5 cm and
 theta <= 75 deg; CORNER is everything else and is gated loosely (see below).
 DET1381 (grid starts at 10 keV), then DET06 (starts at 45 keV):

     band        MAIN      | CORNER
     10- 16 keV  0.29317   | 0.83400   <- gated on absolute error
     16- 22 keV  0.05809   | 0.11464   <- gated on absolute error
     22- 32 keV  0.03470   | 0.16366
     32- 45 keV  0.01185   | 0.05617
     45- 60 keV  0.00534   | 0.02489
     60-200 keV  0.00282   | 0.01279
     200-1k keV  0.00366   | 0.00599
     1k -7k keV  0.00636   | 0.00512

     DET06:      MAIN      | CORNER
     45- 60 keV  0.00329   | 0.03385
     60-200 keV  0.00354   | 0.02593
     200-1k keV  0.00544   | 0.01999
     1k -7k keV  0.00767   | 0.04750

 Read the 10-16 keV row with its magnitude: the band's PEAK efficiency is just
 over eight decades below the crystal's peak, and the worst ABSOLUTE error is
 about a quarter of that band peak.  A 29% relative error on a quantity that small
 cannot affect a spectrum, which is why that band and 16-22 keV are gated
 absolutely.

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
   * CORNER - KERNEL-limited near contact at grazing incidence, and partly the
     deferred past-90-deg work.  Gated separately so a known-deferred regime
     cannot mask an on-axis regression.  Two things meet there:
       - The V == 0 no-data sentinel, which matters here and ONLY here: it is the
         endcap volume, entirely behind the face plane (zero such cells at
         theta <= 90 deg, onset at 92.5 deg), and grazing crystal-frame rungs map
         up to the endcap-front offset behind that plane.  Where a rung's stencil
         touches the zone, neither the grid (a sentinel, or a blend renormalized
         around one) nor the kernel (a source inside the modelled front window)
         describes it, yet it anchors real queries on the face plane, so `makeDrf`
         reads it at the face-plane point in front of it (119 of 2516 rungs on
         DET1381, 151 of 2590 on DET06).  Filling partly-touching rungs as data,
         as was done before, put DET1381's worst near-contact error at -21%
         (32 keV, on the face plane 3.5 cm out); holding them, as all-sentinel
         rungs were held, made the band worsts 2-8x worse.
       - The KERNEL's can corner.  The modelled can is the crystal plus its side
         shells, narrower than the real one (D4 is not modelled), so a source a
         few cm out near the face plane sees the crystal past the model's corner
         rather than through the real can's front.  At low energy that is a deep
         trough in K that the grid does not share - at 10 keV K falls three
         decades between adjacent rungs while the grid's eps barely moves - so
         ln(grid/K) carries a ridge ~6 in ln high and about one node wide
         (DET1381, 10 keV: from 4.9 cm, 75 deg to 6.2 cm, 70 deg, crystal frame).
         On the face plane just beside the modelled corner it is a dip at every
         energy, which the ladder resolves once its subdivided window opens at
         2 cm (see `distance_ladder`).  The ridge's crest is not resolvable by any
         affordable lattice: EccMatch's grazing 10 keV probe (50 mm, 78.7 deg)
         sits on it at -23.7%; half-step columns from 65 deg plus a window from the
         ladder floor (+74% bytes) left -19.5%, and 3.1x the bytes -11.6%.  Where
         the ridge falls moves with the geometry: sizing the can from D4 moves it
         into MAIN 16-22 keV (0.058 -> 0.098 worst over every grid row to 7 cm,
         against a 0.065 gate), so that is not done; D4 and D7's side role remain
         open questions.
     Moving the crystal to D7 (2 mm forward on DET1381) changed which rungs touch
     the zone and where the ridge falls against the lattice, and grew the CORNER
     worst figures 1.2-2.5x over every grid row from 0.8 to 7 cm, almost all of
     it on the face plane 3.6 cm out.  The two fixes above more than recover that
     (the table above); EccMatch's grazing low probe is bounded as CORNER is.

 MONTE-CARLO SPOT-CHECKS AT THE KNEE, and what they can and cannot settle.
 The band table above measures FIDELITY TO THE FILE.  That cannot tell us whether
 the file is right, so CeeLo transport through the SAME descriptor was run at the
 knee, on-axis at 25 cm in vacuum, seeded, asking 1% relative FEP precision but
 capped at 2e8 events - so the achieved precision is whatever the cap allowed,
 which is the "MC sigma" note on the 16 keV row below.  Percentages
 are against `ParEfficiency` at the same point; the energies are FILE NODES and
 on-axis is column 0 exactly, but 25 cm is NOT a grid row (the rows are log-spaced
 at 3.034%, and ln(250 mm)/r_step = 184.72), so the truth side still goes through
 the reader's own interpolation in d.  By the taxonomy above this table is
 therefore metric 3, not metric 1 (DET1381, then DET06):

  Only the RATIOS are recorded - the absolute stored efficiencies they were formed
  from are the vendor's characterization output, and this repository does not carry
  those files or transcribe their content.  Re-measure with a throwaway against a
  local corpus if the absolute numbers are ever needed.

     E keV    | resp/grid   MC/grid  (MC sigma)
      16.0    | +0.0058%    (unmeasurable, see below)
      22.0    | -0.0280%    -0.81%   (0.99%)
      32.0    | -0.0025%    +2.30%   (0.76%)
      45.0    | -0.0166%    +1.03%   (0.47%)
     122.0    | +0.0299%    +0.16%   (0.30%)
     DET06  45.0  | +0.0057%   -20.30%   (0.41%)
     DET06 122.0  | +0.0454%    -9.09%   (0.40%)

 Two separate readings, and conflating them is the trap:
   * `resp/grid` is our error against the reader at a point that is node-exact in
     ENERGY and in ANGLE but off-node in DISTANCE by 0.72 of a row: 0.003-0.046%,
     at or below the +-1 LSB storage floor.  It is a useful corroboration that
     "we add essentially nothing" survives a query a user would actually make -
     but it is NOT the on-lattice claim, which is the gated band table above, and
     it is not a tighter figure than that table because the two measure different
     things.  What it does show is that the distance ladder is dense enough here
     that being 0.72 of a row off-node costs less than storage rounding.
   * `MC/grid` is NOT our error.  `eps_fep_at` is grounded to the grid, so the
     response tracks the file by construction; MC is ungrounded absolute
     transport, and this column compares the vendor's characterization against
     first-principles physics through the descriptor the files describe - housing
     included, and with the crystal at D7 (see #geometryFromDetectorDef).  Within
     ~2% at 22-122 keV on DET1381 shows no gross error at the knee, and is itself
     evidence that the grid contains the housing: were it absent from the grid,
     the descriptor's 2.0 mm Al front would put MC at about a quarter of the grid
     at 22 keV, and 22% below it at 45 keV.  DET06's -20% at
     45 keV and -9% at 122 keV are a real descriptor-versus-vendor discrepancy on
     a thin-window record.  The likely cause, unconfirmed: the vendor does not put
     the slot-2 Al on the path through the window - with it removed, Monte Carlo
     of DET06 comes to 0.98 / 0.99 / 1.01 of the grid at 45 / 60 / 122 keV, but on
     other thin-window grids removing it closes only about half the 60 keV
     deficit.  It is an absolute-grounding question, NOT something node density
     affects, and it does not move the grid-backed efficiency (only the kernel
     and any Monte-Carlo total).
   * 16 keV is genuinely unmeasurable this way, and it is the run that hit the cap
     rather than the 1% request: at an efficiency ~8 decades below the crystal's
     peak the whole 2e8-event budget yields only ~200 FEP counts, 7.0%
     statistical (-6.7%), so no comparison is reportable at the level of the
     others.
     The other rows did reach ~1%, which is why they are quoted and this one is
     not.  Adjudicating that band against truth needs position-biased MC, not more
     nodes.
 So the MC confirms we converged toward the file AND that the file is not grossly
 wrong at the knee; it is not tight enough to adjudicate at the 0.1% level, and
 the gated numbers remain fidelity-to-file by design.

 THE FLOOR, AND WHY IT IS NOT WHAT LIMITS US.  The file stores
 eff = 10^(-V/1000) as a uint16, so one raw-V step is 0.002305 in ln (0.2305%)
 and +-1 LSB is 0.001152 (0.1152%).  That is the tightest claim definable against
 the file's content, and every MAIN band above sits ABOVE it - 2.4x (DET1381
 60-200 keV) to 6.7x (DET06 1k-7k) over the bands >= 45 keV, and 10x and 30x for
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
 inside 10-16 keV, whose peak efficiency is ~8 decades below the crystal's, so the
 residue is negligible in absolute terms - which is why that band is gated
 absolutely.  Adding nodes
 does not remove it (both sides refine, the bases still differ); only putting eta
 on the kernel's basis would, which means changing the interpolator for that.
 Outside those windows the
 response is monotone across 10-45 keV to 0.1%, and the edge RATIO is pinned at
 0.999499x against the grid's 1.000000x.  `KEdgeSegmentsHaveBothFlanks` measures
 both figures, and excludes the window explicitly rather than by a loose gate.

 ANOTHER RESULT WORTH RECORDING, because it is counter-intuitive and will
 otherwise be re-discovered: the achieved worst error is NOT monotone in ladder
 density.
 DET06's 1k-7k MAIN band measures 0.03111 / 0.00955 / 0.02110 at 2 / 3 / 4
 subdivisions, at three different loci.  Rung PHASE against the file's fixed
 3.03%/row grid matters as much as rung spacing, and row alignment is
 unachievable for the same reason the angular axis is skewed: a crystal-frame
 rung sitting on a grid row does not map to a grid row in the face frame.  Gates
 therefore carry phase headroom instead of being set tight against one run.

 COSTS.  The serialized response grows from 398,932 bytes to 4,037,378 for
 DET1381 (10.1x) and 2,028,694 for DET06 (5.09x) - far more than the ~2.5x first
 predicted, because `ln_n` is the dominant block and scales with the ladder (the
 window's opening at 2 cm is +15% of that, on both).
 `makeDrf` fill cost is ne*nc*nd.  The gated corpus suite takes ~5.6 min
 (Release, 8 cores, developer checks on) against its 7200 s timeout, ~72% of it
 in `NodeRoundTripIdentity`, which evaluates every node of both responses.

 KNOCK-ON, deliberate: `CeeLoUtils::setLegacyEfficiencyFromResponse` samples 48
 points from `valid_e_min_keV` (10 keV) through `eps_fep_at`, so its lowest
 points stop being garbage and the exported legacy curve - hence
 `intrinsicEfficiency` below ~30 keV - changes.  That is this fix working on a
 second user-visible path.

 Two candidate causes were checked and RULED OUT; do not re-test them:
   * DET1381's two lowest nodes (10 and 12 keV) are BIT-IDENTICAL across all
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
namespace ceelo
{
  class DetectorResponse;
  struct GroundingPoint;
  struct GeometryDescriptor;
}

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
   Every layer - dead layers and housing alike - is already contained in the
   `.par` grid; the stack only describes the geometry (see the header).
   */
  struct LayerEntry
  {
    std::string material;    ///< elemental symbol as written, lower-case ("ge","al","mg",...); "" = empty slot
    double thickness_mm = 0.0;
    double density = 0.0;    ///< g/cm3 as written; 0 if absent
  };//struct LayerEntry

  /** A parsed DETECTOR.txt record: the friendly name, the linked `.par`
   filename, the crystal/endcap dimensions, and the 9-slot layer stack.

   Field roles, all in millimetres and derived from the files, the reference
   outputs and Monte Carlo (see the header): D1 crystal diameter, D2 crystal
   length, D3 thin-entrance-window diameter (0 for a standard coax with a thick
   Li-diffused front contact - p-type), D4 endcap outer diameter and D5 endcap
   length (the grid's no-data zone, i.e. the can, reaches D4/2 laterally),
   D7 the depth of the crystal's front face behind the endcap's OUTER face,
   front window and endcap included, and D6 - equal to D7 on every D3 = 0 record
   seen, smaller on thin-window records, where its role is unconfirmed.
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
    double d6_mm = 0.0;      ///< D6: = D7 on standard coax; role on thin-window records unconfirmed
    double d7_mm = 0.0;      ///< D7: endcap outer face to crystal front face (front window/endcap included)

    int typeCode = 0;        ///< the constant that precedes the .par field (role unknown; kept verbatim)
    int kCode = 0;           ///< the constant that follows the .par field (role unknown; NOT a layer count)

    std::array<LayerEntry,9> layers;   ///< slots 0/3/6 = front/side/back Ge dead layers; rest = housing

    /** Coarse family, from the definition line, window slot and crystal shape; it decides whether
     the importer guesses a bore.  The names describe the record, not the crystal's doping: a thin
     window is found on p-type (BEGe, extended-range) and n-type (reverse-electrode) detectors alike.
     */
    enum class Kind
    {
      ThickWindowCoax,   ///< D3==0, empty window slot, thick Ge front contact: a standard (p-type) coax
      ThinWindowCoax,    ///< D3>0, crystal length >= 0.7 x diameter: extended-range, reverse-electrode, ...
      ThinWindowPlanar,  ///< D3>0, shorter crystal: BEGe, planar, Falcon - no bore
      Generic            ///< does not fit the standard patterns
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

  /** Cheap test of whether text (e.g., the first ~1 kB of a file) looks like a DETECTOR.txt: it
   must hold a definition line naming a `.par` file, with positive crystal diameter and length.
   */
  bool isCandidateDetectorTxt( const std::string &headerText );


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

  /** Cheap test of whether the start of a file looks like a .par grid, without reading the whole
   file: a plausible energy range, energy count, and ascending energy list in the header.
   `parseParFile` remains the authority; this only gates trying it.
   */
  bool isCandidateParFile( const uint8_t *header, size_t headerLen, size_t fileSize );


  //===========================================================================
  //  Efficiency evaluator (a direct evaluation of the grid, no CeeLo)
  //===========================================================================

  /** Evaluates FEP efficiency from a decoded `.par`, matching the reference
   suite's SPATIAL interpolation: bilinear in the log-efficiency domain over
   (ln distance, theta).  Across energy we use a monotone cubic (PCHIP) over
   ln(energy), which deliberately differs from the reference between nodes (by
   up to ~1.7% at 200 keV) - see the file header comment.

   This is the ground-truth evaluator the assembled CeeLo response is built to
   reproduce, and what the validation test compares against `.ecc` output.
   */
  class ParEfficiency
  {
  public:
    explicit ParEfficiency( ParFile par );

    /** Absolute FEP efficiency at a point, source in vacuum - housing and dead layers included, so
     nothing but air (#airTransmission) belongs on top.
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

    /** Whether a position is behind the endcap face plane with any of its interpolation stencil in
     the V == 0 no-data zone - inside the endcap, or against it.  The file does not characterize
     such a position: #efficiency there is either flagged `no_data` or a blend renormalized around
     the sentinel cells.  #makeDrf reads near-field nodes there at the endcap face in front of them.
     */
    bool inNoDataZone( double dist_from_face_mm, double theta_rad ) const;

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
   `ceelo::TotEffTier::NotCharacterized`, so cascade summing declines to run until a Monte-Carlo
   total is attached (#attachTotalEfficiency).  The DRF's setback is the descriptor's endcap-front
   offset (the GADRAS convention) - the record's D7, see #geometryFromDetectorDef - so flat-disk
   paths and displays put the crystal face where the grid's near field says it is.
   A layer whose material token is not a recognizable element cannot be modeled;
   it is omitted from the geometry and noted in the DRF description rather than
   dropped silently (the grid still holds the true efficiency at its own nodes,
   but the kernel that shapes the between-node and near-field model degrades).

   Throws std::runtime_error on failure.
   */
  std::shared_ptr<DetectorPeakResponse> makeDrf( const ParFile &par,
                                                 const DetectorDef &def );

  /** Whether `resp` is a response #makeDrf built from a grid - an "imported efficiency grid".

   CeeLo has no category for one, so it is labelled a curve transfer (no Monte Carlo); unlike a
   genuine curve transfer it carries its own near-field table, and no transfer sigma model.  Its
   full-energy-peak efficiency is tied to the geometry it was built with; it carries a total
   efficiency only once #attachTotalEfficiency has given it one.
   */
  bool isGridResponse( const std::shared_ptr<const ceelo::DetectorResponse> &resp );

  /** The geometry `makeDrf` models a DETECTOR.txt record as.  Field by field:
     - D1, D2: crystal radius (D1/2) and length.  D3: only classifies the record (#DetectorDef::Kind).
     - Slots 0/3/6: the front/side/back Ge dead layers (`ceelo::DeadLayerConfig`).
     - Slots 1, 2: front shells; slots 4, 5: side shells; slots 7, 8: back plates, inside the side
       shells (`ceelo::GeometryDescriptor::set_layer_back_thickness_cm`).
     - D7: the crystal face sits D7 behind the endcap's outer face, so the front vacuum shell is
       max(0, D7 - slots 1 and 2), and the endcap-front offset (the DRF's setback) is D7, or the
       front stack when that is thicker.  Measured: CeeLo near-field Monte Carlo against the grid
       (0.5-2 cm vs 25 cm, on axis) over 65 grids fits to a median 0.9% (D3 = 0) and 1.3-2.5%
       (D3 > 0) this way.  Reading D6 as a vacuum gap in front of the stack - the former rule - put
       a standard coax ~2 mm too deep, 5-11% low at 0.5-2 cm; reading D6 as the depth on a
       thin-window record is 9-17% off (median).
     - Side vacuum shell: max(0, D7 - D6), which is zero on every D3 = 0 record.  Unconfirmed: side-on
       Monte Carlo cannot tell candidate side gaps apart.
     - D5: the endcap length, from its outer face.  Whatever it leaves behind the crystal and the
       back plates is vacuum, so the back of the modelled can is where the record puts it.
     - D4: not used, so the modelled can is the crystal plus its side shells - narrower than
       the real one (71.6 vs 87.5 mm on one corpus detector), which the grid's no-data zone puts at
       D4.
     The record has no bore or bullet, so coaxial families get a conservative guessed bore - which
   is why users may edit the geometry afterwards (see DetectorPeakResponse::monteCarloGeometry).
   Unmodelable layers are noted in `warnings`; an unmodelable front layer's thickness becomes vacuum,
   so the crystal stays at D7.
   Throws std::runtime_error without a crystal diameter and length.
   */
  ceelo::GeometryDescriptor geometryFromDetectorDef( const DetectorDef &def,
                                                     std::vector<std::string> &warnings );

  /** The one distance (cm, from the endcap face) at which an imported grid stands in for a
   measurement - transfer anchors, grounding a Monte Carlo to it, and pinning an attached total:
   max( 50 cm, 10 x the transverse half-extent ), a calibration-like distance clear of the near
   field.  (The legacy far-field curve a grid DRF also carries is sampled ~1000 half-extents out -
   tens of metres - which is no place to anchor anything.)
   */
  double gridReferenceDistanceCm( const ceelo::GeometryDescriptor &gd );

  /** On-axis absolute FEP efficiencies of a grid response at `dist_from_face_cm` (<= 0 for
   #gridReferenceDistanceCm): 24 log-spaced energies over its validated range, plus flanks either
   side of each crystal K-edge.  Grounding a Monte-Carlo response to these makes it agree with the
   grid; they are also the grid's transfer anchor (CeeLoUtils::gridTransferAnchor).
   */
  std::vector<ceelo::GroundingPoint> groundingPoints( const ceelo::DetectorResponse &grid,
                                                      double dist_from_face_cm = -1.0 );

  /** A copy of grid response `grid` carrying a total efficiency from `donor` - or none, for a null
   `donor`.  The full-energy-peak efficiency is untouched (bit-identical).

   The total is the grid's FEP times the donor's own total/FEP at the same physical point.  The
   donor - typically a Monte Carlo of the user's (possibly edited) geometry - supplies only that
   total-to-peak ratio, in which whatever its geometry gets wrong largely cancels; the grid stays the
   authority for the FEP, and any grounding the donor carries is ignored.  It is stored on the
   GRID's response:
     - an `EtaTotTable` on the grid's own energy x angle nodes, at the far-field pin
       (max(1000 a, 100 cm) from the crystal face); and
     - everywhere closer, a distance table (`ceelo::TotEffPayload::near_field`) holding the same
       product on the GRID's own close-in lattice - its energies, angle columns and distance
       ladder - so neither geometry's kernel decides the total, the grid's close-in FEP and the
       donor's close-in total/FEP both carry over, and the grid's close-in structure, interpolated
       the same way for both, cancels from total/FEP.  A rung makeDrf reads at the endcap face in
       front of it takes the donor's ratio there too.
   Exact at the lattice nodes; between them, total/FEP is the donor's ratio to <=0.05% from 1 cm
   out, and 0.13% at worst at the grazing near-contact corner, on the mock pair
   (test_ImportedGridDrf's `TotalOnEditedGeometry` and `TotalToPeakBetweenNodes`).  That table is
   as large as the grid's own, so the response serializes to ~2.5x the bare grid (1.1 -> 2.8 MB on
   the mock), and attaching it takes ~5 s CPU there.  Where
   nothing reaches a crystal - a kernel underflowing to 0, in either geometry - the ratio comes from
   the donor's kernel-free prefactors, so the total is 0 wherever the grid's FEP is.

   End to end on the mock pair, the app's own run (General profile, scans stopped on total
   precision - ceelo::GenerationOptions::scans_stop_on_total - at 0.3%) against direct point-source
   Monte Carlo of the geometry the grid was made from (1-27 cm, 0-75 deg, 60-2614 keV, 0.2% MC),
   measured when the distance table was a coarser lattice of its own (which added up to ~1% between
   its nodes):
     - that geometry: worst 1.65%, RMS 0.65% (carrying the grounded donor's total instead: worst
       1.58%, RMS 0.55% - the total/FEP takes on the off-axis and close-in FEP's MC noise, 0.3-1.2%
       at those scans);
     - the importer's guessed core: worst 2.43%, RMS 0.93% (the grounded donor's total: worst 7.6%,
       RMS 2.7%).
   And the grounded total failed outright for a geometry far from the grid's - 1 cm of lead in front
   put totals above 1 between grounding energies - where total/FEP stays 1.05-1.26 from 20 to
   120 keV (2-25 cm, on axis and at 75 deg).

   Throws std::runtime_error if `donor` has no total efficiency, or `grid` carries a grounding
   (which would scale its FEP too).
   */
  std::shared_ptr<ceelo::DetectorResponse> attachTotalEfficiency(
                                                    const ceelo::DetectorResponse &grid,
                                                    const ceelo::DetectorResponse *donor );

  /** How closely a Monte Carlo of the user's geometry reproduces the imported grid's FEP close in,
   where grounding on axis at #gridReferenceDistanceCm does not tie them together.
   */
  struct GridFepConsistency
  {
    size_t num_nodes = 0;
    double mean_ln = 0.0;    //mean of ln( Monte-Carlo FEP / grid FEP )
    double rms_ln = 0.0;     //RMS of the same
    double worst_ln = 0.0;   //largest-magnitude value, signed
    double worst_energy_keV = 0.0, worst_dist_cm = 0.0, worst_cos_theta = 1.0;  //from the endcap-face centre
  };//struct GridFepConsistency

  /** Compares `donor` (grounded to `grid` with #groundingPoints) with the grid's FEP at the donor's
   near-field Monte-Carlo nodes - where its FEP is the node's own MC - that lie in front of the
   endcap face and inside the grid's validated energy range.  A large value means the geometry does
   not reproduce the vendor's characterization there, and the Monte-Carlo total close in deserves
   the same doubt.  Empty (`num_nodes` 0) when the donor has no near-field table.

   On the mock pair, from the app's run (whose near-field FEP is only 0.3-1.2% precise): RMS 1.0%,
   worst 4.5%, for the geometry the grid was made from; RMS 3.9%, worst 23%, for the importer's
   guessed core - whose total was then 3.7-8.1% off at worst, against <=1.6%.
   */
  GridFepConsistency gridFepConsistency( const ceelo::DetectorResponse &grid,
                                         const ceelo::DetectorResponse &donor );

  /** Convenience: parse both files from disk, select the matching geometry
   record, and build the DRF.  Throws std::runtime_error on failure.
   */
  std::shared_ptr<DetectorPeakResponse> makeDrfFromFiles(
                                          const std::string &parPath,
                                          const std::string &detectorTxtPath );
}//namespace DetEffG2kPar

#endif //DetectorEffG2kPar_h
