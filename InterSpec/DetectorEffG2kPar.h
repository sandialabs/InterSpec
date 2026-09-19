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

   * AT a node (45/60/80/100/300/500 here) we reproduce ".ecc" to ~6 digits, so
     both tools demonstrably hold the same stored values.
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

 -----------------------------------------------------------------------------
  LOW-ENERGY VALIDITY FLOOR: TRUST THE ASSEMBLED RESPONSE ABOVE ~45 keV
  (measured; a property of the files, not of our code)
 -----------------------------------------------------------------------------
 The grid ITSELF is exact at its own nodes at every energy - the round-trip test
 reproduces stored values to ~6 digits down to the lowest node.  What degrades at
 low energy is the ASSEMBLED CeeLo response, which does not store the grid: it
 factorizes efficiency into an analytic kernel times a PCHIP-in-ln(E) correction
 (see step 7 of `makeDrf`).  Below the absorption knee the true efficiency falls
 by ~5 DECADES over a handful of nodes (on 18211381, on-axis at 25 cm: 8.2e-11 at
 12 keV, 1.4e-05 at 22 keV, 1.6e-03 at 45 keV).  Twenty log-spaced nodes cannot
 resolve a cliff that steep, and the analytic kernel's endcap cutoff does not
 have the same shape as the file's, so the correction term carries the entire
 mismatch: its on-axis node values swing over ~4.6 in ln (0.74 -> -3.81 -> -1.31
 -> -0.81 -> -0.09) across 10-32 keV before settling to ~0.1 for the whole rest
 of the range.  PCHIP through that whipsaws between nodes.

 Measured worst |error| of the assembled response against the grid it was built
 from, over theta 0-88 deg and d in {2, 5, 10, 25, 100} cm (detector 18211381,
 whose grid starts at 10 keV):

     energy band     worst error      efficiency there (rel. to its peak)
     10 -  16 keV    astronomical     ~1e-8   (efficiency is ~1e-10; meaningless)
     16 -  22 keV    ~530%            ~1e-5 - 4e-3
     22 -  32 keV    ~91%             0.004 - 0.10
     32 -  45 keV    ~17%             0.10 - 0.52
     45 -  60 keV    ~5.2%            0.52 - 0.80
     60 - 2000 keV   ~2.8%            0.80 - 1.00   <-- the usable range

 The error is confined to where the absolute efficiency is a small fraction of
 its peak, i.e. it is large in RELATIVE terms exactly where the response is
 negligible in ABSOLUTE terms.  On LAB06 (whose grid starts at 45 keV) there is
 no such region at all: worst ~2.2% across every band.

 Two things were checked and ruled out as causes, so this is not a latent bug:
   * 18211381's two lowest nodes (10 and 12 keV) are BIT-IDENTICAL across all
     32120 grid cells - the file duplicates its first node.  Dropping the
     duplicate changes the errors by exactly nothing, so the degeneracy is not
     the mechanism.
   * The behaviour is not near-field or grazing-angle specific: it is just as
     large on-axis at 25 cm and 100 cm as at 2 cm.
 It is the energy factorization alone.

 CONSEQUENCE, and why this is acceptable: 20 keV is the lower end of trustworthy
 nuclear data in InterSpec generally, and reference-tool output ("*.ecc") is never
 tabulated below 45 keV anywhere in the validation corpus - there is no ground
 truth below that to validate against even in principle.  So:

     >= 60 keV     trust to ~3% (the between-node figure above)
     45 - 60 keV   trust to ~5%
     32 - 45 keV   indicative only (~17%)
     <  32 keV     DO NOT USE the assembled response

 Callers needing values below ~45 keV should query `ParEfficiency::efficiency`
 directly, which reads the grid and is exact at its nodes.  Fixing the assembled
 path would mean changing the energy representation (more nodes through the knee,
 or a kernel whose low-energy cutoff matches the file's) - deliberately NOT done
 here, since it would trade a documented limit outside InterSpec's quantitative
 range for churn in the well-validated 60-7000 keV region.
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
