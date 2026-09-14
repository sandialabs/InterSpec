#ifndef DrfModifyCalc_h
#define DrfModifyCalc_h
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

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class MeasuredDrfPoints;
class DetectorPeakResponse;

/** The non-GUI core of the "Modify Detector Response" tool (`DrfModifyWidget`): turning what the
 Anchor tab holds into a DRF, and checking that the DRF that comes out is one self-consistent
 description of a detector.  Kept free of Wt so it can be unit tested - the defects this code exists
 to prevent were all "the emitted DRF's pieces disagree", which no GUI test can see and every
 apply-path test can.

 ## Which store is authoritative

 A DRF can say what its efficiency is in three places and what its uncertainty is in three.  This is
 decided once (see the table on `DetectorPeakResponse::efficiencyUncert`) and everything here obeys
 it:

  - the value comes from an attached `ceeloResponse()`, else the `efficiencyCurve()`;
  - the uncertainty comes from that same response, else - for an equation curve - the coefficient
    covariance, else the node covariance;
  - `measuredPoints()` are the raw measurements: the input a curve is fit from and a Monte-Carlo
    response is grounded to.  They are never read as a value, and never edited apart from the curve
    they produced (see #AnchorEditor::RefitPoints).

 ## What must hold after an apply

 #checkDrfSelfConsistent is the executable form of this list, and every apply path here is expected
 to leave it satisfied:

  1. the curve is valid and its energy units are positive;
  2. a stored coefficient covariance is usable (symmetric, positive semi-definite) and its rank
     matches the curve's coefficient count - otherwise `DetectorEfficiencyCurve::fracCovariance`
     silently ignores it and the user's edit does nothing;
  3. a stored node covariance is usable;
  4. if the curve is an equation and raw points are stored, the curve is still what those points fit
     to - i.e. the refit source has not stopped describing the curve;
  5. `measuredPoints()` distances make sense for the geometry type.

 Staleness of an attached response is NOT part of this: it depends on what the response was built
 from, which only the editing session knows - see #seedFingerprint.
 */
namespace DrfModifyCalc
{
  /** Which editor the Anchor tab offers, decided by what the DRF actually carries - never by the
   Flat-Disk / Geometry-Modeled toggle, which says how the detector answers off-axis questions and
   has nothing to do with how its efficiency is represented.

   The rule that matters: **the rows are seeded from whatever the apply writes back.**  Getting that
   wrong is what let a Create-DRF detector's absolute measured efficiencies be installed as its
   intrinsic curve, low by one solid angle.
   */
  enum class AnchorEditor
  {
    /** An equation curve that carries the raw points it was fit from (every "Create Detector
     Response Function" detector).  Rows are those points - ABSOLUTE efficiency at each point's own
     source distance - and applying re-fits the equation from them, so the curve, its coefficient
     covariance, the node covariance and the points all come out of one fit.
     */
    RefitPoints,

    /** A `kEnergyEfficiencyPairs` curve of ABSOLUTE efficiencies at one reference distance, which
     also carries those pairs as raw points (an ANGLE .outx import).  Rows are the pairs, which here
     are the points; applying rewrites both, plus the block covariance the points imply.
     */
    AbsolutePoints,

    /** A `kEnergyEfficiencyPairs` curve with no raw points (a GADRAS Efficiency.csv, an ISOCS .ecc).
     Rows are the curve's own numbers, in its own geometry convention; applying replaces the curve
     and the node covariance built from the rows' uncertainty columns.
     */
    CurvePairs,

    /** A `kExpOfLogPowerSeries` curve with no raw points: the coefficients and their covariance,
     which for this form IS the uncertainty every fit propagates.
     */
    Coefficients,

    /** A `kFunctialEfficienyForm` curve: the formula text.  Its point rows are covariance nodes
     only - the efficiency comes from the formula - so the efficiency column is not shown.
     */
    Formula
  };//enum class AnchorEditor

  AnchorEditor editorForDrf( const DetectorPeakResponse &drf );

  /** Whether this editor's rows come from (and are written back to) `measuredPoints()`. */
  bool editorUsesMeasuredPoints( const AnchorEditor editor );

  /** Whether this editor shows the point table at all. */
  bool editorUsesPointTable( const AnchorEditor editor );


  /** Something an apply could not do, or did in a way the user needs to know about.

   Deliberately an id rather than text: the calc layer has no Wt dependency, and the wording lives in
   `InterSpec_resources/app_text/DrfModifyWidget.xml` with everything else the user reads.  `blocking`
   distinguishes "your edit was not applied, here is why" from "it was applied, and here is what else
   happened" - the first refuses the apply, rather than closing the dialog looking successful.
   */
  struct Problem
  {
    std::string messageId;
    std::string arg;              //substituted for {1} when non-empty
    bool blocking = true;

    Problem( const std::string &id, const std::string &argument = std::string(),
             const bool is_blocking = true )
      : messageId( id ), arg( argument ), blocking( is_blocking ){}
  };//struct Problem

  /** Whether any of `problems` refuses the apply. */
  bool anyBlocking( const std::vector<Problem> &problems );


  /** One row of the Anchor tab's point table, parsed.  The widget owns the text-to-number step (it
   knows which row a cell is in, for the message); this is what the apply paths see.
   */
  struct PointRow
  {
    /** 1-based row number as shown, for messages. */
    int rowNumber = 0;

    /** Energy, in `AnchorOptions::energyUnits`. */
    double energy = 0.0;

    /** Efficiency as the editor means it: absolute at #distance for the measured-point editors, the
     curve's own value for #AnchorEditor::CurvePairs, ignored for #AnchorEditor::Formula.
     */
    double efficiency = 0.0;

    /** Fractional 1-sigma, independent between rows; 0 when the cell was left blank. */
    double fracStat = 0.0;

    /** Fractional 1-sigma, correlated: within a source for the measured-point editors, across
     energy (by `AnchorOptions::corrLength`) otherwise. */
    double fracCert = 0.0;

    /** Source-to-face distance, PhysicalUnits; <= 0 means "the reference distance". */
    double distance = -1.0;

    std::string sourceKey;

    /** Index of the `measuredPoints()` entry this row was seeded from, or -1 for a row the user
     added.  Matched by position, NOT by energy: matching on energy loses a point's provenance (peak
     area, live time, file name, distance uncertainty) the moment its energy is edited.
     */
    int seedIndex = -1;
  };//struct PointRow


  /** The Anchor tab's settings that are not per-row. */
  struct AnchorOptions
  {
    /** PhysicalUnits energy unit the rows' energies are written in. */
    float energyUnits = 1.0f;

    /** Distance an absolute curve is anchored at, PhysicalUnits; <= 0 when not applicable. */
    double refDistance = -1.0;

    /** Log-energy correlation length the "Corr. %" column is correlated across energy with - the
     `EccUncertOptions::effectiveCorrLength` semantics, so one number carries all three modes:
     <= 0 is "Uncorrelated" (that column becomes diagonal too),
     `DetectorEfficiencyUncert::sm_fullyCorrelatedLength` is "Fully correlated", and anything between
     is the Gaussian log-energy kernel of that width.  Used as given - the two ends of the range are
     opposite statements, so there is no defaulting here.  Does not affect the "Stat. %" column,
     which is independent between rows by definition. */
    double corrLength = -1.0;

    /** Number of equation terms for a refit; <= 0 keeps the current count. */
    int equationTerms = 0;

    /** Whether the editor shows a row table at all.  Distinguishes "the user emptied every row"
     (which states no uncertainty) from "this caller has no rows to offer" (which keeps whatever the
     curve had) - a distinction the row count cannot make, since emptying the rows produces none. */
    bool hasRowTable = false;
  };//struct AnchorOptions


  /** Applies the point table to `working` per `editor`, which must be one of the three editors that
   uses it.  Returns whether anything was written; `problems` always says why not.

   `seedPoints` is the DRF's points as the dialog was opened with, which rows carry provenance from.
   */
  bool applyPointRows( DetectorPeakResponse &working,
                       const AnchorEditor editor,
                       const std::vector<PointRow> &rows,
                       const AnchorOptions &options,
                       const std::shared_ptr<const MeasuredDrfPoints> &seedPoints,
                       std::vector<Problem> &problems );

  /** Applies the coefficient editor: the equation, and (when `writeCovariance`) the covariance
   composed from `sigmas` and `rho` - refusing an impossible one rather than storing it.

   @param sigmas  1-sigma of each coefficient (absolute; these coefficients are logs).
   @param rho     Row-major N*N correlation matrix, unit diagonal.
   @param writeCovariance  False leaves whatever uncertainty the DRF already had, for the case where
          the matrix shown was only ever a placeholder the user did not touch (see the note on
          #sigmaRhoFromLegacyUncerts).
   */
  bool applyCoefficients( DetectorPeakResponse &working,
                          const std::vector<float> &coefficients,
                          const std::vector<double> &sigmas,
                          const std::vector<double> &rho,
                          const float energyUnits,
                          const bool writeCovariance,
                          std::vector<Problem> &problems );

  /** Applies the formula editor: the formula, plus the node covariance its rows describe (or the
   existing uncertainty when there are no rows and no default).
   */
  bool applyFormula( DetectorPeakResponse &working,
                     const std::string &formula,
                     const float energyUnits,
                     const std::vector<PointRow> &rows,
                     const AnchorOptions &options,
                     std::vector<Problem> &problems );


  /** Row-major N*N covariance `C[i][j] = sigma_i * sigma_j * rho[i][j]`, with `rho`'s diagonal taken
   as 1 whatever it holds.  Empty when the sizes do not agree.

   The editor keeps sigma and rho separately, rather than a covariance, so that typing a 0 sigma does
   not destroy that coefficient's correlations (you cannot recover rho from a row of zeros) and so
   the correlations the user can see are exactly the ones stored.
   */
  std::vector<double> covarianceFromSigmaRho( const std::vector<double> &sigmas,
                                              const std::vector<double> &rho );

  /** The inverse split, for seeding the editor from a stored covariance: sigma from the diagonal,
   rho from the off-diagonal (0 where either sigma is 0).
   */
  void sigmaRhoFromCovariance( const std::vector<float> &covRowMajor,
                               std::vector<double> &sigmas,
                               std::vector<double> &rho );

  /** Sigmas from the legacy per-coefficient uncertainties, with rho all zero.

   This is a PLACEHOLDER, not a covariance the DRF has: fitted log-power-series coefficients are
   strongly correlated, so assuming independence overstates the efficiency band badly.  It exists so
   the editor can show something editable instead of a grid of zeros, and the caller must only store
   it if the user actually edited it.  Returns false when there is nothing to seed from.
   */
  bool sigmaRhoFromLegacyUncerts( const DetectorPeakResponse &drf,
                                  std::vector<double> &sigmas,
                                  std::vector<double> &rho );


  /** A fingerprint of everything an attached CeeLo response is built from: the efficiency curve, both
   covariances, the raw measured points, and the geometry.

   Staleness is derived from this rather than tracked with a "something changed" flag.  A flag has to
   be cleared by whoever regenerates, and every path that forgets to - an automatic rebuild that
   anchored on a pre-edit seed, a redo that restored the edits but not the flag - leaves an edited
   curve behind a response that still answers every query.  Comparing content cannot forget.
   */
  std::size_t seedFingerprint( const DetectorPeakResponse &drf );


  /** Checks the editing invariants listed at the top of this file.  Returns true when the DRF is
   self-consistent; otherwise `why` gets a one-line explanation.

   Called on what the Modify tool emits (as a developer check), and directly by the unit tests.
   */
  bool checkDrfSelfConsistent( const DetectorPeakResponse &drf, std::string &why );


  /** The efficiency uncertainty a DRF reports at one energy, split into what its own data supports
   and what is an ad hoc model envelope (`ceelo::model_sigma` - a transfer model's off-axis
   allowance, a regime floor; things no measurement of this detector constrains).

   Comes from `DetectorPeakResponse::efficiencyFracCovariance`, i.e. the same numbers the
   activity/shielding fit propagates and the efficiency chart draws.
   */
  struct UncertSummary
  {
    bool valid = false;
    float energy = 0.0f;    //keV the summary is for
    double total = 0.0;     //fractional 1-sigma
    double model = 0.0;     //the part that is an ad hoc envelope; <= total
    double data = 0.0;      //sqrt(total^2 - model^2)

    /** True when the DRF states no efficiency uncertainty of its own at all, yet reports one.

     Every shipped GADRAS `Detector.dat` is in this position: it says nothing about how well its
     efficiency is known, InterSpec attaches a curve-transfer response to it at load, and
     `CeeLoUtils::sm_default_anchor_frac_sigma` then supplies a flat 5% that the response carries in
     the slot meant for measured anchor uncertainty - so it arrives in #data, where it looks like a
     measurement.  It is not one: it will not change if the underlying efficiency is excellent or
     terrible.  Callers must say so rather than presenting #data as what the detector's data supports.
     */
    bool dataIsAssumed = false;
  };//struct UncertSummary

  /** @param energy keV; <= 0 picks a representative energy inside the DRF's range. */
  UncertSummary uncertSummary( const DetectorPeakResponse &drf, const float energy = -1.0f );
}//namespace DrfModifyCalc

#endif //DrfModifyCalc_h
