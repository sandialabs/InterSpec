#ifndef DetectorEfficiency_h
#define DetectorEfficiency_h
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

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <functional>

#include "InterSpec/DetectorPeakResponse.h"

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml


/** A single evaluator-specified efficiency uncertainty band.

 Semantics: the efficiency uncertainty is 100% correlated for energies within
 [lowerEnergy, upperEnergy), and uncorrelated with energies in other bands.
 Energies are in keV (PhysicalUnits); the uncertainty is fractional
 (i.e., 0.05 == 5%).

 This historical, piecewise way of specifying DRF uncertainty comes from
 evaluations/specifications used by activity-fitting applications that did
 not fit shielding.
 */
struct EffUncertBand
{
  float lowerEnergy;
  float upperEnergy;
  float fractionalUncert;

  bool operator==( const EffUncertBand &rhs ) const;
};//struct EffUncertBand


/** Optional uncertainty description for a detector efficiency curve.

 Holds two independent, optional components:
 - Piecewise bands (see #EffUncertBand): user/evaluator specified.
 - Node covariance: the canonical, rigorous form.  A set of energy nodes
   {E_j} plus a covariance matrix C of the *fractional* efficiency error at
   those nodes.  The fractional error at any other energy is *defined* as the
   linear (in log-energy) interpolation of the node errors, with constant
   extrapolation beyond the first/last node.  So the covariance among any
   requested energies is L*C*L^T (L being the interpolation operator), which
   is positive-semidefinite by construction whenever C is.

 The node covariance is independent of the functional form of the efficiency
 curve, so it can be populated from a MakeDrf fit (via the fit-coefficient
 covariance, which may optionally be retained alongside as provenance), from
 ANGLE/ISOCS style per-point uncertainties, or from external evaluations.

 Instances are held as `shared_ptr<const DetectorEfficiencyUncert>` and treated as
 immutable once attached to a DetectorPeakResponse - replace, dont mutate.

 All energies passed to/from this class are in keV (PhysicalUnits keV = 1.0),
 regardless of the energy units of the efficiency curve it describes.
 */
class DetectorEfficiencyUncert
{
public:
  DetectorEfficiencyUncert();

  bool hasBands() const;
  bool hasNodeCovariance() const;
  bool isEmpty() const;

  /** Covariance of the fractional efficiency error among `energies` (keV),
   from the node-covariance component only (no piecewise band contribution).

   Computed as L*C*L^T where L linearly interpolates, in log-energy, the node
   errors to the requested energies (constant extrapolation outside the node
   range).  Returned row-major, with size energies.size()*energies.size().

   Returns all zeros if no node covariance is defined.
   */
  std::vector<double> nodeFracCovariance( const std::vector<double> &energies ) const;

  /** Total fractional-efficiency-error covariance among `energies` (keV):
   nodeFracCovariance(...) plus the piecewise-band contribution (u_b^2 added
   to element (i,j) when both energies fall within the same band).
   */
  std::vector<double> efficiencyFracCovariance( const std::vector<double> &energies ) const;

  /** Square root of the diagonal of efficiencyFracCovariance(...) - the
   1-sigma fractional uncertainty envelope at each requested energy (keV).
   */
  std::vector<double> fracUncertainties( const std::vector<double> &energies ) const;

  /** Builds a node covariance from per-point fractional uncertainties, using
   a Gaussian log-energy correlation model:
     rho(E1,E2) = exp( -0.5 * pow( (log(E1) - log(E2))/corrLength, 2.0 ) )
   so C[j][k] = u_j * u_k * rho(E_j,E_k); this kernel is positive-semidefinite
   for any corrLength > 0.

   @param energies Node energies in keV; will be sorted (and exact duplicates
          removed, keeping the first).  All must be > 0.
   @param fracUncerts Fractional 1-sigma uncertainties (>= 0), one per energy.
   @param corrLength Correlation length in natural-log-energy units.  A value
          <= 0 gives an uncorrelated (diagonal) covariance; very large values
          approach 100% correlation.

   Throws std::runtime_error on invalid input.
   */
  static std::shared_ptr<DetectorEfficiencyUncert> fromPointUncerts(
                              const std::vector<float> &energies,
                              const std::vector<float> &fracUncerts,
                              const double corrLength = sm_defaultLogEnergyCorrLength );

  /** Sets the piecewise bands; validates bands are sorted by energy,
   non-overlapping, with upperEnergy > lowerEnergy and fractionalUncert >= 0.
   Throws std::runtime_error on invalid input.
   */
  void setBands( const std::vector<EffUncertBand> &bands );

  /** Sets the node covariance.

   @param energies Node energies in keV, strictly increasing, all > 0;
          at most #sm_maxCovarianceNodes entries.
   @param covRowMajor Row-major N*N covariance of fractional efficiency error
          at the nodes.  Must be symmetric (to within small numerical error -
          it is symmetrized on storage), with non-negative diagonal.

   Throws std::runtime_error on invalid input.
   */
  void setNodeCovariance( const std::vector<float> &energies,
                          const std::vector<float> &covRowMajor );

  /** Sets the (optional) covariance of the efficiency-equation fit
   coefficients - retained for provenance only; not used in any calculation
   here.  Must be a square (M*M) row-major matrix, or empty.
   */
  void setCoefficientCovariance( const std::vector<float> &covRowMajor );

  const std::vector<EffUncertBand> &bands() const;
  const std::vector<float> &covarianceEnergies() const;

  /** Row-major N*N node covariance matrix; empty if not defined. */
  const std::vector<float> &covarianceMatrix() const;

  /** Row-major M*M fit-coefficient covariance; empty if not defined. */
  const std::vector<float> &coefficientCovariance() const;

  /** Log-energy correlation length used to construct the node covariance from
   per-point uncertainties; <= 0 if the covariance was not built that way.
   */
  double correlationLength() const;

  /** Appends a "EfficiencyUncert" node under `parent`. */
  void toXml( ::rapidxml::xml_node<char> *parent, ::rapidxml::xml_document<char> *doc ) const;

  /** Parses a "EfficiencyUncert" node; throws std::runtime_error on error. */
  void fromXml( const ::rapidxml::xml_node<char> *node );

  /** Appends url query-string entries (keys prefix+"EFUB", prefix+"EFUE",
   prefix+"EFUC", prefix+"EFUL") to `parts`.  The covariance matrix is encoded
   as its upper triangle (including diagonal), N*(N+1)/2 values.
   The coefficient covariance is never written to URLs.
   */
  void toUrlParts( std::map<std::string,std::string> &parts, const std::string &prefix ) const;

  /** Creates an instance from url query-string entries; returns nullptr if no
   relevant keys present.  Throws std::runtime_error on malformed entries.
   */
  static std::shared_ptr<DetectorEfficiencyUncert> fromUrlParts(
                              const std::map<std::string,std::string> &parts,
                              const std::string &prefix );

  /** Combines all defined content into `seed` via boost::hash_combine. */
  void appendToHash( std::size_t &seed ) const;

  /** Exact (bitwise) comparison of all members. */
  bool operator==( const DetectorEfficiencyUncert &rhs ) const;

#if( PERFORM_DEVELOPER_CHECKS )
  /** Tests whether the two objects are equal to within small numerical
   rounding (e.g., after a serialization round-trip, where float parsing may
   be off in the last ULP).  Throws std::runtime_error with a brief
   explanation when an issue is found.
   */
  static void equalEnough( const DetectorEfficiencyUncert &lhs,
                           const DetectorEfficiencyUncert &rhs );
#endif

  /** Default log-energy correlation length used when ingesting per-point
   uncertainties (ANGLE .outx, efficiency CSVs).  A value of 0.5 in ln(E)
   means: energies 20% apart have correlation ~0.94, while energies a factor
   of e (~2.72) apart have correlation exp(-2) ~ 0.14.

   These per-point uncertainties typically come from transfer/MC calculations
   whose errors are dominated by smooth, systematic model/geometry effects, so
   nearby energies are strongly correlated - treating them as independent
   would badly overestimate the freedom of the efficiency curve.

   The value actually used is serialized with the covariance, so changing this
   default will not alter previously stored DRFs.
   */
  static const double sm_defaultLogEnergyCorrLength;

  /** Maximum number of covariance node energies accepted. */
  static const size_t sm_maxCovarianceNodes;

private:
  /** Piecewise bands; sorted by energy, non-overlapping. */
  std::vector<EffUncertBand> m_bands;

  /** Node energies in keV; sorted ascending. */
  std::vector<float> m_covEnergies;

  /** Row-major N*N covariance of fractional efficiency error at the nodes. */
  std::vector<float> m_covMatrix;

  /** Optional provenance: row-major M*M covariance of efficiency-equation fit
   coefficients.
   */
  std::vector<float> m_coefCovMatrix;

  /** Log-energy correlation length used by #fromPointUncerts; <= 0 if the
   node covariance was set directly.
   */
  double m_corrLength;
};//class DetectorEfficiencyUncert


/** A reusable efficiency-curve, mirroring the three efficiency
 representations DetectorPeakResponse supports (energy/efficiency pairs,
 formula string, exp-of-log-power-series), plus an optional uncertainty.

 Currently used to optionally store the *total* efficiency of a detector
 (the probability an incident gamma deposits *any* energy - not just its
 full energy), to allow accounting for cascade-summing effects.

 Evaluation delegates to the same static helpers the full-energy efficiency
 of DetectorPeakResponse uses, so no logic is duplicated.

 Instances are held as `shared_ptr<const DetectorEfficiencyCurve>` and treated as
 immutable once attached to a DetectorPeakResponse - replace, dont mutate.
 */
class DetectorEfficiencyCurve
{
public:
  DetectorEfficiencyCurve();

  /** Returns true if an efficiency representation has been set. */
  bool isValid() const;

  /** Evaluates the efficiency at `energy` (in keV, PhysicalUnits).

   The interpretation (intrinsic vs per-decay) follows the geometry type of
   the owning DetectorPeakResponse; no absolute-to-intrinsic or solid angle
   corrections are applied here.

   Throws std::runtime_error if !isValid().
   */
  float efficiency( const float energy ) const;

  /** Sets from energy/efficiency pairs.
   @param pairs At least two entries; will be sorted by energy.  Energies in
          units of `energyUnits`; efficiencies must be >= 0.
   @param energyUnits The PhysicalUnits energy unit of the pair energies
          (e.g. 1.0 for keV, 1000.0 for MeV).
   Throws std::runtime_error on invalid input.
   */
  void setFromPairs( const std::vector<DetectorPeakResponse::EnergyEfficiencyPair> &pairs,
                     const float energyUnits );

  /** Sets from a formula string with energy as 'x'; see
   DetectorPeakResponse::setIntrinsicEfficiencyFormula for syntax.
   Throws std::runtime_error if the formula is invalid.
   */
  void setFromFormula( const std::string &formula, const float energyUnits );

  /** Sets from exp(A0 + A1*log(x) + A2*log(x)^2 + ...) coefficients, with
   energy evaluated in units of `energyUnits`.

   @param uncerts Optional per-coefficient 1-sigma uncertainties; must be
          empty or the same size as `coefs`.

   Throws std::runtime_error on empty coefficients or mismatched uncerts.
   */
  void setFromExpOfLogPowerSeries( const std::vector<float> &coefs,
                                   const std::vector<float> &uncerts,
                                   const float energyUnits );

  /** Clears to an invalid/empty curve. */
  void reset();

  /** Sets the representation directly from the raw fields used by
   DetectorPeakResponse serialization (database columns, XML, app-url).

   Tolerates `form == kNumEfficiencyFnctForms` (yields an invalid curve);
   rebuilds the parsed formula function when `form == kFunctialEfficienyForm`.
   Does NOT touch the (separately serialized) rich uncertainty.

   Throws std::runtime_error on inconsistent input (e.g. a non-empty formula
   that fails to parse, or uncerts whose size doesnt match the coefficients).
   */
  void setRawFields( const DetectorPeakResponse::EfficiencyFnctForm form,
                     const float energyUnits,
                     const std::vector<DetectorPeakResponse::EnergyEfficiencyPair> &pairs,
                     const std::string &formula,
                     const std::vector<float> &coefs,
                     const std::vector<float> &coefUncerts );

  DetectorPeakResponse::EfficiencyFnctForm form() const;
  float energyUnits() const;
  const std::vector<DetectorPeakResponse::EnergyEfficiencyPair> &energyEfficiencies() const;
  const std::string &formula() const;
  const std::function<float(float)> &efficiencyFcn() const;
  const std::vector<float> &expOfLogPowerSeriesCoeffs() const;

  /** Per-coefficient 1-sigma uncertainties for the exp-of-log form; empty if
   not set.  (This is the legacy diagonal uncertainty, distinct from the rich
   #DetectorEfficiencyUncert returned by #uncertainty.)
   */
  const std::vector<float> &expOfLogPowerSeriesUncerts() const;

  std::shared_ptr<const DetectorEfficiencyUncert> uncertainty() const;
  void setUncertainty( std::shared_ptr<const DetectorEfficiencyUncert> uncert );

  /** Returns a copy of this curve with every efficiency value multiplied by
   the energy-independent constant `factor` (> 0).

   For pairs the efficiencies are scaled; for exp-of-log the constant term is
   shifted by log(factor); for a formula the parsed function is wrapped (the
   stored formula string is left unchanged - matching legacy behavior).  The
   (fractional) uncertainty is carried across unchanged.

   Throws std::runtime_error if the curve is invalid or factor is not > 0.
   */
  DetectorEfficiencyCurve scaledByConstant( const double factor ) const;

  /** Appends a node named `node_name` under `parent`. */
  void toXml( ::rapidxml::xml_node<char> *parent, ::rapidxml::xml_document<char> *doc,
              const char *node_name ) const;

  /** Parses a node previously written by #toXml; throws on error. */
  void fromXml( const ::rapidxml::xml_node<char> *node );

  /** Appends url query-string entries, keys prefixed by `prefix` (e.g.,
   prefix "T" gives "TEFT", "TEFX", "TEFY", "TEFE", "TEFC", "TEUNIT", and
   "TEFUB"/"TEFUE"/"TEFUC"/"TEFUL" for the uncertainty).
   */
  void toUrlParts( std::map<std::string,std::string> &parts, const std::string &prefix ) const;

  /** Creates an instance from url query-string entries; returns nullptr if no
   relevant keys present.  Throws std::runtime_error on malformed entries.
   */
  static std::shared_ptr<DetectorEfficiencyCurve> fromUrlParts(
                              const std::map<std::string,std::string> &parts,
                              const std::string &prefix );

  /** Combines all defined content into `seed` via boost::hash_combine. */
  void appendToHash( std::size_t &seed ) const;

  /** Exact (bitwise) comparison of all members (availability-only for the
   parsed formula function).
   */
  bool operator==( const DetectorEfficiencyCurve &rhs ) const;

#if( PERFORM_DEVELOPER_CHECKS )
  /** Tests whether the two curves are equal to within small numerical
   rounding; throws std::runtime_error with a brief explanation when an
   issue is found.
   */
  static void equalEnough( const DetectorEfficiencyCurve &lhs,
                           const DetectorEfficiencyCurve &rhs );
#endif

private:
  DetectorPeakResponse::EfficiencyFnctForm m_form;

  /** PhysicalUnits energy unit the representation expects (keV == 1.0). */
  float m_energyUnits;

  /** Only filled for kEnergyEfficiencyPairs. */
  std::vector<DetectorPeakResponse::EnergyEfficiencyPair> m_energyEfficiencies;

  /** Only filled for kFunctialEfficienyForm. */
  std::string m_formula;

  /** Parsed formula; rebuilt on deserialization, not compared or hashed
   beyond availability.
   */
  std::function<float(float)> m_fcn;

  /** Only filled for kExpOfLogPowerSeries. */
  std::vector<float> m_expOfLogCoeffs;

  /** Legacy per-coefficient 1-sigma uncertainties for kExpOfLogPowerSeries;
   empty or same size as m_expOfLogCoeffs.
   */
  std::vector<float> m_expOfLogCoeffUncerts;

  /** Optional uncertainty of this curve. */
  std::shared_ptr<const DetectorEfficiencyUncert> m_uncert;
};//class DetectorEfficiencyCurve

#endif  //DetectorEfficiency_h
