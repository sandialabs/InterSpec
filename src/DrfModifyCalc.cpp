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

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <boost/functional/hash.hpp>

#include "io/DetectorResponse.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
/** The rows that describe a real point: a positive energy, and (unless the efficiency column is not
 in play) a positive efficiency.
 */
bool row_is_usable( const DrfModifyCalc::PointRow &row, const bool need_efficiency )
{
  if( (row.energy <= 0.0) || std::isnan(row.energy) || std::isinf(row.energy) )
    return false;

  if( !need_efficiency )
    return true;

  return (row.efficiency > 0.0) && !std::isnan(row.efficiency) && !std::isinf(row.efficiency);
}//row_is_usable(...)


/** The node covariance the rows' two uncertainty columns describe: the correlated column correlated
 across energy by `corrLength`, the independent column on the diagonal.
 */
shared_ptr<DetectorEfficiencyUncert> uncert_from_rows( const vector<DrfModifyCalc::PointRow> &rows,
                                                      const DrfModifyCalc::AnchorOptions &options,
                                                      const bool need_efficiency,
                                                      vector<DrfModifyCalc::Problem> &problems )
{
  vector<float> energies, corr_fracs, stat_fracs;
  bool any_uncert = false;
  for( const DrfModifyCalc::PointRow &row : rows )
  {
    if( !row_is_usable( row, need_efficiency ) )
      continue;

    // Covariance node energies are keV by contract, whatever units the curve (and hence the column)
    //  is written in - see DetectorEfficiencyUncert.
    energies.push_back( static_cast<float>( row.energy * options.energyUnits ) );
    corr_fracs.push_back( static_cast<float>( std::max( 0.0, row.fracCert ) ) );
    stat_fracs.push_back( static_cast<float>( std::max( 0.0, row.fracStat ) ) );
    any_uncert = (any_uncert || (corr_fracs.back() > 0.0f) || (stat_fracs.back() > 0.0f));
  }//for( const PointRow &row : rows )

  // Blank uncertainty columns mean the source states none - not a covariance of zeros, which is a
  //  claim of perfect knowledge and, worse, reads downstream as "this detector states an
  //  uncertainty" while contributing nothing.
  if( !any_uncert )
    energies.clear();

  // No rows with an uncertainty in them means the user is stating none.
  if( energies.empty() )
    return nullptr;

  try
  {
    // Passed through exactly as the correlation control gives it - including a value <= 0, which is
    //  how "Uncorrelated" is expressed (`EccUncertOptions::effectiveCorrLength` returns -1 for it)
    //  and which `fromCorrelatedPlusDiagonal` documents as making the correlated part diagonal.
    //  Substituting a default here would silently turn the user's "these errors are independent" into
    //  "these errors are one common mode" - the two extremes of the same control.
    return DetectorEfficiencyUncert::fromCorrelatedPlusDiagonal( energies, corr_fracs, stat_fracs,
                                                                 options.corrLength );
  }catch( std::exception &e )
  {
    problems.push_back( DrfModifyCalc::Problem( "dmw-err-uncert-not-built", e.what() ) );
    return nullptr;
  }
}//uncert_from_rows(...)


/** The `MeasuredEffPoint`s the rows describe, each starting from the point it was seeded with so an
 edit of one number does not erase the provenance (peak area, live time, file, distance uncertainty).
 */
vector<MeasuredEffPoint> measured_points_from_rows( const vector<DrfModifyCalc::PointRow> &rows,
                                    const shared_ptr<const MeasuredDrfPoints> &seedPoints,
                                    const DrfModifyCalc::AnchorOptions &options,
                                    const bool want_distance,
                                    const bool ref_distance_applies )
{
  const vector<MeasuredEffPoint> empty;
  const vector<MeasuredEffPoint> &seed = seedPoints ? seedPoints->points() : empty;

  vector<MeasuredEffPoint> answer;
  for( const DrfModifyCalc::PointRow &row : rows )
  {
    if( !row_is_usable( row, true ) )
      continue;

    MeasuredEffPoint point;
    if( (row.seedIndex >= 0) && (static_cast<size_t>(row.seedIndex) < seed.size()) )
      point = seed[row.seedIndex];

    // MeasuredEffPoint energies are keV by contract, whatever units the column is in.
    point.energy = static_cast<float>( row.energy * options.energyUnits );
    point.efficiency = static_cast<float>( row.efficiency );
    point.fracStatUncert = static_cast<float>( std::max( 0.0, row.fracStat ) );
    point.fracCertUncert = static_cast<float>( std::max( 0.0, row.fracCert ) );
    point.sourceKey = row.sourceKey;

    if( !want_distance )
    {
      point.distance = -1.0f;
      point.distanceUncert = 0.0f;
    }else if( row.distance > 0.0 )
    {
      point.distance = static_cast<float>( row.distance );
    }else if( ref_distance_applies && (options.refDistance > 0.0) )
    {
      // Only an absolute reference curve has a reference distance for a blank cell to mean; for a
      //  re-fit each point carries its own, and substituting another point's would silently move it.
      point.distance = static_cast<float>( options.refDistance );
    }

    answer.push_back( std::move(point) );
  }//for( const PointRow &row : rows )

  return answer;
}//measured_points_from_rows(...)


/** The geometric factor a point at `distance` sees: the ray-traced kernel when the detector states a
 geometry (right in the near field too), else the flat-disk solid angle.  Used to move a point
 measured at one distance onto a curve anchored at another.
 */
/** The ratio of geometric factors that moves a point measured at `from_dist` onto a curve anchored
 at `to_dist`; 0 when it cannot be computed.

 The kernel-or-disk choice is made ONCE for the pair - a ratio of a ray-traced factor to a flat-disk
 one is not a transfer of anything.
 */
double transfer_ratio( const DetectorPeakResponse &drf,
                       CeeLoUtils::GeometryKernel * const kernel,
                       const double energy_keV, const double from_dist, const double to_dist )
{
  if( kernel )
  {
    try
    {
      const double g_from = kernel->intrinsicFactor( energy_keV, from_dist / PhysicalUnits::cm );
      const double g_to = kernel->intrinsicFactor( energy_keV, to_dist / PhysicalUnits::cm );
      if( (g_from > 0.0) && (g_to > 0.0) )
        return g_to / g_from;
    }catch( std::exception & )
    {
      //fall through to the disk, for BOTH distances
    }
  }//if( kernel )

  const double diameter = drf.detectorDiameter();
  const double setback = drf.detectorSetback();
  const double g_from = DetectorPeakResponse::fractionalSolidAngle( diameter, from_dist + setback );
  const double g_to = DetectorPeakResponse::fractionalSolidAngle( diameter, to_dist + setback );

  return ((g_from > 0.0) && (g_to > 0.0)) ? (g_to / g_from) : 0.0;
}//transfer_ratio(...)


/** The kernel for `drf`, or null when it states no usable geometry. */
std::unique_ptr<CeeLoUtils::GeometryKernel> kernel_for_drf( const DetectorPeakResponse &drf )
{
  std::unique_ptr<CeeLoUtils::GeometryKernel> kernel;
  const shared_ptr<const ceelo::GeometryDescriptor> geom = drf.geometry();
  if( geom )
  {
    try
    {
      kernel.reset( new CeeLoUtils::GeometryKernel( *geom ) );
    }catch( std::exception & )
    {
    }
  }//if( geom )

  return kernel;
}//kernel_for_drf(...)


/** The source table that goes with `points`: the seed entries the points still reference, plus a stub
 for a source key the user typed that the table never had.  A table listing sources no point uses (or
 missing one a point names) is how a re-fit ends up using a certificate uncertainty that is not on
 screen.
 */
vector<MeasuredSourceInfo> sources_for_points( const vector<MeasuredEffPoint> &points,
                                    const shared_ptr<const MeasuredDrfPoints> &seedPoints,
                                    vector<DrfModifyCalc::Problem> &problems )
{
  auto referenced = [&points]( const string &key ) -> bool {
    for( const MeasuredEffPoint &point : points )
    {
      if( point.sourceKey == key )
        return true;
    }
    return false;
  };//referenced lambda

  // The sources the seed table still has points for, in the order it had them - so an edit that
  // changes no source leaves the table bit-identical.
  vector<MeasuredSourceInfo> answer;
  if( seedPoints )
  {
    for( const MeasuredSourceInfo &src : seedPoints->sources() )
    {
      if( referenced( src.sourceKey ) )
        answer.push_back( src );
    }
  }//if( seedPoints )

  // Then any key the user typed that the table never had: kept (the points' own certificate
  // uncertainty is what a fit uses), but reported - the assay information a re-fit would want for it
  // is not there.
  string unknown_keys;
  for( const MeasuredEffPoint &point : points )
  {
    if( point.sourceKey.empty() )
      continue;

    bool already_have = false;
    for( const MeasuredSourceInfo &src : answer )
      already_have = (already_have || (src.sourceKey == point.sourceKey));
    if( already_have )
      continue;

    MeasuredSourceInfo stub;
    stub.sourceKey = point.sourceKey;
    stub.fracActivityUncert = point.fracCertUncert;
    stub.distance = point.distance;
    answer.push_back( stub );

    unknown_keys += (unknown_keys.empty() ? "" : ", ") + point.sourceKey;
  }//for( const MeasuredEffPoint &point : points )

  if( !unknown_keys.empty() )
    problems.push_back( DrfModifyCalc::Problem( "dmw-warn-unknown-source", unknown_keys, false ) );

  return answer;
}//sources_for_points(...)

}//namespace


namespace DrfModifyCalc
{

AnchorEditor editorForDrf( const DetectorPeakResponse &drf )
{
  const shared_ptr<const DetectorEfficiencyCurve> curve = drf.efficiencyCurve();
  const shared_ptr<const MeasuredDrfPoints> points = drf.measuredPoints();
  const bool have_points = (points && !points->empty());

  if( curve && curve->isValid() )
  {
    switch( curve->form() )
    {
      case DetectorPeakResponse::kExpOfLogPowerSeries:
        // The points are what the equation was fit from, so they can be edited - but only through a
        //  re-fit, since they are absolute efficiencies and the equation is not.
        //
        //  Except when the equation ITSELF states absolute efficiency: a re-fit divides the source
        //  geometry out, so it would replace an absolute equation with an intrinsic one and leave
        //  the detector wrong by a solid angle (the very defect this editor exists to prevent).
        //  Such a DRF edits its equation directly, and keeps the points as provenance.
        if( have_points
            && (drf.geometryType() != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute) )
        {
          return AnchorEditor::RefitPoints;
        }
        return AnchorEditor::Coefficients;

      case DetectorPeakResponse::kFunctialEfficienyForm:
        return AnchorEditor::Formula;

      case DetectorPeakResponse::kEnergyEfficiencyPairs:
        // An absolute reference curve (ANGLE .outx) stores the same numbers twice - as the curve and
        //  as the raw points - so one edit writes both.  Any other pairs curve edits only the curve.
        if( have_points
            && (drf.geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute) )
        {
          return AnchorEditor::AbsolutePoints;
        }
        return AnchorEditor::CurvePairs;

      case DetectorPeakResponse::kNumEfficiencyFnctForms:
        break;
    }//switch( curve->form() )
  }//if( curve && curve->isValid() )

  return AnchorEditor::CurvePairs;
}//editorForDrf(...)


bool editorUsesMeasuredPoints( const AnchorEditor editor )
{
  return ((editor == AnchorEditor::RefitPoints) || (editor == AnchorEditor::AbsolutePoints));
}


bool editorUsesPointTable( const AnchorEditor editor )
{
  return (editor != AnchorEditor::Coefficients);
}


bool anyBlocking( const std::vector<Problem> &problems )
{
  for( const Problem &problem : problems )
  {
    if( problem.blocking )
      return true;
  }
  return false;
}//anyBlocking(...)


bool applyPointRows( DetectorPeakResponse &working,
                     const AnchorEditor editor,
                     const std::vector<PointRow> &rows,
                     const AnchorOptions &options,
                     const std::shared_ptr<const MeasuredDrfPoints> &seedPoints,
                     std::vector<Problem> &problems )
{
  // The three editors whose rows ARE the efficiency; a formula curve's rows are covariance nodes and
  //  go through #applyFormula, which knows the efficiency does not come from them.
  const bool is_point_editor = ((editor == AnchorEditor::RefitPoints)
                                || (editor == AnchorEditor::AbsolutePoints)
                                || (editor == AnchorEditor::CurvePairs));
  assert( is_point_editor );
  if( !is_point_editor )
  {
    problems.push_back( Problem( "dmw-err-wrong-editor" ) );
    return false;
  }

  // ---- the two editors whose rows are the raw measured points -------------
  if( editorUsesMeasuredPoints( editor ) )
  {
    const bool want_distance = !working.isFixedGeometry();
    const bool ref_distance_applies = (editor == AnchorEditor::AbsolutePoints);
    const vector<MeasuredEffPoint> points
        = measured_points_from_rows( rows, seedPoints, options, want_distance, ref_distance_applies );

    if( points.size() < 2 )
    {
      problems.push_back( Problem( "dmw-err-need-two-points" ) );
      return false;
    }

    if( want_distance )
    {
      for( const MeasuredEffPoint &point : points )
      {
        if( point.distance <= 0.0f )
        {
          problems.push_back( Problem( "dmw-err-need-distance" ) );
          return false;
        }
      }
    }//if( want_distance )

    // Collected separately so a note about the sources is not left standing next to a refusal that
    //  happens further down.
    vector<Problem> source_notes;
    auto measured = make_shared<MeasuredDrfPoints>();
    try
    {
      measured->setPoints( points );
      measured->setSources( sources_for_points( points, seedPoints, source_notes ) );
    }catch( std::exception &e )
    {
      problems.push_back( Problem( "dmw-err-points-invalid", e.what() ) );
      return false;
    }

    if( editor == AnchorEditor::RefitPoints )
    {
      // One fit produces the equation, its covariance, the node covariance and the points - so the
      //  refit source cannot stop describing the curve.
      string warnings;
      try
      {
        MakeDrfCalc::refitEfficiencyFromPoints( working, *measured, options.equationTerms, warnings );
      }catch( std::exception &e )
      {
        problems.push_back( Problem( "dmw-err-refit-failed", e.what() ) );
        return false;
      }

      if( !warnings.empty() )
        problems.push_back( Problem( "dmw-note-refit-warning", warnings, false ) );

      problems.insert( end(problems), begin(source_notes), end(source_notes) );
      return true;
    }//if( editor == AnchorEditor::RefitPoints )

    // AbsolutePoints: the curve IS these numbers, at the reference distance.
    if( options.refDistance <= 0.0 )
    {
      problems.push_back( Problem( "dmw-err-need-ref-distance" ) );
      return false;
    }

    const float diameter = working.detectorDiameter();
    if( diameter <= 0.0f )
    {
      problems.push_back( Problem( "dmw-err-no-diameter" ) );
      return false;
    }

    // The curve states absolute efficiency AT the reference distance, so a point measured somewhere
    //  else is moved there by the ratio of geometric factors first.  Without this a point taken at
    //  50 cm on a curve anchored at 25 cm is written in about four times too small - the same class
    //  of mistake as reading an absolute efficiency as an intrinsic one.
    const std::unique_ptr<CeeLoUtils::GeometryKernel> kernel = kernel_for_drf( working );

    vector<DetectorPeakResponse::EnergyEffPoint> effpts;
    for( const MeasuredEffPoint &point : points )
    {
      double eff_at_ref = point.efficiency;
      if( (point.distance > 0.0f)
          && (fabs(point.distance - options.refDistance) > (1.0E-6*options.refDistance)) )
      {
        const double ratio = transfer_ratio( working, kernel.get(), point.energy,
                                             point.distance, options.refDistance );
        if( ratio <= 0.0 )
        {
          problems.push_back( Problem( "dmw-err-no-distance-transfer",
                                       std::to_string( static_cast<int>(point.energy) ) ) );
          return false;
        }

        eff_at_ref = point.efficiency * ratio;
      }//if( this point was measured somewhere other than the reference distance )

      DetectorPeakResponse::EnergyEffPoint e;
      e.energy = point.energy;
      e.efficiency = static_cast<float>( eff_at_ref );
      const double combo = std::sqrt( point.fracStatUncert*point.fracStatUncert
                                      + point.fracCertUncert*point.fracCertUncert );
      if( combo > 0.0 )
        e.efficiencyUncert = static_cast<float>( eff_at_ref * combo );
      effpts.push_back( e );
    }//for( const MeasuredEffPoint &point : points )

    // setEfficiencyPoints is a re-characterization: it resets the setback, the flags, the air
    //  attenuation choice, the energy range and the creation time, none of which this edit touched.
    //  Carry them across rather than letting one edited cell silently change every absolute
    //  efficiency the detector reports.
    const double setback = working.detectorSetback();
    const bool air_atten = working.absEffCorrectForAirAtten();
    try
    {
      working.setEfficiencyPoints( effpts, diameter, options.refDistance,
                                   DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    }catch( std::exception &e )
    {
      problems.push_back( Problem( "dmw-err-points-invalid", e.what() ) );
      return false;
    }

    if( setback > 0.0 )
      working.setDetectorSetback( setback );
    working.setAbsEffCorrectForAirAtten( air_atten );

    working.setMeasuredPoints( measured );

    // The uncertainty the edited stat/cert/source structure implies:
    //  C[i][j] = delta_ij*stat_i^2 + cert_i*cert_j*[same source].
    const shared_ptr<DetectorEfficiencyUncert> uncert = measured->toEfficiencyUncert();
    if( uncert && !uncert->isEmpty() )
      working.setEfficiencyUncert( uncert );

    problems.insert( end(problems), begin(source_notes), end(source_notes) );
    return true;
  }//if( editorUsesMeasuredPoints( editor ) )

  // ---- CurvePairs: the rows are the curve's own numbers -------------------
  vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
  for( const PointRow &row : rows )
  {
    if( !row_is_usable( row, true ) )
      continue;

    DetectorPeakResponse::EnergyEfficiencyPair pair;
    pair.energy = static_cast<float>( row.energy );
    pair.efficiency = static_cast<float>( row.efficiency );
    pairs.push_back( pair );
  }//for( const PointRow &row : rows )

  if( pairs.size() < 2 )
  {
    problems.push_back( Problem( "dmw-err-need-two-points" ) );
    return false;
  }

  try
  {
    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromPairs( pairs, options.energyUnits );
    curve->setUncertainty( uncert_from_rows( rows, options, true, problems ) );
    if( anyBlocking( problems ) )
      return false;

    working.replaceEfficiencyCurve( curve );
  }catch( std::exception &e )
  {
    problems.push_back( Problem( "dmw-err-points-invalid", e.what() ) );
    return false;
  }

  // The raw points (if any) were the source of the curve that was just replaced, so they no longer
  //  describe it.  Say so rather than leaving a refit source that would quietly undo the edit.
  if( working.measuredPoints() && !working.measuredPoints()->empty() )
    problems.push_back( Problem( "dmw-note-points-stale", string(), false ) );

  return true;
}//applyPointRows(...)


bool applyCoefficients( DetectorPeakResponse &working,
                        const std::vector<float> &coefficients,
                        const std::vector<double> &sigmas,
                        const std::vector<double> &rho,
                        const float energyUnits,
                        const bool writeCovariance,
                        std::vector<Problem> &problems )
{
  if( coefficients.empty() )
  {
    problems.push_back( Problem( "dmw-err-no-coefs" ) );
    return false;
  }

  const size_t n = coefficients.size();

  vector<float> coefCov;
  bool clear_covariance = false;   //the user zeroed it all: "I do not know", not "keep the old one"
  if( writeCovariance )
  {
    const vector<double> cov = covarianceFromSigmaRho( sigmas, rho );
    if( cov.size() != (n*n) )
    {
      problems.push_back( Problem( "dmw-err-cov-size" ) );
      return false;
    }

    std::string why;
    if( !DetectorEfficiencyUncert::covarianceIsUsable( cov, &why ) )
    {
      problems.push_back( Problem( "dmw-err-cov-not-psd" ) );
      return false;
    }

    bool any_nonzero = false;
    coefCov.resize( n*n );
    for( size_t i = 0; i < (n*n); ++i )
    {
      coefCov[i] = static_cast<float>( cov[i] );
      any_nonzero = (any_nonzero || (coefCov[i] != 0.0f));
    }

    if( !any_nonzero )
    {
      coefCov.clear();  //all-zero is "no covariance given", not a covariance of zero
      clear_covariance = true;
    }
  }//if( writeCovariance )

  // The equation's diagonal doubles as the legacy per-coefficient uncertainty, which is what the
  //  older DB fields and (when the covariance has to be dropped to fit) the app-URL carry.  When no
  //  covariance is being written, whatever the curve already stated is kept - erasing it would throw
  //  away the only uncertainty some DRFs have, and it is what the placeholder matrix is seeded from.
  vector<float> uncerts;
  if( !coefCov.empty() )
  {
    uncerts.resize( n );
    for( size_t i = 0; i < n; ++i )
      uncerts[i] = std::sqrt( std::max( 0.0f, coefCov[i*n + i] ) );
  }else if( !clear_covariance )
  {
    const shared_ptr<const DetectorEfficiencyCurve> existing_curve = working.efficiencyCurve();
    const vector<float> &legacy = existing_curve ? existing_curve->expOfLogPowerSeriesUncerts()
                                                 : vector<float>{};
    if( legacy.size() == n )
      uncerts = legacy;
  }

  // Whether the equation itself moved, as opposed to only its uncertainty.
  bool curve_changed = true;
  {
    const shared_ptr<const DetectorEfficiencyCurve> existing_curve = working.efficiencyCurve();
    if( existing_curve && (existing_curve->form() == DetectorPeakResponse::kExpOfLogPowerSeries)
        && (existing_curve->energyUnits() == energyUnits)
        && (existing_curve->expOfLogPowerSeriesCoeffs() == coefficients) )
    {
      curve_changed = false;
    }
  }

  try
  {
    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromExpOfLogPowerSeries( coefficients, uncerts, energyUnits );

    if( !coefCov.empty() )
    {
      // For an equation curve this IS the uncertainty the fits see (DetectorEfficiencyCurve::
      //  fracCovariance propagates it as J*Sigma*J^T).  Any node covariance is left untouched
      //  beside it as provenance; the two are never combined.
      const shared_ptr<const DetectorEfficiencyUncert> existing = working.efficiencyUncert();
      auto uncert = existing ? make_shared<DetectorEfficiencyUncert>( *existing )
                             : make_shared<DetectorEfficiencyUncert>();
      uncert->setCoefficientCovariance( coefCov );
      curve->setUncertainty( uncert->isEmpty() ? nullptr : uncert );
    }else
    {
      // Carrying the existing uncertainty across - except for a coefficient covariance that no
      //  longer applies.  Two ways that happens: the user zeroed every sigma (they are saying they
      //  do not know it, so keeping the old matrix would be putting words in their mouth), or the
      //  term count changed and its rank no longer matches, in which case
      //  `DetectorEfficiencyCurve::fracCovariance` would silently ignore it while the detector still
      //  looked like it stated one.  Either way it goes, and the user is told.
      const shared_ptr<const DetectorEfficiencyUncert> existing = working.efficiencyUncert();
      const bool rank_mismatch = (existing && !existing->coefficientCovariance().empty()
                                  && (existing->coefficientCovariance().size() != (n*n)));

      if( existing && (clear_covariance || rank_mismatch) )
      {
        auto trimmed = make_shared<DetectorEfficiencyUncert>( *existing );
        trimmed->setCoefficientCovariance( {} );
        curve->setUncertainty( trimmed->isEmpty() ? nullptr : trimmed );
        if( rank_mismatch )
          problems.push_back( Problem( "dmw-note-cov-dropped", string(), false ) );
      }else
      {
        curve->setUncertainty( existing );
      }
    }

    working.replaceEfficiencyCurve( curve );
  }catch( std::exception &e )
  {
    problems.push_back( Problem( "dmw-err-coefs-invalid", e.what() ) );
    return false;
  }

  // An equation edited by hand is no longer the fit of the stored points - but only if it actually
  //  changed.  Saying so for an uncertainty-only edit is both false and, because the apply-time
  //  invariant check treats this note as "the user meant to", quietly disarming.
  if( working.measuredPoints() && !working.measuredPoints()->empty() && curve_changed )
    problems.push_back( Problem( "dmw-note-points-stale", string(), false ) );

  return true;
}//applyCoefficients(...)


bool applyFormula( DetectorPeakResponse &working,
                   const std::string &formula,
                   const float energyUnits,
                   const std::vector<PointRow> &rows,
                   const AnchorOptions &options,
                   std::vector<Problem> &problems )
{
  if( formula.empty() )
  {
    problems.push_back( Problem( "dmw-err-no-formula" ) );
    return false;
  }

  try
  {
    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromFormula( formula, energyUnits );

    // The rows here are covariance nodes - the formula supplies the efficiency - so they are read
    //  without an efficiency column.
    const shared_ptr<DetectorEfficiencyUncert> uncert
        = uncert_from_rows( rows, options, false, problems );
    if( anyBlocking( problems ) )
      return false;

    // With rows on screen, they ARE the uncertainty: emptying their columns says "none", the same as
    //  it does for a pairs curve.  Only a formula with no row table at all keeps what it had.
    // With a row table on screen the rows ARE the uncertainty, so emptying them says "none" - the
    //  same as it does for a pairs curve.  Only a caller with no row table at all (not the dialog,
    //  which always shows one for a formula) keeps whatever the curve had.
    curve->setUncertainty( options.hasRowTable ? uncert : working.efficiencyUncert() );
    working.replaceEfficiencyCurve( curve );
  }catch( std::exception &e )
  {
    problems.push_back( Problem( "dmw-err-formula-invalid", e.what() ) );
    return false;
  }

  // The stored points were the source of the curve this just replaced.
  if( working.measuredPoints() && !working.measuredPoints()->empty() )
    problems.push_back( Problem( "dmw-note-points-stale", string(), false ) );

  return true;
}//applyFormula(...)


std::vector<double> covarianceFromSigmaRho( const std::vector<double> &sigmas,
                                            const std::vector<double> &rho )
{
  const size_t n = sigmas.size();
  if( !n || (rho.size() != (n*n)) )
    return {};

  vector<double> cov( n*n, 0.0 );
  for( size_t i = 0; i < n; ++i )
  {
    cov[i*n + i] = sigmas[i] * sigmas[i];
    for( size_t j = i + 1; j < n; ++j )
    {
      const double r = 0.5*(rho[i*n + j] + rho[j*n + i]);
      cov[i*n + j] = cov[j*n + i] = r * sigmas[i] * sigmas[j];
    }
  }//for( size_t i = 0; i < n; ++i )

  return cov;
}//covarianceFromSigmaRho(...)


void sigmaRhoFromCovariance( const std::vector<float> &covRowMajor,
                             std::vector<double> &sigmas,
                             std::vector<double> &rho )
{
  sigmas.clear();
  rho.clear();

  if( covRowMajor.empty() )
    return;

  const size_t n = static_cast<size_t>( std::lround( std::sqrt( static_cast<double>(covRowMajor.size()) ) ) );
  if( (n*n) != covRowMajor.size() )
    return;

  sigmas.resize( n, 0.0 );
  for( size_t i = 0; i < n; ++i )
  {
    const double var = covRowMajor[i*n + i];
    sigmas[i] = (var > 0.0) ? std::sqrt(var) : 0.0;
  }

  rho.assign( n*n, 0.0 );
  for( size_t i = 0; i < n; ++i )
  {
    rho[i*n + i] = 1.0;
    for( size_t j = i + 1; j < n; ++j )
    {
      double r = 0.0;
      if( (sigmas[i] > 0.0) && (sigmas[j] > 0.0) )
        r = covRowMajor[i*n + j] / (sigmas[i] * sigmas[j]);
      r = std::max( -1.0, std::min( 1.0, r ) );
      rho[i*n + j] = rho[j*n + i] = r;
    }
  }//for( size_t i = 0; i < n; ++i )
}//sigmaRhoFromCovariance(...)


bool sigmaRhoFromLegacyUncerts( const DetectorPeakResponse &drf,
                                std::vector<double> &sigmas,
                                std::vector<double> &rho )
{
  sigmas.clear();
  rho.clear();

  const shared_ptr<const DetectorEfficiencyCurve> curve = drf.efficiencyCurve();
  if( !curve )
    return false;

  const vector<float> &legacy = curve->expOfLogPowerSeriesUncerts();
  const size_t n = curve->expOfLogPowerSeriesCoeffs().size();
  if( !n || (legacy.size() != n) )
    return false;

  bool any = false;
  sigmas.assign( n, 0.0 );
  for( size_t i = 0; i < n; ++i )
  {
    sigmas[i] = std::max( 0.0f, legacy[i] );
    any = (any || (sigmas[i] > 0.0));
  }

  if( !any )
  {
    sigmas.clear();
    return false;
  }

  rho.assign( n*n, 0.0 );
  for( size_t i = 0; i < n; ++i )
    rho[i*n + i] = 1.0;

  return true;
}//sigmaRhoFromLegacyUncerts(...)


std::size_t seedFingerprint( const DetectorPeakResponse &drf )
{
  std::size_t seed = 0;

  const shared_ptr<const DetectorEfficiencyCurve> curve = drf.efficiencyCurve();
  if( curve )
    curve->appendToHash( seed );  //covers the representation AND its uncertainty, both covariances

  const shared_ptr<const MeasuredDrfPoints> points = drf.measuredPoints();
  if( points && !points->empty() )
    points->appendToHash( seed );

  // The geometry a response is ray-traced for.  Read the member-or-response geometry, so that
  //  attaching a response generated for the same geometry does not by itself look like a change.
  const shared_ptr<const ceelo::GeometryDescriptor> geom = drf.geometry();
  if( geom )
    boost::hash_combine( seed, geom->to_xml_string() );

  // Things a transfer or grounding anchor depends on beyond the curve itself.  The energy range and
  //  the air-attenuation flag are in here because the anchor builders read them: the range picks the
  //  energies a curve-derived anchor is sampled at (CeeLoUtils::transferAnchorForDrf,
  //  curveAnchorWithCovarianceForDrf, MakeMcResponseForDrf::groundingPointsForDrf), and the flag
  //  changes what `farFieldIntrinsicEfficiency` returns for an absolute curve, which is what they sample.
  boost::hash_combine( seed, drf.detectorDiameter() );
  boost::hash_combine( seed, drf.detectorSetback() );
  boost::hash_combine( seed, static_cast<int>(drf.geometryType()) );
  boost::hash_combine( seed, drf.absoluteEfficiencyDistance() );
  boost::hash_combine( seed, drf.absEffCorrectForAirAtten() );
  boost::hash_combine( seed, drf.lowerEnergy() );
  boost::hash_combine( seed, drf.upperEnergy() );

  return seed;
}//seedFingerprint(...)


bool checkDrfSelfConsistent( const DetectorPeakResponse &drf, std::string &why )
{
  why.clear();

  auto fail = [&why]( const string &msg ) -> bool {
    why = msg;
    return false;
  };//fail lambda

  const shared_ptr<const DetectorEfficiencyCurve> curve = drf.efficiencyCurve();
  if( !curve || !curve->isValid() )
    return fail( "the efficiency curve is not valid" );

  if( curve->energyUnits() <= 0.0f )
    return fail( "the efficiency curve has non-positive energy units" );

  const shared_ptr<const DetectorEfficiencyUncert> uncert = drf.efficiencyUncert();
  if( uncert )
  {
    const vector<float> &coef_cov = uncert->coefficientCovariance();
    if( !coef_cov.empty() )
    {
      const vector<double> as_dbl( begin(coef_cov), end(coef_cov) );
      std::string cov_why;
      if( !DetectorEfficiencyUncert::covarianceIsUsable( as_dbl, &cov_why ) )
        return fail( "the stored coefficient covariance is not usable: " + cov_why );

      // A rank mismatch is not an error the DRF can detect at read time, but it means the covariance
      //  is silently ignored by every query - so an edit that produced one did nothing.
      const size_t ncoef = curve->expOfLogPowerSeriesCoeffs().size();
      if( curve->form() == DetectorPeakResponse::kExpOfLogPowerSeries )
      {
        if( coef_cov.size() != (ncoef*ncoef) )
          return fail( "the coefficient covariance does not match the number of equation terms,"
                       " so nothing will use it" );
      }
    }//if( !coef_cov.empty() )

    if( uncert->hasNodeCovariance() )
    {
      const vector<float> &node_cov = uncert->covarianceMatrix();
      const vector<double> as_dbl( begin(node_cov), end(node_cov) );
      std::string cov_why;
      if( !DetectorEfficiencyUncert::covarianceIsUsable( as_dbl, &cov_why ) )
        return fail( "the stored node covariance is not usable: " + cov_why );
    }
  }//if( uncert )

  const shared_ptr<const MeasuredDrfPoints> points = drf.measuredPoints();
  if( points && !points->empty() )
  {
    const bool far_field = !drf.isFixedGeometry();
    for( const MeasuredEffPoint &point : points->points() )
    {
      if( point.energy <= 0.0f )
        return fail( "a measured point has a non-positive energy" );
      if( far_field && (point.efficiency > 0.0f) && (point.distance <= 0.0f) )
        return fail( "a measured point of a far-field detector has no source distance, so the"
                     " efficiency it states cannot be interpreted" );
      if( !far_field && (point.distance > 0.0f) )
        return fail( "a measured point of a fixed-geometry detector states a source distance, which"
                     " that geometry type has no meaning for" );
    }

    // An absolute curve states its efficiency AT a distance; without one, nothing below can check
    //  the points against it, and the curve cannot be interpreted either.
    if( (drf.geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute)
        && (drf.absoluteEfficiencyDistance() <= 0.0) )
    {
      return fail( "an absolute-efficiency curve has no distance to state it at" );
    }

    // Invariant 4: the raw points are what the curve was made from, so the curve has to still say
    //  what they say.  Checked loosely - a fit does not pass through its points - by comparing the
    //  curve against the efficiencies the points imply in the curve's own convention, so a merely
    //  poor fit passes and a convention mix-up does not.
    //
    //  This is the check that catches the defect the whole editing rework is about: absolute
    //  efficiencies (what a measured point is) installed as the intrinsic curve differ by one solid
    //  angle - a factor of about 200 for a 7 cm crystal at 25 cm.
    //
    //  A FarFieldAbsolute curve states absolute efficiency at one distance, so for it the comparison
    //  is against the points measured at that distance (there is no geometry to divide out).
    const bool abs_curve = (drf.geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute);
    const bool in_mev = (curve->energyUnits() > 10.0f);

    // Per-point |ln(curve/points)|, judged by the MEDIAN rather than the worst.  What this invariant
    //  is for is a convention error - absolute read as intrinsic, the wrong distance, the wrong
    //  units - which moves EVERY point by a similar factor.  One point far off is a mis-typed cell:
    //  a real thing to tell the user about (the apply does), but a user error, not a broken program,
    //  and a developer-check assert on it would abort the build the moment somebody fat-fingers a
    //  number.
    vector<double> residuals;
    if( abs_curve )
    {
      // The curve is absolute efficiency at one distance; a point measured elsewhere is moved there
      //  the same way the apply does, so "this point is at the wrong geometry" is caught rather than
      //  skipped (skipping it was exactly the hole that let a mis-scaled point through).
      const double ref_dist = drf.absoluteEfficiencyDistance();
      const std::unique_ptr<CeeLoUtils::GeometryKernel> kernel = kernel_for_drf( drf );

      for( const MeasuredEffPoint &point : points->points() )
      {
        if( (point.efficiency <= 0.0f) || (point.energy <= 0.0f) || (ref_dist <= 0.0) )
          continue;

        double point_at_ref = point.efficiency;
        if( (point.distance > 0.0f) && (fabs(point.distance - ref_dist) > (1.0E-3*ref_dist)) )
        {
          const double ratio = transfer_ratio( drf, kernel.get(), point.energy,
                                               point.distance, ref_dist );
          if( ratio <= 0.0 )
            continue;
          point_at_ref = point.efficiency * ratio;
        }//if( measured somewhere other than the reference distance )

        const double curve_eff = curve->efficiency( point.energy );
        if( curve_eff <= 0.0 )
          continue;

        residuals.push_back( fabs( log( curve_eff / point_at_ref ) ) );
      }//for( const MeasuredEffPoint &point : points->points() )
    }else
    {
      const MakeDrfCalc::GeometryChoice geom = MakeDrfCalc::geometryChoiceForDrf( drf );

      vector<MakeDrfFit::EffFitPoint> fitpts;
      try
      {
        fitpts = MakeDrfCalc::intrinsicFitPoints( *points, geom, in_mev );
      }catch( std::exception & )
      {
        fitpts.clear();  //no geometry to divide out; nothing to check against
      }

      for( const MakeDrfFit::EffFitPoint &fitpt : fitpts )
      {
        if( (fitpt.efficiency <= 0.0f) || (fitpt.energy <= 0.0f) )
          continue;

        const float energy_kev = in_mev ? (1000.0f*fitpt.energy) : fitpt.energy;
        const double curve_eff = curve->efficiency( energy_kev );
        if( curve_eff <= 0.0 )
          continue;

        const double ratio = curve_eff / static_cast<double>( fitpt.efficiency );
        residuals.push_back( fabs( log( ratio ) ) );
      }//for( const MakeDrfFit::EffFitPoint &fitpt : fitpts )
    }//if( abs_curve ) / else

    // A factor of two at the median is far beyond any fit quality.
    if( !residuals.empty() )
    {
      std::sort( begin(residuals), end(residuals) );
      const double median = residuals[ residuals.size()/2 ];

      if( median > 0.7 )
      {
        char buffer[128];
        snprintf( buffer, sizeof(buffer), "%.3g", exp(median) );
        return fail( "the efficiency curve disagrees with the measured points it is supposed to come"
                     " from, by a factor of " + string(buffer)
                     + " at half its points - one of the two is not what it claims to be" );
      }
    }//if( !residuals.empty() )
  }//if( points )

  return true;
}//checkDrfSelfConsistent(...)


UncertSummary uncertSummary( const DetectorPeakResponse &drf, const float energy )
{
  UncertSummary answer;

  if( !drf.isValid() )
    return answer;

  float at_energy = energy;
  if( at_energy <= 0.0f )
  {
    // Somewhere representative: Cs-137 when the DRF covers it, else the middle of its range.
    at_energy = 661.7f;
    const double lower = drf.lowerEnergy(), upper = drf.upperEnergy();
    if( (upper > lower) && (lower > 0.0) && ((at_energy < lower) || (at_energy > upper)) )
      at_energy = static_cast<float>( 0.5*(lower + upper) );
  }//if( at_energy <= 0.0f )

  vector<double> model_cov;
  const vector<double> cov = drf.efficiencyFracCovariance( { static_cast<double>(at_energy) },
                                                           &model_cov );
  if( cov.size() != 1 )
    return answer;

  answer.valid = true;
  answer.energy = at_energy;
  answer.total = std::sqrt( std::max( 0.0, cov[0] ) );
  answer.model = (model_cov.size() == 1) ? std::sqrt( std::max( 0.0, model_cov[0] ) ) : 0.0;
  answer.model = std::min( answer.model, answer.total );
  answer.data = std::sqrt( std::max( 0.0, answer.total*answer.total - answer.model*answer.model ) );

  // Whether anything the DRF itself states backs up the data-derived part - see the field's comment.
  //  One implementation, shared with DetectorPeakResponse::toJSON, so the chart cannot call a number
  //  "from data" that the chip calls assumed.
  answer.dataIsAssumed = (answer.data > 0.0) && !drf.statesOwnEfficiencyUncert();

  return answer;
}//uncertSummary(...)

}//namespace DrfModifyCalc
