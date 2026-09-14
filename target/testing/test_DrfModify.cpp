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

/** The apply paths of the "Modify Detector Response" tool (`DrfModifyCalc`, which
 `DrfModifyWidget` is the Wt plumbing around).

 These exist because every defect they cover is of the form "the DRF that comes out has pieces that
 disagree with each other", which no click-through can see: the tool's own chart, and the attached
 CeeLo response, hid the worst of them until someone exported the detector or switched it back to the
 flat-disk model.

 The one to read first is `uncertainty_only_edit_keeps_the_efficiency`.  A Create-DRF detector's
 measured points are ABSOLUTE efficiencies at their own source distances, while its curve is
 INTRINSIC; the apply used to install those rows as the curve, which made the detector about 200x
 too insensitive for a 7 cm crystal at 25 cm - triggered by editing any cell, including an
 uncertainty.  `inconsistent_drf_is_detected` builds that exact wrong DRF and shows the invariant
 check catching it.
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE TestDrfModify
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "io/DetectorResponse.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  const int argc = boost::unit_test::framework::master_test_suite().argc;
  char ** const argv = boost::unit_test::framework::master_test_suite().argv;

  string datadir;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }
  SpecUtils::ireplace_all( datadir, "%20", " " );

  if( datadir.empty() )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.reactiongamma.xml") ) )
      {
        datadir = d;
        break;
      }
    }
  }//if( datadir.empty() )

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(datadir, "sandia.decay.xml") ),
                         "sandia.decay.xml not in '" << datadir << "'; pass --datadir=" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path( datadir, "sandia.decay.xml" ) );
}//set_data_dir()


/** The 3"x3" NaI in a 0.5 mm Al can, as test_MakeDrfEndToEnd and test_ShieldingSourceFitCalc use. */
ceelo::GeometryDescriptor synthetic_nai_descriptor()
{
  ceelo::GeometryDescriptor gd;
  gd.shape = ceelo::DetectorShape::Cylinder;
  gd.dimensions_cm = { 3.81, 7.62 };
  gd.materials = { ceelo::MaterialSpec::from( ceelo::make_NaI() ),
                   ceelo::MaterialSpec::from( ceelo::make_Aluminum() ) };
  gd.crystal_material_index = 0;
  ceelo::LayerSpec can;
  can.material_index = 1;
  can.front_thickness_cm = 0.05;
  can.side_thickness_cm = 0.05;
  can.z_end_cm = 7.62;
  gd.layers.push_back( can );
  gd.reference_point = ceelo::ReferencePoint::EndcapFront;
  return gd;
}//synthetic_nai_descriptor()


const vector<double> sm_energies = { 59.5, 81.0, 121.8, 276.4, 302.9, 356.0, 383.8,
                                     661.7, 1173.2, 1332.5, 1408.0 };

/** A plausible intrinsic full-energy efficiency: exp of a power series in ln(E/MeV), so the fit form
 can reproduce it and the test is about the editing, not about fit quality. */
double intrinsic_truth( const double energy_keV )
{
  const double lx = log( energy_keV / 1000.0 );
  return exp( -0.45 - 0.62*lx - 0.20*lx*lx - 0.02*lx*lx*lx );
}


MakeDrfCalc::GeometryChoice nai_geometry()
{
  MakeDrfCalc::GeometryChoice geom;
  geom.geometry = make_shared<const ceelo::GeometryDescriptor>( synthetic_nai_descriptor() );
  geom.diameter = 2.0 * geom.geometry->transverse_half_extent() * PhysicalUnits::cm;
  geom.setback = geom.geometry->endcap_front_offset_cm() * PhysicalUnits::cm;
  return geom;
}


/** The measured points the Create DRF tool records: ABSOLUTE efficiency at each source's own
 distance, with a per-source certificate uncertainty and the peak provenance. */
MeasuredDrfPoints make_measured_points()
{
  const MakeDrfCalc::GeometryChoice geom = nai_geometry();
  CeeLoUtils::GeometryKernel kernel( *geom.geometry );

  struct Src { const char *name; vector<double> energies; double d_cm; };
  const vector<Src> sources = {
    { "Am241/SRS-1", { 59.5 },                             25.0 },
    { "Ba133/SRS-2", { 81.0, 276.4, 302.9, 356.0, 383.8 }, 50.0 },
    { "Cs137/SRS-3", { 661.7 },                            100.0 },
    { "Co60/SRS-4",  { 1173.2, 1332.5 },                   25.0 },
    { "Eu152/SRS-5", { 121.8, 1408.0 },                    100.0 },
  };

  vector<MeasuredEffPoint> points;
  vector<MeasuredSourceInfo> srcs;
  for( const Src &src : sources )
  {
    for( const double energy : src.energies )
    {
      MeasuredEffPoint p;
      p.energy = static_cast<float>( energy );
      p.efficiency = static_cast<float>( intrinsic_truth(energy)
                                         * kernel.intrinsicFactor( energy, src.d_cm ) );
      p.fracStatUncert = 0.004f;
      p.fracCertUncert = 0.03f;
      p.sourceKey = src.name;
      p.distance = static_cast<float>( src.d_cm * PhysicalUnits::cm );
      p.distanceUncert = static_cast<float>( 0.5 * PhysicalUnits::cm );
      p.peakArea = 12345.0f;
      p.peakAreaUncert = 111.0f;
      p.liveTime = 600.0f;
      p.fileName = "characterization.n42";
      p.sampleNumbers = "1";
      points.push_back( p );
    }//for( const double energy : src.energies )

    MeasuredSourceInfo info;
    info.sourceKey = src.name;
    info.nuclide = string(src.name).substr( 0, string(src.name).find('/') );
    info.activity = 10.0 * PhysicalUnits::microCi;
    info.fracActivityUncert = 0.03f;
    info.distance = static_cast<float>( src.d_cm * PhysicalUnits::cm );
    info.assayInfo = "certificate 2026-01-01";
    srcs.push_back( info );
  }//for( const Src &src : sources )

  MeasuredDrfPoints answer;
  answer.setPoints( points );
  answer.setSources( srcs );
  return answer;
}//make_measured_points()


/** A DRF exactly as the "Create Detector Response Function" tool makes one: a fitted equation with
 its coefficient covariance, the raw points beside it, a geometry, and (since the geometry is stated)
 an attached curve-transfer response. */
shared_ptr<DetectorPeakResponse> make_created_drf( const MeasuredDrfPoints &points )
{
  const MakeDrfCalc::GeometryChoice geom = nai_geometry();
  const vector<MakeDrfFit::EffFitPoint> effpts = MakeDrfCalc::intrinsicFitPoints( points, geom, true );
  BOOST_REQUIRE_EQUAL( effpts.size(), points.points().size() );

  MakeDrfCalc::FitResults fit;
  fit.eff = MakeDrfFit::performEfficiencyFit( effpts, 4 );
  fit.effInMeV = true;
  fit.lowerEnergy = static_cast<float>( sm_energies.front() );
  fit.upperEnergy = static_cast<float>( sm_energies.back() );

  shared_ptr<DetectorPeakResponse> drf
      = MakeDrfCalc::assembleDrf( "created", "unit test", points, geom, fit, nullptr );
  BOOST_REQUIRE( drf && drf->isValid() );
  return drf;
}//make_created_drf(...)


/** The point table as the dialog seeds it: one row per measured point, remembering which point it
 came from. */
vector<DrfModifyCalc::PointRow> rows_from_points( const MeasuredDrfPoints &points )
{
  vector<DrfModifyCalc::PointRow> rows;
  const vector<MeasuredEffPoint> &pts = points.points();
  for( size_t i = 0; i < pts.size(); ++i )
  {
    DrfModifyCalc::PointRow row;
    row.rowNumber = static_cast<int>(i) + 1;
    row.seedIndex = static_cast<int>(i);
    row.energy = pts[i].energy;              //column is keV for measured points
    row.efficiency = pts[i].efficiency;
    row.fracStat = pts[i].fracStatUncert;
    row.fracCert = pts[i].fracCertUncert;
    row.distance = pts[i].distance;
    row.sourceKey = pts[i].sourceKey;
    rows.push_back( row );
  }//for( size_t i = 0; i < pts.size(); ++i )

  return rows;
}//rows_from_points(...)


DrfModifyCalc::AnchorOptions measured_options()
{
  DrfModifyCalc::AnchorOptions options;
  options.energyUnits = static_cast<float>( PhysicalUnits::keV );
  options.equationTerms = 0;   //keep whatever the curve has
  return options;
}


string problems_to_string( const vector<DrfModifyCalc::Problem> &problems )
{
  string answer;
  for( const DrfModifyCalc::Problem &p : problems )
    answer += (answer.empty() ? "" : "; ") + p.messageId + (p.blocking ? "" : " (note)");
  return answer;
}


/** Largest fractional difference between the two DRFs' intrinsic efficiencies over the
 characterization energies. */
double worst_intrinsic_diff( const DetectorPeakResponse &lhs, const DetectorPeakResponse &rhs )
{
  double worst = 0.0;
  for( const double energy : sm_energies )
  {
    const double a = lhs.intrinsicEfficiency( static_cast<float>(energy) );
    const double b = rhs.intrinsicEfficiency( static_cast<float>(energy) );
    if( (a > 0.0) && (b > 0.0) )
      worst = std::max( worst, fabs( a/b - 1.0 ) );
  }
  return worst;
}//worst_intrinsic_diff(...)
}//namespace


BOOST_AUTO_TEST_CASE( editor_follows_what_the_drf_carries )
{
  set_data_dir();

  // A Create-DRF detector: an equation that carries the points it was fit from.
  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> created = make_created_drf( points );
  BOOST_CHECK( DrfModifyCalc::editorForDrf( *created ) == DrfModifyCalc::AnchorEditor::RefitPoints );
  BOOST_CHECK( DrfModifyCalc::editorUsesMeasuredPoints( DrfModifyCalc::AnchorEditor::RefitPoints ) );

  // Which editor is shown must not depend on anything the Flat-Disk / Geometry-Modeled toggle
  //  changes; detaching the response is exactly that change.
  shared_ptr<DetectorPeakResponse> detached = make_shared<DetectorPeakResponse>( *created );
  detached->setCeeloResponse( nullptr );
  BOOST_CHECK( DrfModifyCalc::editorForDrf( *detached ) == DrfModifyCalc::AnchorEditor::RefitPoints );

  // The same equation without its points is the coefficient editor.
  shared_ptr<DetectorPeakResponse> no_points = make_shared<DetectorPeakResponse>( *created );
  no_points->setMeasuredPoints( nullptr );
  BOOST_CHECK( DrfModifyCalc::editorForDrf( *no_points ) == DrfModifyCalc::AnchorEditor::Coefficients );

  // An absolute reference curve (ANGLE .outx): the pairs ARE the points.
  {
    vector<DetectorPeakResponse::EnergyEffPoint> effpts;
    vector<MeasuredEffPoint> measpts;
    for( const double energy : sm_energies )
    {
      const double abs_eff = 1.0E-3 * intrinsic_truth( energy );
      DetectorPeakResponse::EnergyEffPoint e;
      e.energy = static_cast<float>( energy );
      e.efficiency = static_cast<float>( abs_eff );
      effpts.push_back( e );

      MeasuredEffPoint m;
      m.energy = static_cast<float>( energy );
      m.efficiency = static_cast<float>( abs_eff );
      m.distance = static_cast<float>( 25.0 * PhysicalUnits::cm );
      measpts.push_back( m );
    }

    auto angle = make_shared<DetectorPeakResponse>( "angle", "reference curve" );
    angle->setEfficiencyPoints( effpts, static_cast<float>(7.62*PhysicalUnits::cm),
                                25.0*PhysicalUnits::cm,
                                DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    auto meas = make_shared<MeasuredDrfPoints>();
    meas->setPoints( measpts );
    angle->setMeasuredPoints( meas );

    BOOST_CHECK( DrfModifyCalc::editorForDrf( *angle ) == DrfModifyCalc::AnchorEditor::AbsolutePoints );

    // The same curve without raw points only edits the curve.
    shared_ptr<DetectorPeakResponse> csv = make_shared<DetectorPeakResponse>( *angle );
    csv->setMeasuredPoints( nullptr );
    BOOST_CHECK( DrfModifyCalc::editorForDrf( *csv ) == DrfModifyCalc::AnchorEditor::CurvePairs );
  }

  // A formula curve.
  {
    auto formula = make_shared<DetectorPeakResponse>( "formula", "formula curve" );
    formula->setIntrinsicEfficiencyFormula( "exp(-0.45 - 0.62*log(x))",
                                static_cast<float>(7.62*PhysicalUnits::cm),
                                static_cast<float>(PhysicalUnits::MeV), 50.0f, 3000.0f,
                                DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    BOOST_CHECK( DrfModifyCalc::editorForDrf( *formula ) == DrfModifyCalc::AnchorEditor::Formula );
  }
}//BOOST_AUTO_TEST_CASE( editor_follows_what_the_drf_carries )


BOOST_AUTO_TEST_CASE( uncertainty_only_edit_keeps_the_efficiency )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  const DrfModifyCalc::UncertSummary before = DrfModifyCalc::uncertSummary( *orig );
  BOOST_REQUIRE( before.valid );
  BOOST_CHECK( before.total > 0.0 );

  // The user changes ONE uncertainty cell: the certificate uncertainty of the Cs137 source, 3% -> 9%.
  vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
  size_t edited = rows.size();
  for( size_t i = 0; i < rows.size(); ++i )
  {
    if( SpecUtils::istarts_with( rows[i].sourceKey, "Cs137" ) )
    {
      rows[i].fracCert = 0.09;
      edited = i;
    }
  }
  BOOST_REQUIRE( edited < rows.size() );

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
  vector<DrfModifyCalc::Problem> problems;
  const bool applied = DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::RefitPoints,
                                        rows, measured_options(), orig->measuredPoints(), problems );
  BOOST_REQUIRE_MESSAGE( applied, "apply refused: " << problems_to_string(problems) );
  BOOST_CHECK_MESSAGE( !DrfModifyCalc::anyBlocking(problems),
                       "unexpected problem: " << problems_to_string(problems) );

  // The efficiency is a fit result, so changing a weight moves it a little - but only a little.  The
  //  defect this pins made it move by one solid angle (a factor of ~200 for this detector at 25 cm),
  //  and it did so for an edit that touched no efficiency at all.
  const double worst = worst_intrinsic_diff( *working, *orig );
  BOOST_CHECK_MESSAGE( worst < 0.02, "an uncertainty-only edit moved the efficiency by "
                       << 100.0*worst << "% - it must stay the same detector" );

  // Still an equation, with a covariance that the queries will actually use (a rank mismatch is
  //  silently ignored by DetectorEfficiencyCurve::fracCovariance).
  const shared_ptr<const DetectorEfficiencyCurve> curve = working->efficiencyCurve();
  BOOST_REQUIRE( curve && curve->isValid() );
  BOOST_CHECK( curve->form() == DetectorPeakResponse::kExpOfLogPowerSeries );
  const size_t ncoef = curve->expOfLogPowerSeriesCoeffs().size();
  BOOST_CHECK_EQUAL( ncoef, orig->efficiencyCurve()->expOfLogPowerSeriesCoeffs().size() );
  BOOST_REQUIRE( working->efficiencyUncert() );
  BOOST_CHECK_EQUAL( working->efficiencyUncert()->coefficientCovariance().size(), ncoef*ncoef );

  // The edit was not a no-op: a source whose activity is three times less well known has to show up
  //  in what the curve says about itself.  Compared with the response detached, because while one is
  //  attached IT answers every query - the response is built from the curve, and rebuilding it is the
  //  dialog's job on "Use" (see DrfModifyWidget::responseStale).
  shared_ptr<DetectorPeakResponse> curve_only_before = make_shared<DetectorPeakResponse>( *orig );
  shared_ptr<DetectorPeakResponse> curve_only_after = make_shared<DetectorPeakResponse>( *working );
  curve_only_before->setCeeloResponse( nullptr );
  curve_only_after->setCeeloResponse( nullptr );

  const DrfModifyCalc::UncertSummary curve_before = DrfModifyCalc::uncertSummary( *curve_only_before );
  const DrfModifyCalc::UncertSummary curve_after = DrfModifyCalc::uncertSummary( *curve_only_after );
  BOOST_REQUIRE( curve_before.valid && curve_after.valid );
  BOOST_CHECK_MESSAGE( curve_after.total > curve_before.total,
                       "tripling a source's certificate uncertainty did not raise the reported"
                       " uncertainty (" << 100.0*curve_before.total << "% -> "
                       << 100.0*curve_after.total << "%)" );

  // And the stale response really does hide it, which is why the dialog refuses to hand over a DRF
  //  whose response predates the edit.
  const DrfModifyCalc::UncertSummary after = DrfModifyCalc::uncertSummary( *working );
  BOOST_REQUIRE( after.valid );
  BOOST_CHECK_CLOSE( after.total, before.total, 1.0E-6 );

  string why;
  BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *working, why ), why );
}//BOOST_AUTO_TEST_CASE( uncertainty_only_edit_keeps_the_efficiency )


BOOST_AUTO_TEST_CASE( refit_of_unedited_rows_reproduces_the_curve )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  // Applying the rows exactly as they were seeded must give back the curve they came from - the fit
  //  is deterministic and its inputs are unchanged.  (The dialog skips the apply entirely when
  //  nothing was edited; this is the check that the round trip itself is faithful.)
  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::RefitPoints,
                              rows_from_points(points), measured_options(),
                              orig->measuredPoints(), problems ) );

  BOOST_CHECK_SMALL( worst_intrinsic_diff( *working, *orig ), 1.0E-4 );
  BOOST_CHECK( (*working->measuredPoints()) == (*orig->measuredPoints()) );

  // And it is idempotent: a second apply changes nothing further.
  shared_ptr<DetectorPeakResponse> twice = make_shared<DetectorPeakResponse>( *working );
  problems.clear();
  BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *twice, DrfModifyCalc::AnchorEditor::RefitPoints,
                              rows_from_points(*working->measuredPoints()), measured_options(),
                              working->measuredPoints(), problems ) );
  BOOST_CHECK_SMALL( worst_intrinsic_diff( *twice, *working ), 1.0E-5 );
}//BOOST_AUTO_TEST_CASE( refit_of_unedited_rows_reproduces_the_curve )


BOOST_AUTO_TEST_CASE( refit_keeps_everything_it_did_not_edit )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  // The Create-DRF geometry sets a setback from the endcap offset; every absolute efficiency depends
  //  on it, and an edit of an efficiency cell has no business changing it.
  BOOST_REQUIRE( orig->detectorSetback() > 0.0 );

  vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
  rows[3].efficiency *= 1.10;   //a corrected peak area

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::RefitPoints,
                              rows, measured_options(), orig->measuredPoints(), problems ) );

  BOOST_CHECK_CLOSE( working->detectorSetback(), orig->detectorSetback(), 1.0E-6 );
  BOOST_CHECK_CLOSE( working->detectorDiameter(), orig->detectorDiameter(), 1.0E-6 );
  BOOST_CHECK( working->geometryType() == orig->geometryType() );
  BOOST_CHECK( working->geometry() != nullptr );
  BOOST_CHECK( working->drfSource() == orig->drfSource() );
  BOOST_CHECK_EQUAL( working->efficiencyCurve()->energyUnits(), orig->efficiencyCurve()->energyUnits() );

  // The edit did land.
  BOOST_CHECK( worst_intrinsic_diff( *working, *orig ) > 1.0E-3 );

  string why;
  BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *working, why ), why );
}//BOOST_AUTO_TEST_CASE( refit_keeps_everything_it_did_not_edit )


BOOST_AUTO_TEST_CASE( editing_an_energy_keeps_the_points_provenance )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
  const size_t row = 2;
  const float old_energy = static_cast<float>( rows[row].energy );
  rows[row].energy = old_energy + 1.5;    //a mis-identified line, corrected

  // And a source key that the source table does not know about.
  rows[row].sourceKey = "Made-up/SRS-9";

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::RefitPoints,
                              rows, measured_options(), orig->measuredPoints(), problems ) );

  // The unknown key is reported - not silently accepted, and not silently dropped.
  bool warned = false;
  for( const DrfModifyCalc::Problem &p : problems )
    warned = (warned || (p.messageId == "dmw-warn-unknown-source"));
  BOOST_CHECK_MESSAGE( warned, "no note about the unknown source key: " << problems_to_string(problems) );

  const shared_ptr<const MeasuredDrfPoints> after = working->measuredPoints();
  BOOST_REQUIRE( after );
  BOOST_CHECK_EQUAL( after->points().size(), points.points().size() );

  const MeasuredEffPoint *moved = nullptr;
  for( const MeasuredEffPoint &p : after->points() )
  {
    if( fabs( p.energy - (old_energy + 1.5f) ) < 1.0E-3f )
      moved = &p;
  }
  BOOST_REQUIRE_MESSAGE( moved, "the edited point is not in the DRF" );

  // Matching the row to its point by energy - as this used to - rebuilds it from a default-constructed
  //  point the moment the energy is the thing being edited, losing all of this.
  BOOST_CHECK_CLOSE( moved->peakArea, 12345.0f, 1.0E-3 );
  BOOST_CHECK_CLOSE( moved->liveTime, 600.0f, 1.0E-3 );
  BOOST_CHECK_CLOSE( moved->distanceUncert, 0.5f*PhysicalUnits::cm, 1.0E-3 );
  BOOST_CHECK_EQUAL( moved->fileName, string("characterization.n42") );

  // The source table lists exactly the keys the points reference - no orphans, no missing entries.
  for( const MeasuredSourceInfo &src : after->sources() )
  {
    bool referenced = false;
    for( const MeasuredEffPoint &p : after->points() )
      referenced = (referenced || (p.sourceKey == src.sourceKey));
    BOOST_CHECK_MESSAGE( referenced, "source table lists an unreferenced source: " << src.sourceKey );
  }
  for( const MeasuredEffPoint &p : after->points() )
    BOOST_CHECK_MESSAGE( after->sourceForKey( p.sourceKey ),
                         "a point references a source the table does not have: " << p.sourceKey );

  // The assay information of a source that was already known survives.
  const MeasuredSourceInfo * const ba133 = after->sourceForKey( "Ba133/SRS-2" );
  BOOST_REQUIRE( ba133 );
  BOOST_CHECK_EQUAL( ba133->assayInfo, string("certificate 2026-01-01") );
}//BOOST_AUTO_TEST_CASE( editing_an_energy_keeps_the_points_provenance )


BOOST_AUTO_TEST_CASE( absolute_points_edit_keeps_setback_and_air_attenuation )
{
  set_data_dir();

  // An ANGLE-style absolute reference curve, with a setback and air attenuation turned off.
  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  vector<MeasuredEffPoint> measpts;
  for( const double energy : sm_energies )
  {
    const double abs_eff = 1.0E-3 * intrinsic_truth( energy );
    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = static_cast<float>( energy );
    e.efficiency = static_cast<float>( abs_eff );
    e.efficiencyUncert = static_cast<float>( 0.02*abs_eff );
    effpts.push_back( e );

    MeasuredEffPoint m;
    m.energy = static_cast<float>( energy );
    m.efficiency = static_cast<float>( abs_eff );
    m.fracStatUncert = 0.01f;
    m.fracCertUncert = 0.02f;
    m.sourceKey = "ANGLE";
    m.distance = static_cast<float>( 25.0 * PhysicalUnits::cm );
    measpts.push_back( m );
  }

  auto drf = make_shared<DetectorPeakResponse>( "angle", "reference curve" );
  drf->setEfficiencyPoints( effpts, static_cast<float>(7.62*PhysicalUnits::cm),
                            25.0*PhysicalUnits::cm,
                            DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  auto meas = make_shared<MeasuredDrfPoints>();
  meas->setPoints( measpts );
  drf->setMeasuredPoints( meas );
  drf->setDetectorSetback( 1.25 * PhysicalUnits::cm );
  drf->setAbsEffCorrectForAirAtten( false );

  BOOST_REQUIRE( DrfModifyCalc::editorForDrf( *drf ) == DrfModifyCalc::AnchorEditor::AbsolutePoints );

  vector<DrfModifyCalc::PointRow> rows = rows_from_points( *drf->measuredPoints() );
  rows[5].fracStat = 0.05;   //an uncertainty-only edit

  DrfModifyCalc::AnchorOptions options = measured_options();
  options.refDistance = 25.0 * PhysicalUnits::cm;

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *drf );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE_MESSAGE( DrfModifyCalc::applyPointRows( *working,
                              DrfModifyCalc::AnchorEditor::AbsolutePoints, rows, options,
                              drf->measuredPoints(), problems ),
                         problems_to_string(problems) );

  // setEfficiencyPoints is a re-characterization; it resets these, and this edit did not touch them.
  BOOST_CHECK_CLOSE( working->detectorSetback(), 1.25*PhysicalUnits::cm, 1.0E-4 );
  BOOST_CHECK_EQUAL( working->absEffCorrectForAirAtten(), false );
  BOOST_CHECK( working->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  BOOST_CHECK_CLOSE( working->absoluteEfficiencyDistance(), 25.0*PhysicalUnits::cm, 1.0E-4 );

  // The efficiency itself is unchanged (the rows were the curve's own numbers).
  BOOST_CHECK_SMALL( worst_intrinsic_diff( *working, *drf ), 1.0E-5 );

  string why;
  BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *working, why ), why );
}//BOOST_AUTO_TEST_CASE( absolute_points_edit_keeps_setback_and_air_attenuation )


BOOST_AUTO_TEST_CASE( absolute_points_are_moved_to_the_reference_distance )
{
  set_data_dir();

  // An ANGLE-style absolute reference curve anchored at 25 cm.
  const double ref_cm = 25.0;
  const MakeDrfCalc::GeometryChoice geom = nai_geometry();
  CeeLoUtils::GeometryKernel kernel( *geom.geometry );

  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  vector<MeasuredEffPoint> measpts;
  for( const double energy : sm_energies )
  {
    const double abs_eff = intrinsic_truth(energy) * kernel.intrinsicFactor( energy, ref_cm );
    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = static_cast<float>( energy );
    e.efficiency = static_cast<float>( abs_eff );
    effpts.push_back( e );

    MeasuredEffPoint m;
    m.energy = static_cast<float>( energy );
    m.efficiency = static_cast<float>( abs_eff );
    m.distance = static_cast<float>( ref_cm * PhysicalUnits::cm );
    measpts.push_back( m );
  }

  auto drf = make_shared<DetectorPeakResponse>( "angle", "reference curve" );
  drf->setEfficiencyPoints( effpts, static_cast<float>(2.0*geom.geometry->transverse_half_extent()*PhysicalUnits::cm),
                            ref_cm*PhysicalUnits::cm,
                            DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  drf->setGeometry( geom.geometry );
  auto meas = make_shared<MeasuredDrfPoints>();
  meas->setPoints( measpts );
  drf->setMeasuredPoints( meas );

  BOOST_REQUIRE( DrfModifyCalc::editorForDrf( *drf ) == DrfModifyCalc::AnchorEditor::AbsolutePoints );

  // The user re-measures one line at 50 cm and types that distance in.  The curve states absolute
  //  efficiency AT 25 cm, so the point has to be moved there - written in as measured it would be
  //  about four times too small, which is the same class of mistake as reading an absolute
  //  efficiency as an intrinsic one.
  const size_t moved = 7;   //661.7 keV
  const double new_dist_cm = 50.0;
  vector<DrfModifyCalc::PointRow> rows = rows_from_points( *drf->measuredPoints() );
  const double truth_at_50 = intrinsic_truth( rows[moved].energy )
                             * kernel.intrinsicFactor( rows[moved].energy, new_dist_cm );
  rows[moved].distance = new_dist_cm * PhysicalUnits::cm;
  rows[moved].efficiency = truth_at_50;

  DrfModifyCalc::AnchorOptions options = measured_options();
  options.refDistance = ref_cm * PhysicalUnits::cm;

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *drf );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE_MESSAGE( DrfModifyCalc::applyPointRows( *working,
                              DrfModifyCalc::AnchorEditor::AbsolutePoints, rows, options,
                              drf->measuredPoints(), problems ),
                         problems_to_string(problems) );

  // The same physical detector, so the curve must come out where it was.
  BOOST_CHECK_SMALL( worst_intrinsic_diff( *working, *drf ), 1.0E-3 );

  const double curve_at_moved = working->efficiencyCurve()->efficiency(
                                              static_cast<float>(rows[moved].energy) );
  const double expected = intrinsic_truth( rows[moved].energy )
                          * kernel.intrinsicFactor( rows[moved].energy, ref_cm );
  BOOST_CHECK_MESSAGE( fabs( curve_at_moved/expected - 1.0 ) < 1.0E-3,
                       "a point measured at " << new_dist_cm << " cm was written into a curve"
                       " anchored at " << ref_cm << " cm without being transferred: got "
                       << curve_at_moved << ", expected " << expected );

  // The invariant check divides each point's own distance out, so a SYSTEMATIC distance error - the
  //  failure mode it is there for, since it is what "absolute efficiency installed as the intrinsic
  //  curve" looks like - has to be caught.
  {
    shared_ptr<DetectorPeakResponse> wrong = make_shared<DetectorPeakResponse>( *drf );
    vector<MeasuredEffPoint> pts = wrong->measuredPoints()->points();
    for( MeasuredEffPoint &p : pts )
      p.distance = static_cast<float>( new_dist_cm * PhysicalUnits::cm );  //same efficiencies
    auto bad = make_shared<MeasuredDrfPoints>();
    bad->setPoints( pts );
    wrong->setMeasuredPoints( bad );

    string why;
    BOOST_CHECK_MESSAGE( !DrfModifyCalc::checkDrfSelfConsistent( *wrong, why ),
                         "every point moved to the wrong distance was not caught" );
  }

  // One point at a wrong distance, though, is deliberately NOT flagged: the check is on the median
  //  residual, because a single point that fits badly is an ordinary measurement, not a broken DRF,
  //  and a check that fired on one would fire on real detectors and get ignored.
  {
    shared_ptr<DetectorPeakResponse> one_off = make_shared<DetectorPeakResponse>( *drf );
    vector<MeasuredEffPoint> pts = one_off->measuredPoints()->points();
    pts[moved].distance = static_cast<float>( 4.0 * new_dist_cm * PhysicalUnits::cm );
    auto altered = make_shared<MeasuredDrfPoints>();
    altered->setPoints( pts );
    one_off->setMeasuredPoints( altered );

    string why;
    BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *one_off, why ),
                         "one outlying point tripped the self-consistency check: " << why );
  }
}//BOOST_AUTO_TEST_CASE( absolute_points_are_moved_to_the_reference_distance )


BOOST_AUTO_TEST_CASE( an_absolute_equation_is_never_re_fit )
{
  set_data_dir();

  // A Create-DRF detector re-interpreted as absolute efficiency (DrfSelect offers exactly this).
  //  Its equation states absolute efficiency; a re-fit produces an INTRINSIC one, so routing it
  //  through the re-fit editor would leave the detector wrong by a solid angle.
  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  shared_ptr<DetectorPeakResponse> absolute
      = orig->reinterpretAsFarFieldAbsEfficiency( orig->detectorDiameter(),
                                                  25.0*PhysicalUnits::cm, false );
  BOOST_REQUIRE( absolute );
  BOOST_REQUIRE( absolute->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  BOOST_REQUIRE( absolute->measuredPoints() && !absolute->measuredPoints()->empty() );

  BOOST_CHECK_MESSAGE( DrfModifyCalc::editorForDrf( *absolute )
                       == DrfModifyCalc::AnchorEditor::Coefficients,
                       "an absolute-efficiency equation must not be offered the intrinsic re-fit" );

  string warnings;
  BOOST_CHECK_THROW( MakeDrfCalc::refitEfficiencyFromPoints( *absolute, *absolute->measuredPoints(),
                                                             0, warnings ), std::exception );
}//BOOST_AUTO_TEST_CASE( an_absolute_equation_is_never_re_fit )


BOOST_AUTO_TEST_CASE( blanking_the_uncertainty_columns_clears_it )
{
  set_data_dir();

  // A pairs curve that carries a node covariance; the user empties both uncertainty columns.
  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  for( const double energy : sm_energies )
  {
    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = static_cast<float>( energy );
    e.efficiency = static_cast<float>( intrinsic_truth(energy) );
    e.efficiencyUncert = static_cast<float>( 0.03*intrinsic_truth(energy) );
    effpts.push_back( e );
  }

  auto drf = make_shared<DetectorPeakResponse>( "csv", "efficiency csv" );
  drf->setEfficiencyPoints( effpts, static_cast<float>(7.62*PhysicalUnits::cm), 0.0,
                            DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  BOOST_REQUIRE( drf->efficiencyUncert() && drf->efficiencyUncert()->hasNodeCovariance() );

  vector<DrfModifyCalc::PointRow> rows;
  for( size_t i = 0; i < sm_energies.size(); ++i )
  {
    DrfModifyCalc::PointRow row;
    row.rowNumber = static_cast<int>(i) + 1;
    row.energy = sm_energies[i];
    row.efficiency = intrinsic_truth( sm_energies[i] );
    rows.push_back( row );   //both uncertainty cells left at zero
  }

  DrfModifyCalc::AnchorOptions options;
  options.energyUnits = static_cast<float>( PhysicalUnits::keV );

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *drf );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::CurvePairs,
                              rows, options, nullptr, problems ) );

  // Check the store, not just the summary: "no uncertainty reported" would also be true if the apply
  //  had damaged the curve, and the efficiency must not have moved either.
  BOOST_CHECK_MESSAGE( !working->efficiencyUncert() || !working->efficiencyUncert()->hasNodeCovariance(),
                       "emptying every uncertainty cell left the node covariance in place" );
  for( const double energy : sm_energies )
  {
    const double before_eff = drf->intrinsicEfficiency( static_cast<float>(energy) );
    const double after_eff = working->intrinsicEfficiency( static_cast<float>(energy) );
    BOOST_CHECK_MESSAGE( fabs(after_eff - before_eff) <= 1.0E-5*before_eff,
                         "clearing the uncertainty moved the efficiency at " << energy << " keV: "
                         << before_eff << " -> " << after_eff );
  }

  const DrfModifyCalc::UncertSummary after = DrfModifyCalc::uncertSummary( *working );
  BOOST_CHECK_MESSAGE( !after.valid || (after.total <= 0.0),
                       "emptying every uncertainty cell left the detector reporting "
                       << 100.0*after.total << "%" );

  // The same for an equation: zeroing every sigma says "I do not know it", not "keep what you had".
  {
    const MeasuredDrfPoints points = make_measured_points();
    const shared_ptr<DetectorPeakResponse> eqn = make_created_drf( points );
    const vector<float> coefs = eqn->efficiencyCurve()->expOfLogPowerSeriesCoeffs();
    const size_t n = coefs.size();
    BOOST_REQUIRE( eqn->efficiencyUncert() );
    BOOST_REQUIRE_EQUAL( eqn->efficiencyUncert()->coefficientCovariance().size(), n*n );

    vector<double> sigmas( n, 0.0 ), rho( n*n, 0.0 );
    for( size_t i = 0; i < n; ++i )
      rho[i*n + i] = 1.0;

    shared_ptr<DetectorPeakResponse> zeroed = make_shared<DetectorPeakResponse>( *eqn );
    vector<DrfModifyCalc::Problem> probs;
    BOOST_REQUIRE( DrfModifyCalc::applyCoefficients( *zeroed, coefs, sigmas, rho,
                                eqn->efficiencyCurve()->energyUnits(), true, probs ) );
    BOOST_CHECK_MESSAGE( !zeroed->efficiencyUncert()
                         || zeroed->efficiencyUncert()->coefficientCovariance().empty(),
                         "zeroing every coefficient sigma put the old covariance back" );
  }
}//BOOST_AUTO_TEST_CASE( blanking_the_uncertainty_columns_clears_it )


BOOST_AUTO_TEST_CASE( nothing_is_applied_without_a_reason )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  // One usable point left: refused, with a reason, and the DRF untouched.
  {
    vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
    rows.resize( 1 );

    shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
    vector<DrfModifyCalc::Problem> problems;
    BOOST_CHECK( !DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::RefitPoints,
                              rows, measured_options(), orig->measuredPoints(), problems ) );
    BOOST_CHECK( DrfModifyCalc::anyBlocking( problems ) );
    BOOST_CHECK_EQUAL( problems.size(), 1 );
    BOOST_CHECK_EQUAL( problems.front().messageId, string("dmw-err-need-two-points") );
    BOOST_CHECK_EQUAL( working->hashValue(), orig->hashValue() );
  }

  // A point the user added and gave no distance: an absolute efficiency nobody can interpret.  (An
  //  existing row whose distance cell is cleared keeps the distance of the point it describes.)
  {
    vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
    DrfModifyCalc::PointRow added;
    added.rowNumber = static_cast<int>( rows.size() ) + 1;
    added.seedIndex = -1;
    added.energy = 511.0;
    added.efficiency = 1.0E-3;
    added.fracStat = 0.01;
    rows.push_back( added );

    shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
    vector<DrfModifyCalc::Problem> problems;
    BOOST_CHECK( !DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::RefitPoints,
                              rows, measured_options(), orig->measuredPoints(), problems ) );
    BOOST_CHECK( DrfModifyCalc::anyBlocking( problems ) );
    BOOST_CHECK_EQUAL( working->hashValue(), orig->hashValue() );
  }

  // An absolute reference curve with no reference distance.
  {
    vector<DetectorPeakResponse::EnergyEffPoint> effpts;
    for( const double energy : sm_energies )
    {
      DetectorPeakResponse::EnergyEffPoint e;
      e.energy = static_cast<float>( energy );
      e.efficiency = static_cast<float>( 1.0E-3*intrinsic_truth(energy) );
      effpts.push_back( e );
    }
    auto drf = make_shared<DetectorPeakResponse>( "angle", "reference curve" );
    drf->setEfficiencyPoints( effpts, static_cast<float>(7.62*PhysicalUnits::cm),
                              25.0*PhysicalUnits::cm,
                              DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
    auto meas = make_shared<MeasuredDrfPoints>();
    vector<MeasuredEffPoint> measpts;
    for( const DetectorPeakResponse::EnergyEffPoint &e : effpts )
    {
      MeasuredEffPoint m;
      m.energy = e.energy;
      m.efficiency = e.efficiency;
      m.distance = static_cast<float>( 25.0*PhysicalUnits::cm );
      measpts.push_back( m );
    }
    meas->setPoints( measpts );
    drf->setMeasuredPoints( meas );

    shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *drf );
    vector<DrfModifyCalc::Problem> problems;
    DrfModifyCalc::AnchorOptions options = measured_options();  //no refDistance
    BOOST_CHECK( !DrfModifyCalc::applyPointRows( *working, DrfModifyCalc::AnchorEditor::AbsolutePoints,
                              rows_from_points(*meas), options, meas, problems ) );
    BOOST_CHECK( DrfModifyCalc::anyBlocking( problems ) );
    BOOST_CHECK_EQUAL( problems.front().messageId, string("dmw-err-need-ref-distance") );
  }
}//BOOST_AUTO_TEST_CASE( nothing_is_applied_without_a_reason )


BOOST_AUTO_TEST_CASE( the_correlation_control_governs_the_correlated_column )
{
  set_data_dir();

  // A pairs curve with no raw points (a GADRAS Efficiency.csv, an ISOCS .ecc): its rows are the
  //  curve's own numbers, and the "Corr. %" column is tied together across energy by the Baseline
  //  correlation control - all three of its modes, which are expressed as one number.
  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  for( const double energy : sm_energies )
  {
    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = static_cast<float>( energy );
    e.efficiency = static_cast<float>( intrinsic_truth(energy) );
    effpts.push_back( e );
  }

  auto drf = make_shared<DetectorPeakResponse>( "csv", "efficiency csv" );
  drf->setEfficiencyPoints( effpts, static_cast<float>(7.62*PhysicalUnits::cm), 0.0,
                            DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  BOOST_REQUIRE( DrfModifyCalc::editorForDrf( *drf ) == DrfModifyCalc::AnchorEditor::CurvePairs );

  // 3% in "Corr. %", 1% in "Stat. %", on every row.
  vector<DrfModifyCalc::PointRow> rows;
  for( size_t i = 0; i < sm_energies.size(); ++i )
  {
    DrfModifyCalc::PointRow row;
    row.rowNumber = static_cast<int>(i) + 1;
    row.energy = sm_energies[i];
    row.efficiency = intrinsic_truth( sm_energies[i] );
    row.fracStat = 0.01;
    row.fracCert = 0.03;
    rows.push_back( row );
  }

  // The correlation between two well-separated energies, as the stored covariance has it.
  auto correlation_between = []( const DetectorPeakResponse &det, const double e1, const double e2 ) -> double {
    const shared_ptr<const DetectorEfficiencyUncert> uncert = det.efficiencyUncert();
    BOOST_REQUIRE( uncert );
    const vector<double> cov = uncert->efficiencyFracCovariance( { e1, e2 } );
    BOOST_REQUIRE_EQUAL( cov.size(), 4 );
    const double denom = std::sqrt( cov[0] * cov[3] );
    return (denom > 0.0) ? (cov[1] / denom) : 0.0;
  };//correlation_between lambda

  struct ModeCase { const char *name; double corrLength; double lo, hi; };
  const vector<ModeCase> cases = {
    // "Uncorrelated" is EccUncertOptions::effectiveCorrLength() == -1: the correlated column has to
    //  become diagonal too, NOT fall back to a default (which would be the opposite extreme).
    { "uncorrelated", -1.0, -0.01, 0.01 },
    { "gaussian (0.35)", 0.35, 0.0, 0.35 },
    { "fully correlated", DetectorEfficiencyUncert::sm_fullyCorrelatedLength, 0.85, 1.0 },
  };

  for( const ModeCase &mode : cases )
  {
    DrfModifyCalc::AnchorOptions options;
    options.energyUnits = static_cast<float>( PhysicalUnits::keV );
    options.corrLength = mode.corrLength;

    shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *drf );
    vector<DrfModifyCalc::Problem> problems;
    BOOST_REQUIRE_MESSAGE( DrfModifyCalc::applyPointRows( *working,
                                DrfModifyCalc::AnchorEditor::CurvePairs, rows, options,
                                nullptr, problems ), problems_to_string(problems) );

    const double rho = correlation_between( *working, 59.5, 1408.0 );
    BOOST_CHECK_MESSAGE( (rho >= mode.lo) && (rho <= mode.hi),
                         "with " << mode.name << " the 59.5 keV / 1408 keV correlation came out "
                         << rho << ", outside [" << mode.lo << ", " << mode.hi << "]" );

    // The total at one energy is the two columns in quadrature however they are correlated, and
    //  "Stat. %" is diagonal in every mode.
    const vector<double> sig = working->efficiencyUncert()->fracUncertainties( { 661.7 } );
    BOOST_REQUIRE_EQUAL( sig.size(), 1 );
    BOOST_CHECK_CLOSE( sig[0], std::sqrt( 0.03*0.03 + 0.01*0.01 ), 1.0 );
  }//for( const ModeCase &mode : cases )
}//BOOST_AUTO_TEST_CASE( the_correlation_control_governs_the_correlated_column )


BOOST_AUTO_TEST_CASE( impossible_correlations_are_refused )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  const shared_ptr<const DetectorEfficiencyCurve> curve = orig->efficiencyCurve();
  const vector<float> coefs = curve->expOfLogPowerSeriesCoeffs();
  const size_t n = coefs.size();
  BOOST_REQUIRE( n >= 3 );

  // Three coefficients cannot all be 0.99 correlated with one anti-correlated that strongly; the
  //  per-pair |rho| <= 1 clamp cannot see it, and downstream it does not throw - the fit's Cholesky
  //  whitening quietly gives up, so the user would get LESS uncertainty for an impossible entry.
  vector<double> sigmas( n, 0.3 ), rho( n*n, 0.0 );
  for( size_t i = 0; i < n; ++i )
    rho[i*n + i] = 1.0;
  rho[0*n + 1] = rho[1*n + 0] = 0.99;
  rho[0*n + 2] = rho[2*n + 0] = 0.99;
  rho[1*n + 2] = rho[2*n + 1] = -0.99;

  BOOST_CHECK( !DetectorEfficiencyUncert::covarianceIsUsable(
                          DrfModifyCalc::covarianceFromSigmaRho( sigmas, rho ) ) );

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_CHECK( !DrfModifyCalc::applyCoefficients( *working, coefs, sigmas, rho,
                              curve->energyUnits(), true, problems ) );
  BOOST_REQUIRE( !problems.empty() );
  BOOST_CHECK_EQUAL( problems.front().messageId, string("dmw-err-cov-not-psd") );
  BOOST_CHECK_EQUAL( working->hashValue(), orig->hashValue() );

  // A possible set is accepted, and is what the queries then use.
  rho[0*n + 1] = rho[1*n + 0] = 0.90;
  rho[0*n + 2] = rho[2*n + 0] = 0.50;
  rho[1*n + 2] = rho[2*n + 1] = 0.40;
  BOOST_REQUIRE( DetectorEfficiencyUncert::covarianceIsUsable(
                          DrfModifyCalc::covarianceFromSigmaRho( sigmas, rho ) ) );
  problems.clear();
  BOOST_REQUIRE_MESSAGE( DrfModifyCalc::applyCoefficients( *working, coefs, sigmas, rho,
                              curve->energyUnits(), true, problems ),
                         problems_to_string(problems) );
  BOOST_REQUIRE( working->efficiencyUncert() );
  BOOST_CHECK_EQUAL( working->efficiencyUncert()->coefficientCovariance().size(), n*n );

  string why;
  BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *working, why ), why );
}//BOOST_AUTO_TEST_CASE( impossible_correlations_are_refused )


BOOST_AUTO_TEST_CASE( changing_the_term_count_drops_a_covariance_nothing_would_use )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  const shared_ptr<const DetectorEfficiencyCurve> curve = orig->efficiencyCurve();
  vector<float> coefs = curve->expOfLogPowerSeriesCoeffs();
  BOOST_REQUIRE( orig->efficiencyUncert() );
  BOOST_REQUIRE_EQUAL( orig->efficiencyUncert()->coefficientCovariance().size(),
                       coefs.size()*coefs.size() );

  // The user adds a term but does not touch the covariance table (which is therefore not written -
  //  `coefCovTouched` false, the state the widget's per-editor dirty check now reports for exactly
  //  this click sequence).  The stored covariance has the wrong rank, and every query silently
  //  ignores it - so it has to go, with a message, rather than sit there looking like an uncertainty.
  coefs.push_back( 0.0f );

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *orig );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE( DrfModifyCalc::applyCoefficients( *working, coefs, {}, {}, curve->energyUnits(),
                                                   false, problems ) );

  bool told = false;
  for( const DrfModifyCalc::Problem &p : problems )
    told = (told || (p.messageId == "dmw-note-cov-dropped"));
  BOOST_CHECK_MESSAGE( told, "no note about the dropped covariance: " << problems_to_string(problems) );

  BOOST_REQUIRE( working->efficiencyCurve() );
  BOOST_CHECK_EQUAL( working->efficiencyCurve()->expOfLogPowerSeriesCoeffs().size(), coefs.size() );
  // The node covariance the points imply must still be there - only the stale coefficient one goes.
  BOOST_REQUIRE_MESSAGE( working->efficiencyUncert(), "the whole uncertainty was dropped" );
  BOOST_CHECK( working->efficiencyUncert()->coefficientCovariance().empty() );
  BOOST_CHECK( working->efficiencyUncert()->hasNodeCovariance() );

  string why;
  BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *working, why ), why );
}//BOOST_AUTO_TEST_CASE( changing_the_term_count_drops_a_covariance_nothing_would_use )


BOOST_AUTO_TEST_CASE( a_zeroed_sigma_does_not_destroy_correlations )
{
  const size_t n = 3;
  vector<double> sigmas = { 0.5, 0.25, 0.125 };
  vector<double> rho( n*n, 0.0 );
  for( size_t i = 0; i < n; ++i )
    rho[i*n + i] = 1.0;
  rho[0*n + 1] = rho[1*n + 0] = -0.8;
  rho[0*n + 2] = rho[2*n + 0] = 0.4;
  rho[1*n + 2] = rho[2*n + 1] = -0.3;

  const vector<double> before = DrfModifyCalc::covarianceFromSigmaRho( sigmas, rho );
  BOOST_REQUIRE_EQUAL( before.size(), n*n );

  // The editor keeps sigma and rho apart, so a transient 0 (a typo, or clearing the cell to retype
  //  it) is recoverable.  With a covariance as the shadow, scaling a row by 0 lost the correlations
  //  for good - a covariance of zeros has none to scale back up.
  const double kept = sigmas[1];
  sigmas[1] = 0.0;
  const vector<double> zeroed = DrfModifyCalc::covarianceFromSigmaRho( sigmas, rho );
  BOOST_CHECK_SMALL( zeroed[1*n + 1], 1.0E-300 );
  BOOST_CHECK_SMALL( zeroed[0*n + 1], 1.0E-300 );

  sigmas[1] = kept;
  const vector<double> restored = DrfModifyCalc::covarianceFromSigmaRho( sigmas, rho );
  for( size_t i = 0; i < (n*n); ++i )
    BOOST_CHECK_CLOSE( restored[i], before[i], 1.0E-9 );

  // And the round trip through a stored covariance gives the same sigmas and correlations back.
  vector<float> as_float( before.begin(), before.end() );
  vector<double> back_sigmas, back_rho;
  DrfModifyCalc::sigmaRhoFromCovariance( as_float, back_sigmas, back_rho );
  BOOST_REQUIRE_EQUAL( back_sigmas.size(), n );
  for( size_t i = 0; i < n; ++i )
    BOOST_CHECK_CLOSE( back_sigmas[i], sigmas[i], 1.0E-4 );
  for( size_t i = 0; i < (n*n); ++i )
    BOOST_CHECK_CLOSE( back_rho[i], rho[i], 1.0E-3 );
}//BOOST_AUTO_TEST_CASE( a_zeroed_sigma_does_not_destroy_correlations )


BOOST_AUTO_TEST_CASE( fingerprint_tracks_what_a_response_is_built_from )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );
  const size_t base = DrfModifyCalc::seedFingerprint( *orig );

  // Renaming the detector does not make a response stale.
  {
    shared_ptr<DetectorPeakResponse> renamed = make_shared<DetectorPeakResponse>( *orig );
    renamed->setName( "a different name" );
    renamed->setDescription( "and description" );
    BOOST_CHECK_EQUAL( DrfModifyCalc::seedFingerprint( *renamed ), base );
  }

  // Editing the curve does.
  {
    shared_ptr<DetectorPeakResponse> edited = make_shared<DetectorPeakResponse>( *orig );
    vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
    rows[2].efficiency *= 1.05;
    vector<DrfModifyCalc::Problem> problems;
    BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *edited, DrfModifyCalc::AnchorEditor::RefitPoints,
                                rows, measured_options(), orig->measuredPoints(), problems ) );
    BOOST_CHECK( DrfModifyCalc::seedFingerprint( *edited ) != base );
  }

  // So does changing the uncertainty alone...
  {
    shared_ptr<DetectorPeakResponse> edited = make_shared<DetectorPeakResponse>( *orig );
    vector<DrfModifyCalc::PointRow> rows = rows_from_points( points );
    rows[2].fracStat = 0.05;
    vector<DrfModifyCalc::Problem> problems;
    BOOST_REQUIRE( DrfModifyCalc::applyPointRows( *edited, DrfModifyCalc::AnchorEditor::RefitPoints,
                                rows, measured_options(), orig->measuredPoints(), problems ) );
    BOOST_CHECK( DrfModifyCalc::seedFingerprint( *edited ) != base );
  }

  // ...and the coefficient covariance alone, which is exactly what the curve-transfer response
  //  propagates - an edit to it that did not force a regenerate would be an edit the program ignores.
  {
    shared_ptr<DetectorPeakResponse> edited = make_shared<DetectorPeakResponse>( *orig );
    const vector<float> coefs = orig->efficiencyCurve()->expOfLogPowerSeriesCoeffs();
    const size_t n = coefs.size();
    vector<double> sigmas( n ), rho( n*n, 0.0 );
    for( size_t i = 0; i < n; ++i )
    {
      sigmas[i] = 0.02*fabs( coefs[i] ) + 1.0E-4;
      rho[i*n + i] = 1.0;
    }
    vector<DrfModifyCalc::Problem> problems;
    BOOST_REQUIRE( DrfModifyCalc::applyCoefficients( *edited, coefs, sigmas, rho,
                                orig->efficiencyCurve()->energyUnits(), true, problems ) );
    BOOST_CHECK( DrfModifyCalc::seedFingerprint( *edited ) != base );
  }

  // ...and changing the geometry, which is what a transfer or grounding is ray-traced through.
  {
    shared_ptr<DetectorPeakResponse> edited = make_shared<DetectorPeakResponse>( *orig );
    edited->setCeeloResponse( nullptr );   //so `geometry()` returns the member we are about to set
    ceelo::GeometryDescriptor gd = synthetic_nai_descriptor();
    gd.dimensions_cm = { 3.81, 5.0 };      //a shorter crystal
    edited->setGeometry( make_shared<const ceelo::GeometryDescriptor>( gd ) );
    BOOST_CHECK( DrfModifyCalc::seedFingerprint( *edited ) != base );
  }
}//BOOST_AUTO_TEST_CASE( fingerprint_tracks_what_a_response_is_built_from )


BOOST_AUTO_TEST_CASE( a_formula_edit_writes_the_formula_and_its_uncertainty )
{
  set_data_dir();

  auto drf = make_shared<DetectorPeakResponse>( "formula", "an efficiency formula" );
  drf->setIntrinsicEfficiencyFormula( "exp(-6.5 + 0.9*log(x) - 0.25*log(x)^2)",
                                      static_cast<float>(7.62*PhysicalUnits::cm),
                                      static_cast<float>(PhysicalUnits::keV), 50.0f, 3000.0f,
                                      DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  BOOST_REQUIRE( drf->isValid() );
  BOOST_REQUIRE( !drf->efficiencyUncert() || drf->efficiencyUncert()->isEmpty() );

  const double before_661 = drf->intrinsicEfficiency( 661.7f );

  // The rows of a formula editor carry ONLY the uncertainty columns - the formula is the efficiency.
  vector<DrfModifyCalc::PointRow> rows;
  for( size_t i = 0; i < sm_energies.size(); ++i )
  {
    DrfModifyCalc::PointRow row;
    row.rowNumber = static_cast<int>(i) + 1;
    row.energy = sm_energies[i];
    row.fracStat = 0.02;
    row.fracCert = 0.03;
    rows.push_back( row );
  }

  DrfModifyCalc::AnchorOptions options;
  options.energyUnits = static_cast<float>( PhysicalUnits::keV );
  options.corrLength = -1.0;   //fully correlated across energy
  options.hasRowTable = true;

  const string edited_formula = "exp(-6.4 + 0.9*log(x) - 0.25*log(x)^2)";

  shared_ptr<DetectorPeakResponse> working = make_shared<DetectorPeakResponse>( *drf );
  vector<DrfModifyCalc::Problem> problems;
  BOOST_REQUIRE_MESSAGE( DrfModifyCalc::applyFormula( *working, edited_formula,
                              static_cast<float>(PhysicalUnits::keV), rows, options, problems ),
                         problems_to_string(problems) );

  // The formula edit took (that coefficient change is e^0.1 = 1.105x) ...
  const double after_661 = working->intrinsicEfficiency( 661.7f );
  BOOST_CHECK_MESSAGE( fabs( after_661/before_661 - std::exp(0.1) ) < 0.01,
                       "the formula edit did not take: " << before_661 << " -> " << after_661 );

  // ... and the rows became the uncertainty, which a formula curve can only hold as node covariance.
  BOOST_REQUIRE( working->efficiencyUncert() );
  BOOST_CHECK( working->efficiencyUncert()->hasNodeCovariance() );

  const DrfModifyCalc::UncertSummary summary = DrfModifyCalc::uncertSummary( *working );
  BOOST_REQUIRE( summary.valid );
  // sqrt(0.02^2 + 0.03^2) = 3.6%, at a node, with both parts fully correlated across energy.
  BOOST_CHECK_MESSAGE( fabs( summary.total - 0.036 ) < 0.004,
                       "a 2% uncorrelated + 3% correlated formula uncertainty came back as "
                       << 100.0*summary.total << "%" );

  string why;
  BOOST_CHECK_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *working, why ), why );

  // A formula that does not parse is refused, and the DRF is left alone.
  {
    shared_ptr<DetectorPeakResponse> bad = make_shared<DetectorPeakResponse>( *drf );
    vector<DrfModifyCalc::Problem> probs;
    BOOST_CHECK( !DrfModifyCalc::applyFormula( *bad, "exp( -6.4 + ", 
                              static_cast<float>(PhysicalUnits::keV), rows, options, probs ) );
    BOOST_CHECK_MESSAGE( !probs.empty(), "a formula that does not parse was refused silently" );
    BOOST_CHECK_CLOSE( bad->intrinsicEfficiency(661.7f), before_661, 1.0E-4 );
  }
}//BOOST_AUTO_TEST_CASE( a_formula_edit_writes_the_formula_and_its_uncertainty )


BOOST_AUTO_TEST_CASE( geometry_survives_a_response_being_detached )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );
  BOOST_REQUIRE( orig->ceeloResponse() );
  BOOST_REQUIRE( orig->geometry() );

  // Switching to Flat Disk must not make the detector forget what it physically is - without the
  //  geometry no response can ever be attached to it again.
  {
    shared_ptr<DetectorPeakResponse> flat = make_shared<DetectorPeakResponse>( *orig );
    flat->setCeeloResponse( nullptr );
    BOOST_REQUIRE_MESSAGE( flat->geometry(), "detaching the response lost the geometry" );
    BOOST_CHECK_EQUAL( flat->geometry()->to_xml_string(), orig->geometry()->to_xml_string() );

    // ...and it can be attached again from that geometry alone.
    BOOST_CHECK( CeeLoUtils::attachCurveTransferResponse( *flat ) );
    BOOST_CHECK( flat->ceeloResponse() );
  }

  // The same through the two serialization paths a stored DRF takes.
  {
    const string blob = orig->drfExtraToXmlString();
    BOOST_CHECK_MESSAGE( blob.find("CeeLoGeometry") != string::npos,
                         "the DB blob of a response-bearing DRF does not carry its geometry" );

    DetectorPeakResponse restored( *orig );
    restored.setCeeloResponse( nullptr );
    restored.setMeasuredPoints( nullptr );
    restored.setDrfExtraFromXmlString( blob );
    BOOST_REQUIRE( restored.ceeloResponse() );
    restored.setCeeloResponse( nullptr );
    BOOST_REQUIRE_MESSAGE( restored.geometry(),
                           "a DRF loaded from the DB blob lost its geometry when detached" );
  }
}//BOOST_AUTO_TEST_CASE( geometry_survives_a_response_being_detached )


BOOST_AUTO_TEST_CASE( inconsistent_drf_is_detected )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  string why;
  BOOST_REQUIRE_MESSAGE( DrfModifyCalc::checkDrfSelfConsistent( *orig, why ), why );

  // The defect the whole rework exists for, built by hand: the ABSOLUTE measured efficiencies
  //  installed as the INTRINSIC curve.  Everything still "works" - the DRF is valid, the response
  //  still answers queries - but the detector is now ~200x too insensitive on the legacy path.
  {
    shared_ptr<DetectorPeakResponse> wrong = make_shared<DetectorPeakResponse>( *orig );
    vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
    for( const MeasuredEffPoint &p : points.points() )
    {
      DetectorPeakResponse::EnergyEfficiencyPair pair;
      pair.energy = p.energy;
      pair.efficiency = p.efficiency;      //absolute, at that point's own distance
      pairs.push_back( pair );
    }
    auto curve = make_shared<DetectorEfficiencyCurve>();
    curve->setFromPairs( pairs, static_cast<float>(PhysicalUnits::keV) );
    wrong->replaceEfficiencyCurve( curve );

    // Caught by the "do the points still describe this curve" check: the absolute points imply an
    //  intrinsic efficiency ~200x larger than the curve now states.  (Not by the covariance rank
    //  check - `replaceEfficiencyCurve` took the new curve's uncertainty, which is none.)
    BOOST_CHECK( !DrfModifyCalc::checkDrfSelfConsistent( *wrong, why ) );
    BOOST_CHECK_MESSAGE( why.find("measured points") != string::npos,
                         "expected the points/curve disagreement to be named; got: " << why );
  }

  // The same mistake with the equation kept: coefficients that no longer describe the points.
  {
    shared_ptr<DetectorPeakResponse> wrong = make_shared<DetectorPeakResponse>( *orig );
    const shared_ptr<const DetectorEfficiencyCurve> curve = orig->efficiencyCurve();
    vector<float> coefs = curve->expOfLogPowerSeriesCoeffs();
    coefs[0] -= 5.3f;    //exp(-5.3) ~ 1/200: one solid angle
    auto shifted = make_shared<DetectorEfficiencyCurve>();
    shifted->setFromExpOfLogPowerSeries( coefs, {}, curve->energyUnits() );
    shifted->setUncertainty( curve->uncertainty() );
    wrong->replaceEfficiencyCurve( shifted );

    BOOST_CHECK( !DrfModifyCalc::checkDrfSelfConsistent( *wrong, why ) );
    BOOST_CHECK_MESSAGE( why.find("measured points") != string::npos,
                         "expected the points/curve disagreement to be named; got: " << why );
  }
}//BOOST_AUTO_TEST_CASE( inconsistent_drf_is_detected )


BOOST_AUTO_TEST_CASE( uncert_summary_splits_data_from_model )
{
  set_data_dir();

  const MeasuredDrfPoints points = make_measured_points();
  const shared_ptr<DetectorPeakResponse> orig = make_created_drf( points );

  // With the transfer response attached, part of what the detector reports is the response's model
  //  envelope - an allowance no measurement of this detector constrains.
  BOOST_REQUIRE( orig->ceeloResponse() );
  const DrfModifyCalc::UncertSummary with_resp = DrfModifyCalc::uncertSummary( *orig );
  BOOST_REQUIRE( with_resp.valid );
  BOOST_CHECK( with_resp.total > 0.0 );
  BOOST_CHECK_MESSAGE( with_resp.model > 0.0,
                       "a curve-transfer response must report a model envelope" );
  BOOST_CHECK( with_resp.model <= with_resp.total );
  BOOST_CHECK_CLOSE( with_resp.total*with_resp.total,
                     with_resp.data*with_resp.data + with_resp.model*with_resp.model, 1.0E-6 );

  // Back on the flat-disk model, the uncertainty is the fit's own covariance - all data-derived.
  shared_ptr<DetectorPeakResponse> flat = make_shared<DetectorPeakResponse>( *orig );
  flat->setCeeloResponse( nullptr );
  const DrfModifyCalc::UncertSummary curve_only = DrfModifyCalc::uncertSummary( *flat );
  BOOST_REQUIRE( curve_only.valid );
  BOOST_CHECK( curve_only.total > 0.0 );
  BOOST_CHECK_SMALL( curve_only.model, 1.0E-12 );
}//BOOST_AUTO_TEST_CASE( uncert_summary_splits_data_from_model )
