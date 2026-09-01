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


/** Validates InterSpec's VOLUMETRIC (trace / self-attenuating) source efficiency against
 extended-source Monte-Carlo truth.

 The question this answers: the volumetric integrand historically used a flat-disk solid angle times
 the intrinsic efficiency - the far-field approximation - and `VolumetricEffMethod` now offers
 near-field/off-axis-correct alternatives.  Which is right, and by how much, is not decidable from
 self-consistency; it needs an independent truth.

 That truth is `ceelo::EfficiencyCalculator`, run as a full extended-source MC with the source
 matrix modeled, over a real detector (the ANGLE GEM35-70 in test_data/det_eff).  Running it is
 minutes of CPU, so it is split:

   VolumetricNearFieldTruth       (disabled) regenerates the bank and prints a paste-ready table.
   VolumetricNearFieldVsTruth     (runs)     checks InterSpec's integration against the recorded one.

 The recorded numbers are absolute full-energy-peak efficiency per emitted gamma, averaged over the
 source volume.  See VolumetricNearFieldScenarios.h for the scenario matrix and why those corners.
 */

#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

#define BOOST_TEST_MODULE VolumetricNearField_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/AngleOutxImport.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

#include "io/SolidAngle.h"
#include "io/ResponseKernel.h"
#include "io/ResponseGenerator.h"
#include "efficiency/EfficiencyCalculator.h"

#include "VolumetricNearFieldScenarios.h"
#include "VolumetricNearFieldTruth.h"

using namespace std;
using namespace VolNearField;

namespace
{
string g_data_dir, g_test_data_dir;

void set_data_dir()
{
  if( !g_data_dir.empty() )
    return;

  const vector<string> args = boost::unit_test::framework::master_test_suite().argv
        ? vector<string>( boost::unit_test::framework::master_test_suite().argv + 1,
                          boost::unit_test::framework::master_test_suite().argv
                            + boost::unit_test::framework::master_test_suite().argc )
        : vector<string>{};

  for( const string &arg : args )
  {
    if( SpecUtils::starts_with( arg, "--datadir=" ) )
      g_data_dir = arg.substr( 10 );
    else if( SpecUtils::starts_with( arg, "--testfiledir=" ) )
      g_test_data_dir = arg.substr( 14 );
  }

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Pass --datadir=..." );
  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "Pass --testfiledir=..." );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
}


/** The ANGLE GEM35-70 from test_data: its physical geometry and its measured reference efficiency
 curve.  The curve is what an EFFTRAN transfer is anchored on, so this one file supplies both halves
 of the comparison - a real detector and a real characterization. */
struct AngleDetector
{
  ceelo::GeometryDescriptor gd;
  double endcap_front_offset_cm = 0.0;
  shared_ptr<ceelo::DetectorResponse> measured_transfer; ///< anchored on the file's MEASURED curve
  shared_ptr<ceelo::DetectorResponse> mc_transfer;      ///< anchored on the recorded MC curve
};


AngleDetector load_angle_detector()
{
  const string angle_path = SpecUtils::append_path(
        SpecUtils::append_path( g_test_data_dir, "det_eff" ), "Angle-example-efficiency.outx" );

  ifstream in( angle_path.c_str(), ios::in | ios::binary );
  BOOST_REQUIRE_MESSAGE( in.is_open(), "Failed to open ANGLE file '" + angle_path + "'" );

  const AngleOutxContents contents = DetectorPeakResponse::parseAngleOutxFileFull( in );
  BOOST_REQUIRE( contents.hasGeometry );
  BOOST_REQUIRE( contents.hasReference );

  AngleDetector det;
  vector<string> warnings;
  det.gd = CeeLoUtils::buildAngleGeometry( contents, warnings );
  det.endcap_front_offset_cm = det.gd.endcap_front_offset_cm();

  // Anchor an EFFTRAN transfer on the file's own measured curve.
  CeeLoUtils::TransferAnchor anchor;
  anchor.ref_distance_cm = contents.referenceDistanceCm;
  anchor.curve_derived = false;
  for( const pair<float,float> &ep : contents.referencePoints )
  {
    if( (ep.first > 0.0f) && (ep.second > 0.0f) )
    {
      anchor.curve.energies_keV.push_back( ep.first );
      anchor.curve.eff.push_back( ep.second );
    }
  }
  BOOST_REQUIRE_GE( anchor.curve.energies_keV.size(), size_t(2) );

  det.measured_transfer = CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{},
                                                            "ANGLE GEM35-70 (EFFTRAN, measured)" );
  BOOST_REQUIRE( det.measured_transfer );

  // The transfer used for the MODEL comparison is anchored on the recorded Monte-Carlo curve
  //  instead.  Anchoring on the measured curve would make every volumetric row inherit the
  //  ANGLE-measured-vs-CeeLo-MC scale difference (up to ~8.5% at high energy - see
  //  TransferVsMcMeasuredAnchor), which has nothing to do with the volumetric model under test.
  //  Same MC that produced the truth => the comparison isolates the volume integral.
  if( !sm_mc_anchor.empty() )
  {
    CeeLoUtils::TransferAnchor mc_anchor;
    mc_anchor.ref_distance_cm = kMcAnchorDistanceCm;
    mc_anchor.curve_derived = false;
    for( const AnchorRow &row : sm_mc_anchor )
    {
      mc_anchor.curve.energies_keV.push_back( row.energy_keV );
      mc_anchor.curve.eff.push_back( row.eff );
      mc_anchor.curve.frac_sigma.push_back( row.frac_sigma );
    }
    det.mc_transfer = CeeLoUtils::makeTransferResponse( det.gd, mc_anchor, ceelo::AnchorCurve{},
                                                        "ANGLE GEM35-70 (EFFTRAN, MC-anchored)" );
    BOOST_REQUIRE( det.mc_transfer );
  }//if( an MC anchor has been recorded )

  return det;
}


/** InterSpec's volumetric efficiency for one scenario at one energy: the volume-integrated
 per-element detector response, divided by the source volume, so it is directly comparable to the
 MC's efficiency-per-emitted-gamma.

 `eff_response` null => the legacy flat-disk model (solid angle x intrinsic efficiency); non-null =>
 that response is evaluated per element and already carries the intrinsic efficiency.
 */
/** Plan 3.4: survival-removal mu for the FEP leg.

 CeeLo: mu_rem = mu_total - mu_Rayleigh - f_win(E, material) * mu_Compton.  Rayleigh is elastic so
 it cannot remove a photon from the FEP window, and forward Compton scatters whose energy loss
 stays inside the window are still counted in the peak.

 CAREFUL: InterSpec's `transmition_length_coefficient` already EXCLUDES Rayleigh
 (`massAttenuationCoefficientElement` returns compton+photoelectric+pair, and the SNL path returns
 0 for RayleighScatter), so subtracting it again would double-count.  In InterSpec terms:

     mu_rem = transmition_length_coefficient(mat, E) - f_win * mu_Compton(mat, E)

 `win_keV` is a HALF-width, matching CeeLo's +-win convention - which is why its 0.75 keV default
 is about FWHM/2 for an HPGe.
 */
double fep_removal_coefficient( const Material &mat, const double energy_keV, const double win_keV )
{
  const double mu_total = GammaInteractionCalc::transmition_length_coefficient(
                                                    &mat, static_cast<float>(energy_keV) );

  double mu_compton = 0.0;
  for( const Material::ElementFractionPair &p : mat.elements )
    mu_compton += p.second * mat.density
                  * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber,
                        static_cast<float>(energy_keV),
                        MassAttenuation::GammaEmProcces::ComptonScatter );
  for( const Material::NuclideFractionPair &p : mat.nuclides )
    mu_compton += p.second * mat.density
                  * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber,
                        static_cast<float>(energy_keV),
                        MassAttenuation::GammaEmProcces::ComptonScatter );

  // Material-aware in-window fraction: S(x,Z) suppresses forward scatter, most strongly at high Z
  //  (@60 keV: water 0.89, Fe 0.77, Pb 0.63), so the free-electron form over-credits exactly where
  //  the shielding is heaviest.  This is a Simpson integration - hoisted per (E, material) by the
  //  caller, never inside a per-ray loop.
  const ceelo::Material cm = CeeLoUtils::to_ceelo_material( mat ).to_material();
  const double f_win = ceelo::kn_in_window_fraction( energy_keV, win_keV, cm );

  return mu_total - f_win*mu_compton;
}


double interspec_volumetric_eff( const AngleDetector &det,
                                 const Scenario &s,
                                 const double energy_keV,
                                 const shared_ptr<const ceelo::DetectorResponse> &eff_response,
                                 const int n_rays = -1,
                                 const double fep_window_keV = -1.0,
                                 const double epsrel = -1.0 )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;

  // Source matrix, as an InterSpec material, so the model attenuates with InterSpec's own tables -
  //  the SAME material the MC truth was generated through (see scenario_matrix_material).
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> matrix = matdb->material( scenario_matrix_material( s.dense ) );
  BOOST_REQUIRE( matrix );

  DistributedSrcCalcT<double> calc;
  calc.m_geometry = GeometryType::CylinderEndOn;
  calc.m_materialIndex = 0;
  calc.m_attenuateForAir = false;
  calc.m_airTransLenCoef = 0.0;
  calc.m_isInSituExponential = false;
  calc.m_inSituRelaxationLength = -1.0;
  calc.m_srcVolumetricActivity = 1.0;
  calc.m_normalizeByVolume = false;
  calc.m_energy = energy_keV;
  calc.m_nuclide = nullptr;
  calc.integral = 0.0;
  calc.m_effResponse = eff_response;
  calc.m_effMethod = eff_response ? ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer
                                  : ShieldingSourceFitCalc::VolumetricEffMethod::FlatDisk;
  if( n_rays > 0 )
    calc.m_effNumRays = n_rays;

  const double det_radius = 0.5 * det.gd.transverse_half_extent() * 2.0 * cm;
  calc.m_detector = detector_geom_from_config<double>( GeometryType::CylinderEndOn,
                                                       scenario_center_distance_cm( s ) * cm,
                                                       det_radius, 0.0 );

  DistributedSrcCalcT<double>::ShellInfo src_shell;
  src_shell.dims = { s.radius_cm*cm, s.half_length_cm*cm, 0.0 };
  src_shell.trans_len_coef = (fep_window_keV > 0.0)
        ? fep_removal_coefficient( *matrix, energy_keV, fep_window_keV )
        : transmition_length_coefficient( matrix.get(), static_cast<float>(energy_keV) );
  src_shell.type = ShellType::Material;
  calc.m_shells.push_back( src_shell );

  if( s.shield_cm > 0.0 )
  {
    const shared_ptr<const Material> iron = matdb->material( scenario_shield_material() );
    BOOST_REQUIRE( iron );
    DistributedSrcCalcT<double>::ShellInfo shield;
    shield.dims = { (s.radius_cm + s.shield_cm)*cm, (s.half_length_cm + s.shield_cm)*cm, 0.0 };
    shield.trans_len_coef = (fep_window_keV > 0.0)
          ? fep_removal_coefficient( *iron, energy_keV, fep_window_keV )
          : transmition_length_coefficient( iron.get(), static_cast<float>(energy_keV) );
    shield.type = ShellType::Material;
    calc.m_shells.push_back( shield );
  }//if( shielded )

  if( epsrel > 0.0 )
    self_shielding_integration_imp<double>( calc, epsrel, 200000000 );
  else
    self_shielding_integration_imp<double>( calc );

  double eff = calc.integral / (scenario_volume_cm3(s) * cm*cm*cm);

  // FlatDisk carries only the geometric solid angle; the intrinsic efficiency is folded in after
  //  the integration, exactly as expected_peak_counts_imp does.
  if( !eff_response )
  {
    // Back the intrinsic efficiency out of the same response the near-field column uses, so the
    //  near-vs-flat difference is purely the geometry model and carries no anchor offset.
    //
    // CAVEAT - do not read the flat-disk column as "how a flat-disk DRF performs".  The intrinsic
    //  efficiency here is taken at 1000 cm, i.e. the model is pinned at infinity, so it carries the
    //  full near-field difference at EVERY distance including the far ones.  A real flat-disk DRF is
    //  characterized AT some distance (50 cm, say) and is ~exact there by construction, degrading
    //  only as the geometry moves away from that point.  Testing that properly means building a
    //  flat-disk efficiency at a stated distance and extrapolating from it; until then this column
    //  shows only that flat-disk-from-infinity is inadequate near-field - it says nothing
    //  quantitative about where a real one should or should not be used.
    const shared_ptr<const ceelo::DetectorResponse> ref = det.mc_transfer ? det.mc_transfer
                                                                          : det.measured_transfer;
    const ceelo::EffResult far = ref->eps_fep( energy_keV, 0.0, 0.0, 1000.0 );
    const double a_cm = ref->transverse_half_extent();
    const double omega = ceelo::disk_solid_angle_fraction( 1000.0, a_cm );
    eff *= (omega > 0.0) ? (far.value / omega) : 0.0;
  }

  return eff;
}//interspec_volumetric_eff(...)

}//namespace


/** Characteristic optical depth (mean free paths) a full-energy photon must survive to leave the
 source and its shield, at `energy_keV`.  A bound on the escape path, not the largest chord.

 This is the single number that says whether a row is even representable: the model transports the
 UNCOLLIDED photon (plus, once the survival-removal mu lands, a credit for small-angle Compton that
 stays in the peak).  Past a few mfp the surviving signal is dominated by scatter the removal kernel
 cannot represent at all - CeeLo puts that limit explicitly at its 9-mfp / 60 keV case, which its
 depth-aware credit improves but "physically cannot" close.  Printing tau next to each comparison is
 what makes a failing row legible as "outside the model" rather than "model is broken".
 */
double scenario_optical_depth( const Scenario &s, const double energy_keV )
{
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  if( !matdb )
    return 0.0;

  const shared_ptr<const Material> matrix = matdb->material( scenario_matrix_material( s.dense ) );
  const shared_ptr<const Material> shield = matdb->material( scenario_shield_material() );
  if( !matrix || !shield )
    return 0.0;

  const float e = static_cast<float>( energy_keV );

  // The escape path for an END-ON detector is axial, not radial: a photon heads toward the detector
  //  face, so what it must cross is at most the source's full length.  Using the radius instead
  //  wildly overstates a pancake - the wide-angle-far scenario (r=12 cm, half-length 0.5 cm) would
  //  read ~101 mfp while the model in fact tracks truth to a few percent, because no photon ever
  //  crosses 12 cm of material on its way out.  Take the smaller of the two half-extents as the
  //  characteristic depth: it is the one that bounds the escape path.
  const double half_extent_cm = std::min( s.radius_cm, s.half_length_cm );
  const double tau_src = GammaInteractionCalc::transmition_length_coefficient( matrix.get(), e )
                         * half_extent_cm * PhysicalUnits::cm;
  const double tau_shield = (s.shield_cm > 0.0)
        ? GammaInteractionCalc::transmition_length_coefficient( shield.get(), e )
          * s.shield_cm * PhysicalUnits::cm
        : 0.0;

  return tau_src + tau_shield;
}//scenario_optical_depth(...)


/** DEVELOPER-ONLY: regenerates the truth bank and prints a paste-ready table.
 Minutes of CPU; enable with `--run_test=VolumetricNearFieldTruth`.
 */
BOOST_AUTO_TEST_CASE( VolumetricNearFieldTruth, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  const AngleDetector det = load_angle_detector();
  const vector<double> energies = scenario_energies();

  cout << "// Regenerated by VolumetricNearFieldTruth.  FEP window: "
       << ceelo::kDefaultFepWindowKeV << " keV.\n";

  {// Monte-Carlo anchor curve: bare point source on axis at kMcAnchorDistanceCm.
    vector<unique_ptr<ceelo::Material>> anchor_owned;
    ceelo::EfficiencyCalculator anchor_calc;
    ceelo::ResponseGenerator::configure_calculator( anchor_calc, det.gd, anchor_owned );
    anchor_calc.set_point_source(
          Eigen::Vector3d( 0.0, 0.0, -( det.endcap_front_offset_cm + kMcAnchorDistanceCm ) ) );

    cout << "static const std::vector<AnchorRow> sm_mc_anchor = {\n";
    for( const double energy : anchor_energies() )
    {
      ceelo::SimulationConfig cfg;
      cfg.energy_keV = energy;
      cfg.termination.target_fep_rel_precision = 0.005;
      cfg.termination.max_events = 40000000;
      cfg.termination.min_events = 100000;
      cfg.seed = 20260830;

      const ceelo::EfficiencyResult r = anchor_calc.compute( cfg );
      const double eff = r.full_energy_peak_efficiency;
      const double frac = (eff > 0.0) ? (r.fep_uncertainty / eff) : 0.0;

      cout << "  { " << setprecision(12) << energy << ", "
           << setprecision(10) << eff << ", " << setprecision(4) << frac << " },\n";
      cout.flush();
    }//for( anchor energies )
    cout << "};\n\n";
  }// end MC anchor curve

  {// Measured-anchor baseline: transfer-vs-MC for a bare point source, no InterSpec code involved.
    //  Recorded so TransferVsMcMeasuredAnchor stays a cheap check.
    vector<unique_ptr<ceelo::Material>> pt_owned;
    ceelo::EfficiencyCalculator pt_calc;
    ceelo::ResponseGenerator::configure_calculator( pt_calc, det.gd, pt_owned );

    cout << "static const std::vector<PointBaselineRow> sm_point_baseline = {\n";
    for( const double d_cm : { 10.0, 25.0, 50.0 } )
    {
      pt_calc.set_point_source( Eigen::Vector3d( 0.0, 0.0, -( det.endcap_front_offset_cm + d_cm ) ) );
      for( const double energy : energies )
      {
        ceelo::SimulationConfig cfg;
        cfg.energy_keV = energy;
        cfg.termination.target_fep_rel_precision = 0.01;
        cfg.termination.max_events = 20000000;
        cfg.termination.min_events = 50000;
        cfg.seed = 20260829;

        const ceelo::EfficiencyResult mc = pt_calc.compute( cfg );
        cout << "  { " << setprecision(12) << d_cm << ", " << energy << ", "
             << setprecision(10) << mc.full_energy_peak_efficiency << " },\n";
        cout.flush();
      }
    }
    cout << "};\n\n";
  }// end measured-anchor baseline

  cout << "static const std::vector<TruthRow> sm_truth = {\n";

  for( const Scenario &s : scenarios() )
  {
    vector<unique_ptr<ceelo::Material>> owned;
    ceelo::EfficiencyCalculator calc;
    ceelo::ResponseGenerator::configure_calculator( calc, det.gd, owned );

    calc.set_cylindrical_source( scenario_center( s, det.endcap_front_offset_cm ),
                                 s.radius_cm, s.half_length_cm );

    const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
    BOOST_REQUIRE( matdb );
    const shared_ptr<const Material> matrix_mat = matdb->material( scenario_matrix_material(s.dense) );
    BOOST_REQUIRE_MESSAGE( matrix_mat, "No material '"
                           << scenario_matrix_material(s.dense) << "'" );

    const ceelo::Material matrix = CeeLoUtils::to_ceelo_material( *matrix_mat ).to_material();
    calc.set_source_material( &matrix );

    const shared_ptr<const Material> shield_mat = matdb->material( scenario_shield_material() );
    BOOST_REQUIRE( shield_mat );
    const ceelo::Material shield = CeeLoUtils::to_ceelo_material( *shield_mat ).to_material();
    if( s.shield_cm > 0.0 )
      calc.add_source_shield( &shield, s.shield_cm, s.shield_cm );

    for( const double energy : energies )
    {
      ceelo::SimulationConfig cfg;
      cfg.energy_keV = energy;
      cfg.termination.target_fep_rel_precision = 0.01;
      cfg.termination.max_events = 20000000;
      cfg.termination.min_events = 50000;
      cfg.seed = 20260829;

      const ceelo::EfficiencyResult r = calc.compute( cfg );

      cout << "  { \"" << s.name << "\", " << setprecision(12) << energy << ", "
           << setprecision(10) << r.full_energy_peak_efficiency << ", "
           << setprecision(4) << r.fep_uncertainty << " },\n";
      cout.flush();
    }//for( energies )
  }//for( scenarios )

  cout << "};\n";
}//BOOST_AUTO_TEST_CASE( VolumetricNearFieldTruth )


/** Checks InterSpec's volumetric integration against the recorded MC truth.
 
 Two claims, and the second is the point of the whole feature:
   - the near-field (transfer) method reproduces truth wherever the physics is representable, and
   - the legacy flat-disk model agrees far-field but VISIBLY FAILS for a large source at contact.
 If flat-disk ever stops failing there, either the truth bank or the model changed and the feature
 has lost its justification.
 */
BOOST_AUTO_TEST_CASE( VolumetricNearFieldVsTruth )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  BOOST_CHECK_MESSAGE( !sm_truth.empty() && !sm_mc_anchor.empty(),
                       "Truth bank and/or MC anchor is empty - a botched regeneration must not"
                       " green the suite silently." );
  if( sm_truth.empty() || sm_mc_anchor.empty() )
  {
    BOOST_TEST_MESSAGE( "Truth bank or MC anchor is empty - run"
                        " --run_test=VolumetricNearFieldTruth and paste both tables into"
                        " VolumetricNearFieldTruth.h.  Skipping." );
    return;
  }

  const AngleDetector det = load_angle_detector();

  size_t n_checked = 0, n_flat_fails_large_near = 0;

  for( const Scenario &s : scenarios() )
  {
    for( const TruthRow &row : sm_truth )
    {
      if( row.scenario != s.name )
        continue;

      const double truth = row.fep_eff;
      if( !(truth > 0.0) )
        continue;

      const double near = interspec_volumetric_eff( det, s, row.energy_keV, det.mc_transfer );
      const double flat = interspec_volumetric_eff( det, s, row.energy_keV, nullptr );

      // Tolerance: the MC's own 1-sigma, plus a model allowance.  The allowance is per-corner -
      //  a dense source at contact at 60 keV is where the removal-mu approximation is weakest, so
      //  holding it to the same budget as a far-field point would just be pinning noise.
      const bool low_energy = (row.energy_keV < 130.0);
      const bool contact = (s.standoff_cm < 5.0);
      double allowance = 0.05;
      if( contact )     allowance = 0.10;
      if( contact && low_energy && s.dense ) allowance = 0.25;

      const double tol = allowance + 3.0*(row.fep_uncert / std::max(truth, 1.0E-30));

      BOOST_CHECK_MESSAGE( fabs(near - truth) <= tol*truth,
                           s.name << " @ " << row.energy_keV << " keV: near-field model "
                           << setprecision(6) << near << " vs MC truth " << truth
                           << " (allowed " << setprecision(3) << 100.0*tol << "%)" );
      ++n_checked;

      // The justification check: flat-disk must be far off for a large source at contact.
      if( (s.name.find("large-near") != string::npos) && (fabs(flat - truth) > 0.15*truth) )
        ++n_flat_fails_large_near;

      BOOST_TEST_MESSAGE( s.name << " @ " << row.energy_keV << " keV: truth=" << setprecision(6)
                          << truth << "  near=" << near << "  flat=" << flat
                          << "  tau=" << setprecision(3)
                          << scenario_optical_depth( s, row.energy_keV ) << " mfp" );
    }//for( truth rows )
  }//for( scenarios )

  BOOST_CHECK_MESSAGE( n_checked > 0, "No truth rows matched the scenario matrix" );
  BOOST_CHECK_MESSAGE( n_flat_fails_large_near > 0,
                       "Flat-disk agreed with MC truth everywhere for a LARGE source at CONTACT."
                       " That is the case the near-field methods exist for - if it now agrees,"
                       " either the truth bank or the model is wrong." );
}//BOOST_AUTO_TEST_CASE( VolumetricNearFieldVsTruth )


/** Separates two things the volumetric comparison above would otherwise conflate.
 
 The "near" column there is InterSpec's volume integral evaluated through an EFFTRAN transfer that is
 anchored on the ANGLE file's MEASURED efficiency curve, while the truth is CeeLo's MC, which gets
 its absolute scale from first principles.  If those two disagree at some energy, every volumetric
 row inherits that offset - and it would read as volumetric model error when it is nothing of the
 kind.
 
 So: compare the transfer against MC for a POINT source, where there is no volume integral, no source
 self-attenuation, and no InterSpec code involved at all.  Whatever shows up here is the shared
 baseline; only the excess beyond it is attributable to the volumetric model.
 */
BOOST_AUTO_TEST_CASE( TransferVsMcMeasuredAnchor )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  if( sm_point_baseline.empty() )
  {
    BOOST_TEST_MESSAGE( "Point-source baseline not recorded - run"
                        " --run_test=VolumetricNearFieldTruth.  Skipping." );
    return;
  }

  const AngleDetector det = load_angle_detector();

  BOOST_TEST_MESSAGE( "Transfer-vs-MC, bare point source on axis (no InterSpec code involved):" );

  double worst = 0.0;
  for( const PointBaselineRow &row : sm_point_baseline )
  {
    // The MC side is recorded (deterministic, seeded); only the cheap transfer evaluation runs
    //  here, so this stays a fast committable check rather than minutes of Monte Carlo per suite.
    const ceelo::EffResult xfer = det.measured_transfer->eps_fep( row.energy_keV, 0.0, 0.0,
                                                                  row.distance_cm );
    BOOST_REQUIRE( row.mc_fep_eff > 0.0 );
    const double ratio = xfer.value / row.mc_fep_eff - 1.0;
    worst = std::max( worst, fabs(ratio) );

    BOOST_TEST_MESSAGE( "  d=" << row.distance_cm << " cm, " << row.energy_keV
                        << " keV: transfer/MC - 1 = " << setprecision(3) << 100.0*ratio << "%" );
  }//for( recorded baseline rows )

  BOOST_TEST_MESSAGE( "Worst |transfer/MC - 1| over the point-source grid: "
                      << setprecision(3) << 100.0*worst << "%" );

  // Not a tight gate - the point is to MEASURE the shared baseline, and to notice if it ever becomes
  //  large enough that the volumetric comparison is really testing the anchor rather than the model.
  BOOST_CHECK_MESSAGE( worst < 0.25,
                       "The EFFTRAN transfer and the MC disagree by " << 100.0*worst
                       << "% for a bare point source.  Every volumetric row inherits this, so the"
                       " volumetric tolerances are measuring the anchor, not the model." );
}//BOOST_AUTO_TEST_CASE( TransferVsMcMeasuredAnchor )


/** The per-element aperture must reproduce eps_fep, and its ray directions must actually point at
 the detector in InterSpec's frame.

 Two independent things can be silently wrong here and neither shows up as a crash: the
 decomposition (prefactor * sum(w) should equal the bare eps_fep) and the crystal-face -> assembly
 frame rotation (if it is wrong, InterSpec traces attenuation along directions that are not toward
 the detector, and every per-ray result is quietly garbage).  Check both.
 */
BOOST_AUTO_TEST_CASE( ElementApertureFrameAndClosure )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  double worst_closure = 0.0, worst_axis = 0.0;

  for( const double dist_cm : { 2.0, 10.0, 50.0 } )
  {
    for( const double cos_theta : { 1.0, 0.8, 0.4 } )
    {
      // An element->detector direction in the assembly frame with this incidence cosine.  The
      //  detector axis points detector->assembly, so the element looks back along -axis.
      const double sin_theta = std::sqrt( std::max(0.0, 1.0 - cos_theta*cos_theta) );
      const double to_det[3] = { sin_theta, 0.0, cos_theta };

      for( const double energy : { 60.0, 344.0, 1332.5 } )
      {
        const GammaInteractionCalc::ElementAperture ap
              = GammaInteractionCalc::build_element_aperture( *det.mc_transfer, energy,
                                                              dist_cm*PhysicalUnits::cm,
                                                              cos_theta, to_det, 2048 );
        BOOST_REQUIRE( !ap.weights.empty() );
        BOOST_REQUIRE_EQUAL( ap.weights.size(), ap.dirs.size() );

        // (1) Closure: prefactor * sum(w) == the bare point-source eps_fep at this geometry.
        double sum_w = 0.0;
        for( const double w : ap.weights ) sum_w += w;
        const double assembled = ap.prefactor * sum_w;
        const double direct = det.mc_transfer->eps_fep( energy, std::acos(cos_theta), 0.0,
                                                        dist_cm ).value;
        if( direct > 0.0 )
          worst_closure = std::max( worst_closure, fabs(assembled/direct - 1.0) );

        // (2) Frame.  The aperture is NOT symmetric about the element->detector line once it is
        //  near and off-axis (the weights favour long chords), so "mean ray points at the detector"
        //  is not a property to test.  What IS defining of a proper rotation is that it preserves
        //  angles: the angle between the mean ray and the element->detector direction must be the
        //  same in CeeLo's frame as in the assembly frame.  A wrong rotation breaks that; an
        //  asymmetric aperture does not.
        const double theta_q = std::acos( cos_theta );
        const Eigen::Vector3d pos = det.mc_transfer->query_position( theta_q, 0.0, dist_cm );
        const ceelo::ApertureQuadrature q =
              ceelo::make_aperture_quadrature( det.mc_transfer->geometry(), pos, 2048 );
        std::vector<double> w_c;
        std::vector<Eigen::Vector3d> dirs_c;
        det.mc_transfer->fep_ray_weights( energy, q, w_c, dirs_c );
        BOOST_REQUIRE_EQUAL( w_c.size(), ap.weights.size() );

        // Mean ray in each frame.
        double mean_i[3] = { 0.0, 0.0, 0.0 };
        Eigen::Vector3d mean_c( 0.0, 0.0, 0.0 );
        for( size_t i = 0; i < ap.dirs.size(); ++i )
        {
          for( int k = 0; k < 3; ++k )
            mean_i[k] += ap.weights[i] * ap.dirs[i][k];
          mean_c += w_c[i] * dirs_c[i];
        }
        const double mi = std::sqrt( mean_i[0]*mean_i[0] + mean_i[1]*mean_i[1] + mean_i[2]*mean_i[2] );
        const double mc = mean_c.norm();
        BOOST_REQUIRE( (mi > 0.0) && (mc > 0.0) );

        // element->detector direction in each frame
        const double pn = pos.norm();
        BOOST_REQUIRE( pn > 0.0 );
        const Eigen::Vector3d u_c( -pos.x()/pn, -pos.y()/pn, -pos.z()/pn );

        const double ang_i = (mean_i[0]*to_det[0] + mean_i[1]*to_det[1] + mean_i[2]*to_det[2]) / mi;
        const double ang_c = mean_c.dot( u_c ) / mc;
        worst_axis = std::max( worst_axis, fabs(ang_i - ang_c) );
      }//for( energies )
    }//for( cos_theta )
  }//for( distances )

  BOOST_TEST_MESSAGE( "element aperture: worst closure rel-diff " << worst_closure
                      << ", worst frame angle mismatch " << worst_axis );

  BOOST_CHECK_MESSAGE( worst_closure < 1.0e-12,
                       "prefactor*sum(w) does not reproduce eps_fep: " << worst_closure );
  // A proper rotation preserves the angle exactly; only floating-point noise should remain.
  BOOST_CHECK_MESSAGE( worst_axis < 1.0e-12,
                       "cos(mean ray, element->detector) differs between CeeLo's frame and the"
                       " assembly frame by " << worst_axis << ": the rotation is not angle-"
                       "preserving, i.e. it is not the rotation it claims to be." );
}//BOOST_AUTO_TEST_CASE( ElementApertureFrameAndClosure )


/** Verification item 7: the ceres::Jet gradient of the per-ray model must match a finite difference
 of the same model.

 This is the check that would have caught dropping `flat_T` from the composition (plan 3.2).  The
 per-ray kernel is `(flat_T/flat_scalar) * prefactor * sum_i w_i * exp(-tau_i(T))`, where `w_i` and
 `prefactor` are built at the SCALAR element position and frozen for the Jet pass.  So the analytic
 gradient carries two of the three parameter dependences exactly - the 1/r^2 solid angle through
 `flat_T`, and the attenuation through `tau_i(T)` - and drops the third, `d w_i / d dims`.

 A finite difference rebuilds `w_i` at each perturbed geometry, so FD-minus-Jet isolates exactly that
 frozen-weight term.  This test therefore does double duty: it fails loudly if the two carried terms
 are wrong, and it MEASURES the omitted one, which is the honest way to size the approximation
 rather than asserting it is small.  The tolerance below is the measured size at contact, not an
 aspiration; if the aperture weights are ever made T-valued it should tighten to ~1e-6.
 */
BOOST_AUTO_TEST_CASE( PerRayGradientVsFiniteDifference )
{
  using namespace GammaInteractionCalc;
  using Jet1 = ceres::Jet<double,1>;
  const double cm = PhysicalUnits::cm;

  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const std::shared_ptr<const Material> steel = matdb->material( scenario_matrix_material(true) );
  const std::shared_ptr<const Material> iron = matdb->material( scenario_shield_material() );
  BOOST_REQUIRE( steel && iron );

  // Contact geometry, shielded: the regime where the near-field terms are largest and where a
  //  dropped Jacobian would do the most damage to a dimension fit.
  const double src_rad = 2.0, src_hz = 1.0, standoff = 1.0;
  const double det_radius = det.gd.transverse_half_extent() * cm;

  // Builds the calculator.  Either the Fe shield thickness or the source radius can carry the Jet
  //  seed; `radius` is the one that MOVES the source elements, and so the one that exercises the
  //  frozen aperture weights.
  const auto make_calc = [&]( const auto &thickness, const auto &radius, const double energy_keV )
  {
    typedef typename std::decay<decltype(thickness+radius)>::type ScalarT;
    DistributedSrcCalcT<ScalarT> calc;
    calc.m_geometry = GeometryType::CylinderEndOn;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_airTransLenCoef = 0.0;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = ScalarT(1.0);
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy_keV;
    calc.m_nuclide = nullptr;
    calc.integral = ScalarT(0.0);
    calc.m_effResponse = det.mc_transfer;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;

    // Source centre distance is fixed; only the shield grows outward, so the source-to-detector
    //  geometry does not move with the parameter and the FD isolates the shield alone.
    calc.m_detector = detector_geom_from_config<ScalarT>( GeometryType::CylinderEndOn,
                                       ScalarT((standoff + src_hz)*cm), det_radius, 0.0 );

    typename DistributedSrcCalcT<ScalarT>::ShellInfo src;
    src.dims = { ScalarT(radius), ScalarT(src_hz*cm), ScalarT(0.0) };
    src.trans_len_coef = ScalarT( transmition_length_coefficient( steel.get(),
                                                       static_cast<float>(energy_keV) ) );
    src.type = ShellType::Material;
    calc.m_shells.push_back( src );

    typename DistributedSrcCalcT<ScalarT>::ShellInfo shield;
    shield.dims = { ScalarT(radius) + ScalarT(thickness), ScalarT(src_hz*cm) + ScalarT(thickness),
                    ScalarT(0.0) };
    shield.trans_len_coef = ScalarT( transmition_length_coefficient( iron.get(),
                                                          static_cast<float>(energy_keV) ) );
    shield.type = ShellType::Material;
    calc.m_shells.push_back( shield );

    self_shielding_integration_imp<ScalarT>( calc );
    return calc.integral;
  };

  const double t0 = 0.5*cm;
  const double h = 1.0e-4*cm;   // ~0.2% of the thickness: above the integrator's 1e-4 epsrel noise

  for( const double energy : { 60.0, 122.0, 661.7 } )
  {
    const Jet1 t_jet( t0, 0 );
    const Jet1 jet_val = make_calc( t_jet, Jet1(src_rad*cm), energy );

    const double v0 = make_calc( t0, src_rad*cm, energy );
    const double vp = make_calc( t0 + h, src_rad*cm, energy );
    const double vm = make_calc( t0 - h, src_rad*cm, energy );
    const double fd = (vp - vm) / (2.0*h);

    // The value must agree to integrator precision - a Jet lane that disagrees in its scalar part
    //  means the two instantiations are not evaluating the same model at all.
    BOOST_CHECK_MESSAGE( std::fabs(jet_val.a - v0) <= 1.0e-6*std::fabs(v0),
                         "E=" << energy << " keV: Jet scalar lane " << jet_val.a
                         << " != double value " << v0 );

    const double rel = (std::fabs(fd) > 0.0) ? std::fabs(jet_val.v[0] - fd)/std::fabs(fd) : 0.0;
    BOOST_TEST_MESSAGE( "  d(integral)/d(shield thickness) @ " << energy << " keV:"
                        << "  Jet=" << jet_val.v[0] << "  FD=" << fd
                        << "  rel.diff=" << 100.0*rel << "%" );

    // Sign and magnitude must be right: a dropped flat_T showed up as a large relative error here.
    BOOST_CHECK_MESSAGE( jet_val.v[0] * fd > 0.0,
                         "E=" << energy << " keV: Jet gradient (" << jet_val.v[0]
                         << ") and finite difference (" << fd << ") disagree in SIGN" );
    BOOST_CHECK_MESSAGE( rel < 0.05,
                         "E=" << energy << " keV: Jet gradient " << jet_val.v[0]
                         << " differs from finite difference " << fd << " by " << 100.0*rel
                         << "% - larger than the frozen-aperture-weight term should account for." );
  }//for( loop over energies )


  // The source radius: this one DOES move the elements, so the frozen aperture weights bite.  The
  //  gap between the analytic gradient and the finite difference is the size of `d w_i / d radius`
  //  that the Jet pass cannot see - reported rather than assumed.
  const double r0 = src_rad*cm;
  double worst_frozen = 0.0;

  // STEP-SIZE SWEEP FIRST.  A single finite difference is only a reference if it is converged;
  //  the model carries per-element aperture quadrature noise, so too small a step measures noise
  //  and too large a step measures curvature.  Report the sweep so the "frozen-weight gap" below
  //  is read against a stable FD, not against whichever step happened to be picked.
  for( const double energy : { 60.0, 122.0, 661.7 } )
  {
    const Jet1 r_jet( r0, 0 );
    const double jet = make_calc( Jet1(t0), r_jet, energy ).v[0];
    std::ostringstream row;
    row << "  FD step sweep @ " << energy << " keV (Jet=" << jet << "):";
    for( const double frac : { 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 3.0e-2 } )
    {
      const double hh = frac*r0;
      const double fd = (make_calc( t0, r0 + hh, energy ) - make_calc( t0, r0 - hh, energy ))
                          / (2.0*hh);
      row << "  h=" << frac << "r0:" << fd;
    }
    BOOST_TEST_MESSAGE( row.str() );
  }

  // 1% of the radius: above the quadrature noise floor, well below the scale on which the
  //  integrand curves.
  const double rh = 1.0e-2*r0;

  for( const double energy : { 60.0, 122.0, 661.7 } )
  {
    const Jet1 r_jet( r0, 0 );
    const Jet1 jet_val = make_calc( Jet1(t0), r_jet, energy );

    const double v0 = make_calc( t0, r0, energy );
    const double vp = make_calc( t0, r0 + rh, energy );
    const double vm = make_calc( t0, r0 - rh, energy );
    const double fd = (vp - vm) / (2.0*rh);

    BOOST_CHECK_MESSAGE( std::fabs(jet_val.a - v0) <= 1.0e-6*std::fabs(v0),
                         "E=" << energy << " keV: Jet scalar lane " << jet_val.a
                         << " != double value " << v0 );

    const double rel = (std::fabs(fd) > 0.0) ? std::fabs(jet_val.v[0] - fd)/std::fabs(fd) : 0.0;
    worst_frozen = std::max( worst_frozen, rel );
    BOOST_TEST_MESSAGE( "  d(integral)/d(source radius) @ " << energy << " keV:"
                        << "  Jet=" << jet_val.v[0] << "  FD=" << fd
                        << "  frozen-weight gap=" << 100.0*rel << "%" );

    BOOST_CHECK_MESSAGE( jet_val.v[0] * fd > 0.0,
                         "E=" << energy << " keV: Jet gradient (" << jet_val.v[0]
                         << ") and finite difference (" << fd << ") disagree in SIGN for the source"
                         " radius - a dimension fit would step the wrong way." );
    // PINNED BASELINE, not an aspiration.  Measured against the CONVERGED finite difference above:
    //  7.3% @ 60 keV, 5.9% @ 122, 2.2% @ 662.  (An earlier reading of 17.5/4.7/0.8% came from a
    //  single un-swept step and was not a converged derivative - hence the sweep, which now runs
    //  first and is reported, so this number is never quoted against an arbitrary step again.)
    //  Budget set just above the worst so the gap cannot silently grow; it should fall to ~1e-6
    //  once the aperture weights are continuous/T-valued.  Do NOT relax it further without a
    //  matching entry in the eval_cylinder comment.
    BOOST_CHECK_MESSAGE( rel < 0.10,
                         "E=" << energy << " keV: source-radius gradient off by " << 100.0*rel
                         << "% (Jet " << jet_val.v[0] << " vs FD " << fd << ") - worse than the"
                         " pinned frozen-aperture-weight baseline." );
  }//for( loop over energies )

  BOOST_TEST_MESSAGE( "  worst frozen-aperture-weight gradient gap: " << 100.0*worst_frozen << "%" );
}//BOOST_AUTO_TEST_CASE( PerRayGradientVsFiniteDifference )


/** Plan 3.9: how many aperture rays does the per-ray kernel actually need?

 `make_aperture_quadrature` samples a cone that subtends the detector's OUTER BOUNDING SPHERE, and
 falls back to the full 4*pi when the source sits inside that sphere
 (ResponseKernel.cpp, `dist > r_sphere*1.0001`).  The bounding sphere is measured from the middle of
 the outermost shell and spans the whole crystal length plus endcap, so an element in CONTACT with
 the endcap is comfortably inside it: sampling flips from a tight cone to all 4*pi, and only the
 ~20% of rays that happen to subtend the crystal contribute.  Instrumentation on the shielded-near
 scenario measured ~100 surviving rays of 512 at d ~ 1.1 cm against ~270 at d ~ 3 cm - i.e. the
 effective quadrature order is LOWEST exactly where the aperture integrand is hardest.

 The estimator stays unbiased (`w_per_ray = cone_omega_frac / n_rays` carries the sampling measure),
 so this is a resolution question, not a correctness one.  This test answers it directly by refining
 the ray count against a high-order run of the same estimator - MC truth's ~1% sigma cannot resolve
 a sub-percent deterministic bias, which is why the reference here is 8192 rays and not the truth
 bank.
 */
BOOST_AUTO_TEST_CASE( ApertureRayConvergence )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const std::vector<int> ray_counts = { 128, 256, 512, 1024, 2048, 8192 };
  const int reference_rays = ray_counts.back();

  // The contact corners - the only place the sampling mode actually changes.
  std::vector<Scenario> probe;
  for( const Scenario &s : scenarios() )
  {
    if( (s.name == "large-near-dense") || (s.name == "shielded-near-dense")
        || (s.name == "small-near-dense") )
      probe.push_back( s );
  }
  BOOST_REQUIRE_EQUAL( probe.size(), 3u );

  double worst_512 = 0.0;
  std::string worst_512_where;

  for( const Scenario &s : probe )
  {
    for( const double energy : { 60.0, 661.7 } )
    {
      const double ref = interspec_volumetric_eff( det, s, energy, det.mc_transfer,
                                                   reference_rays );
      BOOST_REQUIRE( ref > 0.0 );

      std::ostringstream row;
      row << "  " << std::setw(20) << s.name << " @ " << std::setw(6) << energy << " keV:";
      for( const int n : ray_counts )
      {
        const double v = (n == reference_rays)
                            ? ref : interspec_volumetric_eff( det, s, energy, det.mc_transfer, n );
        const double rel = 100.0*(v/ref - 1.0);
        row << "  n=" << n << ":" << std::fixed << std::setprecision(2) << std::showpos << rel
            << "%" << std::noshowpos;

        if( n == 512 )
        {
          if( std::fabs(rel) > worst_512 )
          {
            worst_512 = std::fabs(rel);
            worst_512_where = s.name + " @ " + std::to_string(energy) + " keV";
          }
        }
      }//for( loop over ray counts )
      BOOST_TEST_MESSAGE( row.str() );
    }//for( loop over energies )
  }//for( loop over probe scenarios )

  BOOST_TEST_MESSAGE( "  worst deviation of the shipped 512-ray setting from the "
                      << reference_rays << "-ray reference: " << worst_512 << "% (" << worst_512_where << ")" );

  // MEASURED: 512 rays sits within 0.10% of the 8192-ray reference across all six contact probes
  //  (128 rays is already within 0.76%).  So the low active-ray fraction at contact costs nothing
  //  in the integrated value - per-element aperture error averages out over the volume integral -
  //  and the ~7-9% residual against MC truth on the near rows is MODEL error, not quadrature.
  //  The budget below is a regression guard set well above the measured 0.10%, not a target.
  BOOST_CHECK_MESSAGE( worst_512 < 0.5,
                       "The default aperture ray count is not converged: 512 rays differ from "
                       << reference_rays << " by " << worst_512 << "% at " << worst_512_where
                       << ".  Raise DistributedSrcCalcT::m_effNumRays." );
}//BOOST_AUTO_TEST_CASE( ApertureRayConvergence )


/** Diagnostic for the sampling-efficiency collapse at contact (see ApertureRayConvergence).

 Reports, along the axis and off it, how many of the sampled rays survive to carry the FEP kernel.
 Three counts matter and they are not the same number:
   - `n_rays_total`  : rays actually sampled.
   - `rays.size()`   : rays hitting ANY geometry, endcap and housing included.
   - active          : rays with a non-zero ACTIVE-crystal chord - the only ones that carry FEP
                       weight, and therefore the true order of the aperture quadrature.
 */
BOOST_AUTO_TEST_CASE( ApertureSamplingEfficiency )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const ceelo::Geometry &geom = det.mc_transfer->geometry();
  const double r_out = geom.outer_bounding_radius();
  const std::pair<double,double> zext = geom.outer_z_extent();
  const double half_z = 0.5*(zext.second - zext.first);
  const double r_sphere = std::sqrt( r_out*r_out + half_z*half_z );

  BOOST_TEST_MESSAGE( "  detector bounding sphere: r_out=" << r_out << " cm, half_z=" << half_z
                      << " cm  =>  r_sphere=" << r_sphere << " cm, centered "
                      << -0.5*(zext.first + zext.second) << " cm behind the face" );
  BOOST_TEST_MESSAGE( "  (full-sphere sampling is used whenever the element is inside r_sphere)" );

  // What the ACTIVE volume actually looks like - a bulletized rim or a bore hole gives the scoring
  //  silhouette curvature and the chord length a steep radial gradient near the edge, which is
  //  exactly what a coarse aperture would smear.  Recorded so the convergence result below cannot
  //  be read as covering a feature the test geometry does not actually have.
  BOOST_TEST_MESSAGE( "  active crystal: radius=" << geom.detector_radius()
                      << " cm, length=" << geom.detector_length() << " cm"
                      << ", bulletized=" << (geom.has_bullet_radius() ? "YES" : "no")
                      << " (r=" << geom.bullet_radius() << " cm)"
                      << ", bore hole=" << (geom.has_bore_hole() ? "YES" : "no")
                      << ", dead layer=" << (geom.has_dead_layer() ? "YES" : "no")
                      << ", attenuators=" << geom.num_attenuators() );

  const int n_sampled = 512;
  double worst_active_frac = 1.0;

  for( const double cos_theta : { 1.0, 0.7 } )
  {
    for( const double d : { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 } )
    {
      const Eigen::Vector3d pos = det.mc_transfer->query_position( std::acos(cos_theta), 0.0, d );
      const ceelo::ApertureQuadrature q = ceelo::make_aperture_quadrature( geom, pos, n_sampled );

      size_t n_active = 0;
      for( const ceelo::KernelRay &r : q.rays )
        n_active += (r.active_len > 0.0f) ? 1u : 0u;

      const double active_frac = double(n_active) / double(q.n_rays_total);
      worst_active_frac = std::min( worst_active_frac, active_frac );

      BOOST_TEST_MESSAGE( "  cos_theta=" << cos_theta << " d=" << std::setw(5) << d << " cm:"
                          << "  sampled=" << q.n_rays_total
                          << "  hit-geometry=" << q.rays.size()
                          << "  active=" << n_active
                          << "  cone_omega_frac=" << std::fixed << std::setprecision(4)
                          << q.cone_omega_frac
                          << "  omega_frac_active=" << q.omega_frac_active
                          << (q.cone_omega_frac >= 0.9999 ? "   [FULL SPHERE]" : "") );
    }//for( loop over distances )
  }//for( loop over incidence )

  BOOST_TEST_MESSAGE( "  worst active-ray fraction: " << 100.0*worst_active_frac << "% of sampled" );

  // Purely a report; the assertion just pins that SOME rays score everywhere probed, so a geometry
  //  or reference-point regression cannot make this silently vacuous.
  BOOST_CHECK_MESSAGE( worst_active_frac > 0.0,
                       "No active rays at some probed position - the aperture would be empty." );
}//BOOST_AUTO_TEST_CASE( ApertureSamplingEfficiency )


/** Would a cylinder-exact sampling cone beat the bounding-sphere one?

 `make_aperture_quadrature` circumscribes the geometry in a SPHERE and cones around that, falling
 back to the full 4*pi when the source is inside the sphere.  But the geometry is a cylinder (or
 box), and CeeLo already carries the exact bound elsewhere - `Geometry::ray_hits_outer_boundary()`
 does exact cylinder/box intersection and is even documented "useful for cone-sampling checks", and
 `SourceGeometry::set_detector_bounds()` stores a bounding CYLINDER for electron filtering.  The
 sphere circumscribing a cylinder of radius R and half-length h has radius sqrt(R^2+h^2), so for a
 roughly equilateral detector it is ~40% wider than the object it bounds, and - the part that
 actually bites - it protrudes in FRONT of the crystal face, which is what forces contact sources
 onto the 4*pi fallback.

 This test does not change the sampler; it measures the headroom, so the decision to touch vendored
 code is made on numbers.  The cylinder-exact cone is built by taking the outer shell's two rim
 circles (a cylinder is the convex hull of them, so the enclosing cone is determined by them),
 aiming the axis at the mean direction, then one bisector refinement.
 */
BOOST_AUTO_TEST_CASE( CylinderExactConeHeadroom )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const ceelo::Geometry &geom = det.mc_transfer->geometry();
  const double R = geom.outer_bounding_radius();
  const std::pair<double,double> zext = geom.outer_z_extent();

  // Directions to the outer shell's rim circles, as seen from `src`.
  const auto rim_dirs = [&]( const Eigen::Vector3d &src )
  {
    std::vector<Eigen::Vector3d> dirs;
    for( const double z : { zext.first, zext.second } )
    {
      for( int i = 0; i < 180; ++i )
      {
        const double ph = 2.0*M_PI*i/180.0;
        const Eigen::Vector3d p( R*std::cos(ph), R*std::sin(ph), z );
        dirs.push_back( (p - src).normalized() );
      }
    }
    return dirs;
  };

  // Smallest cone containing `dirs`.  The geometry is axisymmetric, so the optimal cone axis lies
  //  in the plane spanned by the detector axis and the source's radial offset - a 1-D search, which
  //  a scan plus local refinement solves robustly.  (A mean-direction-plus-bisector heuristic was
  //  tried first and produced cones WIDER than the sphere's at 2-25 cm, which is impossible for a
  //  body inside the sphere: it was the optimizer failing, not the idea.)
  const auto cone_frac = [&]( const std::vector<Eigen::Vector3d> &dirs, const Eigen::Vector3d &src )
  {
    Eigen::Vector3d e_r( src.x(), src.y(), 0.0 );
    e_r = (e_r.norm() > 1.0e-9) ? e_r.normalized() : Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d e_z( 0.0, 0.0, 1.0 );

    const auto max_half_angle_cos = [&]( const double psi )
    {
      const Eigen::Vector3d axis = std::cos(psi)*e_z + std::sin(psi)*e_r;
      double min_cos = 2.0;
      for( const Eigen::Vector3d &d : dirs )
        min_cos = std::min( min_cos, axis.dot(d) );
      return min_cos;
    };

    double best_psi = 0.0, best_cos = -2.0;
    for( int i = 0; i < 720; ++i )
    {
      const double psi = 2.0*M_PI*i/720.0;
      const double c = max_half_angle_cos( psi );
      if( c > best_cos ){ best_cos = c; best_psi = psi; }
    }

    // Local refinement around the scan winner.
    double step = 2.0*M_PI/720.0;
    for( int iter = 0; iter < 40; ++iter )
    {
      const double cm_ = max_half_angle_cos( best_psi - step );
      const double cp_ = max_half_angle_cos( best_psi + step );
      if( cm_ > best_cos ){ best_cos = cm_; best_psi -= step; }
      else if( cp_ > best_cos ){ best_cos = cp_; best_psi += step; }
      else step *= 0.5;
    }

    best_cos = std::max( -1.0, std::min( 1.0, best_cos ) );
    return std::make_pair( std::min( 1.0, 0.5*(1.0 - best_cos) ),
                           Eigen::Vector3d(std::cos(best_psi)*e_z + std::sin(best_psi)*e_r) );
  };

  BOOST_TEST_MESSAGE( "  outer shell: R=" << R << " cm, z in [" << zext.first << ", "
                      << zext.second << "] cm  (sphere radius "
                      << std::sqrt(R*R + 0.25*(zext.second-zext.first)*(zext.second-zext.first))
                      << " cm)" );

  double best_gain = 1.0, worst_gain = 1.0e9;

  for( const double cos_theta : { 1.0, 0.7 } )
  {
    for( const double d : { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 } )
    {
      const Eigen::Vector3d pos = det.mc_transfer->query_position( std::acos(cos_theta), 0.0, d );
      const ceelo::ApertureQuadrature q = ceelo::make_aperture_quadrature( geom, pos, 512 );

      const std::pair<double,Eigen::Vector3d> cone = cone_frac( rim_dirs(pos), pos );
      const double exact = cone.first;
      const double gain = (exact > 0.0) ? (q.cone_omega_frac / exact) : 1.0;
      best_gain = std::max( best_gain, gain );
      worst_gain = std::min( worst_gain, gain );

      // Sanity: the "exact" cone must actually contain every ray that hit the geometry, or it is
      //  not a bound and the whole idea is unsafe.
      size_t outside = 0;
      {
        const Eigen::Vector3d axis = cone.second;
        const double cos_lim = 1.0 - 2.0*exact;
        for( const ceelo::KernelRay &r : q.rays )
        {
          if( r.active_len <= 0.0f )
            continue;
          if( axis.dot( r.dir.cast<double>() ) < cos_lim - 1.0e-6 )
            ++outside;
        }
      }

      BOOST_TEST_MESSAGE( "  cos_theta=" << cos_theta << " d=" << std::setw(5) << d
                          << " cm:  sphere-cone=" << std::fixed << std::setprecision(4)
                          << q.cone_omega_frac << "  cylinder-cone=" << exact
                          << "  => " << std::setprecision(2) << gain << "x fewer wasted rays"
                          << (outside ? "   [!! " + std::to_string(outside)
                                        + " scoring rays OUTSIDE the tighter cone]" : "") );

      BOOST_CHECK_MESSAGE( outside == 0,
                           "cos_theta=" << cos_theta << " d=" << d << ": " << outside
                           << " scoring rays fall outside the cylinder-exact cone - it is not a"
                           " valid bound, do not adopt it." );
    }//for( loop over distances )
  }//for( loop over incidence )

  BOOST_TEST_MESSAGE( "  ray-yield gain from a cylinder-exact cone: " << worst_gain << "x to "
                      << best_gain << "x" );
}//BOOST_AUTO_TEST_CASE( CylinderExactConeHeadroom )


/** PER-ELEMENT aperture precision, which is NOT what ApertureRayConvergence measures.

 That test refines the volume-INTEGRATED efficiency, where thousands of elements average; it found
 512 rays within 0.10% of 8192.  That result cannot be read as "512 rays resolve one element well",
 because averaging can hide a large per-element scatter.  Per-element precision is what governs the
 gradient staircase, and it is worst exactly where this geometry is most awkward: at contact, out at
 the BULLETIZED rim (bullet r=0.8 cm on a 2.915 cm crystal) and over the BORE HOLE, where the active
 chord length changes fastest with ray direction.

 So: evaluate one element's FEP efficiency (prefactor * sum w_i, which closure-matches eps_fep) at
 increasing ray counts against a 32768-ray reference, on a grid that deliberately includes the rim
 and the bore shadow.
 */
BOOST_AUTO_TEST_CASE( PerElementAperturePrecision )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const ceelo::Geometry &geom = det.mc_transfer->geometry();
  const double r_cryst = geom.detector_radius();

  const std::vector<int> counts = { 128, 512, 2048, 8192 };
  const int reference = 32768;
  const double energy = 122.0;

  // (axial distance from face, radial offset) - on axis, mid, at the rim, and outside it.
  const std::vector<std::pair<double,double>> probes = {
    { 0.5, 0.0 }, { 0.5, 0.5*r_cryst }, { 0.5, r_cryst }, { 0.5, 1.4*r_cryst },
    { 1.0, 0.0 }, { 1.0, 0.5*r_cryst }, { 1.0, r_cryst }, { 1.0, 1.4*r_cryst },
    { 4.0, 0.0 },                       { 4.0, r_cryst },
  };

  const auto eff_at = [&]( const double d_ax, const double r_off, const int n_rays )
  {
    const double dist = std::sqrt( d_ax*d_ax + r_off*r_off );
    const double cos_theta = (dist > 0.0) ? (d_ax/dist) : 1.0;
    const Eigen::Vector3d pos = det.mc_transfer->query_position( std::acos(cos_theta), 0.0, dist );
    const ceelo::ApertureQuadrature q = ceelo::make_aperture_quadrature( geom, pos, n_rays );

    std::vector<double> w;
    std::vector<Eigen::Vector3d> dirs;
    det.mc_transfer->fep_ray_weights( energy, q, w, dirs );
    double sum_w = 0.0;
    size_t n_active = 0;
    for( size_t i = 0; i < w.size(); ++i )
    {
      sum_w += w[i];
      n_active += (w[i] > 0.0) ? 1u : 0u;
    }
    const ceelo::EffResult pre = det.mc_transfer->fep_prefactor( energy, pos, q );
    return std::make_pair( pre.value * sum_w, n_active );
  };

  double worst_512 = 0.0;
  std::string worst_where;

  BOOST_TEST_MESSAGE( "  crystal radius " << r_cryst << " cm, bulletized r="
                      << geom.bullet_radius() << " cm, bore hole="
                      << (geom.has_bore_hole() ? "yes" : "no") << ", at " << energy << " keV" );

  for( const std::pair<double,double> &p : probes )
  {
    const std::pair<double,size_t> ref = eff_at( p.first, p.second, reference );
    BOOST_REQUIRE( ref.first > 0.0 );

    std::ostringstream row;
    row << "  d_axial=" << std::fixed << std::setprecision(1) << p.first
        << " r_offset=" << std::setprecision(2) << std::setw(5) << p.second << " cm ("
        << std::setw(4) << std::setprecision(0) << 100.0*p.second/r_cryst << "% of R):";
    for( const int n : counts )
    {
      const std::pair<double,size_t> v = eff_at( p.first, p.second, n );
      const double rel = 100.0*(v.first/ref.first - 1.0);
      row << "  n=" << n << ":" << std::showpos << std::fixed << std::setprecision(2) << rel << "%"
          << std::noshowpos << "(" << v.second << " act)";
      if( n == 512 && std::fabs(rel) > worst_512 )
      {
        worst_512 = std::fabs(rel);
        std::ostringstream w; w << "d=" << p.first << " r=" << p.second;
        worst_where = w.str();
      }
    }
    BOOST_TEST_MESSAGE( row.str() );
  }//for( loop over probes )

  BOOST_TEST_MESSAGE( "  worst PER-ELEMENT error of the shipped 512-ray setting: "
                      << worst_512 << "% (" << worst_where << ")" );

  // Pinned as a report, generously: this is the number that should drive any decision to raise
  //  m_effNumRays, and it must not silently drift.
  BOOST_CHECK_MESSAGE( worst_512 < 5.0,
                       "Per-element aperture error at 512 rays is " << worst_512 << "% at "
                       << worst_where << " - large enough that m_effNumRays should be raised." );
}//BOOST_AUTO_TEST_CASE( PerElementAperturePrecision )


/** TEMPORARY A/B (plan 3.4): does the survival-removal mu close the low-energy near-field gap? */
BOOST_AUTO_TEST_CASE( SurvivalRemovalMuSweep )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  if( sm_truth.empty() ) return;

  // The window is NOT a free parameter: CeeLo's convention is a +-half-width, so it is FWHM/2 of
  //  the detector the transfer is anchored on. A GEM35-70 HPGe runs ~0.80 keV FWHM at 122 keV and
  //  ~1.85 keV at 1332 keV; fitting FWHM = sqrt(a + b*E) to those gives a=0.359, b=0.00230, i.e.
  //  half-widths of 0.35 keV at 60 keV rising to 0.93 keV at 1332 keV. A CONSTANT 0.75 keV
  //  therefore roughly DOUBLES the credit at 60 keV, which is exactly where the credit is largest.
  std::map<std::string,Scenario> byname;
  for( const Scenario &s : scenarios() ) byname[s.name] = s;

  for( const TruthRow &row : sm_truth )
  {
    const std::map<std::string,Scenario>::const_iterator it = byname.find( row.scenario );
    if( it == byname.end() ) continue;
    const double fwhm = std::sqrt( 0.359 + 0.00230*row.energy_keV );
    const std::vector<double> windows = { 0.5*fwhm, 0.75, 1.5 };
    std::ostringstream o;
    o << std::setw(22) << row.scenario << " @" << std::setw(7) << row.energy_keV << " keV:  mu_total="
      << std::showpos << std::fixed << std::setprecision(1)
      << 100.0*(interspec_volumetric_eff(det,it->second,row.energy_keV,det.mc_transfer)/row.fep_eff - 1.0) << "%";
    const char * const tags[] = { "FWHM/2", "0.75  ", "1.50  " };
    for( size_t wi = 0; wi < windows.size(); ++wi )
    {
      const double w = windows[wi];
      o << "   " << tags[wi] << "(" << std::noshowpos << std::setprecision(2) << w << "):"
        << std::setprecision(1) << std::showpos
        << 100.0*(interspec_volumetric_eff(det,it->second,row.energy_keV,det.mc_transfer,-1,w)/row.fep_eff - 1.0) << "%";
    }
    BOOST_TEST_MESSAGE( o.str() );
  }
}


/** How much can the survival credit possibly move the answer?  mu_rem = mu_total - f_win*mu_CS, so
 the ENTIRE headroom is f_win * (mu_CS/mu_total). Where photoelectric dominates, there is almost
 nothing to credit no matter how the window is chosen. */
BOOST_AUTO_TEST_CASE( RemovalMuHeadroom )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  for( const bool dense : { true, false } )
  {
    const std::shared_ptr<const Material> m = matdb->material( scenario_matrix_material(dense) );
    BOOST_REQUIRE( m );
    const ceelo::Material cm = CeeLoUtils::to_ceelo_material( *m ).to_material();
    BOOST_TEST_MESSAGE( "  " << (dense ? "steel" : "water") << ":" );
    for( const double e : scenario_energies() )
    {
      const double mu_tot = GammaInteractionCalc::transmition_length_coefficient( m.get(), (float)e );
      double mu_cs = 0.0;
      for( const Material::ElementFractionPair &pr : m->elements )
        mu_cs += pr.second * m->density * MassAttenuation::massAttenuationCoefficientElement(
                   pr.first->atomicNumber, (float)e, MassAttenuation::GammaEmProcces::ComptonScatter );
      const double fwhm = std::sqrt( 0.359 + 0.00230*e );
      const double f_free = ceelo::kn_in_window_fraction( e, 0.5*fwhm );
      const double f_mat  = ceelo::kn_in_window_fraction( e, 0.5*fwhm, cm );
      BOOST_TEST_MESSAGE( "    E=" << std::setw(7) << e << " keV  win=" << std::fixed
        << std::setprecision(2) << 0.5*fwhm
        << "  mu_CS/mu_tot=" << std::setprecision(3) << mu_cs/mu_tot
        << "  f_win(free)=" << f_free << "  f_win(material)=" << f_mat
        << "  => max credit = " << std::setprecision(2) << 100.0*f_mat*mu_cs/mu_tot << "% of mu" );
    }
  }
}


/** DIAGNOSTIC: is the near-field residual in the TRANSFER, or in the volumetric integration?
 *
 * `VolumetricNearFieldVsTruth` compares a volume-integrated model against extended-source MC, so
 * its residual bundles two separable things: how well the transfer response reproduces the real
 * detector at a near-field position at all, and how well InterSpec integrates that response over a
 * source volume. `sm_point_baseline` cannot separate them - it only covers 10/25/50 cm, all
 * far-field.
 *
 * This runs bare POINT sources at the standoffs the near scenarios actually sample and compares
 * MC directly against `eps_fep` from the same transfer response, with no InterSpec integration in
 * the loop at all. Whatever error shows up here is the transfer's own near-field accuracy and is a
 * floor on what the volumetric comparison can achieve; anything beyond it is ours.
 *
 * Developer-only: real MC, minutes of CPU.
 */
BOOST_AUTO_TEST_CASE( PointSourceNearFieldTransferAccuracy, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  std::vector<std::unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator calc;
  ceelo::ResponseGenerator::configure_calculator( calc, det.gd, owned );

  BOOST_TEST_MESSAGE( "  bare POINT source on axis: MC truth vs the transfer response's eps_fep" );
  BOOST_TEST_MESSAGE( "  (no InterSpec volume integration involved)" );

  for( const double d_cm : { 1.0, 1.5, 2.0, 3.0, 5.0, 10.0 } )
  {
    calc.set_point_source( Eigen::Vector3d( 0.0, 0.0, -( det.endcap_front_offset_cm + d_cm ) ) );

    std::ostringstream o;
    o << "    d=" << std::setw(5) << d_cm << " cm: ";
    for( const double energy : scenario_energies() )
    {
      ceelo::SimulationConfig cfg;
      cfg.energy_keV = energy;
      cfg.termination.target_fep_rel_precision = 0.01;
      cfg.termination.max_events = 20000000;
      cfg.termination.min_events = 50000;
      cfg.seed = 20260830;

      const ceelo::EfficiencyResult mc = calc.compute( cfg );
      const ceelo::EffResult tr = det.mc_transfer->eps_fep( energy, 0.0, 0.0, d_cm );
      const double rel = (mc.full_energy_peak_efficiency > 0.0)
              ? 100.0*(tr.value/mc.full_energy_peak_efficiency - 1.0) : 0.0;
      o << " " << std::setw(6) << std::setprecision(0) << std::fixed << energy << "keV:"
        << std::showpos << std::setprecision(1) << std::setw(6) << rel << "%" << std::noshowpos;
    }
    BOOST_TEST_MESSAGE( o.str() );
  }
}


/** DIAGNOSTIC: is the low-energy near-field deficit just integration resolution?
 *
 * At 60 keV in steel mu ~ 9.4/cm, so a 4 cm object emits from a ~1 mm surface skin: the integrand
 * is a thin shell against a nearly-empty interior. An adaptive rule at epsrel=1e-4 may simply not
 * resolve it, and under-resolving a sharply peaked integrand under-estimates. That would be
 * energy-dependent (worse where mu is large), material-dependent (dense worse than light) and
 * roughly geometry-independent - which is exactly the observed pattern.
 */
BOOST_AUTO_TEST_CASE( IntegrationResolutionSweep )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  if( sm_truth.empty() ) return;

  std::map<std::string,Scenario> byname;
  for( const Scenario &s : scenarios() ) byname[s.name] = s;

  for( const TruthRow &row : sm_truth )
  {
    if( row.scenario.find("near") == std::string::npos ) continue;
    if( row.energy_keV > 130.0 ) continue;
    const std::map<std::string,Scenario>::const_iterator it = byname.find( row.scenario );
    if( it == byname.end() ) continue;

    std::ostringstream o;
    o << "  " << std::setw(20) << row.scenario << " @" << std::setw(6) << row.energy_keV << " keV:";
    for( const double eps : { -1.0, 1.0e-5, 1.0e-6, 1.0e-7 } )
    {
      const double v = interspec_volumetric_eff( det, it->second, row.energy_keV,
                                                 det.mc_transfer, -1, -1.0, eps );
      o << (eps < 0 ? "  default:" : "  ") ;
      if( eps > 0 ) o << std::scientific << std::setprecision(0) << eps << ":";
      o << std::fixed << std::showpos << std::setprecision(1)
        << 100.0*(v/row.fep_eff - 1.0) << "%" << std::noshowpos;
    }
    BOOST_TEST_MESSAGE( o.str() );
  }
}


/** STRUCTURED VALIDATION OF THE TRANSFER MECHANISM ALONE (no volume integration anywhere).
 *
 * Anchors a transfer response on a Monte-Carlo POINT-source curve at 25 cm, then asks it for the
 * efficiency at a grid of other point-source positions - on-axis from contact out to 100 cm, and
 * off-axis out to 45 cm lateral at two axial distances - and compares each against Monte Carlo of
 * the SAME position.  Every number on both sides is a point source, so whatever error shows up is
 * the transfer's, and it is a floor under anything the volumetric model can achieve.
 *
 * Both sides are driven from `query_position`, so the endcap-front-vs-crystal-face reference
 * convention cannot drift between them (commit 9c5e9c5f was exactly that bug: sources placed one
 * dead layer too far away, ~3% at contact).
 *
 * Developer-only: this is a LOT of Monte Carlo - 12 anchor energies at 0.2% plus ~19 positions x 12
 * energies at 0.5%.  Achieved precision is printed per point, because the far and wide-angle points
 * have small efficiencies and may hit the event cap rather than the precision target; a row that
 * did not converge must not be read as a transfer error.
 */
BOOST_AUTO_TEST_CASE( TransferValidationStudy, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();

  const std::vector<double> energies = { 50, 60, 80, 100, 150, 200, 350, 661, 1001, 1460, 2000, 2614 };
  const double anchor_d_cm = 25.0;

  std::vector<std::unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator calc;
  ceelo::ResponseGenerator::configure_calculator( calc, det.gd, owned );

  const auto mc_at = [&]( const Eigen::Vector3d &pos, const double energy,
                          const double target_rel, const size_t max_events )
  {
    calc.set_point_source( pos );
    ceelo::SimulationConfig cfg;
    cfg.energy_keV = energy;
    cfg.termination.target_fep_rel_precision = target_rel;
    cfg.termination.max_events = max_events;
    cfg.termination.min_events = 200000;
    cfg.seed = 20260831;
    const ceelo::EfficiencyResult r = calc.compute( cfg );
    const double eff = r.full_energy_peak_efficiency;
    return std::make_pair( eff, (eff > 0.0) ? (r.fep_uncertainty/eff) : 1.0 );
  };

  // ---- 1) the anchor: MC point source at 25 cm, 0.2% ----
  BOOST_TEST_MESSAGE( "== anchor: MC point source at " << anchor_d_cm << " cm, target 0.2% ==" );
  CeeLoUtils::TransferAnchor anchor;
  anchor.ref_distance_cm = anchor_d_cm;
  anchor.curve_derived = false;
  for( const double e : energies )
  {
    const Eigen::Vector3d pos = det.mc_transfer->query_position( 0.0, 0.0, anchor_d_cm );
    const std::pair<double,double> mc = mc_at( pos, e, 0.002, 400000000 );
    anchor.curve.energies_keV.push_back( e );
    anchor.curve.eff.push_back( mc.first );
    anchor.curve.frac_sigma.push_back( mc.second );
    BOOST_TEST_MESSAGE( "   " << std::setw(6) << e << " keV: eff=" << std::scientific
                        << std::setprecision(6) << mc.first << "  achieved sigma="
                        << std::fixed << std::setprecision(3) << 100.0*mc.second << "%" );
  }

  const std::shared_ptr<ceelo::DetectorResponse> transfer =
        CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{},
                                          "ANGLE GEM35-70 (transfer validation)" );
  BOOST_REQUIRE( transfer );

  // ---- 2) evaluation grid ----
  struct EvalPt { std::string label; double axial_cm; double offset_cm; };
  std::vector<EvalPt> pts;
  for( const double d : { 0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 25.0, 30.0, 50.0, 100.0 } )
    pts.push_back( { (d == 0.0 ? std::string("contact") : std::string("on-axis")), d, 0.0 } );
  for( const double ax : { 5.0, 25.0 } )
    for( const double off : { 5.0, 10.0, 22.0, 45.0 } )   //offset 0 is already on-axis above
      pts.push_back( { "off-axis", ax, off } );

  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "== transfer vs MC, point sources.  (transfer/MC - 1), with MC sigma ==" );

  double worst = 0.0; std::string worst_where;
  for( const EvalPt &p : pts )
  {
    const double dist = std::sqrt( p.axial_cm*p.axial_cm + p.offset_cm*p.offset_cm );
    const double theta = (dist > 0.0) ? std::atan2( p.offset_cm, p.axial_cm ) : 0.0;

    std::ostringstream o;
    o << "  " << std::setw(8) << p.label << " axial=" << std::setw(6) << std::fixed
      << std::setprecision(1) << p.axial_cm << " off=" << std::setw(5) << p.offset_cm
      << " (d=" << std::setw(6) << dist << ", theta=" << std::setw(5) << std::setprecision(1)
      << theta*180.0/M_PI << "deg):";

    for( const double e : energies )
    {
      const Eigen::Vector3d pos = transfer->query_position( theta, 0.0, dist );
      const std::pair<double,double> mc = mc_at( pos, e, 0.005, 100000000 );
      const ceelo::EffResult tr = transfer->eps_fep( e, theta, 0.0, dist );
      const double rel = (mc.first > 0.0) ? 100.0*(tr.value/mc.first - 1.0) : 0.0;
      if( (mc.second < 0.02) && (std::fabs(rel) > worst) )
      {
        worst = std::fabs(rel);
        std::ostringstream w; w << p.label << " ax=" << p.axial_cm << " off=" << p.offset_cm
                                << " @" << e << "keV";
        worst_where = w.str();
      }
      o << " " << std::setw(5) << std::setprecision(0) << e << ":" << std::showpos
        << std::setprecision(1) << std::setw(6) << rel << "%" << std::noshowpos
        << "(" << std::setprecision(1) << 100.0*mc.second << ")";
    }
    BOOST_TEST_MESSAGE( o.str() );
  }

  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "  worst |transfer-MC| over converged points: " << worst << "% at " << worst_where );
}


/** Do InterSpec and CeeLo agree on mu for the SAME material?
 *
 * The transfer validation is CeeLo-on-both-sides, but the VOLUMETRIC comparison is a mix: MC truth
 * attenuates the source through CeeLo cross-sections (set_source_material), while InterSpec's model
 * attenuates through MassAttenuationTool (transmition_length_coefficient).  Any disagreement in mu
 * therefore shows up as a model-vs-truth error that has nothing to do with the geometry model.
 *
 * It compounds exponentially: at 60 keV in steel mu ~ 9.4/cm, so a 1% table difference over a 1 cm
 * path is already ~9% in the efficiency.  That is the right magnitude for the unexplained
 * low-energy volumetric deficit, which is why this is worth checking before anything else.
 *
 * Compared like-for-like: CeeLo's mu_total INCLUDES Rayleigh, InterSpec's excludes it, so the
 * comparison uses CeeLo's (mu_total - mu_rs).
 */
BOOST_AUTO_TEST_CASE( CrossSectionConsistency )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  double worst = 0.0; std::string worst_where;

  for( const char * const name : { "Stainless Steel NIST", "Water", "Iron", "Lead" } )
  {
    std::shared_ptr<const Material> m;
    try { m = matdb->material( name ); }catch( std::exception & ){ }
    if( !m ) continue;

    const ceelo::Material cm = CeeLoUtils::to_ceelo_material( *m ).to_material();
    BOOST_TEST_MESSAGE( "  " << name << " (density " << m->density/(PhysicalUnits::g/PhysicalUnits::cm3)
                        << " g/cm3 InterSpec, " << cm.density() << " CeeLo)" );

    for( const double e : { 50.0, 60.0, 80.0, 100.0, 122.0, 150.0, 200.0, 344.0, 661.7, 1332.5, 2614.0 } )
    {
      // InterSpec: per PhysicalUnits length; convert to 1/cm.
      const double mu_is = GammaInteractionCalc::transmition_length_coefficient( m.get(), (float)e )
                             * PhysicalUnits::cm;
      const ceelo::MacroscopicXS xs = cm.macroscopic_xs( e*1.0e-3 );
      const double mu_ce = xs.mu_total() - xs.mu_rs;     //match InterSpec's no-Rayleigh convention

      const double rel = (mu_ce > 0.0) ? 100.0*(mu_is/mu_ce - 1.0) : 0.0;
      if( std::fabs(rel) > worst )
      {
        worst = std::fabs(rel);
        std::ostringstream w; w << name << " @" << e << " keV";
        worst_where = w.str();
      }
      BOOST_TEST_MESSAGE( "    E=" << std::setw(7) << e << " keV   InterSpec=" << std::fixed
          << std::setprecision(5) << std::setw(10) << mu_is << "/cm   CeeLo=" << std::setw(10) << mu_ce
          << "/cm   diff=" << std::showpos << std::setprecision(2) << rel << "%" << std::noshowpos );
    }
  }

  BOOST_TEST_MESSAGE( "  worst mu disagreement: " << worst << "% (" << worst_where << ")" );
}


/** WHICH FACTOR of the transfer carries the near-field high-energy bias?
 *
 * eps_fep = prefactor(E, pos) * sum_i w_i(E, pos), and the response is GROUNDED at the anchor
 * distance - so it is exact there by construction and the whole error is in how the prediction
 * SCALES away from it.  Comparing ratios-to-anchor therefore isolates the mechanism:
 *
 *     predicted ratio = [prefactor(d)/prefactor(ref)] * [sum_w(d)/sum_w(ref)]
 *     true ratio      = MC(d)/MC(ref)
 *
 * `sum_w` is the ray-traced aperture kernel: it knows the real chord-length distribution through
 * the crystal, which is what matters at high energy where the crystal is semi-transparent
 * (mu*L ~ 1.5 in Ge at 2614 keV, so p_int ~ 0.78 and the efficiency is chord-sensitive). The
 * prefactor carries everything else - the angular shape and the near-field correction.
 *
 * Whichever of the two ratios departs from the MC ratio is where the bias lives.
 */
BOOST_AUTO_TEST_CASE( TransferFactorDecomposition, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const double ref_cm = kMcAnchorDistanceCm;
  std::vector<std::unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator calc;
  ceelo::ResponseGenerator::configure_calculator( calc, det.gd, owned );

  const auto mc_at = [&]( const Eigen::Vector3d &pos, const double e )
  {
    calc.set_point_source( pos );
    ceelo::SimulationConfig cfg;
    cfg.energy_keV = e;
    cfg.termination.target_fep_rel_precision = 0.003;
    cfg.termination.max_events = 200000000;
    cfg.termination.min_events = 200000;
    cfg.seed = 20260831;
    return calc.compute( cfg ).full_energy_peak_efficiency;
  };

  const auto decompose = [&]( const double d, const double e )
  {
    const Eigen::Vector3d pos = det.mc_transfer->query_position( 0.0, 0.0, d );
    const ceelo::ApertureQuadrature q =
          ceelo::make_aperture_quadrature( det.mc_transfer->geometry(), pos, 8192 );
    std::vector<double> w; std::vector<Eigen::Vector3d> dirs;
    det.mc_transfer->fep_ray_weights( e, q, w, dirs );
    double sum_w = 0.0;
    for( const double x : w ) sum_w += x;
    const ceelo::EffResult pre = det.mc_transfer->fep_prefactor( e, pos, q );
    return std::make_tuple( pre.value, sum_w, pre.value*sum_w );
  };

  for( const double e : { 60.0, 344.0, 1332.5, 2614.0 } )
  {
    const std::tuple<double,double,double> R = decompose( ref_cm, e );
    const double mc_ref = mc_at( det.mc_transfer->query_position(0.0,0.0,ref_cm), e );
    BOOST_TEST_MESSAGE( "" );
    BOOST_TEST_MESSAGE( "  E=" << e << " keV   (ratios to the " << ref_cm << " cm anchor)" );
    BOOST_TEST_MESSAGE( "      d(cm)   prefactor_ratio   sumw_ratio   predicted   MC_true    error" );

    for( const double d : { 0.0, 1.0, 3.0, 5.0, 10.0, 25.0 } )
    {
      const std::tuple<double,double,double> X = decompose( d, e );
      const double mc = mc_at( det.mc_transfer->query_position(0.0,0.0,d), e );

      const double pre_r  = std::get<0>(X)/std::get<0>(R);
      const double sumw_r = std::get<1>(X)/std::get<1>(R);
      const double pred_r = std::get<2>(X)/std::get<2>(R);
      const double true_r = mc/mc_ref;

      std::ostringstream o;
      o << "    " << std::setw(7) << std::fixed << std::setprecision(1) << d
        << std::setw(16) << std::setprecision(4) << pre_r
        << std::setw(13) << sumw_r
        << std::setw(12) << pred_r
        << std::setw(10) << true_r
        << std::setw(9) << std::showpos << std::setprecision(1)
        << 100.0*(pred_r/true_r - 1.0) << "%";
      BOOST_TEST_MESSAGE( o.str() );
    }
  }
}


/** CONFIRMATION: the near-field high-energy bias is peak-to-total, not solid angle.
 *
 * `TransferFactorDecomposition` localises the entire bias to the aperture kernel
 * K = sum_i omega_w * p_int, where p_int is the probability of INTERACTING somewhere in the active
 * chord. At low energy, interacting means depositing everything (photoelectric, short mfp), so K is
 * a good FEP kernel. At high energy it is not: a first Compton interaction deposits only part of
 * the energy, and whether the rest is captured depends on how much crystal remains along the track.
 *
 * That residual path length is geometry-dependent. At 25 cm the rays are near-parallel to the axis
 * and see the full crystal length; at contact they fan out to steep angles and clip corners, so
 * more first interactions escape. The kernel models the reduced interaction probability but not the
 * reduced capture-given-interaction, and the grounding factor can only absorb the average at the
 * anchor distance - so the discrepancy does not cancel in the ratio.
 *
 * If that is right, the MC peak-to-total ratio must FALL from 25 cm to contact, by about the same
 * amount as the kernel over-predicts, and must be flat at 60 keV where it is ~1 either way.
 */
BOOST_AUTO_TEST_CASE( PeakToTotalVsDistance, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  std::vector<std::unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator calc;
  ceelo::ResponseGenerator::configure_calculator( calc, det.gd, owned );

  for( const double e : { 60.0, 344.0, 1332.5, 2614.0 } )
  {
    BOOST_TEST_MESSAGE( "" );
    BOOST_TEST_MESSAGE( "  E=" << e << " keV" );
    BOOST_TEST_MESSAGE( "     d(cm)   peak-to-total   vs 25cm    mean_chord(cm)  mean|cos|" );
    double pt_ref = -1.0;
    for( const double d : { 0.0, 1.0, 3.0, 5.0, 10.0, 25.0 } )
    {
      const Eigen::Vector3d pos = det.mc_transfer->query_position( 0.0, 0.0, d );
      calc.set_point_source( pos );
      ceelo::SimulationConfig cfg;
      cfg.energy_keV = e;
      cfg.termination.target_fep_rel_precision = 0.003;
      cfg.termination.max_events = 200000000;
      cfg.termination.min_events = 200000;
      cfg.seed = 20260831;
      const ceelo::EfficiencyResult r = calc.compute( cfg );
      const double pt = (r.total_efficiency > 0.0)
                          ? (r.full_energy_peak_efficiency/r.total_efficiency) : 0.0;

      const ceelo::ApertureQuadrature q =
            ceelo::make_aperture_quadrature( det.mc_transfer->geometry(), pos, 8192 );

      if( d >= 25.0 ) pt_ref = pt;   //filled on the last iteration; ratios printed relative below
      std::ostringstream o;
      o << "   " << std::setw(7) << std::fixed << std::setprecision(1) << d
        << std::setw(14) << std::setprecision(4) << pt
        << std::setw(11) << "-"
        << std::setw(16) << std::setprecision(3) << q.mean_chord_cm
        << std::setw(11) << std::setprecision(4) << q.mean_cos_incidence;
      BOOST_TEST_MESSAGE( o.str() );
    }
    BOOST_TEST_MESSAGE( "     (25 cm peak-to-total = " << std::setprecision(4) << pt_ref << ")" );
  }
}


/** Is CeeLo already TELLING us the near-field query is unmodelled?
 *
 * eps_fep raises ResponseFlag::NearFieldUnmodeled and adds 5% fractional sigma whenever the query
 * distance is inside the near-field gate and the response carries no NearFieldModel - which is the
 * case for every curve-anchored transfer, since neither CeeLoUtils nor EfficiencyTransfer ever
 * populates one.  If so, the near-field bias is a DECLARED limitation we are ignoring, not a
 * silent error, and the first fix is to surface the flag and the sigma (plan 3.5).
 */
BOOST_AUTO_TEST_CASE( NearFieldFlagIsRaised )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  BOOST_TEST_MESSAGE( "    d(cm)      eps_fep     frac_sigma   flag" );
  for( const double d : { 0.0, 1.0, 3.0, 5.0, 10.0, 20.0, 25.0, 50.0 } )
  {
    const ceelo::EffResult r = det.mc_transfer->eps_fep( 1332.5, 0.0, 0.0, d );
    const double frac = (r.value > 0.0) ? (r.sigma/r.value) : 0.0;
    BOOST_TEST_MESSAGE( "  " << std::setw(7) << std::fixed << std::setprecision(1) << d
        << std::setw(14) << std::scientific << std::setprecision(4) << r.value
        << std::setw(13) << std::fixed << std::setprecision(2) << 100.0*frac << "%   "
        << ceelo::to_string(r.flag) );
  }
}


/** Can a CHORD-DEPENDENT peak fraction explain the near-field bias?  (design check for the fix)
 *
 * The kernel weights each ray by p_int = 1 - exp(-mu*L): probability of interacting in an active
 * chord L.  CeeLo then factors eps_fep = k(E)*eta(E,theta)*N*K, which assumes the ratio
 * (peak fraction / interaction probability) is the SAME for every ray, so it can be absorbed into
 * eta and the grounding k.  Far-field that is fine - every ray is near-parallel and sees nearly the
 * same chord.  At contact the rays fan out and the chord distribution changes completely, so a
 * single eta cannot represent it.
 *
 * This prints the aperture's chord distribution at each distance.  If the mean and spread move
 * enough to account for the measured peak-to-total drop, then re-weighting each ray by a
 * chord-dependent peak fraction is the right fix; if the distribution barely moves, the bias is
 * something else and the fix would be misdirected.
 */
BOOST_AUTO_TEST_CASE( ApertureChordDistribution )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const ceelo::Geometry &geom = det.mc_transfer->geometry();

  BOOST_TEST_MESSAGE( "    d(cm)   mean_chord  median   p10    p90   frac(L<2cm)  omega_active" );
  for( const double d : { 0.0, 1.0, 3.0, 5.0, 10.0, 25.0 } )
  {
    const Eigen::Vector3d pos = det.mc_transfer->query_position( 0.0, 0.0, d );
    const ceelo::ApertureQuadrature q = ceelo::make_aperture_quadrature( geom, pos, 8192 );

    std::vector<double> chords;
    double wsum = 0.0, short_w = 0.0;
    for( const ceelo::KernelRay &r : q.rays )
    {
      if( r.active_len <= 0.0f ) continue;
      chords.push_back( r.active_len );
      wsum += r.omega_w;
      if( r.active_len < 2.0f ) short_w += r.omega_w;
    }
    std::sort( chords.begin(), chords.end() );
    const auto pct = [&]( const double f ){
      return chords.empty() ? 0.0 : chords[ std::min(chords.size()-1,
                                    size_t(f*(chords.size()-1))) ]; };

    BOOST_TEST_MESSAGE( "  " << std::setw(7) << std::fixed << std::setprecision(1) << d
        << std::setw(12) << std::setprecision(3) << q.mean_chord_cm
        << std::setw(9) << pct(0.5) << std::setw(7) << pct(0.1) << std::setw(7) << pct(0.9)
        << std::setw(12) << std::setprecision(1) << 100.0*(wsum>0.0?short_w/wsum:0.0) << "%"
        << std::setw(14) << std::scientific << std::setprecision(3) << q.omega_frac_active );
  }
}


/** What is the covariance path's sigma ACTUALLY, versus the per-query sigma?
 *
 * frac_covariance() is not just the floor: it is
 *     floor^2 + transfer_i*transfer_j + cov_ln_k + node_i^2
 * so it already carries the grounding TRANSFER term, which grows as the query geometry departs
 * from the anchor.  Only the NearFieldUnmodeled penalty (a hard-coded 0.05 in quadrature) is
 * missing from it.  Measure both so the size of the inflation is a number, not an inference.
 */
BOOST_AUTO_TEST_CASE( CovarianceVsPerQuerySigma )
{
  set_data_dir();
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const std::vector<double> es = { 60.0, 1332.5 };
  BOOST_TEST_MESSAGE( "    d(cm)     E(keV)   frac_cov sigma   eps_fep sigma   inflation   flag" );
  for( const double d : { 0.0, 1.0, 5.0, 10.0, 25.0, 50.0 } )
  {
    const std::vector<double> cov = det.mc_transfer->frac_covariance( es, 0.0, d );
    for( size_t i = 0; i < es.size(); ++i )
    {
      const double cov_sig = std::sqrt( std::max(0.0, cov[i*es.size()+i]) );
      const ceelo::EffResult r = det.mc_transfer->eps_fep( es[i], 0.0, 0.0, d );
      const double q_sig = (r.value > 0.0) ? (r.sigma/r.value) : 0.0;
      BOOST_TEST_MESSAGE( "  " << std::setw(7) << std::fixed << std::setprecision(1) << d
          << std::setw(10) << es[i]
          << std::setw(15) << std::setprecision(2) << 100.0*cov_sig << "%"
          << std::setw(15) << 100.0*q_sig << "%"
          << std::setw(11) << std::setprecision(2) << (cov_sig>0 ? q_sig/cov_sig : 0.0) << "x   "
          << ceelo::to_string(r.flag) );
    }
  }
}
