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

#pragma once

/** The test harness shared by test_VolumetricNearField.cpp and test_VolumetricLadder.cpp: the ANGLE
 GEM35-70 fixture, the InterSpec-side scenario builder/integrator, the CeeLo-side scenario
 configuration, and (at the bottom) an on-disk Monte-Carlo result cache.

 TEST-ONLY.  Everything lives in an anonymous namespace on purpose: each test executable is a single
 translation unit, and the helpers use Boost.Test macros, so the including TU must have defined
 BOOST_TEST_MODULE and included <boost/test/included/unit_test.hpp> BEFORE this header.

 Extracted verbatim from test_VolumetricNearField.cpp (2026-09-02) so a second test TU can drive the
 SAME integrand and the SAME Monte-Carlo configuration as the committed truth-bank comparison.
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <functional>

#include <boost/test/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

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
string g_data_dir, g_test_data_dir, g_cache_dir;

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
    else if( SpecUtils::starts_with( arg, "--cachedir=" ) )
      g_cache_dir = arg.substr( 11 );
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
              + mu_Rayleigh(mat, E) * h(E, mat, mu * normal_thickness_cm)

 `win_keV` is a HALF-width, matching CeeLo's +-win convention - which is why its 0.75 keV default
 is about FWHM/2 for an HPGe.

 The LAST term is the Rayleigh DEFLECTION loss (2026-09-04).  "Elastic so it cannot leave the
 window" is right about energy and wrong about geometry: in Fe at 60 keV half of the coherent
 scatters exceed 20 degrees, and in a THICK layer that deflection lengthens the photon's remaining
 path, so a fraction h of the Rayleigh-scattered photons are absorbed after all.  h is
 ceelo::rayleigh_deflection_loss_fraction, evaluated at the layer's NORMAL non-Rayleigh optical
 depth; it vanishes for a thin layer and reaches 0.25 behind 0.5 cm Fe at 60 keV, where the term is
 worth 9% of the peak - the whole of the deep-shield over-read.  `normal_thickness_cm` <= 0 (the
 default, and what the SOURCE shell passes) turns it off, because a self-attenuating matrix
 self-selects a ~1 mfp skin and its full extent is the wrong depth.
 */
double fep_removal_coefficient( const Material &mat, const double energy_keV, const double win_keV,
                                const double normal_thickness_cm = 0.0 )
{
  const double mu_total = GammaInteractionCalc::transmition_length_coefficient(
                                                    &mat, static_cast<float>(energy_keV) );

  const auto process_mu = [&]( const MassAttenuation::GammaEmProcces process ) -> double {
    double mu = 0.0;
    for( const Material::ElementFractionPair &p : mat.elements )
      mu += p.second * mat.density
            * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber,
                                                static_cast<float>(energy_keV), process );
    for( const Material::NuclideFractionPair &p : mat.nuclides )
      mu += p.second * mat.density
            * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber,
                                                static_cast<float>(energy_keV), process );
    return mu;
  };
  const double mu_compton = process_mu( MassAttenuation::GammaEmProcces::ComptonScatter );

  // Material-aware in-window fraction: S(x,Z) suppresses forward scatter, most strongly at high Z
  //  (@60 keV: water 0.89, Fe 0.77, Pb 0.63), so the free-electron form over-credits exactly where
  //  the shielding is heaviest.  This is a Simpson integration - hoisted per (E, material) by the
  //  caller, never inside a per-ray loop.
  const ceelo::Material cm = CeeLoUtils::to_ceelo_material( mat ).to_material();
  const double f_win = (win_keV > 0.0) ? ceelo::kn_in_window_fraction( energy_keV, win_keV, cm ) : 0.0;

  double rayleigh_loss = 0.0;
  if( normal_thickness_cm > 0.0 )
  {
    const double mu_rayleigh = process_mu( MassAttenuation::GammaEmProcces::RayleighScatter );
    const double tau_nr = mu_total * normal_thickness_cm * PhysicalUnits::cm;
    rayleigh_loss = mu_rayleigh * ceelo::rayleigh_deflection_loss_fraction( energy_keV, tau_nr, cm );
  }

  return mu_total - f_win*mu_compton + rayleigh_loss;
}


/** Builds the DistributedSrcCalcT a scenario maps onto - shells, detector geometry, response and
 ray count - WITHOUT integrating it.

 Factored out of interspec_volumetric_eff so a diagnostic can integrate the very same integrand by
 an independent rule (see BoxOuterQuadratureVsTensorGL) instead of keeping a second copy of this
 setup that could silently drift from the one under test.  Callers that integrate by hand must call
 finalize_shell_coefficients() first: ShellInfo::fep_trans_len_coef defaults to a NEGATIVE sentinel
 and self_shielding_integration_imp is normally what replaces it.
 */
GammaInteractionCalc::DistributedSrcCalcT<double>
build_scenario_calc( const AngleDetector &det,
                     const Scenario &s,
                     const double energy_keV,
                     const shared_ptr<const ceelo::DetectorResponse> &eff_response,
                     const int n_rays = -1,
                     const double fep_window_keV = -1.0,
                     const bool transparent = false )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;

  // Source matrix, as an InterSpec material, so the model attenuates with InterSpec's own tables -
  //  the SAME material the MC truth was generated through (see scenario_matrix_material).
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> matrix = matdb->material( scenario_matrix_material( s.dense ) );
  BOOST_REQUIRE( matrix );

  const GeometryType geom = (s.shape == Shape::Box) ? GeometryType::Rectangular
                          : ((s.shape == Shape::CylinderSideOn) ? GeometryType::CylinderSideOn
                                                                : GeometryType::CylinderEndOn);

  DistributedSrcCalcT<double> calc;
  calc.m_geometry = geom;
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
  // The lateral offset goes in as the detector-side offset; CeeLo displaces the SOURCE by the same
  //  magnitude (configure_scenario_mc, via scenario_center) - mirror images of one geometry.
  calc.m_detector = detector_geom_from_config<double>( geom,
                                                       scenario_center_distance_cm( s ) * cm,
                                                       det_radius, 0.0, s.offset_cm * cm, 0.0 );

  DistributedSrcCalcT<double>::ShellInfo src_shell;
  src_shell.dims = (s.shape == Shape::Box)
        ? std::array<double,3>{ s.half_width_cm*cm, s.half_height_cm*cm, s.half_length_cm*cm }
        : std::array<double,3>{ s.radius_cm*cm, s.half_length_cm*cm, 0.0 };
  // `transparent` keeps the shell GEOMETRY but removes its attenuation, so the multi-shell walk
  //  runs exactly the path the attenuating case does - see TransparentTransferFloor.
  src_shell.trans_len_coef = transparent
        ? 0.0
        : ( (fep_window_keV > 0.0)
              ? fep_removal_coefficient( *matrix, energy_keV, fep_window_keV )
              : transmition_length_coefficient( matrix.get(), static_cast<float>(energy_keV) ) );
  src_shell.type = ShellType::Material;
  calc.m_shells.push_back( src_shell );

  if( s.shield_cm > 0.0 )
  {
    const shared_ptr<const Material> iron = matdb->material( scenario_shield_material() );
    BOOST_REQUIRE( iron );
    DistributedSrcCalcT<double>::ShellInfo shield;
    shield.dims = (s.shape == Shape::Box)
          ? std::array<double,3>{ (s.half_width_cm + s.shield_cm)*cm,
                                  (s.half_height_cm + s.shield_cm)*cm,
                                  (s.half_length_cm + s.shield_cm)*cm }
          : std::array<double,3>{ (s.radius_cm + s.shield_cm)*cm,
                                  (s.half_length_cm + s.shield_cm)*cm, 0.0 };
    // The Rayleigh deflection loss applies to the SHIELD only (see fep_removal_coefficient);
    //  it rides along with the window credit so a `fep_window_keV <= 0` row stays plain mu_total,
    //  which keeps the historical "mu_total" columns of the ladder meaning what they say.
    shield.trans_len_coef = transparent
          ? 0.0
          : ( (fep_window_keV > 0.0)
                ? fep_removal_coefficient( *iron, energy_keV, fep_window_keV, s.shield_cm )
                : transmition_length_coefficient( iron.get(), static_cast<float>(energy_keV) ) );
    shield.type = ShellType::Material;
    calc.m_shells.push_back( shield );
  }//if( shielded )

  return calc;
}//build_scenario_calc(...)


/** The CENTRE-ANCHORED transfer for a scenario: an EFFTRAN response grounded on the recorded
 point-source curve at that scenario's own source-centre distance (`sm_centre_anchor`), so the
 transfer is not extrapolating in distance at all.

 Built from RECORDED Monte Carlo, so this runs no simulation and is safe in a committed test.  Cached
 by distance: a box and its cylinder twin share one response exactly, as does every scenario whose
 centre happens to sit at the same place.  Returns null if the table has no curve for that distance,
 which is what a caller should treat as "regenerate Rung7_CentreAnchorTruth".
 */
shared_ptr<const ceelo::DetectorResponse> centre_anchored_response( const AngleDetector &det,
                                                                    const Scenario &s )
{
  static map<long long,shared_ptr<const ceelo::DetectorResponse>> cache;

  const double d_cm = scenario_center_distance_cm( s );
  const long long key = llround( d_cm * 1.0e6 );
  const map<long long,shared_ptr<const ceelo::DetectorResponse>>::const_iterator it = cache.find( key );
  if( it != end(cache) )
    return it->second;

  CeeLoUtils::TransferAnchor anchor;
  anchor.ref_distance_cm = d_cm;
  anchor.curve_derived = false;
  for( const CentreAnchorRow &row : sm_centre_anchor )
  {
    if( llround( row.distance_cm * 1.0e6 ) != key )
      continue;
    anchor.curve.energies_keV.push_back( row.energy_keV );
    anchor.curve.eff.push_back( row.eff );
    anchor.curve.frac_sigma.push_back( row.frac_sigma );
  }

  shared_ptr<const ceelo::DetectorResponse> resp;
  if( anchor.curve.energies_keV.size() >= 2 )
  {
    ostringstream nm;
    nm << "centre-anchored (d=" << fixed << setprecision(2) << d_cm << " cm)";
    resp = CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{}, nm.str() );
  }

  cache[key] = resp;
  return resp;
}//centre_anchored_response(...)


/** Attaches a detector-side line set (the LINE integration path) to a calculator built by
 #build_scenario_calc, from the calculator's own scalar source dims and detector placement - the
 same construction ShieldingSourceChi2Fcn::volumetricLineCache does inside a fit.
 */
void attach_line_cache( GammaInteractionCalc::DistributedSrcCalcT<double> &calc, const int num_lines )
{
  using namespace GammaInteractionCalc;
  BOOST_REQUIRE( calc.m_effResponse );
  const std::array<double,3> &src = calc.m_shells[calc.m_materialIndex].dims;
  const std::array<double,3> det_pos = { calc.m_detector.position[0], calc.m_detector.position[1],
                                         calc.m_detector.position[2] };
  const std::array<double,3> det_axis = { calc.m_detector.axis[0], calc.m_detector.axis[1],
                                          calc.m_detector.axis[2] };
  calc.m_lineCache = build_volumetric_line_cache( calc.m_effResponse, calc.m_geometry,
                                                  calc.m_materialIndex, src, det_pos, det_axis, 0.0,
                                                  num_lines );
}//attach_line_cache(...)


/** Integrates a calculator on the requested path (Element = the per-element aperture reference;
 Line = the detector-side line set with `num_lines` lines), through the production dispatcher. */
void integrate_on_path( GammaInteractionCalc::DistributedSrcCalcT<double> &calc,
                        const GammaInteractionCalc::VolumetricIntegrator path,
                        const int num_lines = -1 )
{
  using namespace GammaInteractionCalc;
  if( path == VolumetricIntegrator::Line )
    attach_line_cache( calc, (num_lines > 0) ? num_lines : 65536 );
  std::vector<std::unique_ptr<DistributedSrcCalcT<double>>> calcs;
  calcs.push_back( std::make_unique<DistributedSrcCalcT<double>>( calc ) );
  {
    const ScopedVolumetricIntegratorOverride force( path );
    integrate_volumetric_calculators<double>( calcs, true );
  }
  calc.integral = calcs.front()->integral;
  calc.m_num_evals = calcs.front()->m_num_evals;
  calc.m_est_rel_error = calcs.front()->m_est_rel_error;
}//integrate_on_path(...)


/** Harness-wide override: when > 0 (INTERSPEC_HARNESS_LINE_COUNT), every #interspec_volumetric_eff
 call that has a response and was not given an explicit line count runs on the LINE path with this
 many lines, so any existing rung can be re-run against the line quadrature without being rewritten.
 -1 (unset) keeps the element path. */
int harness_line_count()
{
  static const int value = [](){
    const char *env = std::getenv( "INTERSPEC_HARNESS_LINE_COUNT" );
    const int n = env ? std::atoi( env ) : -1;
    if( n > 0 )
      BOOST_TEST_MESSAGE( "  (INTERSPEC_HARNESS_LINE_COUNT=" << n << ": volumetric integrals run on"
                          " the LINE path)" );
    return n;
  }();
  return value;
}


double interspec_volumetric_eff( const AngleDetector &det,
                                 const Scenario &s,
                                 const double energy_keV,
                                 const shared_ptr<const ceelo::DetectorResponse> &eff_response,
                                 const int n_rays = -1,
                                 const double fep_window_keV = -1.0,
                                 const double epsrel = -1.0,
                                 const bool transparent = false,
                                 size_t *num_evals_out = nullptr,
                                 double *est_rel_error_out = nullptr,
                                 const int num_lines = -1 )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;

  DistributedSrcCalcT<double> calc = build_scenario_calc( det, s, energy_keV, eff_response,
                                                          n_rays, fep_window_keV, transparent );

  const int use_lines = (num_lines > 0) ? num_lines
                        : ((eff_response && (epsrel <= 0.0)) ? harness_line_count() : -1);

  if( use_lines > 0 )
  {
    // The LINE path (VolumetricLineIntegration_imp.hpp); epsrel does not apply.
    BOOST_REQUIRE( eff_response );
    integrate_on_path( calc, VolumetricIntegrator::Line, use_lines );
  }else if( epsrel > 0.0 )
  {
    self_shielding_integration_imp<double>( calc, epsrel, 200000000 );
  }else
  {
    self_shielding_integration_imp<double>( calc );
  }

  if( num_evals_out )
    *num_evals_out = calc.m_num_evals;
  if( est_rel_error_out )
    *est_rel_error_out = calc.m_est_rel_error;

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


/** The Scenario named `name`, BY VALUE.

 `scenarios()` returns by value, so a pointer or reference into a range-for over it dangles the
 moment the loop ends - which produced a NaN once.  Copying out is the fix, and it is cheap.
 */
Scenario find_scenario( const string &name )
{
  const vector<Scenario> all = scenarios();
  for( const Scenario &s : all )
  {
    if( s.name == name )
      return s;
  }

  BOOST_REQUIRE_MESSAGE( false, "No scenario named '" << name << "'" );
  return Scenario{};
}//find_scenario(...)


/** An EFFTRAN transfer anchored on a Monte-Carlo POINT source sitting on axis at
 `centre_distance_cm` - i.e. at the CENTRE of a scenario whose source centre is that far out.

 The anchor is a bare point in vacuum: no source matrix and no shield.  It therefore depends only on
 the distance and the energy grid, and NOT on the shape, matrix or shield of the scenario whose
 centre that distance happens to be - so a box and its cylinder twin share one anchor exactly.  That
 is why the callers cache these by distance rather than by scenario.

 Minutes of Monte Carlo per call (one run per `scenario_energies()` entry at 0.5%).
 */
shared_ptr<ceelo::DetectorResponse> make_centre_anchor_response( const AngleDetector &det,
                                                                 ceelo::EfficiencyCalculator &pt,
                                                                 const double centre_distance_cm )
{
  CeeLoUtils::TransferAnchor anchor;
  anchor.ref_distance_cm = centre_distance_cm;
  anchor.curve_derived = false;

  // Matches scenario_center(): the crystal-face frame puts the source in front, at negative z.
  const Eigen::Vector3d centre( 0.0, 0.0, -(det.endcap_front_offset_cm + centre_distance_cm) );

  for( const double e : scenario_energies() )
  {
    pt.set_point_source( centre );
    ceelo::SimulationConfig cfg;
    cfg.energy_keV = e;
    cfg.termination.target_fep_rel_precision = 0.005;
    cfg.termination.max_events = 400000000;
    cfg.termination.min_events = 100000;
    cfg.seed = 20260901;

    const ceelo::EfficiencyResult r = pt.compute( cfg );
    anchor.curve.energies_keV.push_back( e );
    anchor.curve.eff.push_back( r.full_energy_peak_efficiency );
    anchor.curve.frac_sigma.push_back( (r.full_energy_peak_efficiency > 0.0)
                  ? (r.fep_uncertainty/r.full_energy_peak_efficiency) : 0.0 );
  }//for( anchor energies )

  ostringstream nm;
  nm << "near-anchored (centre, d=" << fixed << setprecision(2) << centre_distance_cm << " cm)";

  const shared_ptr<ceelo::DetectorResponse> resp =
        CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{}, nm.str() );
  BOOST_REQUIRE( resp );

  return resp;
}//make_centre_anchor_response(...)


/** One-line description of a scenario's source: shape, outer dimensions, standoff, shield. */
string scenario_description( const Scenario &s )
{
  ostringstream o;
  o << fixed << setprecision(1);
  if( s.shape == Shape::Box )
    o << "box " << 2.0*s.half_width_cm << " x " << 2.0*s.half_height_cm
      << " x " << 2.0*s.half_length_cm << " cm";
  else
    o << ((s.shape == Shape::CylinderSideOn) ? "side-on cylinder r=" : "cylinder r=")
      << s.radius_cm << " cm, len=" << 2.0*s.half_length_cm << " cm";

  o << ", standoff " << s.standoff_cm << " cm";
  if( s.offset_cm != 0.0 )
    o << ", " << s.offset_cm << " cm off axis";
  o << ", " << scenario_matrix_material( s.dense );
  if( s.shield_cm > 0.0 )
    o << ", " << s.shield_cm << " cm " << scenario_shield_material() << " shield";

  return o.str();
}//scenario_description(...)

}//namespace


/** Optical depth (mean free paths) of a scenario's EXTERNAL shield alone, at `energy_keV`; 0 for a
 bare scenario.

 Separate from #scenario_optical_depth because the two depths limit the model in different ways.  A
 deep SOURCE MATRIX self-selects: the escaping signal comes from a surface skin, so getting mu
 slightly wrong there barely moves the integral.  A SHIELD has no such escape - every ray crosses its
 full thickness - so an error in the removal coefficient compounds over the whole of it.  That is why
 the truth-bank comparison groups on this number rather than on the total.
 */
double scenario_shield_optical_depth( const Scenario &s, const double energy_keV )
{
  if( !(s.shield_cm > 0.0) )
    return 0.0;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const shared_ptr<const Material> shield = matdb ? matdb->material( scenario_shield_material() )
                                                  : nullptr;
  if( !shield )
    return 0.0;

  return GammaInteractionCalc::transmition_length_coefficient( shield.get(),
                                                    static_cast<float>(energy_keV) )
         * s.shield_cm * PhysicalUnits::cm;
}//scenario_shield_optical_depth(...)


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
  if( !matrix )
    return 0.0;

  const float e = static_cast<float>( energy_keV );

  // The escape path for an END-ON detector is axial, not radial: a photon heads toward the detector
  //  face, so what it must cross is at most the source's full length.  Using the radius instead
  //  wildly overstates a pancake - the wide-angle-far scenario (r=12 cm, half-length 0.5 cm) would
  //  read ~101 mfp while the model in fact tracks truth to a few percent, because no photon ever
  //  crosses 12 cm of material on its way out.  Take the smaller of the two half-extents as the
  //  characteristic depth: it is the one that bounds the escape path.
  // For a SIDE-ON cylinder the escape path toward the detector is RADIAL, so the radius is the
  //  characteristic depth outright - taking the min would understate a long thin pipe seen sideways.
  const double half_extent_cm = (s.shape == Shape::Box)
        ? std::min( std::min( s.half_width_cm, s.half_height_cm ), s.half_length_cm )
        : ((s.shape == Shape::CylinderSideOn) ? s.radius_cm
                                              : std::min( s.radius_cm, s.half_length_cm ));
  const double tau_src = GammaInteractionCalc::transmition_length_coefficient( matrix.get(), e )
                         * half_extent_cm * PhysicalUnits::cm;

  return tau_src + scenario_shield_optical_depth( s, energy_keV );
}//scenario_optical_depth(...)


/** Runs the MC for each of `scenarios_to_run` at every energy, printing paste-ready TruthRow
 initializers.  Shared by the full and the box-only generators so the two cannot drift apart in how
 they configure the source. */
/** The CeeLo materials a scenario needs, kept alive for as long as the calculator runs. */
struct ScenarioMcMaterials
{
  vector<unique_ptr<ceelo::Material>> owned;   ///< detector materials (configure_calculator)
  // ceelo::Material has no default constructor, so the source materials are heap-held; the
  //  calculator keeps raw pointers to them, so this object must outlive every compute() call.
  unique_ptr<ceelo::Material> matrix;
  unique_ptr<ceelo::Material> shield;
};

/** Configures `calc` as the extended-source Monte-Carlo for scenario `s`: detector geometry, the
 source shape at its scenario centre, the source matrix, and the shield.  `transparent` skips the
 matrix and the shield (an emitting volume in vacuum - CeeLo's null source material is supported).
 Shared by the truth generator and the ladder so the two cannot drift apart. */
void configure_scenario_mc( ceelo::EfficiencyCalculator &calc, const AngleDetector &det,
                            const Scenario &s, ScenarioMcMaterials &mats,
                            const bool transparent = false )
{
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  ceelo::ResponseGenerator::configure_calculator( calc, det.gd, mats.owned );

  // Shape must be configured BEFORE add_source_shield: the 3-thickness overload asserts the
  //  source is rectangular.
  if( s.shape == Shape::Box )
    calc.set_rectangular_source( scenario_center( s, det.endcap_front_offset_cm ),
                                 scenario_box_half_dims( s ) );
  else
    calc.set_cylindrical_source( scenario_center( s, det.endcap_front_offset_cm ),
                                 s.radius_cm, s.half_length_cm,
                                 scenario_source_rotation( s ) );

  if( transparent )
    return;

  const shared_ptr<const Material> matrix_mat = matdb->material( scenario_matrix_material(s.dense) );
  BOOST_REQUIRE_MESSAGE( matrix_mat, "No material '" << scenario_matrix_material(s.dense) << "'" );
  mats.matrix.reset( new ceelo::Material( CeeLoUtils::to_ceelo_material( *matrix_mat ).to_material() ) );
  calc.set_source_material( mats.matrix.get() );

  const shared_ptr<const Material> shield_mat = matdb->material( scenario_shield_material() );
  BOOST_REQUIRE( shield_mat );
  mats.shield.reset( new ceelo::Material( CeeLoUtils::to_ceelo_material( *shield_mat ).to_material() ) );
  if( s.shield_cm > 0.0 )
  {
    if( s.shape == Shape::Box )
      calc.add_source_shield( mats.shield.get(), s.shield_cm, s.shield_cm, s.shield_cm );
    else
      // (radial, end) thicknesses in the source's OWN frame, so this is orientation-independent.
      calc.add_source_shield( mats.shield.get(), s.shield_cm, s.shield_cm );
  }//if( shielded )
}//configure_scenario_mc(...)


/** A textual key describing a scenario's MC configuration, for the on-disk cache. */
string scenario_mc_key( const Scenario &s, const bool transparent = false )
{
  ostringstream o;
  o << setprecision(10);
  if( s.shape == Shape::Box )
    o << "box(hx=" << s.half_width_cm << ",hy=" << s.half_height_cm << ",hz=" << s.half_length_cm << ")";
  else
    o << ((s.shape == Shape::CylinderSideOn) ? "cylSideOn(r=" : "cyl(r=") << s.radius_cm
      << ",hl=" << s.half_length_cm << ")";
  o << ";standoff=" << s.standoff_cm;
  if( s.offset_cm != 0.0 )
    o << ";offset=" << s.offset_cm;
  if( transparent )
    o << ";transparent";
  else
    o << ";matrix=" << scenario_matrix_material( s.dense ) << ";shield=" << s.shield_cm
      << "cm:" << scenario_shield_material();
  return o.str();
}//scenario_mc_key(...)


void print_truth_rows( const AngleDetector &det, const vector<double> &energies,
                       const vector<Scenario> &scenarios_to_run )
{
  for( const Scenario &s : scenarios_to_run )
  {
    ScenarioMcMaterials mats;
    ceelo::EfficiencyCalculator calc;
    configure_scenario_mc( calc, det, s, mats );

    for( const double energy : energies )
    {
      ceelo::SimulationConfig cfg;
      cfg.energy_keV = energy;
      // 0.5% target.  The cap must be large enough that LOW-efficiency rows reach the TARGET
      //  rather than stopping at the CAP: a shielded dense box at 60 keV has eff ~ 1.8e-5, so
      //  0.5% (40k FEP counts) needs ~2.3e9 histories.  A row that stops on the cap carries an
      //  uncertainty larger than the model error measured against it, and reads as model
      //  disagreement rather than as a truth-row limitation.
      cfg.termination.target_fep_rel_precision = 0.005;
      cfg.termination.max_events = 4000000000ULL;
      cfg.termination.min_events = 50000;
      cfg.seed = 20260829;

      const ceelo::EfficiencyResult r = calc.compute( cfg );

      cout << "  { \"" << s.name << "\", " << setprecision(12) << energy << ", "
           << setprecision(10) << r.full_energy_peak_efficiency << ", "
           << setprecision(4) << r.fep_uncertainty << " },\n";
      cout.flush();
    }//for( energies )
  }//for( scenarios_to_run )
}//print_truth_rows(...)


namespace
{
/** A tag identifying the Monte-Carlo PHYSICS that produced a cached result: the vendored CeeLo
 baseline plus a hash of the detector description file.  Part of every cache key, so a later CeeLo
 change or a different detector file can never silently reuse stale truth.

 CeeLo carries no version string of its own; bump the literal here whenever external_libs/CeeLo is
 re-vendored (see the memory notes on the vendoring workflow). */
string mc_physics_tag()
{
  static string tag;
  if( !tag.empty() )
    return tag;

  const string angle_path = SpecUtils::append_path(
        SpecUtils::append_path( g_test_data_dir, "det_eff" ), "Angle-example-efficiency.outx" );
  ifstream in( angle_path.c_str(), ios::in | ios::binary );
  uint64_t h = 1469598103934665603ULL;   //FNV-1a
  char c;
  while( in.get( c ) )
  {
    h ^= static_cast<unsigned char>( c );
    h *= 1099511628211ULL;
  }

  ostringstream o;
  o << "ceelo=0914c4d;det=" << hex << h;
  tag = o.str();
  return tag;
}//mc_physics_tag()


/** One Monte-Carlo answer, as cached. */
struct McResult
{
  double eff = 0.0;
  double frac_sigma = 0.0;
  uint64_t n_events = 0;
  string stop_reason;
  bool from_cache = false;

  /** A row that stopped on a CAP (events or time) rather than on its precision target carries a
   larger sigma than intended; readers must not quote it as model error without saying so. */
  bool hit_cap() const { return (stop_reason != "FepPrecision") && (stop_reason != "TotalPrecision"); }
};//struct McResult


/** Append-only, tab-separated cache of Monte-Carlo results keyed on a full textual description of
 the run.  Lives at `<cachedir>/mc_cache.tsv`; with no `--cachedir=` the cache is in-memory only.

 Columns: key, energy_keV, eff, frac_sigma, n_events, stop_reason, seed, date.  The key already
 contains the energy, precision target and #mc_physics_tag(); the extra columns are for humans.
 Never read the file while a run is appending to it. */
class McCache
{
public:
  McCache()
  {
    if( g_cache_dir.empty() )
      return;

    m_path = SpecUtils::append_path( g_cache_dir, "mc_cache.tsv" );
    ifstream in( m_path.c_str() );
    string line;
    while( getline( in, line ) )
    {
      if( line.empty() || (line[0] == '#') )
        continue;
      vector<string> f;
      SpecUtils::split( f, line, "\t" );
      if( f.size() < 6 )
        continue;
      McResult r;
      r.eff = stod( f[2] );
      r.frac_sigma = stod( f[3] );
      r.n_events = stoull( f[4] );
      r.stop_reason = f[5];
      r.from_cache = true;
      m_rows[f[0]] = r;
    }
    BOOST_TEST_MESSAGE( "MC cache '" << m_path << "': " << m_rows.size() << " rows" );
  }

  bool find( const string &key, McResult &out ) const
  {
    const map<string,McResult>::const_iterator it = m_rows.find( key );
    if( it == end(m_rows) )
      return false;
    out = it->second;
    out.from_cache = true;
    return true;
  }

  void add( const string &key, const double energy_keV, const McResult &r, const unsigned seed )
  {
    m_rows[key] = r;
    if( m_path.empty() )
      return;

    const time_t now = time( nullptr );
    char date[32];
    strftime( date, sizeof(date), "%Y-%m-%dT%H:%M:%S", localtime( &now ) );
    ofstream out( m_path.c_str(), ios::app );
    out << key << '\t' << setprecision(12) << energy_keV << '\t'
        << setprecision(12) << r.eff << '\t' << setprecision(6) << r.frac_sigma << '\t'
        << r.n_events << '\t' << r.stop_reason << '\t' << seed << '\t' << date << '\n';
  }

  size_t size() const { return m_rows.size(); }

private:
  string m_path;
  map<string,McResult> m_rows;
};//class McCache


/** Runs `calc` (already configured with its source) at `energy_keV` to `target_rel` FEP precision,
 through the cache.  `config_key` must describe EVERYTHING about the source/detector configuration
 (shape, dims, centre, materials, shields); energy, precision and the physics tag are appended here.
 */
McResult run_mc( ceelo::EfficiencyCalculator &calc, McCache &cache, const string &config_key,
                 const double energy_keV, const double target_rel = 0.0025,
                 const uint64_t max_events = 16000000000ULL, const unsigned seed = 20260902,
                 const double max_wall_s = 1800.0 )
{
  ostringstream k;
  k << config_key << "|E=" << setprecision(10) << energy_keV << "|prec=" << target_rel
    << "|" << mc_physics_tag();
  const string key = k.str();

  McResult r;
  if( cache.find( key, r ) )
    return r;

  ceelo::SimulationConfig cfg;
  cfg.energy_keV = energy_keV;
  cfg.termination.target_fep_rel_precision = target_rel;
  cfg.termination.max_events = max_events;
  cfg.termination.min_events = 50000;
  // A wall cap so one hopeless row (deep + far) cannot stall the whole rung; a capped row is
  //  flagged by McResult::hit_cap() and must never be quoted as model error.
  cfg.termination.max_wall_seconds = max_wall_s;
  cfg.seed = seed;

  const ceelo::EfficiencyResult res = calc.compute( cfg );
  r.eff = res.full_energy_peak_efficiency;
  r.frac_sigma = (r.eff > 0.0) ? (res.fep_uncertainty / r.eff) : 0.0;
  r.n_events = res.num_events_simulated;
  switch( res.stop_reason )
  {
    case ceelo::StopReason::MaxEvents:      r.stop_reason = "MaxEvents"; break;
    case ceelo::StopReason::FepPrecision:   r.stop_reason = "FepPrecision"; break;
    case ceelo::StopReason::TotalPrecision: r.stop_reason = "TotalPrecision"; break;
    case ceelo::StopReason::WallTime:       r.stop_reason = "WallTime"; break;
    case ceelo::StopReason::CpuTime:        r.stop_reason = "CpuTime"; break;
  }
  r.from_cache = false;

  cache.add( key, energy_keV, r, seed );
  return r;
}//run_mc(...)


/** The per-element kernel at ONE point, with no outer quadrature: `calc` must describe a tiny
 rectangular source (the "point") plus whatever shells surround it; returns eps_fep for a photon
 emitted at its centre, i.e. prefactor * sum_i w_i * exp(-tau_i) evaluated exactly once.

 eval_rect returns the integrand times the full source volume, so dividing by that volume recovers
 the per-photon efficiency.  Use this for the point-at-depth rung instead of integrating the tiny box.
 */
double point_kernel_eff( GammaInteractionCalc::DistributedSrcCalcT<double> &calc )
{
  using namespace GammaInteractionCalc;
  BOOST_REQUIRE( calc.m_geometry == GeometryType::Rectangular );
  BOOST_REQUIRE( !calc.m_shells.empty() );
  calc.finalize_shell_coefficients();

  const std::array<double,3> &d = calc.m_shells[0].dims;
  const double vol = 8.0 * d[0] * d[1] * d[2];
  const double xx[3] = { 0.5, 0.5, 0.5 };
  return calc.eval_rect( xx, 3 ) / vol;
}//point_kernel_eff(...)

}//namespace
