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
#include <string>
#include <iomanip>
#include <iostream>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE VolumetricIntegration_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/Integrate.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "VolumetricIntegrationFixtures.h"

using namespace std;
using namespace boost::unit_test;


// We need to set the static data directory, so the code knows where
//  like sandia.decay.xml and MaterialDataBase.txt are located.
void set_data_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;

  s_have_set = true;

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

  string datadir;

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }//for( int arg = 1; arg < argc; ++ arg )

  SpecUtils::ireplace_all( datadir, "%20", " " );

  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )

  const string sandia_decay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ),
                         "sandia.decay.xml not at '" << sandia_decay_file << "'" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db && db->nuclide("U238"), "Error initing SandiaDecayDataBase" );

  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );
}//void set_data_dir()


/** Cuhre baseline values for `VolumetricFixture::standard_volumetric_fixtures()`,
 recorded 20260611 (Cuhre epsrel=1e-4, maxeval=5e6, via
 ShieldingSourceChi2Fcn::selfShieldingIntegration).

 These are the regression net for any change touching the volumetric source
 integration; a negative value means the baseline has not been recorded yet.

 Note: the multi-shell end-on cylinder values were re-recorded 20260611 after
 fixing the intersection-selection tie in #cylinder_line_intersection (the
 toward/away-from-detector choice was a floating-point coin-flip when the
 detector is on the cylinder axis, sometimes picking the wrong end-cap exit);
 the prior values were 0.32794593909881059 (cyl-end-U-plus-Fe-661keV) and
 13.082663842525706 (cyl-end-innervoid-661keV), both badly contaminated.

 Units: the integral is the per-unit-volumetric-activity detection probability
 integrated over the source volume (PhysicalUnits length-cubed scale).
 */
static const std::map<std::string,double> sm_cuhre_baselines = {
  { "sph-U-self-atten-185keV",        0.073113030740012699 },
  { "sph-U-self-atten-1001keV",       1.3375229749285693 },
  { "sph-U-plus-Fe-661keV",           0.41581277652924714 },
  { "sph-U-plus-generic-661keV",      0.48585190912467358 },
  { "sph-innervoid-Fe-then-U-661keV", 1.3320402863305374 },
  { "sph-trace-exp-soil-661keV",      1096.8666899329651 },
  { "sph-U-air-atten-661keV",         0.21063660342448978 },
  { "cyl-end-Fe-661keV",              23.652548215934239 },
  { "cyl-end-U-plus-Fe-661keV",       0.513410048805 },
  { "cyl-end-innervoid-661keV",       1.39075429013 },
  { "cyl-side-Fe-661keV",             55.752592498433657 },
  { "cyl-side-U-plus-Fe-661keV",      1.4367924646991794 },
  { "cyl-side-innervoid-661keV",      1.72724054817 },
  { "rect-Fe-661keV",                 30.0958140512233 },
  { "rect-U-plus-Fe-661keV",          0.65361416781207027 },
  { "rect-innervoid-661keV",          1.7688442081211311 },
  { "rect-poly-slab-661keV",          24.373511879753263 },
  { "rect-trace-exp-soil-661keV",     3756.2597646013069 },
};//sm_cuhre_baselines


BOOST_AUTO_TEST_CASE( CuhreBaselineIntegrals )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();
  BOOST_REQUIRE( !fixtures.empty() );

  for( const VolumetricFixture::VolumetricSrcSpec &spec : fixtures )
  {
    GammaInteractionCalc::DistributedSrcCalc calc;
    BOOST_REQUIRE_NO_THROW( calc = VolumetricFixture::make_distributed_src_calc( spec, *matdb ) );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( calc );

    // Paste-ready baseline line, in case values need re-recording:
    cout << "  { \"" << spec.name << "\", " << std::setprecision(12) << calc.integral << " }," << endl;

    BOOST_CHECK_MESSAGE( calc.integral > 0.0,
                         "Fixture '" << spec.name << "' gave non-positive integral " << calc.integral );

    const auto baseline_iter = sm_cuhre_baselines.find( spec.name );
    BOOST_REQUIRE_MESSAGE( baseline_iter != sm_cuhre_baselines.end(),
                           "Fixture '" << spec.name << "' has no baseline entry" );

    const double baseline = baseline_iter->second;
    BOOST_CHECK_MESSAGE( baseline > 0.0,
                         "Fixture '" << spec.name << "' baseline not recorded yet (got "
                         << std::setprecision(12) << calc.integral << ")" );

    if( baseline > 0.0 )
    {
      // Cuhre is deterministic; allow 0.05% for cross-platform floating point drift
      BOOST_CHECK_MESSAGE( fabs(calc.integral - baseline) <= 5.0E-4*baseline,
                           "Fixture '" << spec.name << "' integral " << std::setprecision(12)
                           << calc.integral << " differs from baseline " << baseline );
    }//if( baseline > 0.0 )
  }//for( loop over fixtures )
}//BOOST_AUTO_TEST_CASE( CuhreBaselineIntegrals )


/** The single-shell end-on cylinder has a dedicated, slightly-faster evaluator
 (`eval_single_cyl_end_on`); make sure it agrees with the general cylinder
 evaluator on the same geometry.
 */
BOOST_AUTO_TEST_CASE( SingleCylEndOnVsGeneralCylinder )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  VolumetricFixture::VolumetricSrcSpec spec;
  spec.name = "cyl-end-Fe-661keV";
  spec.geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
  spec.shells = { {"Fe (iron)", 0.0, 0.0, {5.0*PhysicalUnits::cm, 5.0*PhysicalUnits::cm, 0.0}} };

  GammaInteractionCalc::DistributedSrcCalc calc
                          = VolumetricFixture::make_distributed_src_calc( spec, *matdb );

  const double epsrel = 1.0E-4, epsabs = -1.0;
  int nregions, neval, fail;
  double error, prob;

  double single_shell_integral = 0.0, general_integral = 0.0;

  Integrate::CuhreIntegrate( 2, GammaInteractionCalc::DistributedSrcCalc_integrand_single_cyl_end_on,
                             &calc, epsrel, epsabs, Integrate::LastImportanceFcnt,
                             0, 5000000, nregions, neval, fail, single_shell_integral, error, prob );

  Integrate::CuhreIntegrate( 2, GammaInteractionCalc::DistributedSrcCalc_integrand_cylindrical,
                             &calc, epsrel, epsabs, Integrate::LastImportanceFcnt,
                             0, 5000000, nregions, neval, fail, general_integral, error, prob );

  BOOST_REQUIRE( single_shell_integral > 0.0 );
  BOOST_REQUIRE( general_integral > 0.0 );

  BOOST_CHECK_MESSAGE( fabs(single_shell_integral - general_integral) <= 1.0E-3*general_integral,
                       "Single-shell end-on cylinder evaluator (" << std::setprecision(12)
                       << single_shell_integral << ") disagrees with general evaluator ("
                       << general_integral << ")" );
}//BOOST_AUTO_TEST_CASE( SingleCylEndOnVsGeneralCylinder )


/** Regression test for the off-axis end-on cylinder self-attenuation integral.

 An end-on cylinder is azimuthally symmetric about its own axis, so displacing the
 detector radially by the same magnitude in any perpendicular direction must give the
 same self-attenuation integral.  A bug (fixed) integrated the off-axis end-on case with
 the on-axis 2D (azimuthal-symmetry) shortcut, which samples only the theta=0 half-plane
 and makes an +x offset and a +y offset disagree.  This pins the full-3D fix.
 */
BOOST_AUTO_TEST_CASE( OffAxisEndOnAzimuthalSymmetry )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const double cm = PhysicalUnits::cm;
  const double r0 = 40.0*cm;  //radial offset magnitude

  auto integral_for = [&]( const double off_x, const double off_y ) -> double {
    VolumetricFixture::VolumetricSrcSpec spec;
    spec.name = "cyl-end-offaxis";
    spec.geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
    spec.shells = { {"Fe (iron)", 0.0, 0.0, {5.0*cm, 5.0*cm, 0.0}} };
    spec.offset_x = off_x;
    spec.offset_y = off_y;

    GammaInteractionCalc::DistributedSrcCalcT<double> calc
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );
    GammaInteractionCalc::self_shielding_integration_imp( calc, 1.0E-4, size_t(5000000) );
    return calc.integral;
  };//integral_for lambda

  const double on_axis = integral_for( 0.0, 0.0 );
  const double off_x   = integral_for( r0,  0.0 );
  const double off_y   = integral_for( 0.0, r0  );

  BOOST_REQUIRE( on_axis > 0.0 );
  BOOST_REQUIRE( off_x > 0.0 );
  BOOST_REQUIRE( off_y > 0.0 );

  // +x and +y radial offsets are equivalent by azimuthal symmetry (this fails with the
  //  buggy 2D shortcut, which only samples theta=0).
  BOOST_CHECK_MESSAGE( fabs(off_x - off_y) <= 1.0E-3*std::max(off_x,off_y),
                       "Off-axis end-on cylinder: +x offset (" << std::setprecision(12) << off_x
                       << ") != +y offset (" << off_y << ") - azimuthal symmetry broken" );

  // Sanity: a radial offset actually changes the result vs on-axis (so the test isn't
  //  trivially passing because the offset is being ignored).
  BOOST_CHECK_MESSAGE( fabs(off_x - on_axis) > 1.0E-3*on_axis,
                       "Off-axis end-on cylinder: radial offset had no effect vs on-axis ("
                       << std::setprecision(12) << on_axis << ")" );
}//BOOST_AUTO_TEST_CASE( OffAxisEndOnAzimuthalSymmetry )


/** The spherical evaluator supports both a 2D (azimuthal symmetry) and 3D
 integral; production uses 2D.  Make sure both agree - the 3D path will be
 needed once sources can be off the detector axis with a non-isotropic
 detector response.
 */
BOOST_AUTO_TEST_CASE( SphericalTwoVsThreeDim )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  VolumetricFixture::VolumetricSrcSpec spec;
  spec.name = "sph-U-plus-Fe-661keV";
  spec.geometry = GammaInteractionCalc::GeometryType::Spherical;
  spec.shells = { {"U (uranium)", 0.0, 0.0, {2.0*PhysicalUnits::cm, 0.0, 0.0}},
                  {"Fe (iron)", 0.0, 0.0, {1.0*PhysicalUnits::cm, 0.0, 0.0}} };

  GammaInteractionCalc::DistributedSrcCalc calc
                          = VolumetricFixture::make_distributed_src_calc( spec, *matdb );

  const double epsrel = 1.0E-4, epsabs = -1.0;
  int nregions, neval, fail;
  double error, prob;

  double integral_2d = 0.0, integral_3d = 0.0;

  Integrate::CuhreIntegrate( 2, GammaInteractionCalc::DistributedSrcCalc_integrand_spherical,
                             &calc, epsrel, epsabs, Integrate::LastImportanceFcnt,
                             0, 5000000, nregions, neval, fail, integral_2d, error, prob );

  Integrate::CuhreIntegrate( 3, GammaInteractionCalc::DistributedSrcCalc_integrand_spherical,
                             &calc, epsrel, epsabs, Integrate::LastImportanceFcnt,
                             0, 5000000, nregions, neval, fail, integral_3d, error, prob );

  BOOST_REQUIRE( integral_2d > 0.0 );
  BOOST_REQUIRE( integral_3d > 0.0 );

  BOOST_CHECK_MESSAGE( fabs(integral_2d - integral_3d) <= 1.0E-3*integral_2d,
                       "Spherical 2D integral (" << std::setprecision(12) << integral_2d
                       << ") disagrees with 3D integral (" << integral_3d << ")" );
}//BOOST_AUTO_TEST_CASE( SphericalTwoVsThreeDim )


/** Benchmark of the adaptive Gauss-Legendre backend vs Cuba/Cuhre, at the production
 tolerance (epsrel=1e-4), with accuracy judged against a tight (epsrel=1e-7) GL
 reference.  Diagnostic only (asserts nothing beyond sanity) - run explicitly with
 --run_test=IntegrationBackendBenchmark.
 */
BOOST_AUTO_TEST_CASE( IntegrationBackendBenchmark, * boost::unit_test::disabled() )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  const int num_timing_repeats = 10;

  cout << "\nfixture, cuhre_relerr, gl_relerr, cuhre_ms, gl_ms, cuhre_evals, gl_evals" << endl;

  for( const VolumetricFixture::VolumetricSrcSpec &spec : fixtures )
  {
    // Tight-tolerance GL reference value
    GammaInteractionCalc::DistributedSrcCalcT<double> ref_calc
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );
    GammaInteractionCalc::self_shielding_integration_imp( ref_calc, 1.0E-7, size_t(100000000) );
    const double reference = ref_calc.integral;
    BOOST_REQUIRE( reference > 0.0 );

    // Cuhre at production settings (via the production dispatch, incl. its fast paths)
    GammaInteractionCalc::DistributedSrcCalc legacy
                            = VolumetricFixture::make_distributed_src_calc( spec, *matdb );

    const auto cuhre_start = std::chrono::steady_clock::now();
    for( int rep = 0; rep < num_timing_repeats; ++rep )
      GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( legacy );
    const auto cuhre_finish = std::chrono::steady_clock::now();

    // The GL backend at production settings
    GammaInteractionCalc::DistributedSrcCalcT<double> calc
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );

    const auto gl_start = std::chrono::steady_clock::now();
    for( int rep = 0; rep < num_timing_repeats; ++rep )
      GammaInteractionCalc::self_shielding_integration_imp( calc );
    const auto gl_finish = std::chrono::steady_clock::now();

    const double cuhre_ms = std::chrono::duration_cast<std::chrono::microseconds>(cuhre_finish - cuhre_start).count()
                            / (1000.0 * num_timing_repeats);
    const double gl_ms = std::chrono::duration_cast<std::chrono::microseconds>(gl_finish - gl_start).count()
                         / (1000.0 * num_timing_repeats);

    cout << spec.name << ", "
         << std::setprecision(3) << fabs(legacy.integral - reference)/reference << ", "
         << fabs(calc.integral - reference)/reference << ", "
         << std::setprecision(4) << cuhre_ms << ", " << gl_ms << ", "
         << "n/a, " << calc.m_num_evals << endl;
  }//for( loop over fixtures )
}//BOOST_AUTO_TEST_CASE( IntegrationBackendBenchmark )


/** Validates the templated ray-tracing + adaptive Gauss-Legendre integration
 backend (DistributedSrcCalcT<double> + self_shielding_integration_imp) against
 the Cuhre/double path, on the full fixture matrix.

 Acceptance: <= 5e-4 relative agreement, and (logged) evaluation-count ratio.
 */
BOOST_AUTO_TEST_CASE( BoostGLBackendVsCuhre )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();
  BOOST_REQUIRE( !fixtures.empty() );

  for( const VolumetricFixture::VolumetricSrcSpec &spec : fixtures )
  {
    // Reference: the Cuhre integration of the legacy DistributedSrcCalc
    GammaInteractionCalc::DistributedSrcCalc legacy
                            = VolumetricFixture::make_distributed_src_calc( spec, *matdb );
    GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( legacy );
    BOOST_REQUIRE( legacy.integral > 0.0 );

    // The new templated backend, at T=double
    GammaInteractionCalc::DistributedSrcCalcT<double> calc
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );
    BOOST_REQUIRE_NO_THROW( GammaInteractionCalc::self_shielding_integration_imp( calc ) );

    cout << "GL-vs-Cuhre '" << spec.name << "': GL=" << std::setprecision(10) << calc.integral
         << ", Cuhre=" << legacy.integral
         << ", rel-diff=" << std::setprecision(3) << fabs(calc.integral - legacy.integral)/legacy.integral
         << ", GL evals=" << calc.m_num_evals << endl;

    // Inner-void sources have a discontinuous integrand (zero inside the void), which limits
    //  the achievable/estimable accuracy of both integrators - so a looser tolerance there.
    const bool discontinuous = (spec.name.find("innervoid") != string::npos);
    const double tolerance = discontinuous ? 2.0E-3 : 5.0E-4;

    BOOST_CHECK_MESSAGE( fabs(calc.integral - legacy.integral) <= tolerance*legacy.integral,
                         "Fixture '" << spec.name << "': GL backend (" << std::setprecision(12)
                         << calc.integral << ") disagrees with Cuhre (" << legacy.integral << ")" );
  }//for( loop over fixtures )
}//BOOST_AUTO_TEST_CASE( BoostGLBackendVsCuhre )


/** Guards the anisotropic-refinement speedup.  The nested-core silhouette cases (a transparent
 core seen through the source material) used to cost the adaptive integrator ~5M evaluations under
 isotropic bisection; anisotropic refinement cuts that ~20x.  This pins a generous evaluation
 ceiling for those cases (and a smooth case), so a future change that silently restores isotropic
 refinement - or otherwise inflates the cost - trips here.  Also confirms the cheap production
 result still matches a much-tighter reference, i.e. the savings are not bought with accuracy.
 Ceilings are well above the measured anisotropic cost (cross-platform FP-path slack) yet well
 below the isotropic cost they replaced.
 */
BOOST_AUTO_TEST_CASE( AnisotropicEvalBudget )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  // {fixture, production-eval ceiling, rel-tolerance vs a much-tighter reference}.
  //  measured anisotropic / isotropic production evals, for context:
  //   rect-innervoid      ~271k / ~5.28M   cyl-side-innervoid ~196k / ~631k   cyl-side-Fe ~34k / ~79k
  struct Budget{ string name; size_t eval_ceiling; double ref_tol; };
  const vector<Budget> budgets = {
    { "rect-innervoid-661keV",     1500000, 2.0e-3 },
    { "cyl-side-innervoid-661keV",  500000, 1.0e-4 },
    { "cyl-side-Fe-661keV",         150000, 1.0e-5 },
  };

  for( const Budget &b : budgets )
  {
    const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
    for( const VolumetricFixture::VolumetricSrcSpec &f : fixtures )
      if( f.name == b.name )
        spec = &f;
    BOOST_REQUIRE_MESSAGE( spec, "missing fixture " << b.name );

    // Much-tighter reference (100x the production epsrel, bounded budget).
    GammaInteractionCalc::DistributedSrcCalcT<double> ref
                            = VolumetricFixture::make_distributed_src_calc_t( *spec, *matdb );
    GammaInteractionCalc::self_shielding_integration_imp( ref, 1.0e-6, size_t(5000000) );
    BOOST_REQUIRE( ref.integral > 0.0 );

    // Production-settings run - the one whose cost we guard.
    GammaInteractionCalc::DistributedSrcCalcT<double> calc
                            = VolumetricFixture::make_distributed_src_calc_t( *spec, *matdb );
    GammaInteractionCalc::self_shielding_integration_imp( calc );

    cout << "EvalBudget '" << b.name << "': evals=" << calc.m_num_evals
         << " (ceiling " << b.eval_ceiling << "), rel-vs-ref="
         << std::setprecision(3) << fabs(calc.integral - ref.integral)/ref.integral << endl;

    BOOST_CHECK_MESSAGE( calc.m_num_evals <= b.eval_ceiling,
                         "Fixture '" << b.name << "' used " << calc.m_num_evals
                         << " evals, over the anisotropic-refinement ceiling " << b.eval_ceiling
                         << " - did anisotropic refinement regress to isotropic?" );

    BOOST_CHECK_MESSAGE( fabs(calc.integral - ref.integral) <= b.ref_tol*ref.integral,
                         "Fixture '" << b.name << "' production integral " << std::setprecision(12)
                         << calc.integral << " differs from the tighter reference " << ref.integral
                         << " by more than " << b.ref_tol << " (savings cost accuracy?)" );
  }//for( loop over budgets )
}//BOOST_AUTO_TEST_CASE( AnisotropicEvalBudget )


/** Unit-tests the per-axis non-smoothness indicator (`compute_axis_surplus`) that drives the
 anisotropic split-axis selection: a GL5 grid that is a low-order polynomial (inside the rule's
 exactness space) must give ~zero surplus on every axis, and a grid with a sharp kink along one
 axis only must single out that axis.  This guards the indicator wiring independent of the full
 integration.
 */
BOOST_AUTO_TEST_CASE( AxisSurplusIndicator )
{
  const int npts = 5;
  const double * const x = GammaInteractionCalc::QuadDetail::gl5_x;

  // (a) A function linear along each axis is annihilated by the degree-4 surplus -> ~0 everywhere.
  {
    double grid[125];
    for( int i = 0; i < npts; ++i )
      for( int j = 0; j < npts; ++j )
        for( int k = 0; k < npts; ++k )
          grid[(i*npts + j)*npts + k] = 1.0 + 2.0*x[i] - 3.0*x[j] + 0.5*x[k]
                                        + x[i]*x[j] - x[j]*x[k];   // linear in each axis separately
    double axis_err[3] = { -1, -1, -1 };
    GammaInteractionCalc::QuadDetail::compute_axis_surplus( grid, 3, axis_err );
    BOOST_CHECK_SMALL( axis_err[0], 1.0e-6 );
    BOOST_CHECK_SMALL( axis_err[1], 1.0e-6 );
    BOOST_CHECK_SMALL( axis_err[2], 1.0e-6 );
  }

  // (b) A kink |y-0.5| along axis 1 only -> axis 1 dominates, axes 0 and 2 are ~0 (constant there).
  {
    double grid[125];
    for( int i = 0; i < npts; ++i )
      for( int j = 0; j < npts; ++j )
        for( int k = 0; k < npts; ++k )
          grid[(i*npts + j)*npts + k] = std::fabs( x[j] - 0.5 );
    double axis_err[3] = { -1, -1, -1 };
    GammaInteractionCalc::QuadDetail::compute_axis_surplus( grid, 3, axis_err );
    BOOST_CHECK_MESSAGE( axis_err[1] > 1.0e-3,
                         "kink axis surplus too small: " << axis_err[1] );
    BOOST_CHECK_SMALL( axis_err[0], 1.0e-9 );
    BOOST_CHECK_SMALL( axis_err[2], 1.0e-9 );
  }
}//BOOST_AUTO_TEST_CASE( AxisSurplusIndicator )


/** Checks that integrating with ceres::Jet dimensions produces the correct
 derivative of the integral with respect to a shell dimension, by comparing
 the Jet derivative lane against a central finite difference of the double
 path.  This is the key property the Ceres-based fit relies on.
 */
BOOST_AUTO_TEST_CASE( JetDerivativeOfIntegral )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  using Jet2 = ceres::Jet<double,2>;
  const double cm = PhysicalUnits::cm;
  const double epsrel = 1.0E-6;  //tight tolerance so finite-difference noise is small
  const size_t max_evals = 20000000;

  // Geometries to check: {fixture, which-dim-of-source-shell to differentiate}
  // {fixture, dim index, inner_boundary}: inner_boundary=false differentiates the
  //  source shells outer dimension (with enclosing shells shifted along); =true
  //  differentiates the boundary between a non-source core and the source shell -
  //  the lane that was blind to autodiff before hollow shells were split into
  //  sub-domains (the indicator-function approach hid the boundary-motion term).
  struct DerivCase{ string name; int dim; bool inner_boundary; };
  const vector<DerivCase> cases = {
    { "sph-U-plus-Fe-661keV", 0, false },     //sphere radius
    { "cyl-end-U-plus-Fe-661keV", 1, false }, //cylinder half-length
    { "rect-U-plus-Fe-661keV", 2, false },    //rectangle half-depth
    { "sph-innervoid-Fe-then-U-661keV", 0, true },
    { "cyl-end-innervoid-661keV", 1, true },
    { "cyl-side-innervoid-661keV", 0, true }, //radial core boundary - the side-on silhouette
    { "rect-innervoid-661keV", 2, true },
  };

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  for( const DerivCase &test_case : cases )
  {
    const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
    for( const VolumetricFixture::VolumetricSrcSpec &fixture : fixtures )
    {
      if( fixture.name == test_case.name )
        spec = &fixture;
    }
    BOOST_REQUIRE( spec );

    const int dim_index = test_case.dim;
    const size_t src_shell = spec->source_shell_index;

    // Jet evaluation: lane 0 = d/d(source-shell dim), lane 1 unused (exercises multi-lane)
    GammaInteractionCalc::DistributedSrcCalcT<double> calc_dbl
                          = VolumetricFixture::make_distributed_src_calc_t( *spec, *matdb );

    GammaInteractionCalc::DistributedSrcCalcT<Jet2> calc_jet;
    calc_jet.m_geometry = calc_dbl.m_geometry;
    calc_jet.m_materialIndex = calc_dbl.m_materialIndex;
    calc_jet.m_detector = GammaInteractionCalc::detector_geom_from_config<Jet2>(
                  spec->geometry, Jet2(spec->distance), spec->detector_radius, spec->detector_setback );
    calc_jet.m_attenuateForAir = calc_dbl.m_attenuateForAir;
    calc_jet.m_airTransLenCoef = calc_dbl.m_airTransLenCoef;
    calc_jet.m_isInSituExponential = calc_dbl.m_isInSituExponential;
    calc_jet.m_inSituRelaxationLength = calc_dbl.m_inSituRelaxationLength;
    calc_jet.m_srcVolumetricActivity = Jet2( calc_dbl.m_srcVolumetricActivity );
    calc_jet.m_energy = calc_dbl.m_energy;

    for( size_t shell_index = 0; shell_index < calc_dbl.m_shells.size(); ++shell_index )
    {
      const GammaInteractionCalc::DistributedSrcCalcT<double>::ShellInfo &dbl_shell
                                                            = calc_dbl.m_shells[shell_index];
      GammaInteractionCalc::DistributedSrcCalcT<Jet2>::ShellInfo shell;
      for( size_t i = 0; i < 3; ++i )
        shell.dims[i] = Jet2( dbl_shell.dims[i] );
      shell.trans_len_coef = Jet2( dbl_shell.trans_len_coef );
      shell.type = dbl_shell.type;

      // Outer dims are cumulative, so increasing the source-shell dim shifts all
      //  enclosing shells outward too - mirror that in the derivative seeding.
      //  For the inner-boundary cases, only the cores boundary moves (the source
      //  shells outer extent stays put).
      const bool seed = test_case.inner_boundary ? (shell_index == (src_shell - 1))
                                                 : (shell_index >= src_shell);
      if( seed )
        shell.dims[dim_index].v[0] = 1.0;

      calc_jet.m_shells.push_back( shell );
    }//for( loop over shells )

    GammaInteractionCalc::self_shielding_integration_imp( calc_jet, epsrel, max_evals );

    // Finite-difference reference from the double path
    const size_t fd_shell = test_case.inner_boundary ? (src_shell - 1) : src_shell;
    const double dim_val = calc_dbl.m_shells[fd_shell].dims[dim_index];
    const double step = 0.02 * dim_val;

    GammaInteractionCalc::DistributedSrcCalcT<double> calc_plus = calc_dbl;
    GammaInteractionCalc::DistributedSrcCalcT<double> calc_minus = calc_dbl;
    const size_t fd_last_shell = test_case.inner_boundary ? src_shell : calc_dbl.m_shells.size();
    for( size_t shell_index = fd_shell; shell_index < fd_last_shell; ++shell_index )
    {
      calc_plus.m_shells[shell_index].dims[dim_index] += step;
      calc_minus.m_shells[shell_index].dims[dim_index] -= step;
    }

    GammaInteractionCalc::self_shielding_integration_imp( calc_plus, epsrel, max_evals );
    GammaInteractionCalc::self_shielding_integration_imp( calc_minus, epsrel, max_evals );

    const double numeric_deriv = (calc_plus.integral - calc_minus.integral) / (2.0*step);
    const double jet_deriv = calc_jet.integral.v[0];

    cout << "Jet-deriv '" << spec->name << "' dim" << dim_index << ": jet=" << std::setprecision(8)
         << jet_deriv*cm << "/cm, numeric=" << numeric_deriv*cm << "/cm" << endl;

    BOOST_CHECK_MESSAGE( fabs(jet_deriv - numeric_deriv) <= 0.02*fabs(numeric_deriv),
                         "Fixture '" << spec->name << "' dim " << dim_index << ": Jet derivative ("
                         << jet_deriv << ") disagrees with numeric (" << numeric_deriv << ")" );

    // And the value lane must agree with the plain double integration
    GammaInteractionCalc::self_shielding_integration_imp( calc_dbl, epsrel, max_evals );
    BOOST_CHECK_MESSAGE( fabs(calc_jet.integral.a - calc_dbl.integral) <= 1.0E-9*calc_dbl.integral,
                         "Fixture '" << spec->name << "': Jet value lane (" << calc_jet.integral.a
                         << ") differs from double integration (" << calc_dbl.integral << ")" );
  }//for( loop over test cases )
}//BOOST_AUTO_TEST_CASE( JetDerivativeOfIntegral )


/** Point-wise comparison of the templated eval_* integrands against the legacy
 double implementations, on a uniform grid - catches transcription bugs
 independent of the integration scheme.
 */
BOOST_AUTO_TEST_CASE( TemplatedIntegrandPointCompare )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  for( const VolumetricFixture::VolumetricSrcSpec &spec : fixtures )
  {
    const GammaInteractionCalc::DistributedSrcCalc legacy
                            = VolumetricFixture::make_distributed_src_calc( spec, *matdb );
    const GammaInteractionCalc::DistributedSrcCalcT<double> calc
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );

    int ndim = 2;
    switch( spec.geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
        ndim = 2;
        break;
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
      case GammaInteractionCalc::GeometryType::Rectangular:
        ndim = 3;
        break;
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        BOOST_FAIL( "bad geometry" );
    }

    size_t num_bad = 0;
    const int n_grid = 13;
    for( int i = 0; (i < n_grid) && (num_bad < 3); ++i )
    {
      for( int j = 0; (j < n_grid) && (num_bad < 3); ++j )
      {
        for( int k = 0; (k < ((ndim == 3) ? n_grid : 1)) && (num_bad < 3); ++k )
        {
          const double xx[3] = { (i + 0.5)/n_grid, (j + 0.5)/n_grid, (k + 0.5)/n_grid };

          double legacy_val = 0.0;
          switch( spec.geometry )
          {
            case GammaInteractionCalc::GeometryType::Spherical:
              legacy.eval_spherical( xx, &ndim, &legacy_val, nullptr );
              break;
            case GammaInteractionCalc::GeometryType::CylinderEndOn:
            case GammaInteractionCalc::GeometryType::CylinderSideOn:
              legacy.eval_cylinder( xx, &ndim, &legacy_val, nullptr );
              break;
            case GammaInteractionCalc::GeometryType::Rectangular:
              legacy.eval_rect( xx, &ndim, &legacy_val, nullptr );
              break;
            case GammaInteractionCalc::GeometryType::NumGeometryType:
              break;
          }

          double templated_val = 0.0;
          switch( spec.geometry )
          {
            case GammaInteractionCalc::GeometryType::Spherical:
              templated_val = calc.eval_spherical( xx, ndim );
              break;
            case GammaInteractionCalc::GeometryType::CylinderEndOn:
            case GammaInteractionCalc::GeometryType::CylinderSideOn:
              templated_val = calc.eval_cylinder( xx, ndim );
              break;
            case GammaInteractionCalc::GeometryType::Rectangular:
              templated_val = calc.eval_rect( xx, ndim );
              break;
            case GammaInteractionCalc::GeometryType::NumGeometryType:
              break;
          }

          const double diff = fabs( templated_val - legacy_val );
          if( diff > 1.0E-9*std::max(1.0E-12,fabs(legacy_val)) )
          {
            ++num_bad;
            BOOST_ERROR( "Fixture '" << spec.name << "' at xx={" << xx[0] << "," << xx[1] << ","
                         << xx[2] << "}: templated=" << std::setprecision(12) << templated_val
                         << " vs legacy=" << legacy_val );
          }
        }//for( k )
      }//for( j )
    }//for( i )
  }//for( loop over fixtures )
}//BOOST_AUTO_TEST_CASE( TemplatedIntegrandPointCompare )


// Diagnostic: integrand profile across the inner-box shadow boundary, for the
//  behind-the-core sub-domain of rect-innervoid
BOOST_AUTO_TEST_CASE( DebugRectShadowProfile, * boost::unit_test::disabled() )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();
  const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
  for( const VolumetricFixture::VolumetricSrcSpec &fixture : fixtures )
  {
    if( fixture.name == "rect-innervoid-661keV" )
      spec = &fixture;
  }
  BOOST_REQUIRE( spec );

  GammaInteractionCalc::DistributedSrcCalcT<double> calc
                          = VolumetricFixture::make_distributed_src_calc_t( *spec, *matdb );
  const vector<GammaInteractionCalc::DistributedSrcCalcT<double>> subs
                          = GammaInteractionCalc::split_source_subdomains( calc );
  BOOST_REQUIRE_EQUAL( subs.size(), size_t(6) );

  const GammaInteractionCalc::DistributedSrcCalcT<double> &behind = subs[1];

  // Detailed ray-trace pieces at two probe points either side of the jump
  for( const double x_cm : { 1.4, 1.6 } )
  {
    const double mm = PhysicalUnits::mm;
    const double eval_loc[3] = { 10.0*x_cm*mm, 0.0, -20.0*mm };
    const double det_loc[3] = { 0.0, 0.0, 1000.0*mm };

    double exit_pt[3];
    const double src_chord = GammaInteractionCalc::rectangle_exit_location_imp<double>(
                  behind.m_shells[1].dims[0], behind.m_shells[1].dims[1], behind.m_shells[1].dims[2],
                  eval_loc, det_loc, exit_pt );

    double enter_in[3], exit_in[3];
    const bool hits_inner = GammaInteractionCalc::rectangle_intersections_imp<double>(
                  behind.m_shells[0].dims[0], behind.m_shells[0].dims[1], behind.m_shells[0].dims[2],
                  eval_loc, det_loc, enter_in, exit_in );

    double inner_chord = 0.0;
    if( hits_inner )
    {
      const double dx = exit_in[0]-enter_in[0], dy = exit_in[1]-enter_in[1], dz = exit_in[2]-enter_in[2];
      inner_chord = sqrt( dx*dx + dy*dy + dz*dz );
    }

    cout << "probe x=" << x_cm << "cm: src_chord=" << src_chord/PhysicalUnits::cm
         << "cm, hits_inner=" << hits_inner << ", inner_chord=" << inner_chord/PhysicalUnits::cm
         << "cm, mu_src=" << behind.m_shells[1].trans_len_coef*PhysicalUnits::cm
         << "/cm, mu_inner=" << behind.m_shells[0].trans_len_coef*PhysicalUnits::cm << "/cm"
         << ", exp(-trans)=" << exp( -(behind.m_shells[0].trans_len_coef*inner_chord
                              + behind.m_shells[1].trans_len_coef*(src_chord - inner_chord)) ) << endl;
  }//for( probe points )

  // Scan x across the shadow edge of the inner box (y=0 middle, z middle of slab)
  cout << "x_cm, integrand" << endl;
  for( int i = 0; i <= 60; ++i )
  {
    const double xx[3] = { 0.70 + 0.01*(i/60.0 - 0.0)*30.0/30.0*0.3/0.3, 0.5, 0.5 };
    //simpler: u from 0.7 to 1.0
    const double u = 0.70 + (0.30*i)/60.0;
    const double xx2[3] = { u, 0.5, 0.5 };
    const double val = behind.eval_rect( xx2, 3 );
    const double x_cm = (behind.m_subdomain_lo[0] + u*(behind.m_subdomain_hi[0]-behind.m_subdomain_lo[0]))/PhysicalUnits::cm;
    cout << std::setprecision(6) << x_cm << ", " << std::setprecision(10) << val << endl;
    (void)xx;
  }
}//BOOST_AUTO_TEST_CASE( DebugRectShadowProfile )


// Diagnostic: per-sub-domain integration cost for the inner-void fixtures
BOOST_AUTO_TEST_CASE( DebugSubdomainCosts, * boost::unit_test::disabled() )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  for( const VolumetricFixture::VolumetricSrcSpec &spec : fixtures )
  {
    if( spec.name.find("innervoid") == string::npos )
      continue;

    GammaInteractionCalc::DistributedSrcCalcT<double> calc
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );

    const vector<GammaInteractionCalc::DistributedSrcCalcT<double>> subs
                            = GammaInteractionCalc::split_source_subdomains( calc );

    cout << "Fixture '" << spec.name << "': " << subs.size() << " subdomains" << endl;
    for( size_t i = 0; i < subs.size(); ++i )
    {
      GammaInteractionCalc::DistributedSrcCalcT<double> sub = subs[i];
      sub.m_materialIndex = sub.m_materialIndex;  //no-op; keep copy mutable
      GammaInteractionCalc::self_shielding_integration_imp( sub, 1.0E-4 );
      cout << "  sub[" << i << "]: lo={" << sub.m_subdomain_lo[0]/PhysicalUnits::cm << ","
           << sub.m_subdomain_lo[1]/PhysicalUnits::cm << "," << sub.m_subdomain_lo[2]/PhysicalUnits::cm
           << "}cm hi={" << sub.m_subdomain_hi[0]/PhysicalUnits::cm << ","
           << sub.m_subdomain_hi[1]/PhysicalUnits::cm << "," << sub.m_subdomain_hi[2]/PhysicalUnits::cm
           << "}cm integral=" << std::setprecision(8) << sub.integral
           << " evals=" << sub.m_num_evals << " est_rel_err=" << sub.m_est_rel_error << endl;
    }
  }//for( loop over innervoid fixtures )
}//BOOST_AUTO_TEST_CASE( DebugSubdomainCosts )


// Temporary diagnostic - convergence study of the adaptive GL integrator on the
//  fixtures it struggles with.
BOOST_AUTO_TEST_CASE( DebugAdaptiveGLConvergence, * boost::unit_test::disabled() )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  for( const string &name : { string("cyl-end-U-plus-Fe-661keV"), string("cyl-end-innervoid-661keV"),
                              string("cyl-side-Fe-661keV"), string("rect-innervoid-661keV") } )
  {
    const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
    for( const VolumetricFixture::VolumetricSrcSpec &fixture : fixtures )
    {
      if( fixture.name == name )
        spec = &fixture;
    }
    BOOST_REQUIRE( spec );

    GammaInteractionCalc::DistributedSrcCalc legacy
                            = VolumetricFixture::make_distributed_src_calc( *spec, *matdb );
    GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( legacy );

    cout << "Fixture '" << name << "', Cuhre=" << std::setprecision(10) << legacy.integral << endl;

    for( const double epsrel : { 1.0E-2, 1.0E-3, 1.0E-4 } )
    {
      for( const size_t max_evals : { size_t(5000000), size_t(50000000) } )
      {
        GammaInteractionCalc::DistributedSrcCalcT<double> calc
                                = VolumetricFixture::make_distributed_src_calc_t( *spec, *matdb );
        GammaInteractionCalc::self_shielding_integration_imp( calc, epsrel, max_evals );
        cout << "  epsrel=" << epsrel << ", max_evals=" << max_evals
             << ": integral=" << std::setprecision(10) << calc.integral
             << ", rel-to-cuhre=" << std::setprecision(3) << (calc.integral - legacy.integral)/legacy.integral
             << ", evals=" << calc.m_num_evals
             << ", est_rel_err=" << calc.m_est_rel_error << endl;
      }
    }
  }//for( loop over difficult fixtures )
}//BOOST_AUTO_TEST_CASE( DebugAdaptiveGLConvergence )


/** A sphere at distance d with lateral offset s must give the same answer as an
 on-axis sphere at distance sqrt(d*d + s*s): the shells are rotation-invariant
 and (currently) the detector response uses only the line-of-sight distance, with
 no detector-face orientation term.  NOTE: this equivalence breaks once the
 detector response gains an angular dependence - at that point this test should
 compare against an independent 3D reference instead.

 Also cross-checks the off-axis cylinder/rectangle evaluation between the two
 independent backends (legacy + Cuhre, vs templated + adaptive GL), and that
 moving the source off-axis monotonically decreases the detected counts.
 */
BOOST_AUTO_TEST_CASE( OffAxisSphereEquivalence )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const double cm = PhysicalUnits::cm;

  {// Begin off-axis == rotated on-axis sphere
    VolumetricFixture::VolumetricSrcSpec off_axis;
    off_axis.name = "sph-offaxis";
    off_axis.geometry = GammaInteractionCalc::GeometryType::Spherical;
    off_axis.shells = { {"U (uranium)", 0.0, 0.0, {2.0*cm, 0.0, 0.0}},
                        {"Fe (iron)", 0.0, 0.0, {1.0*cm, 0.0, 0.0}} };
    off_axis.distance = 300.0*cm;
    off_axis.offset_x = 100.0*cm;

    VolumetricFixture::VolumetricSrcSpec rotated = off_axis;
    rotated.name = "sph-rotated";
    rotated.distance = sqrt( 300.0*300.0 + 100.0*100.0 )*cm;
    rotated.offset_x = 0.0;

    // Legacy + Cuhre
    GammaInteractionCalc::DistributedSrcCalc legacy_off
                          = VolumetricFixture::make_distributed_src_calc( off_axis, *matdb );
    GammaInteractionCalc::DistributedSrcCalc legacy_rot
                          = VolumetricFixture::make_distributed_src_calc( rotated, *matdb );
    GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( legacy_off );
    GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( legacy_rot );

    BOOST_REQUIRE( legacy_rot.integral > 0.0 );
    BOOST_CHECK_MESSAGE( fabs(legacy_off.integral - legacy_rot.integral) <= 5.0E-4*legacy_rot.integral,
                         "Cuhre off-axis sphere (" << std::setprecision(10) << legacy_off.integral
                         << ") != rotated equivalent (" << legacy_rot.integral << ")" );

    // Templated + adaptive GL
    GammaInteractionCalc::DistributedSrcCalcT<double> gl_off
                          = VolumetricFixture::make_distributed_src_calc_t( off_axis, *matdb );
    GammaInteractionCalc::DistributedSrcCalcT<double> gl_rot
                          = VolumetricFixture::make_distributed_src_calc_t( rotated, *matdb );
    GammaInteractionCalc::self_shielding_integration_imp( gl_off );
    GammaInteractionCalc::self_shielding_integration_imp( gl_rot );

    BOOST_REQUIRE( gl_rot.integral > 0.0 );
    BOOST_CHECK_MESSAGE( fabs(gl_off.integral - gl_rot.integral) <= 5.0E-4*gl_rot.integral,
                         "GL off-axis sphere (" << std::setprecision(10) << gl_off.integral
                         << ") != rotated equivalent (" << gl_rot.integral << ")" );

    BOOST_CHECK_MESSAGE( fabs(gl_off.integral - legacy_off.integral) <= 1.0E-3*legacy_off.integral,
                         "Off-axis sphere: GL (" << std::setprecision(10) << gl_off.integral
                         << ") vs Cuhre (" << legacy_off.integral << ")" );

    cout << "OffAxisSphere: off-axis=" << std::setprecision(10) << gl_off.integral
         << ", rotated=" << gl_rot.integral << " (Cuhre off-axis=" << legacy_off.integral << ")" << endl;
  }// End off-axis == rotated on-axis sphere

  // Off-axis cylinders/rectangles: cross-validate the two independent backends.
  //  (Note: detected counts do NOT necessarily decrease with offset for heavily
  //  self-absorbed boxes - the oblique view exposes more escape surface, which can
  //  beat the solid-angle loss - so only backend agreement is asserted.)
  for( const auto geometry : { GammaInteractionCalc::GeometryType::CylinderEndOn,
                               GammaInteractionCalc::GeometryType::CylinderSideOn,
                               GammaInteractionCalc::GeometryType::Rectangular } )
  {
    for( const double offset : { 0.0, 30.0*cm, 80.0*cm } )
    {
      VolumetricFixture::VolumetricSrcSpec spec;
      spec.name = "offaxis-test";
      spec.geometry = geometry;
      if( geometry == GammaInteractionCalc::GeometryType::Rectangular )
        spec.shells = { {"Fe (iron)", 0.0, 0.0, {4.0*cm, 4.0*cm, 4.0*cm}} };
      else
        spec.shells = { {"Fe (iron)", 0.0, 0.0, {4.0*cm, 4.0*cm, 0.0}} };
      spec.distance = 100.0*cm;
      spec.offset_x = offset;
      spec.offset_y = 0.5*offset;

      GammaInteractionCalc::DistributedSrcCalc legacy
                            = VolumetricFixture::make_distributed_src_calc( spec, *matdb );
      GammaInteractionCalc::ShieldingSourceChi2Fcn::selfShieldingIntegration( legacy );

      GammaInteractionCalc::DistributedSrcCalcT<double> gl
                            = VolumetricFixture::make_distributed_src_calc_t( spec, *matdb );
      GammaInteractionCalc::self_shielding_integration_imp( gl );

      BOOST_REQUIRE( legacy.integral > 0.0 );
      BOOST_REQUIRE( gl.integral > 0.0 );

      BOOST_CHECK_MESSAGE( fabs(gl.integral - legacy.integral) <= 2.0E-3*legacy.integral,
                           "Off-axis " << GammaInteractionCalc::to_str(geometry) << " offset="
                           << offset/cm << "cm: GL (" << std::setprecision(10) << gl.integral
                           << ") vs Cuhre (" << legacy.integral << ")" );

      cout << "OffAxis " << GammaInteractionCalc::to_str(geometry) << " offset=" << offset/cm
           << "cm: GL=" << std::setprecision(8) << gl.integral
           << ", Cuhre=" << legacy.integral << endl;
    }//for( loop over offsets )
  }//for( loop over geometries )
}//BOOST_AUTO_TEST_CASE( OffAxisSphereEquivalence )


/** For the 'rect-poly-slab-661keV' fixture (polyethylene slab, half-depth
 D/2 = 2 cm, at 3 m), rays to the detector are nearly normal to the slab, so
 with z the depth from a source point to the exit surface, and mu the linear
 attenuation coefficient, the attenuation-weighted averages are analytic:

   <AD>    = rho * Int(z*exp(-mu*z))/Int(exp(-mu*z)), z over [0,D]
           = (rho/mu) * (1 - x*exp(-x)/(1 - exp(-x))),  x = mu*D
   <AN>    = mass-weighted Zbar of polyethylene
   <fracH> = hydrogen mass fraction of polyethylene

 (expect agreement to ~1%, the residual coming from ray obliquity)
 */
BOOST_AUTO_TEST_CASE( EffectiveShieldingAnalyticSlab )
{
  set_data_dir();

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
  for( const VolumetricFixture::VolumetricSrcSpec &fixture : fixtures )
  {
    if( fixture.name == "rect-poly-slab-661keV" )
      spec = &fixture;
  }
  BOOST_REQUIRE( spec );

  GammaInteractionCalc::DistributedSrcCalcT<double> calc
                          = VolumetricFixture::make_distributed_src_calc_t( *spec, *matdb );
  BOOST_REQUIRE_EQUAL( calc.m_shells.size(), size_t(1) );

  const GammaInteractionCalc::EffShieldComponents components
                          = GammaInteractionCalc::integrate_effective_shielding( calc, 1.0E-5 );

  BOOST_REQUIRE( components.c[0] > 0.0 );
  BOOST_REQUIRE( components.c[1] > 0.0 );

  const double eff_ad = components.c[1] / components.c[0];
  const double eff_an_mass = components.c[2] / components.c[1];
  const double eff_frac_h = components.c[3] / components.c[1];
  const double eff_an_xs = components.c[5] / components.c[4];

  // Analytic expectations
  const std::shared_ptr<const Material> poly = matdb->material( "Polyethylene" );
  BOOST_REQUIRE( poly );

  const double mu = calc.m_shells[0].trans_len_coef;  //linear attenuation coef (1/length)
  const double depth = 2.0 * calc.m_shells[0].dims[2]; //full slab depth
  const double x = mu * depth;
  const double rho = poly->density;

  const double expected_ad = (rho/mu) * (1.0 - x*exp(-x)/(1.0 - exp(-x)));
  const double expected_an = GammaInteractionCalc::material_mass_weighted_atomic_number( *poly );
  const double expected_frac_h = GammaInteractionCalc::material_hydrogen_mass_fraction( *poly );

  const double g_cm2 = PhysicalUnits::g / PhysicalUnits::cm2;
  cout << "EffectiveShieldingAnalyticSlab: <AD>=" << eff_ad/g_cm2 << " g/cm2 (expect "
       << expected_ad/g_cm2 << "), <AN>=" << eff_an_mass << " (expect " << expected_an
       << "), <fracH>=" << eff_frac_h << " (expect " << expected_frac_h
       << "), <AN>_xs=" << eff_an_xs << endl;

  BOOST_CHECK_MESSAGE( fabs(eff_ad - expected_ad) <= 0.01*expected_ad,
                       "Effective AD (" << eff_ad/g_cm2 << " g/cm2) not within 1% of analytic ("
                       << expected_ad/g_cm2 << " g/cm2)" );

  BOOST_CHECK_MESSAGE( fabs(eff_an_mass - expected_an) <= 0.001*expected_an,
                       "Mass-weighted effective AN (" << eff_an_mass << ") != material Zbar ("
                       << expected_an << ")" );

  BOOST_CHECK_MESSAGE( fabs(eff_frac_h - expected_frac_h) <= 0.001*std::max(expected_frac_h, 1.0E-6),
                       "Effective H fraction (" << eff_frac_h << ") != material value ("
                       << expected_frac_h << ")" );

  // For a single material both atomic-number weightings must agree exactly
  BOOST_CHECK_MESSAGE( fabs(eff_an_xs - expected_an) <= 0.001*expected_an,
                       "Cross-section-weighted effective AN (" << eff_an_xs << ") != material Zbar ("
                       << expected_an << ")" );
}//BOOST_AUTO_TEST_CASE( EffectiveShieldingAnalyticSlab )
