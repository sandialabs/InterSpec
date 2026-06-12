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
  const vector<pair<string,int>> cases = {
    { "sph-U-plus-Fe-661keV", 0 },     //sphere radius
    { "cyl-end-U-plus-Fe-661keV", 1 }, //cylinder half-length
    { "rect-U-plus-Fe-661keV", 2 },    //rectangle half-depth
  };

  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();

  for( const pair<string,int> &test_case : cases )
  {
    const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
    for( const VolumetricFixture::VolumetricSrcSpec &fixture : fixtures )
    {
      if( fixture.name == test_case.first )
        spec = &fixture;
    }
    BOOST_REQUIRE( spec );

    const int dim_index = test_case.second;
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
      if( shell_index >= src_shell )
        shell.dims[dim_index].v[0] = 1.0;

      calc_jet.m_shells.push_back( shell );
    }//for( loop over shells )

    GammaInteractionCalc::self_shielding_integration_imp( calc_jet, epsrel, max_evals );

    // Finite-difference reference from the double path
    const double dim_val = calc_dbl.m_shells[src_shell].dims[dim_index];
    const double step = 0.02 * dim_val;

    GammaInteractionCalc::DistributedSrcCalcT<double> calc_plus = calc_dbl;
    GammaInteractionCalc::DistributedSrcCalcT<double> calc_minus = calc_dbl;
    for( size_t shell_index = src_shell; shell_index < calc_dbl.m_shells.size(); ++shell_index )
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


/** TODO: enable once sources can be off the detector axis.

 A sphere at distance d with lateral offset s must give the same answer as an
 on-axis sphere at distance sqrt(d*d + s*s): the shells are rotation-invariant
 and (currently) the detector response uses only the line-of-sight distance
 (`fractionalSolidAngle( 2*detRad, dist + setback )` in eval_spherical), with no
 detector-face orientation term.  NOTE: this equivalence breaks once the
 detector response gains an angular dependence - at that point this test should
 compare against an independent 3D reference instead.

 Planned check: 'sph-U-plus-Fe-661keV' fixture at 300 cm with a 100 cm offset,
 versus the same fixture on-axis at sqrt(10)*100 cm.
 */
BOOST_AUTO_TEST_CASE( OffAxisSphereEquivalence, * boost::unit_test::disabled() )
{
  BOOST_FAIL( "Off-axis source offsets are not implemented yet" );
}//BOOST_AUTO_TEST_CASE( OffAxisSphereEquivalence )


/** TODO: enable once effective-AN/AD/H accumulation during integration exists.

 For the 'rect-poly-slab-661keV' fixture (polyethylene slab, half-depth
 D/2 = 2 cm, at 3 m), rays to the detector are nearly normal to the slab, so
 with z the depth from a source point to the exit surface, and mu the linear
 attenuation coefficient, the attenuation-weighted averages are analytic:

   <AD>    = rho * Int(z*exp(-mu*z))/Int(exp(-mu*z)), z over [0,D]
           = (rho/mu) * (1 - x*exp(-x)/(1 - exp(-x))),  x = mu*D
   <AN>    = mass-weighted Zbar of polyethylene (CH2):
             0.8563*6 + 0.1437*1 = ~5.28
   <fracH> = hydrogen mass fraction of polyethylene = ~0.1437

 (expect agreement to ~1%, the residual coming from ray obliquity)
 */
BOOST_AUTO_TEST_CASE( EffectiveShieldingAnalyticSlab, * boost::unit_test::disabled() )
{
  BOOST_FAIL( "Effective-AN/AD/H accumulation is not implemented yet" );
}//BOOST_AUTO_TEST_CASE( EffectiveShieldingAnalyticSlab )
