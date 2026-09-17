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

/** End-to-end validation of the Create DRF path against CeeLo Monte-Carlo truth.

 A synthetic 3"x3" NaI in an Al can is "characterized": peak efficiencies are generated from CeeLo
 MC (the truth table below, regenerated with INTERSPEC_REGEN_MAKEDRF_MC_TRUTH=1 - it runs minutes
 of MC and prints the table to paste back in), fit through MakeDrfCalc/MakeDrfFit exactly as the
 GUI does, and the resulting DRF is checked at held-out positions and through an activity fit.

 Measured 2026-09-13 (see `report_geometry_vs_flat_disk`, which prints the tables):

 The flat-disk model gets the DISTANCE SCALING wrong, because it puts the whole detector at one
 plane (face + setback) while the real interaction depth grows with energy - fitted from this MC
 truth, an effective depth of 0.1 cm at 59.5 keV rising to 2.8 cm at 1408 keV.  Since the error is
 energy dependent it cannot be absorbed into the fitted curve:
   - 50 -> 200 cm the disk under-scales by 0.2% at 59.5 keV, 7.4% at 662 keV, 8.9% at 1408 keV;
   - 50 -> 25 cm it over-scales by 0.8% / 9.5% / 9.8% at those energies.
 So a flat-disk DRF fit to sources at MIXED distances cannot describe them at once (noise-free
 chi2 11.3 vs 0.48 for the geometry kernel, and the fit residuals split by distance: +4.3% at
 25 cm, -1.8% at 50 cm, -5.2% at 100 cm, against +0.1/-1.0/+0.5% with the kernel).  Even fit to a
 SINGLE distance it is only right at that distance (0.6% mean error at 50 cm, but 5.7% at 25 cm
 and 3.7% at 100 cm), where the kernel holds ~0.7-1.0% everywhere.

 Accuracy of the finished DRF against MC, mean/max over the 11 characterization energies:
   geometry  : 0.7/1.7% at 25 cm, 0.8/1.2% at 50 cm, 0.9/1.8% at 100 cm,
               1.0/2.6% at 15 cm, 0.9/2.2% at 200 cm, 3.0/5.5% at 30 cm 30 degrees off axis
   flat disk : 4.3/8.5%, 1.9/3.5%, 5.1/9.2%, 11.5/19.5%, 6.6/11.3%, 5.8/14.1%
 The off-axis 3.0/5.5% is the angle-flat eta of a measured-curve transfer (see
 external_libs/CeeLo/src/io/EfficiencyTransfer.h), which a Monte-Carlo characterization removes.

 Uncertainty coverage of the geometry DRF against these model errors: no held-out point exceeds
 2 sigma anywhere, worst |error|/sigma is 1.5 (at 200 cm) and 1.3-1.4 off axis - i.e. the reported
 uncertainty is honest, and mildly conservative on axis.

 Activity-fit uncertainty (Ba-133, five peaks, counting statistics negligible; measured 2026-09-13
 by `act_fit_uncertainty_includes_drf` / `act_fit_pulls_calibrated`).  The DRF's covariance is
 now the response's own per-query sigma budget with its model envelopes as common modes
 (ceelo::DetectorResponse::frac_covariance); before, those envelopes were patched onto the
 covariance DIAGONAL and averaged down over the peaks:
   reported activity sigma, before -> after:  15 cm on axis 3.12% -> 3.95%,  30 cm at 30 deg
   2.26% -> 2.36%,  100 cm on axis 2.21% -> 2.24% (far field unchanged, as it must be).
 Pulls over 50 replicas that re-draw the calibration AND the measurement (mean / SD / RMS):
   15 cm 0.16 / 0.41 / 0.44,  30 deg -0.98 / 0.74 / 1.23,  100 cm 0.75 / 0.68 / 1.01;
   without the DRF term the RMS is 30 / 24 / 7.  The transfer's model error is one-signed across
   the fit's energies in every regime (mean DRF/MC-1 over the replicas: 15 cm -0.8..+0.4%, 30 deg
   +1.3..+2.5%, 100 cm -1.2..-2.3%), which is what treating it as a common mode assumes.  The far
   field carries a genuine ~2% curve-form bias covered at 0.75 sigma; 15 cm is over-covered
   (envelope ~2x the error for this detector); 30 deg is covered at ~1 sigma only because the
   curve covariance and floor make up for the 0.75% off-axis envelope term.
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_MakeDrfEndToEnd_suite
#include <boost/test/included/unit_test.hpp>

#include <map>
#include <cmath>
#include <deque>
#include <random>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

#include <Eigen/Core>

#include "Minuit2/MnUserParameters.h"

#include "SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "efficiency/EfficiencyCalculator.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"

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


/** The 3"x3" NaI in a 0.5 mm Al can - the same detector test_ShieldingSourceFitCalc's
 make_synthetic_nai_drf builds. */
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


// The characterization energies (the calibration sources' strong lines) and geometries
const vector<double> sm_energies = { 59.5, 81.0, 121.8, 276.4, 302.9, 356.0, 383.8,
                                     661.7, 1173.2, 1332.5, 1408.0 };
struct Geom { double d_cm, theta_deg; };
const vector<Geom> sm_calib_geoms = { {25.0, 0.0}, {50.0, 0.0}, {100.0, 0.0} };
const vector<Geom> sm_heldout_geoms = { {15.0, 0.0}, {30.0, 30.0}, {200.0, 0.0} };


/** One CeeLo MC truth value: absolute FEP efficiency for a point source `d_cm` from the detector
 face along `theta_deg` from the axis. */
struct TruthRow { double energy, d_cm, theta_deg, eff, sigma; };

/** Generated by running this test with INTERSPEC_REGEN_MAKEDRF_MC_TRUTH=1 (prints rows to paste
 here and to makedrf_mc_truth.txt in the working directory).  Target MC sigma 0.4%; the 66 points
 cost ~550 s of CPU (2026-09-13, arm64 Mac, Debug build).
 */
const vector<TruthRow> sm_truth = {
//MC_TRUTH_BEGIN
  { 59.5, 25.0, 0.0, 4.835929e-03, 1.128e-05 },
  { 81.0, 25.0, 0.0, 5.124911e-03, 1.026e-05 },
  { 121.8, 25.0, 0.0, 5.238093e-03, 9.780e-06 },
  { 276.4, 25.0, 0.0, 4.414799e-03, 1.130e-05 },
  { 302.9, 25.0, 0.0, 4.238055e-03, 1.162e-05 },
  { 356.0, 25.0, 0.0, 3.853358e-03, 1.125e-05 },
  { 383.8, 25.0, 0.0, 3.655719e-03, 1.068e-05 },
  { 661.7, 25.0, 0.0, 2.373879e-03, 8.305e-06 },
  { 1173.2, 25.0, 0.0, 1.475292e-03, 5.477e-06 },
  { 1332.5, 25.0, 0.0, 1.325720e-03, 4.960e-06 },
  { 1408.0, 25.0, 0.0, 1.263978e-03, 4.697e-06 },
  { 59.5, 50.0, 0.0, 1.234243e-03, 2.829e-06 },
  { 81.0, 50.0, 0.0, 1.309077e-03, 2.555e-06 },
  { 121.8, 50.0, 0.0, 1.343922e-03, 2.401e-06 },
  { 276.4, 50.0, 0.0, 1.183933e-03, 2.718e-06 },
  { 302.9, 50.0, 0.0, 1.135052e-03, 2.830e-06 },
  { 356.0, 50.0, 0.0, 1.040380e-03, 2.775e-06 },
  { 383.8, 50.0, 0.0, 9.922699e-04, 2.833e-06 },
  { 661.7, 50.0, 0.0, 6.586915e-04, 2.222e-06 },
  { 1173.2, 50.0, 0.0, 4.102047e-04, 1.498e-06 },
  { 1332.5, 50.0, 0.0, 3.694937e-04, 1.353e-06 },
  { 1408.0, 50.0, 0.0, 3.515910e-04, 1.303e-06 },
  { 59.5, 100.0, 0.0, 3.105217e-04, 7.066e-07 },
  { 81.0, 100.0, 0.0, 3.285661e-04, 6.401e-07 },
  { 121.8, 100.0, 0.0, 3.394743e-04, 6.607e-07 },
  { 276.4, 100.0, 0.0, 3.047285e-04, 7.247e-07 },
  { 302.9, 100.0, 0.0, 2.945623e-04, 6.877e-07 },
  { 356.0, 100.0, 0.0, 2.723522e-04, 7.336e-07 },
  { 383.8, 100.0, 0.0, 2.610652e-04, 6.959e-07 },
  { 661.7, 100.0, 0.0, 1.740252e-04, 5.858e-07 },
  { 1173.2, 100.0, 0.0, 1.094973e-04, 3.925e-07 },
  { 1332.5, 100.0, 0.0, 9.749780e-05, 3.606e-07 },
  { 1408.0, 100.0, 0.0, 9.374146e-05, 3.480e-07 },
  { 59.5, 15.0, 0.0, 1.302269e-02, 3.040e-05 },
  { 81.0, 15.0, 0.0, 1.368615e-02, 2.810e-05 },
  { 121.8, 15.0, 0.0, 1.387969e-02, 2.733e-05 },
  { 276.4, 15.0, 0.0, 1.116865e-02, 2.936e-05 },
  { 302.9, 15.0, 0.0, 1.065995e-02, 3.000e-05 },
  { 356.0, 15.0, 0.0, 9.615983e-03, 2.891e-05 },
  { 383.8, 15.0, 0.0, 9.064341e-03, 2.749e-05 },
  { 661.7, 15.0, 0.0, 5.797731e-03, 2.041e-05 },
  { 1173.2, 15.0, 0.0, 3.572173e-03, 1.336e-05 },
  { 1332.5, 15.0, 0.0, 3.169509e-03, 1.190e-05 },
  { 1408.0, 15.0, 0.0, 3.042251e-03, 1.152e-05 },
  { 59.5, 30.0, 30.0, 3.806234e-03, 1.331e-05 },
  { 81.0, 30.0, 30.0, 4.109175e-03, 1.401e-05 },
  { 121.8, 30.0, 30.0, 4.242952e-03, 1.462e-05 },
  { 276.4, 30.0, 30.0, 3.464972e-03, 1.222e-05 },
  { 302.9, 30.0, 30.0, 3.289982e-03, 1.171e-05 },
  { 356.0, 30.0, 30.0, 2.943284e-03, 1.075e-05 },
  { 383.8, 30.0, 30.0, 2.785062e-03, 1.007e-05 },
  { 661.7, 30.0, 30.0, 1.765134e-03, 6.747e-06 },
  { 1173.2, 30.0, 30.0, 1.093431e-03, 4.240e-06 },
  { 1332.5, 30.0, 30.0, 9.834297e-04, 3.823e-06 },
  { 1408.0, 30.0, 30.0, 9.340867e-04, 3.653e-06 },
  { 59.5, 200.0, 0.0, 7.763911e-05, 1.770e-07 },
  { 81.0, 200.0, 0.0, 8.239264e-05, 1.594e-07 },
  { 121.8, 200.0, 0.0, 8.525552e-05, 1.636e-07 },
  { 276.4, 200.0, 0.0, 7.757259e-05, 1.772e-07 },
  { 302.9, 200.0, 0.0, 7.531320e-05, 1.841e-07 },
  { 356.0, 200.0, 0.0, 6.916914e-05, 1.817e-07 },
  { 383.8, 200.0, 0.0, 6.619737e-05, 1.868e-07 },
  { 661.7, 200.0, 0.0, 4.462952e-05, 1.469e-07 },
  { 1173.2, 200.0, 0.0, 2.793281e-05, 1.014e-07 },
  { 1332.5, 200.0, 0.0, 2.509254e-05, 9.107e-08 },
  { 1408.0, 200.0, 0.0, 2.420909e-05, 8.802e-08 },
//MC_TRUTH_END
};


const TruthRow *find_truth( const double energy, const double d_cm, const double theta_deg )
{
  for( const TruthRow &r : sm_truth )
  {
    if( (fabs(r.energy - energy) < 1.0e-3*energy) && (fabs(r.d_cm - d_cm) < 1.0e-3*d_cm)
        && (fabs(r.theta_deg - theta_deg) < 1.0e-3) )
      return &r;
  }
  return nullptr;
}//find_truth(...)


bool have_truth()
{
  for( const Geom &g : sm_calib_geoms )
    for( const double e : sm_energies )
      if( !find_truth( e, g.d_cm, g.theta_deg ) )
        return false;
  for( const Geom &g : sm_heldout_geoms )
    for( const double e : sm_energies )
      if( !find_truth( e, g.d_cm, g.theta_deg ) )
        return false;
  return true;
}//have_truth()


double truth_eff( const double energy, const double d_cm, const double theta_deg )
{
  const TruthRow * const r = find_truth( energy, d_cm, theta_deg );
  BOOST_REQUIRE_MESSAGE( r, "No MC truth at " << energy << " keV, " << d_cm << " cm, " << theta_deg << " deg" );
  return r->eff;
}


/** A calibration source: nuclide, its lines, and where it was measured. */
struct CalSource
{
  const char *nuclide;
  vector<double> energies;
  double d_cm;
  double dist_uncert_cm;
};

/** Which distances the sources were measured at. */
enum class Scenario { SameDistance, Mixed, MixedWithDistUncert };

vector<CalSource> calibration_sources( const Scenario scenario )
{
  const bool mixed = (scenario != Scenario::SameDistance);
  const double du = (scenario == Scenario::MixedWithDistUncert) ? 1.0 : 0.0;
  return {
    { "Am241", { 59.5 },                              mixed ? 25.0 : 50.0, du },
    { "Ba133", { 81.0, 276.4, 302.9, 356.0, 383.8 },  50.0,                0.0 },
    { "Cs137", { 661.7 },                             mixed ? 100.0 : 50.0, 0.0 },
    { "Co60",  { 1173.2, 1332.5 },                    mixed ? 25.0 : 50.0, du },
    { "Eu152", { 121.8, 1408.0 },                     mixed ? 100.0 : 50.0, 0.0 },
  };
}//calibration_sources(...)


/** Builds the measured points the Create DRF tool would have (absolute efficiency at each
 source's own distance, 3% certificate uncertainty per source, Poisson scatter from `gen` when
 `noise`), exactly as MakeDrf::handleSourcesUpdates records them. */
MeasuredDrfPoints make_measured_points( const Scenario scenario, const bool noise, std::mt19937 &gen )
{
  const double n_emitted = 1.0E7;  //gammas emitted per line during the measurement
  const float cert = 0.03f;
  std::normal_distribution<double> normal( 0.0, 1.0 );

  vector<MeasuredEffPoint> pts;
  vector<MeasuredSourceInfo> srcs;
  int idx = 0;
  for( const CalSource &src : calibration_sources(scenario) )
  {
    const string key = string(src.nuclide) + "#" + std::to_string( idx++ );
    // The source's certificate error: one common factor for all of its peaks
    const double cert_factor = 1.0 + (noise ? cert*normal(gen) : 0.0);
    for( const double energy : src.energies )
    {
      const double eff = truth_eff( energy, src.d_cm, 0.0 );
      const double expected_area = n_emitted * eff * cert_factor;
      const double stat = 1.0 / sqrt( expected_area );
      const double area = expected_area * (1.0 + (noise ? stat*normal(gen) : 0.0));

      MeasuredEffPoint p;
      p.energy = static_cast<float>( energy );
      p.efficiency = static_cast<float>( area / n_emitted );
      p.fracStatUncert = static_cast<float>( stat );
      p.fracCertUncert = cert;
      p.sourceKey = key;
      p.distance = static_cast<float>( src.d_cm * PhysicalUnits::cm );
      p.distanceUncert = static_cast<float>( src.dist_uncert_cm * PhysicalUnits::cm );
      p.peakArea = static_cast<float>( area );
      p.peakAreaUncert = static_cast<float>( sqrt(area) );
      p.liveTime = 600.0f;
      pts.push_back( p );
    }//for( const double energy : src.energies )

    MeasuredSourceInfo info;
    info.sourceKey = key;
    info.nuclide = src.nuclide;
    info.activity = 10.0 * PhysicalUnits::microCi;
    info.fracActivityUncert = cert;
    info.distance = static_cast<float>( src.d_cm * PhysicalUnits::cm );
    info.distanceUncert = static_cast<float>( src.dist_uncert_cm * PhysicalUnits::cm );
    srcs.push_back( info );
  }//for( const CalSource &src : calibration_sources(scenario) )

  MeasuredDrfPoints answer;
  answer.setPoints( pts );
  answer.setSources( srcs );
  return answer;
}//make_measured_points(...)


MakeDrfCalc::GeometryChoice full_geometry()
{
  MakeDrfCalc::GeometryChoice geom;
  geom.geometry = make_shared<const ceelo::GeometryDescriptor>( synthetic_nai_descriptor() );
  geom.diameter = 2.0 * geom.geometry->transverse_half_extent() * PhysicalUnits::cm;
  geom.setback = geom.geometry->endcap_front_offset_cm() * PhysicalUnits::cm;
  return geom;
}

MakeDrfCalc::GeometryChoice flat_disk()
{
  MakeDrfCalc::GeometryChoice geom;
  geom.diameter = 2.0 * 3.81 * PhysicalUnits::cm;
  geom.setback = 0.05 * PhysicalUnits::cm;
  return geom;
}

const int sm_num_coefs = 5;

/** Fits the efficiency equation the way MakeDrf does (MeV, `sm_num_coefs` terms). */
MakeDrfFit::EffFitResult fit_points( const MeasuredDrfPoints &points, const MakeDrfCalc::GeometryChoice &geom )
{
  const vector<MakeDrfFit::EffFitPoint> effpts = MakeDrfCalc::intrinsicFitPoints( points, geom, true );
  BOOST_REQUIRE_EQUAL( effpts.size(), points.points().size() );
  return MakeDrfFit::performEfficiencyFit( effpts, sm_num_coefs );
}

MakeDrfCalc::FitResults fit_results( const MakeDrfFit::EffFitResult &eff )
{
  MakeDrfCalc::FitResults fr;
  fr.eff = eff;
  fr.effInMeV = true;
  fr.lowerEnergy = static_cast<float>( sm_energies.front() );
  fr.upperEnergy = static_cast<float>( sm_energies.back() );
  return fr;
}

shared_ptr<DetectorPeakResponse> make_drf( const MeasuredDrfPoints &points, const MakeDrfCalc::GeometryChoice &geom )
{
  const MakeDrfFit::EffFitResult eff = fit_points( points, geom );
  return MakeDrfCalc::assembleDrf( "e2e", "end-to-end test", points, geom, fit_results(eff), nullptr );
}

/** Fractional 1-sigma of the fitted intrinsic curve at `energy_keV` from the coefficient covariance. */
double curve_frac_sigma( const MakeDrfFit::EffFitResult &fit, const double energy_keV )
{
  const size_t n = fit.coefs.size();
  if( fit.covRowMajor.size() != n*n )
    return 0.0;
  const double x = log( energy_keV / 1000.0 );
  vector<double> j( n );
  for( size_t i = 0; i < n; ++i )
    j[i] = pow( x, static_cast<double>(i) );
  double var = 0.0;
  for( size_t a = 0; a < n; ++a )
    for( size_t b = 0; b < n; ++b )
      var += j[a] * fit.covRowMajor[a*n + b] * j[b];
  return sqrt( std::max( 0.0, var ) );
}

double intrinsic_of_fit( const MakeDrfFit::EffFitResult &fit, const double energy_keV )
{
  return DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy_keV / 1000.0, fit.coefs.data(), fit.coefs.size() );
}


/** Checks a DRF against MC truth at the held-out geometries: within max(2 sigma, tol).

 Tolerances measured 2026-09-13 for this NaI with the measured-curve (angle-flat) transfer:
 far field on axis agrees to well under 1%; 15 cm on axis (3.9 crystal radii) the transfer
 over-predicts by up to ~3.5% at the Ba133 energies; 30 degrees off axis it over-predicts by up to
 ~7% at 662 keV - the angle-flat eta table a curve transfer is built on (see
 external_libs/CeeLo/src/io/EfficiencyTransfer.h), which a Monte-Carlo characterization removes.  These are the limits a user should
 know about, not fit defects, so the checks pin them rather than hide them.
 */
void check_heldout( const shared_ptr<DetectorPeakResponse> &drf, const char *label )
{
  for( const Geom &g : sm_heldout_geoms )
  {
    double tol = 0.015;                 //far field, on axis
    if( g.theta_deg > 0.0 ) tol = 0.08; //30 degrees off axis: angle-flat transfer
    else if( g.d_cm < 20.0 ) tol = 0.05; //15 cm: near-field transfer
    for( const double e : sm_energies )
    {
      const double truth = truth_eff( e, g.d_cm, g.theta_deg );
      const DetectorPeakResponse::EffEval ev = drf->fepEfficiencyEval( static_cast<float>(e),
                            g.theta_deg * M_PI / 180.0, 0.0, g.d_cm * PhysicalUnits::cm );
      BOOST_REQUIRE_GT( ev.value, 0.0 );
      const double frac_sig = (ev.value > 0.0) ? (ev.sigma / ev.value) : 0.0;
      const double dev = ev.value / truth - 1.0;
      BOOST_TEST_MESSAGE( label << ": " << e << " keV at " << g.d_cm << " cm, " << g.theta_deg
                          << " deg: DRF/MC-1 = " << 100.0*dev << "%, reported sigma " << 100.0*frac_sig << "%" );
      BOOST_CHECK_MESSAGE( fabs(dev) < std::max( 2.0*frac_sig, tol ),
        label << ": " << e << " keV at " << g.d_cm << " cm, " << g.theta_deg << " deg: DRF/MC-1 = "
        << 100.0*dev << "%, reported sigma " << 100.0*frac_sig << "%" );
    }
  }
}//check_heldout(...)


/** The Ba133 point-source activity fit scaffold at (d_cm, theta_deg), truth areas from MC. */
struct ActFit
{
  double activity = 0.0, uncert = 0.0;
  bool ok = false;
};

bool find_gamma_transition( const SandiaDecay::Nuclide * const parent, const double energy,
                            const SandiaDecay::Transition *&transition, int &particle_index )
{
  transition = nullptr;
  particle_index = -1;
  for( const SandiaDecay::Nuclide *nuc : parent->descendants() )
  {
    for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
    {
      for( size_t prod_index = 0; prod_index < trans->products.size(); ++prod_index )
      {
        const SandiaDecay::RadParticle &product = trans->products[prod_index];
        if( (product.type == SandiaDecay::GammaParticle) && (fabs(product.energy - energy) < 0.5) )
        {
          transition = trans;
          particle_index = static_cast<int>( prod_index );
          return true;
        }
      }
    }
  }
  return false;
}//find_gamma_transition(...)


ActFit fit_ba133_activity( const shared_ptr<DetectorPeakResponse> &drf, const double d_cm,
                           const double theta_deg, const bool account_for_drf_uncert,
                           std::mt19937 *gen )
{
  ActFit answer;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nuc = db->nuclide( "Ba133" );
  BOOST_REQUIRE( nuc );

  const double true_activity = 10.0 * PhysicalUnits::microCi;
  const double age = 5.0 * PhysicalUnits::year;
  const float live_time = 600.0f;

  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, true_activity, age );
  const vector<SandiaDecay::EnergyRatePair> rates = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );

  std::normal_distribution<double> normal( 0.0, 1.0 );
  deque<shared_ptr<const PeakDef>> peaks;
  for( const double energy : { 81.0, 276.4, 302.9, 356.0, 383.8 } )
  {
    double rate = 0.0;  //gammas per second
    for( const SandiaDecay::EnergyRatePair &r : rates )
      if( fabs(r.energy - energy) < 0.5 )
        rate += r.numPerSecond;
    BOOST_REQUIRE_GT( rate, 0.0 );

    const double eff = truth_eff( energy, d_cm, theta_deg );
    double area = rate * live_time * eff;
    if( gen )
      area *= (1.0 + normal(*gen)/sqrt(area));

    auto peak = make_shared<PeakDef>();
    peak->setMean( energy );
    // Narrow enough that the fit's photopeak clustering cannot fold the 79.6 keV line into the
    //  81 keV peak (the truth area is the one line's), which would bias the activity ~1.5% low
    peak->setSigma( 0.004 * energy );
    peak->setPeakArea( area );
    peak->setPeakAreaUncert( sqrt(area) );
    const SandiaDecay::Transition *transition = nullptr;
    int particle_index = -1;
    BOOST_REQUIRE( find_gamma_transition( nuc, energy, transition, particle_index ) );
    peak->setNuclearTransition( nuc, transition, particle_index, PeakDef::SourceGammaType::NormalGamma );
    peak->useForShieldingSourceFit( true );
    peaks.push_back( peak );
  }//for( energies )

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  ShieldingSourceFitCalc::SourceFitDef src;
  src.nuclide = nuc;
  src.activity = 1.0 * PhysicalUnits::microCi;  //start away from truth
  src.fitActivity = true;
  src.age = age;
  src.fitAge = false;
  src.ageDefiningNuc = nullptr;
  src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;
  options.account_for_drf_uncert = account_for_drf_uncert;

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput input;
  const double theta = theta_deg * M_PI / 180.0;
  input.config.distance = d_cm * cos(theta) * PhysicalUnits::cm;   //axial
  input.config.source_offsets[0] = d_cm * sin(theta) * PhysicalUnits::cm;  //transverse
  input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  input.config.sources = { src };
  input.config.options = options;
  input.detector = drf;
  input.foreground = foreground;
  input.foreground_peaks = peaks;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( input );
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;
  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  auto progress_fcn = [](){};
  auto finished_fcn = [](){};
  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );

  answer.ok = (results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final)
              && (results->fit_src_info.size() == 1);
  if( answer.ok )
  {
    answer.activity = results->fit_src_info[0].activity;
    answer.uncert = results->fit_src_info[0].activityUncertainty.value_or( 0.0 );
  }
  return answer;
}//fit_ba133_activity(...)

/** Field-by-field comparison of two point sets, tolerant of float printing (relative 1e-6), that
 names the first difference. */
std::string points_difference( const MeasuredDrfPoints &a, const MeasuredDrfPoints &b )
{
  if( a.points().size() != b.points().size() )
    return "point count " + std::to_string(a.points().size()) + " vs " + std::to_string(b.points().size());
  auto close = []( const double x, const double y ){ return fabs(x - y) <= 1.0e-6*std::max(fabs(x),fabs(y)) + 1.0e-12; };
  for( size_t i = 0; i < a.points().size(); ++i )
  {
    const MeasuredEffPoint &p = a.points()[i], &q = b.points()[i];
    const std::pair<const char *,std::pair<double,double>> fields[] = {
      {"energy",{p.energy,q.energy}}, {"efficiency",{p.efficiency,q.efficiency}},
      {"fracStatUncert",{p.fracStatUncert,q.fracStatUncert}}, {"fracCertUncert",{p.fracCertUncert,q.fracCertUncert}},
      {"distance",{p.distance,q.distance}}, {"distanceUncert",{p.distanceUncert,q.distanceUncert}},
      {"peakArea",{p.peakArea,q.peakArea}}, {"peakAreaUncert",{p.peakAreaUncert,q.peakAreaUncert}},
      {"liveTime",{p.liveTime,q.liveTime}}, {"bkgPeakArea",{p.bkgPeakArea,q.bkgPeakArea}},
      {"bkgPeakAreaUncert",{p.bkgPeakAreaUncert,q.bkgPeakAreaUncert}} };
    for( const auto &f : fields )
      if( !close( f.second.first, f.second.second ) )
        return std::string("point ") + std::to_string(i) + " " + f.first + ": " + std::to_string(f.second.first) + " vs " + std::to_string(f.second.second);
    if( (p.sourceKey != q.sourceKey) || (p.fileName != q.fileName) || (p.sampleNumbers != q.sampleNumbers) )
      return "point " + std::to_string(i) + " strings differ";
  }
  if( a.sources().size() != b.sources().size() )
    return "source count differs";
  for( size_t i = 0; i < a.sources().size(); ++i )
  {
    const MeasuredSourceInfo &p = a.sources()[i], &q = b.sources()[i];
    if( (p.sourceKey != q.sourceKey) || (p.nuclide != q.nuclide) || !close(p.activity,q.activity)
        || !close(p.fracActivityUncert,q.fracActivityUncert) || !close(p.age,q.age)
        || !close(p.distance,q.distance) || !close(p.distanceUncert,q.distanceUncert) )
      return "source " + std::to_string(i) + " differs";
  }
  return "";
}//points_difference(...)
}//namespace


/** Regenerates the MC truth table when INTERSPEC_REGEN_MAKEDRF_MC_TRUTH is set (minutes of MC);
 otherwise just reports whether the table is complete. */
BOOST_AUTO_TEST_CASE( mc_truth_table )
{
  set_data_dir();

  if( !std::getenv( "INTERSPEC_REGEN_MAKEDRF_MC_TRUTH" ) )
  {
    BOOST_TEST_MESSAGE( "MC truth table " << (have_truth() ? "complete" : "MISSING - run with INTERSPEC_REGEN_MAKEDRF_MC_TRUTH=1") );
    BOOST_CHECK_MESSAGE( have_truth(), "MC truth table is incomplete; regenerate it" );
    return;
  }

  const ceelo::GeometryDescriptor gd = synthetic_nai_descriptor();
  std::vector<std::unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator calc;
  ceelo::ResponseGenerator::configure_calculator( calc, gd, owned );

  vector<Geom> geoms = sm_calib_geoms;
  geoms.insert( end(geoms), begin(sm_heldout_geoms), end(sm_heldout_geoms) );

  ofstream out( "makedrf_mc_truth.txt" );
  for( const Geom &g : geoms )
  {
    calc.set_point_source( CeeLoUtils::sourcePositionFromFace( gd, g.theta_deg*M_PI/180.0, 0.0, g.d_cm ) );
    for( const double e : sm_energies )
    {
      ceelo::SimulationConfig cfg;
      cfg.energy_keV = e;
      cfg.termination.target_fep_rel_precision = 0.004;
      cfg.termination.max_events = 40000000;
      cfg.num_threads = 4;
      cfg.seed = 12345;
      const ceelo::EfficiencyResult r = calc.compute( cfg );
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "  { %.1f, %.1f, %.1f, %.6e, %.3e },", e, g.d_cm, g.theta_deg,
                r.full_energy_peak_efficiency, r.fep_uncertainty );
      cout << buffer << "   // " << r.num_events_simulated << " events, " << r.cpu_time_seconds << " s CPU" << endl;
      out << buffer << endl;
    }
  }//for( const Geom &g : geoms )
}//mc_truth_table


BOOST_AUTO_TEST_CASE( fit_same_distance_recovers_truth )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1001 );
  const MeasuredDrfPoints points = make_measured_points( Scenario::SameDistance, true, gen );
  const shared_ptr<DetectorPeakResponse> drf = make_drf( points, full_geometry() );
  BOOST_REQUIRE( drf && drf->isValid() );
  BOOST_CHECK( drf->ceeloResponse() );
  BOOST_CHECK( drf->geometry() );
  BOOST_REQUIRE( drf->measuredPoints() );
  BOOST_CHECK_EQUAL( drf->measuredPoints()->sources().size(), 5u );

  check_heldout( drf, "same distance" );
}//fit_same_distance_recovers_truth


BOOST_AUTO_TEST_CASE( fit_mixed_distance_recovers_truth )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1002 );
  const MeasuredDrfPoints mixed = make_measured_points( Scenario::Mixed, true, gen );
  const MeasuredDrfPoints same = make_measured_points( Scenario::SameDistance, true, gen );

  const MakeDrfFit::EffFitResult fit_mixed = fit_points( mixed, full_geometry() );
  const MakeDrfFit::EffFitResult fit_same = fit_points( same, full_geometry() );

  // The kernel conversion puts 25/50/100 cm points on one curve: the two fits agree within their
  //  (statistical) uncertainties everywhere in range.
  for( const double e : { 70.0, 150.0, 300.0, 600.0, 1000.0, 1300.0 } )
  {
    const double a = intrinsic_of_fit( fit_mixed, e ), b = intrinsic_of_fit( fit_same, e );
    const double sig = sqrt( pow(curve_frac_sigma(fit_mixed,e),2) + pow(curve_frac_sigma(fit_same,e),2) );
    BOOST_CHECK_MESSAGE( fabs(a/b - 1.0) < std::max( 2.5*sig, 0.01 ),
                         e << " keV: mixed/same-distance intrinsic = " << a/b << " (sigma " << sig << ")" );
  }

  const shared_ptr<DetectorPeakResponse> drf = MakeDrfCalc::assembleDrf( "mixed", "", mixed, full_geometry(),
                                                                          fit_results(fit_mixed), nullptr );
  check_heldout( drf, "mixed distance" );
}//fit_mixed_distance_recovers_truth


BOOST_AUTO_TEST_CASE( fit_flat_disk_mixed_distance_is_worse )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  // Noise-free points: whatever scatter the fit sees is the models, not the datas.
  std::mt19937 gen( 1003 );
  const MeasuredDrfPoints mixed = make_measured_points( Scenario::Mixed, false, gen );
  const MakeDrfFit::EffFitResult with_geom = fit_points( mixed, full_geometry() );
  const MakeDrfFit::EffFitResult with_disk = fit_points( mixed, flat_disk() );

  BOOST_TEST_MESSAGE( "chi2: geometry " << with_geom.chi2 << ", flat disk " << with_disk.chi2 );
  BOOST_CHECK_GT( with_disk.chi2, 1.5*with_geom.chi2 );
}//fit_flat_disk_mixed_distance_is_worse


BOOST_AUTO_TEST_CASE( fit_uncertainty_is_calibrated )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1004 );
  const MakeDrfCalc::GeometryChoice geom = full_geometry();
  const MeasuredDrfPoints exact = make_measured_points( Scenario::Mixed, false, gen );
  const MakeDrfFit::EffFitResult fit_exact = fit_points( exact, geom );

  // Pulls of the fitted curve about the noise-free curve over many replicas, each with fresh
  //  Poisson scatter and fresh per-source certificate common modes - so the reported covariance
  //  (stat diagonal + per-source blocks) must reproduce both
  const vector<double> test_energies = { 200.0, 500.0, 1000.0 };
  const int nrep = 200;
  vector<vector<double>> pulls( test_energies.size() );
  for( int rep = 0; rep < nrep; ++rep )
  {
    const MeasuredDrfPoints pts = make_measured_points( Scenario::Mixed, true, gen );
    const MakeDrfFit::EffFitResult fit = fit_points( pts, geom );
    for( size_t i = 0; i < test_energies.size(); ++i )
    {
      const double e = test_energies[i];
      const double sig = curve_frac_sigma( fit, e );
      BOOST_REQUIRE_GT( sig, 0.0 );
      pulls[i].push_back( (intrinsic_of_fit(fit,e) / intrinsic_of_fit(fit_exact,e) - 1.0) / sig );
    }
  }//for( replicas )

  for( size_t i = 0; i < test_energies.size(); ++i )
  {
    double mean = 0.0, var = 0.0;
    for( const double p : pulls[i] ) mean += p / nrep;
    for( const double p : pulls[i] ) var += (p - mean)*(p - mean) / (nrep - 1);
    const double sd = sqrt( var );
    BOOST_TEST_MESSAGE( test_energies[i] << " keV: pull mean " << mean << ", SD " << sd );
    BOOST_CHECK_MESSAGE( fabs(mean) < 0.3, test_energies[i] << " keV pull mean " << mean );
    BOOST_CHECK_MESSAGE( (sd > 0.75) && (sd < 1.3), test_energies[i] << " keV pull SD " << sd );
  }
}//fit_uncertainty_is_calibrated


BOOST_AUTO_TEST_CASE( distance_uncert_inflates_covariance )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1005 );
  const MakeDrfCalc::GeometryChoice geom = full_geometry();
  const MeasuredDrfPoints without = make_measured_points( Scenario::Mixed, false, gen );
  const MeasuredDrfPoints with = make_measured_points( Scenario::MixedWithDistUncert, false, gen );

  const vector<MakeDrfFit::EffFitPoint> pts_without = MakeDrfCalc::intrinsicFitPoints( without, geom, true );
  const vector<MakeDrfFit::EffFitPoint> pts_with = MakeDrfCalc::intrinsicFitPoints( with, geom, true );
  BOOST_REQUIRE_EQUAL( pts_with.size(), pts_without.size() );
  for( size_t i = 0; i < pts_with.size(); ++i )
  {
    BOOST_CHECK_EQUAL( pts_without[i].fracDistUncert, 0.0f );
    BOOST_CHECK_CLOSE( pts_with[i].efficiency, pts_without[i].efficiency, 1.0e-4 );
    const bool at_25cm = (fabs(with.points()[i].distance/PhysicalUnits::cm - 25.0) < 0.01);
    if( at_25cm )
      BOOST_CHECK_CLOSE( pts_with[i].fracDistUncert, 2.0*1.0/25.0, 15.0 );  //~ inverse square: 2 sigma_d/d
    else
      BOOST_CHECK_EQUAL( pts_with[i].fracDistUncert, 0.0f );
  }

  const MakeDrfFit::EffFitResult fit_without = MakeDrfFit::performEfficiencyFit( pts_without, sm_num_coefs );
  const MakeDrfFit::EffFitResult fit_with = MakeDrfFit::performEfficiencyFit( pts_with, sm_num_coefs );
  for( const double e : { 60.0, 1200.0 } )
    BOOST_CHECK_GT( curve_frac_sigma( fit_with, e ), 1.2*curve_frac_sigma( fit_without, e ) );
}//distance_uncert_inflates_covariance


BOOST_AUTO_TEST_CASE( drf_xml_roundtrip_preserves_refit )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1006 );
  const MakeDrfCalc::GeometryChoice geom = full_geometry();
  const MeasuredDrfPoints points = make_measured_points( Scenario::Mixed, true, gen );
  const MakeDrfFit::EffFitResult fit = fit_points( points, geom );
  const shared_ptr<DetectorPeakResponse> drf = MakeDrfCalc::assembleDrf( "rt", "", points, geom, fit_results(fit), nullptr );

  auto restored = make_shared<DetectorPeakResponse>();
  restored->setDrfExtraFromXmlString( drf->drfExtraToXmlString() );
  BOOST_REQUIRE( restored->measuredPoints() );
  BOOST_CHECK_MESSAGE( points_difference( *restored->measuredPoints(), points ).empty(),
                       points_difference( *restored->measuredPoints(), points ) );
  BOOST_REQUIRE( restored->ceeloResponse() );

  // Refit from what the DRF stored alone
  MakeDrfCalc::GeometryChoice geom2;
  geom2.geometry = make_shared<const ceelo::GeometryDescriptor>( restored->ceeloResponse()->descriptor );
  geom2.diameter = geom.diameter;
  const MakeDrfFit::EffFitResult refit = fit_points( *restored->measuredPoints(), geom2 );
  BOOST_REQUIRE_EQUAL( refit.coefs.size(), fit.coefs.size() );
  for( size_t i = 0; i < fit.coefs.size(); ++i )
  {
    BOOST_CHECK_CLOSE( refit.coefs[i], fit.coefs[i], 1.0e-2 );
    BOOST_CHECK_CLOSE( refit.uncerts[i], fit.uncerts[i], 1.0 );
  }
}//drf_xml_roundtrip_preserves_refit


BOOST_AUTO_TEST_CASE( act_fit_near_field_and_off_axis_recover_truth )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1007 );
  const MeasuredDrfPoints points = make_measured_points( Scenario::Mixed, true, gen );
  const shared_ptr<DetectorPeakResponse> drf = make_drf( points, full_geometry() );
  BOOST_REQUIRE( drf && drf->ceeloResponse() );

  const double true_activity = 10.0 * PhysicalUnits::microCi;
  for( const Geom &g : { Geom{15.0, 0.0}, Geom{30.0, 30.0} } )
  {
    const ActFit fit = fit_ba133_activity( drf, g.d_cm, g.theta_deg, true, nullptr );
    BOOST_REQUIRE_MESSAGE( fit.ok, "activity fit failed at " << g.d_cm << " cm, " << g.theta_deg << " deg" );
    const double dev = fit.activity / true_activity - 1.0;
    BOOST_TEST_MESSAGE( g.d_cm << " cm, " << g.theta_deg << " deg: activity/truth-1 = " << 100.0*dev
                        << "%, reported " << 100.0*fit.uncert/fit.activity << "%" );
    // Same transfer limits as check_heldout: ~3.5% at 15 cm on axis, ~5-7% at 30 degrees
    const double tol = (g.theta_deg > 0.0) ? 0.08 : 0.05;
    BOOST_CHECK_MESSAGE( fabs(dev) < tol, "activity off by " << 100.0*dev << "%" );
    BOOST_CHECK_MESSAGE( fabs(fit.activity - true_activity) < 2.5*fit.uncert + 0.01*true_activity,
                         "activity not within 2.5 sigma of truth" );
  }
}//act_fit_near_field_and_off_axis_recover_truth


BOOST_AUTO_TEST_CASE( act_fit_uncertainty_includes_drf )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 1008 );
  const MeasuredDrfPoints points = make_measured_points( Scenario::Mixed, true, gen );
  const shared_ptr<DetectorPeakResponse> drf = make_drf( points, full_geometry() );
  const MeasuredDrfPoints points_du = make_measured_points( Scenario::MixedWithDistUncert, true, gen );
  const shared_ptr<DetectorPeakResponse> drf_du = make_drf( points_du, full_geometry() );

  // Three regimes: near field on axis, off axis, and far field on axis (which only the curve's
  //  own covariance should touch).  Noise-free areas, so the reported sigma is the whole story.
  double far_drf_part = 0.0, near_drf_part = 0.0, off_drf_part = 0.0;
  for( const Geom &g : { Geom{15.0, 0.0}, Geom{30.0, 30.0}, Geom{100.0, 0.0} } )
  {
    const ActFit with = fit_ba133_activity( drf, g.d_cm, g.theta_deg, true, nullptr );
    const ActFit without = fit_ba133_activity( drf, g.d_cm, g.theta_deg, false, nullptr );
    BOOST_REQUIRE( with.ok && without.ok );

    // The DRF's own fractional uncertainty shows up in the activity uncertainty at roughly its size
    const double drf_part = sqrt( std::max( 0.0, pow(with.uncert/with.activity,2) - pow(without.uncert/without.activity,2) ) );
    BOOST_TEST_MESSAGE( g.d_cm << " cm, " << g.theta_deg << " deg: activity sigma with DRF uncert "
                        << 100.0*with.uncert/with.activity << "%, without "
                        << 100.0*without.uncert/without.activity << "%, DRF part " << 100.0*drf_part << "%" );
    BOOST_CHECK_GT( with.uncert, 1.2*without.uncert );

    if( g.theta_deg > 0.0 )
      off_drf_part = drf_part;
    else if( g.d_cm < 20.0 )
      near_drf_part = drf_part;
    else
      far_drf_part = drf_part;
  }//for( each regime )

  // Far field on axis: the 3% certificate common modes (through the curve covariance)
  //  dominate, and far_drf_part measures ~1.8%.  These bounds guard that the DRF term
  //  stays non-trivial and does not run away; they are not a calibration - that is
  //  act_fit_pulls_calibrated, checked separately.  fep_far_floor is a fully-correlated
  //  common mode, so it feeds this quantity directly.
  BOOST_CHECK_GT( far_drf_part, 0.015 );
  // Keep the upper bound close to the measurement: a bound that cannot fail is not a test.
  BOOST_CHECK_LT( far_drf_part, 0.03 );
  // Near field and off axis add the transfer's model envelope on top
  BOOST_CHECK_GT( near_drf_part, far_drf_part );
  BOOST_CHECK_GT( off_drf_part, far_drf_part );

  const ActFit with = fit_ba133_activity( drf, 15.0, 0.0, true, nullptr );
  const ActFit with_du = fit_ba133_activity( drf_du, 15.0, 0.0, true, nullptr );
  BOOST_REQUIRE( with.ok && with_du.ok );
  BOOST_TEST_MESSAGE( "15 cm: activity sigma with distance uncert in DRF "
                      << 100.0*with_du.uncert/with_du.activity << "% vs " << 100.0*with.uncert/with.activity << "%" );
  BOOST_CHECK_GT( with_du.uncert, with.uncert );
}//act_fit_uncertainty_includes_drf


/** Are the reported activity uncertainties calibrated against the MC truth?

 Each replica re-draws the CALIBRATION (fresh Poisson scatter and certificate common modes on the
 characterization points, hence a fresh DRF) as well as the measurement, so everything the DRF's
 data-derived covariance claims to cover actually varies: in the far field, where the transfer is
 exact, the pull SD must come out near 1.  The transfer's model error (angle-flat eta off axis, the
 kernel-only near field) is the same in every replica: it cannot widen the SD, it shifts the pull
 MEAN.  So the coverage metric is the RMS pull about the truth, sqrt(mean^2 + SD^2), which the
 reported sigma - data plus model envelope - must keep near 1.  Without the DRF term the pulls are
 pure counting statistics and everything goes uncovered.
 */
BOOST_AUTO_TEST_CASE( act_fit_pulls_calibrated )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  const double true_activity = 10.0 * PhysicalUnits::microCi;
  const int nrep = 50;
  const MakeDrfCalc::GeometryChoice geom = full_geometry();

  // The replica DRFs, shared by every regime and both option settings below
  std::mt19937 cal_gen( 1009 );
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  for( int rep = 0; rep < nrep; ++rep )
  {
    const MeasuredDrfPoints points = make_measured_points( Scenario::Mixed, true, cal_gen );
    drfs.push_back( make_drf( points, geom ) );
    BOOST_REQUIRE( drfs.back() && drfs.back()->ceeloResponse() );
  }

  // Is the transfer's model error one-signed across the energies of a fit (a common mode), as
  //  the covariance treats it?  Report it per regime, averaged over the replica DRFs so the
  //  calibration noise averages out, rather than assume it.
  for( const Geom &g : { Geom{15.0, 0.0}, Geom{30.0, 30.0}, Geom{100.0, 0.0} } )
  {
    std::ostringstream msg;
    msg << g.d_cm << " cm, " << g.theta_deg << " deg: mean DRF/MC-1 (%) at";
    for( const double e : { 81.0, 276.4, 302.9, 356.0, 383.8 } )
    {
      double dev = 0.0;
      for( const shared_ptr<DetectorPeakResponse> &drf : drfs )
      {
        const DetectorPeakResponse::EffEval ev = drf->fepEfficiencyEval( static_cast<float>(e),
                                     g.theta_deg*M_PI/180.0, 0.0, g.d_cm*PhysicalUnits::cm );
        dev += (ev.value/truth_eff(e, g.d_cm, g.theta_deg) - 1.0) / nrep;
      }
      msg << " " << e << ":" << std::fixed << std::setprecision(2) << 100.0*dev;
    }
    BOOST_TEST_MESSAGE( msg.str() );
  }

  for( const Geom &g : { Geom{15.0, 0.0}, Geom{30.0, 30.0}, Geom{100.0, 0.0} } )
  {
    for( const bool with_drf : { true, false } )
    {
      std::mt19937 meas_gen( 2000 + static_cast<unsigned>(g.d_cm) + (with_drf ? 1 : 0) );
      vector<double> pulls;
      double mean_sigma = 0.0;
      for( int rep = 0; rep < nrep; ++rep )
      {
        const ActFit fit = fit_ba133_activity( drfs[rep], g.d_cm, g.theta_deg, with_drf, &meas_gen );
        BOOST_REQUIRE( fit.ok && (fit.uncert > 0.0) );
        pulls.push_back( (fit.activity - true_activity) / fit.uncert );
        mean_sigma += (fit.uncert / fit.activity) / nrep;
      }

      double mean = 0.0, var = 0.0;
      for( const double p : pulls ) mean += p / nrep;
      for( const double p : pulls ) var += (p - mean)*(p - mean) / (nrep - 1);
      const double sd = sqrt( var ), rms = sqrt( mean*mean + var );
      BOOST_TEST_MESSAGE( g.d_cm << " cm, " << g.theta_deg << " deg, " << (with_drf ? "with" : "without")
                          << " DRF uncert: pull mean " << mean << ", SD " << sd << ", RMS " << rms
                          << " (mean reported sigma " << 100.0*mean_sigma << "%)" );

      if( with_drf )
      {
        if( (g.theta_deg == 0.0) && (g.d_cm > 50.0) )
        {
          // Far field: the replicas exercise the data-derived part of the reported sigma (the
          //  curve covariance) but not its model floors, so the SD is sigma_data/sigma_total < 1
          //  rather than 1; and the fitted curve form carries a genuine ~1-2% bias at these
          //  energies (see the header comment), which the sigma must cover.
          BOOST_CHECK_MESSAGE( fabs(mean) < 1.0, "far field pull mean " << mean );
          BOOST_CHECK_MESSAGE( (sd > 0.5) && (sd < 1.3), "far field pull SD " << sd );
        }
        BOOST_CHECK_MESSAGE( rms < 1.3, g.d_cm << " cm, " << g.theta_deg << " deg: RMS pull " << rms );
      // AND a lower bound.  An over-covering envelope is a real defect, not a safe
      //  default - it hides disagreement and makes every activity look consistent - and
      //  without this bound it is invisible here, since an envelope ten times too large
      //  passes `rms < 1.3` trivially.  The three regimes read ~0.8 / 1.0 / 1.2, so 0.6
      //  leaves room for honest variation without admitting a 2x-over-covered envelope.
      BOOST_CHECK_MESSAGE( rms > 0.6, g.d_cm << " cm, " << g.theta_deg
                           << " deg: RMS pull " << rms << " - envelope is OVER-covering" );
      }else
      {
        // Counting statistics alone cover neither the calibration nor the transfer's model error
        BOOST_CHECK_MESSAGE( rms > 1.5, g.d_cm << " cm, " << g.theta_deg << " deg: RMS pull " << rms << " without DRF uncert" );
      }
    }//for( with_drf )
  }//for( each regime )
}//act_fit_pulls_calibrated


/** Reports how well each detector description reproduces the MC, and what the geometry kernel buys
 over the flat disk when the calibration sources sit at different distances.

 Noise-free points, so every number below is model error, not counting statistics.
 */
BOOST_AUTO_TEST_CASE( report_geometry_vs_flat_disk )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 2001 );
  const MeasuredDrfPoints mixed = make_measured_points( Scenario::Mixed, false, gen );
  const MeasuredDrfPoints same = make_measured_points( Scenario::SameDistance, false, gen );

  struct Case { const char *label; MeasuredDrfPoints points; MakeDrfCalc::GeometryChoice geom; };
  vector<Case> cases = {
    { "mixed 25/50/100 cm, geometry", mixed, full_geometry() },
    { "mixed 25/50/100 cm, flat disk", mixed, flat_disk() },
    { "all at 50 cm, geometry",        same,  full_geometry() },
    { "all at 50 cm, flat disk",       same,  flat_disk() },
  };

  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "=== Fit quality (noise-free points; chi2 is pure model mis-fit) ===" );
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  vector<MakeDrfFit::EffFitResult> fits;
  for( const Case &c : cases )
  {
    const MakeDrfFit::EffFitResult fit = fit_points( c.points, c.geom );
    fits.push_back( fit );
    drfs.push_back( MakeDrfCalc::assembleDrf( c.label, "", c.points, c.geom, fit_results(fit), nullptr ) );

    // Residuals the fit itself saw: the intrinsic efficiency each point implies, vs the fitted curve
    const vector<MakeDrfFit::EffFitPoint> pts = MakeDrfCalc::intrinsicFitPoints( c.points, c.geom, true );
    double max_res = 0.0, sum_abs = 0.0;
    for( const MakeDrfFit::EffFitPoint &p : pts )
    {
      const double res = intrinsic_of_fit( fit, 1000.0*p.energy ) / p.efficiency - 1.0;
      max_res = std::max( max_res, fabs(res) );
      sum_abs += fabs(res) / pts.size();
    }
    char buffer[256];
    snprintf( buffer, sizeof(buffer), "%-32s chi2 = %7.3f (dof %d)   fit residual: mean |%.2f%%|, max |%.2f%%|",
              c.label, fit.chi2, fit.dof, 100.0*sum_abs, 100.0*max_res );
    BOOST_TEST_MESSAGE( buffer );
  }//for( const Case &c : cases )

  // How the mixed-distance points scatter about the flat-disk fit, by distance: a flat disk cannot
  //  describe 25 and 100 cm at once, so its residuals split by distance (the geometry kernel's do not).
  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "=== Mixed-distance fit residual, grouped by the distance the source was at ===" );
  for( size_t ci = 0; ci < 2; ++ci )  //the two mixed cases
  {
    const vector<MakeDrfFit::EffFitPoint> pts = MakeDrfCalc::intrinsicFitPoints( cases[ci].points, cases[ci].geom, true );
    BOOST_REQUIRE_EQUAL( pts.size(), cases[ci].points.points().size() );
    map<int,pair<double,int>> by_dist;  //cm -> (sum of residual, count)
    for( size_t i = 0; i < pts.size(); ++i )
    {
      const int d_cm = static_cast<int>( std::lround( cases[ci].points.points()[i].distance / PhysicalUnits::cm ) );
      const double res = intrinsic_of_fit( fits[ci], 1000.0*pts[i].energy ) / pts[i].efficiency - 1.0;
      by_dist[d_cm].first += res;
      by_dist[d_cm].second += 1;
    }
    string line = string(cases[ci].label) + ":";
    for( const auto &d : by_dist )
    {
      char buffer[128];
      snprintf( buffer, sizeof(buffer), "  %d cm: %+.2f%% (n=%d)", d.first,
                100.0*d.second.first/d.second.second, d.second.second );
      line += buffer;
    }
    BOOST_TEST_MESSAGE( line );
  }//for( the two mixed cases )

  // Accuracy of the finished DRF against MC, at every energy and distance - the question a user
  //  actually has ("if I measure at X cm, how wrong is the efficiency?")
  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "=== DRF efficiency vs MC truth (DRF/MC - 1), by source position ===" );
  vector<Geom> all_geoms = sm_calib_geoms;
  all_geoms.insert( end(all_geoms), begin(sm_heldout_geoms), end(sm_heldout_geoms) );

  BOOST_TEST_MESSAGE( "                                  25cm     50cm    100cm  |  15cm    30cm@30deg  200cm" );
  for( size_t ci = 0; ci < cases.size(); ++ci )
  {
    string line;
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%-32s", cases[ci].label );
    line = buffer;
    for( const Geom &g : all_geoms )
    {
      double sum_abs = 0.0, max_abs = 0.0;
      for( const double e : sm_energies )
      {
        const DetectorPeakResponse::EffEval ev = drfs[ci]->fepEfficiencyEval( static_cast<float>(e),
                              g.theta_deg * M_PI / 180.0, 0.0, g.d_cm * PhysicalUnits::cm );
        const double dev = ev.value / truth_eff( e, g.d_cm, g.theta_deg ) - 1.0;
        sum_abs += fabs(dev) / sm_energies.size();
        max_abs = std::max( max_abs, fabs(dev) );
      }
      snprintf( buffer, sizeof(buffer), " %5.2f/%5.2f", 100.0*sum_abs, 100.0*max_abs );
      line += buffer;
    }
    BOOST_TEST_MESSAGE( line );
  }//for( cases )
  BOOST_TEST_MESSAGE( "  (each cell is mean/max |DRF/MC-1| in %, over the 11 characterization energies)" );

  // The claims the numbers above must support
  const double flat_mixed_chi2 = fits[1].chi2, geom_mixed_chi2 = fits[0].chi2;
  BOOST_CHECK_MESSAGE( flat_mixed_chi2 > 5.0*geom_mixed_chi2,
                       "flat-disk mixed-distance chi2 " << flat_mixed_chi2
                       << " should be far worse than the geometry's " << geom_mixed_chi2 );

  // ...and at a calibration distance the geometry DRF must actually be more accurate
  for( const Geom &g : sm_calib_geoms )
  {
    double geom_err = 0.0, flat_err = 0.0;
    for( const double e : sm_energies )
    {
      const double truth = truth_eff( e, g.d_cm, 0.0 );
      geom_err += fabs( drfs[0]->fepEfficiencyEval( static_cast<float>(e), 0.0, 0.0, g.d_cm*PhysicalUnits::cm ).value / truth - 1.0 ) / sm_energies.size();
      flat_err += fabs( drfs[1]->fepEfficiencyEval( static_cast<float>(e), 0.0, 0.0, g.d_cm*PhysicalUnits::cm ).value / truth - 1.0 ) / sm_energies.size();
    }
    BOOST_CHECK_MESSAGE( geom_err < flat_err,
                         g.d_cm << " cm: geometry mean error " << 100.0*geom_err
                         << "% is not better than flat disk " << 100.0*flat_err << "%" );
  }
}//report_geometry_vs_flat_disk


/** `MakeDrfCalc::assembleDrf` re-grounds the Monte-Carlo response it is handed.  Two things must
 hold, and both failed silently once:

  1. The response arrives ALREADY grounded (the MC tool grounds before emitting).  k(E) is the ratio
     of measurement to the UNGROUNDED model, so if the old grounding is not cleared before the model
     efficiency is evaluated, the new k(E) comes out ~1 and the measurement anchoring is discarded.
  2. One save calls assembleDrf several times (store, N42 export, reference sheet) on the SAME
     shared response object, so the operation has to be idempotent.

 A curve-transfer response stands in for an MC one here - it exercises the identical code path
 without minutes of Monte Carlo.
 */
BOOST_AUTO_TEST_CASE( assemble_drf_preserves_mc_grounding )
{
  set_data_dir();
  if( !have_truth() ){ BOOST_TEST_MESSAGE( "skipped: no MC truth" ); return; }

  std::mt19937 gen( 3001 );
  const MeasuredDrfPoints points = make_measured_points( Scenario::Mixed, true, gen );
  const MakeDrfCalc::GeometryChoice geom = full_geometry();
  const MakeDrfFit::EffFitResult fit = fit_points( points, geom );

  // Deliberately DEFORMED: the stand-in's model efficiency is 30% high, so grounding it to the real
  //  points has to do real work (a response anchored on those same points would give k == 1, and
  //  the test would pass whether or not the grounding survived).
  MakeDrfCalc::FitResults biased = fit_results( fit );
  biased.eff.coefs[0] += static_cast<float>( log(1.30) );

  auto make_response = [&]() -> shared_ptr<ceelo::DetectorResponse> {
    const shared_ptr<DetectorPeakResponse> seed
          = MakeDrfCalc::assembleDrf( "seed", "", points, flat_disk(), biased, nullptr );
    // ...and the measured points have to go, or transferAnchorForDrf would anchor on THEM and
    //  hand back a response that already reproduces the measurements (k == 1 again).
    seed->setMeasuredPoints( nullptr );
    const CeeLoUtils::TransferAnchor anchor
          = CeeLoUtils::transferAnchorForDrf( seed, *geom.geometry, 50.0 );
    return CeeLoUtils::makeTransferResponse( *geom.geometry, anchor, ceelo::AnchorCurve{}, "stand-in" );
  };

  // A response that arrives GROUNDED must stay grounded, with the same k(E) after repeated saves.
  {
    shared_ptr<ceelo::DetectorResponse> resp = make_response();
    bool curve_derived = false;
    vector<ceelo::GroundingPoint> pts
          = MakeMcResponseForDrf::groundingPointsForDrf(
                MakeDrfCalc::assembleDrf( "seed", "", points, flat_disk(), fit_results(fit), nullptr ),
                *geom.geometry, curve_derived );
    BOOST_REQUIRE( !pts.empty() );
    ceelo::ResponseGenerator::ground_to_points( *resp, pts, curve_derived );
    BOOST_REQUIRE( !resp->grounding.empty() );

    const vector<double> ln_k_at_generation = resp->grounding.ln_k;
    double max_abs_lnk = 0.0;
    for( const double v : ln_k_at_generation )
      max_abs_lnk = std::max( max_abs_lnk, fabs(v) );
    BOOST_REQUIRE_MESSAGE( max_abs_lnk > 1.0e-4,
                           "the stand-in response should need a non-trivial grounding to be a test" );

    const shared_ptr<DetectorPeakResponse> first
          = MakeDrfCalc::assembleDrf( "a", "", points, geom, fit_results(fit), resp );
    BOOST_REQUIRE( first->ceeloResponse() );
    BOOST_REQUIRE( !first->ceeloResponse()->grounding.empty() );
    const vector<double> ln_k_first = first->ceeloResponse()->grounding.ln_k;

    double max_after = 0.0;
    for( const double v : ln_k_first )
      max_after = std::max( max_after, fabs(v) );
    BOOST_CHECK_MESSAGE( max_after > 1.0e-4,
      "assembleDrf cancelled the grounding: max|ln k| went from " << max_abs_lnk << " to " << max_after );

    // ...and doing it again (N42 export, reference sheet) must not drift
    const shared_ptr<DetectorPeakResponse> second
          = MakeDrfCalc::assembleDrf( "b", "", points, geom, fit_results(fit), resp );
    BOOST_REQUIRE( second->ceeloResponse() );
    const vector<double> ln_k_second = second->ceeloResponse()->grounding.ln_k;
    BOOST_REQUIRE_EQUAL( ln_k_second.size(), ln_k_first.size() );
    for( size_t i = 0; i < ln_k_first.size(); ++i )
      BOOST_CHECK_MESSAGE( fabs(ln_k_first[i] - ln_k_second[i]) < 1.0e-9,
        "assembleDrf is not idempotent: knot " << i << " ln k " << ln_k_first[i] << " -> " << ln_k_second[i] );

    // The grounded response reproduces the measurements it was grounded to
    for( const MeasuredEffPoint &p : points.points() )
    {
      const DetectorPeakResponse::EffEval ev = first->fepEfficiencyEval( p.energy, 0.0, 0.0, p.distance );
      BOOST_CHECK_MESSAGE( fabs( ev.value/p.efficiency - 1.0 ) < 0.05,
        p.energy << " keV at " << p.distance/PhysicalUnits::cm << " cm: grounded response gives "
        << ev.value << " vs measured " << p.efficiency );
    }
  }

  // A response that arrives UNGROUNDED was deliberately left so (the MC tool's "correct to this
  //  detector's measured efficiency" checkbox is off) - assembleDrf must not override that.
  {
    shared_ptr<ceelo::DetectorResponse> resp = make_response();
    BOOST_REQUIRE( resp->grounding.empty() );
    const shared_ptr<DetectorPeakResponse> drf
          = MakeDrfCalc::assembleDrf( "c", "", points, geom, fit_results(fit), resp );
    BOOST_REQUIRE( drf->ceeloResponse() );
    BOOST_CHECK_MESSAGE( drf->ceeloResponse()->grounding.empty(),
                         "assembleDrf grounded a response the user asked not to be grounded" );
  }
}//assemble_drf_preserves_mc_grounding
