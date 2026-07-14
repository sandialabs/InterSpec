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
#include <string>
#include <vector>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ShieldSourcePullTrend_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldSourcePullTrend.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/MassAttenuationTool.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
  // The atomic-number discrimination-power check needs the cross-section data directory; find it
  //  from --datadir (or nearby) once.  Returns true if mass-attenuation data is available.
  bool ensure_static_data()
  {
    static int s_state = -1;   //-1 unknown, 0 no, 1 yes
    if( s_state >= 0 )
      return (s_state == 1);
    string datadir;
    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;
    for( int i = 1; i < argc; ++i )
    {
      const string a = argv[i];
      if( SpecUtils::istarts_with(a, "--datadir=") ) datadir = a.substr(10);
    }
    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data",
                             "/Users/wcjohns/coding/InterSpec/data" } )
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) ){ datadir = d; break; }
    }
    try
    {
      if( !datadir.empty() )
      {
        InterSpec::setStaticDataDirectory( datadir );
        s_state = 1;
        return true;
      }
    }catch( std::exception & ){}
    s_state = 0;
    BOOST_WARN_MESSAGE( false, "No datadir found - atomic-number discrimination tests SKIPPED "
                              "(pass --datadir to exercise them)." );
    return false;
  }

  // Build a (PeakResultPlotInfo, PeakDetail) pair for a nuclide at a given energy/pull.
  //  `significance` (observed counts / uncertainty) feeds the discrimination-power check.
  void add_point( vector<GammaInteractionCalc::PeakResultPlotInfo> &comps,
                  vector<GammaInteractionCalc::PeakDetail> &details,
                  const string &nuclide, const double energy, const double pull,
                  const double significance = 0.0, const double attenFactor = 1.0 )
  {
    GammaInteractionCalc::PeakResultPlotInfo c;
    c.energy = energy;
    c.numSigmaOff = pull;
    c.observedOverExpected = 1.0;
    c.observedOverExpectedUncert = 0.1;
    c.observedCounts = significance*significance;  //so counts/uncert == significance
    c.observedUncert = significance;
    comps.push_back( c );

    GammaInteractionCalc::PeakDetail d;
    d.energy = energy;
    d.numSigmaOff = pull;
    d.assignedNuclide = nuclide;
    d.m_totalShieldAttenFactor = attenFactor;  //model optical depth for the discrimination-power check
    details.push_back( d );
  }

  // Realistic total-shield attenuation fraction for a generic (Z, areal-density) shield.
  double atten_factor( const double z, const double ad_gcm2, const double energy )
  {
    const double mu = MassAttenuation::massAttenuationCoefficientFracAN( static_cast<float>(z),
                                                                         static_cast<float>(energy) );
    const double ad = ad_gcm2 * PhysicalUnits::g / PhysicalUnits::cm2;
    return std::exp( -mu * ad );
  }

  // No-shielding-fit config (activity-only fit): both amount and Z diagnoses allowed.
  ShieldSourcePullTrend::FitConfigSummary activity_only_config()
  {
    ShieldSourcePullTrend::FitConfigSummary c;
    return c; //all false
  }

  // Physical thickness was fit: amount optimal, only atomic-number diagnosis allowed.
  ShieldSourcePullTrend::FitConfigSummary thickness_fit_config()
  {
    ShieldSourcePullTrend::FitConfigSummary c;
    c.physicalDimensionFit = true;
    return c;
  }

  // A generic (Z, areal-density) shielding, so the discrimination-power check has a model.
  vector<ShieldingSourceFitCalc::ShieldingInfo> generic_shield( const double z, const double ad_gcm2 )
  {
    ShieldingSourceFitCalc::ShieldingInfo s;
    s.m_isGenericMaterial = true;
    s.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
    s.m_dimensions[0] = z;
    s.m_dimensions[1] = ad_gcm2 * PhysicalUnits::g / PhysicalUnits::cm2;
    s.m_dimensions[2] = 0.0;
    return { s };
  }
}//namespace


// A strong quadratic (low-energy hook) should be detected as Quadratic.
BOOST_AUTO_TEST_CASE( QuadraticDetected )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  // chi(E) with a downward low-energy hook (negative curvature in logE), amplitude ~7.
  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const double emin = 60, emax = 1300;
  for( const double E : energies )
  {
    // map logE to [-1,1]
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double chi = -14.0*(x*x) + 7.0; //strong negative curvature (well above threshold)
    add_point( comps, details, "Ba133", E, chi );
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, activity_only_config(), ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::Quadratic) );
  BOOST_CHECK( std::fabs(tr.curvatureT) > 3.5 );
  BOOST_CHECK( tr.curvature < 0.0 );  //negative curvature
  BOOST_CHECK_EQUAL( tr.numPointsUsed, energies.size() );
  BOOST_REQUIRE_EQUAL( tr.nuclideTrends.size(), 1u );
  BOOST_CHECK( !tr.nuclideTrends[0].curve.empty() );
  // No shielding passed -> not discriminable -> a significant curvature maps to "other", not a
  //  Z direction.  Pin that branch explicitly.
  BOOST_CHECK( !tr.anDiscriminable );
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::OtherModelInconsistency) );
}


// A pure linear trend (no curvature) should be Linear, with correct slope sign.
BOOST_AUTO_TEST_CASE( LinearDetected )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const double emin = 60, emax = 1300;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double chi = 9.0*x; //strong positive slope, no curvature (above threshold)
    add_point( comps, details, "Ba133", E, chi );
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, activity_only_config(), ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::Linear) );
  BOOST_CHECK( tr.slopeT > 6.0 );
  // Positive slope, no shielding fit -> "amount too low".
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::ShieldingAmountTooLow) );
}


// Flat (null) pulls should yield no trend and no conclusion.
BOOST_AUTO_TEST_CASE( NullFlatNoConclusion )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const vector<double> pulls    = { 0.4, -0.3, 0.1, -0.2, 0.25, -0.15, 0.05, -0.1 };
  for( size_t i = 0; i < energies.size(); ++i )
    add_point( comps, details, "Ba133", energies[i], pulls[i] );

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, activity_only_config(), ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::None) );
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::NoConclusion) );
}


// Fewer than three points cannot support a trend.
BOOST_AUTO_TEST_CASE( TooFewPoints )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;
  add_point( comps, details, "Cs137", 300, 2.0 );
  add_point( comps, details, "Cs137", 600, -2.0 );

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, activity_only_config(), ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::None) );
}


// Joint fit: two nuclides with different intercepts but a shared slope should be recovered,
//  with each nuclide getting its own (different) intercept.
BOOST_AUTO_TEST_CASE( JointSharedSlope )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 80, 160, 320, 640, 1000 };
  const double emin = 80, emax = 1000;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    add_point( comps, details, "NucA", E, 9.0*x + 2.0 ); //intercept +2
    add_point( comps, details, "NucB", E, 9.0*x - 2.0 ); //intercept -2, same slope
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, activity_only_config(), ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK( tr.model != ShieldSourcePullTrend::TrendModel::None );
  BOOST_CHECK( tr.slopeT > 6.0 );
  BOOST_REQUIRE_EQUAL( tr.nuclideTrends.size(), 2u );
  // Shared slope recovered with correct sign and magnitude (regressor scaled to [-1,1] so the
  //  fitted slope should be ~ +9); linear model, no significant curvature.
  BOOST_CHECK( tr.slope > 0.0 );
  BOOST_CHECK( std::fabs( tr.slope - 9.0 ) < 1.0 );
  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::Linear) );
  // Each nuclide keeps its own intercept: one near +2, one near -2 (order not guaranteed).
  const double i0 = tr.nuclideTrends[0].intercept, i1 = tr.nuclideTrends[1].intercept;
  const double hi = std::max(i0,i1), lo = std::min(i0,i1);
  BOOST_CHECK( std::fabs( hi - 2.0 ) < 0.8 );
  BOOST_CHECK( std::fabs( lo + 2.0 ) < 0.8 );
}


// When shielding thickness was fit (amount optimal), a curvature should be read as an
//  atomic-number problem, not an amount problem - provided the peaks have the power to resolve Z
//  (low-energy peaks present, and the cross-section data is available).
BOOST_AUTO_TEST_CASE( AmountFitGivesAtomicNumberConclusion )
{
  if( !ensure_static_data() )
    return;  //discrimination-power check needs mass-attenuation data; skip if unavailable

  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const double emin = 60, emax = 1300;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double chi = -14.0*(x*x) + 7.0; //strong negative curvature
    add_point( comps, details, "Ba133", E, chi, 100.0, atten_factor(26.0,20.0,E) );
  }

  // A real generic shield (Z=26, 20 g/cm2) with low-energy peaks -> discriminable.
  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, generic_shield(26.0, 20.0), thickness_fit_config(),
      ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::Quadratic) );
  BOOST_CHECK( tr.anDiscriminable );
  // Amount was fit + concave-down (negative curvature) + power -> effective atomic number too LOW.
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::AtomicNumberTooLow) );
}


// Same curvature, but only high-energy peaks (>=400 keV, Compton-flat): no power to resolve Z, so
//  the atomic-number call must be SUPPRESSED (reported as other/inconclusive, not a Z direction).
BOOST_AUTO_TEST_CASE( NoLowEnergyPeaksSuppressesAtomicNumberCall )
{
  if( !ensure_static_data() )
    return;

  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 700, 850, 1000, 1150, 1300, 1450 };
  const double emin = 700, emax = 1450;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double chi = -14.0*(x*x) + 7.0; //same strong curvature shape
    add_point( comps, details, "Ba133", E, chi, 60.0, atten_factor(26.0,20.0,E) );
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, generic_shield(26.0, 20.0), thickness_fit_config(),
      ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK( !tr.anDiscriminable );  //no low-energy leverage
  // Not asserted as an atomic-number direction (it cannot be a Z signature here).
  BOOST_CHECK( tr.conclusion != ShieldSourcePullTrend::TrendConclusion::AtomicNumberTooLow );
  BOOST_CHECK( tr.conclusion != ShieldSourcePullTrend::TrendConclusion::AtomicNumberTooHigh );
}


// Heavy shield (Pb, Z=82): the reference Z error (82*1.5=123) exceeds the Z<=98 attenuation table,
//  so the power calc must clamp/flip rather than throw and zero the power - otherwise the
//  atomic-number diagnosis would be silently disabled for lead/tungsten/uranium shields.
BOOST_AUTO_TEST_CASE( HeavyShieldAtomicNumberStillDiscriminable )
{
  if( !ensure_static_data() )
    return;

  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 100, 150, 250, 400, 600, 900, 1300 };  //avoid Pb K-edge 88 keV
  const double emin = 100, emax = 1300;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double chi = -14.0*(x*x) + 7.0;  //strong concave-down
    add_point( comps, details, "Ba133", E, chi, 100.0, atten_factor(82.0,20.0,E) );
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, generic_shield(82.0, 20.0), thickness_fit_config(),
      ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK( tr.anDiscriminationT > 0.0 );  //power was computed, not thrown away
  BOOST_CHECK( tr.anDiscriminable );
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::AtomicNumberTooLow) );
}


// Both generic atomic-number and areal-density fit: a residual trend can only be flagged as
//  "other inconsistency" (never amount/Z).
BOOST_AUTO_TEST_CASE( AnAndAdFitGivesOtherOrNone )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const double emin = 60, emax = 1300;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    add_point( comps, details, "Ba133", E, 9.0*x );
  }

  ShieldSourcePullTrend::FitConfigSummary config;
  config.genericAnFit = true;
  config.genericAdFit = true;

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, config, ShieldSourcePullTrend::TrendBasis::LogEnergy );

  const bool ok = (tr.conclusion == ShieldSourcePullTrend::TrendConclusion::OtherModelInconsistency)
               || (tr.conclusion == ShieldSourcePullTrend::TrendConclusion::NoConclusion);
  BOOST_CHECK( ok );
}


// A self-attenuating / trace source present suppresses all conclusions (curves still drawn).
BOOST_AUTO_TEST_CASE( SelfAttenSuppressesConclusion )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const double emin = 60, emax = 1300;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    add_point( comps, details, "U235", E, 9.0*x );
  }

  ShieldSourcePullTrend::FitConfigSummary config;
  config.hasSelfAttenOrTraceSrc = true;

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, config, ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::NoConclusion) );
  // But the trend curve is still computed.
  BOOST_CHECK( tr.model != ShieldSourcePullTrend::TrendModel::None );
}


// Concave-UP (positive curvature) with amount fit + power -> effective atomic number too HIGH.
BOOST_AUTO_TEST_CASE( AmountFitPositiveCurvatureIsAtomicNumberTooHigh )
{
  if( !ensure_static_data() )
    return;

  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  const vector<double> energies = { 60, 100, 150, 250, 400, 600, 900, 1300 };
  const double emin = 60, emax = 1300;
  for( const double E : energies )
  {
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double chi = 14.0*(x*x) - 7.0; //strong POSITIVE curvature (concave up)
    add_point( comps, details, "Ba133", E, chi, 100.0, atten_factor(26.0,20.0,E) );
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, generic_shield(26.0, 20.0), thickness_fit_config(),
      ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK_EQUAL( static_cast<int>(tr.model),
                     static_cast<int>(ShieldSourcePullTrend::TrendModel::Quadratic) );
  BOOST_CHECK( tr.anDiscriminable );
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::AtomicNumberTooHigh) );
}


// Residual-variance scaling: a curvature-shaped signal buried under large over-dispersion
//  (a poorly-converged-fit surrogate) must NOT produce a conclusion, because the residual RMS
//  scaling deflates the t-statistic.  This is what makes the atomic-number direction reliable.
BOOST_AUTO_TEST_CASE( OverDispersionSuppressesConclusion )
{
  vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
  vector<GammaInteractionCalc::PeakDetail> details;

  // Many points, a mild quadratic shape, but with big alternating scatter so the pulls are
  //  strongly over-dispersed about any smooth trend (redChi2 >> 1).
  const vector<double> energies = { 60,80,100,130,170,220,300,400,550,750,1000,1300 };
  const double emin = 60, emax = 1300;
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const double E = energies[i];
    const double x = 2.0*(std::log(E) - std::log(emin))/(std::log(emax)-std::log(emin)) - 1.0;
    const double scatter = (i % 2 == 0) ? 18.0 : -18.0; //huge alternating residual
    add_point( comps, details, "Ba133", E, -3.0*(x*x) + 1.5 + scatter );
  }

  const ShieldSourcePullTrend::TrendResult tr = ShieldSourcePullTrend::fit_pull_trend(
      comps, details, {}, activity_only_config(), ShieldSourcePullTrend::TrendBasis::LogEnergy );

  BOOST_CHECK( tr.redChi2 > 5.0 ); //genuinely over-dispersed
  BOOST_CHECK_EQUAL( static_cast<int>(tr.conclusion),
                     static_cast<int>(ShieldSourcePullTrend::TrendConclusion::NoConclusion) );
}


// Plain-text report sentence (for non-GUI reports): matches the conclusion + confidence, empty
//  when there is no conclusion.
BOOST_AUTO_TEST_CASE( ReportTextMatchesConclusion )
{
  ShieldSourcePullTrend::TrendResult tr;

  tr.conclusion = ShieldSourcePullTrend::TrendConclusion::NoConclusion;
  BOOST_CHECK( ShieldSourcePullTrend::conclusion_report_text(tr).empty() );

  tr.conclusion = ShieldSourcePullTrend::TrendConclusion::AtomicNumberTooLow;
  tr.confidence = ShieldSourcePullTrend::TrendConfidence::Moderate;
  const string an = ShieldSourcePullTrend::conclusion_report_text(tr);
  BOOST_CHECK( an.find("probably indicate") != string::npos );
  BOOST_CHECK( an.find("atomic number is too low") != string::npos );

  tr.conclusion = ShieldSourcePullTrend::TrendConclusion::ShieldingAmountTooHigh;
  tr.confidence = ShieldSourcePullTrend::TrendConfidence::Strong;
  const string amt = ShieldSourcePullTrend::conclusion_report_text(tr);
  BOOST_CHECK( amt.find("strongly indicate") != string::npos );
  BOOST_CHECK( amt.find("too much modeled shielding") != string::npos );

  tr.conclusion = ShieldSourcePullTrend::TrendConclusion::OtherModelInconsistency;
  const string other = ShieldSourcePullTrend::conclusion_report_text(tr);
  BOOST_CHECK( other.find("despite fitting shielding") != string::npos );
}
