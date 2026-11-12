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
#include <utility>
#include <functional>

#ifdef _WIN32
// NOMINMAX must come before Windows.h, or it defines function-like min(a,b)/max(a,b) macros that
//  turn any later `std::max( a, b )` into MSVC error C2589.  This TU includes Windows.h itself,
//  ahead of anything (e.g. a boost/asio chain) that would otherwise define NOMINMAX for us - which
//  is what made it the one test to hit this.  See the "Windows min/max macros" note in CLAUDE.md.
#define NOMINMAX
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ShieldSourcePullTrend_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldSourcePullTrend.h"
#include "InterSpec/ShieldingSourceFitCalc.h"

using namespace std;
using namespace boost::unit_test;

using ShieldSourcePullTrend::TrendResult;
using ShieldSourcePullTrend::TrendConclusion;
using ShieldSourcePullTrend::TrendConfidence;

namespace
{
  // The likelihood-ratio tests need the cross-section data (and the decay database for nuclide
  //  pointers); find the data directory from --datadir (or nearby) once.
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
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) ){ datadir = d; break; }
    }
    try
    {
      if( !datadir.empty() )
      {
        InterSpec::setStaticDataDirectory( datadir );
        s_state = DecayDataBaseServer::database() ? 1 : 0;
        return (s_state == 1);
      }
    }catch( std::exception & ){}
    s_state = 0;
    BOOST_WARN_MESSAGE( false, "No datadir found - pull-trend tests SKIPPED (pass --datadir)." );
    return false;
  }


  const double g_cm2 = PhysicalUnits::g / PhysicalUnits::cm2;

  /** Optical depth of a generic (Z, areal density in g/cm2) shield at an energy (keV). */
  double tau_generic( const double z, const double ad_gcm2, const double energy )
  {
    if( (z <= 0.0) || (ad_gcm2 <= 0.0) )
      return 0.0;
    return GammaInteractionCalc::transmition_coefficient_generic( z, ad_gcm2*g_cm2, energy );
  }


  /** One gamma line of a test source: energy (keV) and detected counts unshielded. */
  struct Line { double energy, counts; };

  /** Eu152-like lines; `scale` sets the statistics (strongest peak ~3e4 counts at scale 1). */
  vector<Line> eu152_lines( const double scale = 1.0, const double min_energy = 0.0 )
  {
    const vector<pair<double,double>> lines = { {121.8, 28.5}, {244.7, 7.6}, {344.3, 26.6}, {411.1, 2.2},
                    {444.0, 3.1}, {778.9, 13.0}, {867.4, 4.2}, {964.1, 14.5}, {1085.8, 10.1},
                    {1112.1, 13.7}, {1408.0, 20.9} };
    vector<Line> out;
    for( const auto &l : lines )
      if( l.first >= min_energy )
        out.push_back( { l.first, scale * 1500.0 * l.second * std::pow( l.first/100.0, -0.8 ) } );
    return out;
  }

  vector<Line> ba133_lines( const double scale = 1.0 )
  {
    vector<Line> out;
    for( const auto &l : vector<pair<double,double>>{ {81.0, 34.0}, {276.4, 7.2}, {302.9, 18.3}, {356.0, 62.1}, {383.8, 8.9} } )
      out.push_back( { l.first, scale * 1200.0 * l.second * std::pow( l.first/100.0, -0.8 ) } );
    return out;
  }


  /** A generic shield (Z, AD in g/cm2); none when z <= 0. */
  struct Shield
  {
    double z = 0.0, ad = 0.0;
  };


  /** A synthetic activity/shielding fit result. */
  struct Case
  {
    vector<GammaInteractionCalc::PeakResultPlotInfo> comps;
    vector<GammaInteractionCalc::PeakDetail> details;
    vector<ShieldingSourceFitCalc::ShieldingInfo> shields;
    ShieldSourcePullTrend::FitConfigSummary config;
    Shield model;   //as fit
  };


  /** Noise-free observed counts through `truth` (optionally distorted by `obs_factor(E)`), and model
   counts through `model` with every nuclide's activity at its least-squares optimum.  With
   `fit_amount`, the model shields areal density is also optimized, as an amount-fit
   activity/shielding fit would.  Uncertainties are sqrt(counts) times `sigma_factor`.
   */
  Case make_case( const vector<pair<string,vector<Line>>> &sources, const Shield &truth, Shield model,
                  const bool fit_amount, const bool fit_an = false,
                  const std::function<double(double)> &obs_factor = nullptr,
                  const double sigma_factor = 1.0 )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE( db );

    struct Pk { string nuc; double energy, unit, obs, sigma, model; };
    vector<Pk> pks;
    for( const auto &src : sources )
      for( const Line &l : src.second )
      {
        Pk p;
        p.nuc = src.first;
        p.energy = l.energy;
        p.unit = l.counts;
        p.obs = l.counts * std::exp( -tau_generic( truth.z, truth.ad, l.energy ) );
        if( obs_factor )
          p.obs *= obs_factor( l.energy );
        p.sigma = sigma_factor * std::sqrt( std::max( 1.0, p.obs ) );
        p.model = 0.0;
        pks.push_back( p );
      }

    // Model counts (activities at their optimum, per nuclide) and chi2, for a model areal density.
    const auto evaluate = [&]( const double ad ) -> double {
      double chi2 = 0.0;
      for( const auto &src : sources )
      {
        double num = 0.0, den = 0.0;
        for( const Pk &p : pks )
        {
          if( p.nuc != src.first )
            continue;
          const double u = p.unit * std::exp( -tau_generic( model.z, ad, p.energy ) );
          num += p.obs*u/(p.sigma*p.sigma);
          den += u*u/(p.sigma*p.sigma);
        }
        const double a = (den > 0.0) ? (num/den) : 0.0;
        for( Pk &p : pks )
        {
          if( p.nuc != src.first )
            continue;
          p.model = a * p.unit * std::exp( -tau_generic( model.z, ad, p.energy ) );
          chi2 += std::pow( (p.obs - p.model)/p.sigma, 2 );
        }
      }
      return chi2;
    };

    if( fit_amount && (model.z > 0.0) )
    {
      double best_ad = 0.0, best = evaluate( 0.0 );
      for( double ad = 0.01; ad < 400.0; ad *= 1.02 )
      {
        const double c2 = evaluate( ad );
        if( c2 < best ){ best = c2; best_ad = ad; }
      }
      model.ad = best_ad;
    }
    evaluate( model.ad );

    Case c;
    c.model = model;
    for( const Pk &p : pks )
    {
      GammaInteractionCalc::PeakResultPlotInfo comp;
      comp.energy = p.energy;
      comp.observedCounts = p.obs;
      comp.observedUncert = p.sigma;
      comp.expectedCounts = p.model;
      comp.numSigmaOff = (p.obs - p.model)/p.sigma;
      comp.observedOverExpected = (p.model > 0.0) ? (p.obs/p.model) : 0.0;
      comp.observedOverExpectedUncert = (p.model > 0.0) ? (p.sigma/p.model) : 0.0;
      c.comps.push_back( comp );

      GammaInteractionCalc::PeakDetail d;
      d.energy = d.decayParticleEnergy = p.energy;
      d.assignedNuclide = p.nuc;
      d.m_totalShieldAttenFactor = std::exp( -tau_generic( model.z, model.ad, p.energy ) );
      GammaInteractionCalc::PeakDetailSrc src;
      src.nuclide = db->nuclide( p.nuc );
      BOOST_REQUIRE( src.nuclide );
      src.energy = p.energy;
      src.countsAtSource = p.unit;
      src.modelContribToPeak = p.model;
      d.m_sources.push_back( src );
      c.details.push_back( d );
    }

    if( model.z > 0.0 )
    {
      ShieldingSourceFitCalc::ShieldingInfo s;
      s.m_isGenericMaterial = true;
      s.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
      s.m_dimensions[0] = model.z;
      s.m_dimensions[1] = model.ad * g_cm2;
      s.m_dimensions[2] = 0.0;
      s.m_fitDimensions[0] = fit_an;
      s.m_fitDimensions[1] = fit_amount;
      s.m_fitDimensions[2] = false;
      c.shields.push_back( s );
      c.config.genericAnFit = fit_an;
      c.config.genericAdFit = fit_amount;
    }
    return c;
  }//make_case(...)


  TrendResult run( const Case &c )
  {
    return ShieldSourcePullTrend::fit_pull_trend( c.comps, c.details, c.shields, c.config,
                                                  ShieldSourcePullTrend::TrendBasis::LogEnergy );
  }

  int conc( const TrendResult &r ){ return static_cast<int>( r.conclusion ); }
  int conc( const TrendConclusion c ){ return static_cast<int>( c ); }

  const Shield iron20{ 26.0, 20.0 };
}//namespace


// A correct model (activity-only fit, shielding fixed at the truth) gives no message.
BOOST_AUTO_TEST_CASE( CorrectModelNoConclusion )
{
  if( !ensure_static_data() )
    return;
  const Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, iron20, false );
  const TrendResult r = run( c );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::NoConclusion) );
  BOOST_CHECK( std::fabs( r.amountT ) < 1.0 );
  BOOST_CHECK( r.modelHasShielding );
  BOOST_CHECK_EQUAL( r.numPeaksTested, eu152_lines().size() );
}


// The most common mistake: the source is shielded, the model has none.  The old polynomial test
//  reported "Other" here (with text claiming shielding was fit); now it is "lacks shielding", with
//  an estimate of how much.
BOOST_AUTO_TEST_CASE( MissingShieldingIsLacksShielding )
{
  if( !ensure_static_data() )
    return;
  const Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{}, false );
  const TrendResult r = run( c );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::ModelLacksShielding) );
  BOOST_CHECK( !r.modelHasShielding );
  BOOST_CHECK( r.amountT > 6.0 );
  BOOST_CHECK_EQUAL( static_cast<int>(r.confidence), static_cast<int>(TrendConfidence::Strong) );
  const double ad = r.suggestedArealDensity / g_cm2;
  BOOST_CHECK_MESSAGE( (ad > 10.0) && (ad < 40.0), "Suggested AD " << ad << " g/cm2 (truth 20 of Fe)" );

  const ShieldSourcePullTrend::MessageHint hint = ShieldSourcePullTrend::message_hint( r );
  BOOST_CHECK( hint.key != nullptr );
  const string txt = ShieldSourcePullTrend::conclusion_report_text( r );
  BOOST_CHECK( txt.find( "lacks shielding" ) != string::npos );
  BOOST_CHECK( txt.find( "g/cm2" ) != string::npos );
}


// Too little / too much of the right material, fixed (activity-only fit).
BOOST_AUTO_TEST_CASE( AmountTooLowAndTooHigh )
{
  if( !ensure_static_data() )
    return;

  const TrendResult low = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{26.0, 10.0}, false ) );
  BOOST_CHECK_EQUAL( conc(low), conc(TrendConclusion::ShieldingAmountTooLow) );
  BOOST_CHECK( low.amountT > 3.5 );
  BOOST_CHECK_MESSAGE( std::fabs( low.suggestedArealDensity/g_cm2 - 20.0 ) < 4.0,
                       "Suggested AD " << low.suggestedArealDensity/g_cm2 );

  const TrendResult high = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{26.0, 40.0}, false ) );
  BOOST_CHECK_EQUAL( conc(high), conc(TrendConclusion::ShieldingAmountTooHigh) );
  BOOST_CHECK( high.amountT < -3.5 );
  BOOST_CHECK_MESSAGE( std::fabs( high.suggestedArealDensity/g_cm2 - 20.0 ) < 4.0,
                       "Suggested AD " << high.suggestedArealDensity/g_cm2 );
}


// A bare source modeled behind 5 g/cm2 of lead (fixed): a huge, step-shaped misfit that the old
//  polynomial test, inflating its own uncertainty by the misfit, reported as "atomic number
//  possibly too low" - or nothing at all for 10 g/cm2 of iron.
BOOST_AUTO_TEST_CASE( BareSourceWithShieldFixedIsTooMuch )
{
  if( !ensure_static_data() )
    return;

  for( const Shield &model : { Shield{82.0, 5.0}, Shield{26.0, 10.0} } )
  {
    const TrendResult r = run( make_case( {{ "Ba133", ba133_lines() }}, Shield{}, model, false ) );
    BOOST_CHECK_MESSAGE( conc(r) == conc(TrendConclusion::ShieldingAmountTooHigh),
                         "Z=" << model.z << " AD=" << model.ad << ": conclusion " << conc(r) );
    BOOST_CHECK( r.amountT < -6.0 );
    const ShieldSourcePullTrend::MessageHint hint = ShieldSourcePullTrend::message_hint( r );
    BOOST_CHECK( hint.key && (string(hint.key) == "x2g-trend-hint-no-shield") );
  }
}


// Wrong material with the areal density fit: the amount is optimal, so a residual must be the
//  atomic number - in the right direction.
BOOST_AUTO_TEST_CASE( WrongAtomicNumberWithAmountFit )
{
  if( !ensure_static_data() )
    return;

  const Case low_case = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{13.0, 20.0}, true );
  const TrendResult low = run( low_case );
  BOOST_CHECK_EQUAL( conc(low), conc(TrendConclusion::AtomicNumberTooLow) );
  BOOST_CHECK( low.atomicNumberT > 3.5 );
  BOOST_CHECK( low.anDiscriminable );
  BOOST_CHECK_MESSAGE( (low.suggestedAtomicNumber > 17.0) && (low.suggestedAtomicNumber < 40.0),
                       "Suggested Z " << low.suggestedAtomicNumber );

  const Case high_case = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{82.0, 20.0}, true );
  const TrendResult high = run( high_case );
  BOOST_CHECK_EQUAL( conc(high), conc(TrendConclusion::AtomicNumberTooHigh) );
  BOOST_CHECK( high.atomicNumberT < -3.5 );
}


// Only high-energy peaks (Compton-dominated, nearly Z-independent shape): no power to resolve Z,
//  so a modest Z error must not be called.
BOOST_AUTO_TEST_CASE( NoLowEnergyPeaksNoAtomicNumberCall )
{
  if( !ensure_static_data() )
    return;
  const Case c = make_case( {{ "Eu152", eu152_lines( 1.0, 700.0 ) }}, iron20, Shield{40.0, 20.0}, true );
  const TrendResult r = run( c );
  BOOST_CHECK( !r.anDiscriminable );
  BOOST_CHECK( (r.conclusion != TrendConclusion::AtomicNumberTooLow)
               && (r.conclusion != TrendConclusion::AtomicNumberTooHigh) );
}


// Atomic number AND areal density fit: neither can be blamed; a smooth residual trend (here a
//  15% detector-efficiency shape error) can only be "other".
BOOST_AUTO_TEST_CASE( AnAndAdFitOnlyOther )
{
  if( !ensure_static_data() )
    return;

  const TrendResult clean = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, iron20, true, true ) );
  BOOST_CHECK_EQUAL( conc(clean), conc(TrendConclusion::NoConclusion) );

  const auto eff_error = []( const double energy ){
    const double x = std::log( energy/400.0 );
    return 1.0 + 0.15*x*x;
  };
  const TrendResult bent = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, iron20, true, true, eff_error ) );
  BOOST_CHECK( (bent.conclusion == TrendConclusion::OtherModelInconsistency)
               || (bent.conclusion == TrendConclusion::NoConclusion) );
}


// Self-attenuating / trace sources: the point-source reasoning does not apply.
BOOST_AUTO_TEST_CASE( SelfAttenSuppressesConclusion )
{
  if( !ensure_static_data() )
    return;
  Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{}, false );
  c.config.hasSelfAttenOrTraceSrc = true;
  BOOST_CHECK_EQUAL( conc( run(c) ), conc(TrendConclusion::NoConclusion) );
}


// Several nuclides behind the same shielding are judged together (each with its own activity).
//  Here each nuclide alone spans a narrow energy range, but together they span 80-1408 keV.
BOOST_AUTO_TEST_CASE( MultipleNuclidesShareTheShielding )
{
  if( !ensure_static_data() )
    return;
  vector<Line> eu_high;
  for( const Line &l : eu152_lines() )
    if( l.energy > 700.0 )
      eu_high.push_back( l );
  const Case c = make_case( {{ "Ba133", ba133_lines() }, { "Eu152", eu_high }}, iron20, Shield{}, false );
  const TrendResult r = run( c );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::ModelLacksShielding) );
}


// Real peak fits scatter beyond their statistical uncertainty, and an occasional poor fit is off by
//  5-10%: one high-statistics peak 8% off, on an otherwise correct model, must not produce a
//  shielding call.
BOOST_AUTO_TEST_CASE( SingleOutlierPeakNoShieldingCall )
{
  if( !ensure_static_data() )
    return;
  const auto outlier = []( const double energy ){ return (std::fabs( energy - 344.3 ) < 0.5) ? 1.08 : 1.0; };
  const TrendResult r = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, iron20, false, false, outlier ) );
  BOOST_CHECK_MESSAGE( (r.conclusion == TrendConclusion::NoConclusion)
                       || (r.conclusion == TrendConclusion::OtherModelInconsistency),
                       "conclusion " << conc(r) << ", amountT " << r.amountT );
}


// Two fit peaks assigned the same gamma share one model entry: the second (which gets no source or
//  attenuation detail) must not be used.
BOOST_AUTO_TEST_CASE( DuplicateEnergyPeakIgnored )
{
  if( !ensure_static_data() )
    return;
  Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{}, false );
  const TrendResult ref = run( c );

  GammaInteractionCalc::PeakResultPlotInfo dup_comp = c.comps[2];
  GammaInteractionCalc::PeakDetail dup_detail;
  dup_detail.energy = dup_detail.decayParticleEnergy = dup_comp.energy;
  dup_detail.assignedNuclide = "Eu152";
  c.comps.push_back( dup_comp );
  c.details.push_back( dup_detail );

  const TrendResult r = run( c );
  BOOST_CHECK_EQUAL( r.numPeaksTested, ref.numPeaksTested );
  BOOST_CHECK_CLOSE( r.amountT, ref.amountT, 1.0E-6 );
}


// Too few peaks to test anything.
BOOST_AUTO_TEST_CASE( TooFewPeaks )
{
  if( !ensure_static_data() )
    return;
  vector<Line> two = eu152_lines();
  two.resize( 2 );
  const TrendResult r = run( make_case( {{ "Eu152", two }}, iron20, Shield{}, false ) );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::NoConclusion) );
}


// The cascade-summing caveat: only for close sources, correction off, and cascade nuclides.
BOOST_AUTO_TEST_CASE( SummingCaveat )
{
  if( !ensure_static_data() )
    return;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  ShieldingSourceFitCalc::SourceFitDef co60, cs137;
  co60.nuclide = db->nuclide( "Co60" );
  cs137.nuclide = db->nuclide( "Cs137" );
  BOOST_REQUIRE( co60.nuclide && cs137.nuclide );

  ShieldSourcePullTrend::FitConfigSummary cfg;
  ShieldSourcePullTrend::set_summing_info( cfg, 5.0*PhysicalUnits::cm, false, nullptr, { co60 } );
  BOOST_CHECK( cfg.possibleCascadeSumming );
  ShieldSourcePullTrend::set_summing_info( cfg, 50.0*PhysicalUnits::cm, false, nullptr, { co60 } );
  BOOST_CHECK( !cfg.possibleCascadeSumming );
  ShieldSourcePullTrend::set_summing_info( cfg, 5.0*PhysicalUnits::cm, true, nullptr, { co60 } );
  BOOST_CHECK( !cfg.possibleCascadeSumming );
  ShieldSourcePullTrend::set_summing_info( cfg, 5.0*PhysicalUnits::cm, false, nullptr, { cs137 } );
  BOOST_CHECK( !cfg.possibleCascadeSumming );

  // Coincident gammas only from progeny actually present (not behind Np237's 2 My half-life, nor
  //  Pa234's 0.16% branch); gamma + annihilation photon summing counts (Na22); a low-yield but
  //  level-linked cascade counts (Lu177 208 -> 113 keV, ~0.2 gammas per decay); single-line or
  //  parallel-branch emitters do not (Co57's 122 keV is followed only by the sub-30 keV 14.4 keV line).
  ShieldingSourceFitCalc::SourceFitDef am241, u238, na22, eu152, lu177, cs134, co57, tl201;
  am241.nuclide = db->nuclide( "Am241" );
  am241.age = 20.0*PhysicalUnits::year;
  u238.nuclide = db->nuclide( "U238" );
  u238.age = 1.0E4*PhysicalUnits::year;
  na22.nuclide = db->nuclide( "Na22" );
  eu152.nuclide = db->nuclide( "Eu152" );
  lu177.nuclide = db->nuclide( "Lu177" );
  cs134.nuclide = db->nuclide( "Cs134" );
  co57.nuclide = db->nuclide( "Co57" );
  tl201.nuclide = db->nuclide( "Tl201" );
  BOOST_REQUIRE( am241.nuclide && u238.nuclide && na22.nuclide && eu152.nuclide
                 && lu177.nuclide && cs134.nuclide && co57.nuclide && tl201.nuclide );
  for( const pair<ShieldingSourceFitCalc::SourceFitDef,bool> &src_expect
        : vector<pair<ShieldingSourceFitCalc::SourceFitDef,bool>>{ {am241, false}, {u238, false}, {na22, true},
            {eu152, true}, {lu177, true}, {cs134, true}, {co57, false}, {tl201, false} } )
  {
    ShieldSourcePullTrend::set_summing_info( cfg, 5.0*PhysicalUnits::cm, false, nullptr, { src_expect.first } );
    BOOST_CHECK_MESSAGE( cfg.possibleCascadeSumming == src_expect.second,
                         src_expect.first.nuclide->symbol << ": caveat " << cfg.possibleCascadeSumming );
  }

  // The caveat is carried into the result and the report text.
  Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{}, false );
  c.config.possibleCascadeSumming = true;
  const TrendResult r = run( c );
  BOOST_CHECK( r.possibleCascadeSumming );
  BOOST_CHECK( ShieldSourcePullTrend::conclusion_report_text( r ).find( "summing" ) != string::npos );
}


// The descriptive trend curve's ambiguity guard: a peak whose model counts come mostly from another
//  nuclide is left out of the curve fit (sm_min_nuclide_fraction).
BOOST_AUTO_TEST_CASE( AmbiguousPeakExcludedFromCurve )
{
  if( !ensure_static_data() )
    return;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  Case c = make_case( {{ "Ba133", ba133_lines() }}, Shield{}, Shield{}, false );
  const size_t nclean = c.comps.size();

  GammaInteractionCalc::PeakResultPlotInfo comp = c.comps[0];
  comp.energy = 661.7;
  comp.numSigmaOff = 9.0;
  GammaInteractionCalc::PeakDetail d;
  d.energy = d.decayParticleEnergy = 661.7;
  d.assignedNuclide = "Ba133";
  GammaInteractionCalc::PeakDetailSrc ba_src, cs_src;
  ba_src.nuclide = db->nuclide( "Ba133" );
  ba_src.countsAtSource = 200.0;
  cs_src.nuclide = db->nuclide( "Cs137" );
  cs_src.countsAtSource = 800.0;
  d.m_sources = { ba_src, cs_src };
  c.comps.push_back( comp );
  c.details.push_back( d );

  const TrendResult r = run( c );
  BOOST_CHECK_EQUAL( r.numPointsUsed, nclean );
}


// Plain-text report sentences.
BOOST_AUTO_TEST_CASE( ReportText )
{
  TrendResult tr;
  tr.conclusion = TrendConclusion::NoConclusion;
  BOOST_CHECK( ShieldSourcePullTrend::conclusion_report_text(tr).empty() );

  tr.conclusion = TrendConclusion::ShieldingAmountTooLow;
  tr.confidence = TrendConfidence::Strong;
  tr.suggestedArealDensity = 24.0 * g_cm2;
  const string low = ShieldSourcePullTrend::conclusion_report_text(tr);
  BOOST_CHECK( low.find("strongly") != string::npos );
  BOOST_CHECK( low.find("too little") != string::npos );
  BOOST_CHECK( low.find("24 g/cm2") != string::npos );

  tr.conclusion = TrendConclusion::AtomicNumberTooHigh;
  tr.confidence = TrendConfidence::Tentative;
  tr.suggestedAtomicNumber = 26.0;
  tr.suggestedAtomicNumberResolved = true;
  const string zhigh = ShieldSourcePullTrend::conclusion_report_text(tr);
  BOOST_CHECK( zhigh.find("possibly") != string::npos );
  BOOST_CHECK( zhigh.find("too high") != string::npos );
  BOOST_CHECK( zhigh.find("near 26") != string::npos );

  tr.conclusion = TrendConclusion::OtherModelInconsistency;
  const string other = ShieldSourcePullTrend::conclusion_report_text(tr);
  BOOST_CHECK( other.find("does not explain") != string::npos );
  BOOST_CHECK( other.find("despite fitting shielding") == string::npos );

  BOOST_CHECK( ShieldSourcePullTrend::areal_density_txt( 24.0*g_cm2 ) == "24 g/cm²" );
  BOOST_CHECK( ShieldSourcePullTrend::areal_density_txt( 2.46*g_cm2 ) == "2.5 g/cm²" );
}


// The amount hint must come from the amount alternative (the model's own material, scaled) - not
//  from the free-Z alternative - so it is on the right side of the model's areal density even when
//  the model material is also wrong.
BOOST_AUTO_TEST_CASE( AmountHintFromAmountAlternative )
{
  if( !ensure_static_data() )
    return;

  const TrendResult low = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{13.0, 5.0}, false ) );
  BOOST_CHECK_EQUAL( conc(low), conc(TrendConclusion::ShieldingAmountTooLow) );
  BOOST_CHECK_MESSAGE( low.suggestedArealDensity/g_cm2 > 5.0, "hint " << low.suggestedArealDensity/g_cm2 );

  const TrendResult high = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{82.0, 40.0}, false ) );
  BOOST_CHECK_EQUAL( conc(high), conc(TrendConclusion::ShieldingAmountTooHigh) );
  BOOST_CHECK_MESSAGE( high.suggestedArealDensity/g_cm2 < 40.0, "hint " << high.suggestedArealDensity/g_cm2 );
}


// A thin model shield must be able to grow to the true amount (more than 25x).
BOOST_AUTO_TEST_CASE( ThinModelShieldReachesTrueAmount )
{
  if( !ensure_static_data() )
    return;
  const TrendResult r = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{26.0, 0.5}, false ) );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::ShieldingAmountTooLow) );
  BOOST_CHECK_MESSAGE( std::fabs( r.suggestedArealDensity/g_cm2 - 20.0 ) < 4.0,
                       "hint " << r.suggestedArealDensity/g_cm2 );
}


// Material shields: the areal density is the thickness along the source-detector line, which depends
//  on the geometry (cylinder end-on: the length; rectangular: the depth) - not always dimension 0.
BOOST_AUTO_TEST_CASE( MaterialShieldGeometry )
{
  if( !ensure_static_data() )
    return;
  MaterialDB::initialize();
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> iron = matdb->material( "Fe" );
  BOOST_REQUIRE( iron );

  const double thick = 10.0*g_cm2 / iron->density;   //10 g/cm2 along the ray
  for( const GammaInteractionCalc::GeometryType geom : { GammaInteractionCalc::GeometryType::Spherical,
         GammaInteractionCalc::GeometryType::CylinderEndOn, GammaInteractionCalc::GeometryType::CylinderSideOn,
         GammaInteractionCalc::GeometryType::Rectangular } )
  {
    Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{26.0, 10.0}, false );
    ShieldingSourceFitCalc::ShieldingInfo &s = c.shields.at(0);
    s.m_isGenericMaterial = false;
    s.m_material = iron;
    s.m_geometry = geom;
    const double big = 10.0*PhysicalUnits::cm;
    switch( geom )
    {
      case GammaInteractionCalc::GeometryType::Spherical:      s.m_dimensions[0] = thick; s.m_dimensions[1] = s.m_dimensions[2] = 0.0; break;
      case GammaInteractionCalc::GeometryType::CylinderEndOn:  s.m_dimensions[0] = big;   s.m_dimensions[1] = thick; s.m_dimensions[2] = 0.0; break;
      case GammaInteractionCalc::GeometryType::CylinderSideOn: s.m_dimensions[0] = thick; s.m_dimensions[1] = big;   s.m_dimensions[2] = 0.0; break;
      case GammaInteractionCalc::GeometryType::Rectangular:    s.m_dimensions[0] = big;   s.m_dimensions[1] = big;   s.m_dimensions[2] = thick; break;
      case GammaInteractionCalc::GeometryType::NumGeometryType: break;
    }
    c.config.genericAnFit = c.config.genericAdFit = false;

    const TrendResult r = run( c );
    BOOST_CHECK_MESSAGE( r.conclusion == TrendConclusion::ShieldingAmountTooLow,
                         "geometry " << static_cast<int>(geom) << ": conclusion " << conc(r) );
    BOOST_CHECK_MESSAGE( std::fabs( r.suggestedArealDensity/g_cm2 - 20.0 ) < 4.0,
                         "geometry " << static_cast<int>(geom) << ": hint " << r.suggestedArealDensity/g_cm2 );
  }
}


// A thickness fit that went to zero is "at bound": the amount is not optimal, so an amount test still
//  applies - silent for a bare source, "lacks shielding" when the source is in fact shielded.
BOOST_AUTO_TEST_CASE( AmountFitAtZeroBound )
{
  if( !ensure_static_data() )
    return;

  const Case bare = make_case( {{ "Ba133", ba133_lines() }}, Shield{}, Shield{26.0, 1.0}, true );
  BOOST_CHECK( bare.model.ad < 0.01 );
  BOOST_CHECK_EQUAL( conc( run(bare) ), conc(TrendConclusion::NoConclusion) );

  Case stuck = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{26.0, 0.0}, false );
  stuck.shields.at(0).m_fitDimensions[1] = true;
  stuck.config.genericAdFit = true;
  const TrendResult r = run( stuck );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::ModelLacksShielding) );
  BOOST_CHECK( r.amountT > 3.5 );
}


// An age fit: no physical (amount / atomic-number) call, and the Z power is not reported - so a tool
//  reading anDiscriminable can never call a Z it was not allowed to judge "consistent".
BOOST_AUTO_TEST_CASE( AgeFitGatesPhysicalCalls )
{
  if( !ensure_static_data() )
    return;
  Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{13.0, 20.0}, true );
  c.config.anyAgeFit = true;
  const TrendResult r = run( c );
  BOOST_CHECK( (r.conclusion == TrendConclusion::NoConclusion)
               || (r.conclusion == TrendConclusion::OtherModelInconsistency) );
  BOOST_CHECK( !r.anDiscriminable );
  BOOST_CHECK( r.testsRun );
}


// A modest atomic-number error on data with the power to resolve Z: somewhere between 2.5 and 3.5
//  sigma it is reported as a soft ("possibly") call.
BOOST_AUTO_TEST_CASE( TentativeAtomicNumberCall )
{
  if( !ensure_static_data() )
    return;

  // A small Z error (model 22 vs iron's 26): smaller than the 1.5x the power is computed for, so as
  //  the uncertainties are scaled up the observed signal reaches 2.5-3.5 sigma while the power is
  //  still above its gate.
  bool found = false;
  for( double sigma_factor = 1.0; (sigma_factor < 200.0) && !found; sigma_factor *= 1.05 )
  {
    const Case c = make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{22.0, 20.0}, true, false,
                              nullptr, sigma_factor );
    const TrendResult r = run( c );
    if( !r.anDiscriminable || (std::fabs(r.atomicNumberT) < 2.6) || (std::fabs(r.atomicNumberT) > 3.4) )
      continue;
    found = true;
    BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::AtomicNumberTooLow) );
    BOOST_CHECK_EQUAL( static_cast<int>(r.confidence), static_cast<int>(TrendConfidence::Tentative) );
    BOOST_CHECK( ShieldSourcePullTrend::conclusion_report_text( r ).find( "possibly" ) != string::npos );
  }
  BOOST_CHECK_MESSAGE( found, "No uncertainty scaling gave a powered 2.6-3.4 sigma atomic-number signal" );
}


// "Lacks shielding" names a Z class only when the reference layers are distinguishable; with
//  low-energy peaks they are.
BOOST_AUTO_TEST_CASE( LacksShieldingZClass )
{
  if( !ensure_static_data() )
    return;

  const TrendResult fe = run( make_case( {{ "Eu152", eu152_lines() }}, iron20, Shield{}, false ) );
  BOOST_CHECK_EQUAL( conc(fe), conc(TrendConclusion::ModelLacksShielding) );
  BOOST_CHECK( fe.suggestedAtomicNumberResolved );
  BOOST_CHECK_CLOSE( fe.suggestedAtomicNumber, 26.0, 1.0E-6 );

  const TrendResult pb = run( make_case( {{ "Ba133", ba133_lines() }}, Shield{82.0, 3.0}, Shield{}, false ) );
  BOOST_CHECK_EQUAL( conc(pb), conc(TrendConclusion::ModelLacksShielding) );
  BOOST_CHECK( pb.suggestedAtomicNumberResolved );
  BOOST_CHECK( pb.suggestedAtomicNumber >= 50.0 );
  const ShieldSourcePullTrend::MessageHint hint = ShieldSourcePullTrend::message_hint( pb );
  BOOST_CHECK( hint.zclass_key != nullptr );
}


// The "other" conclusion does fire for a strong smooth residual that no shielding change explains
//  (here a detector-efficiency shape error, with the shielding's AN and AD both fit).
BOOST_AUTO_TEST_CASE( OtherInconsistencyFires )
{
  if( !ensure_static_data() )
    return;
  const auto eff_error = []( const double energy ){
    const double x = std::log( energy/400.0 );
    return 1.0 + 0.35*x*x;
  };
  const TrendResult r = run( make_case( {{ "Eu152", eu152_lines( 4.0 ) }}, iron20, iron20, true, true, eff_error ) );
  BOOST_CHECK_EQUAL( conc(r), conc(TrendConclusion::OtherModelInconsistency) );
  BOOST_CHECK( r.otherT >= 5.0 );
}

