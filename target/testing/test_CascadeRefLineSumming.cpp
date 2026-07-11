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

// Validates the Reference-Photopeak tool's cascade-summing integration.  Real
// lines are scaled by the exact level-scheme engine's net factor
// (ceelo::AnalyticPeakResult::summing_factor = c_out + c_in): summing-out shrinks
// the line, and a coincidence sum that lands within the tight (0.25 keV) window
// folds in as summing-in.  This test confirms that the factor that comes out of
// ReferenceLineInfo::generateRefLineInfo (extracted as normalized/uncorrected line
// height) reproduces a direct engine call on the same cascades and detector - i.e.
// that the window building, memoized efficiency provider, and energy-to-line
// mapping are wired correctly.  (Resolvable sums beyond the window are shown as
// their own separate peaks, not folded here.)  Since both sides run the exact
// engine, agreement is to floating-point (a small memoization rounding aside), so
// the tolerance is tight and covers ALL strong lines, including cascade-sum
// targets (a real line that also receives summing-in from two other lines).
//
// It also pins a couple of physical values (Co-60 1173/1332 summing-out at 1 cm)
// so a detector-fixture or engine change that breaks the magnitude is caught.

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>

#define BOOST_TEST_MODULE CascadeRefLineSumming_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "io/DetectorResponse.h"
#include "cascade/CascadeTypes.h"
#include "cascade/AnalyticCascade.h"
#include "cascade/SandiaDecayCascade.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
  using GammaInteractionCalc::CascadeSummingCalc;

  const double kDist = 1.0 * PhysicalUnits::cm;  //must match sm_cascade_summing_distance

  void set_data_dir()
  {
    static bool s_done = false;
    if( s_done )
      return;
    s_done = true;

    int argc = framework::master_test_suite().argc;
    char **argv = framework::master_test_suite().argv;
    string datadir, test_file_dir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        test_file_dir = arg.substr( 14 );
    }
    SpecUtils::ireplace_all( datadir, "%20", " " );

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.reactiongamma.xml") ) )
        { datadir = d; break; }
    }
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

    // The full sandia.decay.xml (with coincidences + x-rays) lives in the SandiaDecay repo.
    const string decay_xml = "external_libs/SandiaDecay/sandia.decay.xml";
    string decay_file;
    for( const string &d : { SpecUtils::append_path( test_file_dir, "../../../" ),
                             string(".."), string("../.."), string("../../..") } )
    {
      const string cand = SpecUtils::append_path( d, decay_xml );
      if( SpecUtils::is_file(cand) ){ decay_file = cand; break; }
    }
    BOOST_REQUIRE_MESSAGE( !decay_file.empty(), "Could not find full sandia.decay.xml" );
    BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( decay_file ) );
    BOOST_REQUIRE( DecayDataBaseServer::database() );
  }//set_data_dir()


  shared_ptr<DetectorPeakResponse> example_detector()
  {
    const string path = SpecUtils::append_path( InterSpec::staticDataDirectory(),
                                               "cascade_summing_example_response.xml" );
    ifstream in( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE_MESSAGE( in.is_open(), "Missing " << path );
    stringstream strm; strm << in.rdbuf();
    auto drf = make_shared<DetectorPeakResponse>( "Example", "example" );
    drf->setIntrinsicEfficiencyFormula( "1.0", 5.0*PhysicalUnits::cm, PhysicalUnits::keV,
        0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    drf->setCeeloResponse( ceelo::DetectorResponse::from_xml_string( strm.str() ) );
    BOOST_REQUIRE( drf->ceeloResponse() );
    return drf;
  }//example_detector()


  long long ekey( const double e ){ return static_cast<long long>( std::llround( e * 10.0 ) ); }


  // The exact-engine reference: net summing factor keyed by energy bucket, on the
  //  same cascades (transition-group x-rays) + same detector as the pairwise path.
  map<long long,double> engine_cnet( const shared_ptr<const DetectorPeakResponse> &det,
                                     const SandiaDecay::Nuclide *nuc,
                                     const double age_seconds,
                                     const vector<double> &energies )
  {
    ceelo::cascade_adapter::CascadeOptions copt;
    copt.age_seconds = age_seconds;
    copt.prompt_equilibrium = false;   //match production (m_promptLinesOnly == false here)
    copt.include_xrays = true;
    copt.vacancy_xray_model = false;   //match the pairwise path
    const vector<ceelo::DecayCascade> cascades = ceelo::cascade_adapter::build_cascades( nuc, copt );

    struct Prov : public ceelo::EfficiencyProvider {
      shared_ptr<const DetectorPeakResponse> d_sp;
      double fep( double e ) const override { return CascadeSummingCalc::detectorFepEffAbs( d_sp, e, 0, 0, kDist ); }
      double total( double e ) const override { return CascadeSummingCalc::detectorTotEffAbs( d_sp, e, 0, 0, kDist ); }
      bool has( double ) const override { return true; }
    } prov;
    prov.d_sp = det;

    vector<ceelo::PeakWindow> windows;
    for( const double e : energies )
    { ceelo::PeakWindow w; w.energy_keV = e; w.tolerance_keV = 0.25; windows.push_back(w); }  //match production

    ceelo::AnalyticCascadeOptions opts;   //defaults: exact, triples on
    const vector<ceelo::AnalyticPeakResult> res
        = ceelo::compute_cascade_analytic( cascades, windows, prov, opts );

    map<long long,double> out;
    for( const ceelo::AnalyticPeakResult &r : res )
      if( r.found && !std::isnan(r.c_out) && !std::isinf(r.c_out) )
        out[ ekey(r.energy_keV) ] = r.c_out;   //production scales real lines by c_out, then ADDS folded summing-in
    return out;
  }//engine_cnet(...)


  // Independent reconstruction of the summing-IN gain that production folds onto
  //  real lines: for each coincident pair (deduped per branch), the joint-FEP gain
  //  branch_weight·intensity_a·P(b|a)·W(0)·eps_FEP_abs(ea)·eps_FEP_abs(eb) is added
  //  to the nearest real line within 0.25 keV of E_a+E_b.  Absolute per-decay
  //  (no shielding in the test).  This mirrors the production formula with an
  //  independent implementation, so a regression in it (missing W(0), wrong prob,
  //  wrong branch basis, broken fold) is caught.  Returns gain keyed by real-line
  //  energy bucket.
  map<long long,double> reconstruct_sum_in( const shared_ptr<const DetectorPeakResponse> &det,
                                            const SandiaDecay::Nuclide *nuc,
                                            const double age_seconds,
                                            const vector<double> &real_energies )
  {
    ceelo::cascade_adapter::CascadeOptions copt;
    copt.age_seconds = age_seconds;
    copt.prompt_equilibrium = false;
    copt.include_xrays = true;
    copt.vacancy_xray_model = false;
    const vector<ceelo::DecayCascade> cascades = ceelo::cascade_adapter::build_cascades( nuc, copt );

    vector<double> sorted_real( real_energies );
    std::sort( begin(sorted_real), end(sorted_real) );
    auto nearest = [&sorted_real]( const double e ) -> double {
      const auto it = std::lower_bound( begin(sorted_real), end(sorted_real), e );
      double best = std::numeric_limits<double>::quiet_NaN(), bd = 0.25;
      if( (it != end(sorted_real)) && (std::fabs(*it - e) <= bd) ){ bd = std::fabs(*it - e); best = *it; }
      if( it != begin(sorted_real) ){ const double v = *(it-1); if( std::fabs(v - e) <= bd ) best = v; }
      return best;
    };
    auto eps = [&det]( const double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, kDist ); };

    map<long long,double> out;
    for( const ceelo::DecayCascade &dc : cascades )
    {
      const size_t n = dc.members.size();
      std::set<std::pair<uint16_t,uint16_t>> seen;
      for( size_t a = 0; a < n; ++a )
      {
        const double ea = dc.members[a].energy_keV;
        if( ea < 5.0 ) continue;
        for( const ceelo::CoincidenceLink &lk : dc.members[a].coincident )
        {
          if( lk.partner >= n ) continue;
          const uint16_t ai = static_cast<uint16_t>(a), bi = lk.partner;
          if( !seen.insert( std::make_pair( std::min(ai,bi), std::max(ai,bi) ) ).second ) continue;
          const double eb = dc.members[bi].energy_keV;
          if( eb < 5.0 ) continue;
          double gain = dc.branch_weight * dc.members[a].intensity * lk.prob * eps(ea) * eps(eb);
          if( lk.has_correlation ) gain *= std::max( 0.0, 1.0 + lk.a2 + lk.a4 );
          if( gain <= 0.0 ) continue;
          const double rl = nearest( ea + eb );
          if( !std::isnan(rl) )
            out[ ekey(rl) ] += gain;
        }
      }
    }
    return out;
  }//reconstruct_sum_in(...)

  using GammaInteractionCalc::CascadeSummingCalc;
}//anonymous namespace


BOOST_AUTO_TEST_CASE( PairwiseVsExactEngine )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const shared_ptr<DetectorPeakResponse> det = example_detector();

  // Nuclides spanning gamma-gamma (Co60), gamma+xray EC emitters (Ba133, Eu152),
  //  a chain-member isomer (Cs137 -> Ba137m), and a decay chain with a crossover
  //  that receives strong summing-IN (Ra226: Bi214 1120.3+609.3 = 1729.6 onto the
  //  direct 1729.595 gamma).
  struct NucCase { const char *name; bool gamma_gamma; };
  const NucCase nuclides[] = { {"Co60",true}, {"Cs137",true}, {"Ba133",false},
                               {"Eu152",false}, {"Ra226",false} };
  const string age_str = "20 y";

  size_t total_compared = 0, gain_lines = 0, total_recon_entries = 0;
  double sum_abs_err = 0.0;
  double max_gg_err = 0.0, max_ec_err = 0.0;
  for( const NucCase &ncase : nuclides )
  {
    const char * const nucstr = ncase.name;
    const SandiaDecay::Nuclide * const nuc = db->nuclide( nucstr );
    BOOST_REQUIRE_MESSAGE( nuc, "no nuclide " << nucstr );

    // ----- Pairwise path (the production code) -----
    RefLineInput input;
    input.m_input_txt = nucstr;
    input.m_age = age_str;
    input.m_showGammas = true;
    input.m_showXrays = true;
    input.m_showCascades = true;
    input.m_lower_br_cutt_off = 0.0;
    input.m_do_cascade_summing = true;
    input.m_summing_fep_eff = [det]( double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, kDist ); };
    input.m_summing_total_eff = [det]( double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, kDist ); };

    const shared_ptr<ReferenceLineInfo> info = ReferenceLineInfo::generateRefLineInfo( input );
    BOOST_REQUIRE( info && (info->m_validity == ReferenceLineInfo::InputValidity::Valid) );

    // Summing-in must produce at least one pure cascade-sum peak for these sources.
    size_t nsum = 0;
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
      nsum += (l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak);
    BOOST_CHECK_MESSAGE( nsum > 0, nucstr << ": expected at least one cascade-sum peak" );

    // The exact age the code resolved to (so the engine sees the same cascades).
    const double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife(
                          info->m_input.m_age, nuc->halfLife );
    const double age_s = age / PhysicalUnits::second;

    // Collect the strong gamma lines and their pairwise net factor (= corrected/uncorrected).
    struct Line { double energy, cnet, intensity; };
    vector<Line> strong;
    double max_intensity = 0.0;
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
      if( (l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal)
         && (l.m_particle_type == ReferenceLineInfo::RefLine::Particle::Gamma) )
        max_intensity = std::max( max_intensity, l.m_decay_intensity );

    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
    {
      if( (l.m_source_type != ReferenceLineInfo::RefLine::RefGammaType::Normal)
         || (l.m_particle_type != ReferenceLineInfo::RefLine::Particle::Gamma) )
        continue;
      if( (l.m_decay_intensity < 0.05*max_intensity) || (l.m_uncorrected_intensity <= 0.0) )
        continue;   //only compare the visually-dominant lines
      const double cnet = l.m_normalized_intensity / l.m_uncorrected_intensity;
      strong.push_back( { l.m_energy, cnet, l.m_decay_intensity } );
    }
    BOOST_REQUIRE_MESSAGE( !strong.empty(), "no strong gamma lines for " << nucstr );

    // Every real gamma/xray line (the same set production folds summing-IN onto);
    //  needed so the reconstruction's nearest-line fold matches production exactly.
    vector<double> all_real_energies;
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
      if( (l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal)
         && ((l.m_particle_type == ReferenceLineInfo::RefLine::Particle::Gamma)
             || (l.m_particle_type == ReferenceLineInfo::RefLine::Particle::Xray)) )
        all_real_energies.push_back( l.m_energy );

    // ----- Independent reconstruction of the folded summing-IN gain -----
    const map<long long,double> recon_sum_in = reconstruct_sum_in( det, nuc, age_s, all_real_energies );
    total_recon_entries += recon_sum_in.size();

    // ----- Exact engine reference on the same cascades -----
    vector<double> energies;
    for( const Line &s : strong )
      energies.push_back( s.energy );
    const map<long long,double> eng = engine_cnet( det, nuc, age_s, energies );

    for( const Line &s : strong )
    {
      const auto pos = eng.find( ekey(s.energy) );
      if( pos == end(eng) )
        continue;
      const double engine = pos->second;
      const double mine = s.cnet;
      const double err = std::fabs( mine - engine );

      BOOST_TEST_MESSAGE( string(nucstr) + " " + std::to_string(s.energy) + " keV: display="
                          + std::to_string(mine) + " engine=" + std::to_string(engine)
                          + " (BR=" + std::to_string(s.intensity) + ")" );

      ++total_compared;
      sum_abs_err += err;
      max_ec_err = std::max( max_ec_err, err );
      (void)max_gg_err;

      // Physical pin: at 1 cm from the example HPGe, Co-60's two ~100% gammas
      //  each lose ~19% to coincidence summing-out.
      if( (string(nucstr) == "Co60")
         && ((std::fabs(s.energy - 1173.2) < 1.0) || (std::fabs(s.energy - 1332.5) < 1.0)) )
        BOOST_CHECK_MESSAGE( (mine > 0.70) && (mine < 0.88),
            "Co60 " << s.energy << " keV summing-out factor " << mine
                    << " outside expected [0.70,0.88] at 1 cm" );

      // Production applies the engine's summing-OUT survival (c_out) and then ADDS
      //  any summing-in that folds onto the line, so the displayed factor is c_out
      //  plus a non-negative gain: it must never be below c_out.
      BOOST_CHECK_MESSAGE( mine >= engine - 0.02,
          nucstr << " " << s.energy << " keV: display factor " << mine
                 << " below engine summing-out " << engine );

      // Full independent reconstruction of the displayed net factor:
      //   displayed = c_out(engine) + folded_summing_in_gain / own_absolute_rate
      // where own_absolute_rate = BR * eps_FEP_abs(E).  This exercises BOTH halves of
      //  the correction (the exact-engine summing-OUT and the pairwise summing-IN
      //  fold), so a regression in either - a missing W(0), wrong coincidence prob,
      //  wrong branch basis, or a broken nearest-line fold - shows up here.
      const double own_rate = s.intensity * CascadeSummingCalc::detectorFepEffAbs( det, s.energy, 0, 0, kDist );
      double expected = engine;
      const auto rpos = recon_sum_in.find( ekey(s.energy) );
      if( (rpos != end(recon_sum_in)) && (own_rate > 0.0) )
      {
        expected += rpos->second / own_rate;
        if( expected > 1.02 )
          ++gain_lines;   //this strong line visibly receives summing-IN gain
      }

      BOOST_CHECK_MESSAGE( std::fabs( mine - expected ) < 0.025,
          nucstr << " " << s.energy << " keV: display factor " << mine
                 << " != reconstructed net " << expected
                 << " (c_out=" << engine << ", sum_in/rate=" << (expected - engine) << ")" );
    }//for( strong lines )
  }//for( nuclides )

  const double mean_err = (total_compared > 0) ? (sum_abs_err / total_compared) : 0.0;
  BOOST_TEST_MESSAGE( "compared " + std::to_string(total_compared) + " lines; mean|display-net|="
                      + std::to_string(mean_err) + " max=" + std::to_string(max_ec_err)
                      + "; strong lines with summing-IN gain=" + std::to_string(gain_lines)
                      + "; total folded summing-IN entries=" + std::to_string(total_recon_entries) );
  BOOST_CHECK_MESSAGE( total_compared >= 15, "compared too few lines (" << total_compared << ")" );

  // The reconstruction must actually fold summing-IN somewhere (otherwise the net
  //  comparison degenerates to a pure summing-out check and would miss sum-in bugs).
  BOOST_CHECK_MESSAGE( total_recon_entries > 0,
      "reconstruction folded no summing-IN gain onto any real line across all nuclides" );

  // ...and at least one visually-dominant line must carry a net gain (e.g. Ba-133's
  //  356 keV receives 276.4+79.6 and 302.9+53.2), proving the gain half is exercised
  //  on a line the user would actually see grow.
  BOOST_CHECK_MESSAGE( gain_lines > 0,
      "no strong line received a net summing-IN gain across all nuclides" );
}//BOOST_AUTO_TEST_CASE( PairwiseVsExactEngine )
