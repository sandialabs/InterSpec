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

// Validates the Reference-Photopeak tool's cascade-summing integration.
//
// A displayed line is scaled by the summing-OUT survival the CeeLo engine reports
// for its energy window, and then grown by any coincidence sum that lands on top
// of it; sums that fall clear of every line are drawn as their own peaks.  The
// cases here check each half against an independently written reconstruction on
// the same cascades and detector, so a regression in the window building, the
// memoized efficiency provider, the energy-to-line mapping, or the sum
// enumeration shows up.  A couple of physical values are pinned (Co-60 1173/1332
// summing-out at 1 cm) so a detector-fixture or engine change that breaks the
// magnitude is caught.
//
// CascadeTimingReport is a report, not a check: pass `--timing` after the `--`
// to print the CPU time of a display update.  It is worth re-running in a Release
// build after any change to the summing path - the tool recomputes on every
// option change, so this has to stay well under a second.  Measured on an M-series
// Mac at 2.54 cm, first use / repeat: Co-60 0.014 / 0.000 s, Ba-133 0.023 / 0.000,
// Eu-152 0.10 / 0.002, Ra-226 0.18 / 0.013, U-238 chain 0.22 / 0.020.  The first
// use of a detector and distance fills the efficiency interpolation table
// (CascadeSummingCalc::interpolatedEff); everything after that reuses it.

#include <cmath>
#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#define BOOST_TEST_MODULE CascadeRefLineSumming_suite
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

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

  //A near-contact distance, where summing is large enough to pin values against.
  //  This test builds its own efficiency functors, so it need not match the GUI default.
  const double kDist = 1.0 * PhysicalUnits::cm;

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
    auto tot = [&det]( const double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, kDist ); };

    map<long long,double> out;
    for( const ceelo::DecayCascade &dc : cascades )
    {
      const size_t n = dc.members.size();

      //Survival of the summed pair against the branch's other coincident photons.
      vector<double> msurv( n, 1.0 );
      double bsurv = 1.0;
      for( size_t m = 0; m < n; ++m )
      {
        const double e = dc.members[m].energy_keV;
        if( (e < 5.0) || (dc.members[m].type == ceelo::CascadeParticleType::Annih511) )
          continue;
        const double dep = std::min( 1.0, dc.members[m].intensity ) * tot(e);
        if( dep <= 1.0E-6 )
          continue;
        msurv[m] = std::max( 1.0E-6, 1.0 - dep );
        bsurv *= msurv[m];
      }

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
          gain *= bsurv / (msurv[ai] * msurv[bi]);
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


// The separate (pure) sum peaks now come from the engine as well: each accepted
// candidate energy is handed to compute_cascade_analytic as a window holding no
// emitted line, and the returned absolute rate is the peak amplitude.  Check the
// peaks the user actually looks for are present, labelled, and match a direct
// engine call on the same cascades and detector.
BOOST_AUTO_TEST_CASE( SumPeaksMatchEngine )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const shared_ptr<DetectorPeakResponse> det = example_detector();

  //Sums that are NOT on top of an emitted line, so they are drawn as their own
  //  peaks.  (Co-60's famous 2505.7 keV sum is deliberately not here: Co-60 really
  //  does emit a very weak 2505.7 keV crossover gamma, so the sum folds into that
  //  line instead - checked separately below.)
  struct Expect { const char *nuc; double energy; const char *what; };
  const Expect expected[] = {
    { "Ba133",  357.4, "276.4 + 81.0" },          //resolvable, just above the 356 line
    { "Co60",  1679.6, "347.1 + 1332.5" },        //a weaker branch, no line at the sum
  };

  size_t n_triple_labels = 0;
  for( const Expect &ex : expected )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( ex.nuc );
    BOOST_REQUIRE( nuc );

    RefLineInput input;
    input.m_input_txt = ex.nuc;
    input.m_age = "20 y";
    input.m_showGammas = true;
    input.m_showXrays = true;
    input.m_showCascades = true;
    input.m_lower_br_cutt_off = 0.0;
    input.m_do_cascade_summing = true;
    input.m_summing_fep_eff = [det]( double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, kDist ); };
    input.m_summing_total_eff = [det]( double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, kDist ); };

    const shared_ptr<ReferenceLineInfo> info = ReferenceLineInfo::generateRefLineInfo( input );
    BOOST_REQUIRE( info && (info->m_validity == ReferenceLineInfo::InputValidity::Valid) );

    const ReferenceLineInfo::RefLine *found = nullptr;
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
    {
      if( l.m_source_type != ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak )
        continue;
      if( std::fabs(l.m_energy - ex.energy) < 0.5 )
        found = &l;
      //Count triple-labelled sums seen anywhere in the set.
      if( std::count( begin(l.m_decaystr), end(l.m_decaystr), '+' ) == 2 )
        ++n_triple_labels;
    }

    BOOST_REQUIRE_MESSAGE( found, ex.nuc << ": no cascade-sum peak near " << ex.energy
                          << " keV (" << ex.what << ")" );
    BOOST_CHECK_GT( found->m_decay_intensity, 0.0 );
    BOOST_CHECK_MESSAGE( SpecUtils::icontains( found->m_decaystr, "Cascade sum" ),
                        "Unexpected label: '" << found->m_decaystr << "'" );

    // Recompute the amplitude here from the same cascades and detector: the joint
    //  full-energy rate of the labelled pair, times the survival of the pair
    //  against every other coincident photon of that decay.  This is the same
    //  formula production uses, so it pins the WIRING - that the labelled pair is
    //  really the one that made the peak, that the peak carries that pair's rate,
    //  and that the window/energy bookkeeping in between is right - rather than
    //  independently validating the physics.  The physics itself is checked
    //  against the CeeLo engine in PairwiseVsExactEngine (summing-out) and in
    //  CeeLo's own suite.
    ceelo::cascade_adapter::CascadeOptions copt;
    copt.age_seconds = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife(
                          info->m_input.m_age, nuc->halfLife ) / PhysicalUnits::second;
    copt.prompt_equilibrium = false;
    copt.include_xrays = true;
    copt.vacancy_xray_model = false;   //match production
    const vector<ceelo::DecayCascade> cascades = ceelo::cascade_adapter::build_cascades( nuc, copt );

    auto fep = [&det]( double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, kDist ); };
    auto tot = [&det]( double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, kDist ); };

    double recon = 0.0;
    for( const ceelo::DecayCascade &dc : cascades )
    {
      const size_t nmem = dc.members.size();
      vector<double> msurv( nmem, 1.0 );
      double bsurv = 1.0;
      for( size_t m = 0; m < nmem; ++m )
      {
        const double e = dc.members[m].energy_keV;
        if( (e < 5.0) || (dc.members[m].type == ceelo::CascadeParticleType::Annih511) )
          continue;
        const double dep = std::min( 1.0, dc.members[m].intensity ) * tot(e);
        if( dep <= 1.0E-6 )
          continue;
        msurv[m] = std::max( 1.0E-6, 1.0 - dep );
        bsurv *= msurv[m];
      }

      std::set<std::pair<uint16_t,uint16_t>> seen;
      for( size_t a = 0; a < nmem; ++a )
      {
        const double ea = dc.members[a].energy_keV;
        if( ea < 5.0 ) continue;
        for( const ceelo::CoincidenceLink &lk : dc.members[a].coincident )
        {
          if( lk.partner >= nmem ) continue;
          const uint16_t ai = static_cast<uint16_t>(a), bi = lk.partner;
          if( !seen.insert( make_pair( std::min(ai,bi), std::max(ai,bi) ) ).second ) continue;
          const double eb = dc.members[bi].energy_keV;
          if( eb < 5.0 ) continue;
          if( std::fabs( (ea + eb) - found->m_energy ) > 0.25 ) continue;

          double g = dc.branch_weight * dc.members[a].intensity * lk.prob * fep(ea) * fep(eb);
          if( lk.has_correlation )
            g *= std::max( 0.0, 1.0 + lk.a2 + lk.a4 );
          g *= bsurv / (msurv[ai] * msurv[bi]);
          recon += g;
        }
      }
    }//for( cascades )

    BOOST_CHECK_MESSAGE( recon > 0.0, ex.nuc << ": reconstruction found no pair at "
                        << found->m_energy << " keV" );
    BOOST_CHECK_CLOSE( found->m_decay_intensity, recon, 0.5 );

  }//for( expected )

  // Eu-152 has enough coincident gammas that some three-photon sums survive the
  //  ranking; that is the path the pairwise-only predecessor could not produce.
  {
    RefLineInput input;
    input.m_input_txt = "Eu152";
    input.m_age = "20 y";
    input.m_showGammas = true;
    input.m_showXrays = true;
    input.m_showCascades = true;
    input.m_lower_br_cutt_off = 0.0;
    input.m_do_cascade_summing = true;
    input.m_summing_fep_eff = [det]( double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, kDist ); };
    input.m_summing_total_eff = [det]( double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, kDist ); };

    const shared_ptr<ReferenceLineInfo> info = ReferenceLineInfo::generateRefLineInfo( input );
    BOOST_REQUIRE( info );
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
      if( (l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak)
         && (std::count( begin(l.m_decaystr), end(l.m_decaystr), '+' ) == 2) )
        ++n_triple_labels;
  }
  BOOST_CHECK_MESSAGE( n_triple_labels > 0,
      "no three-photon cascade sum was produced for any test nuclide" );

  // Co-60's 2505.7 keV: the nuclide emits a real but almost immeasurably weak
  //  crossover gamma there, and at contact essentially everything seen at that
  //  energy is the 1173+1332 coincidence sum.  So it must stay a Normal line and
  //  be lifted enormously above its own branching ratio, rather than being
  //  duplicated as a separate sum peak.
  {
    RefLineInput input;
    input.m_input_txt = "Co60";
    input.m_age = "20 y";
    input.m_showGammas = true;
    input.m_showXrays = true;
    input.m_showCascades = true;
    input.m_lower_br_cutt_off = 0.0;
    input.m_do_cascade_summing = true;
    input.m_summing_fep_eff = [det]( double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, kDist ); };
    input.m_summing_total_eff = [det]( double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, kDist ); };

    const shared_ptr<ReferenceLineInfo> info = ReferenceLineInfo::generateRefLineInfo( input );
    BOOST_REQUIRE( info );

    const ReferenceLineInfo::RefLine *crossover = nullptr;
    size_t n_at_energy = 0;
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
    {
      if( std::fabs(l.m_energy - 2505.7) > 0.5 )
        continue;
      ++n_at_energy;
      if( l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal )
        crossover = &l;
    }

    BOOST_REQUIRE_MESSAGE( crossover, "Co60: no 2505.7 keV line" );
    BOOST_CHECK_MESSAGE( n_at_energy == 1,
        "Co60: 2505.7 keV should be one line, but there are " << n_at_energy );
    BOOST_CHECK_MESSAGE( crossover->m_normalized_intensity > 10.0*crossover->m_uncorrected_intensity,
        "Co60 2505.7 keV should be dominated by summing-in, but grew only from "
        << crossover->m_uncorrected_intensity << " to " << crossover->m_normalized_intensity );
  }
}//BOOST_AUTO_TEST_CASE( SumPeaksMatchEngine )


// The assumed summing distance defaults to a near-contact 2.54 cm, and a stored
// state without one (or with an unusable one) comes back to that same default.
BOOST_AUTO_TEST_CASE( DefaultCascadeDistance )
{
  set_data_dir();

  const RefLineInput fresh;
  BOOST_CHECK_EQUAL( fresh.m_cascade_distance, string(RefLineInput::sm_default_cascade_distance) );
  double dist = 0.0;
  BOOST_REQUIRE_NO_THROW( dist = PhysicalUnits::stringToDistance( fresh.m_cascade_distance ) );
  BOOST_CHECK_CLOSE( dist, 2.54*PhysicalUnits::cm, 1.0E-6 );

  //A stored state that predates the field, and one holding a junk value, both come
  //  back as the default rather than as a zero distance.
  {
    RefLineInput input;
    input.m_input_txt = "Co60";
    input.m_cascade_distance = "not a distance";

    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( input.serialize( &doc ) );

    RefLineInput restored;
    BOOST_REQUIRE_NO_THROW( restored.deSerialize( doc.first_node() ) );
    BOOST_CHECK_EQUAL( restored.m_cascade_distance,
                       string(RefLineInput::sm_default_cascade_distance) );
  }

  {
    RefLineInput input;
    input.m_input_txt = "Co60";
    input.m_cascade_distance = "";   //as a state written before the field existed

    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( input.serialize( &doc ) );

    RefLineInput restored;
    BOOST_REQUIRE_NO_THROW( restored.deSerialize( doc.first_node() ) );
    BOOST_CHECK_EQUAL( restored.m_cascade_distance,
                       string(RefLineInput::sm_default_cascade_distance) );
  }

  //A good value must of course survive.
  {
    RefLineInput input;
    input.m_input_txt = "Co60";
    input.m_cascade_distance = "10 cm";

    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( input.serialize( &doc ) );

    RefLineInput restored;
    BOOST_REQUIRE_NO_THROW( restored.deSerialize( doc.first_node() ) );
    BOOST_CHECK_EQUAL( restored.m_cascade_distance, string("10 cm") );
  }
}//BOOST_AUTO_TEST_CASE( DefaultCascadeDistance )


// Timing report for the interactive path: the Reference Photopeak tool re-runs
//  generateRefLineInfo (with cascade summing) on the session thread for every
//  option change, so it must stay well under a second in a Release build.  Only
//  runs when `--timing` is passed after the `--` on the command line; prints CPU
//  time (std::clock, min/median of 5 repeats after a warm-up) per nuclide.
BOOST_AUTO_TEST_CASE( CascadeTimingReport )
{
  bool want_timing = false;
  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
    want_timing |= (string(argv[i]) == "--timing");
  if( !want_timing )
    return;

  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const shared_ptr<DetectorPeakResponse> det = example_detector();
  const double dist = 2.54 * PhysicalUnits::cm;  //the GUI default

  struct TimeCase { const char *nuc; const char *age; };
  const TimeCase cases[] = {
    {"Co60","20 y"}, {"Cs137","20 y"}, {"Ba133","20 y"}, {"Am241","20 y"},
    {"Eu152","20 y"}, {"Ra226","20 y"}, {"Th232","20 y"}, {"U238","20 y"}, {"U238","1 y"}
  };

  cout << "\nCascade ref-line timing (CPU seconds, example HPGe at 2.54 cm):\n"
       << "  nuclide  age    no-casc  1st-use  repeat  real-lines  sum-peaks\n";
  for( const TimeCase &tc : cases )
  {
    BOOST_REQUIRE_MESSAGE( db->nuclide( tc.nuc ), "no nuclide " << tc.nuc );

    RefLineInput input;
    input.m_input_txt = tc.nuc;
    input.m_age = tc.age;
    input.m_showGammas = true;
    input.m_showXrays = true;
    input.m_showCascades = true;
    input.m_lower_br_cutt_off = 0.0;
    input.m_do_cascade_summing = true;
    //Through the same interpolation the GUI uses, so the timing is the real one.
    //  A key unique to this case makes the first call below pay the full
    //  efficiency-table fill, which is what a user sees on their first toggle.
    const string key = string(tc.nuc) + tc.age;
    input.m_summing_fep_eff = CascadeSummingCalc::interpolatedEff( key + "-fep",
        [det,dist]( double e ){ return CascadeSummingCalc::detectorFepEffAbs( det, e, 0, 0, dist ); } );
    input.m_summing_total_eff = CascadeSummingCalc::interpolatedEff( key + "-tot",
        [det,dist]( double e ){ return CascadeSummingCalc::detectorTotEffAbs( det, e, 0, 0, dist ); } );

    const std::clock_t cold_start = std::clock();
    shared_ptr<ReferenceLineInfo> info = ReferenceLineInfo::generateRefLineInfo( input );
    const double cold = double(std::clock() - cold_start) / CLOCKS_PER_SEC;
    BOOST_REQUIRE( info && (info->m_validity == ReferenceLineInfo::InputValidity::Valid) );

    //The same lines with cascade summing off, to separate the summing cost from
    //  the (unavoidable) decay/line enumeration cost.
    RefLineInput plain = input;
    plain.m_showCascades = false;
    plain.m_do_cascade_summing = false;
    plain.m_summing_fep_eff = nullptr;
    plain.m_summing_total_eff = nullptr;
    ReferenceLineInfo::generateRefLineInfo( plain );  //warm-up

    vector<double> secs, plain_secs;
    for( int rep = 0; rep < 5; ++rep )
    {
      std::clock_t start = std::clock();
      ReferenceLineInfo::generateRefLineInfo( plain );
      plain_secs.push_back( double(std::clock() - start) / CLOCKS_PER_SEC );

      start = std::clock();
      info = ReferenceLineInfo::generateRefLineInfo( input );
      secs.push_back( double(std::clock() - start) / CLOCKS_PER_SEC );
    }
    std::sort( begin(secs), end(secs) );
    std::sort( begin(plain_secs), end(plain_secs) );

    size_t nreal = 0, nsum = 0;
    for( const ReferenceLineInfo::RefLine &l : info->m_ref_lines )
    {
      nreal += (l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal);
      nsum += (l.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak);
    }

    char line[256];
    const double casc = secs[secs.size()/2], plainmed = plain_secs[plain_secs.size()/2];
    snprintf( line, sizeof(line), "  %-8s %-5s  %7.3f  %7.3f  %6.3f  %10zu  %9zu\n",
              tc.nuc, tc.age, plainmed, cold, casc, nreal, nsum );
    cout << line;
  }//for( cases )
  cout << endl;
}//BOOST_AUTO_TEST_CASE( CascadeTimingReport )
