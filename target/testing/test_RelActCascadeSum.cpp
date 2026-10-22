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

/* Cascade-summing correction of "Isotopics by nuclides" (`RelEffCurveInput::correct_for_cascade_summing`).

 A noiseless Ba133 + Cs137 spectrum is built with every line's area multiplied by its true-coincidence summing factor,
 computed directly with `CascadeSummingCalc` for an absolute efficiency `s * intrinsic(E)` (no attenuation, so the
 Physical-model relative efficiency is exactly the detector's intrinsic efficiency).  A fit with the correction should
 recover `s`, the per-line factors, and the Ba133/Cs137 activity ratio; a fit without it should not.

 A second, independent check uses Monte-Carlo truth: the line efficiencies (with and without summing) of an Eu152 +
 Cs137 point source near a micro-Detective, from CeeLo's cascade MC for the real detector geometry (with and without
 an iron shield), stored in test_data/rel_act_cascade/.  The fits should recover the Eu152/Cs137 ratio.

 Needs the InterSpec "data" directory (`-- --datadir=...`) and, for the Monte-Carlo case, the test-data directory
 (`--testfiledir=...`).
 */

#include "InterSpec_config.h"

#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <utility>
#include <iostream>
#include <functional>

#define BOOST_TEST_MODULE RelActCascadeSum_suite
#include <boost/test/included/unit_test.hpp>

#include <boost/math/special_functions/erf.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;


namespace
{
  string sm_datadir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = framework::master_test_suite().argc;
    char ** const argv = framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( datadir, "sandia.decay.xml" ) ),
                           "sandia.decay.xml not in '" << datadir << "'; pass -- --datadir=..." );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
    BOOST_REQUIRE( DecayDataBaseServer::database() );
    sm_datadir = datadir;
  }//set_data_dir()


  /** The GADRAS "HPGe 40%" response: it has FWHM and total-efficiency information. */
  shared_ptr<const DetectorPeakResponse> hpge_drf()
  {
    const string dir = SpecUtils::append_path( SpecUtils::append_path( sm_datadir, "GenericGadrasDetectors" ), "HPGe 40%" );
    auto drf = make_shared<DetectorPeakResponse>();
    drf->fromGadrasDirectory( dir );
    return drf;
  }//hpge_drf()


  struct Source
  {
    const SandiaDecay::Nuclide *nuc = nullptr;
    double age = 0.0;
    double rel_act = 0.0;   //what the fit should find: line area = rel_act * live_time * intrinsic(E) * yield * c(E)
  };

  const double sm_live_time = 1000.0;
  const vector<pair<double,double>> sm_rois{ {270.0, 310.0}, {345.0, 395.0}, {650.0, 675.0} };

  double window_tolerance( const shared_ptr<const DetectorPeakResponse> &drf, const double energy )
  {
    return std::max( 0.25, 1.25 * drf->peakResolutionFWHM( static_cast<float>(energy) ) / 2.35482 );
  }


  /** A pure coincidence-sum peak put in a synthetic spectrum. */
  struct TrueSum
  {
    const SandiaDecay::Nuclide *nuc = nullptr;
    double energy = 0.0, counts = 0.0;
  };

  /** Noiseless spectrum of `sources`, with summing at absolute FEP efficiency `s * intrinsic(E)` (no summing if
   `s` is zero): every line's area times its summing factor, plus the pure coincidence-sum peaks inside the ROIs
   (grouped into windows the way the fit groups them).  `true_factors` gets each line's summing factor, keyed by
   (nuclide, energy), and `true_sums` the sum peaks. */
  shared_ptr<SpecUtils::Measurement> make_spectrum( const shared_ptr<const DetectorPeakResponse> &drf,
                                                    const vector<Source> &sources, const double s,
                                                    map<pair<const SandiaDecay::Nuclide *,double>,double> &true_factors,
                                                    vector<TrueSum> *true_sums = nullptr )
  {
    using GammaInteractionCalc::CascadeSummingCalc;

    const size_t nchannel = 8192;
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( nchannel, { 0.0f, 0.25f }, {} );

    auto counts = make_shared<vector<float>>( nchannel, 0.0f );
    for( size_t i = 0; i < nchannel; ++i )
      (*counts)[i] = static_cast<float>( 50.0 + 0.02*(2048.0 - cal->energy_for_channel( i )) );

    vector<pair<const SandiaDecay::Nuclide *,double>> nuc_ages;
    for( const Source &src : sources )
      nuc_ages.emplace_back( src.nuc, src.age );
    const CascadeSummingCalc calc( nuc_ages, {}, 1.25, drf, 0.0, false, sm_datadir );

    const std::function<double(double)> eps_fep = [&]( const double e ) -> double {
      return s * drf->farFieldIntrinsicEfficiency( static_cast<float>(e) );
    };
    const std::function<double(double)> eps_tot = [&]( const double e ) -> double {
      return s * drf->totalIntrinsicEfficiencyAny( static_cast<float>(e) );
    };

    for( const Source &src : sources )
    {
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( src.nuc, 1.0, src.age );
      const vector<SandiaDecay::EnergyRatePair> photons = mix.photons( 0.0 );

      vector<SandiaDecay::EnergyRatePair> lines;
      vector<ceelo::PeakWindow> windows;
      for( const SandiaDecay::EnergyRatePair &p : photons )
      {
        bool in_roi = false;
        for( const pair<double,double> &roi : sm_rois )
          in_roi = (in_roi || ((p.energy > (roi.first - 5.0)) && (p.energy < (roi.second + 5.0))));
        if( !in_roi || (p.numPerSecond <= 0.0) )
          continue;

        lines.push_back( p );
        ceelo::PeakWindow window;
        window.energy_keV = p.energy;
        window.tolerance_keV = window_tolerance( drf, p.energy );
        windows.push_back( window );
      }//for( const SandiaDecay::EnergyRatePair &p : photons )

      vector<double> factors( lines.size(), 1.0 );
      if( s > 0.0 )
      {
        const vector<ceelo::AnalyticPeakResultT<double>> results
                                         = calc.evaluateWindows<double>( src.nuc, src.age, windows, eps_fep, eps_tot );
        BOOST_REQUIRE_EQUAL( results.size(), lines.size() );
        for( size_t i = 0; i < results.size(); ++i )
          factors[i] = results[i].found ? results[i].c_net : 1.0;
      }

      for( size_t i = 0; i < lines.size(); ++i )
      {
        const double energy = lines[i].energy;
        true_factors[{src.nuc, energy}] = factors[i];

        const double area = src.rel_act * sm_live_time * drf->farFieldIntrinsicEfficiency( static_cast<float>(energy) )
                            * lines[i].numPerSecond * factors[i];
        const double sigma = drf->peakResolutionFWHM( static_cast<float>(energy) ) / 2.35482;
        const auto cdf = [&]( const double x ) -> double {
          return 0.5*boost::math::erfc( -(x - energy) / (sigma*std::sqrt(2.0)) );
        };

        const size_t first = static_cast<size_t>( cal->channel_for_energy( energy - 8.0*sigma ) );
        const size_t last = static_cast<size_t>( cal->channel_for_energy( energy + 8.0*sigma ) );
        for( size_t channel = first; (channel <= last) && (channel < nchannel); ++channel )
        {
          const double lower = cal->energy_for_channel( channel ), upper = cal->energy_for_channel( channel + 1 );
          (*counts)[channel] += static_cast<float>( area*(cdf( upper ) - cdf( lower )) );
        }
      }//for( size_t i = 0; i < lines.size(); ++i )

      // Pure coincidence sums in the ROIs: candidate energies grouped into windows of the line windows' width.
      if( s > 0.0 )
      {
        const shared_ptr<const ceelo::AnalyticCascadeModel> model = calc.cascadeModel( src.nuc, src.age );
        vector<ceelo::PeakWindow> sum_windows;
        for( const pair<double,double> &roi : sm_rois )
        {
          const vector<double> energies = ceelo::coincidence_sum_energies( *model, roi.first, roi.second );
          for( size_t i = 0; i < energies.size(); )
          {
            const double tol = window_tolerance( drf, energies[i] );
            size_t j = i;
            double total = 0.0;
            while( (j < energies.size()) && ((energies[j] - energies[i]) <= tol) )
              total += energies[j++];
            ceelo::PeakWindow window;
            window.energy_keV = total / static_cast<double>( j - i );
            window.tolerance_keV = window_tolerance( drf, window.energy_keV );
            sum_windows.push_back( window );
            i = j;
          }
        }//for( const pair<double,double> &roi : sm_rois )

        const vector<ceelo::AnalyticPeakResultT<double>> sum_results
                          = calc.evaluateWindows<double>( src.nuc, src.age, sum_windows, eps_fep, eps_tot, true );
        for( const ceelo::AnalyticPeakResultT<double> &r : sum_results )
        {
          if( r.found || !(r.pure_sum_rate > 0.0) )
            continue;

          // Counts = activity x live time x rate; the relative activity carries `s` (its lines use s*intrinsic(E)).
          const double energy = r.pure_sum_energy_keV;
          const double area = (src.rel_act / s) * sm_live_time * r.pure_sum_rate;
          if( true_sums )
            true_sums->push_back( TrueSum{ src.nuc, energy, area } );
          const double sigma = drf->peakResolutionFWHM( static_cast<float>(energy) ) / 2.35482;
          const size_t first = static_cast<size_t>( cal->channel_for_energy( energy - 8.0*sigma ) );
          const size_t last = static_cast<size_t>( cal->channel_for_energy( energy + 8.0*sigma ) );
          for( size_t channel = first; (channel <= last) && (channel < nchannel); ++channel )
          {
            const double lower = cal->energy_for_channel( channel ), upper = cal->energy_for_channel( channel + 1 );
            (*counts)[channel] += static_cast<float>( area * 0.5 * (boost::math::erf( (upper - energy)/(sigma*sqrt(2.0)) )
                                                                  - boost::math::erf( (lower - energy)/(sigma*sqrt(2.0)) )) );
          }
        }//for( const ceelo::AnalyticPeakResultT<double> &r : sum_results )
      }//if( s > 0.0 )
    }//for( const Source &src : sources )

    auto spectrum = make_shared<SpecUtils::Measurement>();
    spectrum->set_gamma_counts( counts, static_cast<float>(sm_live_time), static_cast<float>(sm_live_time) );
    spectrum->set_energy_calibration( cal );
    return spectrum;
  }//make_spectrum(...)


  /** A Physical-model curve with no attenuation and no correction function, so relative efficiency is exactly the
   detector's intrinsic efficiency; FWHM fixed to the detector's. */
  RelActCalcAuto::Options make_options( const vector<Source> &sources, const bool cascade )
  {
    RelActCalcAuto::Options options;
    options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
    options.fwhm_form = RelActCalcAuto::FwhmForm::NotApplicable;
    options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency;
    options.skew_type = PeakDef::SkewType::NoSkew;
    options.additional_br_uncert = 0.0;
    options.auto_profile_weak_mass_fractions = false;
    options.auto_simplify_model = false;

    RelActCalcAuto::RelEffCurveInput curve;
    curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
    curve.rel_eff_eqn_order = 0;
    curve.phys_model_self_atten.reset();
    curve.phys_model_external_atten.clear();
    curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::None;
    curve.correct_for_cascade_summing = cascade;
    for( const Source &src : sources )
    {
      RelActCalcAuto::NucInputInfo input;
      input.source = src.nuc;
      input.age = src.age;
      input.fit_age = false;
      input.peak_color_css = "rgb(0,0,255)";
      curve.nuclides.push_back( input );
    }
    options.rel_eff_curves.push_back( curve );

    for( const pair<double,double> &range : sm_rois )
    {
      RelActCalcAuto::RoiRange roi;
      roi.lower_energy = range.first;
      roi.upper_energy = range.second;
      roi.continuum_type = PeakContinuum::OffsetType::Linear;
      roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
      options.rois.push_back( roi );
    }

    return options;
  }//make_options(...)


  RelActCalcAuto::RelActAutoSolution fit( const RelActCalcAuto::Options &options,
                                          const shared_ptr<const SpecUtils::Measurement> &spectrum,
                                          const shared_ptr<const DetectorPeakResponse> &drf )
  {
    RelActCalcAuto::RelActAutoSolution sol;
    BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, spectrum, nullptr, drf, {},
                                                         PeakFitUtils::CoarseResolutionType::High, nullptr ) );
    BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                           "solve failed: " << sol.m_error_message );
    for( const string &warning : sol.m_warnings )
      BOOST_TEST_MESSAGE( "  warning (summing correction "
                          << (options.rel_eff_curves[0].correct_for_cascade_summing ? "on" : "off") << "): " << warning );
    return sol;
  }//fit(...)


  double rel_act( const RelActCalcAuto::RelActAutoSolution &sol, const SandiaDecay::Nuclide *nuc )
  {
    for( const RelActCalcAuto::NuclideRelAct &act : sol.m_rel_activities.at(0) )
    {
      if( RelActCalcAuto::nuclide( act.source ) == nuc )
        return act.rel_activity;
    }
    BOOST_REQUIRE_MESSAGE( false, "No relative activity for " << nuc->symbol );
    return 0.0;
  }//rel_act(...)


  vector<Source> ba133_cs137()
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    Source ba133, cs137;
    ba133.nuc = db->nuclide( "Ba133" );
    ba133.age = 5.0*PhysicalUnits::year;
    ba133.rel_act = 5000.0;
    cs137.nuc = db->nuclide( "Cs137" );
    cs137.age = 10.0*PhysicalUnits::year;
    cs137.rel_act = 3000.0;
    BOOST_REQUIRE( ba133.nuc && cs137.nuc );
    return { ba133, cs137 };
  }//ba133_cs137()
}//namespace


// The option is off by default; a curve without it serializes exactly as before (same version, no new element), and
//  one with it round-trips.
BOOST_AUTO_TEST_CASE( option_serialization )
{
  set_data_dir();

  RelActCalcAuto::RelEffCurveInput curve;
  BOOST_CHECK( !curve.correct_for_cascade_summing );

  const auto to_string = []( const RelActCalcAuto::RelEffCurveInput &input ) -> string {
    rapidxml::xml_document<char> doc;
    input.toXml( &doc );
    string answer;
    rapidxml::print( std::back_inserter(answer), doc, 0 );
    return answer;
  };

  const string off_xml = to_string( curve );
  BOOST_CHECK( off_xml.find( "CorrectForCascadeSumming" ) == string::npos );
  BOOST_CHECK( off_xml.find( "version=\"2\"" ) == string::npos );

  RelActCalcAuto::RelEffCurveInput on_curve = curve;
  on_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  on_curve.correct_for_cascade_summing = true;
  BOOST_CHECK( !(on_curve == curve) );

  string on_xml = to_string( on_curve );
  BOOST_CHECK( on_xml.find( "CorrectForCascadeSumming" ) != string::npos );

  rapidxml::xml_document<char> doc;
  doc.parse<rapidxml::parse_trim_whitespace>( &on_xml[0] );
  const rapidxml::xml_node<char> * const node = doc.first_node( "RelEffCurveInput" );
  BOOST_REQUIRE( node );
  RelActCalcAuto::RelEffCurveInput read_back;
  BOOST_REQUIRE_NO_THROW( read_back.fromXml( node ) );
  BOOST_CHECK( read_back.correct_for_cascade_summing );
  BOOST_CHECK( read_back == on_curve );
}//option_serialization


// Summing at a close geometry: the corrected fit recovers s, the line factors, and the activity ratio; the uncorrected
//  fit is biased, and fits worse.
BOOST_AUTO_TEST_CASE( synthetic_summing_is_recovered )
{
  set_data_dir();
  const shared_ptr<const DetectorPeakResponse> drf = hpge_drf();
  BOOST_REQUIRE( drf && drf->hasResolutionInfo() );
  BOOST_REQUIRE( GammaInteractionCalc::CascadeSummingCalc::drfHasNeededInfo( drf ) );

  const vector<Source> sources = ba133_cs137();
  const double s_true = 0.2;

  map<pair<const SandiaDecay::Nuclide *,double>,double> true_factors;
  vector<TrueSum> true_sums;
  const shared_ptr<SpecUtils::Measurement> spectrum = make_spectrum( drf, sources, s_true, true_factors, &true_sums );

  // The 356 keV line should lose several percent, so the test has something to recover.
  double c_356 = 1.0;
  for( const auto &entry : true_factors )
  {
    if( (entry.first.first == sources[0].nuc) && (fabs( entry.first.second - 356.0 ) < 0.1) )
      c_356 = entry.second;
  }
  BOOST_TEST_MESSAGE( "True Ba133 356 keV summing factor: " << c_356 );
  BOOST_REQUIRE( c_356 < 0.97 );

  const double true_ratio = sources[0].rel_act / sources[1].rel_act;

  const RelActCalcAuto::RelActAutoSolution with = fit( make_options( sources, true ), spectrum, drf );
  BOOST_REQUIRE( !with.m_phys_model_results.empty() && with.m_phys_model_results[0] );
  BOOST_REQUIRE( with.m_phys_model_results[0]->cascade_summing.has_value() );
  const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::CascadeSummingInfo &info = *with.m_phys_model_results[0]->cascade_summing;

  BOOST_TEST_MESSAGE( "Fit s = " << info.solid_angle_frac << " +- " << info.solid_angle_frac_uncert.value_or( -1.0 )
                      << " (true " << s_true << "), chi2 = " << with.m_chi2 );
  BOOST_CHECK_MESSAGE( fabs( info.solid_angle_frac - s_true ) < 0.02*s_true,
                       "Fit s = " << info.solid_angle_frac << ", expected " << s_true );

  const double with_ratio = rel_act( with, sources[0].nuc ) / rel_act( with, sources[1].nuc );
  BOOST_CHECK_MESSAGE( fabs( with_ratio/true_ratio - 1.0 ) < 0.003,
                       "Corrected Ba133/Cs137 = " << with_ratio << ", expected " << true_ratio );

  size_t num_compared = 0;
  for( const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::CascadeSummingInfo::LineFactor &line : info.lines )
  {
    const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( line.source );
    const auto pos = true_factors.find( {nuc, line.energy} );
    if( pos == end(true_factors) )
      continue;
    ++num_compared;
    BOOST_CHECK_MESSAGE( fabs( line.factor - pos->second ) < (0.05*fabs( 1.0 - pos->second ) + 2.0E-4),
                         "Line " << line.energy << " keV: fit factor " << line.factor << ", true " << pos->second );
  }
  BOOST_CHECK( num_compared >= 4 );

  // The pure coincidence-sum peaks: the fit models those above a small fraction of their ROI's strongest line, each at
  //  about its true area; in particular the strongest true one.
  const TrueSum *strongest = nullptr;
  for( const TrueSum &t : true_sums )
  {
    BOOST_TEST_MESSAGE( "  true sum peak " << t.energy << " keV: " << t.counts << " counts" );
    if( !strongest || (t.counts > strongest->counts) )
      strongest = &t;
  }
  BOOST_REQUIRE( strongest );

  bool found_strongest = false;
  for( const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::CascadeSummingInfo::SumPeak &p : info.sum_peaks )
  {
    const TrueSum *match = nullptr;
    for( const TrueSum &t : true_sums )
      if( fabs( t.energy - p.energy ) < 0.05 )
        match = &t;
    BOOST_TEST_MESSAGE( "  fit sum peak " << p.energy << " keV: " << p.counts << " counts (true "
                        << (match ? match->counts : -1.0) << ")" );
    BOOST_REQUIRE_MESSAGE( match, "Fit sum peak at " << p.energy << " keV is not a true one" );
    BOOST_CHECK_MESSAGE( fabs( p.counts/match->counts - 1.0 ) < 0.03,
                         "Sum peak at " << p.energy << " keV: fit " << p.counts << " vs true " << match->counts );
    found_strongest = (found_strongest || (match == strongest));
  }
  BOOST_CHECK_MESSAGE( found_strongest, "The fit did not model the strongest sum peak, at " << strongest->energy << " keV" );

  const RelActCalcAuto::RelActAutoSolution without = fit( make_options( sources, false ), spectrum, drf );
  BOOST_CHECK( without.m_phys_model_results.empty() || !without.m_phys_model_results[0]
               || !without.m_phys_model_results[0]->cascade_summing.has_value() );
  for( const string &name : without.m_parameter_names )
    BOOST_CHECK( name.find( "CascS" ) == string::npos );
  // The correction adds `s`, plus one yield nuisance per line whose summing could be large.
  BOOST_REQUIRE( !with.m_phys_model_results.empty() && with.m_phys_model_results[0]
                 && with.m_phys_model_results[0]->cascade_summing.has_value() );
  const size_t num_yield_nuisances = with.m_phys_model_results[0]->cascade_summing->yield_multiples.size();
  BOOST_CHECK_EQUAL( without.m_final_parameters.size() + 1 + num_yield_nuisances, with.m_final_parameters.size() );

  const double without_ratio = rel_act( without, sources[0].nuc ) / rel_act( without, sources[1].nuc );
  BOOST_TEST_MESSAGE( "Uncorrected Ba133/Cs137 = " << without_ratio << " (true " << true_ratio << "), chi2 = "
                      << without.m_chi2 );
  BOOST_CHECK( fabs( without_ratio/true_ratio - 1.0 ) > 0.01 );
  BOOST_CHECK( without.m_chi2 > 10.0*std::max( with.m_chi2, 1.0 ) );
}//synthetic_summing_is_recovered


// No summing in the spectrum: the corrected fit finds s near zero, and the same activities as the uncorrected fit.
BOOST_AUTO_TEST_CASE( no_summing_gives_no_correction )
{
  set_data_dir();
  const shared_ptr<const DetectorPeakResponse> drf = hpge_drf();
  const vector<Source> sources = ba133_cs137();

  map<pair<const SandiaDecay::Nuclide *,double>,double> true_factors;
  const shared_ptr<SpecUtils::Measurement> spectrum = make_spectrum( drf, sources, 0.0, true_factors );

  const RelActCalcAuto::RelActAutoSolution with = fit( make_options( sources, true ), spectrum, drf );
  const RelActCalcAuto::RelActAutoSolution without = fit( make_options( sources, false ), spectrum, drf );

  BOOST_REQUIRE( !with.m_phys_model_results.empty() && with.m_phys_model_results[0] );
  BOOST_REQUIRE( with.m_phys_model_results[0]->cascade_summing.has_value() );
  const double s = with.m_phys_model_results[0]->cascade_summing->solid_angle_frac;
  const double s_uncert = with.m_phys_model_results[0]->cascade_summing->solid_angle_frac_uncert.value_or( -1.0 );
  BOOST_TEST_MESSAGE( "Fit s with no summing: " << s << " +- " << s_uncert );
  BOOST_CHECK( s < 0.002 );

  // At s = 0 the data must still constrain s: its uncertainty should be far below the prior's width (0.5) - which
  //  is what it was when the summing factors lost their derivative at that bound.
  BOOST_CHECK_MESSAGE( (s_uncert > 0.0) && (s_uncert < 0.05), "s uncertainty " << s_uncert );

  for( const Source &src : sources )
  {
    const double a = rel_act( with, src.nuc ), b = rel_act( without, src.nuc );
    BOOST_CHECK_MESSAGE( fabs( a/b - 1.0 ) < 1.0E-3, src.nuc->symbol << ": " << a << " vs " << b );
    BOOST_CHECK_MESSAGE( fabs( b/src.rel_act - 1.0 ) < 1.0E-3, src.nuc->symbol << ": " << b << " vs " << src.rel_act );
  }
}//no_summing_gives_no_correction


namespace
{
  /** The test-data directory: `--testfiledir=...`, else found relative to the working directory. */
  string test_file_dir()
  {
    const int argc = framework::master_test_suite().argc;
    char ** const argv = framework::master_test_suite().argv;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        return arg.substr( 14 );
    }
    for( const char *d : { "test_data", "../test_data", "../../test_data", "target/testing/test_data",
                           "../target/testing/test_data" } )
    {
      if( SpecUtils::is_directory( SpecUtils::append_path( d, "rel_act_cascade" ) ) )
        return d;
    }
    return "test_data";
  }//test_file_dir()


  struct McLine
  {
    string nuclide;
    double energy = 0.0, yield = 0.0, eff_with = 0.0;
  };
}//namespace


// Eu152 + Cs137 point source near a micro-Detective, line efficiencies from CeeLo's cascade MC for the detector's real
//  geometry (rel_act_cascade/eu152_cs137_mc_truth.csv; the response is microdet_mc.drf.xml).  The summing-corrected
//  fit should recover the Eu152/Cs137 activity ratio; the uncorrected one is biased near contact.
BOOST_AUTO_TEST_CASE( mc_truth_eu152_near_contact )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  const string dir = SpecUtils::append_path( test_file_dir(), "rel_act_cascade" );
  const string drf_path = SpecUtils::append_path( dir, "microdet_mc.drf.xml" );
  const string truth_path = SpecUtils::append_path( dir, "eu152_cs137_mc_truth.csv" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( drf_path ) && SpecUtils::is_file( truth_path ),
                         "Missing " << drf_path << " or " << truth_path << "; pass -- --testfiledir=..." );

  vector<char> drf_data;
  SpecUtils::load_file_data( drf_path.c_str(), drf_data );
  rapidxml::xml_document<char> drf_doc;
  drf_doc.parse<rapidxml::parse_trim_whitespace>( &drf_data[0] );
  auto drf_mut = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( drf_mut->fromXml( drf_doc.first_node( "DetectorPeakResponse" ) ) );
  const shared_ptr<const DetectorPeakResponse> drf = drf_mut;
  BOOST_REQUIRE( GammaInteractionCalc::CascadeSummingCalc::drfHasNeededInfo( drf ) );

  // scenario ("dist,shield") -> lines
  map<string,vector<McLine>> scenarios;
  {
    ifstream truth( truth_path.c_str() );
    string line;
    getline( truth, line );  //header
    while( getline( truth, line ) )
    {
      vector<string> fields;
      SpecUtils::split( fields, line, "," );
      if( fields.size() < 7 )
        continue;
      McLine l;
      l.nuclide = fields[2];
      l.energy = stod( fields[3] );
      l.yield = stod( fields[4] );
      l.eff_with = stod( fields[6] );
      scenarios[fields[0] + "," + fields[1]].push_back( l );
    }
  }
  BOOST_REQUIRE_EQUAL( scenarios.size(), 3 );

  // ROIs; the 444 and 964 keV ones are left out: each holds two Eu152 lines, which the MC can only estimate from an
  //  isotropic sub-run, too noisy to be truth here.
  const vector<pair<double,double>> rois{ {117,127}, {240,250}, {339,350}, {406,416}, {655,668}, {773,785},
                                          {862,873}, {1080,1095}, {1106,1118}, {1402,1414} };
  const SandiaDecay::Nuclide * const eu152 = db->nuclide( "Eu152" ), * const cs137 = db->nuclide( "Cs137" );
  const double eu_age = 5.0*PhysicalUnits::year, cs_age = 10.0*PhysicalUnits::year, live_time = 3600.0;

  for( const auto &scenario : scenarios )
  {
    const bool fe_shield = (scenario.first.find( "fe3" ) != string::npos);

    // Noiseless spectrum: the strongest Eu152 line gets ~2e6 counts, Cs137 661.7 keV ~5e5.
    double eu_max = 0.0, cs_rate = 0.0;
    for( const McLine &l : scenario.second )
    {
      if( l.nuclide == "Eu152" )
        eu_max = std::max( eu_max, l.yield*l.eff_with );
      else if( fabs( l.energy - 661.66 ) < 0.1 )
        cs_rate = l.yield*l.eff_with;
    }
    BOOST_REQUIRE( (eu_max > 0.0) && (cs_rate > 0.0) );
    const double k_eu = 2.0E6/(live_time*eu_max), k_cs = 5.0E5/(live_time*cs_rate);

    const size_t nchannel = 8192;
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( nchannel, { 0.0f, 0.25f }, {} );
    auto counts = make_shared<vector<float>>( nchannel, 0.0f );
    for( size_t i = 0; i < nchannel; ++i )
      (*counts)[i] = static_cast<float>( 400.0*exp( -cal->energy_for_channel( i )/400.0 ) + 20.0 );
    for( const McLine &l : scenario.second )
    {
      const double area = ((l.nuclide == "Eu152") ? k_eu : k_cs) * live_time * l.yield * l.eff_with;
      const double sigma = drf->peakResolutionFWHM( static_cast<float>(l.energy) ) / 2.35482;
      const size_t c0 = static_cast<size_t>( cal->channel_for_energy( l.energy - 8.0*sigma ) );
      const size_t c1 = static_cast<size_t>( cal->channel_for_energy( l.energy + 8.0*sigma ) );
      for( size_t c = c0; (c <= c1) && (c < nchannel); ++c )
      {
        const double lo = cal->energy_for_channel( c ), hi = cal->energy_for_channel( c + 1 );
        (*counts)[c] += static_cast<float>( area * 0.5 * (boost::math::erf( (hi - l.energy)/(sigma*sqrt(2.0)) )
                                                          - boost::math::erf( (lo - l.energy)/(sigma*sqrt(2.0)) )) );
      }
    }
    auto spectrum = make_shared<SpecUtils::Measurement>();
    spectrum->set_gamma_counts( counts, static_cast<float>(live_time), static_cast<float>(live_time) );
    spectrum->set_energy_calibration( cal );

    double err[2] = { 0.0, 0.0 }, s = -1.0;
    for( const bool cascade : { false, true } )
    {
      RelActCalcAuto::Options options;
      options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
      options.fwhm_form = RelActCalcAuto::FwhmForm::NotApplicable;
      options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency;
      options.skew_type = PeakDef::SkewType::NoSkew;
      options.additional_br_uncert = 0.0;
      options.auto_profile_weak_mass_fractions = false;
      options.auto_simplify_model = false;

      RelActCalcAuto::RelEffCurveInput curve;
      curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      curve.rel_eff_eqn_order = 0;
      curve.phys_model_self_atten.reset();
      curve.phys_model_external_atten.clear();
      if( fe_shield )
      {
        auto shield = make_shared<RelActCalc::PhysicalModelShieldInput>();
        shield->material = matdb->material( "Fe" );
        shield->areal_density = 1.0*PhysicalUnits::g_per_cm2;
        shield->fit_areal_density = true;
        shield->lower_fit_areal_density = 0.0;
        shield->upper_fit_areal_density = 20.0*PhysicalUnits::g_per_cm2;
        curve.phys_model_external_atten.push_back( shield );
      }
      curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::Hoerl;
      curve.correct_for_cascade_summing = cascade;
      for( const pair<const SandiaDecay::Nuclide *,double> &na : { make_pair(eu152, eu_age), make_pair(cs137, cs_age) } )
      {
        RelActCalcAuto::NucInputInfo input;
        input.source = na.first;
        input.age = na.second;
        input.fit_age = false;
        input.peak_color_css = "rgb(0,0,255)";
        curve.nuclides.push_back( input );
      }
      options.rel_eff_curves.push_back( curve );
      for( const pair<double,double> &r : rois )
      {
        RelActCalcAuto::RoiRange roi;
        roi.lower_energy = r.first;
        roi.upper_energy = r.second;
        roi.continuum_type = PeakContinuum::OffsetType::Linear;
        roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
        options.rois.push_back( roi );
      }

      const RelActCalcAuto::RelActAutoSolution sol = fit( options, spectrum, drf );
      const double ratio = rel_act( sol, eu152 ) / rel_act( sol, cs137 );
      err[cascade ? 1 : 0] = 100.0*(ratio/(k_eu/k_cs) - 1.0);
      if( cascade )
      {
        BOOST_REQUIRE( !sol.m_phys_model_results.empty() && sol.m_phys_model_results[0] );
        BOOST_REQUIRE( sol.m_phys_model_results[0]->cascade_summing.has_value() );
        s = sol.m_phys_model_results[0]->cascade_summing->solid_angle_frac;
        double sum_counts = 0.0;
        for( const auto &p : sol.m_phys_model_results[0]->cascade_summing->sum_peaks )
          sum_counts += p.counts;
        BOOST_TEST_MESSAGE( "  " << scenario.first << ": " << sol.m_phys_model_results[0]->cascade_summing->sum_peaks.size()
                            << " pure-sum peaks fit, " << sum_counts << " counts in all; chi2 " << sol.m_chi2 );
      }else
      {
        BOOST_TEST_MESSAGE( "  " << scenario.first << ": uncorrected chi2 " << sol.m_chi2 );
      }
    }//for( cascade off / on )

    BOOST_TEST_MESSAGE( "MC truth " << scenario.first << ": Eu152/Cs137 error " << err[0] << "% uncorrected, " << err[1]
                        << "% corrected (s = " << s << ")" );

    const bool contact = (scenario.first.rfind( "1,", 0 ) == 0);
    BOOST_CHECK_MESSAGE( fabs( err[1] ) < (contact ? 2.0 : 1.5), scenario.first << ": corrected error " << err[1] << "%" );
    if( contact )
    {
      BOOST_CHECK_MESSAGE( fabs( err[1] ) < 0.5*fabs( err[0] ),
                           scenario.first << ": correction " << err[0] << "% -> " << err[1] << "%" );
      BOOST_CHECK_MESSAGE( (s > 0.1) && (s < 0.3), scenario.first << ": s = " << s );
    }else
    {
      BOOST_CHECK_MESSAGE( (s >= 0.0) && (s < 0.06), scenario.first << ": s = " << s );
    }
  }//for( const auto &scenario : scenarios )
}//mc_truth_eu152_near_contact
