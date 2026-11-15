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

// The cascade-summing kernel of the "Physical" relative-efficiency model (RelActCalcCascade), against
// CeeLo's analytic engine (through CascadeSummingCalc), closed forms, and finite differences.
//
//   ./test_RelActCalcCascade -- --datadir=<InterSpec data> --testfiledir=<target/testing/test_data>

#include "InterSpec_config.h"

#include <map>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <optional>
#include <algorithm>
#include <functional>

#define BOOST_TEST_MODULE RelActCalcCascade_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcCascade.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/RelActCalcCascade_imp.hpp"

using namespace std;
using namespace boost::unit_test;

namespace
{
string g_data_dir, g_test_data_dir;

void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      g_data_dir = arg.substr( 10 );
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_data_dir = arg.substr( 14 );
  }
  SpecUtils::ireplace_all( g_data_dir, "%20", " " );
  SpecUtils::ireplace_all( g_test_data_dir, "%20", " " );

  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty() && !g_test_data_dir.empty(),
                         "Need --datadir=... and --testfiledir=... after the '--'" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
}//set_data_dir()


/** An example HPGe response from a CeeLo Monte-Carlo model: a 50 x 30 mm crystal, with a characterized total
 efficiency and its own near-field efficiencies. */
shared_ptr<const DetectorPeakResponse> example_hpge_mc()
{
  set_data_dir();
  const string path = SpecUtils::append_path( SpecUtils::append_path( g_test_data_dir, "rel_act_cascade" ),
                                              "example_hpge_mc.drf.xml" );
  vector<char> data;
  SpecUtils::load_file_data( path.c_str(), data );
  rapidxml::xml_document<char> doc;
  doc.parse<rapidxml::parse_trim_whitespace>( &data[0] );
  const rapidxml::xml_node<char> * const root = doc.first_node( "DetectorPeakResponse" );
  BOOST_REQUIRE( root );
  shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>();
  drf->fromXml( root );
  BOOST_REQUIRE( drf->isValid() );
  return drf;
}//example_hpge_mc()


const SandiaDecay::Nuclide *nuclide( const string &name )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
  BOOST_REQUIRE_MESSAGE( nuc, "no nuclide " + name );
  return nuc;
}


/** The gamma lines of `nuc` at `age` emitted more often than `min_yield` per decay, above `min_energy`. */
vector<double> strong_lines( const SandiaDecay::Nuclide * const nuc, const double age, const double min_yield,
                             const double min_energy )
{
  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, 1.0, age );
  vector<double> energies;
  for( const SandiaDecay::EnergyRatePair &erp : mix.gammas( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy, false ) )
  {
    if( (erp.numPerSecond >= min_yield) && (erp.energy >= min_energy)
       && (energies.empty() || (erp.energy != energies.back())) )
      energies.push_back( erp.energy );
  }
  return energies;
}//strong_lines(...)


struct TestSource
{
  const char *name;
  double age_years;
  double min_yield;
};

const TestSource ns_sources[] = {
  { "U235", 20.0, 1.0E-5 }, { "U238", 20.0, 1.0E-6 }, { "U232", 20.0, 1.0E-4 }, { "Eu152", 1.0, 1.0E-3 },
  { "Ba133", 1.0, 1.0E-3 }, { "Co60", 1.0, 1.0E-3 }
};


/** Evaluates the model for a point source at a front solid angle, with one external attenuator. */
template<typename T>
RelActCalcCascade::EvalResult<T> evaluate_model( const RelActCalcCascade::SummingModel &model, const T &solid_angle,
                                                 const optional<RelActCalc::PhysModelShield<T>> &self_atten,
                                                 const vector<RelActCalc::PhysModelShield<T>> &ext_atten )
{
  RelActCalcCascade::EvalInput<T> input;
  input.front_solid_angle = solid_angle;
  model.set_attenuation( self_atten, ext_atten, input );
  RelActCalcCascade::EvalResult<T> result;
  model.evaluate( input, result );
  return result;
}//evaluate_model(...)
}//namespace


BOOST_AUTO_TEST_CASE( solid_angle_and_distance )
{
  const double radius = 2.5*PhysicalUnits::cm;
  for( const double dist_cm : { 0.01, 0.5, 2.0, 5.0, 15.0, 100.0, 1.0E4 } )
  {
    const double dist = dist_cm*PhysicalUnits::cm;
    const double g = RelActCalcCascade::on_axis_solid_angle( dist, radius );
    BOOST_CHECK_CLOSE( g, DetectorPeakResponse::fractionalSolidAngle( 2.0*radius, dist ), 1.0E-6 );
    BOOST_CHECK_CLOSE( RelActCalcCascade::on_axis_distance( g, radius ), dist, 1.0E-5 );
  }
}//solid_angle_and_distance


BOOST_AUTO_TEST_CASE( front_solid_angle_of_a_wide_source )
{
  set_data_dir();
  const SandiaDecay::Nuclide * const nuc = nuclide( "Co60" );
  RelActCalcCascade::SourceInput input;
  input.nuclide = nuc;
  input.age = PhysicalUnits::year;
  input.line_energies = { 1173.228, 1332.492 };

  // A library response (a flat disk: 5 cm diameter, 0.5 cm setback), and the geometry-aware Monte-Carlo one.
  shared_ptr<DetectorPeakResponse> flat;
  {
    ifstream tsv( SpecUtils::append_path( g_data_dir, "common_drfs.tsv" ).c_str(), ios::binary );
    vector<string> credits;
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    DetectorPeakResponse::parseMultipleRelEffDrfCsv( tsv, credits, drfs );
    for( const shared_ptr<DetectorPeakResponse> &d : drfs )
    {
      if( d && (d->name().find( "MicroDetective_LANL_025cm" ) != string::npos) && !flat )
        flat = d;
    }
  }
  BOOST_REQUIRE( flat );

  for( const shared_ptr<const DetectorPeakResponse> &drf : { shared_ptr<const DetectorPeakResponse>( flat ), example_hpge_mc() } )
  {
    const RelActCalcCascade::SummingModel model( { input }, drf );
    const double source_radius = 3.5*PhysicalUnits::cm;   //the CBNM / NRM-171 can
    double previous = 1.0;
    for( const double dist_cm : { 0.5, 2.0, 5.0, 15.0, 100.0 } )
    {
      const double dist = dist_cm*PhysicalUnits::cm;
      const double point = model.front_solid_angle( dist ), wide = model.front_solid_angle( dist, source_radius );
      BOOST_CHECK_GT( point, 0.0 );
      BOOST_CHECK_LE( point, model.max_front_solid_angle()*(1.0 + 1.0E-9) );
      BOOST_CHECK_LT( point, previous );    //falls with distance
      BOOST_CHECK_LT( wide, point );        //a wide source sees less of the detector
      previous = point;

      BOOST_CHECK_CLOSE( model.front_distance( point )/PhysicalUnits::cm, dist_cm, 0.1 );
      BOOST_CHECK_CLOSE( model.front_distance( wide, source_radius )/PhysicalUnits::cm, dist_cm, 0.1 );

      BOOST_TEST_MESSAGE( drf->name() << " at " << dist_cm << " cm: point " << point << ", 7 cm disk " << wide
                          << " (" << wide/point << " of it)" );
    }

    // Far away the source's width does not matter.
    BOOST_CHECK_CLOSE( model.front_solid_angle( 100.0*PhysicalUnits::cm, source_radius ),
                       model.front_solid_angle( 100.0*PhysicalUnits::cm ), 0.2 );
  }

  // For the flat disk the 7 cm can at 2 cm sees 0.69 of the on-axis solid angle (the fraction the JRC 2 cm
  //  line trends were measured to need: 0.68 +- 0.13).
  const RelActCalcCascade::SummingModel flat_model( { input }, flat );
  BOOST_CHECK_CLOSE( flat_model.front_solid_angle( 2.0*PhysicalUnits::cm, 3.5*PhysicalUnits::cm )
                     / flat_model.front_solid_angle( 2.0*PhysicalUnits::cm ), 0.692, 1.0 );
}//front_solid_angle_of_a_wide_source


BOOST_AUTO_TEST_CASE( depth_nodes_resolve_attenuation )
{
  // The slab integrals the kernel forms: int_0^1 exp(-x t) dt, and the pair ratio with x1 + x2 over x1.
  const RelActCalcCascade::DepthNodes &nodes = RelActCalcCascade::depth_nodes();
  double weight_sum = 0.0;
  for( const double w : nodes.weight )
    weight_sum += w;
  BOOST_CHECK_CLOSE( weight_sum, 1.0, 1.0E-10 );

  auto quad = [&nodes]( const double x ) -> double {
    double sum = 0.0;
    for( size_t i = 0; i < nodes.fraction.size(); ++i )
      sum += nodes.weight[i] * exp( -x*nodes.fraction[i] );
    return sum;
  };
  auto exact = []( const double x ) -> double { return (x < 1.0E-8) ? (1.0 - 0.5*x) : (-expm1( -x ) / x); };

  double worst = 0.0, worst_ratio = 0.0;
  for( double x = 0.0; x <= 2000.0; x = (x < 1.0) ? (x + 0.1) : (x*1.07) )
  {
    worst = max( worst, fabs( quad( x )/exact( x ) - 1.0 ) );
    for( const double x2 : { 0.05, 1.0, 20.0 } )
      worst_ratio = max( worst_ratio, fabs( (quad( x + x2 )/quad( x )) / (exact( x + x2 )/exact( x )) - 1.0 ) );
  }
  BOOST_TEST_MESSAGE( nodes.fraction.size() << " depth nodes: worst relative error " << worst
                      << ", of a pair ratio " << worst_ratio );
  BOOST_CHECK_LT( worst, 5.0E-4 );
  BOOST_CHECK_LT( worst_ratio, 5.0E-4 );
}//depth_nodes_resolve_attenuation


BOOST_AUTO_TEST_CASE( total_efficiency_policy )
{
  set_data_dir();
  const vector<double> energies{ 20.0, 93.35, 185.71, 766.37, 1001.0, 2614.5 };

  // A response with a characterized total efficiency uses it.
  const shared_ptr<const DetectorPeakResponse> mc = example_hpge_mc();
  const RelActCalcCascade::DetectorEfficiencies mc_eff = RelActCalcCascade::detector_efficiencies( mc, energies );
  BOOST_CHECK( mc_eff.total_source == RelActCalcCascade::TotalEffSource::DetectorResponse );
  BOOST_CHECK( mc_eff.total_source_reason.empty() );
  BOOST_CHECK_GT( mc_eff.detector_diameter, 0.0 );
  for( size_t i = 0; i < energies.size(); ++i )
  {
    BOOST_CHECK_GE( mc_eff.total[i], mc_eff.fep[i] );
    BOOST_CHECK_LE( mc_eff.total[i], 1.0 );
  }
  BOOST_CHECK_GT( mc_eff.fep[2], 0.0 );

  // GADRAS' generic 40% HPGe states a total efficiency above one: not used.
  shared_ptr<DetectorPeakResponse> gadras = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( gadras->fromGadrasDirectory( SpecUtils::append_path(
                             SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" ), "HPGe 40%" ) ) );
  const RelActCalcCascade::DetectorEfficiencies gad_eff = RelActCalcCascade::detector_efficiencies( gadras, energies );
  BOOST_TEST_MESSAGE( "GADRAS HPGe 40%: " << (gad_eff.total_source_reason.empty() ? string("uses its own total")
                                                                                   : gad_eff.total_source_reason) );
  for( size_t i = 0; i < energies.size(); ++i )
  {
    BOOST_CHECK_GE( gad_eff.total[i], gad_eff.fep[i] );
    BOOST_CHECK_LE( gad_eff.total[i], 1.0 );
  }

  // A library response has no total efficiency: the example detector's total-to-peak ratio is used.
  shared_ptr<DetectorPeakResponse> library;
  {
    ifstream tsv( SpecUtils::append_path( g_data_dir, "common_drfs.tsv" ).c_str(), ios::binary );
    BOOST_REQUIRE( tsv.is_open() );
    vector<string> credits;
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    DetectorPeakResponse::parseMultipleRelEffDrfCsv( tsv, credits, drfs );
    for( const shared_ptr<DetectorPeakResponse> &d : drfs )
    {
      if( d && (d->name().find( "Detective_LANL_025cm" ) != string::npos) && !library )
        library = d;
    }
  }
  BOOST_REQUIRE( library );
  const RelActCalcCascade::DetectorEfficiencies lib_eff = RelActCalcCascade::detector_efficiencies( library, energies );
  BOOST_CHECK( lib_eff.total_source == RelActCalcCascade::TotalEffSource::ExampleDetectorRatio );
  BOOST_CHECK( !lib_eff.total_source_reason.empty() );
  BOOST_CHECK_EQUAL( lib_eff.fep[0], 0.0 );    //20 keV is below anything the response is good for
  BOOST_CHECK_EQUAL( lib_eff.total[0], 0.0 );
  for( size_t i = 1; i < energies.size(); ++i )
  {
    BOOST_CHECK_GT( lib_eff.fep[i], 0.0 );
    BOOST_CHECK_GT( lib_eff.total[i], lib_eff.fep[i] );
    BOOST_CHECK_LE( lib_eff.total[i], 1.0 );
  }
  // Total-to-peak grows with energy.
  BOOST_CHECK_GT( lib_eff.total[4]/lib_eff.fep[4], lib_eff.total[2]/lib_eff.fep[2] );

  // The library Detective is a small (5 cm) crystal, and its full-energy efficiency falls with energy like
  //  one: its fall-back total is the small reference crystal's - close to what the Monte-Carlo model of a
  //  crystal that size has - not the example detector's, which is less than half of that at 766 keV.
  const double lib_ratio_766 = lib_eff.total[3]/lib_eff.fep[3], mc_ratio_766 = mc_eff.total[3]/mc_eff.fep[3];
  BOOST_TEST_MESSAGE( "Total-to-peak at 766 keV: library Detective (fall-back) " << lib_ratio_766
                      << ", Monte-Carlo 50 x 30 mm HPGe " << mc_ratio_766 );
  BOOST_CHECK_CLOSE( lib_ratio_766, mc_ratio_766, 15.0 );

  // A large crystal gets a much smaller ratio.
  shared_ptr<DetectorPeakResponse> big;
  {
    ifstream tsv( SpecUtils::append_path( g_data_dir, "common_drfs.tsv" ).c_str(), ios::binary );
    vector<string> credits;
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    DetectorPeakResponse::parseMultipleRelEffDrfCsv( tsv, credits, drfs );
    for( const shared_ptr<DetectorPeakResponse> &d : drfs )
    {
      if( d && (d->name().find( "Cx116%" ) != string::npos) && !big )
        big = d;
    }
  }
  BOOST_REQUIRE( big );
  const RelActCalcCascade::DetectorEfficiencies big_eff = RelActCalcCascade::detector_efficiencies( big, energies );
  BOOST_TEST_MESSAGE( "Total-to-peak at 766 keV: 116% coax (fall-back) " << big_eff.total[3]/big_eff.fep[3] );
  BOOST_CHECK_LT( big_eff.total[3]/big_eff.fep[3], 0.6*lib_ratio_766 );

  BOOST_CHECK_THROW( RelActCalcCascade::detector_efficiencies( nullptr, energies ), std::exception );
}//total_efficiency_policy


BOOST_AUTO_TEST_CASE( point_source_matches_engine )
{
  // No self-attenuation: one node, so the kernel's factor is the engine's summing-out at the same efficiencies.
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  using GammaInteractionCalc::CascadeSummingCalc;
  BOOST_REQUIRE( CascadeSummingCalc::drfHasNeededInfo( drf ) );

  for( const TestSource &src : ns_sources )
  {
    const SandiaDecay::Nuclide * const nuc = nuclide( src.name );
    const double age = src.age_years*PhysicalUnits::year;

    RelActCalcCascade::SourceInput input;
    input.nuclide = nuc;
    input.age = age;
    // Only lines the response is characterized for: below its range the full-energy efficiency is taken as
    //  zero, where the engine reports no correction at all.
    for( const double energy : strong_lines( nuc, age, src.min_yield, 50.0 ) )
    {
      if( RelActCalcCascade::detector_efficiencies( drf, { energy } ).fep[0] > 0.0 )
        input.line_energies.push_back( energy );
    }
    BOOST_REQUIRE( !input.line_energies.empty() );

    const RelActCalcCascade::SummingModel model( { input }, drf );
    const RelActCalcCascade::SummingModel::Source &source = model.sources()[0];

    // The oracle, at exactly the kernel's efficiencies, and separately with nothing pruned.
    const vector<double> &energies = model.energies();
    const RelActCalcCascade::DetectorEfficiencies &det = model.detector();

    vector<pair<double,double>> peak_widths;
    vector<ceelo::PeakWindow> windows;
    for( const double e : input.line_energies )
    {
      peak_widths.emplace_back( e, 1.0 );
      ceelo::PeakWindow w;
      w.energy_keV = e;
      w.tolerance_keV = 1.0E-3;
      windows.push_back( w );
    }
    const CascadeSummingCalc calc( { { nuc, age } }, peak_widths, 1.0, drf, 0.0, false, g_data_dir );

    for( const double dist_cm : { 0.0, 2.0, 10.0 } )
    {
      const double g0 = model.front_solid_angle( dist_cm*PhysicalUnits::cm );
      const RelActCalcCascade::EvalResult<double> result
                                         = evaluate_model<double>( model, g0, nullopt, {} );

      // The kernel's own point efficiencies (a geometry-aware response's near-field ones).
      vector<double> fep_abs, tot_abs;
      model.point_efficiencies( g0, fep_abs, tot_abs );

      const function<double(double)> fep_fcn = [&]( double e ) -> double {
        const int index = CascadeSummingCalc::nearestEnergyIndex( energies, e );
        return ((index >= 0) && (fabs( energies[index] - e ) < 1.0E-6)) ? fep_abs[index] : 0.0;
      };
      const function<double(double)> tot_fcn = [&]( double e ) -> double {
        const int index = CascadeSummingCalc::nearestEnergyIndex( energies, e );
        return ((index >= 0) && (fabs( energies[index] - e ) < 1.0E-6)) ? tot_abs[index] : 0.0;
      };
      const function<double(double)> tot_full_fcn = [&]( double e ) -> double {
        const int index = CascadeSummingCalc::nearestEnergyIndex( energies, e );
        if( (index >= 0) && (fabs( energies[index] - e ) < 1.0E-6) )
          return tot_abs[index];
        return g0*RelActCalcCascade::detector_efficiencies( drf, { e } ).total[0];   //a pruned partner
      };

      const vector<ceelo::AnalyticPeakResultT<double>> same
                                            = calc.evaluateWindows<double>( nuc, age, windows, fep_fcn, tot_fcn );
      const vector<ceelo::AnalyticPeakResultT<double>> full
                                       = calc.evaluateWindows<double>( nuc, age, windows, fep_fcn, tot_full_fcn );

      double worst_same = 0.0, worst_full = 0.0, largest_loss = 0.0;
      size_t nfound = 0;
      for( size_t line = 0; line < input.line_energies.size(); ++line )
      {
        const double kernel = result.summing_out[0][line];
        if( !same[line].found )
        {
          BOOST_CHECK( !source.line_found[line] );
          BOOST_CHECK_EQUAL( kernel, 1.0 );
          continue;
        }
        ++nfound;
        BOOST_CHECK( source.line_found[line] );
        worst_same = max( worst_same, fabs( kernel - same[line].c_out ) );
        worst_full = max( worst_full, fabs( kernel - full[line].c_out ) );
        largest_loss = max( largest_loss, 1.0 - kernel );
        // "Same efficiencies" is not quite the same: a photon too weak to matter to this nuclide's lines is
        //  left out of its program, even where the model tabulates that energy for another reason.
        BOOST_CHECK_MESSAGE( fabs( kernel - same[line].c_out ) < 2.0E-6,
                             src.name << " " << input.line_energies[line] << " keV at " << dist_cm << " cm: kernel "
                             << kernel << ", engine " << same[line].c_out );
      }
      BOOST_TEST_MESSAGE( src.name << " at " << dist_cm << " cm (solid angle " << g0 << "): " << nfound << " of "
                          << input.line_energies.size() << " lines; largest loss " << largest_loss
                          << "; worst |diff| " << worst_same << " (same efficiencies), " << worst_full
                          << " (nothing pruned); " << energies.size() << " energies kept" );
      BOOST_CHECK_GT( nfound, 0u );
      // Pruned photons each change a factor by under 1e-7; there are not many of them.
      BOOST_CHECK_LT( worst_full, 5.0E-6 );
      if( dist_cm == 0.0 )
        BOOST_CHECK_GT( largest_loss, 0.01 );
    }//for( loop over distances )
  }//for( loop over sources )
}//point_source_matches_engine


BOOST_AUTO_TEST_CASE( slab_closed_form_co60 )
{
  // Co-60: each of the two gammas is followed by the other in (nearly) every decay, so for a slab, with only
  //  unscattered photons counted, the loss of line 1 is g0 * eps_tot(2) * S(mu1+mu2)/S(mu1),
  //  S(mu) = (1 - exp(-mu AD))/(mu AD).
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  const SandiaDecay::Nuclide * const nuc = nuclide( "Co60" );
  const double age = PhysicalUnits::year;

  RelActCalcCascade::SourceInput input;
  input.nuclide = nuc;
  input.age = age;
  input.line_energies = strong_lines( nuc, age, 0.5, 1000.0 );
  BOOST_REQUIRE_EQUAL( input.line_energies.size(), 2u );

  const RelActCalcCascade::SummingModel model( { input }, drf, false );
  BOOST_CHECK( !model.models_shield_scatter() );
  const RelActCalcCascade::SummingModel scatter_model( { input }, drf );
  BOOST_REQUIRE_MESSAGE( scatter_model.models_shield_scatter(), "GADRAS shield-scatter data not loaded" );
  const shared_ptr<const Material> lead = MaterialDB::instance()->material( "Pb (lead)" );
  BOOST_REQUIRE( lead );
  const shared_ptr<const vector<double>> mu = model.material_mu( lead );

  const vector<double> &energies = model.energies();
  const int i1 = GammaInteractionCalc::CascadeSummingCalc::nearestEnergyIndex( energies, input.line_energies[0] );
  const int i2 = GammaInteractionCalc::CascadeSummingCalc::nearestEnergyIndex( energies, input.line_energies[1] );
  BOOST_REQUIRE( (i1 >= 0) && (i2 >= 0) );

  auto slab = []( const double x ) -> double { return (x < 1.0E-8) ? (1.0 - 0.5*x) : (-expm1( -x )/x); };

  const double g0 = model.front_solid_angle( 1.0*PhysicalUnits::cm );
  const RelActCalcCascade::EvalResult<double> thin = evaluate_model<double>( model, g0, nullopt, {} );
  const double thin_loss_1 = 1.0 - thin.summing_out[0][0];
  const double thin_loss_2 = 1.0 - thin.summing_out[0][1];
  BOOST_CHECK_GT( thin_loss_1, 0.01 );

  for( const double ad_g_cm2 : { 0.0, 1.0, 20.0, 150.0 } )
  {
    RelActCalc::PhysModelShield<double> self;
    self.material = lead;
    self.areal_density = ad_g_cm2*PhysicalUnits::g_per_cm2;
    const RelActCalcCascade::EvalResult<double> result
        = evaluate_model<double>( model, g0, optional<RelActCalc::PhysModelShield<double>>( self ), {} );

    const double x1 = (*mu)[i1]*self.areal_density, x2 = (*mu)[i2]*self.areal_density;
    const double expect_1 = thin_loss_1 * slab( x1 + x2 )/slab( x1 );
    const double expect_2 = thin_loss_2 * slab( x1 + x2 )/slab( x2 );
    BOOST_CHECK_SMALL( (1.0 - result.summing_out[0][0]) - expect_1, 1.0E-3*thin_loss_1 );
    BOOST_CHECK_SMALL( (1.0 - result.summing_out[0][1]) - expect_2, 1.0E-3*thin_loss_2 );

    // Photons that scatter in the source and still reach the detector are partners too: with them the loss
    //  can only be larger - but they are a correction to it, not the bulk of it.
    const RelActCalcCascade::EvalResult<double> scattered
        = evaluate_model<double>( scatter_model, g0, optional<RelActCalc::PhysModelShield<double>>( self ), {} );
    const double loss = 1.0 - result.summing_out[0][0], scatter_loss = 1.0 - scattered.summing_out[0][0];
    BOOST_CHECK_GE( scatter_loss, loss*(1.0 - 1.0E-9) );
    BOOST_CHECK_LE( scatter_loss, 1.6*loss + 1.0E-12 );
    if( ad_g_cm2 == 0.0 )
      BOOST_CHECK_SMALL( scatter_loss - loss, 1.0E-12 );
    else
      BOOST_CHECK_GT( scatter_loss, loss*1.001 );
    BOOST_TEST_MESSAGE( "Co-60 in " << ad_g_cm2 << " g/cm2 of lead: loss " << loss << " unscattered only, "
                        << scatter_loss << " with scattered photons" );
  }//for( loop over areal densities )
}//slab_closed_form_co60


BOOST_AUTO_TEST_CASE( derivatives_match_finite_differences )
{
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  const SandiaDecay::Nuclide * const u235 = nuclide( "U235" );
  const SandiaDecay::Nuclide * const u238 = nuclide( "U238" );
  const double age = 20.0*PhysicalUnits::year;

  vector<RelActCalcCascade::SourceInput> inputs( 2 );
  inputs[0].nuclide = u235;
  inputs[0].age = age;
  inputs[0].line_energies = strong_lines( u235, age, 1.0E-4, 120.0 );
  inputs[0].sum_ranges.push_back( { 270.0, 400.0, 0.0 } );
  inputs[1].nuclide = u238;
  inputs[1].age = age;
  inputs[1].line_energies = strong_lines( u238, age, 1.0E-5, 120.0 );
  const RelActCalcCascade::SummingModel model( inputs, drf );
  BOOST_REQUIRE( !model.sources()[0].sum_peaks.empty() );

  const shared_ptr<const Material> uo3 = MaterialDB::instance()->material( "Uranium trioxide" );
  const shared_ptr<const Material> iron = MaterialDB::instance()->material( "Fe (iron)" );
  BOOST_REQUIRE( uo3 && iron );

  // Parameters: front solid angle, self areal density, external (iron) areal density, generic shield AN and AD.
  const double x0[5] = { 0.06, 5.0*PhysicalUnits::g_per_cm2, 1.5*PhysicalUnits::g_per_cm2, 48.0,
                         0.7*PhysicalUnits::g_per_cm2 };

  {
    auto run = [&]( const auto &x ) {
      using T = std::decay_t<decltype(x[0])>;
      RelActCalc::PhysModelShield<T> self, ext, generic;
      self.material = uo3;
      self.areal_density = x[1];
      ext.material = iron;
      ext.areal_density = x[2];
      generic.atomic_number = x[3];
      generic.areal_density = x[4];
      return evaluate_model<T>( model, x[0], optional<RelActCalc::PhysModelShield<T>>( self ), { ext, generic } );
    };

    typedef ceres::Jet<double,5> Jet;
    vector<Jet> xj( 5 );
    for( int i = 0; i < 5; ++i )
      xj[i] = Jet( x0[i], i );
    const RelActCalcCascade::EvalResult<Jet> jets = run( xj );

    double worst = 0.0;
    for( int par = 0; par < 5; ++par )
    {
      const double h = 1.0E-5*max( fabs( x0[par] ), 1.0E-3 );
      vector<double> up( x0, x0 + 5 ), down( x0, x0 + 5 );
      up[par] += h;
      down[par] -= h;
      const RelActCalcCascade::EvalResult<double> rup = run( up ), rdown = run( down );

      for( size_t src = 0; src < jets.summing_out.size(); ++src )
      {
        for( size_t line = 0; line < jets.summing_out[src].size(); ++line )
        {
          const double fd = (rup.summing_out[src][line] - rdown.summing_out[src][line]) / (2.0*h);
          const double ad = jets.summing_out[src][line].v[par];
          const double scale = max( fabs( fd ), 1.0E-4/max( fabs( x0[par] ), 1.0E-3 ) );
          worst = max( worst, fabs( fd - ad )/scale );
          // The attenuation of a generic atomic number is interpolated between elements, whose derivative
          //  a finite difference only resolves to a part in a thousand.
          const double tol = (par == 3) ? 2.0E-3 : 2.0E-4;
          BOOST_CHECK_MESSAGE( fabs( fd - ad ) <= tol*scale, "line " << line << " of source " << src
                               << ", parameter " << par << ": finite difference " << fd << ", Jet " << ad );
        }

        for( size_t p = 0; p < jets.sum_pairs[src].size(); ++p )
        {
          for( size_t q = 0; q < jets.sum_pairs[src][p].size(); ++q )
          {
            const double fd = (rup.sum_pairs[src][p][q].b_given_a - rdown.sum_pairs[src][p][q].b_given_a) / (2.0*h);
            const double ad = jets.sum_pairs[src][p][q].b_given_a.v[par];
            const double value = jets.sum_pairs[src][p][q].b_given_a.a;
            const double scale = max( fabs( fd ), 1.0E-3*fabs( value )/max( fabs( x0[par] ), 1.0E-3 ) );
            if( scale > 0.0 )
              BOOST_CHECK_MESSAGE( fabs( fd - ad ) <= 3.0E-3*scale, "sum pair " << p << "/" << q << ", parameter "
                                   << par << ": finite difference " << fd << ", Jet " << ad );
          }
        }
      }//for( loop over sources )
    }//for( loop over parameters )
    BOOST_TEST_MESSAGE( "worst relative derivative difference " << worst );
  }

  // The factors stay finite, with a finite solid-angle derivative, as the solid angle goes to zero.
  typedef ceres::Jet<double,1> Jet1;
  RelActCalc::PhysModelShield<Jet1> self;
  self.material = uo3;
  self.areal_density = Jet1( x0[1] );
  {
    const RelActCalcCascade::EvalResult<Jet1> zero
        = evaluate_model<Jet1>( model, Jet1( 0.0, 0 ), optional<RelActCalc::PhysModelShield<Jet1>>( self ), {} );
    double steepest = 0.0;
    for( const vector<Jet1> &lines : zero.summing_out )
    {
      for( const Jet1 &factor : lines )
      {
        BOOST_CHECK_SMALL( factor.a - 1.0, 1.0E-9 );
        BOOST_CHECK( std::isfinite( factor.v[0] ) );
        BOOST_CHECK_LE( factor.v[0], 1.0E-9 );   //more solid angle never adds counts to a line
        steepest = min( steepest, factor.v[0] );
      }
    }
    BOOST_CHECK_LT( steepest, -0.1 );  //and some line does lose them
  }
}//derivatives_match_finite_differences


BOOST_AUTO_TEST_CASE( sum_peaks_match_engine )
{
  // A point source: the kernel's sum-peak rate is sum over pairs of weight * eps_fep(a) * eps_fep(b) * survival,
  //  with an estimated coincident survival; the engine's pure-sum rate has the exact one (and triples).
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  using GammaInteractionCalc::CascadeSummingCalc;

  for( const char * const name : { "U235", "Eu152", "Ba133", "Co60" } )
  {
    const SandiaDecay::Nuclide * const nuc = nuclide( name );
    const double age = ((string(name) == "U235") ? 20.0 : 1.0)*PhysicalUnits::year;

    RelActCalcCascade::SourceInput input;
    input.nuclide = nuc;
    input.age = age;
    input.line_energies = strong_lines( nuc, age, 1.0E-4, 50.0 );
    input.sum_ranges.push_back( { 100.0, 3000.0, 1.0E-9 } );
    const RelActCalcCascade::SummingModel model( { input }, drf );
    const RelActCalcCascade::SummingModel::Source &source = model.sources()[0];
    BOOST_REQUIRE_MESSAGE( !source.sum_peaks.empty(), name );

    const RelActCalcCascade::DetectorEfficiencies &det = model.detector();
    const CascadeSummingCalc calc( { { nuc, age } }, { { 500.0, 1.0 } }, 1.0, drf, 0.0, false, g_data_dir );

    // The engine's rates are per decay of the parent at age zero; the kernel's per decay now.
    const double per_current_decay = exp( nuc->decayConstant()*age );
    const double window_tol = 0.3;

    for( const double dist_cm : { 0.0, 5.0 } )
    {
      const double g0 = model.front_solid_angle( dist_cm*PhysicalUnits::cm );
      const RelActCalcCascade::EvalResult<double> result = evaluate_model<double>( model, g0, nullopt, {} );

      // The kernel's own point efficiencies, where it has the energy; else the far-field scaling.
      vector<double> fep_abs, tot_abs;
      model.point_efficiencies( g0, fep_abs, tot_abs );
      const vector<double> &model_energies = model.energies();

      const function<double(double)> fep_fcn = [&]( double e ) -> double {
        const int index = CascadeSummingCalc::nearestEnergyIndex( model_energies, e );
        if( (index >= 0) && (fabs( model_energies[index] - e ) < 1.0E-6) )
          return fep_abs[index];
        return g0*RelActCalcCascade::detector_efficiencies( drf, { e } ).fep[0];
      };
      const function<double(double)> tot_fcn = [&]( double e ) -> double {
        const int index = CascadeSummingCalc::nearestEnergyIndex( model_energies, e );
        if( (index >= 0) && (fabs( model_energies[index] - e ) < 1.0E-6) )
          return tot_abs[index];
        return g0*RelActCalcCascade::detector_efficiencies( drf, { e } ).total[0];
      };

      vector<ceelo::PeakWindow> windows;
      vector<size_t> window_peak;
      for( size_t p = 0; p < source.sum_peaks.size(); ++p )
      {
        if( source.sum_peaks[p].line_index >= 0 )
          continue;  //summing-in to a line: the engine reports that with the line, not as a pure sum

        // Only peaks whose window holds their pairs and no other peak's: no neighbour within reach.
        const bool crowded = ((p > 0) && ((source.sum_peaks[p].energy - source.sum_peaks[p-1].energy) < 3.0*window_tol))
                             || (((p + 1) < source.sum_peaks.size())
                                 && ((source.sum_peaks[p+1].energy - source.sum_peaks[p].energy) < 3.0*window_tol));
        if( crowded )
          continue;

        ceelo::PeakWindow w;
        w.energy_keV = source.sum_peaks[p].energy;
        w.tolerance_keV = window_tol;
        windows.push_back( w );
        window_peak.push_back( p );
      }
      if( windows.empty() )
        continue;
      // The engine, without sums of three photons (which the kernel does not model), and with them.
      struct Provider final : public ceelo::EfficiencyProviderT<double>
      {
        const function<double(double)> *fep_f = nullptr, *tot_f = nullptr;
        double fep( double e ) const override { return (*fep_f)( e ); }
        double total( double e ) const override { return (*tot_f)( e ); }
      };
      Provider prov;
      prov.fep_f = &fep_fcn;
      prov.tot_f = &tot_fcn;
      ceelo::AnalyticCascadeOptions opts;
      opts.report_pure_sums = true;
      opts.enumerate_triples = false;
      const shared_ptr<const ceelo::AnalyticCascadeModel> engine_model = calc.cascadeModel( nuc, age );
      const vector<ceelo::AnalyticPeakResultT<double>> ref
                              = ceelo::compute_cascade_analytic( *engine_model, windows, prov, opts );
      opts.enumerate_triples = true;
      const vector<ceelo::AnalyticPeakResultT<double>> ref_triples
                              = ceelo::compute_cascade_analytic( *engine_model, windows, prov, opts );

      double largest_max_rate = 0.0;
      for( const RelActCalcCascade::SumPeak &peak : source.sum_peaks )
        largest_max_rate = max( largest_max_rate, peak.max_rate_per_decay );

      double worst = 0.0, strongest = 0.0, worst_triples = 0.0;
      size_t nchecked = 0, nestimated = 0;
      for( size_t w = 0; w < windows.size(); ++w )
      {
        if( ref[w].found || !(ref[w].pure_sum_rate > 0.0) )
          continue;
        const size_t p = window_peak[w];
        const RelActCalcCascade::SumPeak &peak = source.sum_peaks[p];

        double rate = 0.0, rate_swapped = 0.0;
        for( size_t q = 0; q < peak.pairs.size(); ++q )
        {
          const RelActCalcCascade::SumPair &pair = peak.pairs[q];
          rate += pair.weight * fep_abs[pair.energy_a] * result.sum_pairs[0][p][q].b_given_a;
          rate_swapped += pair.weight * fep_abs[pair.energy_b] * result.sum_pairs[0][p][q].a_given_b;
        }
        // For a point source the two anchorings are the same number.
        BOOST_CHECK_CLOSE( rate, rate_swapped, 1.0E-6 );

        // Only compare peaks whose pairs the window holds exactly (a neighbour's pairs can straddle a window).
        if( peak.max_rate_per_decay < 1.0E-3*largest_max_rate )
          continue;
        ++nchecked;
        strongest = max( strongest, rate );
        const double engine_rate = per_current_decay*ref[w].pure_sum_rate;
        const double rel = rate/engine_rate - 1.0;
        worst = max( worst, fabs( rel ) );
        worst_triples = max( worst_triples, ref_triples[w].pure_sum_rate/ref[w].pure_sum_rate - 1.0 );

        // Exact for pairs of level-scheme transitions, but for the weakest pairs of a peak, which are left out
        //  (up to a percent of its rate); a peak with other pairs (their survival estimated) is allowed a few
        //  percent.
        bool estimated = false;
        for( const RelActCalcCascade::SumPair &pair : peak.pairs )
          estimated = estimated || (pair.program_pair < 0);
        nestimated += estimated;
        const double estimated_tol = (dist_cm == 0.0) ? 0.2 : 0.05;
        const bool agrees = estimated ? (fabs( rel ) < estimated_tol) : ((rel < 2.0E-3) && (rel > -1.2E-2));
        BOOST_CHECK_MESSAGE( agrees, name << " sum at " << peak.energy
                             << " keV, " << dist_cm << " cm: kernel " << rate << ", engine " << engine_rate
                             << (estimated ? " (has estimated-survival pairs)" : "") );
      }//for( loop over windows )
      BOOST_TEST_MESSAGE( name << " at " << dist_cm << " cm: " << nchecked << " sum peaks checked (" << nestimated
                          << " with estimated-survival pairs), worst relative difference from the engine " << worst
                          << "; sums of three photons would add up to " << worst_triples
                          << "; strongest rate per decay " << strongest );
      if( string(name) != "Co60" )  //its one sum, 1173 + 1332, lands on its own (very weak) 2505.7 keV line
        BOOST_CHECK_GT( nchecked, 0u );
    }//for( loop over distances )
  }//for( loop over nuclides )
}//sum_peaks_match_engine


namespace
{
/** The uranium HPGe test spectrum (unshielded, 6000 s), its detector response, and the Physical Model
 isotopics setup embedded in it. */
struct UraniumProblem
{
  shared_ptr<const SpecUtils::Measurement> foreground, background;
  shared_ptr<const DetectorPeakResponse> drf;
  RelActCalcAuto::Options options;
};

UraniumProblem load_uranium_problem()
{
  set_data_dir();

  const string spec_path = SpecUtils::append_path( SpecUtils::append_path( g_test_data_dir, "RelActAutoReport" ),
                                                   "U235_Unshielded_6000.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ), "Test spectrum not at '" << spec_path << "'" );

  const shared_ptr<SpecMeas> meas = make_shared<SpecMeas>();
  BOOST_REQUIRE_MESSAGE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ),
                         "Failed to load '" << spec_path << "'" );

  UraniumProblem problem;
  for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
  {
    if( !m )
      continue;
    if( m->source_type() == SpecUtils::SourceType::Background )
      problem.background = m;
    else if( !problem.foreground )
      problem.foreground = m;
  }
  BOOST_REQUIRE( problem.foreground );

  problem.drf = meas->detector();
  BOOST_REQUIRE_MESSAGE( problem.drf, "Test N42 has no embedded detector response." );

  const unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
  BOOST_REQUIRE_MESSAGE( state && !state->options.rel_eff_curves.empty(), "Test N42 has no isotopics setup." );
  problem.options = state->options;
  BOOST_REQUIRE( problem.options.rel_eff_curves[0].rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel );

  // Keep the solves comparable: one model, no simplification search.
  problem.options.auto_simplify_model = false;

  return problem;
}//load_uranium_problem()


RelActCalcAuto::RelActAutoSolution solve_problem( const UraniumProblem &problem, const RelActCalcAuto::Options &options )
{
  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, problem.foreground, problem.background,
                                 problem.drf, {}, PeakFitUtils::coarse_det_type( problem.foreground, nullptr ),
                                 nullptr ) );
  BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                         "Solve failed: status=" << static_cast<int>(sol.m_status)
                         << ", error='" << sol.m_error_message << "'" );
  return sol;
}//solve_problem(...)


double mass_fraction( const RelActCalcAuto::RelActAutoSolution &sol, const string &nuclide_name )
{
  const pair<double,optional<double>> frac = sol.mass_enrichment_fraction( nuclide( nuclide_name ), 0 );
  return frac.first;
}
}//namespace


// The option survives a trip through XML, and a curve without it is written exactly as before.
BOOST_AUTO_TEST_CASE( solver_options_xml_round_trip )
{
  const UraniumProblem problem = load_uranium_problem();

  auto round_trip = []( const RelActCalcAuto::RelEffCurveInput &curve, string &xml ){
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> * const parent = doc.allocate_node( rapidxml::node_element, "Parent" );
    doc.append_node( parent );
    rapidxml::xml_node<char> * const node = curve.toXml( parent );
    BOOST_REQUIRE( node );
    xml.clear();
    rapidxml::print( std::back_inserter( xml ), doc, 0 );

    RelActCalcAuto::RelEffCurveInput read;
    BOOST_REQUIRE_NO_THROW( read.fromXml( node ) );
    return read;
  };

  string xml;
  RelActCalcAuto::RelEffCurveInput curve = problem.options.rel_eff_curves[0];
  BOOST_REQUIRE( !curve.cascade_summing.enabled );
  RelActCalcAuto::RelEffCurveInput read = round_trip( curve, xml );
  BOOST_CHECK( read == curve );
  BOOST_CHECK( !read.cascade_summing.enabled );
  BOOST_CHECK_MESSAGE( xml.find( "CascadeSumming" ) == string::npos, "option off, but written to XML" );

  curve.cascade_summing.enabled = true;
  curve.cascade_summing.fit_distance = false;
  curve.cascade_summing.distance = 3.5*PhysicalUnits::cm;
  curve.cascade_summing.source_radius = 1.25*PhysicalUnits::cm;
  read = round_trip( curve, xml );
  BOOST_CHECK( read == curve );
  BOOST_CHECK( read.cascade_summing.enabled && !read.cascade_summing.fit_distance );
  BOOST_CHECK_CLOSE( read.cascade_summing.distance, curve.cascade_summing.distance, 1.0E-5 );
  BOOST_CHECK_CLOSE( read.cascade_summing.source_radius, curve.cascade_summing.source_radius, 1.0E-5 );
  BOOST_CHECK( xml.find( "CascadeSumming" ) != string::npos );

  curve.cascade_summing = RelActCalcAuto::RelEffCurveInput::CascadeSummingInput();
  curve.cascade_summing.enabled = true;
  read = round_trip( curve, xml );
  BOOST_CHECK( read == curve );
  BOOST_CHECK( read.cascade_summing.fit_distance );
  BOOST_CHECK_EQUAL( read.cascade_summing.source_radius, 0.0 );

  // Only a Physical Model curve can have the correction.
  RelActCalcAuto::Options options = problem.options;
  options.rel_eff_curves[0].cascade_summing.enabled = true;
  options.rel_eff_curves[0].rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  options.rel_eff_curves[0].rel_eff_eqn_order = 3;
  options.rel_eff_curves[0].phys_model_self_atten.reset();
  options.rel_eff_curves[0].phys_model_external_atten.clear();
  const RelActCalcAuto::RelActAutoSolution bad = RelActCalcAuto::solve( options, problem.foreground,
                  problem.background, problem.drf, {}, PeakFitUtils::coarse_det_type( problem.foreground, nullptr ),
                  nullptr );
  BOOST_CHECK( bad.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success );
}//solver_options_xml_round_trip


// With the source far away the correction does nothing; near, it lowers the lines that have coincident
//  partners, leaves the ones that do not, and moves the answer.
BOOST_AUTO_TEST_CASE( solver_fixed_distance )
{
  const UraniumProblem problem = load_uranium_problem();
  typedef RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::CascadeSummingInfo CascadeInfo;

  const RelActCalcAuto::RelActAutoSolution off = solve_problem( problem, problem.options );
  BOOST_REQUIRE( (off.m_phys_model_results.size() == 1) && off.m_phys_model_results[0] );
  BOOST_CHECK( !off.m_phys_model_results[0]->cascade_summing );
  const double u235_off = mass_fraction( off, "U235" );

  // 5 m: the detector fills ~1e-5 of the directions.
  RelActCalcAuto::Options options = problem.options;
  options.rel_eff_curves[0].cascade_summing.enabled = true;
  options.rel_eff_curves[0].cascade_summing.fit_distance = false;
  options.rel_eff_curves[0].cascade_summing.distance = 5.0*PhysicalUnits::meter;

  const RelActCalcAuto::RelActAutoSolution far_sol = solve_problem( problem, options );
  BOOST_REQUIRE( far_sol.m_phys_model_results[0] && far_sol.m_phys_model_results[0]->cascade_summing );
  {
    const CascadeInfo &info = *far_sol.m_phys_model_results[0]->cascade_summing;
    BOOST_CHECK( !info.distance_was_fit );
    BOOST_CHECK_CLOSE( info.distance, 5.0*PhysicalUnits::meter, 0.1 );
    BOOST_CHECK_LT( info.solid_angle, 1.0E-4 );
    for( const CascadeInfo::LineFactor &line : info.lines )
      BOOST_CHECK_SMALL( line.factor - 1.0, 1.0E-4 );
    // The model is the same to a part in 1e5; the answer only to the half percent this problem's two nearby
    //  minima are apart (chi2 123.1 and 123.7), since which of them a solve ends in turns on its last digits.
    BOOST_CHECK_CLOSE( mass_fraction( far_sol, "U235" ), u235_off, 1.0 );
    BOOST_CHECK_CLOSE( far_sol.m_chi2, off.m_chi2, 1.0 );
    BOOST_CHECK_EQUAL( far_sol.m_final_parameters.size(), off.m_final_parameters.size() + 1 );
  }

  // 1 cm from the detector's face.
  options.rel_eff_curves[0].cascade_summing.distance = 1.0*PhysicalUnits::cm;
  const RelActCalcAuto::RelActAutoSolution near_sol = solve_problem( problem, options );
  BOOST_REQUIRE( near_sol.m_phys_model_results[0] && near_sol.m_phys_model_results[0]->cascade_summing );
  {
    const CascadeInfo &info = *near_sol.m_phys_model_results[0]->cascade_summing;
    BOOST_CHECK( !info.distance_was_fit );
    BOOST_CHECK_CLOSE( info.distance, 1.0*PhysicalUnits::cm, 0.1 );
    BOOST_CHECK_EQUAL( info.distance_lower, info.distance );
    BOOST_CHECK_EQUAL( info.distance_upper, info.distance );
    BOOST_CHECK_GT( info.solid_angle, 0.05 );
    BOOST_CHECK_MESSAGE( !info.total_efficiency_from_response,
                         "the library response has no total efficiency; the generic one should be used" );

    // The factor of the source's line nearest `energy` (the U-238 chain has lines 0.03 keV apart).
    auto factor = [&info]( const string &nuc, const double energy ) -> double {
      double answer = 1.0, nearest = 0.02;   //lines the correction does not change are not listed
      for( const CascadeInfo::LineFactor &line : info.lines )
      {
        if( (RelActCalcAuto::to_name( line.source ) == nuc) && (fabs( line.energy - energy ) < nearest) )
        {
          nearest = fabs( line.energy - energy );
          answer = line.factor;
        }
      }
      return answer;
    };

    // U-235's 185.7 keV and Pa-234m's 766.4 keV have little in coincidence; 202.1 and 258.3 keV a lot.  (This
    //  fit puts 2 cm of iron in front of the source, whose scattered photons add to every line's loss.)
    BOOST_CHECK_SMALL( factor( "U235", 185.71 ) - 1.0, 5.0E-3 );
    BOOST_CHECK_SMALL( factor( "U238", 766.37 ) - 1.0, 2.0E-2 );
    BOOST_CHECK_LT( factor( "U235", 202.11 ), 0.99 );
    BOOST_CHECK_GT( factor( "U235", 202.11 ), 0.90 );
    BOOST_CHECK_LT( factor( "U238", 258.26 ), 0.98 );
    BOOST_CHECK_GT( factor( "U238", 258.26 ), 0.80 );
    BOOST_CHECK_LT( 1.0 - factor( "U238", 766.37 ), 0.25*(1.0 - factor( "U238", 258.26 )) );
    BOOST_CHECK_LT( 1.0 - factor( "U235", 185.71 ), 0.25*(1.0 - factor( "U235", 202.11 )) );

    bool have_sum = false;
    for( const CascadeInfo::SumPeak &peak : info.sum_peaks )
      have_sum = have_sum || (peak.counts > 0.0);
    BOOST_CHECK_MESSAGE( have_sum, "no coincidence-sum peaks at 1 cm" );

    bool warned = false;
    for( const string &warning : near_sol.m_warnings )
      warned = warned || (warning.find( "Cascade summing correction" ) != string::npos);
    BOOST_CHECK_MESSAGE( warned, "no warning that a generic total efficiency was used" );

    BOOST_CHECK_MESSAGE( fabs( mass_fraction( near_sol, "U235" )/u235_off - 1.0 ) > 1.0E-3,
                         "U-235 fraction did not move: " << mass_fraction( near_sol, "U235" ) << " vs " << u235_off );
  }
}//solver_fixed_distance


// A fitted distance: the solve converges, and reports a distance with a range that brackets it.
BOOST_AUTO_TEST_CASE( solver_fit_distance )
{
  const UraniumProblem problem = load_uranium_problem();
  typedef RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::CascadeSummingInfo CascadeInfo;

  RelActCalcAuto::Options options = problem.options;
  options.rel_eff_curves[0].cascade_summing.enabled = true;
  BOOST_REQUIRE( options.rel_eff_curves[0].cascade_summing.fit_distance );

  const RelActCalcAuto::RelActAutoSolution sol = solve_problem( problem, options );
  BOOST_REQUIRE( sol.m_phys_model_results[0] && sol.m_phys_model_results[0]->cascade_summing );
  const CascadeInfo &info = *sol.m_phys_model_results[0]->cascade_summing;
  BOOST_CHECK( info.distance_was_fit );
  BOOST_CHECK_GE( info.distance, 0.0 );
  BOOST_CHECK_LE( info.distance_lower, info.distance );
  BOOST_CHECK_GE( info.distance_upper, info.distance );
  BOOST_CHECK_GE( info.solid_angle, 0.0 );

  bool named = false;
  for( const string &name : sol.m_parameter_names )
    named = named || (name == "CascadeSolidAngle_0");
  BOOST_CHECK_MESSAGE( named, "no CascadeSolidAngle_0 parameter" );

  const double u235 = mass_fraction( sol, "U235" );
  BOOST_CHECK( (u235 > 0.0) && (u235 < 1.0) );
}//solver_fit_distance


// The random-sum (pile-up) model: the option survives XML, adds its two parameters, reports what it fit, and
//  with nothing for it to fit leaves the answer alone.
BOOST_AUTO_TEST_CASE( solver_random_sums )
{
  const UraniumProblem problem = load_uranium_problem();

  RelActCalcAuto::Options options = problem.options;
  BOOST_REQUIRE( !options.random_sum_peaks );

  auto round_trip = []( const RelActCalcAuto::Options &opts, string &xml ){
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> * const parent = doc.allocate_node( rapidxml::node_element, "Parent" );
    doc.append_node( parent );
    rapidxml::xml_node<char> * const node = opts.toXml( parent );
    BOOST_REQUIRE( node );
    xml.clear();
    rapidxml::print( std::back_inserter( xml ), doc, 0 );
    RelActCalcAuto::Options read;
    BOOST_REQUIRE_NO_THROW( read.fromXml( node ) );
    return read;
  };

  string xml;
  BOOST_CHECK( !round_trip( options, xml ).random_sum_peaks );
  BOOST_CHECK_MESSAGE( xml.find( "RandomSumPeaks" ) == string::npos, "option off, but written to XML" );
  options.random_sum_peaks = true;
  BOOST_CHECK( round_trip( options, xml ).random_sum_peaks );
  BOOST_CHECK( xml.find( "RandomSumPeaks" ) != string::npos );

  const RelActCalcAuto::RelActAutoSolution off = solve_problem( problem, problem.options );
  BOOST_CHECK( !off.m_random_sums );

  const RelActCalcAuto::RelActAutoSolution sol = solve_problem( problem, options );
  BOOST_REQUIRE( sol.m_random_sums );
  BOOST_CHECK_EQUAL( sol.m_final_parameters.size(), off.m_final_parameters.size() + 2 );

  size_t num_named = 0;
  for( const string &name : sol.m_parameter_names )
    num_named += ((name == "RandomSumTime") || (name == "RandomSumXray"));
  BOOST_CHECK_EQUAL( num_named, 2u );

  const RelActCalcAuto::RelActAutoSolution::RandomSumInfo &info = *sol.m_random_sums;
  BOOST_CHECK_GE( info.pileup_time, 0.0 );
  BOOST_CHECK_LE( info.pileup_time, 20.0E-6*PhysicalUnits::second );
  BOOST_CHECK_GE( info.xray_fraction, 0.0 );
  BOOST_CHECK_LE( info.xray_fraction, 0.5 );
  BOOST_REQUIRE( info.xray_element );
  BOOST_CHECK_EQUAL( info.xray_element->symbol, string("U") );   //the nuclides are uranium's
  for( const RelActCalcAuto::RelActAutoSolution::RandomSumInfo::Peak &peak : info.peaks )
  {
    BOOST_CHECK_CLOSE( peak.energy, peak.energy_a + peak.energy_b, 1.0E-6 );
    BOOST_CHECK_GE( peak.counts, 0.0 );
  }

  // A 6000 s unshielded spectrum at a modest rate: whatever was fit, the enrichment barely moves.
  const double u235_off = mass_fraction( off, "U235" ), u235 = mass_fraction( sol, "U235" );
  BOOST_CHECK_MESSAGE( fabs( u235/u235_off - 1.0 ) < 0.05, "U-235 fraction " << u235 << " vs " << u235_off );
}//solver_random_sums


namespace
{
/** A condition on the U-235 mass fraction of the first curve. */
RelActCalcAuto::Options::OptionRule::Condition u235_fraction_condition( const optional<double> lower,
                                                                        const optional<double> upper )
{
  RelActCalcAuto::Options::OptionRule::Condition condition;
  condition.quantity.kind = RelActCalcAuto::Options::ProfileTarget::Kind::MassFraction;
  condition.quantity.source = nuclide( "U235" );
  condition.quantity.rel_eff_curve_index = 0;
  condition.lower = lower;
  condition.upper = upper;
  return condition;
}


bool has_warning_with( const RelActCalcAuto::RelActAutoSolution &sol, const string &text )
{
  for( const string &warning : sol.m_warnings )
  {
    if( warning.find( text ) != string::npos )
      return true;
  }
  return false;
}
}//namespace


// What a rule changes, what makes one unusable, and that rules survive a trip through XML.
BOOST_AUTO_TEST_CASE( option_rule_changes_and_xml )
{
  typedef RelActCalcAuto::Options::OptionRule Rule;
  const UraniumProblem problem = load_uranium_problem();
  const RelActCalcAuto::Options &base = problem.options;
  BOOST_REQUIRE_GE( base.rois.size(), 3u );
  BOOST_REQUIRE( base.option_rules.empty() && !base.random_sum_peaks );

  // The fixture's ROIs are not in energy order.
  const vector<RelActCalcAuto::RoiRange>::const_iterator highest = std::max_element( begin(base.rois), end(base.rois),
    []( const RelActCalcAuto::RoiRange &lhs, const RelActCalcAuto::RoiRange &rhs ){
      return lhs.upper_energy < rhs.upper_energy;
  } );
  const double highest_roi_lower = highest->lower_energy, highest_roi_upper = highest->upper_energy;
  const RelActCalcAuto::RoiRange kept = (highest == begin(base.rois)) ? base.rois[1] : base.rois[0];

  Rule rule;
  rule.name = "Test rule";
  rule.conditions.push_back( u235_fraction_condition( 0.2, 0.8 ) );
  rule.remove_rois.push_back( {highest_roi_lower - 0.1, highest_roi_upper + 0.1} );

  RelActCalcAuto::RoiRange above;         //above every ROI: added
  above.lower_energy = highest_roi_upper + 100.0;
  above.upper_energy = highest_roi_upper + 110.0;
  above.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  RelActCalcAuto::RoiRange overlapping = kept;   //overlaps a ROI that stays: not added
  overlapping.lower_energy += 0.25*(overlapping.upper_energy - overlapping.lower_energy);
  overlapping.upper_energy += 1.0;
  rule.add_rois = { above, overlapping };

  rule.random_sum_peaks = true;
  Rule::CurveChange change;
  change.rel_eff_curve_index = 0;
  change.phys_model_corr_fcn = RelActCalc::PhysModelCorrFcn::None;
  change.cascade_summing = RelActCalcAuto::RelEffCurveInput::CascadeSummingInput();
  change.cascade_summing->enabled = true;
  rule.curve_changes.push_back( change );
  Rule::CurveChange no_such_curve = change;
  no_such_curve.rel_eff_curve_index = 7;    //ignored
  rule.curve_changes.push_back( no_such_curve );

  RelActCalcAuto::Options changed = base;
  rule.apply( changed );
  BOOST_CHECK_EQUAL( changed.rois.size(), base.rois.size() );     //one removed, one added
  BOOST_CHECK( std::find( begin(changed.rois), end(changed.rois), above ) != end(changed.rois) );
  BOOST_CHECK( std::find( begin(changed.rois), end(changed.rois), kept ) != end(changed.rois) );
  BOOST_CHECK( std::find( begin(changed.rois), end(changed.rois), overlapping ) == end(changed.rois) );
  for( const RelActCalcAuto::RoiRange &roi : changed.rois )
    BOOST_CHECK( (roi.lower_energy != highest_roi_lower) || (roi.upper_energy != highest_roi_upper) );

  // Added ROIs go in energy order, where the ROIs are in it.
  RelActCalcAuto::Options sorted = base;
  std::sort( begin(sorted.rois), end(sorted.rois), []( const RelActCalcAuto::RoiRange &lhs,
                                                      const RelActCalcAuto::RoiRange &rhs ){
    return lhs.lower_energy < rhs.lower_energy;
  } );
  rule.apply( sorted );
  BOOST_CHECK( sorted.rois.back() == above );
  for( size_t i = 1; i < sorted.rois.size(); ++i )
    BOOST_CHECK_LE( sorted.rois[i-1].lower_energy, sorted.rois[i].lower_energy );
  BOOST_CHECK( changed.random_sum_peaks );
  BOOST_CHECK( changed.rel_eff_curves[0].phys_model_corr.corr_fcn == RelActCalc::PhysModelCorrFcn::None );
  BOOST_CHECK( changed.rel_eff_curves[0].cascade_summing.enabled );
  BOOST_CHECK( changed.option_rules.empty() );

  // Applying a rule again changes nothing more.
  RelActCalcAuto::Options twice = changed;
  rule.apply( twice );
  BOOST_CHECK( twice == changed );

  // A rule that changes nothing leaves the options as they were.
  Rule nothing;
  nothing.name = "Nothing";
  nothing.conditions = rule.conditions;
  RelActCalcAuto::Options same = base;
  nothing.apply( same );
  BOOST_CHECK( same == base );

  // Usable, and what is not.
  RelActCalcAuto::Options options = base;
  options.option_rules = { rule };
  BOOST_CHECK_MESSAGE( options.why_not_usable().empty(), options.why_not_usable() );

  const auto unusable = [&base]( const Rule &bad ) -> bool {
    RelActCalcAuto::Options opts = base;
    opts.option_rules = { bad };
    return !opts.why_not_usable().empty();
  };

  Rule bad = rule;
  bad.name.clear();
  BOOST_CHECK( unusable( bad ) );

  bad = rule;
  bad.conditions.clear();
  BOOST_CHECK( unusable( bad ) );

  bad = rule;
  bad.conditions[0].lower.reset();
  bad.conditions[0].upper.reset();
  BOOST_CHECK( unusable( bad ) );

  bad = rule;
  bad.conditions[0].lower = 0.9;   //not below the upper limit
  BOOST_CHECK( unusable( bad ) );

  bad = rule;
  bad.conditions[0].quantity.kind = RelActCalcAuto::Options::ProfileTarget::Kind::RelativeActivity;
  BOOST_CHECK( unusable( bad ) );

  bad = rule;
  bad.conditions[0].quantity.kind = RelActCalcAuto::Options::ProfileTarget::Kind::ActivityRatio;
  bad.conditions[0].quantity.denominator = nuclide( "U238" );
  bad.conditions[0].quantity.denominator_curve_index = 1;   //another curve's
  BOOST_CHECK( unusable( bad ) );
  bad.conditions[0].quantity.denominator_curve_index = 0;
  BOOST_CHECK( !unusable( bad ) );

  bad = rule;
  bad.remove_rois = { {0.0, 1.0E6} };   //leaves only the added ROI
  BOOST_CHECK( !unusable( bad ) );
  bad.add_rois.clear();                 //leaves no ROI
  BOOST_CHECK( unusable( bad ) );

  bad = rule;
  bad.remove_rois = { {500.0, 400.0} };
  BOOST_CHECK( unusable( bad ) );

  // XML
  auto round_trip = []( const RelActCalcAuto::Options &opts, string &xml ){
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> * const parent = doc.allocate_node( rapidxml::node_element, "Parent" );
    doc.append_node( parent );
    rapidxml::xml_node<char> * const node = opts.toXml( parent );
    BOOST_REQUIRE( node );
    xml.clear();
    rapidxml::print( std::back_inserter( xml ), doc, 0 );
    RelActCalcAuto::Options read;
    BOOST_REQUIRE_NO_THROW( read.fromXml( node ) );
    return read;
  };

  string xml;
  BOOST_CHECK( round_trip( base, xml ).option_rules.empty() );
  BOOST_CHECK_MESSAGE( xml.find( "OptionRule" ) == string::npos, "no rules, but written to XML" );
  BOOST_CHECK_MESSAGE( xml.find( "<Options version=\"8\"" ) == string::npos,
                       "options using nothing new must keep their older version" );

  Rule second = nothing;
  second.name = "Second rule";
  second.conditions = { u235_fraction_condition( 0.8, nullopt ), u235_fraction_condition( nullopt, 0.95 ) };
  second.conditions[1].quantity.kind = RelActCalcAuto::Options::ProfileTarget::Kind::ActivityRatio;
  second.conditions[1].quantity.denominator = nuclide( "U238" );
  second.random_sum_peaks = false;
  Rule::CurveChange summing_off;
  summing_off.cascade_summing = RelActCalcAuto::RelEffCurveInput::CascadeSummingInput();   //disabled
  second.curve_changes = { summing_off };

  options.option_rules = { rule, second };
  const RelActCalcAuto::Options read = round_trip( options, xml );
  BOOST_CHECK( xml.find( "<Options version=\"8\"" ) != string::npos );
  BOOST_REQUIRE_EQUAL( read.option_rules.size(), 2u );
  BOOST_CHECK( read.option_rules[0] == rule );
  BOOST_CHECK( read.option_rules[1] == second );
  BOOST_CHECK( read == options );
  BOOST_REQUIRE_EQUAL( read.option_rules[1].curve_changes.size(), 1u );
  BOOST_REQUIRE( read.option_rules[1].curve_changes[0].cascade_summing.has_value() );
  BOOST_CHECK( !read.option_rules[1].curve_changes[0].cascade_summing->enabled );
  BOOST_CHECK( !read.option_rules[1].curve_changes[0].phys_model_corr_fcn.has_value() );
  BOOST_REQUIRE( read.option_rules[1].random_sum_peaks.has_value() );
  BOOST_CHECK( !read.option_rules[1].random_sum_peaks.value() );

  // Conditions
  const Rule::Condition range = u235_fraction_condition( 0.2, 0.8 );
  BOOST_CHECK( range.holds( 0.2 ) && range.holds( 0.5 ) && !range.holds( 0.8 ) && !range.holds( 0.1 ) );
  BOOST_CHECK( !range.holds( std::numeric_limits<double>::quiet_NaN() ) );
  BOOST_CHECK( u235_fraction_condition( nullopt, 0.8 ).holds( 0.0 ) );
  BOOST_CHECK( u235_fraction_condition( 0.8, nullopt ).holds( 1.0 ) );
}//option_rule_changes_and_xml


// The solver applies the first rule whose conditions hold, and its answer is then that of the changed
//  options; with no rule holding, it is that of the options as given.
BOOST_AUTO_TEST_CASE( solver_option_rules )
{
  typedef RelActCalcAuto::Options::OptionRule Rule;
  const UraniumProblem problem = load_uranium_problem();

  const RelActCalcAuto::RelActAutoSolution plain = solve_problem( problem, problem.options );
  const double u235_plain = mass_fraction( plain, "U235" );
  BOOST_REQUIRE( (u235_plain > 0.0) && (u235_plain < 1.0) );
  BOOST_CHECK( !plain.m_applied_option_rule.has_value() );

  // A change that is easy to see in the result: no correction function.
  BOOST_REQUIRE( problem.options.rel_eff_curves[0].phys_model_corr.corr_fcn != RelActCalc::PhysModelCorrFcn::None );
  Rule rule;
  rule.name = "No correction function";
  Rule::CurveChange change;
  change.phys_model_corr_fcn = RelActCalc::PhysModelCorrFcn::None;
  rule.curve_changes.push_back( change );

  RelActCalcAuto::Options changed_options = problem.options;
  rule.apply( changed_options );
  const RelActCalcAuto::RelActAutoSolution changed = solve_problem( problem, changed_options );
  const double u235_changed = mass_fraction( changed, "U235" );
  BOOST_REQUIRE_MESSAGE( u235_changed != u235_plain, "the test's change must move the answer" );

  // 1. The conditions hold whatever the fit says: the changed options are what is fit.
  rule.conditions = { u235_fraction_condition( 0.0, nullopt ) };
  RelActCalcAuto::Options options = problem.options;
  options.option_rules = { rule };
  RelActCalcAuto::RelActAutoSolution sol = solve_problem( problem, options );
  BOOST_REQUIRE( sol.m_applied_option_rule.has_value() );
  BOOST_CHECK_EQUAL( sol.m_applied_option_rule->rule_index, 0u );
  BOOST_CHECK_EQUAL( sol.m_applied_option_rule->name, rule.name );
  BOOST_REQUIRE_EQUAL( sol.m_applied_option_rule->condition_values.size(), 1u );
  BOOST_CHECK_CLOSE( sol.m_applied_option_rule->condition_values[0], mass_fraction( sol, "U235" ), 1.0E-9 );
  BOOST_CHECK_EQUAL( mass_fraction( sol, "U235" ), u235_changed );
  BOOST_CHECK( sol.m_options.option_rules.empty() );
  BOOST_CHECK( sol.m_options.rel_eff_curves[0].phys_model_corr.corr_fcn == RelActCalc::PhysModelCorrFcn::None );
  BOOST_CHECK( has_warning_with( sol, "'No correction function' option rule was applied" ) );

  // 2. They never hold: the options are fit as given.
  rule.conditions = { u235_fraction_condition( nullopt, 0.0 ) };
  options.option_rules = { rule };
  sol = solve_problem( problem, options );
  BOOST_CHECK( !sol.m_applied_option_rule.has_value() );
  BOOST_CHECK_EQUAL( mass_fraction( sol, "U235" ), u235_plain );
  BOOST_CHECK( sol.m_options.option_rules.empty() );
  BOOST_CHECK( sol.m_options.rel_eff_curves[0].phys_model_corr.corr_fcn
               == problem.options.rel_eff_curves[0].phys_model_corr.corr_fcn );
  BOOST_CHECK( !has_warning_with( sol, "option rule" ) );

  // 3. The first rule that holds is the one applied.
  Rule never = rule;
  never.name = "Never";
  Rule always = rule;
  always.name = "Always";
  always.conditions = { u235_fraction_condition( 0.0, nullopt ) };
  Rule also = always;
  also.name = "Also always";
  also.curve_changes.clear();
  options.option_rules = { never, always, also };
  sol = solve_problem( problem, options );
  BOOST_REQUIRE( sol.m_applied_option_rule.has_value() );
  BOOST_CHECK_EQUAL( sol.m_applied_option_rule->rule_index, 1u );
  BOOST_CHECK_EQUAL( sol.m_applied_option_rule->name, string("Always") );
  BOOST_CHECK_EQUAL( mass_fraction( sol, "U235" ), u235_changed );

  // 4. A condition on something the problem does not have never holds, and the solution says so.
  Rule absent = always;
  absent.name = "Plutonium";
  absent.conditions[0].quantity.source = nuclide( "Pu239" );
  options.option_rules = { absent };
  sol = solve_problem( problem, options );
  BOOST_CHECK( !sol.m_applied_option_rule.has_value() );
  BOOST_CHECK_EQUAL( mass_fraction( sol, "U235" ), u235_plain );
  BOOST_CHECK( has_warning_with( sol, "'Plutonium' was not used" ) );

  // 5. A rule that holds for the answer of the options as given, but not for the answer its own changes
  //  give: whichever of the two the initial estimate picks, the fit's answer picks the other, the fit is
  //  redone once, and the solution says the spectrum is at the boundary.
  const double between = 0.5*(u235_plain + u235_changed);
  rule.name = "Boundary";
  rule.conditions = { (u235_plain > u235_changed) ? u235_fraction_condition( between, nullopt )
                                                  : u235_fraction_condition( nullopt, between ) };
  BOOST_REQUIRE( rule.conditions[0].holds( u235_plain ) && !rule.conditions[0].holds( u235_changed ) );
  options.option_rules = { rule };
  sol = solve_problem( problem, options );
  BOOST_CHECK( has_warning_with( sol, "so the fit was redone that way" ) );
  BOOST_CHECK( has_warning_with( sol, "at the boundary between option rules" ) );
  const double u235_boundary = mass_fraction( sol, "U235" );
  if( sol.m_applied_option_rule.has_value() )
    BOOST_CHECK_EQUAL( u235_boundary, u235_changed );
  else
    BOOST_CHECK_EQUAL( u235_boundary, u235_plain );
}//solver_option_rules


namespace
{
/** A library (flat disk) response: 5 cm diameter. */
shared_ptr<DetectorPeakResponse> library_detective()
{
  set_data_dir();
  ifstream tsv( SpecUtils::append_path( g_data_dir, "common_drfs.tsv" ).c_str(), ios::binary );
  vector<string> credits;
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  DetectorPeakResponse::parseMultipleRelEffDrfCsv( tsv, credits, drfs );
  for( const shared_ptr<DetectorPeakResponse> &d : drfs )
  {
    if( d && (d->name().find( "MicroDetective_LANL_025cm" ) != string::npos) )
      return d;
  }
  BOOST_REQUIRE_MESSAGE( false, "no MicroDetective_LANL_025cm in common_drfs.tsv" );
  return nullptr;
}//library_detective()


/** `<g^2>/<g>` of a disk source face of `source_radius`, over the on-axis solid angle, for a detector disk of
 `det_radius` a `height` away: by brute force over the detector's area and the source's radius. */
double brute_force_wide_source_ratio( const double det_radius, const double height, const double source_radius )
{
  const int num_r = 3000, num_phi = 720, num_rings = 16;
  auto solid_angle = [&]( const double offset ) -> double {
    double sum = 0.0;
    for( int i = 0; i < num_r; ++i )
    {
      const double r = (i + 0.5)*det_radius/num_r;
      for( int j = 0; j < num_phi; ++j )
      {
        const double phi = (j + 0.5)*2.0*PhysicalUnits::pi/num_phi;
        const double dist2 = height*height + r*r + offset*offset - 2.0*r*offset*cos( phi );
        sum += height*r/(dist2*sqrt( dist2 ));
      }
    }
    return sum*(det_radius/num_r)*(2.0*PhysicalUnits::pi/num_phi)/(4.0*PhysicalUnits::pi);
  };

  double sum = 0.0, sum_sq = 0.0;
  for( int ring = 0; ring < num_rings; ++ring )
  {
    const double g = solid_angle( source_radius*sqrt( (ring + 0.5)/num_rings ) );
    sum += g;
    sum_sq += g*g;
  }
  return (sum_sq/sum)/solid_angle( 0.0 );
}//brute_force_wide_source_ratio(...)
}//namespace


// A wide source right at the detector: the solid angle keeps rising smoothly as the source comes in, to the
//  value in the detector's own plane (where every point over the detector sees half of all directions).
BOOST_AUTO_TEST_CASE( wide_source_near_the_detector_plane )
{
  const shared_ptr<DetectorPeakResponse> flat = library_detective();
  // The geometry is what is tested: put the crystal's face where distances are measured from.
  flat->setDetectorSetback( 0.0 );
  BOOST_REQUIRE_EQUAL( flat->detectorSetback(), 0.0 );

  RelActCalcCascade::SourceInput input;
  input.nuclide = nuclide( "Co60" );
  input.age = PhysicalUnits::year;
  input.line_energies = { 1173.228, 1332.492 };
  const RelActCalcCascade::SummingModel model( { input }, flat );

  const double det_radius = 0.5*model.detector().crystal_diameter;
  const double source_radius = 3.5*PhysicalUnits::cm;
  BOOST_REQUIRE_GT( source_radius, det_radius );

  // Against a brute-force integral, where that still converges.
  for( const double height_mm : { 30.0, 3.0, 1.0, 0.5 } )
  {
    const double height = height_mm*PhysicalUnits::mm;
    const double ratio = model.front_solid_angle( height, source_radius )/model.front_solid_angle( height );
    const double expected = brute_force_wide_source_ratio( det_radius, height, source_radius );
    BOOST_CHECK_MESSAGE( fabs( ratio/expected - 1.0 ) < 5.0E-3, "at " << height_mm << " mm: " << ratio
                         << " of the on-axis solid angle, brute force " << expected );
  }

  // Monotone, and continuous into the plane - where the part of the face over the detector sees half of all
  //  directions and the rest nothing, so `<g^2>/<g>` is a half: the same as a point on the axis.
  double previous = model.front_solid_angle( 0.0, source_radius );
  BOOST_CHECK_CLOSE( previous, model.front_solid_angle( 0.0 ), 1.0E-6 );
  BOOST_CHECK_CLOSE( model.front_solid_angle( 0.0 ), 0.5, 1.0E-6 );
  for( const double height_mm : { 1.0E-4, 1.0E-3, 0.01, 0.03, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0 } )
  {
    const double wide = model.front_solid_angle( height_mm*PhysicalUnits::mm, source_radius );
    BOOST_CHECK_MESSAGE( wide <= previous*(1.0 + 1.0E-9), "not falling with distance at " << height_mm << " mm: "
                         << wide << " after " << previous );
    previous = wide;
  }
  BOOST_CHECK_CLOSE( model.front_solid_angle( 1.0E-4*PhysicalUnits::mm, source_radius ), 0.5, 0.5 );

  // So the distance can be recovered from the solid angle all the way in.
  for( const double height_mm : { 0.05, 0.5, 2.0 } )
  {
    const double g = model.front_solid_angle( height_mm*PhysicalUnits::mm, source_radius );
    BOOST_CHECK_SMALL( model.front_distance( g, source_radius )/PhysicalUnits::mm - height_mm, 2.0E-3 );
  }
}//wide_source_near_the_detector_plane


// Co-60's two gammas, and their sum peak, in a self-attenuating slab: with only unscattered photons counted
//  there is a closed form for each, which tells the two members of the pair apart.
BOOST_AUTO_TEST_CASE( sum_pair_in_a_slab_closed_form )
{
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  RelActCalcCascade::SourceInput input;
  input.nuclide = nuclide( "Co60" );
  input.age = PhysicalUnits::year;
  input.line_energies = { 1173.228, 1332.492 };
  input.sum_ranges.push_back( { 2495.0, 2515.0, 0.0 } );

  const RelActCalcCascade::SummingModel model( { input }, drf, false );
  BOOST_REQUIRE_EQUAL( model.sources()[0].sum_peaks.size(), 1u );
  const RelActCalcCascade::SumPeak &peak = model.sources()[0].sum_peaks[0];
  BOOST_CHECK_CLOSE( peak.energy, 1173.228 + 1332.492, 1.0E-3 );
  BOOST_CHECK_EQUAL( peak.line_index, -1 );
  BOOST_REQUIRE_EQUAL( peak.pairs.size(), 1u );
  const RelActCalcCascade::SumPair &pair = peak.pairs[0];
  BOOST_CHECK_LT( model.energies()[pair.energy_a], model.energies()[pair.energy_b] );

  const shared_ptr<const Material> lead = MaterialDB::instance()->material( "Pb (lead)" );
  BOOST_REQUIRE( lead );
  const shared_ptr<const vector<double>> mu = model.material_mu( lead );
  auto slab = []( const double x ) -> double { return (x < 1.0E-8) ? (1.0 - 0.5*x) : (-expm1( -x )/x); };

  const double g0 = model.front_solid_angle( 2.0*PhysicalUnits::cm );
  vector<double> fep, total;
  model.point_efficiencies( g0, fep, total );

  for( const double ad_g_cm2 : { 0.0, 5.0, 40.0 } )
  {
    RelActCalc::PhysModelShield<double> self;
    self.material = lead;
    self.areal_density = ad_g_cm2*PhysicalUnits::g_per_cm2;
    const RelActCalcCascade::EvalResult<double> result
        = evaluate_model<double>( model, g0, optional<RelActCalc::PhysModelShield<double>>( self ), {} );
    const RelActCalcCascade::SumPairValue<double> &value = result.sum_pairs[0][0][0];

    // Nothing else is emitted with the pair, so its coincident survival is one.
    const double xa = (*mu)[pair.energy_a]*self.areal_density, xb = (*mu)[pair.energy_b]*self.areal_density;
    BOOST_CHECK_CLOSE( value.b_given_a, fep[pair.energy_b]*slab( xa + xb )/slab( xa ), 0.05 );
    BOOST_CHECK_CLOSE( value.a_given_b, fep[pair.energy_a]*slab( xa + xb )/slab( xb ), 0.05 );
    if( ad_g_cm2 > 1.0 )   //the two differ: swapping the members would show
      BOOST_CHECK_GT( fabs( slab( xa )/slab( xb ) - 1.0 ), 0.01 );
  }//for( loop over areal densities )
}//sum_pair_in_a_slab_closed_form


// For a geometry-aware response the kernel's absolute efficiencies are the response's own, at distances on
//  and between those it tabulates the response at.
BOOST_AUTO_TEST_CASE( point_efficiencies_follow_the_response )
{
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  BOOST_REQUIRE( drf->ceeloResponse() );

  RelActCalcCascade::SourceInput input;
  input.nuclide = nuclide( "U235" );
  input.age = 20.0*PhysicalUnits::year;
  input.line_energies = strong_lines( input.nuclide, input.age, 5.0E-3, 100.0 );
  BOOST_REQUIRE_GE( input.line_energies.size(), 4u );
  const RelActCalcCascade::SummingModel model( { input }, drf );
  BOOST_REQUIRE( model.detector().total_source == RelActCalcCascade::TotalEffSource::DetectorResponse );

  const vector<double> &energies = model.energies();
  double worst_fep = 0.0, worst_tot = 0.0;
  for( const double dist_cm : { 0.5, 1.0, 1.25, 2.6, 4.5, 7.0, 30.0 } )
  {
    const double dist = dist_cm*PhysicalUnits::cm;
    vector<double> fep, total;
    model.point_efficiencies( model.front_solid_angle( dist ), fep, total );
    for( const double line : input.line_energies )
    {
      const int index = GammaInteractionCalc::CascadeSummingCalc::nearestEnergyIndex( energies, line );
      BOOST_REQUIRE( (index >= 0) && (fabs( energies[index] - line ) < 1.0E-6) );
      const double drf_fep = drf->fepEfficiencyEval( static_cast<float>(line), 0.0, 0.0, dist ).value;
      const double drf_tot = drf->totalEfficiencyEval( static_cast<float>(line), 0.0, 0.0, dist ).value;
      worst_fep = max( worst_fep, fabs( fep[index]/drf_fep - 1.0 ) );
      worst_tot = max( worst_tot, fabs( total[index]/drf_tot - 1.0 ) );
      BOOST_CHECK_MESSAGE( fabs( fep[index]/drf_fep - 1.0 ) < 0.01, line << " keV at " << dist_cm
                           << " cm: full-energy " << fep[index] << ", the response's " << drf_fep );
      BOOST_CHECK_MESSAGE( fabs( total[index]/drf_tot - 1.0 ) < 0.01, line << " keV at " << dist_cm
                           << " cm: total " << total[index] << ", the response's " << drf_tot );
    }
  }//for( loop over distances )
  BOOST_TEST_MESSAGE( "kernel against the response: worst full-energy " << worst_fep << ", total " << worst_tot );

  // A solid angle that is not a number comes back as such, rather than reading outside the tables.
  const RelActCalcCascade::EvalResult<double> nan_result
      = evaluate_model<double>( model, numeric_limits<double>::quiet_NaN(), nullopt, {} );
  BOOST_REQUIRE_EQUAL( nan_result.summing_out.size(), 1u );

  // Whether the response's total efficiency is used does not depend on the energies asked for.
  const RelActCalcCascade::TotalEffSource one = RelActCalcCascade::detector_efficiencies( drf, { 60.0 } ).total_source;
  const RelActCalcCascade::TotalEffSource many
      = RelActCalcCascade::detector_efficiencies( drf, { 20.0, 60.0, 661.7, 2614.5, 9000.0 } ).total_source;
  BOOST_CHECK( one == many );
}//point_efficiencies_follow_the_response


// Rates are per decay of the source now, whatever ages were asked for before; and a source with nothing
//  emitted together has nothing for a geometry to change.
BOOST_AUTO_TEST_CASE( ages_and_inert_sources )
{
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();

  // Cs-134 has no gamma-emitting progeny: per decay of the Cs-134 there now, its lines do not depend on age.
  const SandiaDecay::Nuclide * const cs134 = nuclide( "Cs134" );
  vector<double> emissions;
  for( const double age_years : { 5.0, 5.09, 5.0 } )
  {
    RelActCalcCascade::SourceInput input;
    input.nuclide = cs134;
    input.age = age_years*PhysicalUnits::year;
    input.line_energies = strong_lines( cs134, input.age, 0.5, 500.0 );   //604.7 and 795.9 keV
    BOOST_REQUIRE_EQUAL( input.line_energies.size(), 2u );
    input.sum_ranges.push_back( { 1395.0, 1405.0, 0.0 } );
    const RelActCalcCascade::SummingModel model( { input }, drf );
    BOOST_REQUIRE( model.sources()[0].line_found[0] && model.sources()[0].line_found[1] );
    BOOST_REQUIRE( !model.sources()[0].sum_peaks.empty() );
    emissions.push_back( model.sources()[0].line_emission[0] );
    emissions.push_back( model.sources()[0].sum_peaks[0].pairs[0].weight );
    BOOST_CHECK( model.has_effect() );
  }
  BOOST_CHECK_CLOSE( emissions[2], emissions[0], 1.0E-6 );
  BOOST_CHECK_CLOSE( emissions[3], emissions[1], 1.0E-6 );
  BOOST_CHECK_EQUAL( emissions[4], emissions[0] );
  BOOST_CHECK_EQUAL( emissions[5], emissions[1] );
  BOOST_CHECK_CLOSE( emissions[0], 0.976, 1.0 );   //604.7 keV per decay

  // Cs-137: one gamma (and the Ba K x-rays of its conversion, which the detector response does not reach).
  RelActCalcCascade::SourceInput cs137;
  cs137.nuclide = nuclide( "Cs137" );
  cs137.age = 5.0*PhysicalUnits::year;
  cs137.line_energies = { 661.657 };
  cs137.sum_ranges.push_back( { 600.0, 1400.0, 0.0 } );
  const RelActCalcCascade::SummingModel inert( { cs137 }, drf );
  BOOST_TEST_MESSAGE( "Cs-137: " << inert.sources()[0].sum_peaks.size() << " sum peaks, has an effect: "
                      << inert.has_effect() );
  if( inert.detector().lower_energy > 37.5 )   //the Ba K x-rays are below what the response is characterized for
    BOOST_CHECK( !inert.has_effect() );
}//ages_and_inert_sources


// Behind a steel plate the summing losses fall by less than the plate's transmission would give, because
//  photons scattered in the plate still reach the detector: against Monte Carlo of Eu-152 at 1 cm, bare and
//  behind 3 mm of iron.  (The Monte Carlo is of another, larger, detector, so only the ratio of the loss
//  behind iron to the bare loss is compared.)
BOOST_AUTO_TEST_CASE( shield_scatter_against_monte_carlo )
{
  const shared_ptr<const DetectorPeakResponse> drf = example_hpge_mc();
  const string path = SpecUtils::append_path( SpecUtils::append_path( g_test_data_dir, "rel_act_cascade" ),
                                              "eu152_cs137_mc_truth.csv" );
  ifstream csv( path.c_str() );
  BOOST_REQUIRE_MESSAGE( csv.good(), "no " << path );

  // energy -> the fraction of the line lost to summing, bare and behind iron.
  map<double,pair<double,double>> truth;
  string line;
  getline( csv, line );  //header
  while( getline( csv, line ) )
  {
    vector<string> fields;
    SpecUtils::split( fields, line, "," );
    if( (fields.size() < 7) || (fields[0] != "1") || (fields[2] != "Eu152") )
      continue;
    const double energy = stod( fields[3] ), loss = 1.0 - stod( fields[6] )/stod( fields[5] );
    if( fields[1] == "none" )
      truth[energy].first = loss;
    else if( fields[1] == "fe3" )
      truth[energy].second = loss;
  }
  BOOST_REQUIRE_GE( truth.size(), 10u );

  const SandiaDecay::Nuclide * const eu152 = nuclide( "Eu152" );
  RelActCalcCascade::SourceInput input;
  input.nuclide = eu152;
  input.age = PhysicalUnits::year;
  input.line_energies = strong_lines( eu152, input.age, 1.0E-3, 100.0 );

  const shared_ptr<const Material> iron = MaterialDB::instance()->material( "Fe (iron)" );
  BOOST_REQUIRE( iron );
  RelActCalc::PhysModelShield<double> plate;
  plate.material = iron;
  plate.areal_density = 0.3*PhysicalUnits::cm*iron->density;

  // The Monte Carlo's efficiencies are net of summing-in, which the factors compared here do not have: only
  //  the lines nothing sums into.
  const double clean_lines[] = { 121.78, 244.70, 344.28, 411.12, 778.90, 867.38 };

  size_t num_checked = 0;
  double mean_with = 0.0, mean_without = 0.0, mean_truth = 0.0;
  for( const bool scatter : { true, false } )
  {
    const RelActCalcCascade::SummingModel model( { input }, drf, scatter );
    BOOST_REQUIRE_EQUAL( model.models_shield_scatter(), scatter );
    const double g0 = model.front_solid_angle( 1.0*PhysicalUnits::cm );
    const RelActCalcCascade::EvalResult<double> bare = evaluate_model<double>( model, g0, nullopt, {} );
    const RelActCalcCascade::EvalResult<double> shielded = evaluate_model<double>( model, g0, nullopt, { plate } );

    for( size_t i = 0; i < input.line_energies.size(); ++i )
    {
      const double energy = input.line_energies[i];
      const double loss_bare = 1.0 - bare.summing_out[0][i], loss_fe = 1.0 - shielded.summing_out[0][i];

      const bool clean = std::any_of( begin(clean_lines), end(clean_lines), [energy]( const double clean_energy ){
        return fabs( clean_energy - energy ) < 0.05;
      } );
      const map<double,pair<double,double>>::const_iterator pos = truth.lower_bound( energy - 0.05 );
      if( !clean || (pos == end(truth)) || (fabs( pos->first - energy ) > 0.05) )
        continue;
      BOOST_REQUIRE_GT( pos->second.first, 0.03 );
      BOOST_REQUIRE_GT( loss_bare, 0.02 );

      const double ratio = loss_fe/loss_bare, truth_ratio = pos->second.second/pos->second.first;
      if( scatter )
      {
        num_checked += 1;
        mean_with += ratio;
        mean_truth += truth_ratio;
        BOOST_TEST_MESSAGE( energy << " keV: loss behind 3 mm of iron over the bare loss " << ratio
                            << ", Monte Carlo " << truth_ratio );
      }else
      {
        mean_without += ratio;
      }
    }//for( loop over lines )
  }//for( with and without scattered photons )

  BOOST_REQUIRE_EQUAL( num_checked, sizeof(clean_lines)/sizeof(clean_lines[0]) );
  mean_with /= num_checked;
  mean_without /= num_checked;
  mean_truth /= num_checked;
  BOOST_TEST_MESSAGE( num_checked << " lines: mean ratio " << mean_with << " with scattered photons, "
                      << mean_without << " without, Monte Carlo " << mean_truth );
  // Without the scattered photons the losses behind the plate are too small by more than the scatter
  //  model (good to a fifth or so of what it adds) leaves them too large.
  BOOST_CHECK_GT( mean_with, mean_without + 0.03 );
  BOOST_CHECK_LT( mean_without, mean_truth );
  BOOST_CHECK_LT( fabs( mean_with - mean_truth ), fabs( mean_without - mean_truth ) );
  BOOST_CHECK_SMALL( mean_with - mean_truth, 0.10 );
}//shield_scatter_against_monte_carlo


namespace
{
/** A copy of `data`, scaled by `count_scale`, with the counts in each of the solution's ROIs replaced by what
 its fit predicts there - continuum and every peak, without noise - so that a refit has a known answer.  The
 scaled spectrum is that of a source `count_scale` times as active: random-sum peaks, which go as the square
 of the rate, are scaled by its square.
 */
shared_ptr<SpecUtils::Measurement> model_spectrum( const RelActCalcAuto::RelActAutoSolution &sol,
                                                  const shared_ptr<const SpecUtils::Measurement> &data,
                                                  const double count_scale )
{
  map<const PeakContinuum *,vector<PeakDef>> rois;
  for( const PeakDef &peak : sol.m_fit_peaks_in_spectrums_cal )
  {
    PeakDef copy = peak;
    if( SpecUtils::starts_with( peak.userLabel(), RelActCalcAuto::Options::sm_random_sum_label_prefix ) )
      copy.setAmplitude( count_scale*peak.amplitude() );
    rois[peak.continuum().get()].push_back( copy );
  }

  const vector<float> &energies = *data->channel_energies();
  const shared_ptr<vector<float>> counts = make_shared<vector<float>>( *data->gamma_counts() );
  for( float &count : *counts )
    count *= static_cast<float>( count_scale );

  for( const auto &roi : rois )
  {
    const PeakContinuum &continuum = *roi.first;
    const size_t first = data->find_gamma_channel( static_cast<float>( continuum.lowerEnergy() ) );
    const size_t last = data->find_gamma_channel( static_cast<float>( continuum.upperEnergy() ) );
    BOOST_REQUIRE( (first <= last) && ((last + 1) < energies.size()) );
    const size_t num_channels = last + 1 - first;

    vector<const PeakDef *> peaks;
    for( const PeakDef &peak : roi.second )
      peaks.push_back( &peak );
    vector<double> model( num_channels, 0.0 );
    continuum.offset_integral( &(energies[first]), model.data(), num_channels, data, peaks.data(), peaks.size() );
    for( const PeakDef &peak : roi.second )
      peak.gauss_integral( &(energies[first]), model.data(), num_channels );
    for( size_t i = 0; i < num_channels; ++i )
      (*counts)[first + i] = static_cast<float>( count_scale*std::max( 0.0, model[i] ) );
  }//for( loop over ROIs )

  const shared_ptr<SpecUtils::Measurement> answer = make_shared<SpecUtils::Measurement>( *data );
  answer->set_gamma_counts( counts, data->live_time(), data->real_time() );
  return answer;
}//model_spectrum(...)


RelActCalcAuto::RelActAutoSolution solve_spectrum( const UraniumProblem &problem,
                                                   const shared_ptr<const SpecUtils::Measurement> &spectrum,
                                                   const RelActCalcAuto::Options &options )
{
  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, spectrum, nullptr, problem.drf, {},
                                 PeakFitUtils::coarse_det_type( spectrum, nullptr ), nullptr ) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "Solve failed: status=" << static_cast<int>(sol.m_status)
                         << ", error='" << sol.m_error_message << "'" );
  return sol;
}//solve_spectrum(...)


RelActCalcAuto::RoiRange fixed_roi( const double lower, const double upper )
{
  RelActCalcAuto::RoiRange roi;
  roi.lower_energy = lower;
  roi.upper_energy = upper;
  roi.continuum_type = PeakContinuum::OffsetType::Linear;
  roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  return roi;
}
}//namespace


// A source 1 cm from the detector at a high count: U-235's 345.9 and 387.8 keV peaks are then mostly
//  coincidence sums.  The spectrum is the fit's own model (without noise, so the answer is known), refit with
//  the distance free, starting from far away.  This guards the start of such fits (RelActCalcAuto_PreFit.hpp):
//  when the initial estimate took the sum peaks for lines, this fit ended with the source on the detector's
//  face.  One configuration only; the isotopics benchmark covers the variety.
BOOST_AUTO_TEST_CASE( solver_close_source_with_sum_peaks )
{
  const UraniumProblem problem = load_uranium_problem();
  typedef RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::CascadeSummingInfo CascadeInfo;

  RelActCalcAuto::Options options = problem.options;
  // Continua that do not depend on the data (a step continuum's shape does), so the model is what it predicts.
  for( RelActCalcAuto::RoiRange &roi : options.rois )
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
  // ROIs for the sums of 185.7 keV with the U K x-rays, with itself, and with 202.1 keV (387.8 keV).
  options.rois.push_back( fixed_roi( 276.0, 288.0 ) );
  options.rois.push_back( fixed_roi( 366.5, 376.5 ) );
  options.rois.push_back( fixed_roi( 383.0, 393.0 ) );
  options.random_sum_peaks = true;
  options.rel_eff_curves[0].cascade_summing.enabled = true;
  options.rel_eff_curves[0].cascade_summing.fit_distance = false;
  options.rel_eff_curves[0].cascade_summing.distance = 1.0*PhysicalUnits::cm;

  const RelActCalcAuto::RelActAutoSolution truth = solve_spectrum( problem, problem.foreground, options );
  const double u235_truth = mass_fraction( truth, "U235" );
  BOOST_REQUIRE( (u235_truth > 0.0) && (u235_truth < 1.0) );

  const shared_ptr<SpecUtils::Measurement> spectrum = model_spectrum( truth, problem.foreground, 30.0 );
  options.rel_eff_curves[0].cascade_summing.fit_distance = true;
  const RelActCalcAuto::RelActAutoSolution sol = solve_spectrum( problem, spectrum, options );

  BOOST_REQUIRE( (sol.m_phys_model_results.size() == 1) && sol.m_phys_model_results[0]
                 && sol.m_phys_model_results[0]->cascade_summing );
  const CascadeInfo &info = *sol.m_phys_model_results[0]->cascade_summing;
  BOOST_CHECK( info.distance_was_fit );
  BOOST_CHECK_MESSAGE( (info.distance > 0.7*PhysicalUnits::cm) && (info.distance < 1.5*PhysicalUnits::cm),
                       "fitted distance " << info.distance/PhysicalUnits::cm << " cm, from 1 cm" );
  // Without noise the true answer has a chi2 of about zero; the fit stops in the flat valley short of it.
  BOOST_CHECK_LT( sol.m_chi2/sol.m_dof, 0.5 );
  BOOST_CHECK_CLOSE( mass_fraction( sol, "U235" ), u235_truth, 3.0 );
}//solver_close_source_with_sum_peaks
