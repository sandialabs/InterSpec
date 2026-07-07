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
#include <random>
#include <string>
#include <iostream>

#include <Wt/Utils.h>
#include <Wt/WApplication.h>
#include <Wt/Test/WTestEnvironment.h>

#ifdef _WIN32
// For some reason, we need to include the following includes, before unit_test.hpp,
//  or we get a bunch of errors relating to winsock.k being included multiple times,
//  or something
#include "winsock2.h"
#include "Windows.h"
#endif

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ShieldingSourceFitCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

//Roots Minuit2 includes
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"


#include "ceres/jet.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"


using namespace std;
using namespace boost::unit_test;


std::string g_test_file_dir;

// We need to set the static data directory, so the code knows where
//  like sandia.decay.xml is located.
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
    cout << "Arg " << i << ": '" << argv[i] << "'" << endl;
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
    
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  const string required_data_file = "findCharacteristics/202204_example_problem_1.n42";
  if( g_test_file_dir.empty() )
  {
    for( const auto &d : { "test_data", "../test_data", "../../test_data", "/Users/wcjohns/rad_ana/InterSpec/target/testing/test_data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, required_data_file) ) )
      {
        g_test_file_dir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }
  
  const string sandia_deacay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_deacay_file ), "sandia.decay.xml not at '" << sandia_deacay_file << "'" );
  
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  
  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );
  
  const string required_data_path = SpecUtils::append_path(g_test_file_dir, required_data_file);
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( required_data_path ), "'" << required_data_file << "' not at '" << required_data_path << "'" );
}//void set_data_dir()


// Helper class to manage InterSpec instance for tests
class InterSpecTestFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_update_lock;
  InterSpec *m_interspec;

  InterSpecTestFixture()
  : m_env( nullptr ),
    m_app( nullptr ),
    m_interspec( nullptr )
  {
    set_data_dir();

    string wt_app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), "..");
    wt_app_root = SpecUtils::lexically_normalize_path(wt_app_root);

    // Create a test environment
    const std::string applicationPath = "";
    const std::string configurationFile = "";
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::EntryPointType::Application ) );
    m_env->setAppRoot( wt_app_root );

    // Create the app
    m_app = new InterSpecApp( *m_env );

    m_update_lock.reset( new Wt::WApplication::UpdateLock(m_app) );

    // Get the InterSpec viewer instance
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );
  }

  ~InterSpecTestFixture()
  {
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
  }
};


BOOST_AUTO_TEST_CASE( ShieldingInfoUri )
{
  using namespace ShieldingSourceFitCalc;
  
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const std::shared_ptr<const Material> iron = matdb->material( "Fe (iron)" );
  BOOST_REQUIRE( iron );

  const std::shared_ptr<const Material> uranium = matdb->material( "U (uranium)" );
  BOOST_REQUIRE( uranium );
  
  {// Begin test Trace source serialization
    TraceSourceInfo trace;
    trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
    trace.m_fitActivity = true;
    BOOST_REQUIRE( trace.m_nuclide = db->nuclide( "U238" ) );
    trace.m_activity = 2.0*PhysicalUnits::microCi / PhysicalUnits::cm2;
    trace.m_relaxationDistance = 0;
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( trace.serialize( &doc ) );
      
      TraceSourceInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( TraceSourceInfo::equalEnough( trace, from_xml ) );
    }
    
    trace.m_type = GammaInteractionCalc::TraceActivityType::ExponentialDistribution;
    trace.m_relaxationDistance = 5*PhysicalUnits::cm;
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( trace.serialize( &doc ) );
      
      TraceSourceInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( TraceSourceInfo::equalEnough( trace, from_xml ) );
    }
  }// End test Trace source serialization
  
  
  {// Begin test Trace source serialization, with empty src
    TraceSourceInfo trace;
    trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
    trace.m_fitActivity = false;
    trace.m_nuclide = nullptr;
    trace.m_activity = 0.0f;
    trace.m_relaxationDistance = 0.0f;
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( trace.serialize( &doc ) );
      
      TraceSourceInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( TraceSourceInfo::equalEnough( trace, from_xml ) );
    }
  }// End test Trace source serialization, with empty src
  

  
  
  {// Begin test ShieldingInfo encode/decode from URI/XML
    // URI currently does not encode self-attenuating source, trace-source information, or "truth" values.
    ShieldingInfo info;
    
    info.m_geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
    info.m_isGenericMaterial = false;
    info.m_forFitting = true;
    info.m_material = make_shared<Material>( *iron );
    info.m_dimensions[0] = 1.2*PhysicalUnits::meter;
    info.m_dimensions[1] = 1.1*PhysicalUnits::meter;
    info.m_dimensions[2] = 1.0E-6*PhysicalUnits::meter;
    info.m_fitDimensions[0] = false;
    info.m_fitDimensions[1] = true;
    info.m_fitDimensions[2] = false;
    
    
    {// Begin test URI - doesnt include truth values or trace/self-atten sources, so test before we set those
      const string uri = info.encodeStateToUrl();
      
      ShieldingInfo from_uri;
      BOOST_REQUIRE_NO_THROW( from_uri.handleAppUrl( uri ) );
      BOOST_CHECK_NO_THROW( ShieldingInfo::equalEnough( info, from_uri ) );
    }// End test URI
    
    
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    // Truth info only tested for XML serialization
    info.m_truthDimensions[0] = 1.11E-6f;
    info.m_truthDimensions[1] = boost::optional<float>();
    info.m_truthDimensions[2] = 9.99E-3f;
    
    info.m_truthDimensionsTolerances[0] = 1.11E-7f;
    info.m_truthDimensionsTolerances[1] = 100.0f;
    info.m_truthDimensionsTolerances[2] = boost::optional<float>();
#endif
    
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( info.serialize( &doc ) );
      
      ShieldingInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( ShieldingInfo::equalEnough( info, from_xml ) );
    }
    
    
    // Now modify `info` to include trace and self-attenuating sources
    info.m_material = make_shared<Material>( *uranium );
    
    
    // Self-atten source stuff
    {// Begin add self-atten source
      map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,bool>>> 
        &el_to_nuc_frac = info.m_nuclideFractions_;
      
      const SandiaDecay::Nuclide *nuc = db->nuclide( "U-238" );
      const SandiaDecay::Element *el = db->element( nuc->atomicNumber );
      BOOST_REQUIRE( nuc && el );
      el_to_nuc_frac[el].emplace_back( nuc, 0.8, true );
      
      nuc = db->nuclide( "U-235" );
      el = db->element( nuc->atomicNumber );
      BOOST_REQUIRE( nuc && el );
      el_to_nuc_frac[el].emplace_back( nuc, 0.15, true );
      
      nuc = db->nuclide( "U-232" );
      el = db->element( nuc->atomicNumber );
      BOOST_REQUIRE( nuc && el );
      el_to_nuc_frac[el].emplace_back( nuc, 2.5E-8, false );
    }// End add self-atten source
    
    // Trace-source stuff
    {// Begin add trace source source
      
      TraceSourceInfo trace;
      trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
      trace.m_fitActivity = false;
      trace.m_nuclide = db->nuclide( "Cs137" );
      trace.m_activity = 2.0*PhysicalUnits::microCi / PhysicalUnits::cm2;
      trace.m_relaxationDistance = 0.0f;
      info.m_traceSources.push_back( trace );
      
      
      trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerGram;
      trace.m_fitActivity = true;
      trace.m_nuclide = db->nuclide( "Ba133" );
      trace.m_activity = 5.0*PhysicalUnits::microCi / PhysicalUnits::gram;
      trace.m_relaxationDistance = 0.0f;
      info.m_traceSources.push_back( trace );
      
      trace.m_type = GammaInteractionCalc::TraceActivityType::ExponentialDistribution;
      trace.m_fitActivity = true;
      trace.m_nuclide = db->nuclide( "I131" );
      trace.m_activity = 1.0E5*PhysicalUnits::bq / PhysicalUnits::m2;
      trace.m_relaxationDistance = 1.23*PhysicalUnits::cm;
      info.m_traceSources.push_back( trace );
    }// End add trace source source
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( info.serialize( &doc ) );
      
      ShieldingInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( ShieldingInfo::equalEnough( info, from_xml ) );
    }
  }// End test ShieldingInfo encode/decode from URI
}//BOOST_AUTO_TEST_CASE( ShieldingInfoUri )


BOOST_AUTO_TEST_CASE( SourceFitDefSerialization )
{
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  
  ShieldingSourceFitCalc::SourceFitDef test;
  
  test.nuclide = db->nuclide( "U238" );
  test.activity = 1.1*PhysicalUnits::curie;
  test.fitActivity = true;
  test.age = 20*PhysicalUnits::year;
  test.fitAge = false;
  test.ageDefiningNuc = db->nuclide( "U235" );
  test.sourceType = ShieldingSourceFitCalc::ModelSourceType::Intrinsic;
  test.activityUncertainty = 0.1*test.activity;
  test.ageUncertainty.reset();
  
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  test.truthActivity = 1.0*PhysicalUnits::curie;
  test.truthActivityTolerance = 0.15*PhysicalUnits::curie;
  test.truthAge = 21.5*PhysicalUnits::year;
  test.truthAgeTolerance = 1.05*PhysicalUnits::year;
#endif
  
  rapidxml::xml_document<char> doc;
  BOOST_REQUIRE_NO_THROW( test.serialize( &doc ) );
  
  ShieldingSourceFitCalc::SourceFitDef from_xml;
  BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
  BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::SourceFitDef::equalEnough( test, from_xml ) );
}//BOOST_AUTO_TEST_CASE( SourceFitDefSerialization )



BOOST_AUTO_TEST_CASE( SrcFitOptionsSerialization )
{
  set_data_dir();
  
  for( size_t i = 0; i < 20; ++i )
  {
    ShieldingSourceFitCalc::ShieldingSourceFitOptions test;
    
    test.multiple_nucs_contribute_to_peaks = (rand() % 2);
    test.attenuate_for_air = (rand() % 2);
    test.account_for_decay_during_meas = (rand() % 2);
    test.multithread_self_atten = (rand() % 2);
    test.photopeak_cluster_sigma = ((10.0*rand()) / RAND_MAX);
    test.background_peak_subtract = (rand() % 2);
    test.same_age_isotopes = (rand() % 2);
    
    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( test.serialize( &doc ) );
    
    ShieldingSourceFitCalc::ShieldingSourceFitOptions from_xml;
    BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( &doc ) );
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( test, from_xml ) );
  }//for( size_t i = 0; i < 20; ++i )
}//BOOST_AUTO_TEST_CASE( SourceFitDefSerialization )


// A simple, fully programmatically defined case of fitting a source/shielding model
BOOST_AUTO_TEST_CASE( SimpleSourceFit )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  
  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 100*PhysicalUnits::second;
  const float peak_area = PhysicalUnits::microCi;
  
  const GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::Spherical;
  vector<ShieldingSourceFitCalc::ShieldingInfo> shieldings;
  vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  
  {
    ShieldingSourceFitCalc::SourceFitDef cs137_src;
    cs137_src.nuclide = db->nuclide( "Cs137" );
    cs137_src.activity = 2.0*peak_area;
    cs137_src.fitActivity = true;
    cs137_src.age = 180*PhysicalUnits::day;
    cs137_src.fitAge = false;
    cs137_src.ageDefiningNuc = nullptr;
    cs137_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
    
    src_definitions.push_back(cs137_src);
  }
  
  // Create a 100% efficient detector at our test distance.
  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance,
                                    5*PhysicalUnits::cm, PhysicalUnits::keV,
                                    0, 3000*PhysicalUnits::keV,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  
  auto foreground = make_shared<SpecUtils::Measurement>();
  
  auto spec = make_shared<vector<float>>();
  *spec = {0.0f, 1.0f, 5.0f, 2.0f, 3.5f};
  foreground->set_gamma_counts( spec, live_time, live_time );
  
  std::deque<std::shared_ptr<const PeakDef>> foreground_peaks;
  
  {// Begin create peak for 661 keV in Cs137
    auto peak = make_shared<PeakDef>();
    peak->setMean( 661 );
    peak->setSigma( 10.0 );
    peak->setPeakArea( peak_area );
    peak->setPeakAreaUncert( sqrt(peak_area) );
    
    const auto cs137 = db->nuclide( "Cs137" );
    BOOST_REQUIRE( cs137 );
    const auto ba137m = db->nuclide( "Ba137m" ); //The 661 keV actually comes from Ba137m
    BOOST_REQUIRE( ba137m );
    int radParticle = -1;
    const SandiaDecay::Transition *transition = nullptr;
    
    for( size_t i = 0; !transition && (i < ba137m->decaysToChildren.size()); ++i )
    {
      const SandiaDecay::Transition *trans = ba137m->decaysToChildren[i];
      for( size_t j = 0; !transition && (j < trans->products.size()); ++j )
      {
        if( (trans->products[j].type == SandiaDecay::GammaParticle)
           && fabs(trans->products[j].energy - 661.657) < 1.0)
        {
          transition = trans;
          radParticle = static_cast<int>(j);
        }
      }
    }//
    
    BOOST_REQUIRE( transition && (radParticle >= 0) );
    
    peak->setNuclearTransition( cs137, transition, radParticle, PeakDef::SourceGammaType::NormalGamma );
    peak->useForShieldingSourceFit( true );
    
    foreground_peaks.push_back( peak );
  }// End create peak for 661 keV in Cs137
  
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;
  
  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = geometry;
  chi_input.config.shieldings = shieldings;
  chi_input.config.sources = src_definitions;
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = nullptr;
  chi_input.foreground_peaks = foreground_peaks;
  chi_input.background_peaks = nullptr;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
                GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
 
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;
  
  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  
  
  auto progress_fcn = [progress](){
    // Probably wont get called, since its a simple fit
  };
  
  bool finished_fit_called = false;
  auto finished_fcn = [results, &finished_fit_called](){
    finished_fit_called = true;
  };
  
  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );
  
  BOOST_CHECK( finished_fit_called );
  
  BOOST_REQUIRE( results->fit_src_info.size() == 1 );
  
  const ShieldingSourceFitCalc::SourceFitDef &fit_src_info = results->fit_src_info[0];
  BOOST_REQUIRE( results->fit_src_info.size() == 1 );
  BOOST_CHECK( fit_src_info.nuclide == db->nuclide( "Cs137" ) );
  BOOST_CHECK( !fit_src_info.ageDefiningNuc );
  BOOST_CHECK( fit_src_info.sourceType == ShieldingSourceFitCalc::ModelSourceType::Point );
  
  const double activity = results->fit_src_info[0].activity;
  
  const double br = 0.8533;
  const double drf_eff = detector->efficiency( 661.657, distance );
  
  const double expected_activity = peak_area / (live_time * br * drf_eff);

  BOOST_CHECK_MESSAGE( fabs(activity - expected_activity) < 1.0E-4*expected_activity,
                       "Fit activity (" << PhysicalUnits::printToBestActivityUnits(activity,10)
                       << ") didnt match expected (" << PhysicalUnits::printToBestActivityUnits(expected_activity,10) << ")" );
}//BOOST_AUTO_TEST_CASE( SimpleSourceFit )


namespace
{
/** Finds a transition + product index for a gamma of the given energy, emitted
 anywhere in `parent`s decay chain (i.e., possibly by a descendant nuclide).
 */
bool find_gamma_transition( const SandiaDecay::Nuclide * const parent,
                            const double energy,
                            const SandiaDecay::Transition *&transition,
                            int &particle_index )
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
        if( (product.type == SandiaDecay::GammaParticle)
            && (fabs(product.energy - energy) < 0.25) )
        {
          transition = trans;
          particle_index = static_cast<int>( prod_index );
          return true;
        }
      }//for( loop over transition products )
    }//for( loop over transitions )
  }//for( loop over descendant nuclides )

  return false;
}//find_gamma_transition(...)


/** Creates a Gaussian peak assigned to the gamma of `parent`s decay chain at `energy`. */
std::shared_ptr<PeakDef> make_test_peak( const SandiaDecay::Nuclide * const parent,
                                         const double energy,
                                         const double sigma,
                                         const double area )
{
  auto peak = make_shared<PeakDef>();
  peak->setMean( energy );
  peak->setSigma( sigma );
  peak->setPeakArea( area );
  peak->setPeakAreaUncert( sqrt(area) );

  const SandiaDecay::Transition *transition = nullptr;
  int particle_index = -1;
  BOOST_REQUIRE_MESSAGE( find_gamma_transition( parent, energy, transition, particle_index ),
                         "No gamma at " << energy << " keV in " << parent->symbol << " decay chain" );

  peak->setNuclearTransition( parent, transition, particle_index, PeakDef::SourceGammaType::NormalGamma );
  peak->useForShieldingSourceFit( true );

  return peak;
}//make_test_peak(...)


/** Returns copies of `input.foreground_peaks`, with each peaks area set to the
 forward-model expectation for the configuration in `input` - i.e., the chi2
 minimum of the returned peaks is at `input`s parameter values (the "truth").
 */
deque<shared_ptr<const PeakDef>> peaks_with_model_expected_areas(
                const GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput &input )
{
  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( input );

  const vector<double> truth_pars = fcn_pars.second.Params();

  GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache mix_cache;
  const vector<GammaInteractionCalc::PeakResultPlotInfo> peak_infos
        = fcn_pars.first->energy_chi_contributions( truth_pars, {}, mix_cache, nullptr, nullptr );

  BOOST_REQUIRE_EQUAL( peak_infos.size(), input.foreground_peaks.size() );

  deque<shared_ptr<const PeakDef>> answer;

  for( const shared_ptr<const PeakDef> &orig : input.foreground_peaks )
  {
    const double gamma_energy = orig->gammaParticleEnergy();

    const GammaInteractionCalc::PeakResultPlotInfo *closest = nullptr;
    for( const GammaInteractionCalc::PeakResultPlotInfo &info : peak_infos )
    {
      if( !closest || (fabs(info.energy - gamma_energy) < fabs(closest->energy - gamma_energy)) )
        closest = &info;
    }

    BOOST_REQUIRE( closest );
    BOOST_REQUIRE_MESSAGE( fabs(closest->energy - gamma_energy) < 1.0,
                           "No model prediction near " << gamma_energy << " keV" );
    BOOST_REQUIRE_MESSAGE( (closest->observedOverExpected > 0.0) && !IsInf(closest->observedOverExpected),
                           "Model predicts no counts at " << gamma_energy << " keV" );

    const double expected_area = orig->peakArea() / closest->observedOverExpected;
    BOOST_REQUIRE( expected_area > 0.0 );

    auto peak = make_shared<PeakDef>( *orig );
    peak->setPeakArea( expected_area );
    peak->setPeakAreaUncert( sqrt(expected_area) );
    answer.push_back( peak );
  }//for( loop over input peaks )

  return answer;
}//peaks_with_model_expected_areas(...)
}//namespace


/** Minuit2 baseline: fit generic-shielding atomic number + areal density + source
 activity on synthetic peak areas generated from the forward model, so the true
 minimum is at known parameter values.

 This is the regression fixture the (future) Ceres-based driver gets compared
 against; the same fixture exercises the coarse atomic-number scan in fit_model.
 */
BOOST_AUTO_TEST_CASE( FitGenericShieldingANBaseline )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );

  const SandiaDecay::Nuclide * const ba133 = db->nuclide( "Ba133" );
  BOOST_REQUIRE( ba133 );

  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 600*PhysicalUnits::second;
  const double true_an = 26.0;
  const double true_ad = 10.0*PhysicalUnits::g/PhysicalUnits::cm2;
  const double true_activity = 10.0*PhysicalUnits::microCi;
  const double age = 5.0*PhysicalUnits::year;

  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance,
                                    5*PhysicalUnits::cm, PhysicalUnits::keV,
                                    0, 3000*PhysicalUnits::keV,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  // Ba133s main lines; spanning 81 to 384 keV gives leverage on atomic number
  //  (photoelectric component at low energy).
  const double placeholder_area = 1.0E4;
  deque<shared_ptr<const PeakDef>> truth_peaks;
  for( const double energy : { 80.9979, 276.3989, 302.8508, 356.0129, 383.8485 } )
    truth_peaks.push_back( make_test_peak( ba133, energy, 1.0, placeholder_area ) );

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;

  ShieldingSourceFitCalc::ShieldingInfo generic_shield;
  generic_shield.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  generic_shield.m_isGenericMaterial = true;
  generic_shield.m_forFitting = true;
  generic_shield.m_material = nullptr;
  generic_shield.m_dimensions[0] = true_an;
  generic_shield.m_dimensions[1] = true_ad;
  generic_shield.m_dimensions[2] = 0.0;
  generic_shield.m_fitDimensions[0] = false;
  generic_shield.m_fitDimensions[1] = false;
  generic_shield.m_fitDimensions[2] = false;

  ShieldingSourceFitCalc::SourceFitDef ba133_src;
  ba133_src.nuclide = ba133;
  ba133_src.activity = true_activity;
  ba133_src.fitActivity = false;
  ba133_src.age = age;
  ba133_src.fitAge = false;
  ba133_src.ageDefiningNuc = nullptr;
  ba133_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  chi_input.config.shieldings = { generic_shield };
  chi_input.config.sources = { ba133_src };
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = nullptr;
  chi_input.foreground_peaks = truth_peaks;
  chi_input.background_peaks = nullptr;

  // Set peak areas to the forward-model expectation at the truth values
  const deque<shared_ptr<const PeakDef>> expected_peaks = peaks_with_model_expected_areas( chi_input );

  // Now setup the actual fit, starting away from the truth values
  generic_shield.m_dimensions[0] = 40.0;
  generic_shield.m_dimensions[1] = 2.0*PhysicalUnits::g/PhysicalUnits::cm2;
  generic_shield.m_fitDimensions[0] = true;
  generic_shield.m_fitDimensions[1] = true;
  ba133_src.activity = 1.0*PhysicalUnits::microCi;
  ba133_src.fitActivity = true;

  chi_input.config.shieldings = { generic_shield };
  chi_input.config.sources = { ba133_src };
  chi_input.foreground_peaks = expected_peaks;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;

  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  auto progress_fcn = [](){};
  bool finished_fit_called = false;
  auto finished_fcn = [&finished_fit_called](){ finished_fit_called = true; };

  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );

  BOOST_CHECK( finished_fit_called );
  BOOST_CHECK( results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );
  BOOST_REQUIRE_EQUAL( results->final_shieldings.size(), size_t(1) );
  BOOST_REQUIRE_EQUAL( results->fit_src_info.size(), size_t(1) );

  // Note: for generic materials, FitShieldingInfo::m_dimensions[1] is in g/cm2
  //  (unlike the PhysicalUnits-valued input ShieldingInfo) - see fit_model().
  const double true_ad_g_cm2 = true_ad / (PhysicalUnits::g/PhysicalUnits::cm2);
  const double fit_an = results->final_shieldings[0].m_dimensions[0];
  const double fit_ad_g_cm2 = results->final_shieldings[0].m_dimensions[1];
  const double fit_activity = results->fit_src_info[0].activity;

  cout << "FitGenericShieldingANBaseline: AN=" << fit_an
       << " (truth " << true_an << "), AD=" << fit_ad_g_cm2
       << " g/cm2 (truth " << true_ad_g_cm2
       << "), activity=" << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
       << " (truth " << PhysicalUnits::printToBestActivityUnits(true_activity,6) << ")" << endl;

  BOOST_CHECK_MESSAGE( fabs(fit_an - true_an) < 2.5,
                       "Fit AN (" << fit_an << ") not close to truth (" << true_an << ")" );
  BOOST_CHECK_MESSAGE( fabs(fit_ad_g_cm2 - true_ad_g_cm2) < 0.05*true_ad_g_cm2,
                       "Fit AD (" << fit_ad_g_cm2 << " g/cm2) not close to truth ("
                       << true_ad_g_cm2 << " g/cm2)" );
  BOOST_CHECK_MESSAGE( fabs(fit_activity - true_activity) < 0.05*true_activity,
                       "Fit activity (" << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
                       << ") not close to truth ("
                       << PhysicalUnits::printToBestActivityUnits(true_activity,6) << ")" );

  // Uncertainties should be present and sane
  BOOST_CHECK( results->fit_src_info[0].activityUncertainty.has_value()
               && (results->fit_src_info[0].activityUncertainty.value() > 0.0) );
  BOOST_CHECK( results->final_shieldings[0].m_dimensionUncerts[1] > 0.0 );
}//BOOST_AUTO_TEST_CASE( FitGenericShieldingANBaseline )


/** A non-source shielding whose thickness STARTS at exactly 0, fit to data that requires a real
 (1 cm) thickness.  A shield's effect is pure attenuation (exp(-mu*chord), no volume normalization),
 and the chord is the Stage-1-continuous quantity, so the Jet derivative at thickness 0 is finite,
 nonzero, and correctly signed - Ceres should climb off zero to the truth.  (The variant-2 minimum-
 thickness bound applies only to volume-normalized *source* shells, not to plain shields, so this
 thickness is fit from a true zero.)  Truth data is the forward model itself, so the fit recovers it.
 */
BOOST_AUTO_TEST_CASE( FitNonSourceShieldThicknessFromZero )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const std::shared_ptr<const Material> iron = matdb->material( "Fe (iron)" );
  BOOST_REQUIRE( iron );

  const SandiaDecay::Nuclide * const ba133 = db->nuclide( "Ba133" );
  BOOST_REQUIRE( ba133 );

  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 600*PhysicalUnits::second;
  const double true_thickness = 1.0*PhysicalUnits::cm;     // truth: a 1 cm iron shell
  const double true_activity = 50.0*PhysicalUnits::microCi;
  const double age = 5.0*PhysicalUnits::year;

  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance, 5*PhysicalUnits::cm,
                                     PhysicalUnits::keV, 0, 3000*PhysicalUnits::keV,
                                     DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  // Ba133 81 -> 384 keV: iron attenuates 81 keV far more than 356 keV, so the line-ratio pins the
  //  thickness (and breaks the activity/thickness degeneracy a single line would have).
  const double placeholder_area = 1.0E4;
  deque<shared_ptr<const PeakDef>> truth_peaks;
  for( const double energy : { 80.9979, 276.3989, 302.8508, 356.0129, 383.8485 } )
    truth_peaks.push_back( make_test_peak( ba133, energy, 1.0, placeholder_area ) );

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;

  // Truth shielding: a NON-SOURCE 1 cm iron sphere around a Ba133 point source.
  ShieldingSourceFitCalc::ShieldingInfo iron_shield;
  iron_shield.m_geometry = GammaInteractionCalc::GeometryType::Spherical;
  iron_shield.m_isGenericMaterial = false;
  iron_shield.m_forFitting = true;
  iron_shield.m_material = make_shared<Material>( *iron );
  iron_shield.m_dimensions[0] = true_thickness;
  iron_shield.m_dimensions[1] = iron_shield.m_dimensions[2] = 0.0;
  iron_shield.m_fitDimensions[0] = iron_shield.m_fitDimensions[1] = iron_shield.m_fitDimensions[2] = false;

  ShieldingSourceFitCalc::SourceFitDef ba133_src;
  ba133_src.nuclide = ba133;
  ba133_src.activity = true_activity;
  ba133_src.fitActivity = false;
  ba133_src.age = age;
  ba133_src.fitAge = false;
  ba133_src.ageDefiningNuc = nullptr;
  ba133_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  chi_input.config.shieldings = { iron_shield };
  chi_input.config.sources = { ba133_src };
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = nullptr;
  chi_input.foreground_peaks = truth_peaks;
  chi_input.background_peaks = nullptr;

  // Data = the forward model's peak areas for 1 cm of iron.
  const deque<shared_ptr<const PeakDef>> expected_peaks = peaks_with_model_expected_areas( chi_input );

  // The actual fit: start the iron thickness at EXACTLY 0 (the point of the test) and fit it, plus
  //  the activity from an off guess.
  iron_shield.m_dimensions[0] = 0.0;
  iron_shield.m_fitDimensions[0] = true;
  ba133_src.activity = 0.5*true_activity;
  ba133_src.fitActivity = true;

  chi_input.config.shieldings = { iron_shield };
  chi_input.config.sources = { ba133_src };
  chi_input.foreground_peaks = expected_peaks;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;

  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  auto progress_fcn = [](){};
  bool finished_fit_called = false;
  auto finished_fcn = [&finished_fit_called](){ finished_fit_called = true; };

  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );

  BOOST_CHECK( finished_fit_called );
  BOOST_CHECK( results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );
  BOOST_REQUIRE_EQUAL( results->final_shieldings.size(), size_t(1) );
  BOOST_REQUIRE_EQUAL( results->fit_src_info.size(), size_t(1) );

  const double fit_thickness = results->final_shieldings[0].m_dimensions[0];
  const double fit_activity = results->fit_src_info[0].activity;

  cout << "FitNonSourceShieldThicknessFromZero: thickness=" << fit_thickness/PhysicalUnits::cm
       << " cm (truth " << true_thickness/PhysicalUnits::cm << "), activity="
       << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
       << " (truth " << PhysicalUnits::printToBestActivityUnits(true_activity,6) << ")" << endl;

  // Started at exactly 0; the gradient must have driven it to the truth 1 cm.
  BOOST_CHECK_MESSAGE( fabs(fit_thickness - true_thickness) < 0.02*true_thickness,
    "Fit iron thickness (" << fit_thickness/PhysicalUnits::cm << " cm) not close to truth (1 cm) - "
    "Ceres did not climb off a 0 starting thickness correctly" );
  BOOST_CHECK_MESSAGE( fabs(fit_activity - true_activity) < 0.05*true_activity,
    "Fit activity (" << PhysicalUnits::printToBestActivityUnits(fit_activity,6) << ") not close to truth" );
  BOOST_CHECK_MESSAGE( results->final_shieldings[0].m_dimensionUncerts[0] > 0.0,
    "Fit thickness uncertainty should be positive" );
}//BOOST_AUTO_TEST_CASE( FitNonSourceShieldThicknessFromZero )


/** Minuit2 baseline: fit source age (plus activity) for a Ra226 point source,
 where the in-growth of Rn222 progeny (Pb214/Bi214 lines vs Ra226s own 186 keV
 line) determines the age.  Synthetic peak areas are generated from the forward
 model, so the chi2 minimum is at known truth values.
 */
BOOST_AUTO_TEST_CASE( FitPointSourceAgeBaseline )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );

  const SandiaDecay::Nuclide * const ra226 = db->nuclide( "Ra226" );
  BOOST_REQUIRE( ra226 );

  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 600*PhysicalUnits::second;
  const double true_activity = 10.0*PhysicalUnits::microCi;
  const double true_age = 5.0*PhysicalUnits::day;

  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance,
                                    5*PhysicalUnits::cm, PhysicalUnits::keV,
                                    0, 3000*PhysicalUnits::keV,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  // 186.2 keV is Ra226 itself; 351.9 (Pb214) and 609.3 (Bi214) grow in with
  //  the ~3.8 day half-life of Rn222, so their ratio to 186.2 fixes the age.
  const double placeholder_area = 1.0E4;
  deque<shared_ptr<const PeakDef>> truth_peaks;
  for( const double energy : { 186.211, 351.932, 609.312 } )
    truth_peaks.push_back( make_test_peak( ra226, energy, 1.0, placeholder_area ) );

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;

  ShieldingSourceFitCalc::SourceFitDef ra226_src;
  ra226_src.nuclide = ra226;
  ra226_src.activity = true_activity;
  ra226_src.fitActivity = false;
  ra226_src.age = true_age;
  ra226_src.fitAge = false;
  ra226_src.ageDefiningNuc = nullptr;
  ra226_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  chi_input.config.shieldings = {};
  chi_input.config.sources = { ra226_src };
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = nullptr;
  chi_input.foreground_peaks = truth_peaks;
  chi_input.background_peaks = nullptr;

  const deque<shared_ptr<const PeakDef>> expected_peaks = peaks_with_model_expected_areas( chi_input );

  // Fit starting away from truth
  ra226_src.activity = 1.0*PhysicalUnits::microCi;
  ra226_src.fitActivity = true;
  ra226_src.age = 1.0*PhysicalUnits::day;
  ra226_src.fitAge = true;

  chi_input.config.sources = { ra226_src };
  chi_input.foreground_peaks = expected_peaks;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;

  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  auto progress_fcn = [](){};
  bool finished_fit_called = false;
  auto finished_fcn = [&finished_fit_called](){ finished_fit_called = true; };

  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );

  BOOST_CHECK( finished_fit_called );
  BOOST_CHECK( results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );
  BOOST_REQUIRE_EQUAL( results->fit_src_info.size(), size_t(1) );

  const double fit_activity = results->fit_src_info[0].activity;
  const double fit_age = results->fit_src_info[0].age;

  cout << "FitPointSourceAgeBaseline: activity="
       << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
       << " (truth " << PhysicalUnits::printToBestActivityUnits(true_activity,6)
       << "), age=" << fit_age/PhysicalUnits::day << " days (truth "
       << true_age/PhysicalUnits::day << " days)" << endl;

  BOOST_CHECK_MESSAGE( fabs(fit_activity - true_activity) < 0.05*true_activity,
                       "Fit activity (" << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
                       << ") not close to truth ("
                       << PhysicalUnits::printToBestActivityUnits(true_activity,6) << ")" );
  BOOST_CHECK_MESSAGE( fabs(fit_age - true_age) < 0.15*true_age,
                       "Fit age (" << fit_age/PhysicalUnits::day << " days) not close to truth ("
                       << true_age/PhysicalUnits::day << " days)" );

  BOOST_CHECK( results->fit_src_info[0].activityUncertainty.has_value()
               && (results->fit_src_info[0].activityUncertainty.value() > 0.0) );
  BOOST_CHECK( results->fit_src_info[0].ageUncertainty.has_value()
               && (results->fit_src_info[0].ageUncertainty.value() > 0.0) );
}//BOOST_AUTO_TEST_CASE( FitPointSourceAgeBaseline )


/** End-to-end fit with the source assembly off the detector axis: synthetic peak
 areas are generated from the forward model with offsets set, then activity and
 areal density are fit - the recovered values must match the truth (and would not,
 if the off-axis geometry were ignored: the true distance here is ~12% longer than
 the axial one).
 */
BOOST_AUTO_TEST_CASE( FitOffAxisSourceBaseline )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );

  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 600*PhysicalUnits::second;
  const double true_ad = 10.0*PhysicalUnits::g/PhysicalUnits::cm2;
  const double true_activity = 10.0*PhysicalUnits::microCi;

  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance,
                                    5*PhysicalUnits::cm, PhysicalUnits::keV,
                                    0, 3000*PhysicalUnits::keV,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;

  ShieldingSourceFitCalc::ShieldingInfo generic_shield;
  generic_shield.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  generic_shield.m_isGenericMaterial = true;
  generic_shield.m_forFitting = true;
  generic_shield.m_material = nullptr;
  generic_shield.m_dimensions[0] = 26.0;
  generic_shield.m_dimensions[1] = true_ad;
  generic_shield.m_dimensions[2] = 0.0;
  generic_shield.m_fitDimensions[0] = generic_shield.m_fitDimensions[1] = generic_shield.m_fitDimensions[2] = false;

  ShieldingSourceFitCalc::SourceFitDef ba133_src;
  ba133_src.nuclide = db->nuclide( "Ba133" );
  ba133_src.activity = true_activity;
  ba133_src.fitActivity = false;
  ba133_src.age = 5.0*PhysicalUnits::year;
  ba133_src.fitAge = false;
  ba133_src.ageDefiningNuc = nullptr;
  ba133_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

  deque<shared_ptr<const PeakDef>> peaks;
  for( const double energy : { 80.9979, 276.3989, 302.8508, 356.0129, 383.8485 } )
    peaks.push_back( make_test_peak( ba133_src.nuclide, energy, 1.0, 1.0E4 ) );

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  chi_input.config.source_offsets[0] = 40.0*PhysicalUnits::cm;
  chi_input.config.source_offsets[1] = -25.0*PhysicalUnits::cm;
  chi_input.config.shieldings = { generic_shield };
  chi_input.config.sources = { ba133_src };
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.foreground_peaks = peaks;

  // Truth peak areas, with the offsets in effect
  const deque<shared_ptr<const PeakDef>> expected_peaks = peaks_with_model_expected_areas( chi_input );

  // Fit activity and AD, starting away from truth (offsets stay user-set, never fit)
  generic_shield.m_dimensions[1] = 2.0*PhysicalUnits::g/PhysicalUnits::cm2;
  generic_shield.m_fitDimensions[1] = true;
  ba133_src.activity = 1.0*PhysicalUnits::microCi;
  ba133_src.fitActivity = true;

  chi_input.config.shieldings = { generic_shield };
  chi_input.config.sources = { ba133_src };
  chi_input.foreground_peaks = expected_peaks;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

  BOOST_CHECK_CLOSE( fcn_pars.first->trueSourceToDetectorDistance(),
                     sqrt( 100.0*100.0 + 40.0*40.0 + 25.0*25.0 )*PhysicalUnits::cm, 1.0E-9 );

  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;

  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  auto progress_fcn = [](){};
  auto finished_fcn = [](){};

  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );

  BOOST_REQUIRE( results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );
  BOOST_REQUIRE_EQUAL( results->final_shieldings.size(), size_t(1) );
  BOOST_REQUIRE_EQUAL( results->fit_src_info.size(), size_t(1) );

  const double true_ad_g_cm2 = true_ad / (PhysicalUnits::g/PhysicalUnits::cm2);
  const double fit_ad_g_cm2 = results->final_shieldings[0].m_dimensions[1];
  const double fit_activity = results->fit_src_info[0].activity;

  cout << "FitOffAxisSourceBaseline: AD=" << fit_ad_g_cm2 << " g/cm2 (truth " << true_ad_g_cm2
       << "), activity=" << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
       << " (truth " << PhysicalUnits::printToBestActivityUnits(true_activity,6) << ")" << endl;

  BOOST_CHECK_MESSAGE( fabs(fit_ad_g_cm2 - true_ad_g_cm2) < 0.02*true_ad_g_cm2,
                       "Off-axis fit AD (" << fit_ad_g_cm2 << " g/cm2) not close to truth ("
                       << true_ad_g_cm2 << " g/cm2)" );
  BOOST_CHECK_MESSAGE( fabs(fit_activity - true_activity) < 0.02*true_activity,
                       "Off-axis fit activity (" << PhysicalUnits::printToBestActivityUnits(fit_activity,6)
                       << ") not close to truth ("
                       << PhysicalUnits::printToBestActivityUnits(true_activity,6) << ")" );
}//BOOST_AUTO_TEST_CASE( FitOffAxisSourceBaseline )


/** XML round-trip of the source off-axis offsets in ShieldSourceConfig. */
BOOST_AUTO_TEST_CASE( SourceOffsetsSerialization )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  GammaInteractionCalc::ShieldSourceConfig config;
  config.distance = 150.0*PhysicalUnits::cm;
  config.geometry = GammaInteractionCalc::GeometryType::Rectangular;
  config.source_offsets[0] = 23.5*PhysicalUnits::cm;
  config.source_offsets[1] = -7.25*PhysicalUnits::cm;

  ShieldingSourceFitCalc::SourceFitDef src_def;
  src_def.nuclide = db->nuclide( "Cs137" );
  src_def.activity = 1.0*PhysicalUnits::microCi;
  src_def.fitActivity = true;
  src_def.age = 30.0*PhysicalUnits::day;
  src_def.fitAge = false;
  src_def.ageDefiningNuc = nullptr;
  src_def.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
  config.sources.push_back( src_def );

  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
  doc.append_node( root );
  BOOST_REQUIRE_NO_THROW( config.serialize( root ) );

  GammaInteractionCalc::ShieldSourceConfig reparsed;
  BOOST_REQUIRE_NO_THROW( reparsed.deSerialize( root ) );

  BOOST_CHECK_CLOSE( reparsed.source_offsets[0], config.source_offsets[0], 1.0E-6 );
  BOOST_CHECK_CLOSE( reparsed.source_offsets[1], config.source_offsets[1], 1.0E-6 );

  // And zero offsets dont emit the node (backward-compatible XML)
  GammaInteractionCalc::ShieldSourceConfig zero_config = config;
  zero_config.source_offsets[0] = zero_config.source_offsets[1] = 0.0;

  rapidxml::xml_document<char> doc2;
  rapidxml::xml_node<char> *root2 = doc2.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
  doc2.append_node( root2 );
  BOOST_REQUIRE_NO_THROW( zero_config.serialize( root2 ) );
  BOOST_CHECK( !root2->first_node( "SourceOffsets", 13 ) );

  GammaInteractionCalc::ShieldSourceConfig zero_reparsed;
  BOOST_REQUIRE_NO_THROW( zero_reparsed.deSerialize( root2 ) );
  BOOST_CHECK_EQUAL( zero_reparsed.source_offsets[0], 0.0 );
  BOOST_CHECK_EQUAL( zero_reparsed.source_offsets[1], 0.0 );
}//BOOST_AUTO_TEST_CASE( SourceOffsetsSerialization )


std::tuple<bool,int,int,vector<string>> test_fit_against_truth( const ShieldingSourceFitCalc::ModelFitResults &results,
                                                               const string &filename )
{
  bool successful = true;
  int numCorrect = 0, numTested = 0;
  vector<string> textInfoLines;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();    

  try
  {
    for( int i = 0; i < results.fit_src_info.size(); ++i )
    {
      const ShieldingSourceFitCalc::SourceFitDef &src = results.fit_src_info[i];
      const SandiaDecay::Nuclide *nuc = src.nuclide;
      BOOST_REQUIRE( nuc );
      
      const SandiaDecay::Nuclide *ageNuc = src.ageDefiningNuc;
      const ShieldingSourceFitCalc::ModelSourceType sourceType = src.sourceType;
      const bool selfAttNuc = (sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
      
      // For self-attenuating shieldings, we'll just test the shielding thickness
      // For nuclides whose age is controlled by another nuclide, we dont need to test age.
      if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
        continue;
      
      bool isSelfAtten = false;
      
      for( const auto &shield : results.initial_shieldings )
      {
        const SandiaDecay::Element *el = db->element( nuc->atomicNumber );
        auto el_nucs = shield.m_nuclideFractions_.find(el); 
        if( el_nucs != end(shield.m_nuclideFractions_) )
        {
          for( const auto &nuc_frac : el_nucs->second )
            isSelfAtten |= (std::get<0>(nuc_frac) == nuc);
        }
      }//for( const auto &shield : results.initial_shieldings )

      if( src.fitActivity && !isSelfAtten )
      {
        const boost::optional<double> &truthAct = src.truthActivity;
        const boost::optional<double> &tolerance = src.truthActivityTolerance;
        
        BOOST_CHECK_MESSAGE( truthAct.has_value(), "Did not have truth activity value for " + nuc->symbol + "." );
        BOOST_CHECK_MESSAGE( tolerance.has_value(), "Did not have truth activity tolerance for " + nuc->symbol + "." );
        
        
        if( !truthAct || !tolerance )
        {
          successful = false;
          textInfoLines.push_back( "Did not have truth value for " + nuc->symbol + " activity." );
          continue;
        }
        
        const bool closeEnough = (fabs(*truthAct - src.activity) < *tolerance);
        
        numTested += 1;
        numCorrect += closeEnough;
        
        textInfoLines.push_back( "For " + nuc->symbol + " fit activity "
                                + PhysicalUnits::printToBestActivityUnits(src.activity) + " with the"
                                " truth value of "
                                + PhysicalUnits::printToBestActivityUnits(*truthAct)
                                + " and tolerance "
                                + PhysicalUnits::printToBestActivityUnits(*tolerance)
                                + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                );
        
        BOOST_CHECK_MESSAGE( fabs(*truthAct - src.activity) < *tolerance, textInfoLines.back() );
      }//if( fit activity )
      
      if( src.fitAge )
      {
        const boost::optional<double> &truthAge = src.truthAge;
        const boost::optional<double> &tolerance = src.truthAgeTolerance;
        
        BOOST_CHECK_MESSAGE( truthAge.has_value(), "Did not have truth age value for " + nuc->symbol + " age." );
        BOOST_CHECK_MESSAGE( tolerance.has_value(), "Did not have truth age tolerance for " + nuc->symbol + " age." );
      
        
        if( !truthAge || !tolerance )
        {
          successful = false;
          textInfoLines.push_back( "Did not have truth value for " + nuc->symbol + " age." );
          continue;
        }
        
        const bool closeEnough = (fabs(*truthAge - src.age) < *tolerance);
        
        numTested += 1;
        numCorrect += closeEnough;
        
        textInfoLines.push_back( "For " + nuc->symbol + " fit age "
                                + PhysicalUnits::printToBestTimeUnits(src.age) + " with the"
                                " truth value of "
                                + PhysicalUnits::printToBestTimeUnits(*truthAge)
                                + " and tolerance "
                                + PhysicalUnits::printToBestTimeUnits(*tolerance)
                                + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                );
        
        BOOST_CHECK_MESSAGE( fabs(*truthAge - src.age) < *tolerance, textInfoLines.back() );
      }//if( fit age )
    }//for( int i = 0; i < nnuc; ++i )
    
    for( const ShieldingSourceFitCalc::FitShieldingInfo &shield : results.final_shieldings )
    {
      if( shield.m_isGenericMaterial )
      {
        if( shield.m_fitDimensions[1] )
        {
          const boost::optional<double> &truth = shield.m_truthDimensions[1];
          const boost::optional<double> &tolerance = shield.m_truthDimensionsTolerances[1];
          
          BOOST_CHECK_MESSAGE( truth.has_value(), "Did not have truth value AD for generic shielding" );
          BOOST_CHECK_MESSAGE( tolerance.has_value(), "Did not have truth tolerance AD for generic shielding" );
          
          if( !truth || !tolerance )
          {
            successful = false;
            textInfoLines.push_back( "Did not have truth AD for generic shielding" );
            continue;
          }
          
          const double fitAD = shield.m_dimensions[1];
          const bool closeEnough = (fabs(*truth - fitAD) < *tolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For Generic Shielding fit AN " + std::to_string(fitAD)
                                  + " with the truth value of " + std::to_string(*truth)
                                  + " and tolerance " + std::to_string(*tolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
          
          BOOST_CHECK_MESSAGE( fabs(*truth - fitAD) < *tolerance, textInfoLines.back() );
        }//if( fit AD )
        
        if( shield.m_fitDimensions[0] )
        {
          const boost::optional<double> &truth = shield.m_truthDimensions[0];
          const boost::optional<double> &tolerance = shield.m_truthDimensionsTolerances[0];
          
          BOOST_CHECK_MESSAGE( truth.has_value(), "Did not have truth value AN for generic shielding" );
          BOOST_CHECK_MESSAGE( tolerance.has_value(), "Did not have truth tolerance AN for generic shielding" );
          
          if( !truth || !tolerance )
          {
            successful = false;
            textInfoLines.push_back( "Did not have truth AN for generic shielding" );
            continue;
          }
          
          const double fitAN = shield.m_dimensions[0];
          const bool closeEnough = (fabs(*truth - fitAN) < *tolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For Generic Shielding fit AN " + std::to_string(fitAN)
                                  + " with the truth value of " + std::to_string(*truth)
                                  + " and tolerance " + std::to_string(*tolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
          
          BOOST_CHECK_MESSAGE( fabs(*truth - fitAN) < *tolerance, textInfoLines.back() );
        }//if( fit AN )
      }else
      {
        shared_ptr<const Material> mat = shield.m_material;
        
        BOOST_CHECK( results.geometry == shield.m_geometry );
        
        BOOST_CHECK( mat );
        if( !mat )
        {
          successful = false;
          textInfoLines.push_back( "There was an invalid material." );
          continue;
        }
        
        const auto geom = shield.m_geometry;
        auto checkDimension = [shield, geom, mat, &filename, &successful, &textInfoLines, &numTested, &numCorrect]( const int dim ){
          
          if( !shield.m_fitDimensions[dim] )
            return;
          
          string labeltxt;
          const double fitValue = shield.m_dimensions[dim];
          const boost::optional<double> &thicknessVal = shield.m_truthDimensions[dim];
          const boost::optional<double> &toleranceVal = shield.m_truthDimensionsTolerances[dim];
          
          switch( geom )
          {
            case GammaInteractionCalc::GeometryType::Spherical:
              BOOST_REQUIRE( dim < 1 );
              labeltxt = "thickness";
              break;
              
            case GammaInteractionCalc::GeometryType::CylinderEndOn:
            case GammaInteractionCalc::GeometryType::CylinderSideOn:
              BOOST_REQUIRE( dim < 2 );
              if( dim == 0 )
                labeltxt = "cyl. radius";
              else if( dim == 1 )
                labeltxt = "cyl. length";
              break;
              
            case GammaInteractionCalc::GeometryType::Rectangular:
              BOOST_REQUIRE( dim < 3 );
              if( dim == 0 )
                labeltxt = "rect. width";
              else if( dim == 1 )
                labeltxt = "rect. height";
              else if( dim == 2 )
                labeltxt = "rect. depth";
              
              break;
              
            case GammaInteractionCalc::GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( geometry() )
          
          
          BOOST_CHECK_MESSAGE( thicknessVal.has_value() && toleranceVal.has_value(),
                              "Missing truth " + labeltxt + " for shielding '" + mat->name + "'" );
          
          if( !thicknessVal.has_value() || !toleranceVal.has_value() )
          {
            successful = false;
            textInfoLines.push_back( "Missing truth " + labeltxt + " for shielding '" + mat->name + "'" );
            return;
          }
          
          const bool closeEnough = (fabs((*thicknessVal) - fitValue) < (*toleranceVal));
          if( !closeEnough )
          {
            cerr << "Failing for " << labeltxt << " '" << filename << "'." << endl;
            cerr << "shield.m_dimensions[] = [" << PhysicalUnits::printToBestLengthUnits(shield.m_dimensions[0],4)
            << ", " << PhysicalUnits::printToBestLengthUnits(shield.m_dimensions[1],4)
            << ", " << PhysicalUnits::printToBestLengthUnits(shield.m_dimensions[2],4)
            << "]" << endl;
            cerr <<endl;
          }
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For shielding '" + mat->name + "' fit " + labeltxt + " "
                                  + PhysicalUnits::printToBestLengthUnits(fitValue,4)
                                  + " with the truth value of "
                                  + PhysicalUnits::printToBestLengthUnits(*thicknessVal,4)
                                  + " and tolerance "
                                  + PhysicalUnits::printToBestLengthUnits(*toleranceVal)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
          
          BOOST_CHECK_MESSAGE( fabs((*thicknessVal) - fitValue) < (*toleranceVal), textInfoLines.back() );
        };//auto checkDimension
        
        
        switch( geom )
        {
          case GammaInteractionCalc::GeometryType::Spherical:
            checkDimension(0);
            break;
            
          case GammaInteractionCalc::GeometryType::CylinderEndOn:
          case GammaInteractionCalc::GeometryType::CylinderSideOn:
            checkDimension(0);
            checkDimension(1);
            break;
            
          case GammaInteractionCalc::GeometryType::Rectangular:
            checkDimension(0);
            checkDimension(1);
            checkDimension(2);
            break;
            
          case GammaInteractionCalc::GeometryType::NumGeometryType:
            break;
        }//switch( geom )
        
        auto checkMassFrac = [shield, mat, &successful, &textInfoLines, &numTested, &numCorrect](
                                    const SandiaDecay::Nuclide * const nuc, 
                                    const SandiaDecay::Element * const el, 
                                    double value ){
          BOOST_REQUIRE( el );
          const auto truth_nucs_pos = shield.m_truthFitMassFractions.find(el);
          BOOST_CHECK_MESSAGE( truth_nucs_pos != end(shield.m_truthFitMassFractions),
                              "Missing truth mass-fraction for " + el->symbol );

          if( truth_nucs_pos == end(shield.m_truthFitMassFractions) )
          {
            successful = false;
            textInfoLines.push_back( "Missing truth mass-fractions for element " + el->symbol );
            return;
          }

          const map<const SandiaDecay::Nuclide *,pair<double,double>> &nucs = truth_nucs_pos->second;
          auto truth_pos = nucs.find(nuc);
          BOOST_CHECK_MESSAGE( truth_pos != end(nucs),
                                "Missing truth mass-fraction for " + (nuc ? nuc->symbol : "other nuclides") );
          if( truth_pos == end(nucs) )
          {
            successful = false;
            textInfoLines.push_back( "Missing truth mass-fraction for " + (nuc ? nuc->symbol : "other nuclides") );
            return;
          }                           
          
                                  
          const double truthval = get<0>(truth_pos->second);
          const double truthtol = get<1>(truth_pos->second);
          BOOST_CHECK_MESSAGE( (truthval >= 0.0) && (truthval <= 1.0) && (truthtol >= 0.0) && (truthtol <= 1.0),
                                "Invalid truth mass-fraction for " + (nuc ? nuc->symbol : "other nuclides") );
                                      
          if( (truthval < 0.0) || (truthval > 1.0) || (truthtol < 0.0) || (truthtol > 1.0) )
          {
            successful = false;
            textInfoLines.push_back( "Invalid truth mass-fraction for " + (nuc ? nuc->symbol : "other nuclides") );
            return;
          }
            
          numTested += 1;
          const bool closeEnough = ((value >= (truthval - truthtol)) && (value <= (truthval + truthtol)));
          numCorrect += closeEnough;
            
          textInfoLines.push_back( "For shielding '" + mat->name + "' fit " + (nuc ? nuc->symbol : "other nuclides")
                                     + " to have mass fraction " + SpecUtils::printCompact(value, 5)
                                     + " with the truth value of " + SpecUtils::printCompact(truthval,5)
                                     + " and tolerance " + SpecUtils::printCompact(truthtol,5)
                                     + (closeEnough ? " - within tolerance." : " - out of tolerance." ) );
                                      
          BOOST_CHECK_MESSAGE( (value >= (truthval - truthtol)) && (value <= (truthval + truthtol)),
                                textInfoLines.back() );
        };//checkMassFrac( ... )
          
        for( const auto &el_nucs : shield.m_nuclideFractions_ )
        {
          const SandiaDecay::Element * const el = el_nucs.first;
          for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_frac : el_nucs.second )
          {
            checkMassFrac( std::get<0>(nuc_frac), el, std::get<1>(nuc_frac) );
          }
        }//for( const auto &prev_nuc_frac : shield.m_nuclideFractionUncerts )
      }//if( generic material ) / else
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    successful = (successful && numTested);
  }catch( std::exception &e )
  {
    successful = false;
    textInfoLines.push_back( "Caught exception during testing: " + string(e.what()) );
  }//try / catch
  
  return tuple<bool,int,int,vector<string>>( successful, numCorrect, numTested, textInfoLines );
}//std::tuple<bool,int,int,vector<string>> test_fit_against_truth( const ShieldingSourceFitCalc::ModelFitResults &results )



void set_fit_quantities_to_default_values( vector<ShieldingSourceFitCalc::ShieldingInfo> &shield_definitions,
                                     vector<ShieldingSourceFitCalc::SourceFitDef> &src_definitions )
{
  for( ShieldingSourceFitCalc::SourceFitDef &src : src_definitions )
  {
    const SandiaDecay::Nuclide *nuc = src.nuclide;
    BOOST_REQUIRE( nuc );
    
    const SandiaDecay::Nuclide *ageNuc = src.ageDefiningNuc;
    const ShieldingSourceFitCalc::ModelSourceType sourceType = src.sourceType;
    const bool selfAttNuc = (sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
    
    // For self-attenuating shieldings, we'll just test the shielding thickness
    // For nuclides whose age is controlled by another nuclide, we dont need to test age.
    if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
      continue;
    
    if( src.fitActivity )
      src.activity = 0.001*PhysicalUnits::curie;
    
    if( src.fitAge )
      src.age = PeakDef::defaultDecayTime( nuc, nullptr );
  }//for( int i = 0; i < nnuc; ++i )
  
  
  for( ShieldingSourceFitCalc::ShieldingInfo &shield : shield_definitions )
  {
    if( shield.m_isGenericMaterial )
    {
      if( shield.m_fitDimensions[0] )
        shield.m_dimensions[0] = 26; // Atomic Number
      
      if( shield.m_fitDimensions[1] )
        shield.m_dimensions[1] = 0;  // Areal Density
    }else
    {
      shared_ptr<const Material> mat = shield.m_material;
      if( !mat )
        continue;
      
      // Set mass-fractions to all be even, if we are fitting for them
      for( auto &el_nucs : shield.m_nuclideFractions_ )
      {
        int num_fit = 0;
        double frac_fit_sum = 0.0;
          
        for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_frac : el_nucs.second )
        {
          num_fit += std::get<2>(nuc_frac);
          if( std::get<2>(nuc_frac) )
            frac_fit_sum += std::get<1>(nuc_frac);
        }
        assert( num_fit != 1 );
        for( tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_frac : el_nucs.second )
        {
          if( std::get<2>(nuc_frac) )
            std::get<1>(nuc_frac) = (frac_fit_sum / num_fit);
        }
      }//for( auto &el_nucs : shield.m_nuclideFractions_ )
      
      
      switch( shield.m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          if( shield.m_fitDimensions[0] )
            shield.m_dimensions[0] = 1.0*PhysicalUnits::cm;
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          if( shield.m_fitDimensions[0] )
            shield.m_dimensions[0] = 0.5*PhysicalUnits::cm;
          if( shield.m_fitDimensions[1] )
            shield.m_dimensions[1] = 0.5*PhysicalUnits::cm;
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          if( shield.m_fitDimensions[0] )
            shield.m_dimensions[0] = 0.5*PhysicalUnits::cm;
          if( shield.m_fitDimensions[1] )
            shield.m_dimensions[1] = 0.5*PhysicalUnits::cm;
          if( shield.m_fitDimensions[2] )
            shield.m_dimensions[2] = 0.5*PhysicalUnits::cm;
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( geometry() )
      
      for( ShieldingSourceFitCalc::TraceSourceInfo &trace : shield.m_traceSources )
      {
        if( !trace.m_fitActivity )
          continue;
        
        switch( trace.m_type )
        {
          case GammaInteractionCalc::TraceActivityType::TotalActivity:
            trace.m_activity = 0.001*PhysicalUnits::curie;;
            break;
          case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
          case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
          case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
            trace.m_activity = 1.0E-6*PhysicalUnits::curie;;
            break;
            
          case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
            assert( 0 );
            break;
        }//switch( trace.m_type )
      }//switch( trace.m_type )
    }//if( generic material ) / else
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void set_fit_quantities_to_default_values(...)


BOOST_AUTO_TEST_CASE( FitAnalystTraceSource )
{
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
 
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const std::shared_ptr<const Material> iron = matdb->material( "Fe (iron)" );
  BOOST_REQUIRE( iron );

  const string test_n42_file = SpecUtils::append_path(g_test_file_dir, "../analysis_tests/AEGIS_Eu152_surface_contamination.n42_20230622T113239.276178.n42");
  BOOST_REQUIRE( SpecUtils::is_file(test_n42_file) );
  
  
  
  SpecMeas specfile;
  BOOST_REQUIRE( specfile.load_N42_file( test_n42_file ) );
  BOOST_REQUIRE( specfile.sample_numbers().size() == 1 );
  BOOST_REQUIRE( specfile.num_measurements() == 1 );
  
  shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = specfile.peaks( specfile.sample_numbers() );
  BOOST_REQUIRE( peaks && (peaks->size() > 10) );
  
  shared_ptr<const SpecUtils::Measurement> foreground = specfile.measurements()[0];
  BOOST_REQUIRE( foreground );
  
  shared_ptr<const DetectorPeakResponse> detector = specfile.detector();
  BOOST_REQUIRE( detector && detector->isValid() );
  
  rapidxml::xml_document<char> *model_xml = specfile.shieldingSourceModel();
  BOOST_REQUIRE( model_xml );
  rapidxml::xml_node<char> *base_node = model_xml->first_node();
  BOOST_REQUIRE( base_node );
  
  //string xml;
  //rapidxml::print(std::back_inserter(xml), *model_xml, 0);
  //cout << "XML: " << xml << endl;
  
  const rapidxml::xml_node<char> *geom_node = base_node->first_node( "Geometry" );
  BOOST_REQUIRE( geom_node );
  const rapidxml::xml_node<char> *dist_node = base_node->first_node( "Distance" );
  BOOST_REQUIRE( dist_node );
  
  const string diststr = SpecUtils::xml_value_str( dist_node );
  double distance;
  BOOST_REQUIRE_NO_THROW( distance = PhysicalUnits::stringToDistance( diststr ) );
  
  const string geomstr = SpecUtils::xml_value_str( geom_node );
  GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  for( GammaInteractionCalc::GeometryType type = GammaInteractionCalc::GeometryType(0);
      type != GammaInteractionCalc::GeometryType::NumGeometryType;
      type = GammaInteractionCalc::GeometryType(static_cast<int>(type) + 1) )
  {
    if( SpecUtils::iequals_ascii(geomstr, GammaInteractionCalc::to_str(type)) )
    {
      geometry = type;
      break;
    }
  }
  
  BOOST_REQUIRE( geometry != GammaInteractionCalc::GeometryType::NumGeometryType );
  
  vector<ShieldingSourceFitCalc::ShieldingInfo> shield_definitions;
  vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  BOOST_REQUIRE_NO_THROW( options.deSerialize( base_node ) );
  
  const rapidxml::xml_node<char> *shieldings_node = base_node->first_node( "Shieldings" );
  BOOST_REQUIRE( shieldings_node );
  
  XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
  {
    ShieldingSourceFitCalc::ShieldingInfo info;
    BOOST_REQUIRE_NO_THROW( info.deSerialize( shield_node ) );
    shield_definitions.push_back( info );
  }
  
  BOOST_REQUIRE( shield_definitions.size() == 1 );
  
  const rapidxml::xml_node<char> *srcs_node = base_node->first_node( "Nuclides" );
  BOOST_REQUIRE( srcs_node );
  
  XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
  {
    ShieldingSourceFitCalc::SourceFitDef info;
    BOOST_REQUIRE_NO_THROW( info.deSerialize( src_node ) );
    src_definitions.push_back( info );
  }
  BOOST_REQUIRE( src_definitions.size() == 1 );

  GammaInteractionCalc::ShieldSourceConfig parsed_config;
  BOOST_REQUIRE_NO_THROW( parsed_config.deSerialize( base_node ) );

  BOOST_CHECK_SMALL( fabs(parsed_config.distance - distance), 1e-9 * (std::max)(1.0, fabs(distance)) );
  BOOST_CHECK_EQUAL( static_cast<int>(parsed_config.geometry), static_cast<int>(geometry) );
  BOOST_CHECK_EQUAL( parsed_config.shieldings.size(), shield_definitions.size() );
  std::vector<ShieldingSourceFitCalc::SourceFitDef> parsed_source_defs = parsed_config.sources;
  BOOST_CHECK_EQUAL( parsed_source_defs.size(), src_definitions.size() );

#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  for( size_t i = 0; i < shield_definitions.size(); ++i )
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingInfo::equalEnough( shield_definitions[i], parsed_config.shieldings[i] ) );
  for( size_t i = 0; i < src_definitions.size(); ++i )
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::SourceFitDef::equalEnough( src_definitions[i], parsed_source_defs[i] ) );
  BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( options, parsed_config.options ) );
#endif

  rapidxml::xml_document<char> roundtrip_doc;
  rapidxml::xml_node<char> *round_root = roundtrip_doc.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
  roundtrip_doc.append_node( round_root );
  BOOST_REQUIRE_NO_THROW( parsed_config.serialize( round_root ) );

  GammaInteractionCalc::ShieldSourceConfig reparsed_config;
  BOOST_REQUIRE_NO_THROW( reparsed_config.deSerialize( round_root ) );

  BOOST_CHECK_SMALL( fabs(reparsed_config.distance - parsed_config.distance), 1e-9 * (std::max)(1.0, fabs(parsed_config.distance)) );
  BOOST_CHECK_EQUAL( static_cast<int>(reparsed_config.geometry), static_cast<int>(parsed_config.geometry) );
  BOOST_CHECK_EQUAL( reparsed_config.shieldings.size(), parsed_config.shieldings.size() );
  std::vector<ShieldingSourceFitCalc::SourceFitDef> reparsed_source_defs = reparsed_config.sources;
  BOOST_CHECK_EQUAL( reparsed_source_defs.size(), parsed_source_defs.size() );

#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  for( size_t i = 0; i < reparsed_config.shieldings.size(); ++i )
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingInfo::equalEnough( reparsed_config.shieldings[i], parsed_config.shieldings[i] ) );
  for( size_t i = 0; i < reparsed_source_defs.size(); ++i )
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::SourceFitDef::equalEnough( reparsed_source_defs[i], parsed_source_defs[i] ) );
  BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( reparsed_config.options, parsed_config.options ) );
#endif

  BOOST_CHECK_EQUAL( reparsed_config.sources.size(), parsed_config.sources.size() );
  
  set_fit_quantities_to_default_values( shield_definitions, src_definitions );
  
  // We have all the parts, lets do the computation:
  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = geometry;
  chi_input.config.shieldings = shield_definitions;
  chi_input.config.sources = src_definitions;
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = nullptr;
  chi_input.foreground_peaks.assign( peaks->begin(), peaks->end() );
  chi_input.background_peaks = nullptr;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
                GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
 
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;
  
  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  
  
  auto progress_fcn = [progress](){
    // Probably wont get called, since its a simple fit
  };
  
  bool finished_fit_called = false;
  auto finished_fcn = [results, &finished_fit_called](){
    finished_fit_called = true;
  };
  
  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );
  
  BOOST_CHECK( finished_fit_called );
  
  BOOST_CHECK_MESSAGE( results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final,
                      "Fit status was not successful" );
  
  BOOST_REQUIRE( results->fit_src_info.size() >= 1 );
  
  tuple<bool,int,int,vector<string>> test_results = test_fit_against_truth( *results, test_n42_file );
  const bool successful = get<0>(test_results);
  const int numCorrect = get<1>(test_results);
  const int numTested = get<2>(test_results);
  const vector<string> &textInfoLines = get<3>(test_results);
  
  cout << (successful ? "Successfully" : "Unsuccessfully") << " checked against truth values with "
       << numCorrect << " tests passing, out of " << numTested << ".\nNotes:" << endl;
  for( const auto &msg : textInfoLines )
    cout << "\t" << msg << endl;
}//BOOST_AUTO_TEST_CASE( FitAnalystTraceSource )


namespace
{
/** Loads an analyst-test N42 (single-measurement style, like the AEGIS trace-source file)
 into a ShieldSourceInput; returns false if the file doesnt fit that simple structure.
 */
bool load_simple_analyst_n42( const string &n42_path,
                              GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput &chi_input )
{
  auto specfile = make_shared<SpecMeas>();
  if( !specfile->load_N42_file( n42_path ) || !specfile->num_measurements() )
    return false;

  // Use the first sample-set that has peaks (the analyst files store peaks for the foreground)
  shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks;
  set<int> foreground_samples;
  for( const set<int> &samples : specfile->sampleNumsWithPeaks() )
  {
    shared_ptr<const deque<shared_ptr<const PeakDef>>> these_peaks = specfile->peaks( samples );
    if( these_peaks && !these_peaks->empty() )
    {
      peaks = these_peaks;
      foreground_samples = samples;
      break;
    }
  }//for( loop over sample sets with peaks )

  if( !peaks || peaks->empty() )
    return false;

  shared_ptr<const DetectorPeakResponse> detector = specfile->detector();
  if( !detector || !detector->isValid() )
    return false;

  rapidxml::xml_document<char> *model_xml = specfile->shieldingSourceModel();
  if( !model_xml || !model_xml->first_node() )
    return false;

  GammaInteractionCalc::ShieldSourceConfig config;
  config.deSerialize( model_xml->first_node() );
  set_fit_quantities_to_default_values( config.shieldings, config.sources );

  shared_ptr<const SpecUtils::Measurement> foreground;
  if( !foreground_samples.empty() )
  {
    const vector<string> &det_names = specfile->detector_names();
    foreground = specfile->measurement( *begin(foreground_samples),
                                        det_names.empty() ? string() : det_names[0] );
  }
  if( !foreground )
    foreground = specfile->measurements()[0];

  chi_input.config = config;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = nullptr;
  chi_input.foreground_peaks.assign( peaks->begin(), peaks->end() );
  chi_input.background_peaks = nullptr;

  return true;
}//load_simple_analyst_n42(...)


/** Checks `expected_peak_counts_imp<double>` reproduces the `expectedCounts` of
 the legacy `energy_chi_contributions(...)` at the fits initial parameters.

 `tolerance` is relative; volumetric (trace/self-atten) cases need a looser one
 since the integration backend differs (adaptive GL vs Cuhre).
 */
void check_expected_counts_parity( const GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput &chi_input,
                                   const double tolerance,
                                   const string &label )
{
  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

  const vector<double> params = fcn_pars.second.Params();

  GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache cache_legacy, cache_imp;
  const vector<GammaInteractionCalc::PeakResultPlotInfo> legacy
        = fcn_pars.first->energy_chi_contributions( params, {}, cache_legacy, nullptr, nullptr );
  const vector<double> templated
        = fcn_pars.first->expected_peak_counts_imp<double>( params, cache_imp );

  BOOST_REQUIRE_EQUAL( legacy.size(), templated.size() );

  for( size_t i = 0; i < legacy.size(); ++i )
  {
    BOOST_CHECK_MESSAGE( fabs(templated[i] - legacy[i].expectedCounts)
                            <= tolerance*std::max(legacy[i].expectedCounts, 1.0E-12),
                         label << ": peak " << i << " (" << legacy[i].energy << " keV): templated="
                         << std::setprecision(10) << templated[i] << " vs legacy="
                         << legacy[i].expectedCounts );
  }//for( loop over peaks )
}//check_expected_counts_parity(...)
}//namespace


/** Validates the templated chi2 core (`expected_peak_counts_imp<double>`) against the
 legacy `energy_chi_contributions(...)` on point-source, aged-source, and volumetric
 (trace-source) configurations.
 */
BOOST_AUTO_TEST_CASE( ExpectedPeakCountsImpParity )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 600*PhysicalUnits::second;

  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance,
                                    5*PhysicalUnits::cm, PhysicalUnits::keV,
                                    0, 3000*PhysicalUnits::keV,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;

  {// Case A: Ba133 point source behind a generic shielding
    ShieldingSourceFitCalc::ShieldingInfo generic_shield;
    generic_shield.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
    generic_shield.m_isGenericMaterial = true;
    generic_shield.m_forFitting = true;
    generic_shield.m_material = nullptr;
    generic_shield.m_dimensions[0] = 26.0;
    generic_shield.m_dimensions[1] = 10.0*PhysicalUnits::g/PhysicalUnits::cm2;
    generic_shield.m_dimensions[2] = 0.0;
    generic_shield.m_fitDimensions[0] = generic_shield.m_fitDimensions[1] = generic_shield.m_fitDimensions[2] = false;

    ShieldingSourceFitCalc::SourceFitDef ba133_src;
    ba133_src.nuclide = db->nuclide( "Ba133" );
    ba133_src.activity = 10.0*PhysicalUnits::microCi;
    ba133_src.fitActivity = true;
    ba133_src.age = 5.0*PhysicalUnits::year;
    ba133_src.fitAge = false;
    ba133_src.ageDefiningNuc = nullptr;
    ba133_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

    deque<shared_ptr<const PeakDef>> peaks;
    for( const double energy : { 80.9979, 276.3989, 302.8508, 356.0129, 383.8485 } )
      peaks.push_back( make_test_peak( ba133_src.nuclide, energy, 1.0, 1.0E4 ) );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    chi_input.config.distance = distance;
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
    chi_input.config.shieldings = { generic_shield };
    chi_input.config.sources = { ba133_src };
    chi_input.config.options = options;
    chi_input.detector = detector;
    chi_input.foreground = foreground;
    chi_input.foreground_peaks = peaks;

    check_expected_counts_parity( chi_input, 1.0E-5, "Ba133+generic" );
  }// End Case A

  {// Case B: aged Ra226 point source, no shielding
    ShieldingSourceFitCalc::SourceFitDef ra226_src;
    ra226_src.nuclide = db->nuclide( "Ra226" );
    ra226_src.activity = 10.0*PhysicalUnits::microCi;
    ra226_src.fitActivity = true;
    ra226_src.age = 5.0*PhysicalUnits::day;
    ra226_src.fitAge = true;
    ra226_src.ageDefiningNuc = nullptr;
    ra226_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

    deque<shared_ptr<const PeakDef>> peaks;
    for( const double energy : { 186.211, 351.932, 609.312 } )
      peaks.push_back( make_test_peak( ra226_src.nuclide, energy, 1.0, 1.0E4 ) );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    chi_input.config.distance = distance;
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
    chi_input.config.shieldings = {};
    chi_input.config.sources = { ra226_src };
    chi_input.config.options = options;
    chi_input.detector = detector;
    chi_input.foreground = foreground;
    chi_input.foreground_peaks = peaks;

    check_expected_counts_parity( chi_input, 1.0E-9, "Ra226-aged" );
  }// End Case B

  {// Case C: the AEGIS Eu152 trace-source (volumetric) analyst file
    const string test_n42_file = SpecUtils::append_path( g_test_file_dir,
                  "../analysis_tests/AEGIS_Eu152_surface_contamination.n42_20230622T113239.276178.n42" );
    BOOST_REQUIRE( SpecUtils::is_file(test_n42_file) );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    BOOST_REQUIRE( load_simple_analyst_n42( test_n42_file, chi_input ) );

    // Volumetric: integration backend differs (adaptive GL vs Cuhre), so looser tolerance
    check_expected_counts_parity( chi_input, 2.0E-3, "AEGIS-trace" );
  }// End Case C

  {// Case D: self-attenuating U+Np source with fit mass fractions in two elements
    const string test_n42_file = SpecUtils::append_path( g_test_file_dir,
                                                         "../analysis_tests/u_np_1kg.n42" );
    BOOST_REQUIRE( SpecUtils::is_file(test_n42_file) );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    BOOST_REQUIRE( load_simple_analyst_n42( test_n42_file, chi_input ) );

    check_expected_counts_parity( chi_input, 2.0E-3, "U-Np-self-atten" );
  }// End Case D
}//BOOST_AUTO_TEST_CASE( ExpectedPeakCountsImpParity )


/** Validates the ceres::Jet derivative lanes of `expected_peak_counts_imp` - including
 the numerically-seeded age derivative - against central finite differences of the
 double evaluation.
 */
BOOST_AUTO_TEST_CASE( ExpectedPeakCountsJetDerivatives )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );

  const double distance = 100*PhysicalUnits::cm;
  const float live_time = 600*PhysicalUnits::second;

  auto detector = make_shared<DetectorPeakResponse>();
  detector->fromExpOfLogPowerSeries( {0.0f, 0.0f}, {}, distance,
                                    5*PhysicalUnits::cm, PhysicalUnits::keV,
                                    0, 3000*PhysicalUnits::keV,
                                    DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );

  auto foreground = make_shared<SpecUtils::Measurement>();
  auto spec = make_shared<vector<float>>( 16, 1.0f );
  foreground->set_gamma_counts( spec, live_time, live_time );

  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;

  // Aged Ra226 point source behind a generic shielding: parameters are
  //  {activity, age, AN, AD} - exercising the age numeric-diff seeding and the
  //  AN interpolation derivative in one go.
  ShieldingSourceFitCalc::ShieldingInfo generic_shield;
  generic_shield.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  generic_shield.m_isGenericMaterial = true;
  generic_shield.m_forFitting = true;
  generic_shield.m_material = nullptr;
  generic_shield.m_dimensions[0] = 31.4;  //non-integer, to exercise the AN interpolation
  generic_shield.m_dimensions[1] = 5.0*PhysicalUnits::g/PhysicalUnits::cm2;
  generic_shield.m_dimensions[2] = 0.0;
  generic_shield.m_fitDimensions[0] = true;
  generic_shield.m_fitDimensions[1] = true;
  generic_shield.m_fitDimensions[2] = false;

  ShieldingSourceFitCalc::SourceFitDef ra226_src;
  ra226_src.nuclide = db->nuclide( "Ra226" );
  ra226_src.activity = 10.0*PhysicalUnits::microCi;
  ra226_src.fitActivity = true;
  ra226_src.age = 5.0*PhysicalUnits::day;
  ra226_src.fitAge = true;
  ra226_src.ageDefiningNuc = nullptr;
  ra226_src.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

  deque<shared_ptr<const PeakDef>> peaks;
  for( const double energy : { 186.211, 351.932, 609.312 } )
    peaks.push_back( make_test_peak( ra226_src.nuclide, energy, 1.0, 1.0E4 ) );

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  chi_input.config.distance = distance;
  chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  chi_input.config.shieldings = { generic_shield };
  chi_input.config.sources = { ra226_src };
  chi_input.config.options = options;
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.foreground_peaks = peaks;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
  const shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> &fcn = fcn_pars.first;
  const vector<double> params = fcn_pars.second.Params();

  // Params: [0]=activity (MBq-scaled), [1]=age (seconds), [2]=AN, [3]=AD, [4]=unused
  BOOST_REQUIRE_EQUAL( params.size(), size_t(5) );

  using Jet4 = ceres::Jet<double,4>;

  vector<Jet4> jet_params( params.size() );
  for( size_t i = 0; i < params.size(); ++i )
  {
    jet_params[i] = Jet4( params[i] );
    if( i < 4 )
      jet_params[i].v[i] = 1.0;
  }

  GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache jet_cache;
  const vector<Jet4> jet_counts = fcn->expected_peak_counts_imp<Jet4>( jet_params, jet_cache );

  GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache dbl_cache;
  const vector<double> dbl_counts = fcn->expected_peak_counts_imp<double>( params, dbl_cache );

  BOOST_REQUIRE_EQUAL( jet_counts.size(), dbl_counts.size() );
  BOOST_REQUIRE_EQUAL( jet_counts.size(), size_t(3) );

  // The value lanes must match the double evaluation exactly
  for( size_t i = 0; i < jet_counts.size(); ++i )
    BOOST_CHECK_CLOSE( jet_counts[i].a, dbl_counts[i], 1.0E-9 );

  // Each derivative lane vs a central finite difference.
  //  Step sizes: ~0.1% of each parameter value.
  for( size_t par_index = 0; par_index < 4; ++par_index )
  {
    // The atomic-number lane needs a smaller step: the expected counts have enough
    //  curvature in AN that the h^2 truncation of this check dominates at 1e-3 relative.
    const double rel_step = (par_index == 2) ? 1.0E-4 : 1.0E-3;
    const double step = rel_step * std::max( fabs(params[par_index]), 1.0E-6 );

    vector<double> params_plus = params, params_minus = params;
    params_plus[par_index] += step;
    params_minus[par_index] -= step;

    GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache cache_p, cache_m;
    const vector<double> counts_plus = fcn->expected_peak_counts_imp<double>( params_plus, cache_p );
    const vector<double> counts_minus = fcn->expected_peak_counts_imp<double>( params_minus, cache_m );

    for( size_t i = 0; i < jet_counts.size(); ++i )
    {
      const double numeric = (counts_plus[i] - counts_minus[i]) / (2.0*step);
      const double jet_deriv = jet_counts[i].v[par_index];

      // The age lane uses numeric differencing internally (with its own adaptive step),
      //  so it gets a looser tolerance than the analytically-propagated lanes.
      const double tol = (par_index == 1) ? 0.01 : 1.0E-4;

      BOOST_CHECK_MESSAGE( fabs(jet_deriv - numeric) <= tol*std::max(fabs(numeric), 1.0E-9),
                           "Param " << par_index << ", peak " << i << ": Jet deriv ("
                           << jet_deriv << ") vs numeric (" << numeric << ")" );
    }//for( loop over peaks )
  }//for( loop over parameters )
}//BOOST_AUTO_TEST_CASE( ExpectedPeakCountsJetDerivatives )


/** Runs both fit drivers on a self-attenuating U+Np analyst problem (mass fractions
 in two elements plus a fit thickness - historically the hardest configuration), and
 requires they agree with each other.
 */
BOOST_AUTO_TEST_CASE( CeresVsMinuitDrivers )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  const string n42 = SpecUtils::append_path( g_test_file_dir, "../analysis_tests/u_np_1kg.n42" );
  BOOST_REQUIRE( SpecUtils::is_file(n42) );

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  BOOST_REQUIRE( load_simple_analyst_n42( n42, chi_input ) );

  vector<double> chi2s;
  vector<vector<double>> param_sets;

  for( const bool use_ceres : { false, true } )
  {
    pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                              = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

    auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
    *inputPrams = fcn_pars.second;

    auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
    auto fit_results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
    auto progress_fcn = [](){};
    auto finished_fcn = [](){};

    const auto start_time = std::chrono::steady_clock::now();

    if( use_ceres )
      ShieldingSourceFitCalc::fit_model_ceres( "", fcn_pars.first, inputPrams, progress,
                                               progress_fcn, fit_results, finished_fcn );
    else
      ShieldingSourceFitCalc::fit_model_minuit2( "", fcn_pars.first, inputPrams, progress,
                                                 progress_fcn, fit_results, finished_fcn );

    const auto finish_time = std::chrono::steady_clock::now();
    const double seconds = 1.0E-3 * std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time).count();
    cout << (use_ceres ? "Ceres" : "Minuit") << " driver took " << seconds << " s, "
         << fit_results->num_fcn_calls << " evals" << endl;

    BOOST_REQUIRE( fit_results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final );

    chi2s.push_back( fit_results->chi2 );
    param_sets.push_back( fit_results->paramValues );
  }//for( both drivers )

  BOOST_REQUIRE_EQUAL( chi2s.size(), size_t(2) );
  BOOST_REQUIRE_EQUAL( param_sets[0].size(), param_sets[1].size() );

  cout << "CeresVsMinuitDrivers: Minuit chi2=" << chi2s[0] << ", Ceres chi2=" << chi2s[1] << endl;

  // Both drivers should find (essentially) the same minimum
  BOOST_CHECK_MESSAGE( fabs(chi2s[0] - chi2s[1]) <= 0.01*std::max(chi2s[0], chi2s[1]),
                       "Driver chi2s differ: Minuit=" << chi2s[0] << " vs Ceres=" << chi2s[1] );

  for( size_t i = 0; i < param_sets[0].size(); ++i )
  {
    const double minuit_val = param_sets[0][i];
    const double ceres_val = param_sets[1][i];
    const double scale = std::max( {fabs(minuit_val), fabs(ceres_val), 1.0E-9} );

    BOOST_CHECK_MESSAGE( fabs(minuit_val - ceres_val) <= 0.01*scale,
                         "Parameter " << i << " differs: Minuit=" << minuit_val
                         << " vs Ceres=" << ceres_val );
  }//for( loop over parameters )
}//BOOST_AUTO_TEST_CASE( CeresVsMinuitDrivers )


/** Sanity-checks ShieldingSourceChi2Fcn::computeEffectiveShielding on a volumetric
 (trace-source) and a point-source-with-mass-fractions analyst problem.
 (The quantitative validation of the underlying accumulation is the analytic-slab
 test in test_VolumetricIntegration.cpp.)
 */
BOOST_AUTO_TEST_CASE( ComputeEffectiveShielding )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  for( const string &filename : { string("AEGIS_Eu152_surface_contamination.n42_20230622T113239.276178.n42"),
                                  string("u_np_1kg.n42") } )
  {
    const string n42 = SpecUtils::append_path( g_test_file_dir, "../analysis_tests/" + filename );
    BOOST_REQUIRE( SpecUtils::is_file(n42) );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    BOOST_REQUIRE( load_simple_analyst_n42( n42, chi_input ) );

    pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                              = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

    GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache cache;
    vector<GammaInteractionCalc::EffectiveShieldingInfo> eff_shieldings;
    BOOST_REQUIRE_NO_THROW( eff_shieldings = fcn_pars.first->computeEffectiveShielding(
                                                          fcn_pars.second.Params(), cache ) );

    BOOST_REQUIRE_MESSAGE( !eff_shieldings.empty(), filename << ": no effective shielding entries" );

    size_t num_volumetric = 0;
    for( const GammaInteractionCalc::EffectiveShieldingInfo &info : eff_shieldings )
    {
      num_volumetric += !info.is_point_source;

      BOOST_CHECK_MESSAGE( info.nuclide, filename << ": null nuclide" );
      BOOST_CHECK_MESSAGE( info.energy > 0.0, filename << ": non-positive energy" );
      BOOST_CHECK_MESSAGE( info.effective_ad >= 0.0,
                           filename << ": negative effective AD at " << info.energy << " keV" );

      if( info.effective_ad > 0.0 )
      {
        BOOST_CHECK_MESSAGE( (info.effective_an_mass >= 1.0) && (info.effective_an_mass <= 98.0),
                             filename << ": effective AN (mass) " << info.effective_an_mass
                             << " out of range at " << info.energy << " keV" );
        BOOST_CHECK_MESSAGE( (info.effective_an_xs >= 1.0) && (info.effective_an_xs <= 98.0),
                             filename << ": effective AN (xs) " << info.effective_an_xs
                             << " out of range at " << info.energy << " keV" );
        BOOST_CHECK_MESSAGE( (info.hydrogen_frac_of_ad >= 0.0) && (info.hydrogen_frac_of_ad <= 1.0),
                             filename << ": H fraction " << info.hydrogen_frac_of_ad
                             << " out of range at " << info.energy << " keV" );
      }//if( there was any shielding )
    }//for( loop over entries )

    cout << "ComputeEffectiveShielding '" << filename << "': " << eff_shieldings.size()
         << " entries (" << num_volumetric << " volumetric); first: E="
         << eff_shieldings[0].energy << " keV, AD="
         << eff_shieldings[0].effective_ad/(PhysicalUnits::g/PhysicalUnits::cm2)
         << " g/cm2, AN_mass=" << eff_shieldings[0].effective_an_mass
         << ", AN_xs=" << eff_shieldings[0].effective_an_xs
         << ", fracH=" << eff_shieldings[0].hydrogen_frac_of_ad << endl;

    BOOST_CHECK_MESSAGE( num_volumetric > 0, filename << ": expected volumetric entries" );
  }//for( loop over test files )
}//BOOST_AUTO_TEST_CASE( ComputeEffectiveShielding )


// Diagnostic: Jet-vs-numeric Jacobian of expected_peak_counts_imp on an analyst file
//  (useful when the Ceres driver misbehaves on a particular configuration).
BOOST_AUTO_TEST_CASE( DebugAnalystFileJacobian, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  const string n42 = SpecUtils::append_path( g_test_file_dir, "../analysis_tests/u_np_1kg.n42" );
  BOOST_REQUIRE( SpecUtils::is_file(n42) );

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
  BOOST_REQUIRE( load_simple_analyst_n42( n42, chi_input ) );

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
  const shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> &fcn = fcn_pars.first;
  const vector<double> params = fcn_pars.second.Params();

  cout << "Parameters:" << endl;
  for( size_t i = 0; i < params.size(); ++i )
  {
    const auto &p = fcn_pars.second.Parameters()[i];
    cout << "  [" << i << "] " << p.GetName() << " = " << params[i]
         << (p.IsConst() ? " (const)" : "") << (p.IsFixed() ? " (fixed)" : "")
         << (p.HasLowerLimit() ? (" lo=" + std::to_string(p.LowerLimit())) : string())
         << (p.HasUpperLimit() ? (" up=" + std::to_string(p.UpperLimit())) : string()) << endl;
  }

  using JetN = ceres::Jet<double,16>;
  BOOST_REQUIRE( params.size() <= 16 );

  vector<JetN> jet_params( params.size() );
  for( size_t i = 0; i < params.size(); ++i )
  {
    jet_params[i] = JetN( params[i] );
    jet_params[i].v[i] = 1.0;
  }

  GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache jc, dc;
  const vector<JetN> jets = fcn->expected_peak_counts_imp<JetN>( jet_params, jc );
  const vector<double> dbls = fcn->expected_peak_counts_imp<double>( params, dc );

  BOOST_REQUIRE_EQUAL( jets.size(), dbls.size() );

  for( size_t par = 0; par < params.size(); ++par )
  {
    const auto &p = fcn_pars.second.Parameters()[par];
    if( p.IsConst() || p.IsFixed() || (string(p.GetName()).find("_FIXED") != string::npos) )
      continue;

    const double step = 1.0E-4 * std::max( fabs(params[par]), 1.0E-6 );
    vector<double> pp = params, pm = params;
    pp[par] += step;
    pm[par] -= step;

    GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache cp, cm;
    const vector<double> fp = fcn->expected_peak_counts_imp<double>( pp, cp );
    const vector<double> fm = fcn->expected_peak_counts_imp<double>( pm, cm );

    for( size_t i = 0; i < jets.size(); ++i )
    {
      const double numeric = (fp[i] - fm[i]) / (2.0*step);
      const double jet_d = jets[i].v[par];
      const double denom = std::max( {fabs(numeric), fabs(jet_d), 1.0E-9} );
      if( fabs(jet_d - numeric) > 0.02*denom )
        cout << "  MISMATCH par[" << par << "]=" << p.GetName() << " peak " << i
             << ": jet=" << jet_d << " numeric=" << numeric << endl;
    }
  }//for( loop over params )

  cout << "Jacobian check done" << endl;

  // Now run both fit drivers and compare where they end up
  for( const bool use_ceres : { false, true } )
  {
    auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
    *inputPrams = fcn_pars.second;

    auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
    auto fit_results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
    auto progress_fcn = [](){};
    auto finished_fcn = [](){};

    if( use_ceres )
      ShieldingSourceFitCalc::fit_model_ceres( "", fcn, inputPrams, progress, progress_fcn, fit_results, finished_fcn );
    else
      ShieldingSourceFitCalc::fit_model_minuit2( "", fcn, inputPrams, progress, progress_fcn, fit_results, finished_fcn );

    cout << (use_ceres ? "CERES" : "MINUIT") << " result: chi2=" << fit_results->chi2 << ", params={";
    for( const double v : fit_results->paramValues )
      cout << std::setprecision(6) << v << ", ";
    cout << "}" << endl;

    const double recomputed_chi2 = fcn->DoEval( fit_results->paramValues );
    cout << "  DoEval at solution: " << recomputed_chi2 << endl;
    for( const string &msg : fit_results->errormsgs )
      cout << "  errormsg: " << msg << endl;
  }//for( both drivers )
}//BOOST_AUTO_TEST_CASE( DebugAnalystFileJacobian )


/** A self-attenuating uranium fit started deep in the flat self-attenuation regime (source radius set
 to ~6x its true value, but kept inside the detector distance), combined with a far-low Fe shell.  A
 thick uranium source emits only from a thin surface layer, so the model is nearly insensitive to the
 radius there, and the far U-radius / far-low-Fe combination traps a plain local optimizer near the
 far start with a large chi2.

 This test does two things:
   1) It guards that such a fit COMPLETES rather than crashing.  A far/degenerate start used to drive
      the model past `radius > distance` (NaN residuals) and then segfault in the post-fit covariance
      SVD; both are now handled, and a segfault here would abort the test process.
   2) It verifies the far-start recovery: run_ceres_solve_with_recovery() detects the disastrous chi2
      and re-solves from a bounded multi-start over the fitted shielding dimensions, so the fit now
      climbs back to the true radius.  (This was an EXPECTED FAILURE before that recovery landed.)
 */
BOOST_AUTO_TEST_CASE( SelfAttenUraniumFarStartStuck )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  const string n42 = SpecUtils::append_path( g_test_file_dir, "uranium_40pct_selfatten_HPGe_15cm.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(n42), "missing fixture: " << n42 );

  SpecMeas specfile;
  BOOST_REQUIRE( specfile.load_N42_file( n42 ) );

  shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks;
  set<int> fg_samples;
  for( const set<int> &samples : specfile.sampleNumsWithPeaks() )
  {
    shared_ptr<const deque<shared_ptr<const PeakDef>>> p = specfile.peaks( samples );
    if( p && !p->empty() ){ peaks = p; fg_samples = samples; break; }
  }
  BOOST_REQUIRE( peaks && !peaks->empty() );

  shared_ptr<const DetectorPeakResponse> det = specfile.detector();
  BOOST_REQUIRE( det && det->isValid() );

  rapidxml::xml_document<char> *model_xml = specfile.shieldingSourceModel();
  BOOST_REQUIRE( model_xml && model_xml->first_node() );
  GammaInteractionCalc::ShieldSourceConfig config;
  BOOST_REQUIRE_NO_THROW( config.deSerialize( model_xml->first_node() ) );
  BOOST_REQUIRE( config.geometry == GammaInteractionCalc::GeometryType::Spherical );

  const vector<string> &dn = specfile.detector_names();
  shared_ptr<const SpecUtils::Measurement> fg = specfile.measurement( *fg_samples.begin(),
                                                                      dn.empty() ? string() : dn[0] );
  if( !fg )
    fg = specfile.measurements()[0];

  // The uranium (self-attenuating source) shell; its saved radius is the true/converged answer.
  ShieldingSourceFitCalc::ShieldingInfo *u_shell = nullptr;
  for( ShieldingSourceFitCalc::ShieldingInfo &sh : config.shieldings )
    if( sh.m_material && SpecUtils::icontains(sh.m_material->name, "uranium") )
      u_shell = &sh;
  BOOST_REQUIRE_MESSAGE( u_shell, "fixture changed: no uranium shell found" );
  BOOST_REQUIRE_MESSAGE( u_shell->m_fitDimensions[0], "fixture changed: uranium radius is not fit" );

  const double true_U_radius = u_shell->m_dimensions[0];
  BOOST_REQUIRE( true_U_radius > 0.0 );

  // Reproduce the documented far-start failure: perturb every fitted dimension by a fixed-seed random
  //  0.1-10x factor (deterministic std::mt19937, so portable).  For this fixture this drives the U
  //  source radius to ~6x the truth (deep in the flat self-attenuation regime) while also driving the
  //  thin Fe shell far low - the combination the optimizer cannot recover from.  (Perturbing the U
  //  radius alone, with a correct Fe, does recover - so both must be off to exercise the failure.)
  std::mt19937 rng( 20240601u + 1000u );
  std::uniform_real_distribution<double> logfac( std::log(0.1), std::log(10.0) );
  for( ShieldingSourceFitCalc::ShieldingInfo &sh : config.shieldings )
    for( int d = 0; d < 3; ++d )
      if( sh.m_fitDimensions[d] && (sh.m_dimensions[d] > 0.0) )
        sh.m_dimensions[d] *= std::exp( logfac(rng) );
  const double far_start = u_shell->m_dimensions[0];

  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput ci;
  ci.config = config;
  ci.detector = det;
  ci.foreground = fg;
  ci.background = nullptr;
  ci.foreground_peaks.assign( peaks->begin(), peaks->end() );
  ci.background_peaks = nullptr;

  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars
                            = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( ci );
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;
  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  auto pf = [](){};
  auto ff = [](){};

  // (1) Must complete without crashing / throwing (the regression guard - a segfault would abort here).
  BOOST_REQUIRE_NO_THROW( ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams,
                                                             progress, pf, results, ff ) );
  BOOST_CHECK( std::isfinite( results->chi2 ) );

  double fitted_U_radius = -1.0;
  for( const ShieldingSourceFitCalc::FitShieldingInfo &sh : results->final_shieldings )
    if( sh.m_material && SpecUtils::icontains(sh.m_material->name, "uranium") )
      fitted_U_radius = sh.m_dimensions[0];

  BOOST_TEST_MESSAGE( "Self-atten U far-start (started " << far_start/PhysicalUnits::cm << " cm): fitted "
                      << fitted_U_radius/PhysicalUnits::cm << " cm vs true "
                      << true_U_radius/PhysicalUnits::cm << " cm" );

  // (2) The far-start recovery must climb back to the true radius from the far start.
  BOOST_CHECK_MESSAGE( std::fabs(fitted_U_radius - true_U_radius) <= 0.10*true_U_radius,
                       "far-start self-attenuating fit did not recover the true U radius: got "
                       << fitted_U_radius/PhysicalUnits::cm << " cm, true "
                       << true_U_radius/PhysicalUnits::cm << " cm" );
}//BOOST_AUTO_TEST_CASE( SelfAttenUraniumFarStartStuck )


BOOST_AUTO_TEST_CASE( FitAnalystShieldingSourcecases )
{
  set_data_dir();

  cout << "Starting FitAnalystShieldingSourcecases" << endl;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
 
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const std::shared_ptr<const Material> iron = matdb->material( "Fe (iron)" );
  BOOST_REQUIRE( iron );

  const string base_dir = SpecUtils::append_path(g_test_file_dir, "../analysis_tests");
  BOOST_REQUIRE( SpecUtils::is_directory(base_dir) );
  
  const vector<string> files = SpecUtils::recursive_ls( base_dir, ".n42" );
  
  // Make sure we found some files
  BOOST_REQUIRE( files.size() >= 1 );
  
  for( const string n42_filename : files )
  {
    SpecMeas specfile;
    const bool loaded = specfile.load_N42_file( n42_filename );
    BOOST_CHECK_MESSAGE( loaded, "Analyst file '" << n42_filename << "' couldnt be loaded - skipping testing." );
    if( !loaded )
      continue;
    
    rapidxml::xml_document<char> *model_xml = specfile.shieldingSourceModel();
    BOOST_CHECK_MESSAGE( model_xml, "Analyst file '" << n42_filename << "' doesnt have a ShieldingSourceModel - skipping test." );
    if( !model_xml )
      continue;
    
    rapidxml::xml_node<char> *base_node = model_xml->first_node();
    BOOST_CHECK_MESSAGE( base_node, "Analyst file '" << n42_filename << "' has invalid ShieldingSourceModel - skipping test." );
    if( !base_node )
      continue;
    
    // We will go through and try to figure out what the "foreground" should be.
    //  Note: we could probably just call `specfile.displayedSampleNumbers()`, but we'll
    //        be hard on ourselves
    const vector<string> &detectors = specfile.detector_names();
    const set<set<int>> samplesWithPeaks = specfile.sampleNumsWithPeaks();
    size_t num_foreground = 0, num_background = 0;
    set<int> foreground_samples, background_samples;
    
    for( const set<int> &samples : samplesWithPeaks )
    {
      bool classified_sample_nums = false;
      for( const int sample : samples )
      {
        for( const string &det : detectors )
        {
          auto m = specfile.measurement( sample, det );
          if( !m )
            continue;
          
          switch( m->source_type() )
          {
            case SpecUtils::SourceType::IntrinsicActivity:
            case SpecUtils::SourceType::Calibration:
              break;
              
            case SpecUtils::SourceType::Background:
              classified_sample_nums = true;
              num_background += 1;
              background_samples = samples;
              break;
            
            case SpecUtils::SourceType::Foreground:
            case SpecUtils::SourceType::Unknown:
              classified_sample_nums = true;
              num_foreground += 1;
              foreground_samples = samples;
              break;
          }//switch( m->source_type() )
          
          if( classified_sample_nums )
            break;
        }//for( const string &det : detectors )
        
        if( classified_sample_nums )
          break;
      }//for( const int sample : samples )
      
      assert( !classified_sample_nums || ((num_foreground >= 1) || (num_background >= 1)) );
    }//for( const set<int> &samples : samplesWithPeaks )
    
    if( (num_background == 1) && (num_foreground == 0) )
    {
      std::swap( num_background, num_foreground );
      std::swap( foreground_samples, background_samples );
    }
    
    // Check that it is unambiguous what foreground the model file is to
    BOOST_CHECK_MESSAGE( num_foreground == 1,
                        "Analyst file '" << n42_filename << "' does not have exactly one foreground"
                        " - skipping test." );
    if( num_foreground != 1 )
      continue;
    
    const shared_ptr<const SpecUtils::Measurement> foreground
                 = ((foreground_samples.size() == 1) && (detectors.size() == 1))
                   ? specfile.measurement( *begin(foreground_samples), detectors.front() )
                   : specfile.sum_measurements( foreground_samples, detectors,  nullptr );
    
    BOOST_CHECK_MESSAGE( foreground,
                        "Analyst file '" << n42_filename << "' failed to extract foreground"
                        " - skipping test." );
    if( !foreground )
      continue;
    
    
    // Find foreground that has act/shielding fit associated
    shared_ptr<const deque<shared_ptr<const PeakDef>>> foreground_peaks = specfile.peaks( foreground_samples );
    BOOST_CHECK_MESSAGE( foreground_peaks && (foreground_peaks->size() >= 1),
                         "Analyst file '" << n42_filename << "' didnt have peaks for the"
                         " identified foreground - skipping test." );
    if( !foreground_peaks || foreground_peaks->empty() )
      continue;
        
    shared_ptr<const DetectorPeakResponse> detector = specfile.detector();
    BOOST_CHECK_MESSAGE( detector && detector->isValid(),
                        "Analyst file '" << n42_filename << "' didnt have a DRF - skipping test." );
    
    if( !detector || !detector->isValid() )
      continue;
    
    
    //string xml;
    //rapidxml::print(std::back_inserter(xml), *model_xml, 0);
    //cout << "XML: " << xml << endl;
    
    double distance;
    std::string diststr;
    try
    {
      const rapidxml::xml_node<char> *dist_node = base_node->first_node( "Distance" );
      diststr = SpecUtils::xml_value_str( dist_node );
      distance = PhysicalUnits::stringToDistance( diststr );
    }catch( std::exception &e )
    {
      BOOST_CHECK_MESSAGE( false,
                          "Analyst file '" << n42_filename << "' <Distance> node invalid: "
                          << e.what() << " - skipping test." );
      continue;
    }//try / catch to get distance
    
    GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
    
    const rapidxml::xml_node<char> *geom_node = base_node->first_node( "Geometry" );
    const string geomstr = SpecUtils::xml_value_str( geom_node );
    
    for( GammaInteractionCalc::GeometryType type = GammaInteractionCalc::GeometryType(0);
        type != GammaInteractionCalc::GeometryType::NumGeometryType;
        type = GammaInteractionCalc::GeometryType(static_cast<int>(type) + 1) )
    {
      if( SpecUtils::iequals_ascii(geomstr, GammaInteractionCalc::to_str(type)) )
      {
        geometry = type;
        break;
      }
    }
    
    BOOST_CHECK_MESSAGE( geometry != GammaInteractionCalc::GeometryType::NumGeometryType,
                        "Analyst file '" << n42_filename << "' <Geometry> node missing or invalid"
                        " - skipping test.");
    if( geometry == GammaInteractionCalc::GeometryType::NumGeometryType )
      continue;
    
    
    vector<ShieldingSourceFitCalc::ShieldingInfo> shield_definitions;
    vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
    ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
    try
    {
      options.deSerialize( base_node );
    }catch( std::exception &e )
    {
      BOOST_CHECK_MESSAGE( false,
                          "Analyst file '" << n42_filename << "' Options invalid: "
                          << e.what()
                          << " - skipping test.");
      continue;
    }

    
    BOOST_CHECK_MESSAGE( !options.background_peak_subtract || (num_background <= 1),
                        "Analyst file '" << n42_filename << "' Options indicated"
                        " background-subtract, but background is ambiguous. - skipping test.");
    if( options.background_peak_subtract && (num_background > 1) )
      continue;
    
    // It looks like background subtract _could_ be set to true, even if we didnt have any
    //  background peaks for the fit - at least if things are defined not using the GUI (the
    //  GUI enforces background peaks must be present to have this option checked) so I guess
    //  we'll just be a bit loose around this for now, but we probably dont have to be.
    shared_ptr<const SpecUtils::Measurement> background;
    shared_ptr<const deque<shared_ptr<const PeakDef>>> background_peaks;
    if( options.background_peak_subtract && (num_background == 1) )
    {
      background = specfile.sum_measurements( background_samples, detectors, nullptr );
      BOOST_CHECK_MESSAGE( background && (background->num_gamma_channels() > 16),
                          "Analyst file '" << n42_filename << "' Options indicated"
                          " background-subtract, but background couldnt be extracted - skipping test.");
      if( !background || (background->num_gamma_channels() <= 16) )
        continue;
      
      background_peaks = specfile.peaks( background_samples );
      
      if( background_peaks && background_peaks->empty() )
      {
        background = nullptr;
        background_peaks.reset();
      }
      
      BOOST_WARN_MESSAGE( background_peaks, "Analyst file '" + n42_filename + "' Options indicated"
                   " background-subtract, but no background peaks found - continuing anyway." );
    }//if( options.background_peak_subtract && (num_background == 1) )
    
    
    const rapidxml::xml_node<char> *shieldings_node = base_node->first_node( "Shieldings" );
    BOOST_REQUIRE( shieldings_node );
    
    bool success_parsing = true;
    XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
    {
      try
      {
        ShieldingSourceFitCalc::ShieldingInfo info;
        info.deSerialize( shield_node );
        shield_definitions.push_back( std::move(info) );
      }catch( std::exception &e )
      {
        success_parsing = false;
        BOOST_CHECK_MESSAGE( false,
                            "Analyst file '" << n42_filename << "' <Shielding> node invalid: "
                             << e.what() << " - skipping test.");
      }
    }//XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
    
    if( !success_parsing )
      continue;
    
    
    const rapidxml::xml_node<char> *srcs_node = base_node->first_node( "Nuclides" );
    BOOST_CHECK_MESSAGE( srcs_node,
                        "Analyst file '" << n42_filename << "' missing <Nuclides> node: "
                         " - skipping test.");
    
    XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
    {
      try
      {
        ShieldingSourceFitCalc::SourceFitDef info;
        info.deSerialize( src_node );
        src_definitions.push_back( std::move(info) );
      }catch( std::exception &e )
      {
        success_parsing = false;
        BOOST_CHECK_MESSAGE( false,
                            "Analyst file '" << n42_filename << "' <Nuclide> node invalid: "
                             << e.what() << " - skipping test.");
      }// try / catch
    }//XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
    
    if( !success_parsing )
      continue;
    
    BOOST_CHECK_MESSAGE( src_definitions.size() >= 1,
                        "Analyst file '" << n42_filename << "' no sources defined." );
    
    // Check that the number of fit nuclides in `src_definitions` matches the sources in the peaks
    set<const SandiaDecay::Nuclide *> peak_fit_nucs;
    for( const auto &p : *foreground_peaks )
    {
      if( p->parentNuclide() && p->useForShieldingSourceFit())
        peak_fit_nucs.insert( p->parentNuclide() );
    }
    
    BOOST_CHECK_MESSAGE( peak_fit_nucs.size() == src_definitions.size(),
                        "Analyst file '" << n42_filename << "' has mismatching number of sources in peaks ("
                        << peak_fit_nucs.size() << "), vs src_definitions (" << src_definitions.size() << ")" );
    
    
    GammaInteractionCalc::ShieldSourceConfig parsed_config;
    BOOST_CHECK_NO_THROW( parsed_config.deSerialize( base_node ) );

    BOOST_CHECK_SMALL( fabs(parsed_config.distance - distance), 1e-9 * (std::max)(1.0, fabs(distance)) );
    BOOST_CHECK_EQUAL( static_cast<int>(parsed_config.geometry), static_cast<int>(geometry) );
    BOOST_CHECK_EQUAL( parsed_config.shieldings.size(), shield_definitions.size() );
    std::vector<ShieldingSourceFitCalc::SourceFitDef> parsed_source_defs = parsed_config.sources;
    BOOST_CHECK_EQUAL( parsed_source_defs.size(), src_definitions.size() );

#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    for( size_t i = 0; i < shield_definitions.size(); ++i )
      BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingInfo::equalEnough( shield_definitions[i], parsed_config.shieldings[i] ) );
  for( size_t i = 0; i < src_definitions.size(); ++i )
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::SourceFitDef::equalEnough( src_definitions[i], parsed_source_defs[i] ) );
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( options, parsed_config.options ) );
#endif

    rapidxml::xml_document<char> roundtrip_doc;
    rapidxml::xml_node<char> *round_root = roundtrip_doc.allocate_node( rapidxml::node_element, "ShieldingSourceFit" );
    roundtrip_doc.append_node( round_root );
    BOOST_CHECK_NO_THROW( parsed_config.serialize( round_root ) );

    GammaInteractionCalc::ShieldSourceConfig reparsed_config;
    BOOST_CHECK_NO_THROW( reparsed_config.deSerialize( round_root ) );

    BOOST_CHECK_SMALL( fabs(reparsed_config.distance - parsed_config.distance), 1e-9 * (std::max)(1.0, fabs(parsed_config.distance)) );
    BOOST_CHECK_EQUAL( static_cast<int>(reparsed_config.geometry), static_cast<int>(parsed_config.geometry) );
    BOOST_CHECK_EQUAL( reparsed_config.shieldings.size(), parsed_config.shieldings.size() );
  std::vector<ShieldingSourceFitCalc::SourceFitDef> reparsed_source_defs = reparsed_config.sources;
  BOOST_CHECK_EQUAL( reparsed_source_defs.size(), parsed_source_defs.size() );

#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    for( size_t i = 0; i < reparsed_config.shieldings.size(); ++i )
      BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingInfo::equalEnough( reparsed_config.shieldings[i], parsed_config.shieldings[i] ) );
  for( size_t i = 0; i < reparsed_source_defs.size(); ++i )
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::SourceFitDef::equalEnough( reparsed_source_defs[i], parsed_source_defs[i] ) );
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( reparsed_config.options, parsed_config.options ) );
#endif

    set_fit_quantities_to_default_values( shield_definitions, src_definitions );
    
    // We have all the parts, lets do the computation:
    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    chi_input.config.distance = distance;
    chi_input.config.geometry = geometry;
    chi_input.config.shieldings = shield_definitions;
    chi_input.config.sources = src_definitions;
    chi_input.config.options = options;
    chi_input.detector = detector;
    chi_input.foreground = foreground;
    chi_input.background = background;
    chi_input.foreground_peaks.assign( foreground_peaks->begin(), foreground_peaks->end() );
    chi_input.background_peaks = background_peaks;
    pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
    GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
    
    auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
    *inputPrams = fcn_pars.second;
    
    auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
    auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
    
    
    auto progress_fcn = [progress](){
      // Probably wont get called, since its a simple fit
    };
    
    bool finished_fit_called = false;
    auto finished_fcn = [results, &finished_fit_called](){
      finished_fit_called = true;
    };
    
    ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, results, finished_fcn );
    
    BOOST_CHECK( finished_fit_called );
    
    BOOST_CHECK_MESSAGE( results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final,
                        "Analyst file '" << n42_filename << "': Fit status was not successful" );
    
    BOOST_CHECK_MESSAGE( results->fit_src_info.size() == src_definitions.size(),
                        "Analyst file '" << n42_filename << "': Fit sources size not equal to input." );
    
    tuple<bool,int,int,vector<string>> test_results = test_fit_against_truth( *results, n42_filename );
    const bool successful = get<0>(test_results);
    const int numCorrect = get<1>(test_results);
    const int numTested = get<2>(test_results);
    const vector<string> &textInfoLines = get<3>(test_results);
    
    cout << (successful ? "Successfully" : "Unsuccessfully") << " checked against truth values with "
    << numCorrect << " tests passing, out of " << numTested << ".\nNotes:" << endl;
    for( const auto &msg : textInfoLines )
      cout << "\t" << msg << endl;
  }//for( const string n42_filename : files )
}//BOOST_AUTO_TEST_CASE( FitAnalystShieldingSourcecases )


BOOST_AUTO_TEST_CASE( ShieldingSourceDisplayGuiRoundTrip )
{
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const string base_dir = SpecUtils::append_path(g_test_file_dir, "../analysis_tests");
  BOOST_REQUIRE( SpecUtils::is_directory(base_dir) );
  
  const vector<string> files = SpecUtils::recursive_ls( base_dir, ".n42" );
  
  // Make sure we found some files
  BOOST_REQUIRE( files.size() >= 1 );
  
  for( const string n42_filename : files )
  {
    // Create a fresh fixture for each file to ensure clean state
    InterSpecTestFixture fixture;
    InterSpec *m_interspec = fixture.m_interspec;
    
    SpecMeas specfile;
    const bool loaded = specfile.load_N42_file( n42_filename );
    BOOST_CHECK_MESSAGE( loaded, "Analyst file '" << n42_filename << "' couldnt be loaded - skipping testing." );
    if( !loaded )
      continue;
    
    rapidxml::xml_document<char> *model_xml = specfile.shieldingSourceModel();
    BOOST_CHECK_MESSAGE( model_xml, "Analyst file '" << n42_filename << "' doesnt have a ShieldingSourceModel - skipping test." );
    if( !model_xml )
      continue;
    
    rapidxml::xml_node<char> *base_node = model_xml->first_node();
    BOOST_CHECK_MESSAGE( base_node, "Analyst file '" << n42_filename << "' has invalid ShieldingSourceModel - skipping test." );
    if( !base_node )
      continue;
    
    // Parse original config from XML
    GammaInteractionCalc::ShieldSourceConfig original_config;
    BOOST_CHECK_NO_THROW( original_config.deSerialize( base_node ) );
    
    // Load the file into InterSpec
    m_interspec->userOpenFileFromFilesystem( n42_filename, n42_filename );
    
    shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
    BOOST_CHECK_MESSAGE( meas, "Analyst file '" << n42_filename << "' failed to load into InterSpec - skipping test." );
    if( !meas )
      continue;
    
    // Verify the measurement has the shielding source model
    rapidxml::xml_document<char> *loaded_model_xml = meas->shieldingSourceModel();
    BOOST_CHECK_MESSAGE( loaded_model_xml && loaded_model_xml->first_node(),
                        "Analyst file '" << n42_filename << "' measurement doesn't have shielding source model after loading - skipping test." );
    if( !loaded_model_xml || !loaded_model_xml->first_node() )
      continue;
    
    
    // Create the GUI (it will auto-deserialize from the measurement's shielding source model)
    ShieldingSourceDisplay *gui = m_interspec->shieldingSourceFit();
    BOOST_CHECK_MESSAGE( gui, "Analyst file '" << n42_filename << "' failed to create ShieldingSourceDisplay GUI - skipping test." );
    if( !gui )
      continue;
    
    // Serialize the GUI state back using the new serialize() method
    // This tests that the new serialize() overload works correctly
    ShieldingSourceDisplay::ShieldingSourceDisplayState gui_state;
    BOOST_CHECK_NO_THROW( gui_state = gui->serialize() );
    
    BOOST_CHECK_MESSAGE( gui_state.config, "Analyst file '" << n42_filename << "' GUI state has no config - skipping test." );
    if( !gui_state.config )
      continue;
    
    const GammaInteractionCalc::ShieldSourceConfig &gui_config = *gui_state.config;
    
    // Compare the configs - note that the GUI may have auto-deserialized, so we check if it matches
    // If the GUI has sources/shieldings, compare them; otherwise skip this file
    if( gui_config.sources.empty() && gui_config.shieldings.empty() )
    {
      BOOST_WARN_MESSAGE( false, "Analyst file '" << n42_filename << "' GUI state is empty - GUI may not have auto-deserialized. Skipping comparison." );
      continue;
    }
    
    // Compare basic properties with detailed error messages
    const double distance_diff = fabs(gui_config.distance - original_config.distance);
    const double distance_tolerance = 1e-6 * (std::max)(1.0, fabs(original_config.distance));
    BOOST_CHECK_MESSAGE( distance_diff < distance_tolerance,
                        "Distance mismatch for '" << n42_filename << "': original=" 
                        << PhysicalUnits::printToBestLengthUnits(original_config.distance, 6) << " (" << original_config.distance << ")"
                        << ", GUI=" << PhysicalUnits::printToBestLengthUnits(gui_config.distance, 6) << " (" << gui_config.distance << ")"
                        << ", diff=" << distance_diff << ", tolerance=" << distance_tolerance );
    
    BOOST_CHECK_MESSAGE( static_cast<int>(gui_config.geometry) == static_cast<int>(original_config.geometry),
                        "Geometry mismatch for '" << n42_filename << "': original=" 
                        << static_cast<int>(original_config.geometry) << " (" << GammaInteractionCalc::to_str(original_config.geometry) << ")"
                        << ", GUI=" << static_cast<int>(gui_config.geometry) << " (" << GammaInteractionCalc::to_str(gui_config.geometry) << ")" );
    
    // Compare shieldings
    BOOST_CHECK_MESSAGE( gui_config.shieldings.size() == original_config.shieldings.size(),
                        "Shieldings count mismatch for '" << n42_filename << "': original=" 
                        << original_config.shieldings.size() << ", GUI=" << gui_config.shieldings.size() );
    
    if( gui_config.shieldings.size() == original_config.shieldings.size() )
    {
      for( size_t i = 0; i < original_config.shieldings.size(); ++i )
      {
        try
        {
          ShieldingSourceFitCalc::ShieldingInfo::equalEnough( original_config.shieldings[i], gui_config.shieldings[i] );
        }catch( std::exception &e )
        {
          cout << "Shielding " << i << " mismatch for '" << n42_filename << "': " << e.what() << endl;
          const auto &orig_sh = original_config.shieldings[i];
          const auto &gui_sh = gui_config.shieldings[i];
          cout << "  Original: geometry=" << static_cast<int>(orig_sh.m_geometry) 
               << ", isGeneric=" << orig_sh.m_isGenericMaterial
               << ", material=" << (orig_sh.m_material ? orig_sh.m_material->name : "null") << endl;
          cout << "  GUI: geometry=" << static_cast<int>(gui_sh.m_geometry)
               << ", isGeneric=" << gui_sh.m_isGenericMaterial
               << ", material=" << (gui_sh.m_material ? gui_sh.m_material->name : "null") << endl;
          throw;
        }
      }
    }
    
    // Compare sources
    BOOST_CHECK_MESSAGE( gui_config.sources.size() == original_config.sources.size(),
                        "Sources count mismatch for '" << n42_filename << "': original=" 
                        << original_config.sources.size() << ", GUI=" << gui_config.sources.size() );
    
    if( gui_config.sources.size() == original_config.sources.size() )
    {
      for( size_t i = 0; i < original_config.sources.size(); ++i )
      {
        const auto &orig_src = original_config.sources[i];
        const auto &gui_src = gui_config.sources[i];
        
        // Compare nuclide
        BOOST_CHECK_MESSAGE( orig_src.nuclide == gui_src.nuclide,
                            "Source " << i << " nuclide mismatch for '" << n42_filename << "': original="
                            << (orig_src.nuclide ? orig_src.nuclide->symbol : "null")
                            << ", GUI=" << (gui_src.nuclide ? gui_src.nuclide->symbol : "null") );
        
        // Compare sourceType
        BOOST_CHECK_MESSAGE( static_cast<int>(orig_src.sourceType) == static_cast<int>(gui_src.sourceType),
                            "Source " << i << " sourceType mismatch for '" << n42_filename << "': original="
                            << static_cast<int>(orig_src.sourceType) << ", GUI=" << static_cast<int>(gui_src.sourceType) );
        
        // Compare fitActivity
        BOOST_CHECK_MESSAGE( orig_src.fitActivity == gui_src.fitActivity,
                            "Source " << i << " fitActivity mismatch for '" << n42_filename << "': original="
                            << orig_src.fitActivity << ", GUI=" << gui_src.fitActivity );
        
        // Compare activity - allow 1% difference for self-attenuating sources
        const double act_diff = fabs(orig_src.activity - gui_src.activity);
        const double max_act = (std::max)(fabs(orig_src.activity), fabs(gui_src.activity));
        const bool is_self_attenuating = (orig_src.sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
        const double activity_tolerance = is_self_attenuating ? 0.01 * max_act : 1.0E-6 * max_act;
        
        BOOST_CHECK_MESSAGE( act_diff <= activity_tolerance,
                            "Source " << i << " activity mismatch for '" << n42_filename << "': original="
                            << PhysicalUnits::printToBestActivityUnits(orig_src.activity) << " (" << orig_src.activity << ")"
                            << ", GUI=" << PhysicalUnits::printToBestActivityUnits(gui_src.activity) << " (" << gui_src.activity << ")"
                            << ", diff=" << act_diff << ", tolerance=" << activity_tolerance
                            << (is_self_attenuating ? " (self-attenuating, 1% allowed)" : "") );
        
        // Compare age
        const double age_diff = fabs(orig_src.age - gui_src.age);
        const double max_age = (std::max)(fabs(orig_src.age), fabs(gui_src.age));
        const double age_tolerance = 1.0E-6 * (std::max)(1.0, max_age);
        BOOST_CHECK_MESSAGE( age_diff <= age_tolerance,
                            "Source " << i << " age mismatch for '" << n42_filename << "': original="
                            << orig_src.age << ", GUI=" << gui_src.age
                            << ", diff=" << age_diff << ", tolerance=" << age_tolerance );
        
        // Compare fitAge
        BOOST_CHECK_MESSAGE( orig_src.fitAge == gui_src.fitAge,
                            "Source " << i << " fitAge mismatch for '" << n42_filename << "': original="
                            << orig_src.fitAge << ", GUI=" << gui_src.fitAge );
        
        // Compare ageDefiningNuc
        BOOST_CHECK_MESSAGE( orig_src.ageDefiningNuc == gui_src.ageDefiningNuc,
                            "Source " << i << " ageDefiningNuc mismatch for '" << n42_filename << "': original="
                            << (orig_src.ageDefiningNuc ? orig_src.ageDefiningNuc->symbol : "null")
                            << ", GUI=" << (gui_src.ageDefiningNuc ? gui_src.ageDefiningNuc->symbol : "null") );
        
        // Note: activityUncertainty and ageUncertainty are not preserved by GUI, so we don't compare them
      }
    }
    
    BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::ShieldingSourceFitOptions::equalEnough( original_config.options, gui_config.options ) );
  }//for( const string n42_filename : files )
}//BOOST_FIXTURE_TEST_CASE( ShieldingSourceDisplayGuiRoundTrip )
