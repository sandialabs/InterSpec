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

#include <string>
#include <iostream>

#include <Wt/Utils>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ShieldingSourceFitCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

//Roots Minuit2 includes
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"


#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"


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


BOOST_AUTO_TEST_CASE( ShieldingInfoUri )
{
  using namespace ShieldingSourceFitCalc;
  
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  const Material * const iron = matdb.material("Fe (iron)");
  BOOST_REQUIRE( iron );
  
  const Material * const uranium = matdb.material("U (uranium)");
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
      BOOST_REQUIRE_NO_THROW( from_uri.handleAppUrl( uri, &matdb ) );
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
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node(), &matdb ) );
      BOOST_REQUIRE_NO_THROW( ShieldingInfo::equalEnough( info, from_xml ) );
    }
    
    
    // Now modify `info` to include trace and self-attenuating sources
    info.m_material = make_shared<Material>( *uranium );
    
    
    // Self-atten source stuff
    info.m_fitMassFrac = true;
    
    {// Begin add self-atten source
      map<const SandiaDecay::Nuclide *,double> &el_to_nuc_frac = info.m_nuclideFractions;
      
      const SandiaDecay::Nuclide *nuc = db->nuclide( "U-238" );
      BOOST_REQUIRE( nuc );
      el_to_nuc_frac[nuc] = 0.8;
      
      nuc = db->nuclide( "U-235" );
      BOOST_REQUIRE( nuc );
      el_to_nuc_frac[nuc] = 0.15;
      
      nuc = db->nuclide( "U-232" );
      BOOST_REQUIRE( nuc );
      el_to_nuc_frac[nuc] = 2.5E-8;
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
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node(), &matdb ) );
      BOOST_REQUIRE_NO_THROW( ShieldingInfo::equalEnough( info, from_xml ) );
    }
  }// End test ShieldingInfo encode/decode from URI
  
}//BOOST_AUTO_TEST_CASE( ShieldingInfoUri )


BOOST_AUTO_TEST_CASE( IsoFitStructSerialization )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  
  ShieldingSourceFitCalc::IsoFitStruct test;
  
  test.nuclide = db->nuclide( "U238" );
  test.activity = 1.1*PhysicalUnits::curie;
  test.fitActivity = true;
  test.age = 20*PhysicalUnits::year;
  test.fitAge = false;
  test.ageDefiningNuc = db->nuclide( "U235" );
  test.sourceType = ShieldingSourceFitCalc::ModelSourceType::Intrinsic;
  test.numProgenyPeaksSelected = 2;
  test.ageIsFittable = !PeakDef::ageFitNotAllowed( test.ageDefiningNuc ) ;
  test.activityUncertainty = 0.1*test.activity;
  test.ageUncertainty = 0.0;
  
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  test.truthActivity = 1.0*PhysicalUnits::curie;
  test.truthActivityTolerance = 0.15*PhysicalUnits::curie;
  test.truthAge = 21.5*PhysicalUnits::year;
  test.truthAgeTolerance = 1.05*PhysicalUnits::year;
#endif
  
  rapidxml::xml_document<char> doc;
  BOOST_REQUIRE_NO_THROW( test.serialize( &doc ) );
  
  ShieldingSourceFitCalc::IsoFitStruct from_xml;
  BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
  BOOST_CHECK_NO_THROW( ShieldingSourceFitCalc::IsoFitStruct::equalEnough( test, from_xml ) );
}//BOOST_AUTO_TEST_CASE( IsoFitStructSerialization )


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
  detector->fromExpOfLogPowerSeriesAbsEff( {0.0f, 0.0f}, {}, distance,
                                          5*PhysicalUnits::cm, PhysicalUnits::keV, 0, 3000*PhysicalUnits::keV );
  
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
  
  
  GammaInteractionCalc::ShieldingSourceFitOptions options;
  options.attenuate_for_air = false;
  
  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
                GammaInteractionCalc::ShieldingSourceChi2Fcn::create( distance, geometry,
                                  shieldings, src_definitions, detector,
                                  foreground, nullptr, foreground_peaks, nullptr, options );
 
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
  
  const ShieldingSourceFitCalc::IsoFitStruct &fit_src_info = results->fit_src_info[0];
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
