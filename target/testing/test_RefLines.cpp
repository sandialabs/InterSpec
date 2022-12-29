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
#include <map>
#include <deque>
#include <tuple>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <cstdlib>


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RefLines_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;
using namespace boost::unit_test;

namespace SetDataDirs
{
// The "InterSpec/data" directory that contains all the cross-sections, reactions, etc
std::string sm_static_data_dir = "";

// The sandia.decay.xml file in InterSpec/data is usually the minimized (~6MB) version of the file
//  without coincidences and stuff - we want to use the full ~30MB version of the file that is in
//  the SandiaDecay repository.
std::string sm_sandia_decay_file = "";

void set_data_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;
  
  s_have_set = true;
  
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  std::string datadir, test_file_dir;
  
  for( int i = 1; i < argc; ++i )
  {
    const std::string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
    
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      test_file_dir = arg.substr( 14 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( test_file_dir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.reactiongamma.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  const std::string sandia_reaction_file = SpecUtils::append_path(datadir, "sandia.reactiongamma.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_reaction_file ), "sandia.reactiongamma.xml not at '" << sandia_reaction_file << "'" );
  
  // Set the static data directory so we have cross-sections, reactions, and all that
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  
  // Now find the full-version of sandia.decay.xml
  
  // What we actually want is the InterSpec code base-path - we'll search around a little for it
  const string decay_xml = "external_libs/SandiaDecay/sandia.decay.xml";
  string potential_code_dirs[] = {
    SpecUtils::append_path( test_file_dir, "../../../" ),
    "..", "../..", "../../..", "/Users/wcjohns/rad_ana/InterSpec/"
  };
  
  for( const auto &d : potential_code_dirs )
  {
    const string candidate = SpecUtils::append_path( d, decay_xml );
    if( SpecUtils::is_file(candidate) )
    {
      sm_sandia_decay_file = candidate;
      break;
    }
  }
  
  BOOST_REQUIRE_MESSAGE( !sm_sandia_decay_file.empty() && SpecUtils::is_file(sm_sandia_decay_file),
                        "Error finding the full version of sandia.decay.xml" );
  
  BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( sm_sandia_decay_file ) );
  
  sm_static_data_dir = datadir;
  
  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing full SandiaDecayDataBase" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE_MESSAGE( u238, "Full SandiaDecayDataBase empty?" );
  
  // Now make sure we have the full database with x-rays and coincidences
  bool has_coincidences = false, has_xray = false;
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( u238, 1.0E-3 * SandiaDecay::curie );
  
  for( const SandiaDecay::NuclideActivityPair &nap : mixture.activity( 20*SandiaDecay::year ) )
  {
    for( const SandiaDecay::Transition *transition : nap.nuclide->decaysToChildren )
    {
      for( const SandiaDecay::RadParticle &particle : transition->products )
      {
        has_xray |= (particle.type == SandiaDecay::ProductType::XrayParticle);
        
        if( particle.type == SandiaDecay::GammaParticle )
          has_coincidences |= !particle.coincidences.empty();
      }//for( const SandiaDecay::RadParticle &particle : transition->products )
    }//for( const SandiaDecay::Transition *transition : nap.nuclide->decaysToChildren )
  }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
  
  BOOST_REQUIRE_MESSAGE( has_coincidences, "U238 decay in SandiaDecay didnt have coincidences" );
  BOOST_REQUIRE_MESSAGE( has_xray, "U238 decay in SandiaDecay didnt have x-rays" );
}//void set_data_dir()

}//namespace


struct NumPartSum
{
  size_t ngamma = 0, nxray = 0, nalpha = 0, nbeta = 0, ncoinc = 0, nrandsum = 0, nde = 0, nse = 0, nannih = 0, nnormal = 0;
};

NumPartSum get_part_summary( const std::vector<ReferenceLineInfo::RefLine> &lines )
{
  NumPartSum answer;
  
  for( const ReferenceLineInfo::RefLine &line : lines )
  {
    //line.m_drf_factor
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:  answer.nalpha += 1; break;
      case ReferenceLineInfo::RefLine::Particle::Beta:   answer.nbeta += 1;  break;
      case ReferenceLineInfo::RefLine::Particle::Gamma:  answer.ngamma += 1; break;
      case ReferenceLineInfo::RefLine::Particle::Xray:   answer.nxray += 1;  break;
    }//switch( line.m_particle_type )
    
    switch( line.m_source_type )
    {
      case ReferenceLineInfo::RefLine::RefGammaType::Normal:             answer.nnormal += 1;  break;
      case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:       answer.nannih += 1;   break;
      case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:       answer.nse += 1;      break;
      case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:       answer.nde += 1;      break;
      case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak: answer.ncoinc += 1;   break;
      case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:       answer.nrandsum += 1; break;
    }//switch( line.m_source_type )
  }//for( const ReferenceLineInfo::RefLine &line : ref_lines->m_ref_lines )
  
  return answer;
}//NumPartSum get_part_summary( const std::vector<ReferenceLineInfo::RefLine> &lines )




// Checks the fields that shouldnt be changed by #ReferenceLineInfo::generateRefLineInfo, are
// not changed
void check_most_input_same( const RefLineInput &lhs, const RefLineInput &rhs )
{
  BOOST_CHECK( lhs.m_color == rhs.m_color );
  BOOST_CHECK_EQUAL( lhs.m_promptLinesOnly, rhs.m_promptLinesOnly );
  BOOST_CHECK_EQUAL( lhs.m_lower_br_cutt_off, rhs.m_lower_br_cutt_off );
  BOOST_CHECK_EQUAL( lhs.m_showGammas, rhs.m_showGammas );
  BOOST_CHECK_EQUAL( lhs.m_showXrays, rhs.m_showXrays );
  BOOST_CHECK_EQUAL( lhs.m_showAlphas, rhs.m_showAlphas );
  BOOST_CHECK_EQUAL( lhs.m_showBetas, rhs.m_showBetas );
  BOOST_CHECK_EQUAL( lhs.m_showCascades, rhs.m_showCascades );
}//check_most_input_same(...)


// Checks that the specified energy is the largest amplitude line, of the specified
//  particle and line type
void check_largest_line( const float energy,
                        ReferenceLineInfo::RefLine::Particle particle_type,
                        ReferenceLineInfo::RefLine::RefGammaType line_type,
                        const vector<ReferenceLineInfo::RefLine> &lines )
{
  const ReferenceLineInfo::RefLine *current_line = nullptr;
  for( const ReferenceLineInfo::RefLine &line : lines )
  {
    if( (line.m_source_type == line_type)
       && (line.m_particle_type == particle_type)
       && (!current_line || (line.m_normalized_intensity > current_line->m_normalized_intensity)) )
    {
      current_line = &line;
    }
  }
  
  BOOST_CHECK_MESSAGE( current_line, "There was no reference lines with specified Particle and RefGammaType values " );
  if( !current_line )
    return;
  
  BOOST_CHECK_MESSAGE( fabs(energy - current_line->m_energy) < 0.1,
                      "The largest intensity line is not " << energy << " keV, but " << current_line->m_energy << " keV" );
}//check_largest_line(...)


BOOST_AUTO_TEST_CASE( basicRefLineTest )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "137-Cs";
  input.m_age = "2.5 HL";
  input.m_color = Wt::WColor("red");
  input.m_lower_br_cutt_off = 0.0;
  input.m_promptLinesOnly = false;
  
  input.m_showGammas = true;
  input.m_showXrays = true;
  input.m_showAlphas = false;
  input.m_showBetas = false;
  input.m_showCascades = false;
  
  input.m_detector_name = "";
  input.m_det_intrinsic_eff = nullptr;
  
  input.m_shielding_name = "";
  input.m_shielding_thickness = "";
  input.m_shielding_an = "";
  input.m_shielding_ad = "";
  input.m_shielding_att = nullptr;
    
   
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  input.m_input_txt = "";
  input.m_age = "";
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Blank );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::None );
  check_most_input_same( input, ref_lines->m_input );
  
  input.m_input_txt = "Tl201";
  input.m_age = "An Invalid Age";
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::InvalidAge );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::None );
  check_most_input_same( input, ref_lines->m_input );
  
  input.m_input_txt = "Invalid Source";
  input.m_age = "20 y";
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::InvalidSource );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::None );
  check_most_input_same( input, ref_lines->m_input );
  
  
  input.m_input_txt = "137-Cs";
  input.m_age = "2 HL";
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  
  // Check that the name is actually normalized
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_input_txt, "Cs137" );
  check_most_input_same( input, ref_lines->m_input );
  
  BOOST_CHECK( ref_lines->m_nuclide );
  BOOST_CHECK( !ref_lines->m_element );
  BOOST_CHECK( ref_lines->m_reactions.empty() );
  if( ref_lines->m_nuclide )
  {
    BOOST_CHECK( ref_lines->m_nuclide->symbol == "Cs137" );
  }
  
  NumPartSum summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_GT( summary.ngamma, 1 );
  BOOST_CHECK_GT( summary.nxray, 1 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  check_largest_line( 661.66, ReferenceLineInfo::RefLine::Particle::Gamma,
                          ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  check_largest_line( 32.19, ReferenceLineInfo::RefLine::Particle::Xray,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  for( const ReferenceLineInfo::RefLine &line : ref_lines->m_ref_lines )
  {
    BOOST_CHECK( line.m_parent_nuclide && line.m_parent_nuclide->symbol == ref_lines->m_input.m_input_txt );
    BOOST_CHECK( line.m_energy > 0 );
    BOOST_CHECK( !line.m_element );
    BOOST_CHECK( !line.m_reaction );
    BOOST_CHECK( line.m_transition );
    BOOST_CHECK_EQUAL( line.m_drf_factor, 1.0 );
    BOOST_CHECK_EQUAL( line.m_shield_atten, 1.0 );
    BOOST_CHECK( (line.m_particle_type == ReferenceLineInfo::RefLine::Particle::Xray)
                || (line.m_particle_type == ReferenceLineInfo::RefLine::Particle::Gamma) );
    BOOST_CHECK( line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal );
    BOOST_CHECK_GE( line.m_normalized_intensity, 0.0 );
    BOOST_CHECK_LE( line.m_normalized_intensity, 1.00001 );
  }//for( const ReferenceLineInfo::RefLine &line : ref_lines->m_ref_lines )
  
  
  // Check if we dont want gammas, we dont get them
  input.m_input_txt = "U238";
  input.m_age = "20 y";
  input.m_showGammas = false;
  input.m_showXrays = true;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  
  summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_EQUAL( summary.ngamma, 0 );
  BOOST_CHECK_GT( summary.nxray, 10 );  //Actual value is 194
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  
  // Check if we dont want anything, we dont get anything.
  input.m_input_txt = "U238";
  input.m_age = "20 y";
  input.m_showGammas = false;
  input.m_showXrays = false;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  
  summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_EQUAL( summary.ngamma, 0 );
  BOOST_CHECK_EQUAL( summary.nxray, 0 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
}//BOOST_AUTO_TEST_CASE( basicRefLineTest )


BOOST_AUTO_TEST_CASE( testReaction )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "Fe(n,g)";
  input.m_showGammas = true;
  input.m_showXrays = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Reaction );
  BOOST_CHECK( !ref_lines->m_reactions.empty() );
  BOOST_CHECK( !ref_lines->m_element );
  BOOST_CHECK( !ref_lines->m_nuclide );
  BOOST_CHECK( !ref_lines->m_has_coincidences );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_input_txt, "Fe(n,g)" );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_age, "" );
  check_most_input_same( input, ref_lines->m_input );
  
  
  NumPartSum summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_GT( summary.ngamma, 1 );
  BOOST_CHECK_EQUAL( summary.nxray, 0 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  check_largest_line( 7631.13, ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
}//BOOST_AUTO_TEST_CASE( testReaction )


BOOST_AUTO_TEST_CASE( testXRay )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "Pb";
  input.m_showGammas = true;
  input.m_showXrays = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::FluorescenceXray );
  BOOST_CHECK( ref_lines->m_reactions.empty() );
  BOOST_CHECK( ref_lines->m_element );
  BOOST_CHECK( !ref_lines->m_nuclide );
  BOOST_CHECK( !ref_lines->m_has_coincidences );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_input_txt, "Pb" );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_age, "" );
  check_most_input_same( input, ref_lines->m_input );
  
  
  NumPartSum summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_EQUAL( summary.ngamma, 0 );
  BOOST_CHECK_GT( summary.nxray, 1 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  check_largest_line( 75.0, ReferenceLineInfo::RefLine::Particle::Xray,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  for( const ReferenceLineInfo::RefLine &line : ref_lines->m_ref_lines )
  {
    BOOST_CHECK( line.m_element && line.m_element->symbol == ref_lines->m_input.m_input_txt );
    BOOST_CHECK( (line.m_energy > 0) && (line.m_energy < 200) );
    BOOST_CHECK( !line.m_parent_nuclide );
    BOOST_CHECK( !line.m_reaction );
    BOOST_CHECK( !line.m_transition );
    BOOST_CHECK_EQUAL( line.m_drf_factor, 1.0 );
    BOOST_CHECK_EQUAL( line.m_shield_atten, 1.0 );
    BOOST_CHECK( line.m_particle_type == ReferenceLineInfo::RefLine::Particle::Xray );
    BOOST_CHECK( line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal );
    BOOST_CHECK_GE( line.m_normalized_intensity, 0.0 );
    BOOST_CHECK_LE( line.m_normalized_intensity, 1.00001 );
  }
}//BOOST_AUTO_TEST_CASE( testXRay )


BOOST_AUTO_TEST_CASE( testCustomEnergy )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "511 keV";
  input.m_showGammas = true;
  input.m_showXrays = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::CustomEnergy );
  BOOST_CHECK( ref_lines->m_reactions.empty() );
  BOOST_CHECK( !ref_lines->m_element );
  BOOST_CHECK( !ref_lines->m_nuclide );
  BOOST_CHECK( !ref_lines->m_has_coincidences );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_input_txt, "511 keV" );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_age, "" );
  check_most_input_same( input, ref_lines->m_input );
  
  
  NumPartSum summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_EQUAL( summary.ngamma, 1 );
  BOOST_CHECK_EQUAL( summary.nxray, 0 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  check_largest_line( 511, ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  BOOST_REQUIRE_MESSAGE( ref_lines->m_ref_lines.size() == 1,
                        "We expect exactly 1 ref line from custom energy input" );
  
  BOOST_CHECK( ref_lines->m_ref_lines[0].m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal );
  
  for( const ReferenceLineInfo::RefLine &line : ref_lines->m_ref_lines )
  {
    BOOST_CHECK( !line.m_element );
    BOOST_CHECK( fabs(line.m_energy - 511) < 0.0001 );
    BOOST_CHECK( !line.m_parent_nuclide );
    BOOST_CHECK( !line.m_reaction );
    BOOST_CHECK( !line.m_transition );
    BOOST_CHECK_EQUAL( line.m_drf_factor, 1.0 );
    BOOST_CHECK_EQUAL( line.m_shield_atten, 1.0 );
    BOOST_CHECK( line.m_particle_type == ReferenceLineInfo::RefLine::Particle::Gamma );
    BOOST_CHECK( line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal );
    BOOST_CHECK_GE( line.m_normalized_intensity, 0.0 );
    BOOST_CHECK_LE( line.m_normalized_intensity, 1.00001 );
  }
}//BOOST_AUTO_TEST_CASE( testCustomEnergy )


BOOST_AUTO_TEST_CASE( testBackground )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "Background";
  input.m_showGammas = true;
  input.m_showXrays = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Background );
  BOOST_CHECK( ref_lines->m_reactions.empty() );
  BOOST_CHECK( !ref_lines->m_element );
  BOOST_CHECK( !ref_lines->m_nuclide );
  BOOST_CHECK( !ref_lines->m_has_coincidences );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_input_txt, "background" );
  BOOST_CHECK_EQUAL( ref_lines->m_input.m_age, "" );
  check_most_input_same( input, ref_lines->m_input );
  
  // We will just make sure we get some lines - wont be too detailed.
  NumPartSum summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_GT( summary.ngamma, 50 ); //should be 77
  BOOST_CHECK_GT( summary.nxray, 5 );   //should be 12
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_GE( summary.nse, 1 );  //should be 1
  BOOST_CHECK_GE( summary.nde, 1 );  //should be 1
  BOOST_CHECK_EQUAL( summary.ncoinc, 0 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
}//BOOST_AUTO_TEST_CASE( testBackground )


BOOST_AUTO_TEST_CASE( testDrfEffect )
{
  SetDataDirs::set_data_dir();
  
  // Init a generic DRF, to make sure it gets applied
  const string drf_dir = SpecUtils::append_path( SetDataDirs::sm_static_data_dir, "GenericGadrasDetectors/HPGe 10%/" );
  BOOST_CHECK_MESSAGE( SpecUtils::is_directory(drf_dir), "Default DRF somehow missing? '" << drf_dir << "'" );
  std::shared_ptr<DetectorPeakResponse> drf;
  BOOST_REQUIRE_NO_THROW( drf = DrfSelect::initAGadrasDetectorFromDirectory(drf_dir) );
  BOOST_REQUIRE_MESSAGE( drf, "Default DRF somehow missing? '" << drf_dir << "'" );
  
  RefLineInput input;
  input.m_input_txt = "Co60";
  input.m_showGammas = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  check_most_input_same( input, ref_lines->m_input );
  
  // With no DRF, the 1332 keV line is the largest.
  check_largest_line( 1332.49,
                      ReferenceLineInfo::RefLine::Particle::Gamma,
                      ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );

  // With a DRF, the 1173.23 line becomes largest
  input.m_detector_name = drf->name();
  input.m_det_intrinsic_eff = drf->intrinsicEfficiencyFcn();
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  check_most_input_same( input, ref_lines->m_input );
  
  check_largest_line( 1173.23,
                     ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  // TODO: go through and check the DRF values are getting applied correctly
}//BOOST_AUTO_TEST_CASE( testDrfEffect )


BOOST_AUTO_TEST_CASE( testShieldingEffect )
{
  SetDataDirs::set_data_dir();
  
  // We will do a basic sanity check that shielding is actually getting aplied
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  
  MaterialDB materialDB;
  BOOST_REQUIRE_NO_THROW( materialDB.parseGadrasMaterialFile( materialfile, db, false ) );
  
  
  RefLineInput input;
  input.m_input_txt = "U238";
  input.m_age = "20 y";
  
  input.m_showGammas = true;
  input.m_showXrays = true;

  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  BOOST_CHECK( ref_lines->m_nuclide && (ref_lines->m_nuclide->symbol == "U238") );
  check_most_input_same( input, ref_lines->m_input );
  
  check_largest_line( 63.29,
                    ReferenceLineInfo::RefLine::Particle::Gamma,
                    ReferenceLineInfo::RefLine::RefGammaType::Normal,
                    ref_lines->m_ref_lines );
  
  input.m_shielding_name = "Pb (lead)";
  input.m_shielding_thickness = "1 cm";
  input.m_shielding_an = "";
  input.m_shielding_ad = "";
  input.m_shielding_att = nullptr;
  
  BOOST_REQUIRE_NO_THROW( input.setShieldingAttFcn( &materialDB ) );
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  
  check_largest_line( 1000.99,
                     ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  // TODO: go through and check the attenuation gets applied correctly
}//BOOST_AUTO_TEST_CASE( testShieldingEffect )


BOOST_AUTO_TEST_CASE( testAge )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "U238";
  input.m_age = "0 y";
  
  input.m_showGammas = true;
  input.m_showXrays = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  BOOST_CHECK( ref_lines->m_nuclide && (ref_lines->m_nuclide->symbol == "U238") );
  check_most_input_same( input, ref_lines->m_input );
  
  check_largest_line( 49.55,
                     ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  
  input.m_age = "20 y";
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  check_largest_line( 63.29,
                     ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
  
  input.m_age = "1 HL";
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  check_largest_line( 609.31,
                     ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::Normal,
                     ref_lines->m_ref_lines );
}//BOOST_AUTO_TEST_CASE( testAge )



BOOST_AUTO_TEST_CASE( testCascade )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "U238";
  input.m_age = "20 y";
  
  // Check we can get cascades with gammas/x-rays
  input.m_showGammas = true;
  input.m_showXrays = true;
  input.m_showCascades = true;
  
  std::shared_ptr<ReferenceLineInfo> ref_lines;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  BOOST_CHECK( ref_lines->m_nuclide && (ref_lines->m_nuclide->symbol == "U238") );
  check_most_input_same( input, ref_lines->m_input );
  
  
  NumPartSum summary = get_part_summary( ref_lines->m_ref_lines );
  
  BOOST_CHECK_GT( summary.ngamma, 100 );
  BOOST_CHECK_GT( summary.nxray, 50 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_GT( summary.ncoinc, 100 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  
  // Check we dont get gammas/xrays if we dont want them
  input.m_showGammas = false;
  input.m_showXrays = false;
  input.m_showCascades = true;
  BOOST_REQUIRE_NO_THROW( ref_lines = ReferenceLineInfo::generateRefLineInfo(input) );
  BOOST_REQUIRE_MESSAGE( ref_lines, "Ref lines for '" << input.m_input_txt << "' should not be null");
  BOOST_CHECK( ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid );
  BOOST_CHECK( ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide );
  
  summary = get_part_summary( ref_lines->m_ref_lines );
  BOOST_CHECK_EQUAL( summary.ngamma, 0 );
  BOOST_CHECK_EQUAL( summary.nxray, 0 );
  BOOST_CHECK_EQUAL( summary.nalpha, 0 );
  BOOST_CHECK_EQUAL( summary.nbeta, 0 );
  BOOST_CHECK_EQUAL( summary.nse, 0 );
  BOOST_CHECK_EQUAL( summary.nde, 0 );
  BOOST_CHECK_GT( summary.ncoinc, 100 );
  BOOST_CHECK_EQUAL( summary.nrandsum, 0 );
  
  check_largest_line( 1171.7,
                     ReferenceLineInfo::RefLine::Particle::Gamma,
                     ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak,
                     ref_lines->m_ref_lines );
  
  // TODO: check actual values, as well as check actual effects of DRF and shielding
}//BOOST_AUTO_TEST_CASE( testAge )



BOOST_AUTO_TEST_CASE( testSerialization )
{
  SetDataDirs::set_data_dir();
  
  RefLineInput input;
  input.m_input_txt = "Cs137";
  input.m_age = "4 y";
  input.m_color = Wt::WColor("red");
  input.m_lower_br_cutt_off = 0.1;
  input.m_promptLinesOnly = false;
  
  input.m_showGammas = true;
  input.m_showXrays = true;
  input.m_showAlphas = false;
  input.m_showBetas = false;
  input.m_showCascades = true;
  
  input.m_detector_name = "";
  input.m_det_intrinsic_eff = nullptr;
  
  input.m_shielding_name = "Pb (iron)";
  input.m_shielding_thickness = "1 cm";
  input.m_shielding_an = "";
  input.m_shielding_ad = "";
  input.m_shielding_att = nullptr;
  
  auto check_serialization = []( const RefLineInput &lhs ){
    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( lhs.serialize( &doc ) );
    
    RefLineInput rhs;
    BOOST_REQUIRE_NO_THROW( rhs.deSerialize( doc.first_node() ) );
    
    check_most_input_same( lhs, rhs );
    BOOST_CHECK_EQUAL( lhs.m_input_txt, rhs.m_input_txt );
    BOOST_CHECK_EQUAL( lhs.m_age, rhs.m_age );
  };//check_serialization lambda
  
  check_serialization( input );
  
  
  input.m_lower_br_cutt_off = 0.0;
  input.m_shielding_name = "";
  input.m_shielding_thickness = "";
  input.m_shielding_an = "26";
  input.m_shielding_ad = "10";
  
  check_serialization( input );
}//BOOST_AUTO_TEST_CASE( testSerialization )


