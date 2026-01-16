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
#include <cmath>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PeakModel_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;
using namespace boost::unit_test;


std::string sm_static_data_dir = "";

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

  const std::string sandia_decay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( !sandia_decay_file.empty() && SpecUtils::is_file(sandia_decay_file),
                        "Error finding sandia.decay.xml" );
  
  // Set the static data directory so we have cross-sections, reactions, and all that
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

  sm_static_data_dir = datadir;

  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE_MESSAGE( u238, "Full SandiaDecayDataBase empty?" );
}//void set_data_dir()


BOOST_AUTO_TEST_CASE( testSetNuclideXrayRctn )
{
  set_data_dir();

  auto check_extract_energy = []( std::string testval, const double expected_energy, const std::string expected_str ) {
    const double energy = PeakDef::extract_energy_from_peak_source_string(testval);
    BOOST_CHECK( energy == expected_energy );
    BOOST_CHECK( testval == expected_str );
  };

  check_extract_energy( "fe xray 98.2 kev", 98.2, "fe xray" );
  check_extract_energy( "5.34e+2 kev", 534, "" );
  check_extract_energy( "hf178m 5.34e-3 Mev", 5.34, "hf178m" );
  check_extract_energy( "8.0e+02 kev hf178m", 800, "hf178m" );
  check_extract_energy( "8.0E+02 kev hf178m", 800, "hf178m" );
  check_extract_energy( "hf178m2 574. KEV", 574, "hf178m2" );
  check_extract_energy( "hf178m2 574.", 574, "hf178m2" );
  check_extract_energy( "u232 xray 98.", 98, "u232 xray" );
  check_extract_energy( "u232 xray 98", 98, "u232 xray" );
  check_extract_energy( "u232 98", 98, "u232" );
  check_extract_energy( "98 u232", -1, "98 u232" );
  check_extract_energy( "u-232", -1.0, "u-232" );
  check_extract_energy( "321 u-232", -1.0, "321 u-232" );
  check_extract_energy( "321 keV u-232", 321, "u-232" );
  check_extract_energy( "3.3mev be(a,n)", 3300, "be(a,n)" );
  check_extract_energy( "co60 1173.23", 1173.23, "co60" );
  check_extract_energy( "co60 1173.23 kev", 1173.23, "co60" );
  check_extract_energy( "1173.23 kev co60", 1173.23, "co60" );
  check_extract_energy( "CO60 1173.23", 1173.23, "CO60" );
  check_extract_energy( "CO60 1173", 1173, "CO60" );
  check_extract_energy( "1173 CO60", -1, "1173 CO60" );
  check_extract_energy( "1173.0 CO60", -1, "1173.0 CO60" );
  check_extract_energy( "Pb 98", -1, "Pb 98" );
  check_extract_energy( "Pb 98.2", 98.2, "Pb" );


  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  PeakModel::SetGammaSource result;

  PeakDef peak;
  const SandiaDecay::Nuclide *nuc = nullptr;


  peak = PeakDef( 1001, 1, 1.8E6 );
  nuc = db->nuclide( "U238" );
  result = PeakModel::setNuclide( peak, PeakDef::SourceGammaType::NormalGamma, nuc, 1001, 4.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  BOOST_CHECK( peak.parentNuclide() == nuc );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 1001) < 1.0 );

  nuc = db->nuclide( "Th232" );
  peak = PeakDef( 2614-511, 5, 1.8E6 );
  BOOST_CHECK( !peak.useForShieldingSourceFit() );
  BOOST_CHECK( !peak.useForManualRelEff() );
  result = PeakModel::setNuclide( peak, PeakDef::SourceGammaType::SingleEscapeGamma, nuc, 2614, 4.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( !peak.useForShieldingSourceFit() );
  BOOST_CHECK( !peak.useForManualRelEff() );
  BOOST_CHECK( peak.parentNuclide() == nuc );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2614-511)) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::SingleEscapeGamma );


  nuc = db->nuclide( "Th232" );
  peak = PeakDef( 2614-511-511, 5, 1.8E6 );
  result = PeakModel::setNuclide( peak, PeakDef::SourceGammaType::DoubleEscapeGamma, nuc, 2614, -1 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() == nuc );
  BOOST_CHECK( !peak.useForShieldingSourceFit() );
  BOOST_CHECK( !peak.useForManualRelEff() );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2614 - 511 - 511)) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );

  peak = PeakDef( 2614-511, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Th232 S.E.", -1.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2614-511)) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::SingleEscapeGamma );


  peak = PeakDef( 2614 - 511 - 511, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Th232 D.E.", -1.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2614 - 511 - 511)) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );

  peak = PeakDef( 2614 - 511 - 100, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Th232 2614 keV D.E.", -1.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2614 - 511 - 511)) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );

  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248 keV", -1.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.reaction() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 2223.248) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 100.0, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248", -1.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.reaction() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 2223.248) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g)", 4.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.reaction() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 2223.248) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248 keV S.E.", 4.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.reaction() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2223.248-511)) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::SingleEscapeGamma );

  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248 keV D.E.", 4.0 );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.reaction() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - (2223.248-2*511)) < 2.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );

  peak = PeakDef( 100, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "U xray 98.4340 kev", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.xrayElement() != nullptr );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 98.4340) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 574, 5, 1.8E6 );
  nuc = db->nuclide( "hf178m2" );
  BOOST_REQUIRE( nuc );
  result = PeakModel::setNuclideXrayReaction( peak, "hf178m2 574.219971 kev", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  BOOST_CHECK( peak.parentNuclide() == nuc );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 574.219971) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 574, 5, 1.8E6 );
  nuc = db->nuclide( "hf178m2" );
  BOOST_REQUIRE( nuc );
  result = PeakModel::setNuclideXrayReaction( peak, "hf178m2", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  BOOST_CHECK( peak.parentNuclide() == nuc );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 574.219971) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 100, 5, 1.8E6 );
  nuc = db->nuclide( "hf178m2" );
  BOOST_REQUIRE( nuc );
  result = PeakModel::setNuclideXrayReaction( peak, "hf178m2 574.219971", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  BOOST_CHECK( peak.parentNuclide() == nuc );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 574.219971) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 84.9, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Pb xray 84.9 kev", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.xrayElement() );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 84.9) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 84.9, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Pb 84.9 kev", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.xrayElement() );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 84.9) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );

  peak = PeakDef( 84.9, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Pb212 84.9 kev", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() );
  BOOST_CHECK( peak.parentNuclide()->symbol == "Pb212" );
  BOOST_CHECK( fabs(peak.gammaParticleEnergy() - 84.865) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::XrayGamma );

  peak = PeakDef( 511, 5, 1.8E6 );
  BOOST_CHECK( !peak.hasSourceGammaAssigned() );
  result = PeakModel::setNuclideXrayReaction( peak, "Na22 511 kev", -1. );
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() );
  BOOST_CHECK( peak.parentNuclide()->symbol == "Na22" );
  BOOST_CHECK( fabs( peak.gammaParticleEnergy() - 511 ) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::AnnihilationGamma );
  BOOST_CHECK( peak.hasSourceGammaAssigned() );

  peak = PeakDef( 511, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Na22 511 kev, I=1.8E+02%", -1. ); //e.g., from "Search for Peaks" dialog
  BOOST_CHECK( result == PeakModel::SetGammaSource::SourceChange );
  BOOST_CHECK( peak.parentNuclide() );
  BOOST_CHECK( peak.parentNuclide()->symbol == "Na22" );
  BOOST_CHECK( fabs( peak.gammaParticleEnergy() - 511 ) < 1.0 );
  BOOST_CHECK( peak.sourceGammaType() == PeakDef::SourceGammaType::AnnihilationGamma );
  BOOST_CHECK( peak.hasSourceGammaAssigned() );
}
