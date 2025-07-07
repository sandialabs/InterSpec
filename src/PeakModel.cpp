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

#include <deque>
#include <string>
#include <vector>
#include <memory>
#include <cctype>
#include <algorithm>

#include <Wt/WResource>
#include <Wt/WModelIndex>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStringStream>
#include <Wt/WAbstractItemModel>


// Block out some warnings occurring in boost files.
#pragma warning(disable:4244)  // warning C4244: 'initializing' : conversion from 'std::streamoff' to 'size_t', possible loss of data

#include <boost/any.hpp>
#include <boost/tokenizer.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/PeakInfoDisplay.h"  //Only for ALLOW_PEAK_COLOR_DELEGATE
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;
using namespace Wt;

#if( PERFORM_DEVELOPER_CHECKS )
namespace
{
void testSetNuclideXrayRctn()
{
  auto check_extract_energy = []( std::string testval, const double expected_energy, const std::string expected_str ) {
    const double energy = PeakDef::extract_energy_from_peak_source_string(testval);
    assert( energy == expected_energy );
    assert( testval == expected_str );
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
  assert( db );
  
  PeakModel::SetGammaSource result;
  
  PeakDef peak;
  const SandiaDecay::Nuclide *nuc = nullptr;
  
  
  peak = PeakDef( 1001, 1, 1.8E6 );
  nuc = db->nuclide( "U238" );
  result = PeakModel::setNuclide( peak, PeakDef::SourceGammaType::NormalGamma, nuc, 1001, 4.0 );
  assert( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  assert( peak.parentNuclide() == nuc );
  assert( fabs(peak.gammaParticleEnergy() - 1001) < 1.0 );
  
  nuc = db->nuclide( "Th232" );
  peak = PeakDef( 2614-511, 5, 1.8E6 );
  assert( !peak.useForShieldingSourceFit() );
  assert( !peak.useForManualRelEff() );
  result = PeakModel::setNuclide( peak, PeakDef::SourceGammaType::SingleEscapeGamma, nuc, 2614, 4.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( !peak.useForShieldingSourceFit() );
  assert( !peak.useForManualRelEff() );
  assert( peak.parentNuclide() == nuc );
  assert( fabs(peak.gammaParticleEnergy() - (2614-511)) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::SingleEscapeGamma );
  
  
  nuc = db->nuclide( "Th232" );
  peak = PeakDef( 2614-511-511, 5, 1.8E6 );
  result = PeakModel::setNuclide( peak, PeakDef::SourceGammaType::DoubleEscapeGamma, nuc, 2614, -1 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() == nuc );
  assert( !peak.useForShieldingSourceFit() );
  assert( !peak.useForManualRelEff() );
  assert( fabs(peak.gammaParticleEnergy() - (2614 - 511 - 511)) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );
  
  peak = PeakDef( 2614-511, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Th232 S.E.", -1.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - (2614-511)) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::SingleEscapeGamma );
  
  
  peak = PeakDef( 2614 - 511 - 511, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Th232 D.E.", -1.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - (2614 - 511 - 511)) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );
  
  peak = PeakDef( 2614 - 511 - 100, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Th232 2614 keV D.E.", -1.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - (2614 - 511 - 511)) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );
  
  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248 keV", -1.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.reaction() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - 2223.248) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 100.0, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248", -1.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.reaction() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - 2223.248) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g)", 4.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.reaction() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - 2223.248) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248 keV S.E.", 4.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.reaction() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - (2223.248-511)) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::SingleEscapeGamma );
  
  peak = PeakDef( 2223.248, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "H(n,g) 2223.248 keV D.E.", 4.0 );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.reaction() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - (2223.248-2*511)) < 2.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::DoubleEscapeGamma );
  
  peak = PeakDef( 100, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "U xray 98.4340 kev", -1. );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.xrayElement() != nullptr );
  assert( fabs(peak.gammaParticleEnergy() - 98.4340) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 574, 5, 1.8E6 );
  nuc = db->nuclide( "hf178m2" );
  assert( nuc );
  result = PeakModel::setNuclideXrayReaction( peak, "hf178m2 574.219971 kev", -1. );
  assert( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  assert( peak.parentNuclide() == nuc );
  assert( fabs(peak.gammaParticleEnergy() - 574.219971) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 574, 5, 1.8E6 );
  nuc = db->nuclide( "hf178m2" );
  assert( nuc );
  result = PeakModel::setNuclideXrayReaction( peak, "hf178m2", -1. );
  assert( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  assert( peak.parentNuclide() == nuc );
  assert( fabs(peak.gammaParticleEnergy() - 574.219971) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 100, 5, 1.8E6 );
  nuc = db->nuclide( "hf178m2" );
  assert( nuc );
  result = PeakModel::setNuclideXrayReaction( peak, "hf178m2 574.219971", -1. );
  assert( result == PeakModel::SetGammaSource::SourceAndUseChanged );
  assert( peak.parentNuclide() == nuc );
  assert( fabs(peak.gammaParticleEnergy() - 574.219971) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 84.9, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Pb xray 84.9 kev", -1. );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.xrayElement() );
  assert( fabs(peak.gammaParticleEnergy() - 84.9) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 84.9, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Pb 84.9 kev", -1. );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.xrayElement() );
  assert( fabs(peak.gammaParticleEnergy() - 84.9) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma );
  
  peak = PeakDef( 84.9, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Pb212 84.9 kev", -1. );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() );
  assert( peak.parentNuclide()->symbol == "Pb212" );
  assert( fabs(peak.gammaParticleEnergy() - 84.865) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::XrayGamma );

  peak = PeakDef( 511, 5, 1.8E6 );
  assert( !peak.hasSourceGammaAssigned() );
  result = PeakModel::setNuclideXrayReaction( peak, "Na22 511 kev", -1. );
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() );
  assert( peak.parentNuclide()->symbol == "Na22" );
  assert( fabs( peak.gammaParticleEnergy() - 511 ) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::AnnihilationGamma );
  assert( peak.hasSourceGammaAssigned() );

  peak = PeakDef( 511, 5, 1.8E6 );
  result = PeakModel::setNuclideXrayReaction( peak, "Na22 511 kev, I=1.8E+02%", -1. ); //e.g., from "Search for Peaks" dialog
  assert( result == PeakModel::SetGammaSource::SourceChange );
  assert( peak.parentNuclide() );
  assert( peak.parentNuclide()->symbol == "Na22" );
  assert( fabs( peak.gammaParticleEnergy() - 511 ) < 1.0 );
  assert( peak.sourceGammaType() == PeakDef::SourceGammaType::AnnihilationGamma );
  assert( peak.hasSourceGammaAssigned() );
}//void testSetNuclideXrayRctn()

}//namespace
#endif


bool PeakModel::recommendUseForFit( const SandiaDecay::Nuclide *nuc,
                                   const float energy )
{
  //These are the lines that should be recomended to users to use in a
  //  shielding source fit
  static const float Ag110m[] = {657.762f, 706.682f, 763.944f, 884.685f, 937.493f, 1384.29f, 1505.04f};
  static const float Al26[]   = {1129.67f, 1808.65f, 2938.0f};
  static const float Am241[]  = {26.35f, 59.54f, 125.3f, 208.01f, 335.37f, 376.65f, 662.4f, 722.01f};
  static const float Am243[]  = {74.66f, 103.032f, 106.123f, 209.753f, 228.183f, 277.599f, 315.88f, 334.31f};
  static const float Au198[]  = {411.802f, 675.884f, 1087.68f};
  static const float Ba133[]  = {80.9971f, 276.4f, 302.853f, 356.017f, 383.848f};
  static const float Ba140[]  = {487.021f, 537.3f, 815.8f, 925.2f, 1596.21f, 2521.7f};
  static const float Be7[]    = {477.595f};
  static const float Bi207[]  = {75.0f, 569.702f, 1063.66f, 1442.2f, 1770.24f};
  static const float Cd109[]  = {88.04f};
  static const float Ce139[]  = {165.85f};
  static const float Ce141[]  = {145.443f};
  static const float Cf249[]  = {109.2f, 252.8f, 333.37f, 388.16f};
  static const float Cm243[]  = {103.734f, 106.125f, 209.753f, 228.183f, 277.599f, 315.88f, 334.31f};
  static const float Cm244[]  = {98.86f, 104.23f, 152.63f, 554.5f, 817.8f};
  static const float Cm245[]  = {103.68f, 133.05f, 175.01f, 231.96f};
  static const float Co56[]   = {846.771f, 1238.28f, 1771.35f, 2598.46f, 3253.42f};
  static const float Co57[]   = {122.061f, 136.474f, 692.03f};
  static const float Co58[]   = {810.759f};
  static const float Co60[]   = {1173.23f, 1332.49f};
  static const float Cr51[]   = {320.082f};
  static const float Cs134[]  = {569.331f, 604.721f, 795.8f, 1038.5f, 1167.86f, 1365.13f};
  static const float Cs136[]  = {340.547f, 818.514f, 1048.07f, 1235.36f, 1537.9f};
  static const float Cs137[]  = {661.657f};
  static const float Cs138[]  = {462.796f, 871.8f, 1009.78f, 1435.86f, 2218.0f, 2639.59f};
  static const float Cu61[]   = {282.956f, 656.008f, 1185.23f};
  static const float Cu67[]   = {184.577f, 300.219f};
  static const float Eu152[]  = {121.781f, 244.697f, 344.278f, 411.116f, 443.965f, 778.904f, 867.378f, 964.079f, 1085.86f, 1112.07f, 1408.01f};
  static const float Eu154[]  = {123.14f, 723.305f, 873.19f, 1004.76f, 1274.39f, 1596.48f};
  static const float Eu155[]  = {86.54f, 105.308f};
  static const float Fe55[]   = {126.0f};
  static const float Fe59[]   = {192.35f, 334.8f, 1099.24f, 1291.59f, 1481.7f};
  static const float Ga67[]   = {93.311f, 184.577f, 300.219f, 393.529f, 887.693f};
  static const float Ge68[]   = {805.75f, 1077.35f, 1260.97f, 1883.09f};
  static const float Ge75[]   = {264.6f, 468.8f};
  static const float Hg194[]  = {293.58f, 328.5f, 948.29f, 1468.89f, 2043.67f};
  static const float Hg203[]  = {279.195f};
  static const float Ho166m[] = {184.41f, 280.459f, 410.944f, 711.683f, 810.276f, 1241.4f};
  static const float I123[]   = {158.97f, 528.96f, 624.57f, 1068.12f};
  static const float I125[]   = {27.472f};
  static const float I126[]   = {388.633f, 666.331f, 753.819f, 1420.17f, 2045.17f};
  static const float I131[]   = {284.305f, 364.489f, 503.004f, 636.989f, 722.911f};
  static const float In111[]  = {171.28f, 245.395f};
  static const float Ir192[]  = {295.957f, 316.506f, 468.069f, 604.411f, 884.537f, 1061.48f};
  static const float K40[]    = {1460.75f};
  static const float La140[]  = {487.0f, 815.8f, 925.2f, 1596.21f, 2521.7f};
  static const float Lu172[]  = {181.525f, 810.064f, 900.724f, 1093.63f, 1584.12f};
  static const float Lu177[]  = {112.94f, 208.36f, 321.31f};
  static const float Lu177m[] = {55.15f, 208.366f, 228.484f, 378.503f, 418.539f};
  static const float Mn54[]   = {834.848f};
  static const float Mn56[]   = {846.754f};
  static const float Mo99[]   = {140.511f, 181.063f, 366.421f, 739.5f, 777.921f};
  static const float Na22[]   = {1274.54f};
  static const float Na24[]   = {996.82f, 1368.63f, 2754.03f, 3866.19f};
  static const float Nb95[]   = {765.803f};
  static const float Np237[]  = {86.477f, 194.95f, 212.29f, 300.129f, 311.904f, 340.476f, 375.45f, 398.492f, 415.76f};
  static const float Np239[]  = {103.032f, 106.123f, 209.753f, 228.183f, 277.599f, 315.88f, 334.31f};
  static const float Pb202[]  = {439.56f};
  static const float Pd103[]  = {294.98f, 357.47f, 497.08f};
  static const float Po210[]  = {803.13f};
  static const float Pu238[]  = {43.498f, 99.853f, 152.72f, 742.81f, 766.39f, 1001.03f};
  static const float Pu239[]  = {129.297f, 203.55f, 255.38f, 345.01f, 375.05f, 413.71f, 451.48f, 639.97f, 645.896f, 769.37f};
  static const float Pu240[]  = {45.244f, 104.23f, 160.308f, 642.35f, 687.59f};
  static const float Pu241[]  = {59.54f, 98.97f, 103.68f, 148.567f, 208.0f};
  static const float Pu242[]  = {44.915f};
  static const float Ra226[]  = {74.815f, 77.107f, 186.211f, 241.997f, 295.24f, 351.932f, 609.312f, 768.356f, 1120.29f, 1238.11f, 1377.67f, 1407.98f, 1764.49f, 2204.21f, 2447.86f};
  static const float Th232[]  = {77.34f, 238.632f, 338.32f, 583.191f, 727.33f, 860.564f, 911.204f, 968.971f, 1588.2f, 2614.53f};
  static const float U233[]   = {146.345f, 164.522f, 291.354f, 317.16f, 1567.09f};
  static const float U235[]   = {106.587f, 143.76f, 163.38f, 185.71f, 205.309f};
  static const float U238[]   = {63.24f, 92.35f, 258.26f, 742.77f, 766.37f, 1000.99f, 1737.73f, 1831.36f};
  static const std::string symbols[]
         = { "Ag110m", "Al26", "Am241", "Am243", "Au198", "Ba133", "Ba140",
             "Be7", "Bi207", "Cd109", "Ce139", "Ce141", "Cf249", "Cm243",
             "Cm244", "Cm245", "Co56", "Co57", "Co58", "Co60", "Cr51", "Cs134",
             "Cs136", "Cs137", "Cs138", "Cu61", "Cu67", "Eu152", "Eu154",
             "Eu155", "Fe55", "Fe59", "Ga67", "Ge68", "Ge75", "Hg194", "Hg203",
             "Ho166m", "I123", "I125", "I126", "I131", "In111", "Ir192", "K40",
             "La140", "Lu172", "Lu177", "Lu177m", "Mn54", "Mn56", "Mo99",
             "Na22", "Na24", "Nb95", "Np237", "Np239", "Pb202", "Pd103",
             "Po210", "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", "Ra226",
             "Th232", "U233", "U235", "U238"
           };
  static const size_t num_sumbols = sizeof(symbols) / sizeof(symbols[0]);
  static const size_t num_energies[]
           = { 7, 3, 8, 8, 3, 5, 6, 1, 5, 1, 1, 1, 4, 7, 5, 4, 5, 3, 1, 2, 1,
               6, 5, 1, 6, 3, 2, 11, 6, 2, 1, 5, 5, 4, 2, 5, 1, 6, 4, 1, 5, 5,
               2, 6, 1, 5, 5, 3, 5, 1, 1, 5, 1, 4, 1, 9, 7, 1, 3, 1, 6, 10, 5,
               5, 1, 15, 10, 5, 5, 8
             };
  static const float *isotope_energies[]
         = { &(Ag110m[0]), &(Al26[0]), &(Am241[0]), &(Am243[0]), &(Au198[0]),
             &(Ba133[0]), &(Ba140[0]), &(Be7[0]), &(Bi207[0]), &(Cd109[0]),
             &(Ce139[0]), &(Ce141[0]), &(Cf249[0]), &(Cm243[0]), &(Cm244[0]),
             &(Cm245[0]), &(Co56[0]), &(Co57[0]), &(Co58[0]), &(Co60[0]),
             &(Cr51[0]), &(Cs134[0]), &(Cs136[0]), &(Cs137[0]), &(Cs138[0]),
             &(Cu61[0]), &(Cu67[0]), &(Eu152[0]), &(Eu154[0]), &(Eu155[0]),
             &(Fe55[0]), &(Fe59[0]), &(Ga67[0]), &(Ge68[0]), &(Ge75[0]),
             &(Hg194[0]), &(Hg203[0]), &(Ho166m[0]), &(I123[0]), &(I125[0]),
             &(I126[0]), &(I131[0]), &(In111[0]), &(Ir192[0]), &(K40[0]),
             &(La140[0]), &(Lu172[0]), &(Lu177[0]), &(Lu177m[0]), &(Mn54[0]),
             &(Mn56[0]), &(Mo99[0]), &(Na22[0]), &(Na24[0]), &(Nb95[0]),
             &(Np237[0]), &(Np239[0]), &(Pb202[0]), &(Pd103[0]), &(Po210[0]),
             &(Pu238[0]), &(Pu239[0]), &(Pu240[0]), &(Pu241[0]), &(Pu242[0]),
             &(Ra226[0]), &(Th232[0]), &(U233[0]), &(U235[0]), &(U238[0])
           };

  if( !nuc )
    return false;
  
  const string *pos = lower_bound( symbols, symbols+num_sumbols, nuc->symbol );
  if( pos && (pos!=(symbols+num_sumbols)) && (*pos)==nuc->symbol )
  {
    const float *energies = isotope_energies[pos-symbols];
    const size_t nenergies = num_energies[pos-symbols];
    const float *epos = lower_bound( energies, energies+nenergies, energy );
    
    //*epos will be greater or equal to energy
    if( epos == (energies+nenergies) )
      return (fabs((*(epos-1))-energy)<0.1);  //assumes at least one element in array
    
    return (((*epos)==energy) || (fabs((*epos)-energy)<0.1)
             || ((epos!=energies) && (fabs((*(epos-1))-energy)<0.1)));
  }//if( we have info on this elemement )
  
  //todo - put in logic here as how to guess if a photopeak should be used or not
  
//  if( nuc && nuc->)
//  if( fabs(511.0-energy) < 1.0 )
//    return false;
  
  return true;
}//recommendUseForFit( const SandiaDecay::Nuclide *nuc )


bool PeakModel::recommendUseForManualRelEff( const SandiaDecay::Nuclide *n, const float energy )
{
  // We currently only have customized recommendations for U, Pu, and Am
  if( !n || ((n->atomicNumber != 92) && (n->atomicNumber != 94) && (n->atomicNumber != 95)) )
    return (energy > 90.0f);
  
  const auto use_gamma = [energy]( const float * const start, const float * const end ) -> bool {
    for( auto iter = start; iter != end; ++iter )
    {
      if( fabs(energy - *iter) < 0.1 ) // 0.1 is arbitrary
        return true;
    }
    return false;
  };// use_gamma lamda
  
  
  if( n->atomicNumber == 92 )
  {
    // Recommended values for uranium taken from chapter 14 of FRMAC Gamma Spectroscopist Knowledge Guide
    const float u232_energies[] = { 238.625f, 583.187f, 727.3f, 860.56f };
    const float u234_energies[] = { 120.905f };
    const float u235_energies[] = { 143.76f, 163.36f, 185.715f, 202.11f, 205.311f, 221.38f, 246.84f, 345.9f };
    const float u238_energies[] = { 258.26f, 569.3173913f, 742.83f, 766.4f, 880.47f, 883.24f, 945.95f, 1001.03f };
    
    switch( n->massNumber )
    {
      case 232: return use_gamma( begin(u232_energies), end(u232_energies) );
      case 234: return use_gamma( begin(u234_energies), end(u234_energies) );
      case 235: return use_gamma( begin(u235_energies), end(u235_energies) );
      case 238: return use_gamma( begin(u238_energies), end(u238_energies) );
      default:
        break;
    }//switch( n->massNumber )
  }//if( Uranium )
  
  
  if( n->atomicNumber == 94 )
  {
    // Recommended value taken from FRAM - see LA-UR-20-21287 Duc T. Vo, and Thomas E. Sampson.
    //  https://www.osti.gov/servlets/purl/1599022
    const float pu238_energies[] = { 152.72 };
    const float pu239_energies[] = { 129.3, 144.2, 161.45, 203.55, 255.38, 345.01, 375.05, 413.71, 451.48, 645.9, 658.86 };
    const float pu240_energies[] = { 160.31 };
    const float pu241_energies[] = { 146.55, 164.61, 208, 267.54, 619.01, 722.01 };
    
    switch( n->massNumber )
    {
      case 238: return use_gamma( begin(pu238_energies), end(pu238_energies) );
      case 239: return use_gamma( begin(pu239_energies), end(pu239_energies) );
      case 240: return use_gamma( begin(pu240_energies), end(pu240_energies) );
      case 241: return use_gamma( begin(pu241_energies), end(pu241_energies) );
      default:
        break;
    }//switch( n->massNumber )
  }//if( Plutonium )
  
  
  if( n->atomicNumber == 95 )
  {
    // Recommended value taken from FRAM - see LA-UR-20-21287 Duc T. Vo, and Thomas E. Sampson.
    //  https://www.osti.gov/servlets/purl/1599022
    const float am241_energies[] = { 125.3, 335.37, 368.65 };
    
    switch( n->massNumber )
    {
      case 241: return use_gamma( begin(am241_energies), end(am241_energies) );
      default:
        break;
    }//switch( n->massNumber )
  }//if( Americium )
  
  // Some other U/Pu/Am isotope - use it if its above the x-ray absorption edge.
  return (energy > 122.0f);
}//bool recommendUseForManualRelEff( const SandiaDecay::Nuclide *n, const float energy )


std::vector<PeakDef> PeakModel::csv_to_candidate_fit_peaks(
                                      std::shared_ptr<const SpecUtils::Measurement> meas,
                                      std::istream &csv )
{
  //Info that will be parsed is based on PeakModel::PeakCsvResource::handleRequest(...)
  //  with a few small accommodations for how other programs may save peak CSVs.
  
  using SpecUtils::trim_copy;
  using SpecUtils::to_lower_ascii_copy;
  
  typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
  boost::escaped_list_separator<char> separator("\\",",\t", "\"");
  
  if( !meas || !meas->gamma_counts() || meas->gamma_counts()->size() < 7 )
    throw runtime_error( "input data invalid" );
  
  const float minenergy = meas->gamma_energy_min();
  const float maxenergy = meas->gamma_energy_max();
  
  
  string line;
  
  //Get first non-empty, non-comment ('#') line.
  while( SpecUtils::safe_get_line(csv, line, 2048) )
  {
    SpecUtils::trim(line);
    if( !line.empty() && line[0]!='#' )
      break;
  }
  
  if( line.empty() || !csv )
    throw runtime_error( "Failed to get first line" );
  
  //Columns guaranteed to be in file, or we'll throw an exception.
  int mean_index = -1, area_index = -1, fwhm_index = -1;
  
  //Columns that may or not be in file, in whcih case will be >= 0.
  int roi_lower_index = -1, roi_upper_index = -1, nuc_index = -1, nuc_energy_index = -1;
  int color_index = -1, label_index = -1, cont_type_index = -1, skew_type_index = -1;
  int cont_coef_index = -1, skew_coef_index = -1, area_uncert_index = -1, peak_type_index = -1;
  
  {//begin to get field_pos
    vector<string> headers;
    Tokeniser t( line, separator );
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      headers.push_back( to_lower_ascii_copy( trim_copy(*it) ) );
    
    const auto centroid_pos = std::find( begin(headers), end(headers), "centroid");
    const auto net_area_pos = std::find( begin(headers), end(headers), "net_area");
    //Note First FWHM is in keV; there is a second FWHM that is %
    const auto fwhm_pos = std::find( begin(headers), end(headers), "fwhm");
    const auto nuc_pos = std::find( begin(headers), end(headers), "nuclide");
    const auto nuc_energy_pos = std::find( begin(headers), end(headers), "photopeak_energy");
    const auto roi_start_pos = std::find( begin(headers), end(headers), "roi_lower_energy");
    const auto roi_end_pos = std::find( begin(headers), end(headers), "roi_upper_energy");
    const auto color_pos = std::find( begin(headers), end(headers), "color");
    const auto label_pos = std::find( begin(headers), end(headers), "user_label");
    const auto cont_type_pos = std::find( begin(headers), end(headers), "continuum_type");
    const auto skew_type_pos = std::find( begin(headers), end(headers), "skew_type");
    const auto cont_coef_pos = std::find( begin(headers), end(headers), "continuum_coefficients");
    const auto skew_coef_pos = std::find( begin(headers), end(headers), "skew_coefficients");
    const auto peak_type_pos = std::find( begin(headers), end(headers), "peak_type");
    
    if( centroid_pos == end(headers) )
      throw runtime_error( "Header did not contain 'Centroid'" );
    mean_index = static_cast<int>( centroid_pos - begin(headers) );
    
    if( net_area_pos == end(headers) )
      throw runtime_error( "Header did not contain 'Net_Area'" );
    area_index = static_cast<int>( net_area_pos - begin(headers) );
    
    // The second "Net_Area" column is uncertainty
    const auto net_area_uncert_pos = std::find( net_area_pos + 1, end(headers), "net_area");
    if( net_area_uncert_pos != end(headers) )
      area_uncert_index = static_cast<int>( net_area_uncert_pos - begin(headers) );
    
    if( fwhm_pos == end(headers) )
      throw runtime_error( "Header did not contain 'FWHM'" );
    fwhm_index = static_cast<int>( fwhm_pos - begin(headers) );
    
    if( nuc_pos != end(headers) )
      nuc_index = static_cast<int>( nuc_pos - begin(headers) );
    
    if( nuc_energy_pos != end(headers) )
      nuc_energy_index = static_cast<int>( nuc_energy_pos - begin(headers) );
    
    if( roi_start_pos != end(headers) )
      roi_lower_index = static_cast<int>( roi_start_pos - begin(headers) );
    
    if( roi_end_pos != end(headers) )
      roi_upper_index = static_cast<int>( roi_end_pos - begin(headers) );
    
    if( color_pos != end(headers) )
      color_index = static_cast<int>( color_pos - begin(headers) );
    
    if( label_pos != end(headers) )
      label_index = static_cast<int>( label_pos - begin(headers) );
    
    if( cont_type_pos != end(headers) )
      cont_type_index = static_cast<int>( cont_type_pos - begin(headers) );
    
    if( skew_type_pos != end(headers) )
      skew_type_index = static_cast<int>( skew_type_pos - begin(headers) );
    
    if( cont_coef_pos != end(headers) )
      cont_coef_index = static_cast<int>( cont_coef_pos - begin(headers) );
    
    if( skew_coef_pos != end(headers) )
      skew_coef_index = static_cast<int>( skew_coef_pos - begin(headers) );
    
    if( peak_type_pos != end(headers) )
      peak_type_index = static_cast<int>( peak_type_pos - begin(headers) );
  }//end to get field_pos
  
  
  vector<PeakDef> answer;
  set<std::shared_ptr<PeakContinuum>> continuums_with_type_set;
  
  while( SpecUtils::safe_get_line(csv, line, 2048) )
  {
    SpecUtils::trim(line);
    
    if( SpecUtils::istarts_with(line, "#END ") || SpecUtils::istarts_with(line, "# END ") )
      break;
    
    if( line.empty() || line[0]=='#' || (!isdigit(line[0]) && line[0]!='+' && line[0]!='-') )
      continue;
    
    vector<string> fields;
    Tokeniser t( line, separator );
    int token_col = 0;
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it, ++token_col )
    {
      string val = trim_copy(*it);
      if( token_col != label_index )
        val = to_lower_ascii_copy( val );
      fields.push_back( val );
    }
    
    const int nfields = static_cast<int>( fields.size() );
    
    if( nfields < mean_index || nfields < area_index || nfields < fwhm_index )
      continue;
    
    try
    {
      const float centroid = std::stof( fields[mean_index] );
      const float fwhm = std::stof( fields[fwhm_index] );
      const float area = std::stof( fields[area_index] );
      
      if( centroid <= minenergy || centroid >= maxenergy || fwhm <= 0.0 || area <= 0.0 )
        continue;
      
      PeakDef peak( centroid, fwhm/2.35482, area );
      
      
      if( (area_uncert_index >= 0) && (area_uncert_index < nfields) )
      {
        const float area_uncert = std::stof( fields[area_uncert_index] );
        if( area_uncert > 0.0 )
          peak.setPeakAreaUncert( area_uncert );
      }//if( (area_uncert_index >= 0) && (area_uncert_index < nfields) )
      
      if( roi_lower_index >= 0 && roi_lower_index < nfields
         && roi_upper_index >= 0 && roi_upper_index < nfields )
      {
        const float roi_lower = std::max( minenergy, std::stof( fields[roi_lower_index] ) );
        const float roi_upper = std::min( maxenergy, std::stof( fields[roi_upper_index] ) );
        if( roi_lower >= roi_upper || centroid < roi_lower || centroid > roi_upper )
          throw runtime_error( "ROI range invalid." );
        
        peak.continuum()->setRange( roi_lower, roi_upper );
        peak.continuum()->calc_linear_continuum_eqn( meas, centroid, roi_lower, roi_upper, 3, 3 );
      }else
      {
        vector<std::shared_ptr<const PeakDef>> peakv( 1, make_shared<const PeakDef>(peak) );
        const auto resolution_type = PeakFitUtils::coarse_resolution_from_peaks( peakv );
        const bool isHPGe = (resolution_type == PeakFitUtils::CoarseResolutionType::High);
        
        double lowerEnengy, upperEnergy;
        findROIEnergyLimits( lowerEnengy, upperEnergy, peak, meas, isHPGe );
        
        peak.continuum()->setRange( lowerEnengy, upperEnergy );
        peak.continuum()->calc_linear_continuum_eqn( meas, centroid, lowerEnengy, upperEnergy, 3, 3 );
      }//if( CSV five ROI extent ) / else( find from data )
      
      if( cont_type_index >= 0 && cont_type_index < nfields )
      {
        const string &strval = fields[cont_type_index];
        
        try
        {
          const PeakContinuum::OffsetType type
                          = PeakContinuum::str_to_offset_type_str( strval.c_str(), strval.size() );
          peak.continuum()->setType( type );
          continuums_with_type_set.insert( peak.continuum() );
          
          if( (cont_coef_index >= 0) && (cont_coef_index < nfields) )
          {
            const string &flt_list_str = fields[cont_coef_index];
            vector<float> values;
            // The delimiter should be a space, but for the moment we'll be loose with this
            //  incase we switch things up.
            SpecUtils::split_to_floats( flt_list_str.c_str(), values, " ,\r\n\t;", false );
            
            if( values.size() == (1 + PeakContinuum::num_parameters(type)) )
            {
              const float ref_energy = values[0];
              vector<double> dvalues( begin(values) + 1, end(values) );
              peak.continuum()->setParameters( ref_energy, dvalues, {} );
            }else
            {
              const string msg = "For peak at " + std::to_string(peak.mean()) + " keV, read in "
                      + std::to_string(values.size()) + " continuum coefficients, but expected 1+"
                      + std::to_string(PeakContinuum::num_parameters(type));
              cerr << msg << endl;
              
#if( PERFORM_DEVELOPER_CHECKS )
              log_developer_error( __func__, msg.c_str() );
              assert( 0 ); //just for development.
#endif
            }//if( we have the reference energy plus correct number of continuum parameters ) / else
          }//if( (cont_coef_index >= 0) && (cont_coef_index < nfields) )
        }catch( std::exception & )
        {
          const string msg = "Failed to convert '" + strval + "' to a PeakContinuum::OffsetType.";
          cerr << msg << endl;
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
        }//try / catch
      }//if( cont_type_index >= 0 && cont_type_index < nfields )
      
      
      if( (skew_type_index >= 0) && (skew_type_index < nfields) )
      {
        const string &strval = fields[skew_type_index];
        
        try
        {
          const PeakDef::SkewType type = PeakDef::skew_from_string( strval );
          peak.setSkewType( type );
          
          if( (skew_coef_index >= 0) && (skew_coef_index < nfields) )
          {
            const string &flt_list_str = fields[skew_coef_index];
            vector<float> values;
            // Again, being loos with allowed delimiters
            SpecUtils::split_to_floats( flt_list_str.c_str(), values, " ,\r\n\t;", false );
            if( values.size() == PeakDef::num_skew_parameters(type) )
            {
              for( size_t i = 0; i < values.size(); ++i )
              {
                const auto coef = PeakDef::CoefficientType(PeakDef::SkewPar0 + i);
                peak.set_coefficient( values[i], coef );
              }//for( set skew coeficient values )
            }else
            {
              const string msg = "For peak at " + std::to_string(peak.mean()) + " keV, read in "
                              + std::to_string(values.size()) + " skew coefficients, but expected "
                              + std::to_string(PeakDef::num_skew_parameters(type));
              
              cerr << msg << endl;
#if( PERFORM_DEVELOPER_CHECKS )
              log_developer_error( __func__, msg.c_str() );
              assert( 0 ); //just for development.
#endif
            }//if( we have correct number of coefficients ) / else
          }//if( (skew_coef_index >= 0) && (skew_coef_index < nfields) )
        }catch( std::exception & )
        {
          const string msg = "Failed to convert '" + strval + "' to a PeakDef::OffsetType.";
          cerr << msg << endl;
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
        }//try / catch
      }//if( (skew_type_index >= 0) && (skew_type_index < nfields) )
      
      if( nuc_index >= 0 && nuc_energy_index >= 0
         && nuc_index < nfields && nuc_energy_index < nfields
         && !fields[nuc_index].empty() && !fields[nuc_energy_index].empty() )
      {
        const string nuctxt = fields[nuc_index] + " " + fields[nuc_energy_index] + " keV";
        const SetGammaSource result = setNuclideXrayReaction( peak, nuctxt, 4.0 );
        if( result == SetGammaSource::NoSourceChange )
          cerr << "csv_to_candidate_fit_peaks: could not assign src txt '"
          << nuctxt << "' as a nuc/xray/rctn" << endl;
      }//if( nuc_index >= 0 || nuc_energy_index >= 0 )
      
      if( (color_index >= 0) && (color_index < nfields) && !fields[color_index].empty() )
      {
        try
        {
          // I dont think WColor will throw, but we'll wrap in try/catch, just in case
          //  Also, we need to remove leading/trailing quotes
          peak.setLineColor( WColor( WString::fromUTF8( fields[color_index] ) ) );
        }catch( std::exception & )
        {
        }
      }//if( we have color index )
      
      if( (label_index >= 0) && (label_index < nfields) && !fields[label_index].empty() )
      {
        // TODO: it looks like all double quote characters never make it here, even if they are in the file correctly
        peak.setUserLabel( fields[label_index] );
      }
      
      if( (peak_type_index >= 0) && (peak_type_index < nfields) )
      {
        const string &strval = fields[peak_type_index];
        
        try
        {
          const PeakDef::DefintionType peak_type = PeakDef::peak_type_from_str( strval.c_str() );
          
          switch( peak_type )
          {
            case PeakDef::DefintionType::GaussianDefined:
              // Nothing to do here.
              break;
              
            case PeakDef::DefintionType::DataDefined:
            {
              peak.m_type = peak_type;
              const auto energycal = meas ? meas->energy_calibration() : nullptr;
              const double lx = peak.lowerX(), ux = peak.upperX();
              
              // Adjust the linear continuum to match the data at the ROI edges.  Right now
              //  we are only doing this for linear continua - should we do for others?
              if( energycal && (peak.continuum()->type() == PeakContinuum::OffsetType::Linear) )
              {
                try
                {
                  const double ref_energy = 0.5*(lx + ux);
                  const size_t start_channel      = meas->find_gamma_channel( lx );
                  const size_t end_channel        = meas->find_gamma_channel( ux );
                  
                  const size_t num_side_bins = 3;
                  double coefficients[2] = { 0.0, 0.0 };
                  
                  PeakContinuum::eqn_from_offsets( start_channel, end_channel, ref_energy,
                                                  meas, num_side_bins, num_side_bins,
                                                  coefficients[1], coefficients[0] );
                  
                  peak.continuum()->setParameters( ref_energy, coefficients, nullptr );
                }catch( std::exception &e )
                {
                  //
                }//try /catch
              }//if( linear continuum )
              
              
              double continuumsum = 0.0, datasum = 0.0;
              if( meas && energycal && energycal->valid() && meas->channel_energies() )
              {
                datasum = meas->gamma_integral( lx, ux );
                continuumsum = peak.continuum()->offset_integral( lx, ux, meas );
              }
              double peaksum = datasum - continuumsum;
              peaksum = ((peaksum >= 0.0) && !IsNan(peaksum) && !IsInf(peaksum)) ? peaksum : 0.0;
              datasum = ((datasum >= 0.0) && !IsNan(datasum) && !IsInf(datasum)) ? datasum : 0.0;
              
              peak.set_coefficient( peaksum, PeakDef::GaussAmplitude );
              peak.set_uncertainty( sqrt(datasum), PeakDef::GaussAmplitude );
              
              break;
            }//case PeakDef::DefintionType::DataDefined:
          }//switch( peak_type )
        }catch( std::exception & )
        {
          const string msg = "Failed to convert '" + strval + "' to a peak type.";
          cerr << msg << endl;
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
        }//try / catch
      }//if( (skew_type_index >= 0) && (skew_type_index < nfields) )
      
      
      //Go through existing peaks and if the new peak should share a ROI, do that here
      if( peak.continuum()->energyRangeDefined() )
      {
        for( auto &p : answer )
        {
          if( !p.continuum()->energyRangeDefined() )
            continue;
          if( fabs(p.continuum()->lowerEnergy() - peak.continuum()->lowerEnergy()) < 0.0001
             && fabs(p.continuum()->upperEnergy() - peak.continuum()->upperEnergy()) < 0.0001 )
          {
            peak.setContinuum( p.continuum() );
            
            if( !continuums_with_type_set.count(p.continuum()) )
            {
              //ToDo: determine if we need to make the continuum quadratic or higher... before InterSpec v1.0.11, this info wasnt in the CSV, and from other programs it is never.  We also need to consider that that we may not yet have seen all the peaks sharing a continuum.
              //peak.setType( OffsetType::Quadratic );
            }//if( !continuums_with_type_set.count(p.continuum()) )
          }//if( ROI bounds are about same as another continuum )
        }//for( loop over peaks previously found to see if we should share continuum )
      }//if( new peak has a defined energy range )
      
      answer.push_back( peak );
    }catch( std::exception &e )
    {
      throw runtime_error( "Invalid value on line '" + line + "', " + string(e.what()) );
    }//try / catch to parse a line into a peak
  }//while( SpecUtils::safe_get_line(csv, line, 2048) )
  
  
  if( answer.empty() )
    throw runtime_error( "No peak rows found in file." );
  
  return answer;
}//csv_to_candidate_fit_peaks(...)


vector<PeakDef> PeakModel::gadras_peak_csv_to_peaks( std::shared_ptr<const SpecUtils::Measurement> meas,
                                                     std::istream &csv )
{
  if( !meas || !meas->gamma_counts() || (meas->gamma_counts()->size() < 7) || !csv )
    throw runtime_error( "input data invalid" );
  
  // TODO: switch to using escaped fields - e.g.:
  //typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
  //boost::escaped_list_separator<char> separator("\\",",\t", "\"");
  //vector<string> fields;
  //Tokeniser t( line, separator );
  //for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
  //  fields.push_back( to_lower_ascii_copy( trim_copy(*it) ) );
  
  enum class PeakCsvFormat{ PeakEasy, Gadras, Unknown };
 
  using SpecUtils::trim;
  using SpecUtils::trim_copy;
  using SpecUtils::to_lower_ascii_copy;
  
  PeakCsvFormat csv_format = PeakCsvFormat::Unknown;
  
  string line;
  while( std::getline( csv, line ) )
  {
    trim( line );
    if( line.empty() || line[0] == '#' )
      continue;
    
    // Line should either be "Energy(keV),sigma,Rate(cps)...", or "Centroid,  Net_Area,   Net_Area"
    vector<string> fields;
    SpecUtils::split_no_delim_compress( fields, line, "," );
    if( fields.empty() )
      continue;
    
    if( fields.size() < 9 )
      throw runtime_error( "Invalid Peak CSV header line: '" + line + "'" );
    
    if( (fields[0] == "Energy(keV)") && (fields[1] == "sigma")
       && (fields[2] == "Rate(cps)") && (fields[3] == "sigma") )
    {
      csv_format = PeakCsvFormat::Gadras;
      break;
    }else if( (fields[0] == "Centroid") && (SpecUtils::trim_copy(fields[1]) == "Net_Area")
             && (SpecUtils::trim_copy(fields[2]) == "Net_Area") && (SpecUtils::trim_copy(fields[3]) == "Peak") )
    {
      if( !std::getline( csv, line ) )
        throw runtime_error( "Failed to get second line of PeakEasy CSV" );
      
      SpecUtils::split_no_delim_compress( fields, line, "," );
      
      if( (fields.size() < 9) 
         || (SpecUtils::trim_copy(fields[0]) != "keV")
         || (SpecUtils::trim_copy(fields[1]) != "Counts")
         || (SpecUtils::trim_copy(fields[2]) != "Uncertainty")
         || (SpecUtils::trim_copy(fields[3]) != "CPS") )
      {
        throw runtime_error( "Second line of PeakEasy CSV file is not correct: '" + line + "'" );
      }
      
      csv_format = PeakCsvFormat::PeakEasy;
      break;
    }else
    {
      throw runtime_error( "Invalid peak CSV line: '" + line + "'" );
    }
  }//while( std::getline( csv, line ) )
  
  
  
  vector<PeakDef> answer;
  
  while( std::getline(csv, line) )
  {
    trim(line);
    if( line.empty() || line[0]=='#' || (!isdigit(line[0]) && line[0]!='+' && line[0]!='-') )
      continue;
    
    vector<string> fields;
    SpecUtils::split_no_delim_compress(fields, line, ",");
    
    const size_t nfields = fields.size();
    if( nfields == 0 )
      continue;
    
    if( nfields < 9 )
      throw runtime_error( "Encountered line in GADRAS CSV file with only "
                          + std::to_string(nfields) + " fields.\n\tLine: \"" + line + "\"" );
    
    const double meas_live_time = (meas && (meas->live_time() > 0.0)) ? meas->live_time() : 1.0f;
    
    try
    {
      double energy, energy_uncert, counts, counts_uncert, fwhm, fwhm_uncert;
      
      switch( csv_format )
      {
        case PeakCsvFormat::PeakEasy:
        {
          for( string &v : fields )
            SpecUtils::trim( v );
          
          energy = std::stod( fields[0] );
          energy_uncert = -1.0;
          counts = std::stod( fields[1] );
          counts_uncert = std::stod( fields[2] );
          fwhm = std::stod( fields[4] );
          fwhm_uncert = -1.0;
          
          const double live_time = ((fields.size() > 10) ? std::stod(fields[10]) : meas_live_time );
          counts /= live_time;
          counts_uncert /= live_time;
          
          //const size_t peak_cps_index = 3, fwhm_percent_index = 5;
          //const size_t roi_total_counts_index = 6, roi_id_index = 7, filename_index = 8;
          //const size_t live_time_index = 9, date_index = 10, time_index = 11;
          
          break;
        }//case PeakCsvFormat::PeakEasy:
          
        case PeakCsvFormat::Gadras:
        {
          energy = std::stod( fields[0] );
          energy_uncert = std::stod( fields[1] );
          counts = std::stod( fields[2] );
          counts_uncert = std::stod( fields[3] );
          fwhm = std::stod( fields[4] );
          fwhm_uncert = std::stod( fields[5] );
          
          //const size_t leakage_per_second_index = 6, centroid_index = 7, filename_index = 8;
          //const size_t record_idx_index = 9, title_index = 10, date_time_index = 11;
          break;
        }
          
        case PeakCsvFormat::Unknown:
          assert( 0 );
          throw runtime_error( "Not a peak CSV file." );
          break;
      }//switch( csv_format )
      
      counts *= meas_live_time;
      counts_uncert *= meas_live_time;
      const double sigma = fwhm / 2.35482f;
      
      PeakDef info( energy, sigma, counts );
      
      if( energy_uncert > 0.0 )
        info.setMeanUncert( energy_uncert );
      
      if( counts_uncert > 0.0 )
        info.setAmplitudeUncert( counts_uncert );
      
      if( fwhm_uncert > 0.0 )
        info.setSigmaUncert( fwhm_uncert  / 2.35482f );
      
      answer.push_back( info );
    }catch( std::exception &e )
    {
      throw runtime_error( "Invalid value on line '" + line + "', " + string(e.what()) );
    }//try / catch to parse a line into a peak
  }//while( SpecUtils::safe_get_line(csv, line, 2048) )
  
  if( answer.empty() )
    throw runtime_error( "No peak rows found in file." );
  
  return answer;
}//gadras_peak_csv_to_peaks(...)



PeakModel::PeakCsvResource::PeakCsvResource( PeakModel *parent )
  : WResource( parent ),
    m_model( parent ),
    m_app( WApplication::instance() )
{
  assert( m_app );
}


PeakModel::PeakCsvResource::~PeakCsvResource()
{
  beingDeleted();
}




void PeakModel::PeakCsvResource::handleRequest( const Wt::Http::Request &/*request*/,
                                                Wt::Http::Response& response )
{
  WApplication::UpdateLock lock( m_app );
  
  if( !lock )
  {
    log("error") << "Failed to WApplication::UpdateLock in PeakCsvResource.";
    response.out() << "Error grabbing application lock to form PeakCsvResource resource; please report to InterSpec@sandia.gov.";
    response.setStatus(500);
    assert( 0 );
    return;
  }//if( !lock )
  
  //If you update the fields or headers of this function, you should also update
  //  PeakModel::csv_to_candidate_fit_peaks(...)
  
  std::shared_ptr<SpecMeas> meas = m_model->m_measurment.lock();
  string filename = "peaks.CSV", specfilename = "Unknown";
  if( meas && !meas->filename().empty() )
  {
    specfilename = meas->filename();
    filename = "peaks_" + specfilename;
    const string::size_type pos = filename.find_last_of( "." );
    if( pos != string::npos )
      filename = filename.substr( 0, pos );
    filename += ".CSV";
  }//if( meas && !meas->filename().empty() )
    
  suggestFileName( filename, WResource::Attachment ); //WResource::NoDisposition

  response.setMimeType( "text/csv" );

  if( !m_model || !m_model->m_peaks )
    return;

  const shared_ptr<const SpecUtils::Measurement> &data = m_model->m_foreground;
  
  string backfilename;
  shared_ptr<const SpecUtils::Measurement> background;
  shared_ptr<const deque<shared_ptr<const PeakDef>>> background_peaks;
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( interspec )
  {
    shared_ptr<const SpecMeas> bmeas = interspec->measurment(SpecUtils::SpectrumType::Background);
    background = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
    if( background && bmeas )
    {
      backfilename = bmeas->filename();
      const set<int> &samples = interspec->displayedSamples(SpecUtils::SpectrumType::Background);
      background_peaks = bmeas->peaks( samples );
    }//if( bmeas )
  }//if( interspec )
  
  
  if( background && background_peaks && !background_peaks->empty() )
  {
    PeakModel::write_for_and_back_peak_csv( response.out(), specfilename,
                                        PeakModel::PeakCsvType::Full, m_model->m_sortedPeaks, data,
                                        backfilename, background_peaks.get(), background );
  }else
  {
    PeakModel::write_peak_csv( response.out(), specfilename, PeakModel::PeakCsvType::Full,
                              m_model->m_sortedPeaks, data );
  }
}//void handleRequest(...)






PeakModel::PeakModel( Wt::WObject *parent )
  : WAbstractItemModel( parent ),
    m_foreground( nullptr ),
    m_sortColumn( kMean ),
    m_sortOrder( Wt::AscendingOrder ),
    m_csvResource( NULL )
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( "PeakModel" );
  
  m_csvResource = new PeakCsvResource( this );
  
#if( PERFORM_DEVELOPER_CHECKS )
  testSetNuclideXrayRctn();
#endif
}//PeakModel constructor

PeakModel::~PeakModel()
{
}//~PeakModel()


void PeakModel::setForeground( shared_ptr<const SpecUtils::Measurement> spec )
{
  m_foreground = std::move(spec);
}//void setForeground( std::shared_ptr<const SpecUtils::Measurement> spec )


void PeakModel::setPeakFromSpecMeas( std::shared_ptr<SpecMeas> meas,
                                     const std::set<int> &samplenums )
{
  std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > > peaks;

  if( !!meas )
    peaks = meas->peaks( samplenums );
  
  m_measurment = meas;
  
  if( peaks == m_peaks )
    return;
  
  if( !!m_peaks && !m_peaks->empty() )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks->size()-1) );
    m_peaks = std::make_shared< deque<std::shared_ptr<const PeakDef> > >();
    m_sortedPeaks.clear();
    endRemoveRows();
  }//if( !m_peaks->empty() )
  
  const size_t npeaksadd = !!peaks ? peaks->size() : size_t(0);

  if( !npeaksadd )
  {
    m_peaks = peaks;
  }else
  {
    boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn, meansort;
    
    const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
    sortfcn = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                          m_sortColumn, m_sortOrder, data );
    meansort = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                           kMean, Wt::AscendingOrder, data );

    //Make sure no null peaks
    peaks->erase( std::remove_if( peaks->begin(), peaks->end(), [](std::shared_ptr<const PeakDef> a){return !a;}), peaks->end() );
    
    
    std::sort( peaks->begin(), peaks->end(), meansort );
    
    beginInsertRows( WModelIndex(), 0, static_cast<int>(npeaksadd-1) );

    m_sortedPeaks = *peaks;
    std::stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );

    m_peaks = peaks;
    endInsertRows();
  }//if( peaks.size() )
  
  // There is an apparent bug in Wt 3.7.1 (at least), that once you remove all the rows and add
  //  new ones back in (at least in the same event loop), then if you scroll down the table, it
  //  will get to indicating its loading more rows, but just never will; this is triggered when
  //  you change spectrum files.
  //  To trigger, I suspect it requires the removal of the rows, and re-adding them back in, since
  //  this issue isnt seen in #PeakModel::addNewPeakInternal.
  //  Since we have to re-render the table anyway, this isnt really that heavy handed.
  //  \sa PeakModel::setPeaks
  //  \sa IsotopeSearchByEnergyModel::updateSearchResults
  //reset();
  layoutAboutToBeChanged().emit();
  layoutChanged().emit();
}//void setPeakFromSpecMeas(...)


void PeakModel::setNoSpecMeasBacking()
{
  if( m_peaks && !m_peaks->empty() )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks->size()-1) );
    m_peaks = std::make_shared< deque<std::shared_ptr<const PeakDef> > >();
    m_sortedPeaks.clear();
    endRemoveRows();
  }else
  {
    m_peaks = std::make_shared< deque<std::shared_ptr<const PeakDef> > >();
  }//if( !m_peaks->empty() )
  
  m_measurment.reset();
  
  layoutAboutToBeChanged().emit();
  layoutChanged().emit();
}//void setNoSpecMeasBacking()


size_t PeakModel::npeaks() const
{
  if( !m_peaks )
    return 0;

  return m_peaks->size();
}//size_t npeaks() const


PeakModel::PeakShrdPtr PeakModel::peakPtr( const size_t peakn ) const
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );
  const PeakShrdPtr &peak = m_peaks->at( peakn );
  assert( peak );
  return peak;
}//PeakModel::PeakShrdPtr peakPtr( const size_t peakn ) const


PeakModel::PeakShrdPtr PeakModel::nearestPeak( double energy ) const
{
  if( !m_peaks || m_peaks->empty() )
    return nullptr;
  
  boost::function<bool( const PeakShrdPtr &, const PeakShrdPtr &)> meansort;
  const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
  meansort = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                         kMean, Wt::AscendingOrder, data );
  PeakDef *new_peak_ptr = new PeakDef();
  PeakShrdPtr peak_ptr( new_peak_ptr );
  new_peak_ptr->setMean( energy );
  deque< PeakShrdPtr >::const_iterator mean_pos
          = lower_bound( m_peaks->begin(), m_peaks->end(), peak_ptr, meansort );
  
  if( mean_pos == m_peaks->end() )
    return *(mean_pos-1);
  
  const double ahead_mean = ((mean_pos+1) != m_peaks->end())
                            ? (*(mean_pos+1))->mean()
                            : DBL_MAX;
  const double behind_mean = (mean_pos != m_peaks->begin())
                            ? (*(mean_pos-1))->mean()
                            : DBL_MAX;
  const double this_mean = (*mean_pos)->mean();
  
  const double ahead_diff = fabs( ahead_mean - energy );
  const double behind_diff = fabs( behind_mean - energy );
  const double this_diff = fabs( this_mean - energy );
  
  if( (ahead_diff < behind_diff) && (ahead_diff < this_diff) )
    return (*(mean_pos+1));
  if( behind_diff < this_diff )
    return (*(mean_pos-1));
  return (*mean_pos);
}//PeakShrdPtr nearestPeak( double energy ) const


Wt::WModelIndex PeakModel::indexOfPeak( const PeakShrdPtr &peak ) const
{
  int row = -1;  
  for( std::deque<PeakShrdPtr>::const_iterator i = m_sortedPeaks.begin();
      i != m_sortedPeaks.end(); ++i )
  {
    ++row;
    if( *i == peak )
      return PeakModel::index(row, 0);
  }
  
  return WModelIndex();
}//Wt::WModelIndex indexOfPeak( const PeakShrdPtr &peak ) const


const PeakDef &PeakModel::peak( const size_t peakn ) const
{
  return  *peakPtr(peakn);
}//const PeakDef &peak( const size_t peakn ) const


std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > PeakModel::peaks() const
{
  return m_peaks;
}//const deque<PeakDef> &peaks() const


std::vector<PeakDef> PeakModel::peakVec() const
{
  vector<PeakDef> answer;
  
  if( !m_peaks )
    return answer;
  
  deque<PeakModel::PeakShrdPtr>::const_iterator pos;
  for( pos = m_peaks->begin(); pos != m_peaks->end(); ++pos )
    answer.push_back( *(*pos) );
  
  return answer;
}//std::vector<PeakDef> peakVec() const


const std::deque<std::shared_ptr<const PeakDef>> &PeakModel::sortedPeaks() const
{
  return m_sortedPeaks;
}


bool PeakModel::isWithinRange( const PeakDef &peak ) const
{
  assert( m_foreground );
  
  if( !m_foreground )
  {
    cerr << "PeakModel::isWithinRange(...)\n\tThere is no xaxis!" << endl;
    return false;
  }

  const auto energycal = m_foreground->energy_calibration();
  assert( energycal );
  
  const double lowerx = m_foreground->gamma_energy_min();
  const double upperx = m_foreground->gamma_energy_max();

  switch( peak.type() )
  {
    case PeakDef::GaussianDefined:
      if( peak.mean() < lowerx || peak.mean() > upperx )
        return false;
    break;
    
    case PeakDef::DataDefined:
      if( peak.lowerX() > upperx || peak.upperX() < lowerx )
        return false;
    break;
  }//switch( peak.type() )

  return true;
} // bool PeakModel::isWithinRange( const PeakDef &peak )


bool PeakModel::isOutOfRange( const PeakDef &peak ) const
{
  return !isWithinRange( peak );
}


WResource *PeakModel::peakCsvResource()
{
  return m_csvResource;
}


void PeakModel::notifySpecMeasOfPeakChange()
{
  std::shared_ptr<SpecMeas> meas = m_measurment.lock();
  if( meas )
    meas->setModified();
}//void notifySpeakMeasOfPeakChange();


WModelIndex PeakModel::addNewPeak( const PeakDef &peak )
{
  return addNewPeakInternal( peak ).second;
}


std::pair<std::shared_ptr<const PeakDef>,Wt::WModelIndex> PeakModel::addNewPeakInternal( const PeakDef &peak )
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );

  if( !isWithinRange( peak ) )
    return { nullptr, WModelIndex{} };

  notifySpecMeasOfPeakChange();
  
  PeakDef *new_peak_ptr = new PeakDef( peak );
  PeakShrdPtr peak_ptr( new_peak_ptr );

  definePeakXRange( *new_peak_ptr );
  //Need to go through and estimate the Chi2Dof here.
  
  boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn, meansort;
  const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
  sortfcn = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                        m_sortColumn, m_sortOrder, data );
  meansort = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                         kMean, Wt::AscendingOrder, data );

  deque< PeakShrdPtr >::iterator mean_pos
         = lower_bound( m_peaks->begin(), m_peaks->end(), peak_ptr, meansort );
  deque< PeakShrdPtr >::iterator sort_pos
       = lower_bound( m_sortedPeaks.begin(), m_sortedPeaks.end(), peak_ptr, sortfcn );
  
  const int indexpos = static_cast<int>( sort_pos - m_sortedPeaks.begin() );

  beginInsertRows( WModelIndex(), indexpos, indexpos );
  m_peaks->insert( mean_pos, peak_ptr );
  m_sortedPeaks.insert( sort_pos, peak_ptr );
  endInsertRows();
  
  return {peak_ptr, index( indexpos, 0 )};
  //return index( indexpos, 0 );
}//void addNewPeakInternal( const PeakDef &peak )



void PeakModel::addPeaks( const vector<PeakDef> &peaks )
{
  for( size_t i = 0; i < peaks.size(); ++i )
    addNewPeak( peaks[i] );
}//void addPeaks( const vector<PeakDef> &peaks )


void PeakModel::addPeaks( const std::vector<std::shared_ptr<const PeakDef>> &peaks )
{
  for( const auto &p : peaks )
  {
    assert( p );
    if( p )
      addNewPeak( *p );
  }//for( const auto &p : peaks )
}//void addPeaks( const std::vector<std::shared_ptr<const PeakDef>> peaks &peaks )


void PeakModel::definePeakXRange( PeakDef &peak )
{
  const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
 
  shared_ptr<PeakContinuum> peakcont = peak.continuum();
  assert( peakcont );
  
  if( peakcont->energyRangeDefined() )
    return;
  
  vector<shared_ptr<const PeakDef>> all_peaks(begin(m_sortedPeaks), end(m_sortedPeaks));
  all_peaks.push_back( make_shared<const PeakDef>(peak) );
  
  auto res_type = PeakFitUtils::coarse_resolution_from_peaks( all_peaks );
  if( res_type == PeakFitUtils::CoarseResolutionType::Unknown )
  {
    InterSpec *interspec = InterSpec::instance();
    if( interspec )
    {
      if( PeakFitUtils::is_likely_high_res(interspec) )
        res_type = PeakFitUtils::CoarseResolutionType::High;
      else
        res_type = PeakFitUtils::CoarseResolutionType::Low;
    }
  }//if( res_type == PeakFitUtils::CoarseResolutionType::Unknown )
  
  const bool isHPGe = (res_type == PeakFitUtils::CoarseResolutionType::High);
  
  {
    // We'll set the peak limits to save cpu (or rather memory-time) later
    
    // TODO: We could do a better job setting the range if external continuum defined.
    //       i.e. if( peakcont->type() == PeakContinuum::External ){}
    
    double lowerEnengy, upperEnergy;
    findROIEnergyLimits( lowerEnengy, upperEnergy, peak, data, isHPGe );
    peakcont->setRange( lowerEnengy, upperEnergy );
  }//if( !peakcont->energyRangeDefined() )
}//void definePeakXRange( PeakDef &peak )


/*
 // This version of function is such that a copy of the peaks are set to the model
void PeakModel::setPeaks( const std::vector<std::shared_ptr<const PeakDef> > &peaks )
{
  vector<PeakDef> peakcopies;
  peakcopies.reserve( peaks.size() );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    if( peaks[i] )
      peakcopies.push_back( *peaks[i] );
  }
  
  setPeaks( peakcopies );
}//void PeakModel::setPeaks()
*/

void PeakModel::setPeaks( vector<shared_ptr<const PeakDef>> peaks )
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );
  
  const bool had_peaks = !m_peaks->empty();
  if( had_peaks )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks->size()-1) );
    m_peaks->clear();
    m_sortedPeaks.clear();
    endRemoveRows();
  }//if( had_peaks )
  
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const shared_ptr<const PeakDef> orig_peak = peaks[i];
    assert( orig_peak );
    if( !orig_peak )
      continue;
    
    const auto c = orig_peak->continuum();
    if( c->energyRangeDefined() )
      continue;
    
    Wt::log("warn") << "PeakModel::setPeaksNoCopy: found peak not energy range defined.";
    
    shared_ptr<PeakDef> new_peak = make_shared<PeakDef>( *orig_peak );
    definePeakXRange( *new_peak );
    peaks[i] = new_peak;
    
    // Go through and replace continuum of any peaks that shared a conitnuum with `orig_peak`
    for( size_t j = i + 1; j < peaks.size(); ++j )
    {
      if( peaks[j] && (peaks[j]->continuum() == c) )
      {
        shared_ptr<PeakDef> new_peak_j = make_shared<PeakDef>( *peaks[j] );
        new_peak_j->setContinuum( new_peak->continuum() );
        peaks[j] = new_peak_j;
      }
    }
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  // Remove any nullptr peaks, or peaks with means larger/smaller than data
  const auto peaks_end = remove_if( begin(peaks), end(peaks),
    [this]( const shared_ptr<const PeakDef> &p ) -> bool {
       return !p || !isWithinRange( *p );
    } );
  
  const size_t npeaksadd = peaks_end - begin(peaks);
  
  if( npeaksadd )
  {
    const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
     
    auto sortfcn = [this, data]( const PeakShrdPtr &lhs, const PeakShrdPtr &rhs ) -> bool {
      return PeakModel::compare( lhs, rhs, m_sortColumn, m_sortOrder, data );
    };
    
    auto meansort = [this, data]( const PeakShrdPtr &lhs, const PeakShrdPtr &rhs ) -> bool {
      return PeakModel::compare( lhs, rhs, kMean, Wt::AscendingOrder, data );
    };

    beginInsertRows( WModelIndex(), 0, int(npeaksadd - 1) );
    
    m_peaks->insert( end(*m_peaks), begin(peaks), peaks_end );
    m_sortedPeaks.insert( end(m_sortedPeaks), begin(peaks), peaks_end );

    std::sort( m_peaks->begin(), m_peaks->end(), meansort );
    std::stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );
    
    endInsertRows();
    
    // Trigger a layout change to keep issue with scrolling in table not working from happening.
    layoutAboutToBeChanged().emit();
    layoutChanged().emit();
  }//if( npeaksadd )
  
  
  if( had_peaks || npeaksadd )
    notifySpecMeasOfPeakChange();
}//void setPeaks( const std::vector<std::shared_ptr<const PeakDef> > &peaks )


void PeakModel::setPeaks( vector<PeakDef> peaks )
{
  vector<shared_ptr<const PeakDef>> peak_ptrs;
  for( const PeakDef &p : peaks )
    peak_ptrs.push_back( make_shared<PeakDef>( p ) );
  setPeaks( peak_ptrs );
  
  /*
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before adding peak" );

  if( !m_peaks->empty() )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks->size()-1) );
    m_peaks->clear();
    m_sortedPeaks.clear();
    endRemoveRows();
  }//if( !m_peaks->empty() )

  // Remove any peaks that aren't within the range.
  const vector<PeakDef>::const_iterator peaks_end
              = remove_if( peaks.begin(), peaks.end(),
                           boost::bind( &PeakModel::isOutOfRange, this,
                                       boost::placeholders::_1 ) );
  
  
  const size_t npeaksadd = peaks_end - peaks.begin();
  
  if( npeaksadd )
  {
    for( size_t i = 0; i < npeaksadd; ++i )
      definePeakXRange( peaks[i] );
//    need to go through and estimate the chi2DOF here
    
    boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn, meansort;
    const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
    sortfcn = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                          m_sortColumn, m_sortOrder, data );
    meansort = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                           kMean, Wt::AscendingOrder, data );

    beginInsertRows( WModelIndex(), 0, int(npeaksadd - 1) );

    for( size_t i = 0; i < npeaksadd; ++i )
    {
      const PeakDef &inpeak = peaks[i];
      std::shared_ptr<PeakDef> peak = std::make_shared<PeakDef>( inpeak );
      m_peaks->push_back( peak );
      m_sortedPeaks.push_back( peak );
    }// for( loop over peaks )

    std::sort( m_peaks->begin(), m_peaks->end(), meansort );
    std::stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );
    
    endInsertRows();
    
    // Trigger a layout change to keep issue with scrolling in table not working from happening.
    //reset();
    layoutAboutToBeChanged().emit();
    layoutChanged().emit();
  }//if( peaks.size() )
  
  notifySpecMeasOfPeakChange();
   */
}//void setPeaks( const vector<PeakDef> &peaks )


void PeakModel::removePeak( const size_t peakn )  //throws if invalid peak numbers
{
  if( peakn >= m_peaks->size() )
    throw std::runtime_error( "PeakModel::removePeak(): invalid index" );

  PeakShrdPtr peak = (*m_peaks)[peakn];

  std::deque< PeakShrdPtr >::iterator sort_pos = find( m_sortedPeaks.begin(), m_sortedPeaks.end(), peak );
  assert( sort_pos != m_sortedPeaks.end() );
  const int index = static_cast<int>( sort_pos - m_sortedPeaks.begin() );

  beginRemoveRows( WModelIndex(), index, index );
  m_peaks->erase( m_peaks->begin() + peakn );
  m_sortedPeaks.erase( sort_pos );
  endRemoveRows();
  
  notifySpecMeasOfPeakChange();
}//void removePeak( const size_t peakn )


void PeakModel::removePeak( Wt::WModelIndex index )
{
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::removePeak(): non-valid index" );

  const int row = index.row();
  if( row < 0 || (static_cast<size_t>(row) >= m_sortedPeaks.size()) )
    throw std::runtime_error( "PeakModel::removePeak(): invalid index row." );
  
  PeakShrdPtr peak = m_sortedPeaks[row];
  
  removePeakInternal( peak );
}//void removePeak( Wt::WModelIndex index )


void PeakModel::removePeak( PeakModel::PeakShrdPtr peak )
{
  removePeakInternal( peak );
}//void removePeak( PeakModel::PeakShrdPtr peak )


void PeakModel::removePeakInternal( PeakModel::PeakShrdPtr peak )
{
  if( !peak )
    throw std::runtime_error( "PeakModel::removePeakInternal: invalid peak passed in." );
  
  if( !m_peaks )
    throw std::runtime_error( "PeakModel::removePeakInternal: data not set." );
  
  const WModelIndex index = indexOfPeak( peak );
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::removePeakInternal: peak passed in doesnt belong to model." );
  
  std::deque< PeakShrdPtr >::iterator energy_pos = find( begin(*m_peaks), end(*m_peaks), peak );
  if( energy_pos == m_peaks->end() )
    throw std::runtime_error( "PeakModel::removePeakInternal: peak passed in doesnt belong to model (or logic error sorting m_peaks)." );
  
  beginRemoveRows( WModelIndex(), index.row(), index.row() );
  m_peaks->erase( energy_pos );
  m_sortedPeaks.erase( m_sortedPeaks.begin() + index.row() );
  endRemoveRows();
  
  notifySpecMeasOfPeakChange();
}//void removePeakInternal( PeakModel::PeakShrdPtr peak )


void PeakModel::removeAllPeaks()
{
  setPeaks( vector<PeakDef>{} );
}//void removeAllPeaks()


void PeakModel::removePeaks( const std::vector<PeakModel::PeakShrdPtr> &peaks )
{
  for( const auto &p : peaks )
    removePeakInternal( p );
}//void removePeaks( const std::vector<PeakModel::PeakShrdPtr> &peak )



vector<shared_ptr<const PeakDef>> PeakModel::peaksSharingRoi( const std::shared_ptr<const PeakDef> &peak )
{
  vector<shared_ptr<const PeakDef>> answer;
  
  if( !peak || !m_peaks )
    return answer;
  
  bool found = false;
  for( auto &p : *m_peaks )
  {
    found |= (p == peak);
    if( p->continuum() == peak->continuum() )
      answer.push_back( p );
  }
  
  if( !found )
    throw std::runtime_error( "PeakModel::peaksSharingRoi: passed in peak not owned by this model." );
  
  return answer;
}//peaksSharingRoi( peak )


std::vector<std::shared_ptr<const PeakDef>> PeakModel::peaksNotSharingRoi( const std::shared_ptr<const PeakDef> &peak )
{
  vector<shared_ptr<const PeakDef>> answer;
  
  if( !m_peaks )
    return answer;
  
  if( !peak )
  {
    for( auto &p : *m_peaks )
      answer.push_back( p );
    return answer;
  }
  
  bool found = false;
  for( auto &p : *m_peaks )
  {
    found |= (p == peak);
    if( p->continuum() != peak->continuum() )
    {
      assert( p != peak );
      answer.push_back( p );
    }//if( p->continuum() != peak->continuum() )
  }//for( auto &p : *m_peaks )
  
  if( !found )
    throw std::runtime_error( "PeakModel::peaksNotSharingRoi: passed in peak not owned by this model." );
  
  return answer;
}//peaksNotSharingRoi( peak )


void PeakModel::updatePeak( const std::shared_ptr<const PeakDef> &originalPeak, const PeakDef &newPeak )
{
  removePeakInternal( originalPeak );
  addNewPeak( newPeak );
}//updatePeak( originalPeak, newPeak )


void PeakModel::updatePeaks( const std::vector<std::shared_ptr<const PeakDef>> &originalPeaks,
                 const std::vector<PeakDef> &newPeaks )
{
  for( const auto &p : originalPeaks )
    removePeakInternal( p );
  for( const auto &p : newPeaks )
    addNewPeakInternal( p );
}//void updatePeaks(...)


void PeakModel::setPeakFitFor( const Wt::WModelIndex index,
                               const PeakDef::CoefficientType coef,
                               const bool fitfor )
{
  if( !m_peaks )
    throw std::runtime_error( "PeakModel::setPeakFitFor(...): no foreground"
                              " spectrum loaded" );
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::setPeakFitFor(...): invalid index" );
  
  const int row = index.row();
  PeakShrdPtr old_peak = m_sortedPeaks.at( row );
  std::deque< PeakShrdPtr >::iterator energy_pos
                          = find( m_peaks->begin(), m_peaks->end(), old_peak );
  if( energy_pos == m_peaks->end() || (*energy_pos) != old_peak )
    throw runtime_error( "Serious logic error in setPeakFitFor(...)" );
  const size_t energy_index = energy_pos - m_peaks->begin();
  
  PeakDef newPeak = *old_peak;
  newPeak.setFitFor( coef, fitfor );
  
  m_sortedPeaks[row] = std::make_shared<PeakDef>( newPeak );
  (*m_peaks)[energy_index] = m_sortedPeaks[row];
  
  // We dont need to emit data changed, because the "Fit For" values
  //  are not displayed as part of the MVC values
  
  if( coef == PeakDef::CoefficientType::Mean )
  {
    // If we are fixing the mean, we wont allow using this peak for energy calibration,
    //  so we need to emit that this value may have changed.
    WModelIndex peak_index = PeakModel::index(row, kUseForCalibration);
    dataChanged().emit( peak_index, peak_index );
  }//if( coef == PeakDef::CoefficientType::Mean )
  
  notifySpecMeasOfPeakChange();
}//setPeakFitFor(...)


void PeakModel::setContinuumPolynomialFitFor( const Wt::WModelIndex index,
                                              size_t polyCoefNum,
                                              const bool fitfor )
{
  if( !m_peaks )
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...): no"
                              " foreground spectrum loaded" );
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...):"
                              " invalid index" );
  
  const int row = index.row();
  PeakShrdPtr old_peak = m_sortedPeaks.at( row );
  
  if( !old_peak )
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...):"
                             " invalid peak" );
  
  std::shared_ptr<const PeakContinuum> old_cont = old_peak->continuum();
  
  assert( old_cont );
  if( !old_cont )// shouldnt happen, but JIC
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...):"
                             " invalid continuum" );
  
  std::shared_ptr<PeakContinuum> new_cont = make_shared<PeakContinuum>( *old_cont );
  
  const bool madeChange = new_cont->setPolynomialCoefFitFor( polyCoefNum, fitfor );
  if( !madeChange )
    throw std::runtime_error( "PeakModel::setContinuumPolynomialFitFor(...):"
                             " invalid coefficient number" );
  
  // If we create a new peak, using the copy constructor of the peak to be modified, the _old_
  //  continuum will be used.  This is a dangling issue of a poor design.  So we will manually
  //  create a new continuum, and make new peaks to replace all the old peaks that shared this
  //  continuum.
  
  for( PeakShrdPtr &p : m_sortedPeaks )
  {
    if( p && (p->continuum() == old_cont) )
    {
      auto new_peak = std::make_shared<PeakDef>( *old_peak );
      new_peak->setContinuum( new_cont );
      const auto pos = std::find( std::begin(*m_peaks), std::end(*m_peaks), p );
      
      assert( pos != std::end(*m_peaks) );
      if( pos == std::end(*m_peaks) )
        throw runtime_error( "PeakModel::setContinuumPolynomialFitFor(...):"
                            "Failed to find peak in m_peaks" );
      
      // We will set the shared ptrs in m_sortedPeaks and m_peaks to equal the newly created pointer
      *pos = new_peak;
      p = new_peak;
    }//if( this peaks continuum was old_cont )
  }//for( const auto &p : m_sortedPeaks )
  
  notifySpecMeasOfPeakChange();
  
  //dataChanged().emit( index, PeakModel::index(...) );
}//void setContinuumPolynomialFitFor(...)


//Functions for the Wt::WAbstractItemModel interface
int PeakModel::rowCount( const WModelIndex & index) const
{
  //Returns # of children.  Root = size(), children = 0
  if (!index.isValid())
  { //root
    if( !m_peaks )
        return 0;
    else
      return static_cast<int>( m_peaks->size() );
  } //root
  else
  { //children
    return 0;
  } //children

}//int rowCount( const WModelIndex & ) const


int PeakModel::columnCount( const WModelIndex & ) const
{
  return PeakModel::kNumColumns;
}//int columnCount( const WModelIndex & ) const


WModelIndex PeakModel::parent( const Wt::WModelIndex & ) const
{
  return WModelIndex();
}//WModelIndex parent( const Wt::WModelIndex & ) const


const PeakModel::PeakShrdPtr &PeakModel::peak( const Wt::WModelIndex &index ) const
{
  if( !index.isValid() )
    throw std::runtime_error( "PeakModel::peak(WModelIndex): invalid input index" );

  if( index.row() < 0 || index.row() >= static_cast<int>(m_sortedPeaks.size()) )
    throw runtime_error( "PeakModel::peak(WModelIndex): requested index out of range" );
  
  return m_sortedPeaks[index.row()];
}//PeakShrdPtr peak( const Wt::WModelIndex &index ) const

boost::any PeakModel::data( const WModelIndex &index, int role ) const
{
  using boost::any;
  
  if( !m_peaks )
    return any();

  if( m_sortedPeaks.size() != m_peaks->size() )
  {
    string msg = "PeakModel::data(...)\n\tm_sortedPeaks.size()="
                  + std::to_string(m_sortedPeaks.size())
                  + ", m_peaks->size()=" + std::to_string(m_peaks->size());
    cerr << endl << msg << endl << endl;
    throw std::runtime_error( msg );
  }//if( m_sortedPeaks.size() != m_peaks->size() )

  // TODO: implement ToolTipRole
  if( role == ToolTipRole )
  {
    if( index.column() == kUseForCalibration )
      return WString::tr("pm-tt-use-for-cal");
    return any();
  }//if( role == ToolTipRole )
  
  if( (role != Wt::DisplayRole)
     && (role != Wt::EditRole)
     && !((role==Wt::CheckStateRole)
            && ((index.column()==kUseForCalibration)
            || (index.column()==kUseForShieldingSourceFit)
            || (index.column()==kUseForManualRelEff))
          )
      )
  {
    return any();
  }


  assert( m_sortedPeaks.size() == m_peaks->size() );

  const int row = index.row();
  const int column = index.column();
  const int nrows = static_cast<int>( m_peaks->size() );

  if( row<0 || column<0 || column>=kNumColumns || row>=nrows )
    return any();


  const PeakShrdPtr &peak = m_sortedPeaks[row];

  
  auto getPeakArea = [&peak,this]() -> boost::any {
    switch( peak->type() )
    {
      case PeakDef::GaussianDefined:
        return peak->peakArea();
        
      case PeakDef::DataDefined:
      {
        if( !m_foreground )
          return boost::any();
        
        double contArea = 0.0;
        const std::shared_ptr<const SpecUtils::Measurement> &dataH = m_foreground;
        
        if( !dataH )
          return boost::any();
        
        assert( peak->continuum()->energyRangeDefined() );
        const double lowx = peak->lowerX();
        const double upperx = peak->upperX();
        
        assert( peak->continuum()->parametersProbablySet() );
        if( peak->continuum()->parametersProbablySet() )
          contArea = peak->offset_integral( lowx, upperx, dataH );
        
        
        const size_t lower_channel = dataH->find_gamma_channel( lowx );
        const size_t upper_channel = dataH->find_gamma_channel( upperx - 0.00001 );
        const double dataArea = dataH->gamma_channels_sum( lower_channel, upper_channel );
        
        return (dataArea - contArea);
      }//case PeakDef::DataDefined:
    }//switch( peak->type() )

    return boost::any();
  };//auto getPeakArea lambda
  
  switch( column )
  {
    case kMean:
      return peak->mean();
      
    case kFwhm:
      switch( peak->type() )
      {
        case PeakDef::GaussianDefined:
          return 2.3548201*peak->sigma();
          
        case PeakDef::DataDefined:
          return boost::any();
      }//switch( peak->type() )
      
    case kAmplitude:
    {
      boost::any areaAny = getPeakArea();
      if( areaAny.empty() )
        return areaAny;
      
      const double area = boost::any_cast<double>(areaAny);
      double uncert = peak->amplitudeUncert();
      
      switch( peak->type() )
      {
        case PeakDef::GaussianDefined:
          break;
        case PeakDef::DataDefined:
          uncert = -1.0; // JIC
          break;
      }//switch( peak->type() )
      
      if( uncert <= 0.0 )
      {
        char text[64];
        snprintf( text, sizeof(text), "%.2f", area );
        return WString::fromUTF8(text);
      }
      
      // TODO: Figure out how many significant figures to show - this is kinda a guess for the moment
      const int numValLeft = 1 + static_cast<int>( std::floor( std::log10(area) ) );  //1.2 will give 1, 10.1 will give 2
      //const int numUncertLeft = 1 + static_cast<int>( std::floor( std::log10(uncert) ) ); //0.11 will give 0, 0.011 will give -1
        
      const int nsigfig = 1 + ((area > 1.0) ? std::min(std::max(3,numValLeft), 6) : 4); //we'll print out one place past decimal point
      const string txt = PhysicalUnits::printValueWithUncertainty( area, uncert, nsigfig );
      return WString::fromUTF8(txt);
    }//case kAmplitude:
      
    case PeakModel::kCps:
    {
      boost::any areaAny = getPeakArea();
      
      if( areaAny.empty() )
        return areaAny;
      
      const double area = boost::any_cast<double>(areaAny);
      const shared_ptr<const SpecUtils::Measurement> &dataH = m_foreground;
      const float liveTime = dataH ? dataH->live_time() : -1.0f;
      if( liveTime <= 0.0f )
        return boost::any();
      
      const double cps = area / liveTime;
      const double uncert = peak->amplitudeUncert();
      
      char text[64];
      if( (peak->type() != PeakDef::GaussianDefined) || (uncert <= 0.0) )
      {
        snprintf( text, sizeof(text), "%.4g", cps );
        return WString(text);
      }
      const double cpsUncert = uncert / liveTime;
      const string txt = PhysicalUnits::printValueWithUncertainty( cps, cpsUncert, 4 );
      return WString::fromUTF8(txt);
    }//case PeakModel::kCps:
      
    case kIsotope:
    {
      const SandiaDecay::Nuclide *parentNuclide = peak->parentNuclide();

      if( parentNuclide )
        return WString( parentNuclide->symbol );
      else if( peak->xrayElement() )
        return WString( peak->xrayElement()->symbol + " xray" );
      else if( peak->reaction() )
        return WString( peak->reaction()->name() );
      
      return boost::any();
    }//case kIsotope:

    case kDifference:
    {
      if( !peak->hasSourceGammaAssigned() )
        return boost::any();
      
      try
      {
        const float energy = peak->gammaParticleEnergy() - peak->mean();
        char text[32];
        snprintf( text, sizeof(text), "%.2f keV", (floor(100.0*energy + 0.5)/100.0) );
        return WString( text );
      }catch( std::exception & )
      {
        return boost::any();
      }
    }
      
    case kPhotoPeakEnergy:
    {
      if( !peak->hasSourceGammaAssigned() )
        return boost::any();
      
      try
      {
        float energy = peak->gammaParticleEnergy();
        const char *prefix = "", *postfix = "";
        
        switch( peak->sourceGammaType() )
        {
          case PeakDef::NormalGamma:
          case PeakDef::AnnihilationGamma:
            break;
          case PeakDef::SingleEscapeGamma:
            prefix = "S.E. ";
            if( peak->decayParticle() )
              energy = peak->decayParticle()->energy;
          break;
          case PeakDef::DoubleEscapeGamma:
            prefix = "D.E. ";
            if( peak->decayParticle() )
              energy = peak->decayParticle()->energy;
          break;
          case PeakDef::XrayGamma:
            postfix = " xray";
          break;
        }//switch( peak->sourceGammaType() )
        
        char text[32];
        snprintf( text, sizeof(text), "%s%.2f keV%s", prefix, (floor(100.0*energy + 0.5)/100.0), postfix );
        return WString( text );
      }catch( std::exception & )
      {
        return boost::any();
      }
    }//case kPhotoPeakEnergy:

    case kUseForShieldingSourceFit:
    {
      // Make so we only "use for shielding/source fit" show checkbox for decay gammas and x-rays,
      //  and not for peaks with no source associated, or with florescence x-rays, or for reactions
      
      switch( peak->sourceGammaType() )
      {
        case PeakDef::NormalGamma:
        case PeakDef::XrayGamma:
          if( !peak->nuclearTransition() || !peak->parentNuclide() || (peak->decayParticleIndex() < 0) )
            return boost::any();
          break;
          
        case PeakDef::AnnihilationGamma:
          // Annihilation gammas wont have a nuclearTransition or decay particle index associated with them
          if( !peak->parentNuclide() )
            return boost::any();
          break;
        
        case PeakDef::SingleEscapeGamma:
        case PeakDef::DoubleEscapeGamma:
          // Don't show "use for shielding/source fit" checkbox for single and double escape peaks
          return boost::any();
          break;
      }//switch( peak->sourceGammaType() )
      
      return peak->useForShieldingSourceFit();
    }//case kUseForShieldingSourceFit:

    case kUseForCalibration:
    {
      const bool fixed_mean = !peak->fitFor(PeakDef::CoefficientType::Mean);
      if( fixed_mean )
        return boost::any();
      return peak->useForEnergyCalibration();
    }
      
    case kUseForManualRelEff:
    {
      switch( peak->sourceGammaType() )
      {
        case PeakDef::XrayGamma:
          if( peak->parentNuclide() )
            return peak->useForManualRelEff();
          return boost::any();
        break;

        case PeakDef::NormalGamma:
        case PeakDef::AnnihilationGamma:
          if( !peak->parentNuclide() && !peak->reaction() )
            return boost::any();
          break;

        case PeakDef::SingleEscapeGamma:
        case PeakDef::DoubleEscapeGamma:
          return boost::any();
          break;
      }//switch( peak->sourceGammaType() )
      
      return peak->useForManualRelEff();
    }//case kUseForManualRelEff:
      
    case kPeakLineColor:
    {
      if( peak->lineColor().isDefault() )
        return boost::any();
      return boost::any( WString::fromUTF8(peak->lineColor().cssText(false)) );
    }//case kPeakLineColor:
      
    case kUserLabel:
    {
      if( peak->userLabel().empty() )
        return boost::any();
      return boost::any( WString::fromUTF8(peak->userLabel()) );
    }//case kUserLabel:
      
    case kHasSkew:
    {
      return WString::tr( (peak->skewType() == PeakDef::NoSkew) ? "False" : "True" );
    }//case kHasSkew:
      
    case kSkewAmount:
    {
      char text[64] = { '\0' };
      
      switch( peak->skewType() )
      {
        case PeakDef::NoSkew:
          snprintf( text, sizeof(text), "NA" );
          break;
          
        case PeakDef::Bortel:
        case PeakDef::GaussExp:
          snprintf( text, sizeof(text), "%.3f", peak->coefficient(PeakDef::SkewPar0) );
          break;
          
        case PeakDef::CrystalBall:
        case PeakDef::ExpGaussExp:
          snprintf( text, sizeof(text), "%.3f, %.3f",
                   peak->coefficient(PeakDef::SkewPar0),
                   peak->coefficient(PeakDef::SkewPar1) );
          break;
          
        case PeakDef::DoubleSidedCrystalBall:
          snprintf( text, sizeof(text), "%.3f, %.3f, %.3f, %.4f",
                    peak->coefficient(PeakDef::SkewPar0),
                    peak->coefficient(PeakDef::SkewPar1),
                    peak->coefficient(PeakDef::SkewPar2),
                   peak->coefficient(PeakDef::SkewPar3) );
          break;
      }//switch( peak->skewType() )
      
      return WString::fromUTF8( text );
    }//case kSkewAmount:
      
    case kType:
    {
      if( peak->gausPeak() )
        return WString::tr( "pm-gaussian" );
      else
        return WString::tr( "pm-region" );
    }
      
    case kLowerX:
    case kUpperX:
    {
      shared_ptr<const PeakContinuum> continuum = peak->continuum();
      assert( continuum->energyRangeDefined() );
      return (column == kLowerX) ? peak->lowerX() : peak->upperX();
    }//case kLowerX / case kUpperX:
    
      
    case kRoiCounts:
    {
      if( !m_foreground )
        return boost::any();
      
      assert( peak->continuum()->energyRangeDefined() );
      const double lowx = peak->lowerX();
      const double upperx = peak->upperX();
      
      return m_foreground->gamma_integral( lowx, upperx );
    }//case kRoiCounts:
      
    case kContinuumType:
      return WString::tr( PeakContinuum::offset_type_label_tr(peak->continuum()->type()) );
      
    case kNumColumns:
      return any();
  }//switch( section )

  return any();
}//boost::any PeakModel::data( const WModelIndex &index, int role ) const



PeakModel::SetGammaSource PeakModel::setNuclide( PeakDef &peak,
                                 const PeakDef::SourceGammaType src_type,
                                 const SandiaDecay::Nuclide * const nuclide,
                                 const double ref_energy,
                                 const double nsigma_window )
{
  const bool hadSource = (peak.parentNuclide() || peak.xrayElement() || peak.reaction() );
  assert( !nuclide || (ref_energy > 1.0) );
  
  size_t transition_index = 0;
  const SandiaDecay::Transition *transition = nullptr;
  
  // If width is negative, then nearest photopeak by energy will be used.
  const double width = peak.gausPeak() ? nsigma_window*peak.sigma() : (nsigma_window/8.0)*peak.roiWidth();
  
  const bool xrayOnly = (src_type == PeakDef::SourceGammaType::XrayGamma);
  PeakDef::SourceGammaType sourceGammaType = src_type;
  
  PeakDef::findNearestPhotopeak( nuclide, ref_energy, width,
                                xrayOnly, -1.0, transition, transition_index, sourceGammaType );
  
  //There wasnt any photopeaks within 4 sigma, so instead we'll just use
  //  the closest photpopeak
  if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma) && (nsigma_window>=0.0) )
    PeakDef::findNearestPhotopeak( nuclide, ref_energy, -1.0, xrayOnly, -1.0,
                                  transition, transition_index, sourceGammaType );
  
  switch( src_type )
  {
    case PeakDef::NormalGamma:
    case PeakDef::AnnihilationGamma:
    case PeakDef::XrayGamma:
      break;
      
    case PeakDef::SingleEscapeGamma: case PeakDef::DoubleEscapeGamma:
      if( transition )
        sourceGammaType = src_type;
      else if( sourceGammaType == PeakDef::AnnihilationGamma )
        transition = nullptr;  //dont expect this to ever happen
      break;
  }//switch( src_type )
  
  
  if( !transition && (sourceGammaType != PeakDef::AnnihilationGamma) )
  {
    peak.clearSources();
    
    return (hadSource ? SourceChange : NoSourceChange);
  }
  
  
  bool changedFit = false, shouldFit = false, shouldUseForRe = false;
  
  switch( sourceGammaType )
  {
    case PeakDef::NormalGamma:
      shouldFit = recommendUseForFit( nuclide, transition->products[transition_index].energy );
      shouldUseForRe = recommendUseForManualRelEff( nuclide, transition->products[transition_index].energy );
      break;
      
    case PeakDef::AnnihilationGamma:
      shouldFit = recommendUseForFit( nuclide, 510.99891f );
      break;
      
    case PeakDef::SingleEscapeGamma:
    case PeakDef::DoubleEscapeGamma:
    case PeakDef::XrayGamma:
      shouldFit = false;
      break;
  }//switch( src_type )
  
  
  changedFit |= (shouldFit != peak.useForShieldingSourceFit());
  changedFit |= (shouldFit && !peak.nuclearTransition());
  changedFit |= (shouldUseForRe != peak.useForManualRelEff());
  changedFit |= (shouldUseForRe && !peak.nuclearTransition());
  
  
  peak.useForShieldingSourceFit( shouldFit );
  peak.useForManualRelEff( shouldUseForRe );
  peak.setNuclearTransition( nuclide, transition, int(transition_index), sourceGammaType );
  
  return (changedFit ? SourceAndUseChanged : SourceChange);
}//setNuclide(...)


PeakModel::SetGammaSource PeakModel::setXray( PeakDef &peak,
                              const SandiaDecay::Element *el,
                              const double ref_energy )
{
  const SandiaDecay::Element * const prevEl = peak.xrayElement();
  const bool hadSource = (peak.parentNuclide() || prevEl || peak.reaction() );
  assert( !el || (ref_energy > 1.0) );
  
  if( !el )
  {
    peak.clearSources();
    return (hadSource ? SourceChange : NoSourceChange);
  }//if( !el )
  
  double xray_energy = ref_energy;
  const SandiaDecay::EnergyIntensityPair *nearXray = PeakDef::findNearestXray( el, ref_energy );
  xray_energy = (nearXray ? nearXray->energy : 0.0);
  
  const bool isSame = ((prevEl == el) && (fabs(peak.xrayEnergy() - xray_energy) < 0.0001));
  if( isSame )
    return NoSourceChange;
  
  peak.setXray( el, xray_energy );
  
  if( !nearXray && !hadSource )
    return NoSourceChange;
  
  if( peak.useForShieldingSourceFit() || peak.useForManualRelEff() )
  {
    peak.useForShieldingSourceFit( false );
    peak.useForManualRelEff( false );
    return SourceAndUseChanged;
  }
  
  return SourceChange;
}//setXray(...)


PeakModel::SetGammaSource PeakModel::setReaction( PeakDef &peak,
                                  std::string label,
                                  const PeakDef::SourceGammaType src_type,
                                  const double ref_energy,
                                  const double windowHalfWidth )
{
  const bool hadSource = (peak.parentNuclide() || peak.xrayElement() || peak.reaction() );
  assert( label.empty() || (ref_energy > 1.0) );
  
  const ReactionGamma *rctndb = ReactionGammaServer::database();
  string rctstr = label;
  size_t pos = label.find_first_of( ')' );
  if( pos != string::npos )
    rctstr = label.substr( 0, pos+1 );
  
  vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
  try
  {
    rctndb->gammas( rctstr, possible_rctns );
    if( possible_rctns.empty() )
    {
      peak.clearSources();
      return (hadSource ? SourceChange : NoSourceChange);
    }
  }catch(...)
  {
    peak.clearSources();
    return (hadSource ? SourceChange : NoSourceChange);
  }
  
  // TODO: just taking first reaction, however there may be multiple
  const ReactionGamma::Reaction *rctn = possible_rctns[0].reaction;
  
  if( !rctn )
  {
    peak.clearSources();
    return (hadSource ? SourceChange : NoSourceChange);
  }
  
  double best_delta_e = std::numeric_limits<double>::max(), nearestE = 0.0;
  for( const ReactionGamma::Reaction::EnergyYield &eip : rctn->gammas )
  {
    if( eip.abundance <= 0.0 )
      continue;
    
    const double delta_e = fabs( eip.energy - ref_energy );
    double scaleDeltaE = (0.1*windowHalfWidth + delta_e) / eip.abundance;
    if( windowHalfWidth <= 0.0 )
      scaleDeltaE = delta_e;
    
    if( scaleDeltaE <= best_delta_e )
    {
      nearestE = eip.energy;
      best_delta_e = scaleDeltaE;
    }
  }//for( const ReactionGamma::EnergyAbundance &eip : rctn->gammas )
  
  if( nearestE == 0.0 )
  {
    peak.clearSources();
    return (hadSource ? SourceChange : NoSourceChange);
  }
  
  peak.setReaction( rctn, static_cast<float>(nearestE), src_type );
  
  return (rctn ? SourceChange : NoSourceChange);
}//setReaction


PeakModel::SetGammaSource PeakModel::setNuclideXrayReaction( PeakDef &peak,
                                                            std::string label,
                                                    const double nsigma_window )
{  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  if( !db )
    return NoSourceChange;
  
  SpecUtils::to_lower_ascii( label );
  
  const bool hadSource = (peak.parentNuclide() || peak.xrayElement() || peak.reaction() );
  
  if( label.empty()
     || SpecUtils::contains( label, "none" )
     || SpecUtils::contains( label, "undef" )
     || label == "na" )
  {
    const bool use = (peak.useForShieldingSourceFit() || peak.useForManualRelEff());
    peak.clearSources();
    
    if( !hadSource )
      return NoSourceChange;
    
    if( use )
      return SourceAndUseChanged;
    
    return SourceChange;
  }//if( getting rid of source )
  
  
  PeakDef::SourceGammaType srcType;
  PeakDef::gammaTypeFromUserInput( label, srcType );
  
  // This next line removes the reference energy text from `label`
  double ref_energy = PeakDef::extract_energy_from_peak_source_string( label );
  if( ref_energy < std::numeric_limits<float>::epsilon() )
  {
    ref_energy = peak.mean();
    
    double extraEnergy = 0.0;
    switch( srcType )
    {
      case PeakDef::NormalGamma:
      case PeakDef::AnnihilationGamma:
      case PeakDef::XrayGamma:
        break;
        
      case PeakDef::SingleEscapeGamma:
        ref_energy += 510.99891;
        break;
        
      case PeakDef::DoubleEscapeGamma:
        ref_energy += 2.0*510.99891;
        break;
    }//switch( src_type )
  }//if( input string didnt contain energy )
  
  
   string nuclabel = label;

  {//begin manipulation of nuclabel
    // We are assuming the first set of numbers is atomic number, but
    // then there may be more numbers that indicate meta-level, or 
    // they could just be something we dont care about.
    // Example inputs might be:
    //   "hf178m2 574.219971 kev", or "hf178m2 , I=1.8E-6%", 
    //   "Hf178-meta2", "Hf178-meta-2", "Hf178 meta-2", "Hf178 meta 2", etc
    // 
    //  However, note that the energy string *should* have already been removed
    //  by the call to extract_energy_from_peak_source_string(...)

    //SandiaDecay doesnt care about spacing, so get rid of spaces, and similar
    //SpecUtils::ireplace_all( nuclabel, " \t-\n\r", "" );

    size_t number_start = nuclabel.find_first_of( "0123456789" );
    size_t number_stop = string::npos;
    if( number_start != string::npos )
      number_stop = nuclabel.find_first_not_of( "0123456789", number_start );

    if( number_stop != string::npos )
    {
      vector<string> after_words;
      const string afterstr = nuclabel.substr( number_stop );
      SpecUtils::split( after_words, afterstr, " \t-\n\r" );
      nuclabel = nuclabel.substr( 0, number_stop );
      
      string meta = "";
      if( after_words.size() >= 1 )
      {
        const string &first_word = after_words[0];
        if( SpecUtils::iequals_ascii( first_word, "m" )
          || SpecUtils::iequals_ascii( first_word, "meta" ) )
        {
          nuclabel += "m";
          if( after_words.size() >= 2 )
          {
            const string &second_word = after_words[1];
            if( (second_word == "1") || (second_word == "2") || (second_word == "3") )
              nuclabel += second_word;
          }
        }else if( SpecUtils::iequals_ascii( first_word, "m1" ) 
          || SpecUtils::iequals_ascii( first_word, "m2" )
          || SpecUtils::iequals_ascii( first_word, "m3" ) 
          || SpecUtils::iequals_ascii( first_word, "meta1" )
          || SpecUtils::iequals_ascii( first_word, "meta2" )
          || SpecUtils::iequals_ascii( first_word, "meta3" )
          )
        {
          nuclabel += first_word;
        }
      }//if( after_words.size() >= 1 )
    }//if( number_stop != string::npos )
  }//end manipulation of nuclabel

  const SandiaDecay::Nuclide *nuclide = db->nuclide( nuclabel );
  if( nuclide )
    return setNuclide( peak, srcType, nuclide, ref_energy, nsigma_window );
  
  //lets check for a reaction
  const size_t paren_pos = label.find_first_of( ')' );
  if( paren_pos != string::npos )
  {
    const string rctstr = SpecUtils::trim_copy( label.substr( 0, paren_pos+1 ) );
    return PeakModel::setReaction( peak, rctstr, srcType, ref_energy, nsigma_window );
  }//if( paren_pos != string::npos )
  
  
  bool isXray = (srcType == PeakDef::SourceGammaType::XrayGamma);
  if( !isXray )
  {
    SpecUtils::trim( label ); //jic, but probably already fine
    
    isXray = true;
    for( size_t i = 0; isXray && (i < label.size()); ++i )
      isXray = (label[i]>='a' && label[i]<='z');
  }//if( !isXray )
  
  if( isXray )
  {
    //Example input that might get here is "Pb xray 84.9 kev", which will have been transformed
    //  to "pb 98.2 kev" by PeakDef::gammaTypeFromUserInput, and
    //  PeakDef::extract_energy_from_peak_source_string would have transformed into
    //  "pb"
    const SandiaDecay::Element *el = db->element( label );
    
    // We'll require the label to now _only_ be either the element label (ex Pb), or name (ex Lead).
    if( el && (SpecUtils::iequals_ascii(el->name, label)
              || SpecUtils::iequals_ascii(el->symbol, label)) )
    {
      return setXray( peak, el, ref_energy );
    }
  }//if( srcType == PeakDef::SourceGammaType::XrayGamma )
  
  
  peak.clearSources();
  return (hadSource ? SourceChange : NoSourceChange);
}//bool PeakModel::setNuclideXrayReaction( PeakDef &peak, std::string )


bool PeakModel::setData( const WModelIndex &index,
                         const boost::any &value, int role )
{
  if( !m_peaks )
    throw runtime_error( "Set a primary spectrum before setting peak data" );

  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  notifySpecMeasOfPeakChange();
  
  try
  {
    if( (role != Wt::EditRole) && (role != Wt::CheckStateRole) )
      return false;

    const int row = index.row();
    const int column = index.column();
    const int nrows = static_cast<int>( m_peaks->size() );

    bool changedFit = false;
    
    if( value.empty() || row < 0 || column < 0 || column>=kNumColumns || row >= nrows )
      return false;

    switch( column )
    {
      case kMean:
      case kIsotope:
      case kPhotoPeakEnergy:
      case kCandidateIsotopes:
      case kUseForCalibration:
      case kUseForManualRelEff:
      case kUseForShieldingSourceFit:
      case kUserLabel:
      case kPeakLineColor:
      break;
        
      case kFwhm:
      case kAmplitude:
      {
        switch( m_sortedPeaks[row]->type() )
        {
          case PeakDef::GaussianDefined:
            break;
          case PeakDef::DataDefined:
            cerr << "PeakModel::setData(...)\n\tCant set area or FWHM for non-Gaussian peak" << endl;
            return false;
        }//switch( m_sortedPeaks[row]->type() )
        
        break;
      }//case kFwhm or kAmplitude

      case kCps: case kHasSkew: case kSkewAmount: case kType: case kLowerX: case kUpperX:
      case kRoiCounts: case kContinuumType: case kNumColumns: case kDifference:
      default:
        cerr << "PeakModel::setData(...)\n\tUn Supported column" << endl;
        return false;
    }//switch( section )

    double dbl_val = 0.0, uncert_val = -1.0;
    WString txt_val;

    try
    {
      txt_val = boost::any_cast<WString>( value );
    }catch(...)
    {}

    switch( column )
    {
      case kMean: case kFwhm: case kAmplitude:
        try
        {
          const string strval = txt_val.toUTF8();
          
          const auto parseValWithUncert = []( const string &input, string &value, string &uncert ){
            string::size_type pos = input.find( "\xC2\xB1" );
            if( pos == string::npos )
              pos = input.find( "+-" );
            if( pos == string::npos )
              pos = input.find( "-+" );
            if( pos == string::npos )
            {
              value = input;
              SpecUtils::trim(value);
              return;
            }
            
            value = input.substr(0, pos);
            uncert = input.substr(pos + 2); //All "+-" strings are two bytes long
            SpecUtils::trim(value);
            SpecUtils::trim(uncert);
          };
          
          string valstr, uncertstr;
          parseValWithUncert( strval, valstr, uncertstr );
          
          dbl_val = std::stod( valstr );
          if( !uncertstr.empty() )
            uncert_val = std::stod( uncertstr );
          
          // There could be some rounding in the string representation of the uncertainty, or the
          //  area so lets avoid this if the user hasnt changed that part of it
          const boost::any prevData = data( index, role );
          if( !prevData.empty() )
          {
            WString prevDataWstr = boost::any_cast<WString>( prevData );
            
            assert( row < m_sortedPeaks.size() );
            const PeakShrdPtr &peak = m_sortedPeaks[row];
            
            string prevValStr, prevUncertStr;
            parseValWithUncert( prevDataWstr.toUTF8(), prevValStr, prevUncertStr );
            
            if( prevValStr == valstr )
            {
              if( column == kMean )
                dbl_val = peak->mean();
              else if( column == kFwhm )
                dbl_val = peak->fwhm();
              else if( column == kAmplitude && (peak->type() == PeakDef::GaussianDefined) )
                dbl_val = peak->amplitude();
            }//if( prevValStr == valstr )
            
            if( (uncert_val > 0.0) && (prevUncertStr == uncertstr) )
            {
              if( column == kMean )
                uncert_val = peak->meanUncert();
              else if( column == kFwhm )
                uncert_val = 2.3548201 * peak->sigmaUncert();
              else if( column == kAmplitude && (peak->type() == PeakDef::GaussianDefined) )
                uncert_val = peak->amplitudeUncert();
            }//if( prevValStr == valstr )
          }//if( prevData, so we can potentially avoid rounding values )
        }catch(...)
        {
          cerr << "PeakModel::setData(...)\n\tUnable to convert '" << txt_val
               << "' to a dbl" << endl;
          return false;
        }//try / catch
      break;
      default:
        break;
    }//switch( column )

    
    PeakShrdPtr old_peak = m_sortedPeaks[row];
    std::deque< PeakShrdPtr >::iterator energy_pos
                           = find( m_peaks->begin(), m_peaks->end(), old_peak );
    assert( energy_pos != m_peaks->end() );

    PeakDef new_peak( *old_peak );

    switch( column )
    {
      case kMean:
        switch( new_peak.type() )
        {
          case PeakDef::GaussianDefined:
            new_peak.setMean( dbl_val );
            if( uncert_val > 0.0 )
              new_peak.setMeanUncert( uncert_val );
          break;
            
          case PeakDef::DataDefined:
            return false;
          break;
        }//switch( new_peak.type() )
      break;

      case kFwhm:
        switch( new_peak.type() )
        {
          case PeakDef::GaussianDefined:
            new_peak.setSigma( dbl_val / 2.3548201 );
            if( uncert_val > 0.0 )
              new_peak.setSigmaUncert( uncert_val / 2.3548201 );
          break;
          
          case PeakDef::DataDefined:
            return false;
          break;
        }//switch( new_peak.type() )
      break;

      case kAmplitude:
        switch( new_peak.type() )
        {
          case PeakDef::GaussianDefined:
            new_peak.setAmplitude( dbl_val );
            if( uncert_val > 0.0 )
              new_peak.setAmplitudeUncert( uncert_val );
          break;
          
          case PeakDef::DataDefined:
            return false;
          break;
        }//switch( new_peak.type() )
      break;

      case kIsotope:
      {
        const std::string srctxt = txt_val.toUTF8();
        const SetGammaSource result
                   = setNuclideXrayReaction( new_peak, srctxt, 4.0 );
        changedFit |= (result==SourceAndUseChanged);
        
        //Lambda to see if any peaks, other than the ignorePeak, that has the
        //  same nuc/xray/rctn assigned as srcPeak has a color, and if so
        //  return that color.  If the source has more than one color of peaks
        //  then dont return a color.
        auto sameSrcColor = [this]( const PeakDef &srcPeak, const PeakShrdPtr &ignorePeak ) -> WColor {
          if( !m_peaks || !ignorePeak )
            return WColor();
          
          vector<WColor> src_colors;
          for( const auto &p : *m_peaks )
          {
            if( p
               && p != ignorePeak
               && (srcPeak.parentNuclide() || srcPeak.xrayElement() || srcPeak.reaction())
               && srcPeak.parentNuclide()==p->parentNuclide()
               && srcPeak.xrayElement()==p->xrayElement()
               && srcPeak.reaction()==p->reaction()
               && !p->lineColor().isDefault() )
            {
              if( std::find( begin(src_colors), end(src_colors),p->lineColor()) == end(src_colors) )
                src_colors.push_back( p->lineColor() );
            }
          }//for( const auto &p : *m_peaks )
          
          return src_colors.size()==1 ? *begin(src_colors) : WColor();
        };//sameSrcColor lambda
        
        
        if( !new_peak.hasSourceGammaAssigned() && old_peak->hasSourceGammaAssigned()
            && (result != NoSourceChange)  )
        {
          if( old_peak->lineColor() == sameSrcColor(*old_peak,old_peak) )
            new_peak.setLineColor( WColor() );
        }else if( (result != NoSourceChange) )
        {
          const auto oldsrccolor = sameSrcColor(*old_peak,old_peak);
          const auto newcolor = sameSrcColor(new_peak,old_peak);
          
          //If the old peak had a different color than its source - then dont change the color.
          
          if( (!oldsrccolor.isDefault() && oldsrccolor==old_peak->lineColor()) || !newcolor.isDefault() )
          {
            if( oldsrccolor.isDefault() || old_peak->lineColor().isDefault() || oldsrccolor==old_peak->lineColor() )
              new_peak.setLineColor( newcolor );
          }
        }
        
        break;
      }//case kIsotope:

      case kPhotoPeakEnergy:
      {
        const SandiaDecay::Transition *trans = new_peak.nuclearTransition();
        const SandiaDecay::Element *el = new_peak.xrayElement();
        const ReactionGamma::Reaction *rctn = new_peak.reaction();
        
        string text = SpecUtils::to_lower_ascii_copy( txt_val.narrow() );

        PeakDef::SourceGammaType srcType;
        PeakDef::gammaTypeFromUserInput( text, srcType );
        
        double unit = 1.0;
        if( text.find("mev") != string::npos )
          unit = 1000.0;
        else if( (text.find("kev") == string::npos)
                 && (text.find("ev") != string::npos) )
          unit = 0.001;

        SpecUtils::trim( text );

        double energy = -999.0;
        if( !(stringstream(text) >> energy) )
        {
          switch( srcType )
          {
            case PeakDef::NormalGamma:
            case PeakDef::XrayGamma:
              return false;
              break;
           
            case PeakDef::AnnihilationGamma:
              energy = 511.0;
              break;
              
            case PeakDef::SingleEscapeGamma:
              energy = new_peak.mean() + 511.0;
              break;
              
            case PeakDef::DoubleEscapeGamma:
              energy = new_peak.mean() + 1022.0;
              break;
          }//switch( srcType )
        }//if( !(convertstr >> energy) )
        
        energy *= unit;
        
        if( trans && trans->parent )
        {
          size_t trans_index = 0;
          const SandiaDecay::Nuclide *nuclide = new_peak.parentNuclide();
          const SandiaDecay::Transition *transition = NULL;
        
          //The 0.0 below means find the actual closest in energy, and not the
          //  most likely
          PeakDef::SourceGammaType sourceGammaType;
          const bool xrayOnly = (srcType == PeakDef::SourceGammaType::XrayGamma);
          PeakDef::findNearestPhotopeak( nuclide, energy, 0.0, xrayOnly, -1.0,
                                transition, trans_index, sourceGammaType );
          
          if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma) )
            return false;
          
          switch( srcType )
          {
            case PeakDef::NormalGamma:
            case PeakDef::AnnihilationGamma:
            case PeakDef::XrayGamma:
              break;
            case PeakDef::SingleEscapeGamma:
            case PeakDef::DoubleEscapeGamma:
              if( energy > 1022.0 )
                sourceGammaType = srcType;
              break;
          }//switch( srcType )
          
          new_peak.setNuclearTransition( new_peak.parentNuclide(), transition,
                                         int(trans_index), sourceGammaType );
          
          bool shouldFit = false, shouldUseForRe = false;
          switch( srcType )
          {
            case PeakDef::NormalGamma:
              shouldFit = recommendUseForFit( nuclide, transition->products[trans_index].energy );
              shouldUseForRe = recommendUseForManualRelEff( nuclide, transition->products[trans_index].energy );
              break;
            case PeakDef::AnnihilationGamma:
              shouldFit = recommendUseForFit( nuclide, 510.99891f );
              break;
            case PeakDef::SingleEscapeGamma:
            case PeakDef::DoubleEscapeGamma:
              break;
            case PeakDef::XrayGamma:
              shouldFit = false;
              break;
          }//switch( srcType )
          
          changedFit |= (shouldFit != old_peak->useForShieldingSourceFit());
          changedFit |= (shouldUseForRe != old_peak->useForManualRelEff());
          
          new_peak.useForShieldingSourceFit( shouldFit );
          new_peak.useForManualRelEff( shouldUseForRe );
        }else if( el )
        {
          double nearestE = -999.9;
          for( const SandiaDecay::EnergyIntensityPair &eip : el->xrays )
            if( fabs(eip.energy-energy) < fabs(nearestE-energy) )
              nearestE = eip.energy;
          
          if( fabs(nearestE-energy) < 5.0 )
            energy = nearestE;
          else
            return false;
          
          new_peak.setXray( el, energy );
        }else if( rctn )
        {
          double nearestE = -999.9;
          for( const ReactionGamma::Reaction::EnergyYield &eip : rctn->gammas )
            if( fabs(eip.energy-energy) < fabs(nearestE-energy) )
              nearestE = eip.energy;
          if( nearestE > 0.0 )
            new_peak.setReaction( rctn, static_cast<float>(nearestE), srcType );
        }else
        {
          return false;
        }
        
        break;
      }//case kPhotoPeakEnergy:

      case kCandidateIsotopes:
      {
        try
        {
          const vector<PeakDef::CandidateNuclide> candidates = boost::any_cast< vector<PeakDef::CandidateNuclide> >( value );
          new_peak.setCandidateNuclides( candidates );
        }catch(...)
        {
          return false;
        }

        break;
      }//case kCandidateIsotopes:

      case kUseForShieldingSourceFit:
      {
        try
        {
          const bool use = boost::any_cast<bool>( value );
          
//          cerr << "useForShieldingSourceFit=" << use << " type="
//               << value.type().name() << endl;
//          boost::any myfalsebool( bool(false) ), mytrueebool( bool(true) );
//          cerr << "myfalsebool value=" << boost::any_cast<bool>( myfalsebool)
//               << " mytrueebool value=" << boost::any_cast<bool>( mytrueebool)
//               << " and typename=" << myfalsebool.type().name() << endl;
          if( use && !new_peak.parentNuclide() )
            passMessage( WString::tr("pm-err-use-fit-no-nuc"), WarningWidget::WarningMsgHigh );
          
          if( new_peak.useForShieldingSourceFit() == use )
            return false;
          
          const SandiaDecay::RadParticle *radpart = new_peak.decayParticle();
          if( use && radpart && (radpart->type == SandiaDecay::XrayParticle) )
            passMessage( WString::tr("pm-warn-use-xray-fit"), WarningWidget::WarningMsgLow );
          
          new_peak.useForShieldingSourceFit( use );
        }catch(...)
        {
          cerr << "\n\tFailed in casting checkbox" << endl;
          return false;
        }

        break;
      }//case kUseForShieldingSourceFit:

      case kUseForCalibration:
      {
        try
        {
          const bool use = boost::any_cast<bool>( value );
          if( use && !new_peak.xrayElement() && !new_peak.parentNuclide()
              && !new_peak.reaction() )
            passMessage( WString::tr("pm-err-use-cal-no-nuc"), WarningWidget::WarningMsgHigh );
          
          if( use == new_peak.useForEnergyCalibration() )
            return false;
          new_peak.useForEnergyCalibration( use );
        }catch(...)
        {
          cerr << "\n\tFailed in casting checkbox" << endl;
          return false;
        }

        break;
      }//case kUseForCalibration:

         
      case kUseForManualRelEff:
      {
        try
        {
          const bool use = boost::any_cast<bool>( value );
          
          if( use )
          {
            const bool has_parent = (new_peak.parentNuclide() || new_peak.reaction());
            const bool is_gamma = ((new_peak.sourceGammaType() == PeakDef::SourceGammaType::NormalGamma)
                                   || (new_peak.sourceGammaType() == PeakDef::SourceGammaType::AnnihilationGamma) );
            const bool is_decay_xray = (new_peak.parentNuclide()
                                        && (new_peak.sourceGammaType() == PeakDef::SourceGammaType::XrayGamma));

            if( (!has_parent || !is_gamma) && !is_decay_xray )
               passMessage( WString::tr("pm-err-use-rel-act-no-nuc"), WarningWidget::WarningMsgHigh );
          }//if( use )
          
          if( use == new_peak.useForManualRelEff() )
            return false;
          new_peak.useForManualRelEff( use );
        }catch(...)
        {
          cerr << "\n\tFailed in casting checkbox" << endl;
          return false;
        }
        
        break;
      }//case kUseForManualRelEff:
        
        
      case kPeakLineColor:
      {
        const string css_color = txt_val.toUTF8();
        if( SpecUtils::iequals_ascii(css_color, "none")
            || SpecUtils::iequals_ascii(css_color, "na")
            || css_color.empty() )
          new_peak.setLineColor( Wt::WColor() );
        else
        {
          try{ new_peak.setLineColor( Wt::WColor(css_color) ); }catch(...){ return false; }
        }
        break;
      }//case kPeakLineColor:
        
      case kUserLabel:
      {
        new_peak.setUserLabel( txt_val.toUTF8() );
        break;
      }
        
      case kHasSkew: case kSkewAmount: case kType:
      case kLowerX: case kUpperX: case kRoiCounts:
      case kCps: case kContinuumType:
      case kNumColumns:
        return false;
    }//switch( section )

    m_sortedPeaks[row] = std::make_shared<const PeakDef>( new_peak );
    (*energy_pos) = m_sortedPeaks[row];

    if( column != kIsotope
        && column != kPhotoPeakEnergy
        && column != kMean 
        && column != kAmplitude )
    {
      dataChanged().emit( index, index );
    }else if( (column == kPhotoPeakEnergy) || (column == kAmplitude) )
    {
      if( changedFit )
        dataChanged().emit( PeakModel::index(row, kIsotope), PeakModel::index(row, kNumColumns-1) ); //just update whole row, jic
      else
        dataChanged().emit( index, PeakModel::index(row, kDifference) );
    }else if( column == kMean )
    {
      if( new_peak.decayParticle() || new_peak.parentNuclide() || new_peak.xrayElement() )
        dataChanged().emit( index, PeakModel::index(row, kPhotoPeakEnergy) );
      else
        dataChanged().emit( index, index );
    }else
    {
      dataChanged().emit( index, PeakModel::index(row, kUserLabel) );
    }


    if( column == m_sortColumn )
      sort( m_sortColumn, m_sortOrder );

    return true;
  }catch( std::exception &e )
  {
    cerr << "PeakModel::setData(...)\n\tUnexpected exception: " << e.what() << endl;
  }
  return false;
}//bool setData(...)


WFlags<ItemFlag> PeakModel::flags( const WModelIndex &index ) const
{
  switch( index.column() )
  {
    case kMean: case kUserLabel:
    case kIsotope: case kPhotoPeakEnergy:
#if( ALLOW_PEAK_COLOR_DELEGATE )
    case kPeakLineColor:
#endif
      return ItemIsEditable | ItemIsSelectable; //ItemIsSelectabl

    case kFwhm: case kAmplitude:
    {
      const int row = index.row();
      if( row < static_cast<int>(m_sortedPeaks.size()) )
      {
        switch( m_sortedPeaks[row]->type() )
        {
          case PeakDef::GaussianDefined:
            return ItemIsEditable | ItemIsSelectable;
          case PeakDef::DataDefined:
            return ItemIsSelectable;
        }//switch( peak type )
      }//if( a valid row of a peak )
      
      return ItemIsEditable | ItemIsSelectable;
    }//case kFwhm: case kAmplitude:
      
    case kUseForShieldingSourceFit:
    case kUseForCalibration:
    case kUseForManualRelEff:
      return ItemIsUserCheckable | ItemIsSelectable;

    case kHasSkew: case kSkewAmount: case kType:
    case kLowerX: case kUpperX: case kRoiCounts:
    case kContinuumType: case kCps:
#if( !ALLOW_PEAK_COLOR_DELEGATE )
    case kPeakLineColor:
#endif
    case kNumColumns:
      return ItemIsSelectable;
  }//switch( section )

  return WFlags<ItemFlag>();
}//WFlags<ItemFlag> flags( const WModelIndex &index )


WModelIndex PeakModel::index( int row, int column, const WModelIndex & ) const
{
  return  WAbstractItemModel::createIndex( row, column, (void *)0 );
}//WModelIndex PeakModel::index( int row, int column, const WModelIndex & ) const



boost::any PeakModel::headerData( int section, Orientation orientation, int role ) const
{
  //When orientation is Horizontal, section is a column number,
  //  when orientation is Vertical, section is a row (peak) number.

  if( role == LevelRole )
  {
    return 0;
  } //LevelRole
  else if (role == DisplayRole)
  {
    //If we are here, we want the column title
    switch( section )
    {
      case kMean:           return boost::any( WString::tr("Mean") );
      case kFwhm:           return boost::any( WString::tr("FWHM") ); //\x03C3
      case kAmplitude:      return boost::any( WString::tr("Area") );
      case kCps:            return boost::any( WString::tr("CPS") );
      case kIsotope:        return boost::any( WString::tr("Nuclide") );
      case kPhotoPeakEnergy:return boost::any( WString::tr("Photopeak") );
      case kDifference:     return boost::any( WString::tr("pm-hdr-diff") );
      case kUseForShieldingSourceFit: return boost::any( WString::tr("Use") );
      case kCandidateIsotopes:  return boost::any();
      case kUseForCalibration:  return boost::any( WString::tr("pm-hdr-cal-peak") );
      case kUseForManualRelEff: return boost::any( WString::tr("pm-hdr-rel-act") );
      case kUserLabel:      return boost::any( WString::tr("pm-hdr-label") );
      case kPeakLineColor:  return boost::any( WString::tr("pm-hdr-color") );
      case kHasSkew:        return boost::any( WString::tr("Skew") );
      case kSkewAmount:     return boost::any( WString::tr("pm-hdr-skew-amp") );
      case kType:           return boost::any( WString::tr("pm-hdr-peak-type") );
      case kLowerX:         return boost::any( WString::tr("pm-hdr-low-energy") );
      case kUpperX:         return boost::any( WString::tr("pm-hdr-up-energy") );
      case kRoiCounts:      return boost::any( WString::tr("pm-hdr-roi-counts") );
      case kContinuumType:  return boost::any( WString::tr("cont-type") );
      case kNumColumns:     return boost::any();
    }//switch( section )
  } //DisplayRole
  else if (role == ToolTipRole)
  {
    switch( section )
    {
      case kMean:           return boost::any( WString::tr("pm-hdr-tt-mean") );
      case kFwhm:           return boost::any( WString::tr("pm-hdr-tt-fwhm") ); //\x03C3
      case kAmplitude:      return boost::any( WString::tr("pm-hdr-tt-amp") );
      case kCps:            return boost::any( WString::tr("pm-hdr-tt-cps") );
      case kIsotope:        return boost::any();
      case kPhotoPeakEnergy:return boost::any( WString::tr("pm-hdr-tt-photopeak-energy") );
      case kDifference:     return boost::any( WString::tr("pm-hdr-tt-diff") );
      case kUseForShieldingSourceFit: return boost::any();
      case kCandidateIsotopes:  return boost::any();
      case kUseForCalibration:  return boost::any();
      case kUseForManualRelEff: return boost::any( WString::tr("pm-hdr-tt-man-rel-eff") );
      case kPeakLineColor:     return boost::any( WString::tr("pm-hdr-tt-color") );
      case kUserLabel:         return boost::any( WString::tr("pm-hdr-tt-label") );
      case kRoiCounts:         return boost::any( WString::tr("pm-hdr-tt-roi-counts") );
      case kHasSkew:
      case kSkewAmount:
      case kType:
      case kLowerX:
      case kUpperX:
      case kContinuumType:
      case kNumColumns:     return boost::any();
    }//switch( section )
    
  } //ToolTipRole
  else if( (orientation != Horizontal) || (role != DisplayRole) )
    return WAbstractItemModel::headerData( section, orientation, role );

  return boost::any();
}//any headerData( int section, Orientation orientation, int role ) const


bool PeakModel::removeRows( int row, int count, const WModelIndex & )
{
  if( !m_peaks )
    throw runtime_error( "Can not remove rows without a primary spectrum" );

  const int nrow = static_cast<int>( m_peaks->size() );

  if( !count )
    return false;

  if( (row >= nrow) || ((row+count) > nrow) || (count < 0) )
    throw runtime_error( "PeakModel::removeRows(...): invalid intput indexes" );

  beginRemoveRows( WModelIndex(), row, row+count-1 );

  deque< PeakShrdPtr >::iterator index, mean_index;
  const deque< PeakShrdPtr >::iterator start = m_sortedPeaks.begin() + row;
  const deque< PeakShrdPtr >::iterator end = start + count;

  for( index = start; index != end; ++index )
  {
    const PeakShrdPtr &peak = *index;
    mean_index = find( m_peaks->begin(), m_peaks->end(), peak );
    assert( mean_index != m_peaks->end() );
    m_peaks->erase( mean_index );
  }//for( index = start; index != end; ++index )

  m_sortedPeaks.erase( start, end );

  endRemoveRows();
  
  notifySpecMeasOfPeakChange();

  return true;
}//bool removeRows( int row, int count, const WModelIndex &)


//removeColumns(...) should not be used and will throw std::runtime_error
//  if called.
bool PeakModel::removeColumns( int, int, const Wt::WModelIndex & )
{
  throw std::runtime_error( "PeakModel::removeColumns(...): should not be called" );
  return false;
}//bool removeColumns( int, int, const Wt::WModelIndex & )


void PeakModel::sort( int col, Wt::SortOrder order )
{
  using boost::bind;

  if( m_sortedPeaks.empty() )
    return;

  m_sortOrder = order;
  m_sortColumn = Columns(col);

  if( m_sortColumn==kNumColumns )
    m_sortColumn = kMean;

  boost::function<bool(const PeakShrdPtr &, const PeakShrdPtr &)> sortfcn;
  const shared_ptr<const SpecUtils::Measurement> &data = m_foreground;
  sortfcn = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2,
                        m_sortColumn, order, data );

  layoutAboutToBeChanged().emit();
  stable_sort( m_sortedPeaks.begin(), m_sortedPeaks.end(), sortfcn );
  layoutChanged().emit();
  
  // I think this next part about saying the data changed is probably unnecessary; but didnt fully
  //  test, so leaving it in for now
  const int nrow = rowCount();
  const int ncol = columnCount();
  
  WModelIndex ul = index(0,0);
  WModelIndex lr = index(nrow-1,ncol-1);
  
  dataChanged().emit(ul,lr);
  
//XXX - I wouldnt think the next few lines should be necessary, but it appears
//      that it is
//  layoutChanged();
  
//  beginRemoveRows( WModelIndex(), 0, nrow-1 );
//  std::deque< PeakShrdPtr > peaks, sortedPeaks;
//  m_peaks->swap( peaks );
//  sortedPeaks.swap( m_sortedPeaks );
//  endRemoveRows();
//  beginInsertRows( WModelIndex(), 0, nrow-1 );
//  peaks.swap( *m_peaks );
//  m_sortedPeaks.swap( sortedPeaks );
//  endInsertRows();
}//void sort( int column, Wt::SortOrder order )


bool PeakModel::compare( const PeakShrdPtr &lhs, const PeakShrdPtr &rhs,
                         Columns column, Wt::SortOrder order,
                         const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  const bool asscend = (order==AscendingOrder);

  if( lhs == rhs )
    return false;

  if( !lhs || !rhs )
    return (asscend ? (lhs.get() < rhs.get()) : (lhs.get() > rhs.get()));

  switch( column )
  {
    case kMean:
      return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

    case kFwhm:
    {
      const double lw = lhs->gausPeak() ? lhs->sigma() : 0.5*lhs->roiWidth();
      const double rw = rhs->gausPeak() ? rhs->sigma() : 0.5*rhs->roiWidth();
      return (asscend ? (lw < rw) : (lw > rw));
    }//case kFwhm:
      
    case kCps:  //We dont have live time, but if peaks are for same spectrum, it doesnt matter
    case kAmplitude:
      return (asscend ? (lhs->peakArea() < rhs->peakArea()) : (lhs->peakArea() > rhs->peakArea()));

    case kUserLabel:
    {
      if( lhs->userLabel() == lhs->userLabel() )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));
      return (asscend ? (lhs->userLabel() < rhs->userLabel()) : (lhs->userLabel() > rhs->userLabel()));
    }

    case kPeakLineColor:
    {
      const std::string lhsColor = lhs->lineColor().cssText();
      const std::string rhsColor = rhs->lineColor().cssText();

      if( lhsColor == rhsColor )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));
      return (asscend ? (lhsColor < rhsColor) : (lhsColor > rhsColor));
    }

    case kHasSkew:
    {
      if( lhs->skewType() == rhs->skewType() )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));
      return (asscend ? (lhs->skewType() < rhs->skewType()) : (lhs->skewType() > rhs->skewType()));
    }
    case kSkewAmount:
      return true;

    case kType:
    {
      if( lhs->gausPeak() == rhs->gausPeak() )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

      return (asscend ? (lhs->gausPeak() < rhs->gausPeak()) : (lhs->gausPeak() > rhs->gausPeak()));
    }

    case kLowerX:
      return (asscend ? (lhs->lowerX() < rhs->lowerX()) : (lhs->lowerX() > rhs->lowerX()));

    case kUpperX:
      return (asscend ? (lhs->upperX() < rhs->upperX()) : (lhs->upperX() > rhs->upperX()));

    case kRoiCounts:
    {
      double lhs_area( 0.0 ), rhs_area( 0.0 );
      if( data )
      {
        rhs_area = data->gamma_integral( rhs->lowerX(), rhs->upperX() );
        lhs_area = data->gamma_integral( lhs->lowerX(), lhs->upperX() );
      }else
      {
        try
        {
          rhs_area = rhs->offset_integral( rhs->lowerX(), rhs->upperX(), data );
          lhs_area = lhs->offset_integral( lhs->lowerX(), lhs->upperX(), data );
        }catch(...)
        {
          //Will only fail for FlatStep, LinearStep, and BiLinearStep continuum - which I doubt we will ever get here anyway.
        }
      }
      
      return (asscend ? (lhs_area < rhs_area) : (lhs_area > rhs_area));
    }//case kRoiCounts:

    case kContinuumType:
    {
      const PeakContinuum::OffsetType lhsType = lhs->continuum()->type();
      const PeakContinuum::OffsetType rhsType = rhs->continuum()->type();

      if( lhsType == rhsType )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));
      return (asscend ? (lhsType < rhsType) : (lhsType > rhsType) );
    }

    case kIsotope:
    {
      const SandiaDecay::Nuclide *lhsParent = lhs->parentNuclide();
      const SandiaDecay::Nuclide *rhsParent = rhs->parentNuclide();

     if( lhsParent == rhsParent )
       return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

      bool thisorder = false;
      if( !lhsParent )
        thisorder = false;
      else if( !rhsParent )
        thisorder = true;
      else
        thisorder = (lhsParent->symbol < rhsParent->symbol);
      
      return (asscend ? thisorder : (!thisorder));
    }//case kIsotope

    case kPhotoPeakEnergy:
    {
      float lhsEnergy = 0.0f, rhsEnergy = 0.0f;
      try
      {
        if( lhs->hasSourceGammaAssigned() )
          lhsEnergy = lhs->gammaParticleEnergy();
      }catch(std::exception &){}
      try
      {
        if( rhs->hasSourceGammaAssigned() )
          rhsEnergy = rhs->gammaParticleEnergy();
      }catch(std::exception &){}

      if( lhs == rhs )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

      return (asscend ? (lhsEnergy < rhsEnergy) : (lhsEnergy > rhsEnergy));
    }//case kPhotoPeakEnergy:

    case kDifference:
    {
      double lhsdiff = DBL_MAX, rhsdiff = DBL_MAX;
      
      try
      {
        if( lhs->hasSourceGammaAssigned() )
          lhsdiff = lhs->gammaParticleEnergy() - lhs->mean();
      }catch(std::exception &){}
      try
      {
        if( rhs->hasSourceGammaAssigned() )
          rhsdiff = rhs->gammaParticleEnergy() - rhs->mean();
      }catch(std::exception &){}

      return (asscend ? (lhsdiff < rhsdiff) : (lhsdiff > rhsdiff));
    }//case kDifference:
      
    case kUseForShieldingSourceFit:
    {
      const bool lhsUse = lhs->useForShieldingSourceFit();
      const bool rhsUse = rhs->useForShieldingSourceFit();

      if( lhsUse == rhsUse )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

      return (asscend ? (lhsUse < rhsUse) : (rhsUse < lhsUse));
    }//case kUseForShieldingSourceFit:

    case kUseForCalibration:
    {
      const bool lhsUse = lhs->useForEnergyCalibration();
      const bool rhsUse = rhs->useForEnergyCalibration();

      if( lhsUse == rhsUse )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

      return (asscend ? (lhsUse < rhsUse) : (rhsUse < lhsUse));
    }//case kUseForCalibration:
      
    case kUseForManualRelEff:
    {
      const bool lhsUse = lhs->useForManualRelEff();
      const bool rhsUse = rhs->useForManualRelEff();

      if( lhsUse == rhsUse )
        return (asscend ? (lhs->mean() < rhs->mean()) : (lhs->mean() > rhs->mean()));

      return (asscend ? (lhsUse < rhsUse) : (rhsUse < lhsUse));
    }//case kUseForManualRelEff:

    case kCandidateIsotopes:
    case kNumColumns:
    break;
  }//switch( section )

  cerr << "PeakModel::compare(...): invalid sorting column" << endl;

  return false;
}//bool compare(...)


void PeakModel::write_peak_csv( std::ostream &outstrm,
                               std::string specfilename,
                               const PeakModel::PeakCsvType type,
                               const std::deque<std::shared_ptr<const PeakDef>> &peaks,
                               const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  bool write_html = false, write_compact = false, write_header = false;
  switch( type )
  {
    case PeakCsvType::Full:
      write_html = false;
      write_compact = false;
      write_header = true;
      break;
      
    case PeakCsvType::NoHeader:
      write_html = false;
      write_compact = false;
      write_header = false;
      break;
      
    case PeakCsvType::Compact:
      write_html = false;
      write_compact = true;
      write_header = false;
      break;
      
    case PeakCsvType::FullHtml:
      write_html = true;
      write_compact = false;
      write_header = true;
      break;
      
    case PeakCsvType::NoHeaderHtml:
      write_html = true;
      write_compact = false;
      write_header = false;
      break;
      
    case PeakCsvType::CompactHtml:
      write_html = true;
      write_compact = true;
      write_header = false;
      break;
  }//switch( type )
  
  const size_t npeaks = peaks.size();
  const string eol_char = "\r\n"; //for windows - could potentially customize this for the users operating system
  
  if( write_html )
  {
    outstrm << "<table>" << eol_char;
  }
  
  if( write_header )
  {
    if( write_html )
    {
      outstrm << "  <thead>" << eol_char
      << "    <tr>" << eol_char
      << "      <th>Centroid (keV)</th>"
      "<th>Net Area Counts</th>"
      "<th>Net Area Uncertainty</th>"
      "<th>Peak CPS</th>"
      "<th>FWHM (keV)</th>"
      "<th>FWHM (%)</th>"
      "<th>Reduced Chi2</th>"
      "<th>ROI_Total Counts</th>"
      "<th>ROI ID#</th>"
      "<th>File Name</th>";
      if( !write_compact )
      {
        outstrm << "<th>LiveTime (s)</th>"
        "<th>Date</th>"
        "<th>Time</th>"
        "<th>Nuclide</th>"
        "<th>Photopeak Energy (keV)</th>"
        "<th>ROI Lower Energy</th>"
        "<th>ROI_Upper_Energy</th>"
        "<th>Color</th>"
        "<th>User_Label</th>"
        "<th>Continuum_Type</th>"
        "<th>Skew_Type</th>"
        "<th>Continuum_Coefficients</th>"
        "<th>Skew_Coefficients</th>"
        "<th>RealTime (s)</th>"
        "<th>Peak_Type</th>";
      }
      outstrm << eol_char << "    </tr>" << eol_char
      << "  </thead>" << eol_char;
    }else
    {
      outstrm <<
      "Centroid,  Net_Area,   Net_Area,      Peak, FWHM,   FWHM,Reduced, ROI_Total,ROI, "
      "File";
      if( !write_compact )
      {
        outstrm <<
        ",         ,     ,     , Nuclide, Photopeak_Energy, ROI_Lower_Energy, ROI_Upper_Energy, Color, User_Label, Continuum_Type, "
        "Skew_Type, Continuum_Coefficients, Skew_Coefficients,         , Peak_Type";
      }
      outstrm << eol_char
      <<
      "     keV,    Counts,Uncertainty,       CPS,  keV,Percent,Chi_Sqr,    Counts,ID#, "
      "Name";
      
      if( !write_compact )
      {
        outstrm <<
        ", LiveTime, Date, Time,        ,              keV,              keV,              keV, (css),           ,               , "
        "         ,                       ,                  , RealTime,          ";
      }
      outstrm << eol_char;
    }//if( write_html ) / else
  }//if( write_header )
  
  if( write_html )
    outstrm << " <tbody>" << eol_char;
  
  for( size_t peakn = 0; peakn < npeaks; ++peakn )
  {
    const PeakDef &peak = *peaks[peakn];
    
    const double xlow = peak.lowerX();
    const double xhigh = peak.upperX();
    
    float live_time = 1.0f, real_time = 0.0f;
    double region_area = 0.0;
    SpecUtils::time_point_t meastime{};
    
    if( data )
    {
      live_time = data->live_time();
      real_time = data->real_time();
      meastime = data->start_time();
      region_area = gamma_integral( data, xlow, xhigh );
    }//if( data )
    
    string nuclide;
    char energy[32];
    energy[0] = '\0';
    
    if( peak.hasSourceGammaAssigned() )
    {
      try
      {
        const float gammaEnergy = peak.gammaParticleEnergy();
        snprintf( energy, sizeof(energy), "%.2f", gammaEnergy );
        
        if( peak.parentNuclide()
           && (peak.decayParticle() || peak.sourceGammaType()==PeakDef::AnnihilationGamma) )
        {
          nuclide = peak.parentNuclide()->symbol;
        }else if( peak.reaction() )
        {
          nuclide = peak.reaction()->name();
          SpecUtils::ireplace_all( nuclide, ",", " " );
        }else if( peak.xrayElement() )
        {
          nuclide = peak.xrayElement()->symbol + "-xray";
        }
        
        if( peak.parentNuclide() || peak.reaction() )
        {
          switch( peak.sourceGammaType() )
          {
            case PeakDef::SourceGammaType::NormalGamma:
            case PeakDef::SourceGammaType::AnnihilationGamma:
              break;
              
            case PeakDef::SourceGammaType::XrayGamma:
              nuclide += " (x-ray)";
              break;
              
            case PeakDef::SourceGammaType::SingleEscapeGamma:
              nuclide += " (S.E.)";
              break;
              
            case PeakDef::SourceGammaType::DoubleEscapeGamma:
              nuclide += " (D.E.)";
              break;
          }//switch( peak.sourceGammaType() )
        }//if( peak.parentNuclide() || peak.reaction() )
      }catch( std::exception & )
      {
      }//try / catch
    }//if( peak.hasSourceGammaAssigned() )
    
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.2f", peak.mean() );
    string meanstr = buffer;
    while( meanstr.size() < 8 )
      meanstr = " " + meanstr;
    
    snprintf( buffer, sizeof(buffer), "%.1f", peak.peakArea() );
    string areastr = buffer;
    while( areastr.size() < 10 )
      areastr = " " + areastr;
    
    snprintf( buffer, sizeof(buffer), "  %.1f", peak.peakAreaUncert() );
    string areauncertstr = buffer;
    while( areauncertstr.size() < 11 )
      areauncertstr = areauncertstr + " ";
    
    string cpststr;
    if( (live_time > 0.0f) && data )
    {
      snprintf( buffer, sizeof(buffer), "%1.4e", (peak.peakArea()/live_time) );
      cpststr = buffer;
      size_t epos = cpststr.find( "e" );
      if( epos != string::npos )
      {
        if( cpststr[epos+1]!='-' && cpststr[epos+1]!='+' )
          cpststr.insert( cpststr.begin() + epos + 1, '+' );
        while( (cpststr.size()-epos) < 5 )
          cpststr.insert( cpststr.begin() + epos + 2, '0' );
      }
    }//if( data && data->live_time() > 0.0f )
    
    
    const double width = peak.gausPeak() ? (2.35482*peak.sigma()) : 0.5*peak.roiWidth();
    snprintf( buffer, sizeof(buffer), "%.2f", width );
    string widthstr = buffer;
    while( widthstr.size() < 5 )
      widthstr = " " + widthstr;
    
    snprintf( buffer, sizeof(buffer), "%.2f", (100.0*width/peak.mean()) );
    string widthprecentstr = buffer;
    widthprecentstr += "%";
    while( widthprecentstr.size() < 7 )
      widthprecentstr = " " + widthprecentstr;
    
    snprintf( buffer, sizeof(buffer), "  %.2f", peak.chi2dof() );
    string chi2str = buffer;
    while( chi2str.size() < 7 )
      chi2str = " " + chi2str;
    
    snprintf( buffer, sizeof(buffer), "  %.1f", region_area );
    string roiareastr = buffer;
    while( roiareastr.size() < 10 )
      roiareastr = " " + roiareastr;
    
    snprintf( buffer, sizeof(buffer), "  %i", int(peakn+1) );
    string numstr = buffer;
    while( numstr.size() < 3 )
      numstr = " " + numstr;
    
    specfilename = "  ";
    SpecUtils::ireplace_all( specfilename, ",", "-" );
    
    string live_time_str;
    if( live_time > 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%.3f", live_time );
      live_time_str = buffer;
    }
    
    string real_time_str;
    if( real_time > 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%.3f", real_time );
      real_time_str = buffer;
    }
    
    string datestr, timestr;
    if( !SpecUtils::is_special(meastime) )
    {
      const string tstr = SpecUtils::to_common_string( meastime, true );
      const auto pos = tstr.find(' ');
      if( pos != string::npos )
      {
        datestr = tstr.substr(0,pos);
        timestr = tstr.substr(pos+1);
      }
    }//if( !meastime.is_special() )
    
    
    auto csvEscape = []( string &s ){
      if( s.empty() )
        return;
      
      SpecUtils::ireplace_all( s, "\"", "\"\"" );
      if( (s.find_first_of(",\n\r\"\t") != string::npos)
         || (s.front() == ' ')
         || (s.back() == ' ') )
        s = "\"" + s + "\"";
    };
    
    
    const Wt::WColor &color = peak.lineColor();
    string color_str = color.isDefault() ? string("") : color.cssText(false);
    string user_label = peak.userLabel();
    csvEscape( color_str );
    csvEscape( user_label );

    std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
    assert( continuum );
    const PeakContinuum::OffsetType cont_type = continuum->type();
    const string continuum_type = PeakContinuum::offset_type_str( cont_type );
    const string skew_type = PeakDef::to_string( peak.skewType() );
    const char *peak_type = PeakDef::to_str( peak.type() );
    
    string cont_coefs;
    switch( cont_type )
    {
      case PeakContinuum::NoOffset:
      case PeakContinuum::External:
        break;
        
      case PeakContinuum::Constant:  case PeakContinuum::Linear:
      case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
      case PeakContinuum::FlatStep:  case PeakContinuum::LinearStep:
      case PeakContinuum::BiLinearStep:
      {
        const size_t num_cont_par = PeakContinuum::num_parameters( cont_type );
        const double ref_energy = continuum->referenceEnergy();
        const vector<double> &pars = continuum->parameters();
        const vector<double> &uncerts = continuum->uncertainties();
        const vector<bool> fit_for = continuum->fitForParameter();
        
        // Print reference energy, in keV
        cont_coefs += SpecUtils::printCompact( ref_energy, 6 );
         
        for( size_t i = 0; i < num_cont_par; ++i )
        {
          cont_coefs += " " + SpecUtils::printCompact(pars[i],7);
          
          // We could print triplicate [Coefficient,Uncertainty,FitFor], but I'm unsure about doing
          //   this at the moment.
          //cont_coefs += "[";
          //if( pars.size() > i )
          //  cont_coefs += SpecUtils::printCompact( pars[i], 7 );
          //cont_coefs += ";";
          //if( uncerts.size() > i )
          //  cont_coefs += SpecUtils::printCompact( uncerts[i], 7 );
          //cont_coefs += ";";
          //if( fit_for.size() > i )
          //  cont_coefs += (fit_for[i] ? "1" : "0");
          //cont_coefs += "]";
        }
        break;
      }//case any continuum type that has parameters
    }//switch( cont_type )
    
    
    string skew_coefs;
    const size_t num_skew = PeakDef::num_skew_parameters( peak.skewType() );
    for( size_t i = 0; i < num_skew; ++i )
    {
      const auto par = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
      const double val = peak.coefficient(par);
      //const double uncert = peak.uncertainty(par);
      //skew_coefs += "[" + SpecUtils::printCompact( val, 7 )
      //              + ";" + SpecUtils::printCompact( uncert, 7 )
      //              + (peak.fitFor(par) ? ";1]" : ";1]");
      skew_coefs += (skew_coefs.empty() ? "" : " ")
                    + SpecUtils::printCompact(val,7);
    }//for( loop over skew parameters )
    
    if( write_html )
    {
      const string field_sep = "</td><td>";
      
      outstrm << "    <tr><td>"
      << meanstr
      << field_sep << areastr
      << field_sep << areauncertstr
      << field_sep << cpststr
      << field_sep << widthstr
      << field_sep << widthprecentstr
      << field_sep << chi2str
      << field_sep << roiareastr
      << field_sep << numstr
      << field_sep << specfilename;
      
      if( !write_compact )
      {
        outstrm << field_sep << live_time_str
        << field_sep << datestr
        << field_sep << timestr
        << field_sep << nuclide
        << field_sep << energy
        << field_sep << xlow
        << field_sep << xhigh
        << field_sep << color_str
        << field_sep << user_label
        << field_sep << continuum_type
        << field_sep << skew_type
        << field_sep << cont_coefs
        << field_sep << skew_coefs
        << field_sep << real_time_str
        << field_sep << peak_type;
      }
      outstrm << "</td></tr>" << eol_char;
    }else
    {
      outstrm << meanstr
      << ',' << areastr
      << ',' << areauncertstr
      << ',' << cpststr
      << ',' << widthstr
      << ',' << widthprecentstr
      << ',' << chi2str
      << ',' << roiareastr
      << ',' << numstr
      << ',' << specfilename;
      
      if( !write_compact )
      {
        outstrm << ',' << live_time_str
        << ',' << datestr
        << ',' << timestr
        << ',' << nuclide
        << ',' << energy
        << ',' << xlow
        << ',' << xhigh
        << ',' << color_str
        << ',' << user_label
        << ',' << continuum_type
        << ',' << skew_type
        << ',' << cont_coefs
        << ',' << skew_coefs
        << ',' << real_time_str
        << ',' << peak_type;
      }
      
      outstrm << eol_char;
    }//if( write_html ) / else
  }//for( loop over peaks, peakn )
  
  if( write_html )
    outstrm << " </tbody>" << eol_char << "</table>" << eol_char;
}//void PeakModel::write_peak_csv(...)


void PeakModel::write_for_and_back_peak_csv( std::ostream &outstrm,
                           std::string specfilename,
                           const PeakCsvType type,
                           const std::deque<std::shared_ptr<const PeakDef>> &peaks,
                           const std::shared_ptr<const SpecUtils::Measurement> &data,
                           std::string background_specfilename,
                           const std::deque<std::shared_ptr<const PeakDef>> *background_peaks,
                           const std::shared_ptr<const SpecUtils::Measurement> &background )
{
  write_peak_csv( outstrm, specfilename, type, peaks, data );
  if( !background_peaks || background_peaks->empty() || !background || (background->live_time() <= 0.0f) )
    return;
  
  const string eol_char = "\r\n"; //for windows - could potentially customize this for the users operating system
  
  const double scale = data->live_time() / background->live_time();
  
  switch( type )
  {
    case PeakCsvType::Full:
    case PeakCsvType::NoHeader:
    case PeakCsvType::Compact:
      outstrm << eol_char
      << "#END FOREGROUND PEAKS"
      << eol_char
      << eol_char
      << "#Background Spectrum Peaks (LiveTime " << background->live_time() << " s - scale by "
      << scale << " to make comparable):" << eol_char;
      break;
      
    case PeakCsvType::FullHtml:
    case PeakCsvType::NoHeaderHtml:
    case PeakCsvType::CompactHtml:
      outstrm << eol_char
      << "<!-- END FOREGROUND PEAKS -->" << eol_char
      << "<br />" << eol_char << "<br />" << eol_char << "<div>Background Spectrum Peaks (LiveTime "
      << background->live_time() << " s - scale by "
      << scale << " to make comparable):</div>" << eol_char;
      break;
  }//switch( type )
  
  write_peak_csv( outstrm, background_specfilename, type, *background_peaks, background );
  
  
  // Now we need to perform background subtraction
  const double nsigmaNear = 1.0;
  size_t num_affected_peaks = 0;
  deque<shared_ptr<const PeakDef>> back_sub_peaks;
  
  for( const shared_ptr<const PeakDef> &orig_peak : peaks )
  {
    assert( orig_peak );
    // TODO: we should maybe handle data-defined peaks...
    if( !orig_peak || !orig_peak->gausPeak() )
      continue;
    
    double backCounts = 0.0, backUncert2 = 0.0;
    for( const shared_ptr<const PeakDef> &backPeak : *background_peaks )
    {
      assert( backPeak );
      if( !backPeak || !backPeak->gausPeak() )
        continue;
      
      const double sigma = orig_peak->gausPeak() ? orig_peak->sigma() : 0.25*orig_peak->roiWidth();
      if( fabs(backPeak->mean() - orig_peak->mean()) < (nsigmaNear*sigma) )
      {
        backCounts += scale * backPeak->peakArea();
        const double uncert = scale * std::max( 0.0, backPeak->peakAreaUncert() );
        backUncert2 += uncert * uncert;
      }//if( fabs(backPeak.mean()-peak.mean()) < sigma )
    }//for( const PeakDef &peak : backPeaks )

    auto updated_peak = make_shared<PeakDef>( *orig_peak );
    
    if( backCounts > 0.0 )
    {
      num_affected_peaks += 1;
      const double counts = orig_peak->peakArea() - backCounts;
      const double orig_uncert = std::max( orig_peak->peakAreaUncert(), 0.0 );
      const double uncert = sqrt( orig_uncert*orig_uncert + backUncert2 );
      
      updated_peak->setPeakArea( counts );
      updated_peak->setPeakAreaUncert( uncert );
    }//if( backCounts > 0.0 )
    
    back_sub_peaks.push_back( updated_peak );
  }//for( const shared_ptr<const PeakDef> &peak : peaks )
  
  
  
  switch( type )
  {
    case PeakCsvType::Full:
    case PeakCsvType::NoHeader:
    case PeakCsvType::Compact:
      outstrm << eol_char
      << "#END BACKGROUND PEAKS"
      << eol_char
      << eol_char
      << "#Background Subtracted"
      << " (after live time normalization and whose means are within " << nsigmaNear
      << " sigma) foreground peaks (" << num_affected_peaks << " peaks adjusted):" << eol_char;
      break;
      
    case PeakCsvType::FullHtml:
    case PeakCsvType::NoHeaderHtml:
    case PeakCsvType::CompactHtml:
      outstrm << eol_char << "<!-- END BACKGROUND PEAKS -->" 
      << eol_char << "<br />" << eol_char << "<br />"
      << eol_char << "<div>Background Subtracted"
      << " (after live time normalization and whose means are within " << nsigmaNear
      << " sigma) foreground peaks (" << num_affected_peaks << " peaks adjusted):</div>" 
      << eol_char;
      break;
  }//switch( type )
  
  
  write_peak_csv( outstrm, specfilename, type, back_sub_peaks, data );
  
  
  switch( type )
  {
    case PeakCsvType::Full:
    case PeakCsvType::NoHeader:
    case PeakCsvType::Compact:
      outstrm << eol_char << "#END BACKGROUND-SUBTRACTED PEAKS" << eol_char;
      break;
      
    case PeakCsvType::FullHtml:
    case PeakCsvType::NoHeaderHtml:
    case PeakCsvType::CompactHtml:
      outstrm << eol_char << "<!-- END BACKGROUND-SUBTRACTED PEAKS -->" << eol_char;
      break;
  }//switch( type )
}//PeakModel::write_for_and_back_peak_csv(...)
