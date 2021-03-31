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

#include <map>
#include <mutex>
#include <vector>
#include <string>

#include <Wt/WServer>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/IsotopeId.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace Wt;
using namespace std;



namespace
{
  static std::mutex sm_data_dir_mutex;
  static std::string sm_data_dir = "data";
  
  string dataDirectory()
  {
    std::lock_guard<mutex> lock( sm_data_dir_mutex );
    return sm_data_dir;
  }
  
  bool less_than_by_weight( const IsotopeId::NuclideStatWeightPair &lhs,
                            const IsotopeId::NuclideStatWeightPair &rhs )
  {
    return lhs.weight < rhs.weight;
  }
  
  bool more_than_by_weight( const IsotopeId::NuclideStatWeightPair &lhs,
                           const  IsotopeId::NuclideStatWeightPair &rhs )
  {
    return lhs.weight > rhs.weight;
  }
  
  template <class T>
  bool less_than_by_first( const T &a, const T &b)
  {
    return a.first < b.first;
  }
  
}//namespace


namespace IsotopeId
{
  void setDataDirectory( const std::string &dir )
  {
    std::lock_guard<mutex> lock( sm_data_dir_mutex );
    sm_data_dir = dir;
  }

double minDetectableCounts( double energy, double sigma, std::shared_ptr<const SpecUtils::Measurement> data )
{
  //If we are calling this function, we are assuming no peak was detected,
  //  and that the entire contribution to data is background.  We will return
  //  the number of counts, from which we would be able to detect a peak 95%
  //  of the time, assuming simple coutning statistics (for now).
  //See eqns 5.52 through 5.56 in Practical Gamma-ray Spectrometry, by Gilmore.
  const double background = gamma_integral( data, energy-3.0*sigma, energy+3.0*sigma );
  return 2.33 * sqrt(background);
}//double minDetectableCounts(...)


double minDetectableCounts( std::shared_ptr<const PeakDef> peak, std::shared_ptr<const SpecUtils::Measurement> data )
{
  double contArea = 0.0;
  if( peak->continuum()->defined() )
  {
    double lowx(0.0), upperx(0.0);
    findROIEnergyLimits( lowx, upperx, *peak, data );
    contArea = peak->offset_integral( lowx, upperx );
  }else
  {
    const double lowerx = (peak->gausPeak() ? (peak->mean()-3.0*peak->sigma()) : peak->lowerX());
    const double upperx = (peak->gausPeak() ? (peak->mean()+3.0*peak->sigma()) : peak->upperX());
    
    contArea = gamma_integral( data, lowerx, upperx );
    contArea -= peak->peakArea();
  }//if( peak->continuum()->defined() ) / else
  
  contArea = max( contArea, 0.0 );

  return 2.33 * sqrt( contArea );
}//double minDetectableCounts(...)



map<const SandiaDecay::Nuclide *, int> characteristics(
      std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > all_peaks )
{
  //Now go through CharacteristicGammas.txt, and possibly augment weights
  typedef pair<double,const SandiaDecay::Nuclide*> EnergyNucPair;
  vector<EnergyNucPair > characlines;
  
  {//begin codeblock to read characteristic gammas
    string line;
    const string filename = SpecUtils::append_path( dataDirectory(), "CharacteristicGammas.txt" );
    
#ifdef _WIN32
    const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
    ifstream characteristicFile( filename.c_str(), ios_base::binary|ios_base::in );
#else
    ifstream characteristicFile( filename.c_str(), ios_base::binary|ios_base::in );
#endif
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    
    while( SpecUtils::safe_get_line( characteristicFile, line ) )
    {
      vector<string> fields;
      SpecUtils::trim( line );
      SpecUtils::split( fields, line, " \t" );
      if( fields.size() != 2 )
        continue;
      
      const double testenergy = std::stod( fields[1] );
      
      const SandiaDecay::Nuclide *testnuc = db->nuclide( fields[0] );
      
      if( testnuc )
        characlines.push_back( pair<double,const SandiaDecay::Nuclide*>(testenergy,testnuc) );
    }//while( getline( characteristicFile, line ) )
    
    std::sort( characlines.begin(), characlines.end(), &less_than_by_first< pair<double,const SandiaDecay::Nuclide*> >  );
  }//end codeblock to read characteristic gammas
  
  
  map<const SandiaDecay::Nuclide *, int> answer;
  
  for( std::shared_ptr<const PeakDef> peak : *all_peaks )
  {
    const double lowe = (peak->gausPeak() ? (peak->mean()-1.5*peak->sigma()) : peak->lowerX());
    const double highe = (peak->gausPeak() ? (peak->mean() + 1.5*peak->sigma()) : peak->upperX());
    for( const EnergyNucPair &enp : characlines )
    {
      if( (enp.first>=lowe) && (enp.first<=highe) )
      {
        if( answer.find(enp.second) == answer.end() )
          answer[enp.second] = 0;
        answer[enp.second] = answer[enp.second] + 1;
      }
    }//for( const EnergyNucPair &enp : characlines )
  }//for( std::shared_ptr<const PeakDef> peak : all_peaks )
  
  return answer;
}//characteristics(...)










double fractionDetectedWeight( const std::vector<SandiaDecay::EnergyRatePair> &source_gammas, //normailization doent matter
                           std::shared_ptr<const DetectorPeakResponse> response,
                           double shielding_an,
                           double shielding_ad,
                           std::shared_ptr<const SpecUtils::Measurement> data,
                           std::shared_ptr<const PeakDef> test_peak,
                           std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > all_peaks
                          )
{
  //loop over all expected gamma lines for nuclide, and count the fraction of
  //  counts expected to be detected, that actually where
  //  -if photopeak exists under a peak, then assume photopeak was detectable if
  //   expected counts is greater than 2.33*continuum_area.
  //  -if photopeak is not under a detected peak, assume photopeak was
  //   detectable if its counts is greater than 2.33*data_area in region.
  
  //get all gammas between +-1.5 sigma of the mean;
  const double mean  = test_peak->mean();
  const double lowE  = (test_peak->gausPeak() ? (mean-1.5*fabs(test_peak->sigma())) : test_peak->lowerX());
  const double highE = (test_peak->gausPeak() ? (mean+1.5*fabs(test_peak->sigma())) : test_peak->upperX());
  const double distance = 1.0*PhysicalUnits::m;
  const bool hasResolutionResponse = (!!response && response->hasResolutionInfo());
  
  
  vector<SandiaDecay::EnergyRatePair>::const_iterator start, end, iter;
  
  SandiaDecay::EnergyRatePair aep( 0.0, lowE );
  start = lower_bound( source_gammas.begin(), source_gammas.end(), aep,
                         &SandiaDecay::EnergyRatePair::lessThanByEnergy );
  aep.energy = highE;
  end = upper_bound( source_gammas.begin(), source_gammas.end(), aep,
                         &SandiaDecay::EnergyRatePair::lessThanByEnergy );
  
  double expectedAbund = 0.0;
//  for( iter = source_gammas.begin(); iter != source_gammas.end(); ++iter )
  for( iter = start; iter != end; ++iter )
    if( iter->energy >= lowE && iter->energy <= highE )
      expectedAbund += iter->numPerSecond;
  
//  cerr << "In fractionDetectedWeight for peak->mean()=" << test_peak->mean()
//       << " and AN=" << shielding_an << "\texpectedAbund=" << expectedAbund
//       << endl;
  
  if( expectedAbund == 0.0 )
    return 0.0;

  if( expectedAbund == 0.0 )
    throw runtime_error( "fractionDetectedWeight(...): Peak with no candiates" );
  
  const double det_sf = (!!response ? response->efficiency( mean, distance ) : 1.0);
  const double xs = MassAttenuation::massAttenuationCoeficient( shielding_an, mean );
  const double shielding_sf = exp( -shielding_ad * xs );
  const double sf = test_peak->peakArea() / shielding_sf / det_sf /expectedAbund;
  
  double total_expected = 0.0;
  double accounted_for = 0.0;
  int num_accounted_for = 0;
  int num_total = 0;
  
  for( size_t i = 0; i < source_gammas.size(); ++i )
  {
    //check to see if there is a peak cooresponding to this 'aep'
    //check to see if this 'aep' should have been detectable.
    //  if it was, increment total_expected.  If it was actually detected,
    //  increment accounted_for
    
    const double energy = source_gammas[i].energy;
    const double exp_resolution = (hasResolutionResponse ? response->peakResolutionSigma( energy ) : float((highE-lowE)/3.0) );
    const double det_eff = (!!response ? response->efficiency( energy, distance ) : 1.0);
    const double xs = MassAttenuation::massAttenuationCoeficient( shielding_an, energy );
    const double transmition = exp( -shielding_ad * xs );
    
    typedef deque< std::shared_ptr<const PeakDef> >::const_iterator Iter_t;
    Iter_t iter, nearest = all_peaks->end();
    
    //Inefficiently loop over all the 
    double min_dx = 10000000.0*PhysicalUnits::keV;
    for( iter = all_peaks->begin(); iter != all_peaks->end(); ++iter )
    {
      const std::shared_ptr<const PeakDef> &peak = *iter;
      const double mean = peak->mean();
      const double sigma = (peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth());
      const double dx = fabs(mean - energy);
      if( (distance < min_dx) && (dx < 1.5*sigma) )
      {
        min_dx = dx;
        nearest = iter;
      }
    }//for( iter = all_peaks->begin(); iter != all_peaks->end(); ++iter )
    
    double minDetectable = 0.0;
    if( nearest != all_peaks->end() )
      minDetectable = minDetectableCounts( *nearest, data );
    else
      minDetectable = minDetectableCounts( energy, exp_resolution, data );
    
    double expected = sf * det_eff * transmition * source_gammas[i].numPerSecond;
    
    if( expected > minDetectable )
    {
      ++num_total;
      total_expected += source_gammas[i].numPerSecond;
      if( nearest != all_peaks->end() )
      {
        ++num_accounted_for;
        const std::shared_ptr<const PeakDef> &nearest_peak = *nearest;
        const double peakmean = nearest_peak->mean();
        const double peaksigma = (nearest_peak->gausPeak() ? nearest_peak->sigma() : 0.25*nearest_peak->roiWidth());
        const double calib_sf = 1.0 - 0.5*fabs(energy - peakmean) / peaksigma;
        accounted_for += calib_sf * source_gammas[i].numPerSecond;
      }//if( nearest != all_peaks->end() )
    }//if( expected > minDetectable )
  }//for( size_t i = 0; i < source_gammas.size(); ++i )
  
  const double fracDetected = accounted_for / total_expected;
//  cerr << "\tsf=" << fracDetected << endl << endl;
  
  const double nphotopeak_sf = pow(1.0*num_accounted_for, 3.0) / num_total;
  
  return nphotopeak_sf * fracDetected;
}//fractionDetectedWeight(...)



vector<NuclideStatWeightPair> gammasNearInEnergy( float x, double dx )
{
  using namespace SandiaDecay;
  
  dx = fabs( dx );
  const double min_halflife = 6000.0*second;
  const double min_br = 1.0E-6;
  
  static std::shared_ptr< const EnergyToNuclideServer::EnergyNuclidePairVec > energyToNuc;
  static std::mutex energyToNucMutex;
  {//begin codeblock to ensure thread safety
    std::lock_guard<std::mutex> lock( energyToNucMutex );
    if( !energyToNuc )
    {
      EnergyToNuclideServer::setLowerLimits( min_halflife, min_br );
      energyToNuc = EnergyToNuclideServer::energyToNuclide();
    }
  }//end codeblock to ensure thread safety
  
  vector<const Nuclide *> nucs = EnergyToNuclideServer::nuclidesWithGammaInRange( x-dx, x+dx, *energyToNuc );
  
  const EnergyToNuclideServer::EnergyNuclidePairVec::const_iterator begin_nuc = energyToNuc->begin();
  const EnergyToNuclideServer::EnergyNuclidePairVec::const_iterator end_nuc   = energyToNuc->end();
  
  EnergyToNuclideServer::EnergyNuclidePairVec::const_iterator lower, upper, iter;
  const EnergyToNuclideServer::EnergyNuclidePair lowerE( x-dx, NULL ), upperE( x+dx, NULL );
  
  lower = lower_bound( begin_nuc, end_nuc, lowerE );
  upper = upper_bound( begin_nuc, end_nuc, upperE );

  vector<NuclideStatWeightPair> answer;
  
  for( iter = lower; iter != upper; ++iter )
  {
    NuclideStatWeightPair thisnuc( iter->nuclide, fabs( iter->energy - x ) );
    answer.push_back( thisnuc );
    
    //XXX - should check that the BR to this nuclide is large enough to satisfy
    //      min_br
    for( const SandiaDecay::Nuclide *parent : iter->nuclide->forebearers() )
    {
      if( parent->halfLife > min_halflife )
      {
        thisnuc.nuclide = parent;
        answer.push_back( thisnuc );
      }//if( parent->halfLife > min_halflife )
    }//for(...)
  }//for( iter = lower; iter != upper; ++iter )
  
  std::sort( answer.begin(), answer.end(), &less_than_by_weight );
  
  return answer;
}//vector<PeakToNuclideMatch> gammasNearInEnergy( float energy, double window )



void suggestNuclides(
  PeakToNuclideMatch &answer,
  std::shared_ptr<const PeakDef> peak,
  std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
  std::shared_ptr<const SpecUtils::Measurement> data,
  std::shared_ptr<const DetectorPeakResponse> response )
{
  //1) check input:
  //     1a) if 'peak' is empty, return empty results
  //     1b) if 'peaks', or 'data' is empty, return simple search
  //         by energy, sorted by proximity
  //     1c) if 'response' doesnt have resolution info, approximate this and
  //         continue on with below
  //2) get nuclides with gamma under 'peak'
  //3) find fractionDetectedWeight for each nuclide under a 'no shielding',
  //   'medium shielding', and 'heavy shielding' sceneriou, and assign best
  //   fractionDetectedWeight to be the nuclides weight
  //   -could also add in trying at different times of nuclides decay...
  //4) (optional) add in extra weight nuclides matching in file
  //   CharacteristicGammas.txt
  //5) Order results according to the determined weight, and return results

  if( !peak )
    throw runtime_error( "suggestNuclides(...): invalid input" );
  
  const double energy = peak->mean();
  const double sigma = fabs( (peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth()) );
  vector<NuclideStatWeightPair> nearNucs
                                      = gammasNearInEnergy( energy, 1.0*sigma );
  
  answer.energy = energy;
  
  if( !peaks || !data )
  {
    answer.nuclideWeightPairs = nearNucs;
    return;
  }//if( !peaks || !data )
  
  const double gcm2 = PhysicalUnits::g / PhysicalUnits::cm2;
  const double atomic_nums[]   = { 1.0, 26.0, 74.0 };
  const double areal_density[] = { 0.0*gcm2, 10.0*gcm2, 25.0*gcm2 };
  
  double max_hl = 0.0;
  
  vector<NuclideStatWeightPair> candidates;
  
  //We may get duplicate isotopes back from gammasNearInEnergy(...) if there are
  //  multiple photopeaks
  set<const SandiaDecay::Nuclide *> nucschecked;
  
  for( size_t i = 0; i < nearNucs.size(); ++i )
  {
    const SandiaDecay::Nuclide *nuc = nearNucs[i].nuclide;
    const double age = PeakDef::defaultDecayTime( nuc );
    
    if( nucschecked.count(nuc) )
      continue;
    nucschecked.insert( nuc );
    
    max_hl = max( max_hl, nuc->halfLife );
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuc, 1.0E6*SandiaDecay::becquerel );
    const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas( age,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, true );
    nearNucs[i].weight = 0.0;
    for( int j = 0; j < 3; ++j )
    {
      const double val = fractionDetectedWeight( gammas, response,
                                           atomic_nums[j], areal_density[j],
                                           data, peak, peaks );
      nearNucs[i].weight = max( nearNucs[i].weight, val );
    }//for( int j = 0; j < 3; ++j )
    
    candidates.push_back( nearNucs[i] );
  }//for( size_t i = 0; i < nearNucs.size(); ++i )
  
  map<const SandiaDecay::Nuclide *, int> charclines = characteristics( peaks );
  for( NuclideStatWeightPair &n : candidates )
  {
    if( charclines.count(n.nuclide) )
      n.weight *= 2.0*charclines[n.nuclide];
  }
  
  //Preffer isotopes with longer half lives
  //I want some thing like:
  //  n.weight *= (log(1.0 + n.nuclide->halfLife) / log(max_hl));
  //but inorder to get Th232 instead of U232 when only the 2614 peak is
  //  identified, I have to add a slight weight in....
  for( NuclideStatWeightPair &n : candidates )
    n.weight *= (1.0 - 1.2*(1.0-(log(1.0 + n.nuclide->halfLife) / log(max_hl))));
  
  
  //Slightly preffer nuclides that obtain secular equilibrium
  for( NuclideStatWeightPair &n : candidates )
  {
    if( n.nuclide->canObtainSecularEquilibrium() )
      n.weight *= 1.2;
    if( n.nuclide->canObtainPromptEquilibrium() )
      n.weight *= 1.2;
//    if( n.nuclide->decaysToStableChildren() )
//      n.weight *= 1.2;
  }//for( NuclideStatWeightPair &n : nearNucs )
  
  
  std::sort( candidates.begin(), candidates.end(), &more_than_by_weight );
  
  cerr << "For mean=" << peak->mean() << "keV, top candiates are:" << endl;
  for( size_t i = 0; i < candidates.size() && i < 10; ++i  )
    cerr << "\t" << i << "\t" << candidates[i].nuclide->symbol
         << "\tw=" << candidates[i].weight << endl;
  cerr << endl << endl;
  
  answer.nuclideWeightPairs.swap( candidates );
}//suggestNuclides(...)

//////--------------------------------------------------------------------//////




void findCandidates( vector<string> &suggestednucs,
                    std::shared_ptr<const PeakDef> peak,
                    std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > allpeaks,
                    std::shared_ptr<const DetectorPeakResponse> detector,
                    std::shared_ptr<const SpecUtils::Measurement> data )
{
  if( !peak || !allpeaks || !data )
  {
    suggestednucs.clear();
    return;
  }
  
  
  try
  {
    if( !!detector && (!detector->isValid() || !detector->hasResolutionInfo()) )
      detector.reset();
    
    if( !detector || !detector->isValid() )
    {
      std::shared_ptr<DetectorPeakResponse> detPtr
                                  = std::make_shared<DetectorPeakResponse>();
      detector = detPtr;
      
      const string basename = SpecUtils::append_path( dataDirectory(), "GenericGadrasDetectors" );
      
      string csvfilename = SpecUtils::append_path( basename, "HPGe 40%/Efficiency.csv" );
      string datFilename = SpecUtils::append_path( basename, "HPGe 40%/Detector.dat" );
      
      if( data->num_gamma_channels() < HIGH_RES_NUM_CHANNELS )
      {
        csvfilename = SpecUtils::append_path( basename, "NaI 3x3/Efficiency.csv" );
        datFilename = SpecUtils::append_path( basename, "NaI 3x3/Detector.dat" );
      }//if( this is a NaI or other low resolution detector )
      
#ifdef _WIN32
      const std::wstring wcsvfilename = SpecUtils::convert_from_utf8_to_utf16(csvfilename);
      const std::wstring wdatFilename = SpecUtils::convert_from_utf8_to_utf16(datFilename);
      ifstream csv( wcsvfilename.c_str(), ios_base::binary|ios_base::in );
      ifstream datFile( wdatFilename.c_str(), ios_base::binary|ios_base::in );
#else
      ifstream csv( csvfilename.c_str(), ios_base::binary|ios_base::in );
      ifstream datFile( datFilename.c_str(), ios_base::binary|ios_base::in );
#endif
      if( csv.good() && datFile.good() )
      {
        detPtr->fromGadrasDefinition( csv, datFile );
      }else
      {
        cerr << "guessIsotopesForPeaks(...): error opening default detector file" << endl;
        //ToDo: get approximate intrinsic efficiency formula for a HPGe and NaI detector here, and just use those for this function.
        detPtr->setIntrinsicEfficiencyFormula( "1.0", 3.0*PhysicalUnits::cm, PhysicalUnits::keV, 0.0f, 0.0f );
      }
      //      try{ detPtr->fitResolution( allpeaks, DetectorPeakResponse::kGadrasResolutionFcn ); }catch(...){}
    }//if( !detector || !detector->isValid() )
    
    
    PeakToNuclideMatch suggestedNucs;
    suggestNuclides( suggestedNucs, peak, allpeaks, data, detector );
    
    const vector<NuclideStatWeightPair> &sugestions
                                             = suggestedNucs.nuclideWeightPairs;
    char buffer[64];
    
    for( size_t j = 0; j < sugestions.size() && j < 10; ++j )
    {
      const NuclideStatWeightPair &p = sugestions[j];
      
      PeakDef::SourceGammaType sourceGammaType;
      size_t radparticleIndex;
      const SandiaDecay::Transition *transition = NULL;
      const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
      PeakDef::findNearestPhotopeak( p.nuclide, suggestedNucs.energy,
                                     4.0*sigma, false, transition,
                                     radparticleIndex, sourceGammaType );
      
      if( p.nuclide && (transition || (sourceGammaType==PeakDef::AnnihilationGamma)) )
      {
        switch( sourceGammaType )
        {
          case PeakDef::NormalGamma:
            snprintf( buffer, sizeof(buffer), "%s %.2f keV",
                     p.nuclide->symbol.c_str(),
                     transition->products[radparticleIndex].energy );
            break;
          case PeakDef::AnnihilationGamma:
            snprintf( buffer, sizeof(buffer), "%s 511 keV",
                     p.nuclide->symbol.c_str() );
            break;
          case PeakDef::SingleEscapeGamma:
            snprintf( buffer, sizeof(buffer), "S.E. %s %.2f keV",
                     p.nuclide->symbol.c_str(),
                     transition->products[radparticleIndex].energy );
            break;
          case PeakDef::DoubleEscapeGamma:
            snprintf( buffer, sizeof(buffer), "D.E. %s %.2f keV",
                     p.nuclide->symbol.c_str(),
                     transition->products[radparticleIndex].energy );
            break;
          case PeakDef::XrayGamma:
            snprintf( buffer, sizeof(buffer), "%s xray %.2f keV",
                     p.nuclide->symbol.c_str(),
                     transition->products[radparticleIndex].energy );
            break;
        }//switch( sourceGammaType )
        
        suggestednucs.push_back( buffer );
      }//if( transition )
    }//sugestions.size()
  }catch( std::exception &e )
  {
    cerr << "findCandidates(for right click menu) caught exception: "
    << e.what() << endl;
  }//try / catch
}//void findCandidates( )


void findCharacteristics( vector<string> &characteristicnucs,
                         std::shared_ptr<const PeakDef> peak )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  if( !peak || !peak->gausPeak() || !db )
    return;
  
  const ReactionGamma *reactionDb = NULL;
  vector<const ReactionGamma::Reaction *> suggest_reactions;
  try
  {
    reactionDb = ReactionGammaServer::database();
  }catch(...){ cerr << "Failed to open gamma reactions XML file" << endl; }
  
  
  //We could probably cache the parsed file into memorry, to save CPU, but
  //  whatever for now.
  const string filename = SpecUtils::append_path( dataDirectory(), "PhotoPeak.lis" );

#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream input( wfilename.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream input( filename.c_str(), ios_base::binary | ios_base::in );
#endif
  
  if( !input.good() )
  {
    cerr << "Failed to open " << filename << endl;
    return;
  }//if( !input.good() )
  
  //PhotoPeak.lis actually looks to be structured such that we could speed up
  //  access to the line with the starting energy we care about, by fitting a
  //  polynomial that gives the offset, or some offset below, the start of
  //  the line for that energy.
  //This would of course also require implementing a unit test to check this,
  //  and possibly also a mechanism to account for if the user makes a change
  //  such as deleting a line in the file.
  
  
  const double mean = peak->mean();
  const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
  const double nsigma = (0.03 < (sigma/mean)) ? 1.25 : 2.5;
  double width = nsigma*sigma;
  if( mean > 3200.0 )
    width = std::max( width, 20.0 );
  if( mean > 4500.0 )
    width = std::max( width, 30.0 );
  
  const double minenergy = peak->mean() - width;
  const double maxenergy = peak->mean() + width;
  
  //First place the results into a map, keyed off of difference in energy from
  //  the peak mean, so this way results will be returned, ordered by how close
  //  the gamma is to the peak mean.
  map<double,vector<string> > results;
  
  string line;
  char buff[64];
  
  while( SpecUtils::safe_get_line( input, line ) )
  {
    //Example lines:
    //    " 2836.35  Al n-gamma          1.851	1.548"
    //    " 2002.15  125Sn               1.920	26.299"
    //    " 1082.60  62Co de             0.164	0.225"
    //    "  883.24  238Pu             7.7E-07	2.7E-02"
    //    "  706.40  57Co                0.025	0.077"
    //    "  152.00  Ta n-gamma         14.964	7.482"
    //    "  184.40  166Hom             73.200	166.667"
    //    "  185.80  245Cm               0.010	0.100"
    if( line.size() < 25 )
      continue;
    
    try
    {
      string energystr = line.substr( 0, 9 );
      SpecUtils::trim( energystr );
      double energy = std::stod( energystr.c_str() );
      
      if( energy < minenergy )
        continue;
      
      if( energy > maxenergy )
        break;
      
      const double diff = fabs( mean - energy );
      
      string nucstr = line.substr( 10, 16 );
      SpecUtils::trim( nucstr );
      
      if( nucstr.size() < 2 )
        continue;
      
      const bool isReaction = (nucstr.find("n-gamma") != string::npos);
      const bool isXRay = (nucstr.find("x-ray") != string::npos);
      const bool isDoubleEscape = (nucstr.find(" de") != string::npos);
      const bool isSingleEscape = (nucstr.find(" se") != string::npos);
      
      if( isXRay )
      {
        results[diff].push_back( nucstr );
      }else if( !isDoubleEscape && !isSingleEscape && !isReaction )
      {
        //166Hom
        string::size_type pos = nucstr.find_first_not_of( "0123456789" );
        if( pos == string::npos )
          continue;
        
        string anstr = nucstr.substr( 0, pos );
        nucstr = nucstr.substr( pos );
        if( nucstr.size() >= 3 && nucstr[2]=='m' )
        {
          nucstr = nucstr.substr( 0, 2 );
          anstr += 'm';
        }
        nucstr += anstr;
        
        
        const SandiaDecay::Nuclide *nuc = db->nuclide( nucstr );
        if( !nuc )
        {
          cerr << "Couldnt find '" << nucstr << "' in SandiaDecay db" << endl;
          continue;
        }
        
        //Energies in PhotoPeak.lis may not exactly line up with the energies
        //  in SandiaDecay, and since we are relying on text matching later on
        //  to keep from giving duplicates to users, we have to match this gamma
        //  energy from PhotoPeak.lis up to the one in SandiaDecay
        size_t index = 0;
        const SandiaDecay::Transition *trans = NULL;
        PeakDef::SourceGammaType type;
        PeakDef::findNearestPhotopeak( nuc, energy, -1.0, false, trans, index, type );
        if( trans )
        {
          energy = trans->products[index].energy;
        }else
        {
          cerr << "Couldnt find gamma at " << energy <<" keV for " << nucstr
               << " in SandiaDecay database" << endl;
          continue;
        }
        
        if( type == PeakDef::NormalGamma )
        {
          snprintf( buff, sizeof(buff), "%s %.2f keV", nucstr.c_str(), energy );
          results[diff].push_back( buff );
        }
      }else if( isReaction && reactionDb && !isSingleEscape && !isDoubleEscape )
      {
//        "  152.00  Ta n-gamma         14.964	7.482"
        const string::size_type pos = nucstr.find( " " );
        if( pos == string::npos )
          continue;
        
        const SandiaDecay::Element *el = db->element( nucstr.substr(0,pos) );
        if( !el )
        {
          cerr << "Couldnt get element for element '" << nucstr.substr(0,pos)
               << "'" << endl;
          continue;
        }//if( !el )
        
        vector<const ReactionGamma::Reaction *> reactions;
        reactionDb->reactions( el, reactions );
        
        for( const ReactionGamma::Reaction *rctn : reactions )
        {
          const vector<ReactionGamma::EnergyAbundance> &gammas = rctn->gammas;
          for( const ReactionGamma::EnergyAbundance &g : gammas )
          {
            if( fabs(g.energy - energy) < 2.0 )
            {
              snprintf( buff, sizeof(buff), "%s %.2f keV",
                        rctn->name().c_str(), g.energy );
              
              //due to the nested loops, and the non-uniquness of only matching
              //  down to 2 keV, the same reaction/gamma may be inserted into
              //  the result more than once; this is protected against the user
              //  seing it when all results are collated together.
              results[diff].push_back( buff );
            }//if( close enough energy )
          }//for( const ReactionGamma::EnergyAbundance &g : gammas )
        }//for( const ReactionGamma::Reaction *rctn : reactions )
      }else if( !isReaction && (isSingleEscape || isDoubleEscape) )
      {
        //Not implemented
      }//if( isXRay ) / if else / else
      
    }catch( std::exception & )
    {
      
    }//try / catch
  }//while( SpecUtils::safe_get_line( input, line ) )
  
  
  for( map<double,vector<string> >::const_iterator i = results.begin();
      i != results.end(); ++i )
  {
    for( const string &nuc : i->second )
      characteristicnucs.push_back( nuc );
  }
  
}//void findCharacteristics(...)


void isotopesFromOtherPeaks( vector<string> &otherpeaknucs,
                            std::shared_ptr<const PeakDef> peak,
                            std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > allpeaks )
{
  if( !peak || !allpeaks )
    return;
  
  const ReactionGamma *reactionDb = NULL;
  vector<const ReactionGamma::Reaction *> suggest_reactions;
  try
  {
    reactionDb = ReactionGammaServer::database();
  }catch(...){ cerr << "Failed to open gamma reactions XML file" << endl; }

  
  const double mean = peak->mean();
  const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
  const double nsigma = (0.03 < (sigma/mean)) ? 1.25 : 2.5;
  double width = nsigma*sigma;
  if( mean > 3200.0 )
    width = std::max( width, 20.0 );
  if( mean > 4500.0 )
    width = std::max( width, 30.0 );
  
  const double minenergy = peak->mean() - width;
  const double maxenergy = peak->mean() + width;

  char buff[64];
  typedef std::shared_ptr<const PeakDef> ConstPeakDefPtr;
  
  for( const ConstPeakDefPtr &p : *allpeaks )
  {
    
    if( p == peak )
      continue;
    
    if( p->parentNuclide() )
    {
      const SandiaDecay::Nuclide *nuc = p->parentNuclide();
      size_t index = 0;
      const SandiaDecay::Transition *trans = NULL;
      PeakDef::SourceGammaType type;
      PeakDef::findNearestPhotopeak( nuc, mean, sigma, false, trans, index, type );
     
      if( trans )
      {
        const float candidateE = trans->products[index].energy;
        if( candidateE > minenergy && candidateE < maxenergy )
        {
          switch( type )
          {
            case PeakDef::NormalGamma:
              snprintf( buff, sizeof(buff), "%s %.2f keV",
                       nuc->symbol.c_str(), candidateE );
              otherpeaknucs.push_back( buff );
            break;
              
            case PeakDef::AnnihilationGamma:
            break;
              
            case PeakDef::SingleEscapeGamma:
            break;
              
            case PeakDef::DoubleEscapeGamma:
            break;
              
            case PeakDef::XrayGamma:
              snprintf( buff, sizeof(buff), "%s xray %.2f keV",
                       nuc->symbol.c_str(), candidateE );
              otherpeaknucs.push_back( buff );
            break;
          }//switch( type )
        }//if( candidateE > minenergy && candidateE < maxenergy )
      }//if( trans )
      
      if( mean > 510.99
          && (p->sourceGammaType() == PeakDef::NormalGamma
               || p->sourceGammaType() == PeakDef::SingleEscapeGamma) )
      {
        PeakDef::findNearestPhotopeak( nuc, mean+510.99, sigma, false, trans, index, type );
        if( trans )
        {
          const float candidateE = trans->products[index].energy;
          
          if( candidateE > (minenergy+510.99) && candidateE < (maxenergy+510.99) )
          {
            //Should also require the BR of this candidate to be detectable.
/*
            const float se_frac = (p->sourceGammaType() == PeakDef::NormalGamma) ? 0.2 : 1.0;
            const float br = se_frac * trans->branchRatio * trans->products[index].intensity;

            const SandiaDecay::Transition *parent_trans = p->nuclearTransition();
            const SandiaDecay::RadParticle *parent_particle = p->decayParticle();
            const float parent_energy = p->gammaParticleEnergy();
            const float parent_br = parent_trans->branchRatio * parent_particle->intensity;
*/
            
            snprintf( buff, sizeof(buff), "%s S.E. %.2f keV",
                      nuc->symbol.c_str(), candidateE );
            otherpeaknucs.push_back( buff );
          }//if( candidateE > minenergy && candidateE < maxenergy )
        }//if( trans )
      }//if( mean > 510.99 )
      
      if( p->sourceGammaType() == PeakDef::NormalGamma )
      {
        PeakDef::findNearestPhotopeak( nuc, mean+2.0*510.99, sigma, false, trans, index, type );
        if( trans )
        {
          const float candidateE = trans->products[index].energy;
          if( candidateE > (minenergy+2.0*510.99) && candidateE < (maxenergy+2.0*510.99) )
          {
             snprintf( buff, sizeof(buff), "%s D.E. %.2f keV",
                    nuc->symbol.c_str(), candidateE );
            otherpeaknucs.push_back( buff );
          }//if( candidateE > minenergy && candidateE < maxenergy )
        }//if( trans )
      }
    }else if( p->xrayElement() )
    {
      const SandiaDecay::Element *el = p->xrayElement();

      for( const SandiaDecay::EnergyIntensityPair &xray : el->xrays )
      {
        if( (xray.energy > minenergy) && (xray.energy < maxenergy) )
        {
          otherpeaknucs.push_back( el->symbol + " x-ray" );
          break;
        }//
      }//for( const SandiaDecay::EnergyIntensityPair &xray : el->xrays )
    }else if( p->reaction() && reactionDb )
    {
      for( const ReactionGamma::EnergyAbundance &g : p->reaction()->gammas )
      {
        if( g.energy > minenergy && g.energy < maxenergy )
        {
          snprintf( buff, sizeof(buff), "%s %.2f keV",
                    p->reaction()->name().c_str(), g.energy );
        }
      }//for( const ReactionGamma::EnergyAbundance &g : gammas )
    }//if( p->parentNuclide() ) / else / else
  }//for( const ConstPeakDefPtr &p : *allpeaks )
  
}//void findCharacteristics(...)


void peakCandidateSourceFromRefLines( std::shared_ptr<const PeakDef> peak, const std::vector<ReferenceLineInfo> &showingRefLines,
                                       std::shared_ptr< vector<string> > candidates )
{
  if( !peak || !candidates )
    return;
  
  char buffer[128];
  
  const auto initialsize = candidates->size();
  const double mean = peak->mean();
  const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
  
  for( const ReferenceLineInfo &ref : showingRefLines )
  {
    //Right now just look for nuclides
    if( !ref.nuclide )
      continue;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( ref.nuclide, 1.0E6*SandiaDecay::becquerel );
    const double age = ((ref.age >= 0.0) ? ref.age : PeakDef::defaultDecayTime(ref.nuclide) );
    
    const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas( age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
    const vector<SandiaDecay::EnergyRatePair> xrays = (ref.showXrays ? mix.xrays( age ) : vector<SandiaDecay::EnergyRatePair>{});
    
    vector<pair<string,double>> trials;
    
    //We will take the simeple approach and assume transmision and detection efficiency
    //  are about the same for all lines we care about, and create a simple
    //  weight of sf/(0.25*sigma + distance)
    for( const auto &gamma : gammas )
    {
      const double distance = fabs( gamma.energy - mean );
      
      if( distance > 4.0*sigma )
        continue;
      
      const double amp = gamma.numPerSecond;
      const double weight = amp / (0.25*sigma + distance);
      
      snprintf( buffer, sizeof(buffer), "%s %.2f keV",
                ref.nuclide->symbol.c_str(), gamma.energy );
      
      trials.emplace_back( buffer, weight );
    }//for( const auto &gammas : gammas )
    
    std::sort( begin(trials), end(trials),
               []( const pair<string,double> &lhs, const pair<string,double> &rhs ){
                 return lhs.second > rhs.second;
    } );
    
    //Add a maximum of three trial RefLines
    for( size_t i = 0; i < 3 && i < trials.size(); ++i )
      candidates->push_back( trials[i].first );
  }//for( const ReferenceLineInfo &ref : showingRefLines )
  
  
  if( initialsize != candidates->size() )
    candidates->push_back( "" );  //create the spacer
}//void peakCandidateSourceFromRefLines(...)
  

void populateCandidateNuclides( std::shared_ptr<const SpecUtils::Measurement> data,
                               std::shared_ptr<const PeakDef> peak,
                               std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > hintpeaks,
                               std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > userpeaks,
                               const std::vector<ReferenceLineInfo> showingRefLines,
                               std::shared_ptr<const DetectorPeakResponse> detector,
                               const std::string sessionid,
                               std::shared_ptr< vector<string> > candidates,
                               boost::function<void(void)> doupdate )
{
  //TODO: should probably add in some error logging or something
  if( !data || !peak || !candidates || doupdate.empty() )
    return;
  
  char buffer[128];
  
  //we will use 'entries' to make sure we dont have duplicates
  set<string> entries;
  vector<string> suggestednucs, characteristicnucs, otherpeaksnucs;
  
  SpecUtilsAsync::ThreadPool pool;
  
  /*
  std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > allpeaks
                                                                    = hintpeaks;
  
  if( !allpeaks || (!!userpeaks && userpeaks->size() >= allpeaks->size()) )
    allpeaks = userpeaks;
  */
  
  //Construct a set of peaks that has all user peaks, but also the hint peaks
  //  the user hasnt identified yet.
  std::shared_ptr<deque< std::shared_ptr<const PeakDef> > > allpeaks;
  if( userpeaks )
    allpeaks = make_shared< deque<std::shared_ptr<const PeakDef>> >( *userpeaks );
  else
    allpeaks = make_shared< deque<std::shared_ptr<const PeakDef>> >();
  
  auto &peaks = *allpeaks;
  std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMeanShrdPtr );
  
  if( hintpeaks )
  {
    for( const auto &hintp : *hintpeaks )
    {
      if( !hintp || !hintp->gausPeak() )  //shouldnt ever happen
        continue;
      const float mean = hintp->mean();
      const float sigma = hintp->sigma();
      
      auto pos = std::lower_bound( begin(peaks), end(peaks), hintp, &PeakDef::lessThanByMeanShrdPtr );
      bool already_have = false;
      if( pos != begin(peaks) && fabs((*(pos-1))->mean()-mean)<sigma )  //I dont think we need this one
        already_have = true;
      if( pos != end(peaks) && fabs((*pos)->mean()-mean)<sigma )
        already_have = true;
      if( !already_have )
        peaks.insert( pos, hintp );
    }//for( const auto &hintp : *hintpeaks )
  }//if( hintpeaks )
  
  
  pool.post( boost::bind( &findCandidates,
                         boost::ref(suggestednucs),
                         peak, allpeaks, detector, data ) );
  
  pool.post( boost::bind( &findCharacteristics,
                         boost::ref(characteristicnucs), peak ) );
  
  if( userpeaks )
    pool.post( boost::bind( &isotopesFromOtherPeaks,
                           boost::ref(otherpeaksnucs), peak, userpeaks ) );
  
  const char *modifier = "";
  switch( peak->sourceGammaType() )
  {
    case PeakDef::NormalGamma:       break;
    case PeakDef::AnnihilationGamma: break;
    case PeakDef::SingleEscapeGamma: modifier = "S.E. "; break;
    case PeakDef::DoubleEscapeGamma: modifier = "D.E. "; break;
    case PeakDef::XrayGamma:         modifier = "xray "; break;
  }//switch( peak->sourceGammaType() )
  
  if( peak->parentNuclide() && peak->decayParticle() )
  {
    snprintf( buffer, sizeof(buffer), "%s %s%.2f keV",
             peak->parentNuclide()->symbol.c_str(), modifier,
             peak->decayParticle()->energy );
    
    candidates->push_back( buffer );
    entries.insert( buffer );
  }//if( peak->parentNuclide() && peak->decayParticle() )
  
  if( peak->xrayElement() )
  {
    snprintf( buffer, sizeof(buffer), "%s %.2f keV",
             peak->xrayElement()->symbol.c_str(), peak->xrayEnergy() );
    candidates->push_back( buffer );
    entries.insert( buffer );
  }//if( peak->xrayElement() )
  
  if( peak->reaction() )
  {
    snprintf( buffer, sizeof(buffer), "%s %s%.2f keV",
              peak->reaction()->name().c_str(), modifier,
              peak->reactionEnergy() );
    candidates->push_back( buffer );
    entries.insert( buffer );
  }//if( peak->reaction() )
  
  pool.join();
  
  
  //We want to filter out candidate lines for a nuclide, where there is clearly
  //  a better candidate line for that nuclide.  This shows up mostly for NaI
  //  spectra where a peak will cover a large B.R. as well as many small B.R.
  //  lines
  //
//  const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
//  const double delta_e = fabs( product.energy - energy);
//  double scaleDeltaE = (0.1*windowHalfWidth + delta_e) / fracIntensity;
//  if( windowHalfWidth <= 0.0 )
//    scaleDeltaE = delta_e;
  //characteristicnucs have the most entries.  Should sort result be energy, or closeness or something
  
  //for( const std::string &nuc : suggestednucs )
  //  cerr << "suggestednucs: " << nuc << endl;
  //for( const std::string &nuc : characteristicnucs )
  //  cerr << "characteristicnucs: " << nuc << endl;
  //for( const std::string &nuc : otherpeaksnucs )
  //  cerr << "otherpeaksnucs: " << nuc << endl;
  
  
  //Go through and add found sources from each source to 'candidates', but
  //  first add in sources found in the highest number of places (which has max
  //  of 3 right now).
  map<string,int> nsources;
  
  for( const std::string &nuc : suggestednucs )
    nsources[nuc] = nsources.count(nuc) ? nsources[nuc] : 1;
  for( const std::string &nuc : characteristicnucs )
    nsources[nuc] = nsources.count(nuc) ? nsources[nuc] : 1;
  for( const std::string &nuc : otherpeaksnucs )
    nsources[nuc] = nsources.count(nuc) ? nsources[nuc] : 1;
  
  //Aritifically bump up the rating for background and common sources; note that
  //  this is put in just to play around with, we may not want to keep it, I
  //  really dont like hard coding any nuclides anywhere in InterSpec (20141120)
  for( map<string,int>::iterator i = nsources.begin(); i != nsources.end(); ++i )
  {
    //Not a complete list...
    if( SpecUtils::starts_with(i->first, "K40" )
        || SpecUtils::starts_with(i->first, "U235" )
        || SpecUtils::starts_with(i->first, "Th232" )
        || SpecUtils::starts_with(i->first, "Ra226" )
        || SpecUtils::starts_with(i->first, "U238" )
        || SpecUtils::starts_with(i->first, "F18" )
        || SpecUtils::starts_with(i->first, "Cs137" )
        || SpecUtils::starts_with(i->first, "Co60" )
        || SpecUtils::starts_with(i->first, "Na22" )
        || SpecUtils::starts_with(i->first, "Ba133 " )
        || SpecUtils::starts_with(i->first, "Na22" )
        || SpecUtils::starts_with(i->first, "Pu239" )
        || SpecUtils::starts_with(i->first, "Np239" )
        || SpecUtils::starts_with(i->first, "Pu241" ) )
      i->second += 1;
  }
  
  int occurance = 0;
  for( map<string,int>::const_iterator i = nsources.begin();
       i != nsources.end(); ++i )
  {
    occurance = std::max( occurance, i->second );
  }
  
  
  //Add a blank entry so we'll get a spacer in the right click popup menu
  if( candidates->size() && nsources.size() )
    candidates->push_back( "" );
  
  const double mean = peak->mean();
  const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
  
  vector<string> escapes;
  for( deque< std::shared_ptr<const PeakDef> >::const_iterator i = allpeaks->begin();
      i != allpeaks->end(); ++i )
  {
    std::shared_ptr<const PeakDef> p = *i;
    const double pmean = p->mean();
    const double tolerance = 1.5*(p->gausPeak() ? p->sigma() : 0.125*p->roiWidth());
    
    const bool singleEscape = (fabs((pmean-mean) - 510.99) < tolerance);
    const bool doubleEscape = (fabs((pmean-mean) - 2.*510.99) < tolerance);
    const bool normal = (p->sourceGammaType() == PeakDef::NormalGamma);
    
    if( normal && (singleEscape || doubleEscape) )
    {
      char buffer[128];
      const char *prefix = singleEscape ? "S.E." : "D.E.";
      
      if( p->reaction() )
      {
        snprintf( buffer, sizeof(buffer), "%s %s %.2f keV",
                  p->reaction()->name().c_str(), prefix,
                  p->gammaParticleEnergy() );
        
      }else if( p->parentNuclide() )
      {
        snprintf( buffer, sizeof(buffer), "%s %s %.2f keV",
                  p->parentNuclide()->symbol.c_str(), prefix,
                  p->gammaParticleEnergy() );
      }else
      {
        //The opening and closing paranthesis as first/last character will cause
        //  InterSpec::updateRightClickNuclidesMenu(...) to disable the
        //  menu item in the right click menu
        snprintf( buffer, sizeof(buffer), "(%s of %.2f keV Peak)",
                  prefix, pmean );
      }
      
      
      if( !std::count( candidates->begin(), candidates->end(), buffer ) )
      {
        escapes.push_back( buffer );
        entries.insert( buffer );
      }
    }//f( normal && (singleEscape || doubleEscape) )
  }//for( loop over allpeaks looking for potential escape peaks )
  
  if( escapes.size() )
  {
    candidates->insert( candidates->end(), escapes.begin(), escapes.end() );
    candidates->push_back( "" );
  }//if( escapes.size() )
  
  
  //Now look for possible sum peak combinations, and let the user know of the
  //  possiblity (but selecting them wont do any good)
  //TODO: Extend this to check for peaks of identified nuclides who are
  //  higher in energy than the detector when up to
  vector<string> summs;
  for( size_t i = 0; i < allpeaks->size(); ++i )
  {
    const double peakAmp = peak->gausPeak() ? peak->amplitude() : peak->areaFromData(data);
    if( !peak->gausPeak() )
      continue;  //be lazy for the moment and not deal with this edge case
      
    char buffer[128];
    std::shared_ptr<const PeakDef> rpeak = (*allpeaks)[i];
    
    const double rpeakAmp = rpeak->gausPeak() ? rpeak->amplitude() : rpeak->areaFromData(data);
    
    if( fabs(2.0*rpeak->mean() - mean) < 2.0*sigma && (peakAmp < rpeakAmp) )
    {
      cerr << "Have sum for peak " << rpeak.get() << endl;
      snprintf( buffer, sizeof(buffer), "(Sum from %.2f keV peak)", rpeak->mean() );
      summs.push_back( buffer );
    }
    
    for( size_t j = 0; j < i; ++j )
    {
      std::shared_ptr<const PeakDef> lpeak = (*allpeaks)[j];
      
      const double lpeakAmp = lpeak->gausPeak() ? lpeak->amplitude() : lpeak->areaFromData(data);
      
      if( !lpeak->gausPeak() )
        continue; //be lazy for the moment and not deal with this edge case
      
      if( fabs(lpeak->mean() + rpeak->mean() - mean) < 2.0*sigma
          && ((rpeakAmp+lpeakAmp) > peakAmp) )
      {
        cerr << "Have double sum for peak " << rpeak.get() << " and " << lpeak.get() << endl;
        snprintf( buffer, sizeof(buffer), "(Sum %.2f and %.2f peaks)",
                  lpeak->mean(), rpeak->mean() );
        summs.push_back( buffer );
        entries.insert( buffer );
      }
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 0; i < allpeaks->size(); ++i )
  
  if( summs.size() )
  {
    candidates->insert( candidates->end(), summs.begin(), summs.end() );
    candidates->push_back( "" );
  }//if( escapes.size() )
  
  //Look at the currently showing reference lines, and add those in
  peakCandidateSourceFromRefLines( peak, showingRefLines, candidates );
  
  //These nested loops make me cringe at the inefficiency
  for( ; occurance > 0; --occurance )
  {
    for( const std::string &nuc : suggestednucs )
    {
      if( nsources[nuc] != occurance )
        continue;
      if( !entries.count(nuc) )
        candidates->push_back( nuc );
      entries.insert( nuc );
    }//for( const std::string &nuc : suggestednucs )
    
    for( const std::string &nuc : otherpeaksnucs )
    {
      if( nsources[nuc] != occurance )
        continue;
      if( !entries.count(nuc) )
        candidates->push_back( nuc );
      entries.insert( nuc );
    }//for( const std::string &nuc : otherpeaksnucs )
    
    for( const std::string &nuc : characteristicnucs )
    {
      if( nsources[nuc] != occurance )
        continue;
      if( !entries.count(nuc) )
        candidates->push_back( nuc );
      entries.insert( nuc );
    }//for( const std::string &nuc : characteristicnucs )
  }//for( ; occurance > 0; --occurance )
  
  
//  //Add in suggestednucs
//  if( candidates->size() && suggestednucs.size() )
//    candidates->push_back( "" );
//  for( const std::string &nuc : suggestednucs )
//  {
//    if( !entries.count(nuc) )
//      candidates->push_back( nuc );
//    entries.insert( nuc );
//  }//for( const std::string &nuc : suggestednucs )
//  
//  //Add in characteristicnucs
//  if( candidates->size() && characteristicnucs.size() )
//    candidates->push_back( "" );
//  for( const std::string &nuc : characteristicnucs )
//  {
//    if( !entries.count(nuc) )
//      candidates->push_back( nuc );
//    entries.insert( nuc );
//  }//for( const std::string &nuc : characteristicnucs )
//  
//  //Add in otherpeaksnucs
//  if( candidates->size() && otherpeaksnucs.size() )
//    candidates->push_back( "" );
//  for( const std::string &nuc : otherpeaksnucs )
//  {
//    if( !entries.count(nuc) )
//      candidates->push_back( nuc );
//    entries.insert( nuc );
//  }//for( const std::string &nuc : otherpeaksnucs )
  
  Wt::WServer *server = Wt::WServer::instance();
  if( server )
    server->post( sessionid, doupdate );
  else
    cerr << "populateCandidateNuclides(): error getting pointer to server"
         << endl;
}//void populateCandidateNuclides(...)

}//namespace IsotopeId
