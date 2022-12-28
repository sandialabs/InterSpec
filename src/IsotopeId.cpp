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
#include "InterSpec/PeakFitUtils.h"
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
  

// We'll read CharacteristicGammas.txt in just once; if you want to access this files contents,
//  just call the #characteristicGammas() function, and it will take care of everything for you.
//  There are about 110 entries in "CharacteristicGammas.txt"
std::mutex sm_CharacteristicGammas_mutex;
bool sm_CharacteristicGammas_inited = false;
vector<tuple<string, const SandiaDecay::Nuclide *, float>> sm_CharacteristicGammas;

const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &characteristicGammas()
{
  std::lock_guard<std::mutex> lock( sm_CharacteristicGammas_mutex );
  if( sm_CharacteristicGammas_inited )
    return sm_CharacteristicGammas;
  
  string line;
  const string filename = SpecUtils::append_path( dataDirectory(), "CharacteristicGammas.txt" );
  
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream characteristicFile( filename.c_str(), ios_base::binary|ios_base::in );
#else
  ifstream characteristicFile( filename.c_str(), ios_base::binary|ios_base::in );
#endif
  
  sm_CharacteristicGammas_inited = true;
  
  if( !characteristicFile.is_open() )
  {
    cerr << "Failed to open '" << filename << "' - nuclide suggestions will not use it." << endl;
    return sm_CharacteristicGammas;
  }
  
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    assert( db );
    cerr << "Failed to get SandiaDecayDataBase - nuclide suggestions will not use it." << endl;
    return sm_CharacteristicGammas;
  }
  
  
  while( SpecUtils::safe_get_line( characteristicFile, line ) )
  {
    vector<string> fields;
    SpecUtils::trim( line );
    if( line.empty() )
      continue;
    
    SpecUtils::split( fields, line, " \t" );
    
    double testenergy;
    if( (fields.size() != 2) || !(stringstream(fields[1]) >> testenergy) )
    {
      cerr << "CharacteristicGammas.txt: Failed to read line '" << line << "' - skipping." << endl;
      continue;
    }
    
    const SandiaDecay::Nuclide *testnuc = db->nuclide( fields[0] );
    
    sm_CharacteristicGammas.push_back( {fields[0], testnuc, testenergy} );
  }//while( getline( characteristicFile, line ) )
  
  std::sort( sm_CharacteristicGammas.begin(), sm_CharacteristicGammas.end(),
            []( const tuple<string, const SandiaDecay::Nuclide *, float> &lhs,
               const tuple<string, const SandiaDecay::Nuclide *, float> &rhs ) -> bool {
    return get<2>(lhs) < get<2>(rhs);
  } );
  
  return sm_CharacteristicGammas;
}//const vector<...> &characteristicGammas()



std::mutex sm_PhotPeakLis_mutex;
bool sm_PhotPeakLis_inited = false;

/** Similar to CharacteristicGammas.txt, we will parse PhotoPeak.lis a single time into memory.
 
 Previously, we were looping over the file each time we wanted data from it, and stopping once we
 found the energy we wanted; at an energy of ~1400 keV, it would take like 15ms to find the entry.
 
 Now the upfront parsing into memory takes ~11 ms, but subsequent searches take no time.
 We are not finding the Nuclide, Element, or Reaction pointers, because this then costs ~100ms
 to initialize (although I'm not sure SandiaDecay was being compiled with -O3).
 */
enum class PhotopeakLisType
{
  Xray, Reaction, Nuclide
};//

vector<tuple<float, //Energy, in keV
             float, //BR
             string, //symbol, ex. "Co60", "Fe".  Xrays and reactions give just symbol (it looks like short string optimization always kicks in, so std::string is no cost).
             PhotopeakLisType,
             ReferenceLineInfo::RefLine::RefGammaType>> sm_PhotPeakLis;


template<size_t n>
const char *str_pos( const char *line, const size_t length, const char (&label)[n] )
{
  static_assert( n > 1, "You must call with a non-empty search string" );
  
  const char * const end = line + length;
  const char * const pos = std::search( line, end, label, label + (n-1) );
  
  return (pos == end) ? nullptr : pos;
}


const vector<tuple<float,
                   float,
                   string,
                   PhotopeakLisType,
                   ReferenceLineInfo::RefLine::RefGammaType>> &photopeak_lis()
{
  // Single threaded: Took 10820.5 microseconds to parse PhotoPeak.lis; 367.693 microseconds of this was loading data
                   
  // There were 10567 total entries.
  //const auto start_time = std::chrono::steady_clock::now();
                     
  std::lock_guard<std::mutex> lock( sm_PhotPeakLis_mutex );
  if( sm_PhotPeakLis_inited )
    return sm_PhotPeakLis;
  
  sm_PhotPeakLis_inited = true;
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return sm_PhotPeakLis;
  
  const ReactionGamma *reactionDb = nullptr;
  try
  {
    reactionDb = ReactionGammaServer::database();
  }catch(...)
  {
    cerr << "photopeak_lis(): Failed to open gamma reactions XML file" << endl;
  }
  
  
  const string filename = SpecUtils::append_path( dataDirectory(), "PhotoPeak.lis" );
  
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream input( wfilename.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream input( filename.c_str(), ios_base::binary | ios_base::in );
#endif
  
  if( !input.good() )
  {
    cerr << "photopeak_lis(): Failed to open " << filename << endl;
    return sm_PhotPeakLis;
  }//if( !input.good() )
  

  input.unsetf(ios::skipws);
  
  // Determine stream size
  input.seekg(0, ios::end);
  const size_t size = static_cast<size_t>( input.tellg() );
  input.seekg(0);
  
  if( size < 100*1024 )
  {
    cerr << "photopeak_lis(): Photopeak.lis file smaller than expected - only " << size << " bytes" << endl;
    return sm_PhotPeakLis;
  }
                     
  // Load data
  string data;
  data.resize(size);
  if( !input.read(&data.front(), static_cast<streamsize>(size)) )
  {
    cerr << "photopeak_lis(): Failed to read " << filename << endl;
    return sm_PhotPeakLis;
  }
  
  // We expect 10974 entries, so we'll reserve this space upfront, plus a little
  size_t line_num = 0;
  sm_PhotPeakLis.resize( 10974 + 26 );
                    
#define VALIDATE_PHOTOPEAK_LIS_NUCS 0
                     
  //const auto data_load_time = std::chrono::steady_clock::now();
                     
  bool is_sorted = true;
  float last_energy = -1000;
                     
  // TODO: we could break parsing this file up into multiple threads without much difficulty.
  size_t line_start = 0, line_end = 0;
  for( ; (line_start + 1) < size; line_start = line_end + 1 )
  {
    //Example lines:
    //    "  955.34  AlCapture           0.476  9.2e-05"
    //    "  954.55  I132               17.569  0.066"
    //    " 7463.00  CrCapture de        8.664  3.3e-04"
    //    " 7847.30  N16 de              0.076  1.2e-04"
    //    " 1332.49  Co60               99.983  0.516"
    //    "   57.98  W x-ray            57.838  0.263"
    
    line_end = data.find('\n', line_start + 1);
    if( line_end == string::npos )
      line_end = size;
    
    assert( line_end > line_start );
    assert( (line_end - line_start) < 50 );
    assert( ((line_end - line_start) < 5) || ((line_end - line_start) > 38) );
    
    const size_t line_len = (line_end > line_start) ? (line_end - line_start) : size_t(0);
    if( line_len < 38 )
    {
      cerr << "photopeak_lis(): Found line of length " << line_len << endl;
      continue;
    }
    
    const char * const begin_energy = data.c_str() + line_start;
    const char * const begin_nuc = begin_energy + 10;
    const char * const begin_importance = begin_energy + 28;
    const char * const begin_br = begin_energy + 36;
    const char * const end_line = data.c_str() + line_end - 1; //-1 is for \r
    
    //cout << "Energy='" << string(begin_energy, begin_nuc) << "', nuc='" << string(begin_nuc, begin_importance)
    //<< "', importance='" << string(begin_importance,begin_br) << "', br='" << string(begin_br,end_line) << "'" << endl;
    
    if( line_num >= sm_PhotPeakLis.size() )
      sm_PhotPeakLis.resize( sm_PhotPeakLis.size() + 100 );
    
    auto &line_data = sm_PhotPeakLis[line_num];
    
    if( !SpecUtils::parse_float( begin_energy, 9, get<0>(line_data) ) )
    {
      cerr << "photopeak_lis(): Failed to parse float for line '" << std::string(begin_energy,end_line) << "'" << endl;
      continue;
    }
    
    if( !SpecUtils::parse_float( begin_br, line_len - 36, get<1>(line_data) ) )
    {
      cerr << "photopeak_lis(): Failed to parse BR for line '" << std::string(begin_energy,end_line) << "'" << endl;
      continue;
    }
  
    if( get<0>(line_data) < last_energy )
      is_sorted = false;
    last_energy = get<0>(line_data);
    
    const char * const xray_pos = str_pos( begin_nuc, 16, "x-ray" );
    const char * const capture_pos = str_pos( begin_nuc, 16, "Capture" );
    // There is one 'AlInelastic' and one 'FeInelastic' in the file, but no other 'Inelastic' - so we'll ignore them
    const char * const se_pos = str_pos( begin_nuc, 16, " se" );
    const char * const de_pos = str_pos( begin_nuc, 16, " de" );
    
    
    if( se_pos )
      get<4>(line_data) = ReferenceLineInfo::RefLine::RefGammaType::SingleEscape;
    else if( de_pos )
      get<4>(line_data) = ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape;
    else
      get<4>(line_data) = ReferenceLineInfo::RefLine::RefGammaType::Normal;
    
    if( xray_pos )
    {
      assert( xray_pos >= begin_nuc );
      
      get<2>(line_data) = string(begin_nuc, xray_pos);
      get<3>(line_data) = PhotopeakLisType::Xray;
      
#if( VALIDATE_PHOTOPEAK_LIS_NUCS )
      
      // For debug builds, check we do actually find the element
      const SandiaDecay::Element *el = db->element( get<2>(line_data) );
      // We could check the energy, but its not that important.
      
      if( !el )
        cerr << "photopeak_lis(): Failed to get element for line '" << get<2>(line_data) << "'" << endl;
       
#endif
    }else if( capture_pos )
    {
      if( !reactionDb )
        continue;
      
      get<2>(line_data) = string(begin_nuc, capture_pos);
      get<3>(line_data) = PhotopeakLisType::Reaction;
      
#if( VALIDATE_PHOTOPEAK_LIS_NUCS )
      const SandiaDecay::Element * const el = db->element( string(begin_nuc, capture_pos) );
      if( !el )
      {
        cerr << "photopeak_lis(): Failed to get element for line '" << std::string(begin_energy,end_line) << "'" << endl;
        continue;
      }

      // Implementation to match the gamma up to the reaction
      vector<const ReactionGamma::Reaction *> reactions;
      reactionDb->reactions( el, reactions );
      
      const float energy = get<0>(line_data) + (se_pos ? 510.9989f : (de_pos ? 1021.998f : 0.0f));
      const ReactionGamma::Reaction *rctn = nullptr;
      
      // TODO: Find the reaction with largest BR at this energy - there could be multiple reactions with this energy (e.g., slow and fast share an energy)
      float nearest_de = 1000000.0f;
      for( size_t rctn_index = 0; !rctn && (rctn_index < reactions.size()); ++rctn_index )
      {
        const ReactionGamma::Reaction * const r = reactions[rctn_index];
        
        for( size_t gamma_index = 0; gamma_index < r->gammas.size(); ++gamma_index )
        {
          const ReactionGamma::EnergyAbundance &ea = r->gammas[gamma_index];
          const float de = fabs(ea.energy - energy);
          if( (de < 1.0f) && (de < nearest_de) )
          {
            nearest_de = de;
            rctn = r;
          }
        }
      }//for( const ReactionGamma::Reaction *r : reactions )
      
      if( !rctn  || (nearest_de > 2.0f) )
      {
        // About 405 of 2073 'Capture' lines dont get matches... InterSpec really needs reaction data updated.
        //cerr << "Failed to get reaction for line '" << std::string(begin_energy,end_line) << "'" << endl;
        continue;
      }
      
      //get<4>(line_data) = rctn;
#endif
    }else //if( is probably a nuclide )
    {
      const char * const end_nuc = se_pos ? se_pos : (de_pos ? de_pos : (begin_nuc + 9));
      get<2>(line_data) = string(begin_nuc, end_nuc);
      get<3>(line_data) = PhotopeakLisType::Nuclide;
      
#if( VALIDATE_PHOTOPEAK_LIS_NUCS )
      string nuclide(begin_nuc, end_nuc);
      const SandiaDecay::Nuclide *nuc = db->nuclide( nuclide );
      
      if( !nuc )
      {
        if( nuclide.find("Inelast") == string::npos ) //Dont print error message for the 'AlInelastic' and 'FeInelastic' lines
          cerr << "photopeak_lis(): Failed to get nuclide for line '" << std::string(begin_energy,end_line) << "' - nuclide='" << nuclide << "'" << endl;
      }//
#endif
    }//
    
    line_num += 1;
  }//for( ; (line_start + 1) < size; line_start = line_end + 1 )
  
  sm_PhotPeakLis.resize( line_num );
                     
  assert( is_sorted );
  if( !is_sorted )
  {
    cerr << "PhotoPeak.lis not sorted" << endl;
    std::sort( begin(sm_PhotPeakLis), end(sm_PhotPeakLis),
              []( const tuple<float,float,string,PhotopeakLisType,ReferenceLineInfo::RefLine::RefGammaType> &lhs,
                 const tuple<float,float,string,PhotopeakLisType,ReferenceLineInfo::RefLine::RefGammaType> &rhs) -> bool {
      
    } );
  }//if( !is_sorted )
                     
  //const auto end_time = std::chrono::steady_clock::now();
  //cout << "Took " << std::chrono::duration<double, std::micro>(end_time - start_time).count()
  //     << " microseconds to parse PhotoPeak.lis; "
  //     << std::chrono::duration<double, std::micro>(data_load_time - start_time).count()
  //     << " microseconds of this was loading data"
  //     << endl;
  //cout << "There were " << sm_PhotPeakLis.size() << " total entries." << endl;
  
  return sm_PhotPeakLis;
}//photopeak_lis(...)

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
  if( peak->continuum()->parametersProbablySet() )
  {
    double lowx(0.0), upperx(0.0);
    findROIEnergyLimits( lowx, upperx, *peak, data );
    contArea = peak->offset_integral( lowx, upperx, data );
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
  const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &characlines = characteristicGammas();
  
  map<const SandiaDecay::Nuclide *, int> answer;
  
  for( std::shared_ptr<const PeakDef> peak : *all_peaks )
  {
    const double lowe = (peak->gausPeak() ? (peak->mean()-1.5*peak->sigma()) : peak->lowerX());
    const double highe = (peak->gausPeak() ? (peak->mean() + 1.5*peak->sigma()) : peak->upperX());
    for( const tuple<string, const SandiaDecay::Nuclide *, float> &enp : characlines )
    {
      const float energy = get<2>(enp);
      const SandiaDecay::Nuclide *nuc = get<1>(enp);
      
      if( nuc && (energy >= lowe) && (energy <= highe) )
      {
        if( answer.find(nuc) == answer.end() )
          answer[nuc] = 0;
        answer[nuc] = answer[nuc] + 1;
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
    if( energy < 1.0 )
      continue;
    
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
  
  // The characteristics(...) gives weight based on number of peaks that
  // nuclides in CharacteristicGammas.txt - it doesnt single out the single
  // peak we are interested in.
  map<const SandiaDecay::Nuclide *, int> charclines = characteristics( peaks );
  for( NuclideStatWeightPair &n : candidates )
  {
    if( charclines.count(n.nuclide) )
      n.weight *= 2.0*charclines[n.nuclide];
  }
  
  // Here we will check if our peak of interest matches one of the lines in CharacteristicGammas.txt
  //  And if so, we'll double its weight, or at least make sure the nuclide is a candidate
  //  TODO: move this logic to #characteristics(...), and maybe improve the logic to have better fractional weighting.
  {// Begin check if `peak` corresponds to nuclide in CharacteristicGammas.txt, if so
    double max_weight = 0.0;
    for( NuclideStatWeightPair &n : candidates )
      max_weight = std::max( max_weight, n.weight );
    
    const float lower_e = static_cast<float>(energy - 2.0*sigma);
    const float upper_e = static_cast<float>(energy + 2.0*sigma);
   
    // We could use lower/upper_bound, since characteristicGammas() is sorted by energy,
    //  but its few enough elements I wont bother with this at the moment.
    map<const SandiaDecay::Nuclide *, size_t> n_times_nuc;
    for( const tuple<string, const SandiaDecay::Nuclide *, float> &line : characteristicGammas() )
    {
      const SandiaDecay::Nuclide * const nuc = get<1>(line);
      if( n_times_nuc.count(nuc) )
        n_times_nuc[nuc] = 0;
      n_times_nuc[nuc] += 1;
      
      const float line_energy = get<2>(line);
      if( nuc && (line_energy >= lower_e) && (line_energy <= upper_e) )
      {
        bool found_candidate = false;
        for( NuclideStatWeightPair &n : candidates )
        {
          if( n.nuclide == nuc )
          {
            found_candidate = true;
            n.weight += 0.5*max_weight; //the 0.5 and max_weight are arbitrarily decided, without testing
            break;
          }
        }
        
        if( !found_candidate )
        {
          NuclideStatWeightPair nuswp;
          nuswp.nuclide = nuc;
          nuswp.weight = 0.5*max_weight / n_times_nuc[nuc]; //Weight is arbitrarily decided, without testing
          candidates.push_back( std::move(nuswp) );
        }
      }//if( our peak matches this characteristic line )
    }//for( loop over CharacteristicGammas.txt )
  }// End check if `peak` corresponds to nuclide in CharacteristicGammas.txt, if so
  
  
  //Prefer isotopes with longer half lives
  //I want some thing like:
  //  n.weight *= (log(1.0 + n.nuclide->halfLife) / log(max_hl));
  //but in order to get Th232 instead of U232 when only the 2614 peak is
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
  
  // cerr << "For mean=" << peak->mean() << "keV, top candiates are:" << endl;
  // for( size_t i = 0; i < candidates.size() && i < 10; ++i  )
  //   cerr << "\t" << i << "\t" << candidates[i].nuclide->symbol
  //        << "\tw=" << candidates[i].weight << endl;
  // cerr << endl << endl;
  
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
      
      if( !PeakFitUtils::is_high_res(data) )
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
  const vector<tuple<float,float,string,PhotopeakLisType,ReferenceLineInfo::RefLine::RefGammaType>> &lines = photopeak_lis();
 
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
  
  const double mean = peak->mean();
  const double sigma = peak->gausPeak() ? peak->sigma() : 0.125*peak->roiWidth();
  const double nsigma = (0.03 < (sigma/mean)) ? 1.25 : 2.5;
  double width = nsigma*sigma;
  if( mean > 3200.0 )
    width = std::max( width, 20.0 );
  if( mean > 4500.0 )
    width = std::max( width, 30.0 );
  
  const float minenergy = static_cast<float>( peak->mean() - width );
  const float maxenergy = static_cast<float>( peak->mean() + width );
  
  //First place the results into a map, keyed off of difference in energy from
  //  the peak mean, so this way results will be returned, ordered by how close
  //  the gamma is to the peak mean.
  map<double,vector<string> > results;
  
  string line;
  char buff[64];
  
  auto compare = []( const tuple<float,float,string,PhotopeakLisType,ReferenceLineInfo::RefLine::RefGammaType> &line, const float energy ) ->bool {
    return (get<0>(line) < energy);
  };
  
  auto iter = std::lower_bound( begin(lines), end(lines), minenergy, compare );
  
  //cout << "Start val=" << std::get<0>(*iter) << " keV for " << minenergy << " keV - peak->mean()=" << peak->mean() << endl;
  
  for( ; (iter != end(lines)) && (get<0>(*iter) <= maxenergy); ++iter )
  {
    const tuple<float,float,string,PhotopeakLisType,ReferenceLineInfo::RefLine::RefGammaType> &line = *iter;
    
    const string nucstr = SpecUtils::trim_copy( get<2>(line) ); // I dont think the trim is really necassary, but whatever
    const double diff = fabs( mean - get<0>(line) );
    
    double energy = get<0>(line);
    const char *escape_str = "";
    switch( get<4>(line) )
    {
      case ReferenceLineInfo::RefLine::RefGammaType::Normal:
      case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
      case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
      case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
        break;
      case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
        energy += 510.9989;
        escape_str = "S.E. ";
        break;
        
      case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
        energy += 1021.998f;
        escape_str = "D.E. ";
        break;
    }//switch( get<4>(line) )
    
    
    switch( get<3>(line) )
    {
      case PhotopeakLisType::Xray:
        results[diff].push_back( nucstr );
        break;
        
      case PhotopeakLisType::Reaction:
      {
        if( !reactionDb )
          continue;
      
        const SandiaDecay::Element *el = db->element( nucstr );
        if( !el )
        {
          cerr << "Couldnt get element for element '" << nucstr << "'" << endl;
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
              snprintf( buff, sizeof(buff), "%s %s%.2f keV",
                       rctn->name().c_str(), escape_str, g.energy );
              
              //due to the nested loops, and the non-uniquness of only matching
              //  down to 2 keV, the same reaction/gamma may be inserted into
              //  the result more than once; this is protected against the user
              //  seeing it when all results are collated together.
              results[diff].push_back( buff );
            }//if( close enough energy )
          }//for( const ReactionGamma::EnergyAbundance &g : gammas )
        }//for( const ReactionGamma::Reaction *rctn : reactions )
      
        break;
      }//case PhotopeakLisType::Reaction:
        
      case PhotopeakLisType::Nuclide:
      {
        const SandiaDecay::Nuclide *nuc = db->nuclide( nucstr );
        if( nuc )
        {
          //Energies in PhotoPeak.lis may not exactly line up with the energies in SandiaDecay
          size_t index = 0;
          const SandiaDecay::Transition *trans = NULL;
          PeakDef::SourceGammaType type;
          PeakDef::findNearestPhotopeak( nuc, energy, -1.0, false, trans, index, type );
          if( trans )
          {
            energy = trans->products[index].energy;
          }else
          {
#if( PERFORM_DEVELOPER_CHECKS )
            //log_developer_error( __func__, "..." ); this is a bit more than is necassary - we know there will be missing reactions
            cerr << "Couldnt find gamma at " << energy <<" keV for " << nucstr
            << " in SandiaDecay database" << endl;
#endif
            continue;
          }
          
          snprintf( buff, sizeof(buff), "%s %s%.2f keV", nucstr.c_str(), escape_str, energy );
          results[diff].push_back( buff );
        }
        break;
      }//case PhotopeakLisType::Nuclide:
    }//switch( get<3>(line) )
  }//for( ; (iter != end(lines)) && (get<0>(*iter) <= maxenergy); ++iter )
  
  // We can easily get over 100 results, which after the first few, are mostly useless - so we
  //  will arbitrary limit to about 25 results.
  //  Alternatively we could limit by weight-value, but I havent investigated reasonable ways
  //  to do this.
  const size_t max_characteristic_nucs = 25;
  
  for( const auto &val : results )
  {
    for( const string &nuc : val.second )
      characteristicnucs.push_back( nuc );
    
    // Note that if a nuclide had exact same weight, we could have a few more results
    //  than max_characteristic_nucs
    if( characteristicnucs.size() > max_characteristic_nucs )
      break;
  }//for( const auto &val : results )
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
    switch( ref.m_source_type )
    {
      case ReferenceLineInfo::SourceType::Nuclide:
      case ReferenceLineInfo::SourceType::FluorescenceXray:
      case ReferenceLineInfo::SourceType::Reaction:
      case ReferenceLineInfo::SourceType::Background:
        break;
        
      case ReferenceLineInfo::SourceType::CustomEnergy:
      case ReferenceLineInfo::SourceType::None:
        continue;
        break;
    }//switch( ref.m_source_type )
    
  
    vector<pair<string,double>> trials;
    
    //We will take the simeple approach and assume transmision and detection efficiency
    //  are about the same for all lines we care about, and create a simple
    //  weight of sf/(0.25*sigma + distance)
    for( const ReferenceLineInfo::RefLine &line : ref.m_ref_lines )
    {
      double energy = line.m_energy;
      const double distance = fabs( energy - mean );
      
      if( distance > 4.0*sigma )
        continue;
      
      const double amp = line.m_normalized_intensity;
      const double weight = amp / (0.25*sigma + distance);
      
      string source = ref.m_input.m_input_txt;
      
      //if( ref.m_source_type == ReferenceLineInfo::SourceType::Background )
      //{
      if( line.m_parent_nuclide )
        source = line.m_parent_nuclide->symbol;
      else if( line.m_element )
        source = line.m_element->symbol;
      else if( line.m_reaction )
        source = line.m_reaction->name();
      else
        continue; //We want to be able to assign peaks to a source
      //}//if( ref.m_source_type == ReferenceLineInfo::SourceType::Background )
      
      string typestr;
      
      switch( line.m_particle_type )
      {
        case ReferenceLineInfo::RefLine::Particle::Alpha:
        case ReferenceLineInfo::RefLine::Particle::Beta:
          // We dont want alpha or beta particles
          continue;
          break;
          
        case ReferenceLineInfo::RefLine::Particle::Gamma:
          break;
          
        case ReferenceLineInfo::RefLine::Particle::Xray:
          typestr = " xray";
          break;
      }//switch( line.m_particle_type )
      
      
      switch( line.m_source_type )
      {
        case ReferenceLineInfo::RefLine::RefGammaType::Normal:
        case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
          break;
          
        case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
          typestr = " S.E.";
          energy += 510.99891;
          break;
        case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
          typestr = " D.E.";
          energy += 2.0*510.99891;
          break;
          
        case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
        case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
          continue;
      }//switch( line.m_source_type )
      
      snprintf( buffer, sizeof(buffer), "%s%s %.2f keV",
               source.c_str(), typestr.c_str(), energy );
      
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
  
  //Artificially bump up the rating for background and common sources; note that
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
        //The opening and closing parenthesis as first/last character will cause
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
      //cerr << "Have sum for peak " << rpeak.get() << endl;
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
        //cerr << "Have double sum for peak " << rpeak.get() << " and " << lpeak.get() << endl;
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
