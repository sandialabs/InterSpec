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

#include <set>
#include <deque>
#include <tuple>
#include <thread>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <type_traits>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "Wt/WDateTime"
#include "Wt/WApplication"

#include "ceres/ceres.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
const double ns_decay_act_mult = SandiaDecay::MBq;

struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct


// We'll define some XML helper functions for serializing to/from XML
//  (we should probably put these in a header somewhere and use through the code to minimize
//   things a little)
void append_float_node( rapidxml::xml_node<char> *base_node, const char *node_name, const double value )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  
  char buffer[128];
  snprintf( buffer, sizeof(buffer), "%1.8e", value );
  const char *strvalue = doc->allocate_string( buffer );
  xml_node<char> *node = doc->allocate_node( node_element, node_name, strvalue );
  base_node->append_node( node );
}//append_float_node(...)


void append_bool_node( rapidxml::xml_node<char> *base_node, const char *node_name, const bool value )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  xml_node<char> *node = doc->allocate_node( node_element, node_name, (value ? "true" : "false") );
  base_node->append_node( node );
}//append_bool_node(...)


template <class T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr >
void append_int_node( rapidxml::xml_node<char> *base_node, const char *node_name, const T value )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  const char *strvalue = doc->allocate_string( std::to_string(value).c_str() );
  xml_node<char> *node = doc->allocate_node( node_element, node_name, strvalue );
  base_node->append_node( node );
}//append_bool_node(...)


rapidxml::xml_node<char> *append_string_node( rapidxml::xml_node<char> *base_node, const char *node_name, const string &value )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  const char *strvalue = doc->allocate_string( value.c_str(), value.size() + 1 );
  xml_node<char> *node = doc->allocate_node( node_element, node_name, strvalue );
  base_node->append_node( node );
  
  return node;
}//append_string_node(...)


void append_attrib( rapidxml::xml_node<char> *base_node, const string &name, const string &value )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  
  const char *name_str = doc->allocate_string( name.c_str() );
  const char *value_str = doc->allocate_string( value.c_str() );
  
  xml_attribute<char> *attrib = doc->allocate_attribute( name_str, value_str );
  base_node->append_attribute( attrib );
}//void append_attrib(...)


void append_version_attrib( rapidxml::xml_node<char> *base_node, const int version )
{
  char buffer[32];
  snprintf( buffer, sizeof(buffer), "%i", version );
  append_attrib( base_node, "version", buffer );
}//append_version_attrib(...)


template<size_t n>
double get_float_node_value( const rapidxml::xml_node<char> * const parent_node, const char (&name)[n] )
{
  assert( parent_node );
  assert( name );
  
  if( !parent_node )
    throw runtime_error( "null parent node." );
  
  const rapidxml::xml_node<char> *node = XML_FIRST_NODE(parent_node, name);
  if( !node )
    throw runtime_error( "Missing node '" + std::string(name) + "'" );
  
  const string value = SpecUtils::xml_value_str(node);
  double answer;
  if( !(stringstream(value) >> answer) )
    throw runtime_error( "Value ('" + value + "') of node '"
                        + string(name) + "' was a valid float." );
  
  return answer;
}//double get_float_node_value(...)


template<size_t n>
bool get_bool_node_value( const rapidxml::xml_node<char> * const parent_node, const char (&name)[n] )
{
  assert( parent_node );
  assert( name );
  
  if( !parent_node )
    throw runtime_error( "null parent node." );
  
  const rapidxml::xml_node<char> *node = XML_FIRST_NODE(parent_node, name);
  if( !node )
    throw runtime_error( "Missing node '" + std::string(name) + "'" );
  
  if( XML_VALUE_ICOMPARE(node,"yes") || XML_VALUE_ICOMPARE(node,"true") || XML_VALUE_ICOMPARE(node,"1") )
    return true;
  
  if( !XML_VALUE_ICOMPARE(node,"no") && !XML_VALUE_ICOMPARE(node,"false") && !XML_VALUE_ICOMPARE(node,"0") )
    throw runtime_error( "Invalid boolean value in node '" + string(name) + "' with value '"
                        + SpecUtils::xml_value_str(node) );
  
  return false;
}//bool get_bool_node_value(...)


void check_xml_version( const rapidxml::xml_node<char> * const node, const int required_version )
{
  assert( node );
  const rapidxml::xml_attribute<char> *att = XML_FIRST_ATTRIB(node, "version");
  
  int version;
  if( !att || !att->value()
     || (sscanf(att->value(), "%i", &version) != 1) )
    throw runtime_error( "invalid or missing version" );
  
  if( (version < 0) || (version > required_version) )
    throw runtime_error( "Invalid version: " + std::to_string(version) + ".  "
                        + "Only up to version " + to_string(required_version)
                        + " supported." );
}//check_xml_version(...)

template<size_t n>
const rapidxml::xml_node<char> *get_required_node( const rapidxml::xml_node<char> *parent, const char (&name)[n] )
{
  assert( parent );
  const auto child_node = XML_FIRST_INODE(parent, name);
  if( !child_node )
    throw runtime_error( "No <" + string(name) + "> node" );
  
  return child_node;
}//get_required_node(...)


struct RoiRangeChannels : public RelActCalcAuto::RoiRange
{
  // TODO: currently if the energy calibration is adjusted - we have to keep the total number of channels (i.e. the number of residuals constant) we have to move the first/last channel the same number of channels
  size_t num_channels;
  
  
  /** Returns the the channel range corresponding to the desired energy range.
   
   If an energy extends more than halfway into a channel, the entire channel will be used.
   
   @param wanted_nchan The number of channels wanted in the range; if the nominal channel
          range doesnt match this, the channels will be adjusted so they will.  If zero, no
          adjustment will be made.
   
   @returns The lower and upper channels for the energy range; its guaranteed the \c .first element
            will be less than or equal to the \c .second element.
   */
  static pair<size_t,size_t> channel_range( const double lower_energy, const double upper_energy,
                                          const size_t wanted_nchan,
                                          const shared_ptr<const SpecUtils::EnergyCalibration> &cal )
  {
    if( !cal || !cal->valid() || (cal->num_channels() < 16) )
      throw runtime_error( "RoiRangeChannels: invalid energy calibration." );
    
    if( upper_energy <= lower_energy )
      throw runtime_error( "RoiRangeChannels: energy range invalid." );
    
    const size_t nchannel = cal->num_channels();
    
    if( wanted_nchan > nchannel )
      throw runtime_error( "RoiRangeChannels: invalid wanted number of channels." );
    
    const double first_channel_f = cal->channel_for_energy( lower_energy );
    const double last_channel_f = cal->channel_for_energy( upper_energy );
    
    double first_channel_i, last_channel_i;
    double first_fractional = modf( first_channel_f, &first_channel_i );
    double last_fractional = modf( last_channel_f, &last_channel_i );
    
    if( first_fractional >= 0.5 )
    {
      first_channel_i += 1.0;
      first_fractional -= 1.0;
    }
    
    if( last_fractional >= 0.5 )
    {
      last_channel_i += 1.0;
      last_fractional -= 1.0;
    }
    
    first_channel_i = std::max( first_channel_i, 0.0 );
    last_channel_i = std::max( last_channel_i, 0.0 );
    
    size_t first_channel = std::min( static_cast<size_t>(first_channel_i), nchannel - 1 );
    size_t last_channel = std::min(static_cast<size_t>(last_channel_i), nchannel - 1 );
    
    if( last_channel < first_channel )
      throw runtime_error( "Starting and ending channels for "
                          + std::to_string(lower_energy) + " to "
                          + std::to_string(upper_energy)
                          + " keV energy range are inconsistent [" + std::to_string(first_channel)
                          + ", " + std::to_string(last_channel) + "]." );
    
    const size_t niave_nchan = 1 + last_channel - first_channel;
    if( wanted_nchan && (niave_nchan != wanted_nchan) )
    {
      // TODO: right now we are just adding/subtracting the difference of channels from one side
      //       which should be fine if we are only off by one channel - which is probably all we
      //       expect to encounter, but
      const size_t diff = (niave_nchan > wanted_nchan) ? (niave_nchan - wanted_nchan)
                                                       : (wanted_nchan - niave_nchan);
      if( niave_nchan < wanted_nchan )
      {
        // We need to add `diff` channels to our channel range
        if( last_fractional > first_fractional )
        {
          last_channel += diff;
        }else if( first_channel >= diff )
        {
          first_channel -= diff;
        }else
        {
          assert( diff > first_channel );
          last_channel += diff - first_channel;
          first_channel = 0;
        }
      }else
      {
        // We need to subtract `diff` channels to our channel range
        if( first_fractional < last_fractional )
        {
          last_channel -= diff;
        }else if( first_channel >= diff )
        {
          first_channel -= diff;
        }else
        {
          assert( diff > first_channel );
          last_channel += diff - first_channel;
          first_channel = 0;
        }
      }// if( niave_nchan < wanted_nchan ) / else
      
      if( last_channel >= nchannel )
      {
        first_channel -= ((1 + last_channel) - nchannel);
        last_channel = nchannel - 1;
      }
    }//if( we didnt have the number of channels we wanted )
    
    // Some sanity checks on the above logic
    assert( first_channel <= last_channel );
    assert( !wanted_nchan || (((1 + last_channel) - first_channel) == wanted_nchan) );
    
    return { first_channel, last_channel };
  }//channel_range(...)
  
  RoiRangeChannels( const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                    const RelActCalcAuto::RoiRange &roi_range )
  : RelActCalcAuto::RoiRange( roi_range )
  {
    auto range = channel_range( roi_range.lower_energy, roi_range.upper_energy, 0, energy_cal );
    
    num_channels =  1 + range.second - range.first;
  }//RoiRangeChannels constructor
};//struct RoiRangeChannels


struct NucInputGamma : public RelActCalcAuto::NucInputInfo
{
  struct EnergyYield
  {
    double energy;
    double yield;
  };
  
  vector<EnergyYield> nominal_gammas;
  
  static size_t remove_gamma( const double energy, vector<SandiaDecay::EnergyRatePair> &gammas )
  {
    const size_t ninitial = gammas.size();
    
    gammas.erase( std::remove_if(begin(gammas), end(gammas),
      [energy](const SandiaDecay::EnergyRatePair &e ){
        return fabs(e.energy - energy) < 0.001;
    }), end(gammas) );
    
    return gammas.size() - ninitial;
  }
  
  NucInputGamma( const RelActCalcAuto::NucInputInfo &info )
   : RelActCalcAuto::NucInputInfo( info )
  {
    if( !nuclide )
      throw runtime_error( "NucInputGamma: null Nuclide." );
    
    if( age < 0.0 )
      throw runtime_error( "NucInputGamma: age may not be negative (" + nuclide->symbol + ": "
                           + PhysicalUnits::printToBestTimeUnits(age)  + ")" );
    
    double nominal_age = info.age;
    
    //if( !fit_age )
    //  nominal_age = PeakDef::defaultDecayTime( nuclide, nullptr );
    
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuclide, ns_decay_act_mult, nominal_age );
    auto gammas = mix.gammas( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
    
    for( double energy : gammas_to_exclude )
    {
      const size_t did_remove = remove_gamma( energy, gammas );
      assert( did_remove );
    }
    
    nominal_gammas.resize( gammas.size() );
    for( size_t i = 0; i < gammas.size(); ++i )
    {
      nominal_gammas[i].energy = gammas[i].energy;
      nominal_gammas[i].yield = gammas[i].numPerSecond / ns_decay_act_mult;
    }
  }//NucInputGamma constructor
};//struct NucInputGamma



struct RelActAutoCostFcn /* : ROOT::Minuit2::FCNBase() */
{
  RelActCalcAuto::Options m_options;
  std::vector<NucInputGamma> m_nuclides;
  std::vector<RoiRangeChannels> m_energy_ranges;
  std::vector<RelActCalcAuto::FloatingPeak> m_extra_peaks;
  
  const std::shared_ptr<const SpecUtils::Measurement> m_spectrum;
  
  const float m_live_time;
  const vector<float> m_channel_counts; //background subtracted channel counts
  const vector<float> m_channel_count_uncerts; //e.g. sqrt( m_channel_counts[i] )
  const std::shared_ptr<const SpecUtils::EnergyCalibration> m_energy_cal;
  
  /** We want to anchor the relative efficiency curve to be equal to 1.0 at the lowest energy; to
   do this we will add an extra residual term that is this value, multiplied by the difference
   in the relative efficiency curve, from 1.0, at the lowest energy ranges lowest energy.
   
   In the end, I dont think the answer, or any of the uncertainties, will be dependent on the value
   of this value chosen, since this is just eliminating an ambiguity between the normalization of
   the relative efficiency curve, and the relative activities, but I do need to double check on
   this.
   */
  double m_rel_eff_anchor_enhancement;
  
  /** Will either be null, or have FWHM info. */
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  
  /** just for debug purposes, we'll keep track of how many times the eval function gets called. */
  mutable std::atomic<size_t> m_ncalls;
  
  
  RelActAutoCostFcn( RelActCalcAuto::Options options,
                     vector<RelActCalcAuto::RoiRange> energy_ranges,
                     vector<RelActCalcAuto::NucInputInfo> nuclides,
                     vector<RelActCalcAuto::FloatingPeak> extra_peaks,
                     shared_ptr<const SpecUtils::Measurement> spectrum,
                     const vector<float> &channel_counts_uncert,
                     std::shared_ptr<const DetectorPeakResponse> drf,
                     std::vector<std::shared_ptr<const PeakDef>> all_peaks
                           )
  : m_options( options ),
  m_nuclides{},
  m_energy_ranges{},
  m_extra_peaks( extra_peaks ),
  m_spectrum( spectrum ),
  m_live_time( spectrum ? spectrum->live_time() : 0.0f ),
  m_channel_counts( (spectrum && spectrum->gamma_counts()) ? *spectrum->gamma_counts() : vector<float>() ),
  m_channel_count_uncerts( channel_counts_uncert ),
  m_energy_cal( spectrum ? spectrum->energy_calibration() : nullptr ),
  m_rel_eff_anchor_enhancement( 1000.0 ),
  m_drf( nullptr ),
  m_ncalls( 0 )
  {
    if( !spectrum || (spectrum->num_gamma_channels() < 128) )
      throw runtime_error( "RelActAutoCostFcn: invalid spectrum." );
    
    if( !m_energy_cal || !m_energy_cal->valid() )
      throw runtime_error( "RelActAutoCostFcn: invalid energy calibration." );
    
    if( spectrum->num_gamma_channels() != channel_counts_uncert.size() )
      throw runtime_error( "RelActAutoCostFcn: Diff number of spectrum channels and channel counts uncert." );
    
    if( nuclides.empty() )
      throw runtime_error( "RelActAutoCostFcn: no nuclides specified." );
    
    for( const auto &n : nuclides )
    {
      if( !n.nuclide )
        throw runtime_error( "RelActAutoCostFcn: null Nuclide." );
      
      for( const auto &pn : m_nuclides )
      {
        if( n.nuclide == pn.nuclide )
          throw runtime_error( "RelActAutoCostFcn: duplicate nuclide (" + n.nuclide->symbol + ")." );
      }
      
      m_nuclides.emplace_back( n );
    }//for( const auto &n : nuclides )
    
    
    if( !drf || !drf->hasResolutionInfo() )
    {
      std::shared_ptr<DetectorPeakResponse> new_drf;
      if( drf && drf->isValid() )
        new_drf = make_shared<DetectorPeakResponse>( *drf );
      else
        new_drf = make_shared<DetectorPeakResponse>( "FLAT", "FLAT" );
      
      try
      {
        const std::vector<float> drf_coefs{ 0.0f, 0.0f, 0.0f, 0.0f }, uncerts;
        new_drf->fromExpOfLogPowerSeriesAbsEff( drf_coefs, uncerts,
                                               25*PhysicalUnits::cm,
                                               2*PhysicalUnits::cm,
                                               PhysicalUnits::keV,
                                               spectrum->gamma_energy_min(),
                                               spectrum->gamma_energy_max() );
        
        if( all_peaks.empty() )
          throw runtime_error( "No peaks provided to fit for FWHM parameters." );
        
        auto peaks_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(all_peaks), end(all_peaks) );
        
        new_drf->fitResolution( peaks_deque, spectrum, DetectorPeakResponse::kGadrasResolutionFcn );
        
        drf = new_drf;
      }catch( std::exception &e )
      {
        cerr << "RelActAutoCostFcn: error fitting FWHM for DRF: " << e.what() << endl;
        
        drf.reset();
      }//try / catch (setup drf FWHM info )
        
    }//if( we dont have fwhm info )
    
    m_drf = drf;
    
    
    //Need to initialize m_energy_ranges
    if( energy_ranges.empty() )
      throw runtime_error( "RelActAutoCostFcn: no energy ranges defined." );
    
    const bool highres = PeakFitUtils::is_high_res( spectrum );
    
    for( size_t roi_index = 0; roi_index < energy_ranges.size(); ++roi_index )
    {
      const RelActCalcAuto::RoiRange &roi_range = energy_ranges[roi_index];
      
      if( (roi_range.lower_energy >= roi_range.upper_energy)
         || (roi_range.lower_energy < 0.0) )
        throw runtime_error( "RelActAutoCostFcn: mal-formed energy range" );
      
      // We'll check the ranges dont overlap (but they can touch)
      for( size_t other_roi = roi_index + 1; other_roi < energy_ranges.size(); ++other_roi )
      {
        const RelActCalcAuto::RoiRange &other_roi_range = energy_ranges[other_roi];
        
        if( (roi_range.lower_energy < other_roi_range.upper_energy)
           && (other_roi_range.lower_energy < roi_range.upper_energy) )
        {
          throw runtime_error( "RelActAutoCostFcn: input energy ranges are overlapping ["
                              + std::to_string(roi_range.lower_energy) + ", "
                              + std::to_string(roi_range.upper_energy) + "] and ["
                              + std::to_string(other_roi_range.lower_energy) + ", "
                              + std::to_string(other_roi_range.upper_energy) + "]."
                              );
        }
      }//for( loop over ROIs that come after roi_range )

      if( roi_range.force_full_range )
      {
        assert( !roi_range.allow_expand_for_peak_width );
        if( roi_range.allow_expand_for_peak_width )
          throw runtime_error( "RelActAutoCostFcn: RoiRange::force_full_range and RoiRange::allow_expand_for_peak_width can not both be true." );
        
        m_energy_ranges.emplace_back( m_energy_cal, roi_range );
        continue;
      }//if( roi_range.force_full_range )
      
      //We'll try to limit/break-up energy ranges
      const double min_br = numeric_limits<double>::min();  //arbitrary
      const double num_fwhm_roi = 2.5; // arbitrary...
      
      vector<pair<double,double>> gammas_in_range;
      
      // Define a helper function to add a gamma at an energy into \c gammas_in_range
      auto add_peak_to_range = [&gammas_in_range, this, highres, &roi_range, num_fwhm_roi]( const double energy ){
        double energy_sigma;
        float min_sigma, max_sigma;
        expected_peak_width_limits( energy, highres, min_sigma, max_sigma );
        
        if( m_drf && m_drf->hasResolutionInfo() )
        {
          energy_sigma = m_drf->peakResolutionSigma(energy);
          
          // A sanity check... maybe we dont want this?
          if( energy_sigma < min_sigma )
            energy_sigma = min_sigma;
          if( energy_sigma > max_sigma )
            energy_sigma = max_sigma;
        }else
        {
          energy_sigma = max_sigma;
        }
        
        double gamma_row_lower = energy - num_fwhm_roi*energy_sigma;
        double gamma_row_upper = energy + num_fwhm_roi*energy_sigma;
        
        if( !roi_range.allow_expand_for_peak_width )
        {
          gamma_row_lower = std::max( gamma_row_lower, roi_range.lower_energy );
          gamma_row_upper = std::min( gamma_row_upper, roi_range.upper_energy );
        }
        
        gammas_in_range.push_back( {gamma_row_lower,gamma_row_upper} );
      };//add_peak_to_range(...) lamda
      
      for( const auto &n : m_nuclides )
      {
        for( const auto &g : n.nominal_gammas )
        {
          const double energy = g.energy;
          const double yield = g.yield;
          
          if( (yield > min_br)
             && (energy >= roi_range.lower_energy)
             && (energy <= roi_range.upper_energy) )
          {
            add_peak_to_range( g.energy );
          }
        }//for( const auto &g : n.nominal_gammas )
        
        for( const auto &peak : m_extra_peaks )
        {
          if( (peak.energy >= roi_range.lower_energy) && (peak.energy <= roi_range.upper_energy) )
            add_peak_to_range( peak.energy );
        }//for( const auto &peak : m_extra_peaks )
      }//for( const auto &n : m_nuclides )
      
      std::sort( begin(gammas_in_range), end(gammas_in_range) );
      
      for( size_t index = 0; index < gammas_in_range.size();  )
      {
        const size_t start_index = index;
        for( index += 1; index < gammas_in_range.size(); ++index )
        {
          if( gammas_in_range[index - 1].second < gammas_in_range[index].first )
            break;
        }
        
        RelActCalcAuto::RoiRange this_range = roi_range;
        this_range.lower_energy = gammas_in_range[start_index].first;
        this_range.upper_energy = gammas_in_range[index - 1].second;
        
        // TODO: its possible the the channels of the ranges could overlap - right now the overlapping channel will be double counted; should fix.
        
        m_energy_ranges.emplace_back( m_energy_cal, this_range );
      }//for( loop over gammas_in_range )
    }//for( const RelActCalcAuto::RoiRange &input : energy_ranges )
    
    
    // Check to make sure the floating peak is in a ROI range
    for( const RelActCalcAuto::FloatingPeak &peak : m_extra_peaks )
    {
      bool in_a_range = false;
      for( const RoiRangeChannels &r : m_energy_ranges )
      {
        if( (peak.energy >= r.lower_energy) && (peak.energy <= r.upper_energy) )
        {
          in_a_range = true;
          break;
        }//if( `peak` is in a energy range )
      }//for( loop over energy ranges )
      
      if( !in_a_range )
        throw runtime_error( "Free floating peak at " + std::to_string(peak.energy) + " is not in a ROI." );
    }//for( const auto &peak : m_extra_peaks )
    
    
    // TODO: Figure out a proper value to set m_rel_eff_anchor_enhancement to, like maybe largest peak area divided m_live_time, or something like that
    m_rel_eff_anchor_enhancement = 0.0;
    for( const auto &r : m_energy_ranges )
      m_rel_eff_anchor_enhancement += spectrum->gamma_integral( r.lower_energy, r.upper_energy );
  }//RelActAutoCostFcn constructor.
  
  
  
  size_t number_parameters() const
  {
    size_t num_pars = 0;
    
    // Energy calibration; we will always have these, even if fixed values
    num_pars += 2;
    
    // The FWHM equation
    num_pars += num_parameters( m_options.fwhm_form );
    
    // The Relative Eff coefficients
    num_pars += (m_options.rel_eff_eqn_order + 1);
    
    // The Activities; one parameter for activity, one for age (which may be fixed)
    num_pars += 2*m_nuclides.size();
    
    // Floating peaks; one parameter for amplitude, one for FWHM (which will usually be unused)
    num_pars += 2*m_extra_peaks.size();
    
    // Anything else?
    
    return num_pars;
  }//number_parameters()
  
  size_t number_residuals() const
  {
    // Number of gamma channels, plus one to anchor relative eff
    size_t num_resids = 1;
    for( const auto &r : m_energy_ranges )
      num_resids += r.num_channels;
    
    return num_resids;
  }//size_t number_residuals() const
  
  /** Solve the problem, using the Ceres optimizer. */
  static RelActCalcAuto::RelActAutoSolution solve_ceres( RelActCalcAuto::Options options,
                                                        std::vector<RelActCalcAuto::RoiRange> energy_ranges,
                                                        std::vector<RelActCalcAuto::NucInputInfo> nuclides,
                                                        std::vector<RelActCalcAuto::FloatingPeak> extra_peaks,
                                                        std::shared_ptr<const SpecUtils::Measurement> foreground,
                                                        std::shared_ptr<const SpecUtils::Measurement> background,
                                                        const std::shared_ptr<const DetectorPeakResponse> input_drf,
                                                        std::vector<std::shared_ptr<const PeakDef>> all_peaks
                                                        )
  {
    const auto start_time = std::chrono::high_resolution_clock::now();
    
    const bool highres = PeakFitUtils::is_high_res( foreground );
    
    RelActCalcAuto::RelActAutoSolution solution;
    
    DoWorkOnDestruct setFinalTime( [&solution,start_time](){
      const auto end_time = std::chrono::high_resolution_clock::now();
      solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
    });
    
    solution.m_foreground       = foreground;
    solution.m_background       = background;
    solution.m_options          = options;
    solution.m_rel_eff_form     = options.rel_eff_eqn_type;
    solution.m_fwhm_form        = options.fwhm_form;
    solution.m_input_roi_ranges = energy_ranges;
    if( input_drf && input_drf->isValid() && input_drf->hasResolutionInfo() )
      solution.m_drf = input_drf;
    
    
    // We will make a (potentially background subtracted) spectrum - we also need to track the
    //  uncertainties in each channel
    shared_ptr<SpecUtils::Measurement> spectrum;
    vector<float> channel_counts, channel_count_uncerts;
    
    channel_counts = *foreground->gamma_counts();
    channel_count_uncerts.resize( channel_counts.size(), 0.0 );
    
    if( !background )
    {
      for( size_t i = 0; i < channel_counts.size(); ++i )
      {
        const double counts = channel_counts[i];
        channel_count_uncerts[i] = (counts <= 1.0) ? 1.0 : sqrt( counts );
      }
      
      spectrum = make_shared<SpecUtils::Measurement>( *foreground );
    }else
    {
      if( !background->energy_calibration() || !background->energy_calibration()->valid() )
        throw runtime_error( "RelActAutoCostFcn: invalid background spectrum." );
      
      if( (background->live_time() <= 0.0) || (foreground->live_time() <= 0.0) )
        throw runtime_error( "RelActAutoCostFcn: live-time missing from spectrum." );
      
      const double lt_sf = foreground->live_time() / background->live_time();
      const vector<float> &orig_back_counts = *background->gamma_counts();
      const vector<float> &back_energies = *background->energy_calibration()->channel_energies();
      const vector<float> &fore_energies = *foreground->energy_calibration()->channel_energies();
      
      vector<float> background_counts;
      SpecUtils::rebin_by_lower_edge( back_energies, orig_back_counts, fore_energies, background_counts );
      
      assert( background_counts.size() == channel_counts.size() );
      for( size_t i = 0; i < channel_counts.size(); ++i )
      {
        const double fore_counts = (channel_counts[i] < 0.0f) ? 0.0 : channel_counts[i];
        const double back_counts = (background_counts[i] < 0.0f) ? 0.0f : background_counts[i];
        const double uncert_2 = fore_counts*fore_counts + lt_sf*lt_sf*back_counts*back_counts;
        const double sub_val = fore_counts - lt_sf*back_counts;
        
        channel_counts[i] = static_cast<float>( std::max(sub_val,0.0) );
        channel_count_uncerts[i] = static_cast<float>( std::max( 1.0, sqrt(uncert_2) ) );
      }//for( loop over and set channel counts and uncertainties )
      
      spectrum = make_shared<SpecUtils::Measurement>( *foreground );
      spectrum->set_gamma_counts( make_shared<vector<float>>(channel_counts), foreground->live_time(), foreground->real_time() );
    }//if( !background ) / else
    
    solution.m_channel_counts = channel_counts;
    solution.m_channel_counts_uncerts = channel_count_uncerts;
    
    
    if( all_peaks.empty() )
      all_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( spectrum, nullptr, {}, false );
    solution.m_spectrum_peaks = all_peaks;
    
    vector<shared_ptr<const PeakDef>> peaks_in_rois;
    for( const shared_ptr<const PeakDef> &p : all_peaks )
    {
      for( const RelActCalcAuto::RoiRange &range : energy_ranges )
      {
        if( (p->mean() >= range.lower_energy) && (p->mean() <= range.upper_energy) )
        {
          peaks_in_rois.push_back( p );
          break;
        }
      }//for( const RelActCalcAuto::RoiRange &range : energy_ranges )
    }//for( const shared_ptr<const PeakDef> &p : peaks_in_roi )
    
    // Now use peaks_in_rois and manual stuff to try and estimate initial RelEff and activities
    
    
    
    
    RelActAutoCostFcn *cost_functor = nullptr;
    try
    {
      cost_functor = new RelActAutoCostFcn( options, energy_ranges, nuclides,
                                           extra_peaks, spectrum, channel_count_uncerts,
                                           input_drf, all_peaks );
    }catch( std::exception &e )
    {
      solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      solution.m_error_message = "Error initializing problem: " + string(e.what());
      
      return solution;
    }//try / catch
  
    auto cost_function = new ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn>( cost_functor );
    cost_function->SetNumResiduals( cost_functor->number_residuals() );
    
    
    solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    solution.m_drf = cost_functor->m_drf;
    solution.m_spectrum = make_shared<SpecUtils::Measurement>( *spectrum );
    
    const size_t num_pars = cost_functor->number_parameters();
    
    for( size_t i = 0; i < num_pars; ++i )
      cost_function->AddParameterBlock( 1 );
    
    vector<double> parameters( num_pars, 0.0 );
    double *pars = &parameters[0];
    
    vector<double *> parameter_blocks( num_pars );
    for( size_t i = 0; i < num_pars; ++i )
      parameter_blocks[i] = pars + i;
    
    
    ceres::Problem problem;
    
    // TODO: investigate using a LossFunction - probably really need it
    ceres::LossFunction *lossfcn = nullptr;
    problem.AddResidualBlock( cost_function, lossfcn, parameter_blocks );
    
    
    if( options.fit_energy_cal )
    {
      assert( !cost_functor->m_energy_ranges.empty() );
      
      const double lowest_energy = cost_functor->m_energy_ranges.front().lower_energy;
      const double highest_energy = cost_functor->m_energy_ranges.back().upper_energy;
      
      // We'll allow changing the gain by 1% (limit chosen fairly arbitrarily)
      problem.SetParameterLowerBound(pars + 1, 0, -0.01 );
      problem.SetParameterUpperBound(pars + 1, 0, +0.01 );
      
      if( (lowest_energy < 200) && (highest_energy > 600) )
      {
        //We'll allow changing the offset by 5 keV (limit chosen fairly arbitrarily)
        solution.m_fit_energy_cal[0] = true;
        solution.m_fit_energy_cal[1] = true;
        problem.SetParameterLowerBound(pars + 0, 0, -5.0 );
        problem.SetParameterUpperBound(pars + 0, 0, +5.0 );
      }else
      {
        //We'll only fit offset
        solution.m_fit_energy_cal[0] = true;
        problem.SetParameterBlockConstant( pars + 1 );
      }
    }else
    {
      problem.SetParameterBlockConstant( pars + 0 );
      problem.SetParameterBlockConstant( pars + 1 );
    }
    
    auto res_drf = cost_functor->m_drf;
    const size_t fwhm_start = 2;
    const size_t num_fwhm_pars = num_parameters( options.fwhm_form );
    const size_t rel_eff_start = fwhm_start + num_fwhm_pars;
    const size_t num_rel_eff_par = options.rel_eff_eqn_order + 1;
    const size_t acts_start = rel_eff_start + num_rel_eff_par;
    const size_t num_acts_par = 2*cost_functor->m_nuclides.size();
    const size_t free_peak_start = acts_start + num_acts_par;
    const size_t num_free_peak_par = 2*extra_peaks.size();
    
    assert( (free_peak_start + num_free_peak_par) == cost_functor->number_parameters() );
    
    try
    {
      if( !res_drf || !res_drf->hasResolutionInfo() )
        throw runtime_error( "No starting FWHM info was passed in, or derived from the spectrum;"
                             " please select a DRF with FWHM information." );
    
      //assert( cost_functor->m_options.fwhm_form ==
      vector<float> drfpars = res_drf->resolutionFcnCoefficients();
      bool needToFitOtherType = false;
      
      switch( cost_functor->m_options.fwhm_form )
      {
        case RelActCalcAuto::FwhmForm::Gadras:
        {
          // We are fitting to the GADRAS functional form
          assert( num_parameters( options.fwhm_form ) ==3 );
          
          switch( res_drf->resolutionFcnType() )
          {
            case DetectorPeakResponse::kGadrasResolutionFcn:
            {
              assert( drfpars.size() == 3 );
              needToFitOtherType = false;
              break;
            }//case drf FWHM is of form DetectorPeakResponse::kGadrasResolutionFcn:
              
            case DetectorPeakResponse::kSqrtPolynomial:
            {
              needToFitOtherType = true;
              break;
            }//case DetectorPeakResponse::kSqrtPolynomial:
              
            case DetectorPeakResponse::kNumResolutionFnctForm:
              assert( 0 );
              break;
          }//switch( cost_functor->m_drf->resolutionFcnType() )
          
          break;
        }//case we want GADRAS formGadras:
          
        case RelActCalcAuto::FwhmForm::Polynomial_2:
        case RelActCalcAuto::FwhmForm::Polynomial_3:
        case RelActCalcAuto::FwhmForm::Polynomial_4:
        case RelActCalcAuto::FwhmForm::Polynomial_5:
        case RelActCalcAuto::FwhmForm::Polynomial_6:
        {
          assert( num_parameters( cost_functor->m_options.fwhm_form ) == (1 + static_cast<int>(cost_functor->m_options.fwhm_form)) );
          
          switch( res_drf->resolutionFcnType() )
          {
            case DetectorPeakResponse::kGadrasResolutionFcn:
            {
              assert( drfpars.size() == 3 );
              needToFitOtherType = true;
              break;
            }//kGadrasResolutionFcn
              
            case DetectorPeakResponse::kSqrtPolynomial:
            {
              needToFitOtherType = false;
              break;
            }
              
            case DetectorPeakResponse::kNumResolutionFnctForm:
              assert( 0 );
              break;
          }//switch( cost_functor->m_drf->resolutionFcnType() )
          
          break;
        }//case: we want polynomial form
      }//switch( m_options.fwhm_form )
    
      if( needToFitOtherType )
      {
        // Make a vector of ~10 equally spaced peaks, with uncert 10% that peaks width -
        //  fairly arbitrary
        const double lowest_energy = cost_functor->m_energy_ranges.front().lower_energy;
        const double highest_energy = cost_functor->m_energy_ranges.back().upper_energy;
        const double delta_energy = 0.1*(highest_energy - lowest_energy);
        
        auto fake_peaks = make_shared<std::deque< std::shared_ptr<const PeakDef> > >();
        for( double ene = lowest_energy; ene <=(1.001*highest_energy); ene += delta_energy )
        {
          const auto sigma = res_drf->peakResolutionSigma(ene);
          auto p = make_shared<PeakDef>( ene, sigma, 1000.0 );
          p->setSigmaUncert( 0.1*sigma );
          fake_peaks->push_back( p );
        }
      
        const bool haveGadras = (res_drf->resolutionFcnType() == DetectorPeakResponse::kGadrasResolutionFcn);
        const bool wantGadras = (cost_functor->m_options.fwhm_form == RelActCalcAuto::FwhmForm::Gadras);
        assert( haveGadras != wantGadras );
        assert( !haveGadras || (drfpars.size() == 3) );
        
        const auto formToFit = wantGadras
                               ? DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn
                               : DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
                                
        
        vector<float> new_sigma_coefs, sigma_coef_uncerts;
        MakeDrfFit::performResolutionFit( fake_peaks, formToFit,
                                       highres, num_fwhm_pars, new_sigma_coefs, sigma_coef_uncerts );
      
        assert( new_sigma_coefs.size() == num_fwhm_pars );
        assert( haveGadras || (new_sigma_coefs.size() == 3) );
        assert( new_sigma_coefs.size() == num_fwhm_pars );
        
        drfpars = new_sigma_coefs;
      }//if( needToFitOtherType )
      
      if( drfpars.size() != num_fwhm_pars )
        throw logic_error( "Unexpected num parameters from fit to FWHM function (logic error)." );
      
      for( size_t i = 0; i < drfpars.size(); ++i )
        parameters[fwhm_start + i] = drfpars[i];
    }catch( std::exception &e )
    {
      solution.m_warnings.push_back( "Failed to create initial FWHM estimation, but will continue anyway: " + string(e.what()) );
     
      if( cost_functor->m_options.fwhm_form == RelActCalcAuto::FwhmForm::Gadras )
      {
        // These values are from an arbitrary HPGe and NaI hand-held detector
        if( highres )
        {
          parameters[fwhm_start + 0] = 1.54;
          parameters[fwhm_start + 1] = 0.264;
          parameters[fwhm_start + 2] = 0.33;
        }else
        {
          parameters[fwhm_start + 0] = -6.5;
          parameters[fwhm_start + 1] = 7.5;
          parameters[fwhm_start + 2] = 0.55;
        }
      }else
      {
        // The below coefficient values are from a simple fit to the GADRAS coefficients above
        vector<double> fwhm_pars( num_fwhm_pars, 0.0 );
        if( highres )
        {
          if( num_fwhm_pars >= 1 )
            parameters[fwhm_start + 0] = 2.4083;
          if( num_fwhm_pars >= 2 )
            parameters[fwhm_start + 1] = 0.00113706;
          if( num_fwhm_pars >= 3 )
            parameters[fwhm_start + 2] = 3.36758e-07;
        }else
        {
          // Valid above ~20 keV, which low res would never go down this far anyway
          if( num_fwhm_pars >= 1 )
            parameters[fwhm_start + 0] = -19.7642;
          if( num_fwhm_pars >= 2 )
            parameters[fwhm_start + 1] = 1.73526;
          if( num_fwhm_pars >= 3 )
            parameters[fwhm_start + 2] = 0.00283492;
          if( num_fwhm_pars >= 4 )
            parameters[fwhm_start + 3] = -7.18936e-07;
        }//if( highres ) / else
      }//if( we want GADRAS FWHM function ) / else
    }//if( has DRF resolution info ) / else
    
    assert( options.rel_eff_eqn_order != 0 );
    if( options.rel_eff_eqn_order == 0 )
      throw runtime_error( "Relative efficiency order must be at least 1 (but should probably be at least 2)" );
    
    // Note that \p rel_eff_order is the number of energy-dependent terms we will fit for in the
    //  relative efficiency equation (i.e., we will also have 1 non-energy dependent term, so we
    //  will fit for (1 + rel_eff_order) parameters).
 //
    //  We'll start with all values for rel eff equation at zero, except maybe the first one
    //  (e.g., we'll start with rel eff line == 1.0 for all energies), and then after we estimate
    //  starting activities, we'll do a little better job estimating things.
    for( size_t rel_eff_index = 0; rel_eff_index <= options.rel_eff_eqn_order; ++rel_eff_index )
    {
      parameters[rel_eff_start + rel_eff_index] = 0.0;
    }
    
    switch( options.rel_eff_eqn_type )
    {
      case RelActCalc::RelEffEqnForm::LnX:
      case RelActCalc::RelEffEqnForm::LnY:
        parameters[rel_eff_start + 0] = 1.0;
        break;
      case RelActCalc::RelEffEqnForm::LnXLnY:
      case RelActCalc::RelEffEqnForm::FramEmpirical:
        break;
    }//switch( eqn_form )
    
    
    
    
    // Estimate the initial rel-eff equation
    //  Some things to consider are that we have a really poor estimate of activities right now
    //  and that we want the Rel Eff curve to be 1.0 at the lower energy of the lowest ROI...
    //  Maybe the thing to do is to use the auto-fit peaks, implement the manual fitting back-end
    //   and assigning of gammas to peaks, then use that to first solve for relative activity and
    //   relative efficiency
    bool succesfully_estimated_re_and_ra = false;
    try
    {
      const double base_rel_eff_uncert = 1.0;
      
      // Filter peaks to those in ranges we want
      vector<RelActCalcManual::GenericPeakInfo> peaks_in_range;
      vector<std::shared_ptr<const PeakDef> > debug_manual_display_peaks;
      for( const shared_ptr<const PeakDef> &p : all_peaks )
      {
        bool use_peak = energy_ranges.empty();
        for( const auto &r : energy_ranges )
        {
          if( (p->mean() >= r.lower_energy) && (p->mean() <= r.upper_energy) )
            use_peak = true;
        }
        
        if( use_peak )
        {
          RelActCalcManual::GenericPeakInfo peak;
          peak.m_energy = p->mean();
          peak.m_fwhm = p->gausPeak() ? p->fwhm() : 2.634*0.25*p->roiWidth();
          peak.m_counts = p->amplitude();
          peak.m_counts_uncert = p->amplitudeUncert();
          peak.m_base_rel_eff_uncert = 0.5;
          peaks_in_range.push_back( peak );
          
          debug_manual_display_peaks.push_back( p );
        }//if( use_peak )
      }//for( const shared_ptr<const PeakDef> &p : all_peaks )
      
      
      vector<RelActCalcManual::SandiaDecayNuc> nuc_sources;
      for( const RelActCalcAuto::NucInputInfo &info : nuclides )
      {
        RelActCalcManual::SandiaDecayNuc nucinfo;
        nucinfo.nuclide = info.nuclide;
        nucinfo.age = info.age;
        nuc_sources.push_back( nucinfo );
      }
      
      const double cluster_sigma = 1.5;
      const auto peaks_with_nucs = add_nuclides_to_peaks( peaks_in_range, nuc_sources, cluster_sigma );
      
      vector<RelActCalcManual::GenericPeakInfo> peaks_with_sources;
      for( const auto &p : peaks_with_nucs )
      {
        bool is_floater_peak = false;
        // TODO: - Need to deal with extra floater peaks in getting the initial "manual" solution
        //         Right now were just removing the "fit" peak if its within 1 sigma of a floating
        //         peak (since if the peak is floating, it wont add any info to the manual solution)
        //         but maybe there is a better way of dealing with things?
        for( const RelActCalcAuto::FloatingPeak &floater : extra_peaks )
        {
          if( fabs(floater.energy - p.m_energy) < 1.0*(p.m_fwhm/2.634) )
            is_floater_peak = true;
        }//for( const RelActCalcAuto::FloatingPeak &floater : extra_peaks )
        
        
        if( !is_floater_peak && !p.m_source_gammas.empty() )
          peaks_with_sources.push_back( p );
      }//for( const auto &p : peaks_with_nucs )
      
      
      
      RelActCalcManual::RelEffSolution manual_solution
                 = RelActCalcManual::solve_relative_efficiency( peaks_with_sources,
                                              options.rel_eff_eqn_type, options.rel_eff_eqn_order );
      
      cout << "Initial estimates:" << endl;
      manual_solution.print_summary( cout );
      
      ofstream debug_manual_html( "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/initial_manual_estimate.html" );
      manual_solution.print_html_report( debug_manual_html, options.spectrum_title, spectrum, debug_manual_display_peaks );
      
      //Need to fill out rel eff starting values and rel activities starting values
      
      if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
        throw runtime_error( manual_solution.m_error_message );
      
      assert( manual_solution.m_rel_eff_eqn_coefficients.size() == (options.rel_eff_eqn_order + 1) );
      
      const string rel_eff_eqn_str = RelActCalc::rel_eff_eqn_text( manual_solution.m_rel_eff_eqn_form, manual_solution.m_rel_eff_eqn_coefficients );
      cout << "Starting with initial rel. eff. eqn = " << rel_eff_eqn_str << endl;
      
      for( size_t i = 0; i <= options.rel_eff_eqn_order; ++i )
        parameters[rel_eff_start + i] = manual_solution.m_rel_eff_eqn_coefficients[i];
            
      // Manual "relative_activity" assumes a measurement of 1-second (or rather peaks are in CPS),
      //  but this "auto" relative activity takes into account live_time
      const double live_time = spectrum->live_time();
      
      for( size_t nuc_num = 0; nuc_num < nuclides.size(); ++nuc_num )
      {
        const RelActCalcAuto::NucInputInfo &nuc = nuclides[nuc_num];
        const double rel_act = manual_solution.relative_activity( nuc.nuclide->symbol ) / live_time;
        
        const size_t act_index = acts_start + 2*nuc_num;
        cout << "Updating initial activity estimate for " << nuc.nuclide->symbol << " from "
        << parameters[act_index] << " to " << rel_act << endl;
        
        parameters[act_index] = rel_act;
      }
      
      succesfully_estimated_re_and_ra = true;
    }catch( std::exception &e )
    {
      cerr << "Failed to do initial estimate ov RelEff curve: " << e.what() << endl;
      
      solution.m_warnings.push_back( "Initial estimate of relative efficiency curve failed ('"
                                     + string(e.what())
                                     + "'), using a flat line as a starting point" );
      
 
    }//try / catch
    
    
/*
    cerr << "\n\n\n\nSetting RelEff eqn to canned initial values!\n\n\n";
#warning "Setting RelEff eqn to preset values for debugging!"
    assert( options.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramEmpirical );
    assert( options.rel_eff_eqn_order == 3 );
    
    parameters[rel_eff_start + 0] = -39.549338;
    parameters[rel_eff_start + 1] = 5208.317658;
    parameters[rel_eff_start + 2] = 12.714142;
    parameters[rel_eff_start + 3] = -0.970901;
*/
    
    
    for( size_t nuc_num = 0; nuc_num < nuclides.size(); ++nuc_num )
    {
      const RelActCalcAuto::NucInputInfo &nuc = nuclides[nuc_num];
      const size_t act_index = acts_start + 2*nuc_num;
      const size_t age_index = act_index + 1;
      
      assert( nuc.nuclide );
      assert( nuc.age >= 0.0 );
      if( !nuc.nuclide || nuc.nuclide->isStable() )
        throw runtime_error( "Invalid nuclide" );
      
      if( nuc.age < 0.0 )
        throw runtime_error( "Invalid nuclide (" + nuc.nuclide->symbol + ") age" );
      
      parameters[age_index] = nuc.age;
      bool fit_age = nuc.fit_age;
      
      // If ages of nuclides of an element are set to all be the same, we will use the age slot
      //  of the first nuclide to control the age of all nuclides for that element; we'll check
      //  for this here, and if we have a nuclide of the same element preceding \c nuc_num, we'll
      //  fix the age here
      if( (nuc_num > 0) && options.nucs_of_el_same_age )
      {
        bool found_control_age = false;
        for( size_t i = 0; i < nuc_num; ++i )
        {
          if( nuclides[i].nuclide->atomicNumber == nuc.nuclide->atomicNumber )
          {
            if( (nuclides[i].age != nuc.age) || (nuclides[i].fit_age != nuc.fit_age) )
              throw runtime_error( "When its specified that all nuclides of the element"
                                   " must have the same age, and same (initial) age value and"
                                   " wether or not to be fit must be specified." );
            
            fit_age = false;
            parameters[age_index] = -1.0; //so we will catch logic error as an assert later.
            
            break;
          }//if( we have found the controlling nuclide age )
        }//for( loop over previously seen nuclides )
      }//if( (nuc_num > 0) && options.nucs_of_el_same_age )
      
      
      if( fit_age )
      {
        // TODO: maybe figure out a better method of setting upper age bound
        
        // If we are tying the ages of all nuclides in an element together, we dont want to
        //  accidentally limit the max age based on the shortest nuclides half-life (think U-237 for
        //  uranium).
        double half_life = nuc.nuclide->halfLife;
        if( options.nucs_of_el_same_age )
        {
          for( size_t i = 0; i < nuclides.size(); ++i )
          {
            if( nuclides[i].nuclide->atomicNumber == nuc.nuclide->atomicNumber )
              half_life = std::max(half_life, nuclides[i].nuclide->halfLife);
          }
        }//if( options.nucs_of_el_same_age )
        
        double max_age = std::max( 5.0*nuc.age, 15.0*half_life );
        
        // We'll clamp the upper age, to the range where humans may have done the seperation
        // TODO: is there any problem with a nuclide >100 years, that isnt effectively infinite?
        max_age = std::min( max_age, 120*PhysicalUnits::year );
        
        problem.SetParameterUpperBound(pars + age_index, 0, max_age );
        problem.SetParameterLowerBound(pars + age_index, 0, 0.0 );
      }else
      {
        problem.SetParameterBlockConstant( pars + age_index );
      }
      
      if( !succesfully_estimated_re_and_ra )
      {
        // Now do a very rough initial activity estimate; there are two obvious ideas
        //  1) Either take in user peaks, or auto-search peaks, and try to match up.  However, we may
        //     run into issues that nuc.nuclide may not have any fit-peaks (in which case I guess we
        //     could start with an activity of 0 bq).  We could use
        //     #RelActCalcManual::fit_act_to_rel_eff to fit the initial activities.
        //  2) Take the gammas 10% of the highest yield in each energy range, assume they make peaks,
        //     and estimate the peak area, neglecting interferences and continuum, then take the
        //     median estimated relative activity.  This is all fairly arbitrary.
        //
        //  To avoid fitting peaks, and hassles that come along, lets try (2) out, and also, we will
        //  only do the simplest of peak area estimations
        vector<tuple<double,double,double>> top_energy_to_rel_act; //{energy, br, rel. act.}
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc.nuclide, ns_decay_act_mult, nuc.age );
        const vector<SandiaDecay::EnergyRatePair> gammas = mix.photons(0.0);
        for( const auto &erange : energy_ranges )
        {
          vector<tuple<double,double,double>> energy_to_rel_act; //{energy, br, rel. act.}
          for( const SandiaDecay::EnergyRatePair &er : gammas )
          {
            if( er.numPerSecond <= std::numeric_limits<float>::min() ) //1.17549e-38
              continue;
            
            if( (er.energy < erange.lower_energy) || (er.energy > erange.upper_energy) )
              continue;
            
            const double det_fwhm = cost_functor->fwhm( er.energy, parameters );
            double data_count = foreground->gamma_integral(er.energy - 0.5*det_fwhm, er.energy + 0.5*det_fwhm);
            if( background )
            {
              const double backarea = background->gamma_integral(er.energy - 0.5*det_fwhm, er.energy + 0.5*det_fwhm);
              data_count -= backarea * foreground->live_time() / background->live_time();
              data_count = std::max( 0.0, data_count );
            }//
            
            const double yield = er.numPerSecond / ns_decay_act_mult;
            
            // The FWHM area covers about 76% of gaussian area
            const double rel_act = data_count / yield / foreground->live_time() / 0.76;
            
            energy_to_rel_act.emplace_back( er.energy, yield, rel_act );
          }//for( const SandiaDecay::EnergyRatePair &er : gammas )
          
          std::sort( begin(energy_to_rel_act), end(energy_to_rel_act),
                    []( const auto &lhs, const auto &rhs) -> bool{
            return get<1>(lhs) > get<1>(rhs);
          } );
          
          for( const auto &i : energy_to_rel_act )
          {
            if( get<1>(i) >= 0.1*get<1>(energy_to_rel_act[0]) )
              top_energy_to_rel_act.push_back( i );
          }//for( const auto &l : energy_to_rel_act )
        }//for( const auto &erange : energy_ranges )
        
        std::sort( begin(top_energy_to_rel_act), end(top_energy_to_rel_act),
                  []( const auto &lhs, const auto &rhs) -> bool{
          return get<1>(lhs) > get<1>(rhs);
        } );
        
        if( !top_energy_to_rel_act.empty() )
        {
          parameters[act_index] = std::get<2>( top_energy_to_rel_act[top_energy_to_rel_act.size()/2] );
          cout << "Setting initial relative activity for " << nuc.nuclide->symbol << " to " << parameters[act_index] << " bq" << endl;
        }else
        {
          parameters[act_index] = 100*PhysicalUnits::bq;
          solution.m_warnings.push_back( "Could not estimate a starting activity for "
                                        + nuc.nuclide->symbol + ".  This may be because there are no"
                                        " significant gammas for the nuclide in the selected energy"
                                        " ranges." );
        }
      }//if( !succesfully_estimated_re_and_ra )
      
      problem.SetParameterLowerBound(pars + act_index, 0, 0.0 );
    }//for( size_t nuc_num = 0; nuc_num < nuclides.size(); ++nuc_num )
    
    
    // Floating peaks; one parameter for amplitude, one for FWHM (which will usually be unused)
    for( size_t extra_peak_index = 0; extra_peak_index < extra_peaks.size(); ++extra_peak_index )
    {
      const RelActCalcAuto::FloatingPeak &peak = extra_peaks[extra_peak_index];
      assert( peak.energy > 0.0 );
      if( peak.energy <= 0.0 )
        throw runtime_error( "Invalid floating peak energy." );
      
      const size_t amp_index = free_peak_start + 2*extra_peak_index;
      const size_t fwhm_index = amp_index + 1;
      
      const double det_fwhm = cost_functor->fwhm( peak.energy, parameters );
      double data_count = foreground->gamma_integral(peak.energy - 0.5*det_fwhm, peak.energy + 0.5*det_fwhm);
      if( background )
      {
        const double backarea = background->gamma_integral(peak.energy - 0.5*det_fwhm, peak.energy + 0.5*det_fwhm);
        data_count -= backarea * foreground->live_time() / background->live_time();
        data_count = std::max( 0.0, data_count );
      }//
      
      // Assume data is entirely due to peak; mean +- 0.5*fwhm is ~76% gaussian area
      parameters[amp_index] = data_count / 0.76;
      problem.SetParameterLowerBound(pars + amp_index, 0, 0.0 );
      
      parameters[fwhm_index] = 1.0;
      
      if( peak.release_fwhm )
      {
        // We'll say the width of this peak can not be less than 25% of what is expected
        //  TODO: we probably also want to put a lower bound of 3 or 4 channels as well
        problem.SetParameterLowerBound(pars + fwhm_index, 0, 0.25 );
        
        // TODO: set an upper bound on peak width - this needs to be a convolution of Q value and detector resolution - maybe an input to this function
        //  Right now we'll just say a factor of 4 - which is plenty, unless its caused by a large Q reaction.
        problem.SetParameterUpperBound(pars + fwhm_index, 0, 4.0 );
      }else
      {
        problem.SetParameterBlockConstant(pars + fwhm_index);
      }
    }//for( size_t extra_peak_index = 0; extra_peak_index < extra_peaks.size(); ++extra_peak_index )
    
    
    // Okay - we've set our problem up
    ceres::Solver::Options ceres_options;
    ceres_options.linear_solver_type = ceres::DENSE_QR;
    ceres_options.minimizer_progress_to_stdout = true; //true;
    
    // TODO: there are a ton of ceres::Solver::Options that might be useful for us to set
    
    /*
     // If we to implement cancelling a calculation, we could do something like:
    class CheckTerminateCallback : public ceres::IterationCallback
    {
      std::function<bool()> m_check_terminate;
    public:
      CheckTerminateCallback( std::function<bool()> fcn ) : m_check_terminate(fcn){};
      
      virtual CallbackReturnType operator()(const IterationSummary& summary)
      {
        return (!m_check_terminate || !m_check_terminate())
               ? ceres::CallbackReturnType::SOLVER_CONTINUE
               : ceres::CallbackReturnType::SOLVER_ABORT; //CallbackReturnType::SOLVER_TERMINATE_SUCCESSFULLY
      }
    };//CheckTerminateCallback
     
     CheckTerminateCallback callback( [](){ return false; } );
     ceres_options.callbacks.push_back( &callback );
     */
    
    // Setting ceres_options.num_threads >1 doesnt seem to do much (any?) good - so instead we do
    //  some multiple threaded computations in RelActAutoSolution::eval(...)
    ceres_options.num_threads = std::thread::hardware_concurrency();
    if( !ceres_options.num_threads )
    {
      assert( 0 );
      ceres_options.num_threads = 4;
    }
    
    
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_options, &problem, &summary);
    //std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
    
    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
        break;
        
      case ceres::NO_CONVERGENCE:
      case ceres::FAILURE:
      case ceres::USER_FAILURE:
        solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
        solution.m_error_message += "The L-M solving failed.";
        break;
    }//switch( summary.termination_type )
    
    solution.m_num_function_eval_solution = cost_functor->m_ncalls;
    
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //
    cov_options.num_threads = ceres_options.num_threads;
    
    vector<double> uncertainties( num_pars, 0.0 ), uncerts_squared( num_pars, 0.0 );
    
    ceres::Covariance covariance(cov_options);
    vector<pair<const double*, const double*> > covariance_blocks;
    for( size_t i = 0; i < num_pars; ++i )
      covariance_blocks.push_back( make_pair( pars + i, pars + i ) );
    
    auto add_cov_block = [&covariance_blocks, &pars]( size_t start, size_t num ){
      for( size_t i = start; i < num; ++i )
        for( size_t j = start; j < i; ++j )
          covariance_blocks.push_back( make_pair( pars + i, pars + j ) );
    };
    
    add_cov_block( rel_eff_start, num_rel_eff_par );
    add_cov_block( acts_start, num_acts_par );
    add_cov_block( fwhm_start, num_fwhm_pars );
    
    if( !covariance.Compute(covariance_blocks, &problem) )
    {
      cerr << "Failed to compute final covariances!" << endl;
      solution.m_warnings.push_back( "Failed to compute final covariances." );
    }else
    {
      for( size_t i = 0; i < num_pars; ++i )
      {
        covariance.GetCovarianceBlock( pars + i, pars + i, &(uncerts_squared[i]) );
        if( uncerts_squared[i] >= 0.0 )
          uncertainties[i] = sqrt( uncerts_squared[i] );
        else
          uncertainties[i] = std::numeric_limits<double>::quiet_NaN();
      }
      
      // TODO: check the covariance is actually defined right (like not swapping row/col, calling the right places, etc).
      auto get_cov_block = [&pars,&covariance]( const size_t start, const size_t num, vector<vector<double>> &cov ){
        cov.clear();
        cov.resize( num, vector<double>(num,0.0) );
        
        for( size_t i = start; i < num; ++i )
        {
          for( size_t j = start; j < i; ++j )
          {
            covariance.GetCovarianceBlock( pars + i, pars + j, &(cov[i][j]) );
            cov[j][i] = cov[i][j];
          }
        }
      };
      
      get_cov_block( rel_eff_start, num_rel_eff_par, solution.m_rel_eff_covariance );
      get_cov_block( acts_start, num_acts_par, solution.m_rel_act_covariance );
      get_cov_block( fwhm_start, num_fwhm_pars, solution.m_fwhm_covariance );
    }//if( we failed to get covariance ) / else
    
    
    solution.m_num_function_eval_total = cost_functor->m_ncalls;
  
    solution.m_final_parameters = parameters;
    
    solution.m_energy_cal_adjustments[0] = parameters[0];
    solution.m_energy_cal_adjustments[1] = parameters[1];
    
    shared_ptr<const SpecUtils::EnergyCalibration> new_cal = cost_functor->m_energy_cal;
    if( options.fit_energy_cal )
    {
      auto new_cal = make_shared<SpecUtils::EnergyCalibration>( *cost_functor->m_energy_cal );
      
      const size_t num_channel = cost_functor->m_energy_cal->num_channels();
      const SpecUtils::EnergyCalType energy_cal_type = cost_functor->m_energy_cal->type();
      
      switch( energy_cal_type )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::FullRangeFraction:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        {
          vector<float> coefs = cost_functor->m_energy_cal->coefficients();
          assert( coefs.size() >= 2 );
          coefs[0] += parameters[0];
          coefs[1] *= (1.0 + parameters[1]);
          
          const auto &dev_pairs = cost_functor->m_energy_cal->deviation_pairs();
          
          if( energy_cal_type == SpecUtils::EnergyCalType::FullRangeFraction )
            new_cal->set_full_range_fraction( num_channel, coefs, dev_pairs );
          else
            new_cal->set_polynomial( num_channel, coefs, dev_pairs );
          
          break;
        }//case polynomial or FRF
          
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        {
          vector<float> lower_energies = *cost_functor->m_energy_cal->channel_energies();
          for( float &energy : lower_energies )
            energy = parameters[0] + ((1.0 + parameters[1]) * energy);
          
          new_cal->set_lower_channel_energy( num_channel, std::move(lower_energies) );
          
          break;
        }//case LowerChannelEdge
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
          break;
      }//switch( m_energy_cal->type() )
      
      new_cal = new_cal;
      solution.m_spectrum->set_energy_calibration( new_cal );
    }//if( options.fit_energy_cal )
    
  
    vector<PeakDef> fit_peaks;
    for( const auto &range : cost_functor->m_energy_ranges )
    {
      PeaksForEnergyRange these_peaks = cost_functor->peaks_for_energy_range( range, parameters );
      fit_peaks.insert( end(fit_peaks), begin(these_peaks.peaks), end(these_peaks.peaks) );
    }
    
    std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMean );
    
    // \c fit_peaks are in the original energy calibration of the spectrum, we may need to adjust
    //  them to match the new energy calibration
    if( new_cal != cost_functor->m_energy_cal )
    {
      deque<shared_ptr<const PeakDef>> tmp_peaks;
      for( const auto &p : fit_peaks )
        tmp_peaks.push_back( make_shared<const PeakDef>( p ) );
      
      auto adjusted_peaks = EnergyCal::translatePeaksForCalibrationChange( tmp_peaks,
                                       cost_functor->m_energy_cal, new_cal );
      
      fit_peaks.clear();
      for( const auto &p : adjusted_peaks )
        fit_peaks.push_back( *p );
    }//if( new_cal != cost_functor->m_energy_cal )
    
    solution.m_fit_peaks = fit_peaks;
    
    
    const auto rel_eff_iter = begin(parameters) + rel_eff_start;
    solution.m_rel_eff_coefficients.clear();
    solution.m_rel_eff_coefficients.insert( end(solution.m_rel_eff_coefficients),
                                      rel_eff_iter, rel_eff_iter + options.rel_eff_eqn_order + 1 );
    
    const auto rel_act_iter = begin(parameters) + acts_start;
    for( size_t act_index = 0; act_index < cost_functor->m_nuclides.size(); ++ act_index )
    {
      const NucInputGamma &nuc_input = cost_functor->m_nuclides[act_index];
      
      const size_t par_nuc_index = cost_functor->nuclide_index( nuc_input.nuclide );
      const size_t par_act_index = acts_start + 2*par_nuc_index;
      const size_t par_age_index = par_act_index + 1;
      
      RelActCalcAuto::NuclideRelAct nuc_output;
      nuc_output.nuclide = nuc_input.nuclide;
      nuc_output.age = cost_functor->age( nuc_input.nuclide, parameters );
      nuc_output.age_was_fit = nuc_input.fit_age;
      nuc_output.rel_activity = cost_functor->relative_activity( nuc_input.nuclide, parameters );
      
      nuc_output.age_uncertainty = uncertainties[par_age_index];
      nuc_output.rel_activity_uncertainty = uncertainties[par_act_index];
      
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc_input.nuclide, ns_decay_act_mult, nuc_output.age );
      const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
      nuc_output.gamma_energy_br.clear();
      for( const auto &gamma : gammas )
      {
        const double yield = gamma.numPerSecond / ns_decay_act_mult;
        nuc_output.gamma_energy_br.push_back( {gamma.energy, yield} );
      }
      
      if( IsNan(nuc_output.age_uncertainty) )
      {
        solution.m_warnings.push_back( "Variance for nuclide " + nuc_input.nuclide->symbol
                                      + " age, is invalid ("
                                      + std::to_string(uncerts_squared[par_age_index]) + ")" );
      }
      
      if( IsNan(nuc_output.rel_activity_uncertainty) )
      {
        solution.m_warnings.push_back( "Variance for nuclide " + nuc_input.nuclide->symbol
                                      + " activity, is invalid ("
                                      + std::to_string(uncerts_squared[par_act_index]) + ")" );
      }

      solution.m_rel_activities.push_back( nuc_output );
    }//for( size_t act_index = 0; act_index < cost_functor->m_nuclides.size(); ++ act_index )
    
    solution.m_fwhm_form = options.fwhm_form;
    solution.m_fwhm_coefficients.clear();
    for( size_t i = 0; i < num_fwhm_pars; ++i )
      solution.m_fwhm_coefficients.push_back( parameters[fwhm_start + i] );
    
    for( size_t i = 0; i < cost_functor->m_extra_peaks.size(); ++i )
    {
      const size_t amp_index = free_peak_start + 2*i + 0;
      const size_t fwhm_index = free_peak_start + 2*i + 1;
      
      RelActCalcAuto::FloatingPeakResult peak;
      peak.energy = cost_functor->m_extra_peaks[i].energy;
      peak.amplitude = parameters[amp_index];
      peak.fwhm = parameters[fwhm_index];
      
      if( !cost_functor->m_extra_peaks[i].release_fwhm )
      {
        const double true_energy = cost_functor->un_apply_energy_cal_adjustment( peak.energy, parameters );
        peak.fwhm = cost_functor->fwhm( true_energy, parameters );
        
        // TODO: implement evaluating uncertainty of FWHM, given covariance.
        peak.fwhm_uncert = -1;
      }else
      {
        peak.fwhm_uncert = uncertainties[fwhm_index];
        if( IsNan(peak.fwhm_uncert) )
        {
          solution.m_warnings.push_back( "Variance for floating peak at " + std::to_string(peak.energy)
                                        + " FWHM is invalid ("
                                        + std::to_string(uncerts_squared[fwhm_index]) + ")" );
        }
      }//if( we didnt let FWHM float ) / else
      
      peak.amplitude_uncert = uncertainties[amp_index];
      
      
      if( IsNan(peak.amplitude_uncert) )
      {
        solution.m_warnings.push_back( "Variance for floating peak at " + std::to_string(peak.energy)
                                      + " amplitude is invalid ("
                                      + std::to_string(uncerts_squared[amp_index]) + ")" );
      }
      
      solution.m_floating_peaks.push_back( peak );
    }//for( size_t i = 0; i < cost_functor->m_extra_peaks.size(); ++i )
    
    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    cost_functor->eval( parameters, residuals.data() );
    solution.m_chi2 = 0.0;
    for( const double v : residuals )
      solution.m_chi2 += v*v;
    
    // TODO: need to setup the DOF
    solution.m_dof = -1.0;
    solution.m_warnings.push_back( "Not currently calculating DOF - need to implement" );
    //solution.m_dof = residuals.size() - ;
    
    solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
    
    return solution;
  }//RelActCalcAuto::RelActAutoSolution solve_ceres( ceres::Problem &problem )
  
  
  float fwhm( const float energy, const std::vector<double> &x ) const
  {
    const auto drf_start = begin(x) + 2;
    const size_t num_drf_par = num_parameters(m_options.fwhm_form);
    
    const vector<float> drfx( drf_start, drf_start + num_drf_par );
    
    const auto fctntype = (m_options.fwhm_form == RelActCalcAuto::FwhmForm::Gadras)
                          ? DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn
                          : DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
    
    return DetectorPeakResponse::peakResolutionFWHM( energy, fctntype, drfx );
  }//float fwhm(...)
  
  
  size_t nuclide_index( const SandiaDecay::Nuclide * const nuc ) const
  {
    assert( nuc );
    if( !nuc )
      throw runtime_error( "Null nuclide" );
    
    size_t nuc_index = m_nuclides.size();
    for( size_t i = 0; i < nuc_index; ++i )
    {
      if( m_nuclides[i].nuclide == nuc )
        nuc_index = i;
    }
    
    assert( nuc_index < m_nuclides.size() );
    if( nuc_index >= m_nuclides.size() )
      throw std::logic_error( "Invalid nuclide." );
    
    return nuc_index;
  }
  
  double relative_activity( const SandiaDecay::Nuclide * const nuc, const std::vector<double> &x ) const
  {
    const size_t nuc_index = nuclide_index( nuc );
    const size_t act_start_index = 2  //energy adjustments
                                   + num_parameters(m_options.fwhm_form)
                                   + m_options.rel_eff_eqn_order + 1;
    assert( (act_start_index + 2*m_nuclides.size() + 2*m_extra_peaks.size()) == number_parameters() );
    
    assert( x.size() == number_parameters() );
    return x[act_start_index + 2*nuc_index];
  }//double relative_activity(...)
  
  /** Returns the nuclide that is responsible for setting the passed in nuclides age.
   Will return the input nuclide if it controls its own age.
   */
  const SandiaDecay::Nuclide *age_controlling_nuc( const SandiaDecay::Nuclide * const nuc ) const
  {
    if( !m_options.nucs_of_el_same_age )
      return nuc;
    
    const size_t nuc_index = nuclide_index( nuc );
    
    for( size_t i = 0; i < nuc_index; ++i )
    {
      if( m_nuclides[i].nuclide->atomicNumber == m_nuclides[nuc_index].nuclide->atomicNumber )
        return m_nuclides[i].nuclide;
    }
    
    return nuc;
  }//age_controlling_nuc(...)
  
  
  double age( const SandiaDecay::Nuclide * const nuc, const std::vector<double> &x ) const
  {
    const size_t act_start_index = 2  //energy adjustments
                                   + num_parameters(m_options.fwhm_form)
                                   + m_options.rel_eff_eqn_order + 1;
    
    // The `parent_nuc` will often times be `nuc`.
    const SandiaDecay::Nuclide *parent_nuc = age_controlling_nuc( nuc );
    const size_t parent_nuc_index = nuclide_index( parent_nuc );
    
    const double age = x[act_start_index + 2*parent_nuc_index + 1];
    
    // We'll allow a little bit of slack on the age...
    assert( age >= -std::numeric_limits<float>::epsilon() );
    if( age < -std::numeric_limits<float>::epsilon() )
      throw runtime_error( "Negative age for " + m_nuclides[parent_nuc_index].nuclide->symbol + " found." );
    
    return std::max( age, 0.0 );
  }//age( nuclide )
  
  
  bool is_fixed_age( const SandiaDecay::Nuclide * const nuc ) const
  {
    const SandiaDecay::Nuclide *parent_nuc = age_controlling_nuc( nuc );
    const size_t nuc_index = nuclide_index( parent_nuc );
    return m_nuclides[nuc_index].fit_age;
  }//
  
  
  double relative_eff( const double energy, const std::vector<double> &x ) const
  {
    const size_t rel_eff_start_index = 2 + num_parameters(m_options.fwhm_form);
    assert( (rel_eff_start_index + m_options.rel_eff_eqn_order + 1) < x.size() );
    
    return RelActCalc::eval_eqn( energy, m_options.rel_eff_eqn_type,
                                 &(x[rel_eff_start_index]), m_options.rel_eff_eqn_order + 1 );
  }//
  
  
  /** Translates from a "true" energy (e.g., that of a gamma, or ROI bounds), to the energy
   of m_energy_cal.
   */
  double apply_energy_cal_adjustment( double energy, const std::vector<double> &x ) const
  {
    assert( x.size() > 2 );
    
    if( !m_options.fit_energy_cal )
    {
      assert( fabs(x[0]) < std::numeric_limits<float>::epsilon() );
      assert( fabs(x[1]) < std::numeric_limits<float>::epsilon() );
      
      return energy;
    }//if( we arent fitting energy cal )
    
    // Check adjustments are near the limits we placed (which was [-5,5], and [-0.01,0.01])
    assert( fabs(x[0]) <= 5.1 );
    assert( fabs(x[1]) <= 0.011 );
      
    // TODO: implement and try out optimization so that if there is no adjustment to be made, skip doing the energy cal adjustment.
    //       Note: we do need to be careful we dont make the cuts so large that the steps of the
    //       numerical differentiation around zero will fail (I did have trouble with this for
    //       relative activities and using FLT_EPSILON as a cutoff).
    //if( (fabs(x[0]) < 1.0E-9 ) && (fabs(x[1]) < 1.0E-12) )
    //  return energy;
    
    switch( m_energy_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      {
        vector<float> coefs = m_energy_cal->coefficients();
        assert( coefs.size() >= 2 );
        coefs[0] += x[0];
        coefs[1] *= (1.0 + x[1]);
        
        const auto &dev_pairs = m_energy_cal->deviation_pairs();
        const double channel = m_energy_cal->channel_for_energy( energy );
        
        if( m_energy_cal->type() == SpecUtils::EnergyCalType::FullRangeFraction )
          return SpecUtils::fullrangefraction_energy( channel, coefs, m_energy_cal->num_channels(), dev_pairs );
        
        return SpecUtils::polynomial_energy( channel, coefs, dev_pairs );
        break;
      }//case polynomial or FRF
      
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        return x[0] + ((1.0 + x[1]) * energy);
        
      case SpecUtils::EnergyCalType::InvalidEquationType:
        break;
    }//switch( m_energy_cal->type() )
    
    assert( 0 );
    throw runtime_error( "Energy cal must be valid" );
    
    return energy;
  }//double apply_energy_cal_adjustment( double energy, const std::vector<double> &x ) const
  
  
  /** Translates from energy of m_energy_cal, to the "true" energy.  E.g., maps from energy observed
   on the spectrum, to truth energy.
   */
  double un_apply_energy_cal_adjustment( double adjusted_energy, const std::vector<double> &x ) const
  {
    assert( x.size() > 2 );
    
    if( !m_options.fit_energy_cal )
    {
      assert( fabs(x[0]) < std::numeric_limits<float>::epsilon() );
      assert( fabs(x[1]) < std::numeric_limits<float>::epsilon() );
      
      return adjusted_energy;
    }//if( we arent fitting energy cal )
    
    // Check adjustments are near the limits we placed (which was [-5,5], and [-0.01,0.01])
    assert( fabs(x[0]) <= 5.1 );
    assert( fabs(x[1]) <= 0.011 );
    
    // TODO: is this the best way to apply corrections?
    return (adjusted_energy / (1.0 + x[1])) - x[0];
  }//double un_apply_energy_cal_adjustment( double energy, const std::vector<double> &x ) const
  
  
  /** We could just return the vector of peaks from #peaks_for_energy_range, and then from their
   continuum energy range get the channel numbers of the ROI (which may vary as energy calibration
   varies), but being explicit about channel numbers used will help keep things consistent
   (currently, converting from the ROI energy range to channel numbers is a little inconsistent, or
   at least confusing, because we havent defined how to round the energy to channel number, and
   probably arent being consistent throughout the app...)
   */
  struct PeaksForEnergyRange
  {
    std::vector<PeakDef> peaks;
    size_t first_channel;
    size_t last_channel;
    bool no_gammas_in_range;
  };//struct PeaksForEnergyRange
  
  
  /** Creates the peaks for a ROI.
   
   All peaks will share a single continuum, and have amplitudes and FWHM according to input
   paramaters \c x.
   
   If energy calibration is being adjusted, peaks are returned with energy range and mean w.r.t.,
   the original energy calibration - i.e., the mean will not be that of the gamma, but the energy
   that gamma would be observed in the spectrum, with its original energy calibration.
   
   @param range The energy range to create peaks for.
   @param x The vector of parameters that specify relative activities, efficiencies, energy
          calibration, etc. of the problem.
   
   @returns A peak for each gamma, of each nuclide, as well as each free-floating peak, in the
            problem, that is within the specified energy range.
   */
  PeaksForEnergyRange peaks_for_energy_range( const RoiRangeChannels &range, const std::vector<double> &x ) const
  {
    const size_t num_channels = range.num_channels;
    
    // We will use "adjusted" to refer to energies that have been mapped into the spectrums original
    //  energy calibrations
    const double adjusted_lower_energy = apply_energy_cal_adjustment( range.lower_energy, x );
    const double adjusted_upper_energy = apply_energy_cal_adjustment( range.upper_energy, x );
    
    pair<size_t,size_t> channel_range = range.channel_range( adjusted_lower_energy, adjusted_upper_energy, num_channels, m_energy_cal );
    const size_t first_channel = channel_range.first;
    const size_t last_channel = channel_range.second;
    
    //const float adjusted_first_channel_lower = m_energy_cal->energy_for_channel(first_channel);
    //const float adjusted_last_channel_upper = m_energy_cal->energy_for_channel(last_channel + 0.999 );
    
    const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters( range.continuum_type ) );
    const bool is_step_continuum = PeakContinuum::is_step_continuum( range.continuum_type );
    
    size_t num_free_peak_pars = 0;
    set<const SandiaDecay::Nuclide *> nuclides_used;
    
    PeaksForEnergyRange answer;
    answer.first_channel = first_channel;
    answer.last_channel = last_channel;
    answer.no_gammas_in_range = false;
    
    vector<PeakDef> &peaks = answer.peaks;
    
    // Go through and create peaks based on rel act, eff, etc
    for( const NucInputGamma &nucinfo : m_nuclides )
    {
      const vector<NucInputGamma::EnergyYield> *gammas = nullptr;
      std::unique_ptr<vector<NucInputGamma::EnergyYield>> aged_gammas_cache;
      
      const double rel_act = relative_activity( nucinfo.nuclide, x );
      
      //cout << "peaks_for_energy_range: Relative activity of " << nucinfo.nuclide->symbol
      //     << " is " << PhysicalUnits::printToBestActivityUnits(rel_act) << endl;
      
      if( is_fixed_age(nucinfo.nuclide) )
      {
        gammas = &(nucinfo.nominal_gammas);
      }else
      {
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nucinfo.nuclide, ns_decay_act_mult, age(nucinfo.nuclide,x) );
        auto aged_gammas = mix.gammas( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
        
        for( double energy : nucinfo.gammas_to_exclude )
        {
          const size_t did_remove = nucinfo.remove_gamma( energy, aged_gammas );
          assert( did_remove );
        }
        
        aged_gammas_cache.reset( new vector<NucInputGamma::EnergyYield>(aged_gammas.size()) );
        gammas = aged_gammas_cache.get();
        vector<NucInputGamma::EnergyYield> &aged_gammas_ref = *aged_gammas_cache;
        for( size_t i = 0; i < aged_gammas.size(); ++i )
        {
          aged_gammas_ref[i].energy = aged_gammas[i].energy;
          aged_gammas_ref[i].yield = aged_gammas[i].numPerSecond / ns_decay_act_mult;
        }
      }//if( age is fixed ) / else( age may vary )
      
      assert( gammas );
      
      for( const NucInputGamma::EnergyYield &gamma : *gammas )
      {
        const double energy = gamma.energy;
        const double yield = gamma.yield;
        
        //if( nucinfo.nuclide->symbol == "Pu238" )
        //{
//#warning "Hacking Pu238 to a single gamma at 152 keV"
          //if( fabs(gamma.energy - 152.72) > 0.1 )
          //  continue;
             
        //  if( (fabs(gamma.energy - 152.72) < 0.1) || (fabs(gamma.energy - 368.65) < 0.1) )
        //    cout << "\t" << nucinfo.nuclide->symbol << ": {" << gamma.energy << "," << yield << "}\n";
        //}

        // Filter out zero-amplitude gammas
        // TODO: - come up with more intelligent lower bound of gamma rate to bother wit
        if( yield < std::numeric_limits<float>::min() )
          continue;
        
        // Filter the energy range based on true energies (i.e., not adjusted energies)
        if( (gamma.energy < range.lower_energy) || (gamma.energy > range.upper_energy) )
          continue;
        
        nuclides_used.insert( nucinfo.nuclide );
        
        // We compute the relative efficiency and FWHM based off of "true" energy
        const double rel_eff = relative_eff( gamma.energy, x );
        if( IsInf(rel_eff) || IsNan(rel_eff) )
          throw runtime_error( "peaks_for_energy_range: inf or NaN rel. eff for "
                              + std::to_string(gamma.energy) + " keV."  );

        
        const double peak_fwhm = fwhm( gamma.energy, x );
        
        if( IsInf(peak_fwhm) || IsNan(peak_fwhm) )
        {
          stringstream msg;
          msg << "peaks_for_energy_range: " << peak_fwhm << " FWHM for "
          << std::setprecision(2) << gamma.energy << " keV, from pars={";
          const size_t num_drf_par = num_parameters(m_options.fwhm_form);
          for( size_t i = 0; i < num_drf_par; ++i )
            msg << (i ? ", " : "") << x[2 + i];
          msg << "}";
          
          throw runtime_error( msg.str() );
        }//if( IsInf(peak_fwhm) || IsNan(peak_fwhm) )
        
        if( peak_fwhm < 0.001 )
          throw runtime_error( "peaks_for_energy_range: for peak at " + std::to_string(gamma.energy)
                            + " keV, FWHM=" + std::to_string(peak_fwhm) + " which is too small." );
        
        const double peak_amplitude = rel_act * m_live_time * rel_eff * yield;
        if( IsInf(peak_amplitude) || IsNan(peak_amplitude) )
          throw runtime_error( "peaks_for_energy_range: inf or NaN peak amplitude for "
                              + std::to_string(gamma.energy) + " keV.");
        
        const double peak_mean = apply_energy_cal_adjustment( gamma.energy, x );
        if( IsInf(peak_mean) || IsNan(peak_mean) )
          throw runtime_error( "peaks_for_energy_range: inf or NaN peak mean for "
                              + std::to_string(gamma.energy) + " keV.");
        
        if( peak_amplitude < std::numeric_limits<float>::min() )
        {
          //cout << "peaks_for_energy_range: Peak at " << gamma.energy << " keV for " << nucinfo.nuclide->symbol
          //<< " has a mean of " << peak_mean << " keV, FWHM=" << peak_fwhm << ", RelEff=" << rel_eff
          //<< ", and AMP=" << peak_amplitude << ", rel_act=" << rel_act << ", yield="
          //<< yield << ", m_live_time=" << m_live_time << endl;
          continue;
        }
        
        peaks.emplace_back( peak_mean, peak_fwhm/2.35482, peak_amplitude );
        
        PeakDef &new_peak = peaks.back();

        // TODO: (low priority) Instead of calling #PeakDef::findNearestPhotopeak, should instead keep gamma source information from decay... and should actually implement this in SandiaDecay...
        const bool xray_only = false;
        size_t transition_index = 0;
        const SandiaDecay::Transition *transition = nullptr;
        PeakDef::SourceGammaType gamma_type;
        const double window_half_width = -1.0;
        PeakDef::findNearestPhotopeak( nucinfo.nuclide, gamma.energy, window_half_width, xray_only,
                                       transition, transition_index, gamma_type );
        assert( transition );
        
        const double matched_energy = transition->products[transition_index].energy;
        assert( fabs(matched_energy - gamma.energy) < 0.0001 );
        
        if( transition )
          new_peak.setNuclearTransition( nucinfo.nuclide, transition,
                                           static_cast<int>(transition_index), gamma_type );
        
        if( !nucinfo.peak_color_css.empty() )
          new_peak.setLineColor( Wt::WColor( Wt::WString::fromUTF8(nucinfo.peak_color_css) ) );
      }//for( const SandiaDecay::EnergyRatePair &gamma : gammas )
    }//for( const NucInputGamma &nucinfo : m_nuclides )
    
    const size_t free_peaks_start_index = 2
                                          + num_parameters( m_options.fwhm_form )
                                          + (m_options.rel_eff_eqn_order + 1)
                                          + 2*m_nuclides.size();
    
    for( size_t index = 0; index < m_extra_peaks.size(); ++index )
    {
      const RelActCalcAuto::FloatingPeak &peak = m_extra_peaks[index];
      
      if( (peak.energy < range.lower_energy) || (peak.energy > range.upper_energy) )
        continue;
      
      const size_t amp_index = free_peaks_start_index + 2*index + 0;
      const size_t fwhm_index = free_peaks_start_index + 2*index + 1;
      assert( fwhm_index < x.size() );
      
      num_free_peak_pars += 1;
      const double peak_amp = x[amp_index];
      double peak_fwhm = x[fwhm_index];
      
      if( !peak.release_fwhm )
      {
        const double true_energy = un_apply_energy_cal_adjustment( peak.energy, x );
        peak_fwhm = fwhm( true_energy, x );
      }else
      {
        num_free_peak_pars += 1;
      }
      
      //cout << "peaks_for_energy_range: free peak at " << peak.energy << " has a FWHM=" << peak_fwhm << " and AMP=" << peak_amp << endl;
      
      peaks.emplace_back( peak.energy, peak_fwhm/2.35482, peak_amp );
    }//for( const RelActCalcAuto::FloatingPeak &peak : m_extra_peaks )

    if( peaks.empty() )
    {
      // We will add a zero-amplitude peak, and fit the continuum, so this way we account for this
      //  region, even if age or something has drove all the gammas out of this regions
      cerr << "peaks_for_energy_range: no peaks in range [" << range.lower_energy << ", "
           << range.upper_energy << "] keV." << endl;
      
      answer.no_gammas_in_range = true;
      
      const double middle_energy = 0.5*(adjusted_lower_energy + adjusted_upper_energy);
      const double middle_fwhm = fwhm( middle_energy, x );
      const double amplitude = 0.0;
      peaks.emplace_back( middle_energy, middle_fwhm/2.35482, amplitude );
    }//if( peaks.empty() )
    
    
    shared_ptr<PeakContinuum> continuum = peaks.front().continuum();
    continuum->setType( range.continuum_type );
    continuum->setRange( adjusted_lower_energy, adjusted_upper_energy );
    
    for( size_t i = 0; i < peaks.size(); ++i )
      peaks[i].setContinuum( continuum );
    
    assert( !peaks.empty() );
    assert( num_channels == ((1 + last_channel) - first_channel) );
    assert( first_channel < m_channel_counts.size() );
    assert( last_channel < m_channel_counts.size() );
    assert( m_energy_cal && m_energy_cal->channel_energies() );
    assert( m_energy_cal->channel_energies()->size() >= m_channel_counts.size() );
    
    const double ref_energy = adjusted_lower_energy;
    const float * const data = &(m_channel_counts[first_channel]);
    const float * const energies = &((*m_energy_cal->channel_energies())[first_channel]);
    
    vector<double> dummy_amps, continuum_coeffs, dummy_amp_uncert, continuum_uncerts;
    
    const double chi2 = fit_amp_and_offset( energies, data, num_channels, num_polynomial_terms,
                                           is_step_continuum, ref_energy, {}, {}, peaks, dummy_amps,
                                           continuum_coeffs, dummy_amp_uncert, continuum_uncerts );
    
    // TODO: - currently not defining degrees of freedom well - not using number of relative efficiency terms, or FWHM terms at all, and just blindly using all activity and free peak terms.
    const double approx_dof = 1.0*range.num_channels - nuclides_used.size() - num_polynomial_terms - num_free_peak_pars;
    
    const double chi2Dof = chi2 / approx_dof;
    
    peaks[0].continuum()->setType( range.continuum_type );
    peaks[0].continuum()->setParameters( ref_energy, continuum_coeffs, continuum_uncerts );
    peaks[0].continuum()->setRange( adjusted_lower_energy, adjusted_upper_energy );
        
    for( size_t j = 0; j < peaks.size(); ++j )
      peaks[j].set_coefficient( chi2Dof, PeakDef::Chi2DOF );

    std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMean );
    
    
    //cerr << "peaks_for_energy_range: returning " << peaks.size() << " in range [" << range.lower_energy << ", "
    //<< range.upper_energy << "] keV." << endl;
    
    return answer;
  }//PeaksForEnergyRange peaks_for_energy_range( const double lower_energy, const double upper_energy ) const
  
  
  
  void eval( const std::vector<double> &x, double *residuals ) const
  {
    m_ncalls += 1;
    
    assert( x.size() == number_parameters() );
    assert( residuals );
    assert( !m_energy_ranges.empty() );
    
    // We'll do the simplest parallelization we can by computing peaks in multiple threads; this
    //  actually seems to be reasonably effective in filling up the CPU cores during fitting.
    //  (setting ceres_options.num_threads >1 doesnt seem to do much (any?) good)
    SpecUtilsAsync::ThreadPool pool;
    vector<PeaksForEnergyRange> peaks_in_ranges( m_energy_ranges.size() );
    for( size_t i = 0; i < m_energy_ranges.size(); ++i )
    {
      pool.post( [i,&peaks_in_ranges,this,&x](){
        peaks_in_ranges[i] = peaks_for_energy_range( m_energy_ranges[i], x );
      } );
    }//
    pool.join();
    
    
    size_t residual_index = 0;
    for( size_t i = 0; i < m_energy_ranges.size(); ++i )
    {
      const RoiRangeChannels &energy_range = m_energy_ranges[i];
      const PeaksForEnergyRange &info = peaks_in_ranges[i];
      
      assert( info.peaks.size() );
      assert( info.first_channel < m_channel_counts.size() );
      assert( info.last_channel < m_channel_counts.size() );
      assert( m_channel_counts.size() == m_channel_count_uncerts.size() );
      
      shared_ptr<const PeakContinuum> continuum = info.peaks.at(0).continuum();
      
      assert( continuum );
      
      for( size_t channel = info.first_channel; channel <= info.last_channel; ++channel )
      {
        const double data_counts = m_channel_counts[channel];
        const double data_uncert = m_channel_count_uncerts[channel];
        const double lower_energy = m_energy_cal->energy_for_channel(channel);
        const double upper_energy = m_energy_cal->energy_for_channel(channel + 1);
        const double continuum_counts = continuum->offset_integral(lower_energy, upper_energy, m_spectrum);
        
        double gaussian_area = 0.0;
        for( const PeakDef &peak : info.peaks )
          gaussian_area += peak.gauss_integral( lower_energy, upper_energy );
          
        residuals[residual_index] = (data_counts - continuum_counts - gaussian_area) / data_uncert;
        
        //cout << "eval: For channel " << channel << " (" << lower_energy << " to " << upper_energy
        //<< " keV) there are " << data_counts << " data counts, and currently fitting for "
        //<< gaussian_area << " gaussian counts, and " << continuum_counts
        //<< " continuum counts (total fit: " << (gaussian_area + continuum_counts) << ")" << endl;
        
        ++residual_index;
      }//for( loop over channels )
    }//for( loop over m_energy_ranges )
    
    
    // Now make sure the relative efficiency curve is anchored to 1.0 (this removes the degeneracy
    //  between the relative efficiency amplitude, and the relative activity amplitudes).
    const double lowest_energy = m_energy_ranges.front().lower_energy;
    const double lowest_energy_rel_eff = relative_eff( lowest_energy, x );
    residuals[residual_index] = m_rel_eff_anchor_enhancement * (1.0 - lowest_energy_rel_eff);
    ++residual_index;
    
    assert( residual_index == number_residuals() );
  }//void eval( const std::vector<double> &x, double *residuals ) const
  
  
  virtual double operator()( const std::vector<double> &x ) const
  {
    vector<double> residuals( number_residuals(), 0.0 );
    try
    {
      eval( x, residuals.data() );
    }catch( std::exception &e )
    {
      cerr << "RelActAutoCostFcn::operator() caught: " << e.what() << endl;
      return std::numeric_limits<double>::max();
    }
    
    double chi2 = 0.0;
    for( const double &d : residuals )
      chi2 += d*d;
    
    return chi2;
  }//operator() - for minuit
  
  
  // For Minuit2
  virtual double Up() const
  {
    return 1.0;
  }

  
  // The return value indicates whether the computation of the
  // residuals and/or jacobians was successful or not.
  bool operator()( double const *const *parameters, double *residuals ) const
  {
    try
    {
      const size_t num_pars = number_parameters();
      
      vector<double> pars( num_pars, 0.0 );
      for( size_t i = 0; i < num_pars; ++i )
        pars[i] = parameters[i][0];
        
      eval( pars, residuals );
    }catch( std::exception &e )
    {
      cerr << "RelActAutoCostFcn::operator() caught: " << e.what() << endl;
      return false;
    }
    
    return true;
  };//bool operator() - for Ceres
};//struct RelActAutoCostFcn

}


namespace RelActCalcAuto
{





int run_test()
{
  try
  {
    //const char *xml_file_path = "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/simple_pu_test.xml";
    const char *xml_file_path = "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/thor_core_614_668_kev_test.xml";
    
    rapidxml::file<char> input_file( xml_file_path );
    
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( input_file.data() );
    
    
    const auto base_node = get_required_node(&doc, "RelActCalcAuto");
  
    auto open_file = [xml_file_path]( const string &filename ) -> shared_ptr<SpecMeas> {
      
      if( filename.empty() )
        return nullptr;
      
      string file_path = filename;
      if( !SpecUtils::is_file(file_path) )
      {
        // We'll let the test XML specify paths, relative to itself.
        string base_path = SpecUtils::parent_path( xml_file_path );
        file_path = SpecUtils::append_path( base_path, filename );
      }
      
      if( !SpecUtils::is_file(file_path) )
        throw runtime_error( "Spec file at '" + filename + "', or '" + file_path + "' is not a file" );
      
      auto spec = make_shared<SpecMeas>();
      if( !spec->load_file( file_path, SpecUtils::ParserType::Auto ) )
        throw runtime_error( "File at '" + file_path + "' could not be parsed" );
      
      if( spec->num_measurements() != 1 )
        throw runtime_error( "File at '" + file_path + "' had "
                             + to_string(spec->num_measurements()) + " measurements; requires 1.");
      
      return spec;
    };//open_file() lamda
    
    const auto foreground_node = get_required_node(base_node, "ForegroundFileName");
    const string foreground_file = SpecUtils::xml_value_str(foreground_node);
    
    const auto background_node = XML_FIRST_INODE(base_node, "BackgroundFileName");
    const string background_file = SpecUtils::xml_value_str(background_node);
    
    auto foreground = open_file( foreground_file );
    auto background = open_file( background_file );
    
    const string output_html_name = SpecUtils::xml_value_str( XML_FIRST_NODE(base_node, "OutputHtmlFileName") );
    
    const auto options_node = get_required_node(base_node, "Options");
    Options options;
    options.fromXml( options_node );

    bool extract_info_from_n42 = false;
    try
    {
      extract_info_from_n42 = get_bool_node_value(base_node, "RoiAndNucsFromFile");
    }catch(std::exception &)
    {
      //<RoiAndNucsFromFile> is optional, so will get here if it doesnt exist (or invalid value in it)
    }
    
    
    vector<RoiRange> energy_ranges;
    vector<NucInputInfo> nuclides;
    vector<FloatingPeak> extra_peaks;
    vector<string> input_warnings;
    
    
    if( extract_info_from_n42 )
    {
      // Do a quick sanity check that we arent accidentally over-specifying things or something
      if( XML_FIRST_INODE(base_node, "RoiRangeList")
          || XML_FIRST_INODE(base_node, "NucInputInfoList")
          || XML_FIRST_INODE(base_node, "FloatingPeakList") )
      {
        throw runtime_error( "When RoiAndNucsFromFile is true, <RoiRangeList>, <NucInputInfoList>,"
                              " or <FloatingPeakList> may not also be specified" );
      }
      
      shared_ptr<deque<shared_ptr<const PeakDef>>> peaks = foreground->peaks( foreground->sample_numbers() );
      if( !peaks || peaks->empty() )
        throw runtime_error( "No peaks specified in file even though <RoiAndNucsFromFile> option is set true." );
      
      set<const SandiaDecay::Nuclide *> nucs;
      map<int,set<const SandiaDecay::Nuclide *>> nuc_per_el;
      set<shared_ptr<const PeakContinuum>> continuums;
      for( const auto &peak : *peaks )
      {
        continuums.insert( peak->continuum() );
        const auto n = peak->parentNuclide();
        if( n )
        {
          if( !nucs.count(n) )
          {
            nucs.insert( n );
            nuc_per_el[n->atomicNumber].insert( n );
            
            NucInputInfo nucinfo;
            // We dont actually know the age the user wanted; we could look for the shielding/source
            //  model, but thats to much to bother with for our purposes right now.
            nucinfo.age = PeakDef::defaultDecayTime( n, nullptr );
            nucinfo.fit_age = false;
            //nucinfo.gammas_to_exclude
            nucinfo.nuclide = n;
            nucinfo.peak_color_css = peak->lineColor().cssText();
            
            nuclides.push_back( nucinfo );
          }//if( we havnet seen this nuclide yet )
        }else
        {
          // If peak doesnt have a nuclide, it will be a floating peak
          FloatingPeak fp;
          fp.energy = peak->mean();
          fp.release_fwhm = false;
          
          extra_peaks.push_back( fp );
        }//if( peak has a nuclide ) / else ( a floating nuclide )
      }//for( loop over user peaks )
      
      // Make sure isotopes of the same element all have the same age, if this option is selected
      if( options.nucs_of_el_same_age )
      {
        for( const auto &en_nuc : nuc_per_el )
        {
          if( en_nuc.second.size() < 2 )
            continue;
          
          // We'll use the min-default-age of all the nuclides for this element, and only up to
          //   20 years - with not much good justification for this choice, other than we're
          //   probably interested in human-made nuclides (hence the 20 years - which is a
          //   reasonable default for like U and Pu), and we want to make sure to include all the
          //   nuclides, so we dont want a short-half-lived nuclides to be excluded because it is a
          //   ton of half lives (hence the std::min).
          //   This isnt perfect, but good enough for code development, which is probably the only
          //   time we'll be here anyway.
          double common_age = 20*PhysicalUnits::year;
          for( const auto &n : nuclides )
          {
            if( en_nuc.second.count(n.nuclide) )
              common_age = std::max( common_age, n.age );
          }
          
          for( auto &n : nuclides )
          {
            if( en_nuc.second.count(n.nuclide) )
              n.age = common_age;
          }
        }//for( const auto &en_nuc : nuc_per_el )
      }//if( options.nucs_of_el_same_age )
      
      for( const auto &cont : continuums )
      {
        RoiRange range;
        range.lower_energy = cont->lowerEnergy();
        range.upper_energy = cont->upperEnergy();
        range.continuum_type = cont->type();
        range.force_full_range = true;
        range.allow_expand_for_peak_width = false;
        
        energy_ranges.push_back( range );
      }//for( loop over ROIs )
    }else
    {
      const auto roi_ranges_node = get_required_node(base_node, "RoiRangeList");
      XML_FOREACH_DAUGHTER(roi_range_node, roi_ranges_node, "RoiRange")
      {
        RoiRange range;
        range.fromXml( roi_range_node );
        energy_ranges.push_back( range );
      }
      
      const auto nucs_node = get_required_node(base_node, "NucInputInfoList");
      XML_FOREACH_DAUGHTER(nuc_node, nucs_node, "NucInputInfo")
      {
        NucInputInfo nuc;
        nuc.fromXml( nuc_node );
        nuclides.push_back( nuc );
      }
      
      auto float_peak_node = XML_FIRST_NODE(base_node, "FloatingPeakList");
      XML_FOREACH_DAUGHTER(float_peak_node, nucs_node, "FloatingPeak")
      {
        FloatingPeak peak;
        peak.fromXml( float_peak_node );
        extra_peaks.push_back( peak );
      }
    }//if( extract_info_from_n42 )
    
    
    std::sort( begin(energy_ranges), end(energy_ranges), []( const RoiRange &lhs, const RoiRange &rhs ) -> bool {
      return lhs.lower_energy < rhs.lower_energy;
    });
    
    std::sort( begin(extra_peaks), end(extra_peaks), []( const FloatingPeak &lhs, const FloatingPeak &rhs ) -> bool {
      return lhs.energy < rhs.energy;
    });
    
    
    // Make sure ranges dont overlap - if they do, split the difference.
    //  This only catches simple small
    for( size_t i = 1; i < energy_ranges.size(); ++i )
    {
      RoiRange &prev_range = energy_ranges[i-1];
      RoiRange &this_range = energy_ranges[i];
      if( prev_range.upper_energy > this_range.lower_energy )
      {
        const double diff = prev_range.upper_energy - this_range.lower_energy;
        
        const double prev_kev = prev_range.upper_energy - prev_range.lower_energy;
        const double this_kev = this_range.upper_energy - this_range.lower_energy;
        
        char buffer[512] = { '\0' };
        snprintf( buffer, sizeof(buffer), "Energy range [%.1f, %.1f] keV, overlaps with range"
                 " [%.1f, %.1f] keV",
                 prev_range.lower_energy, prev_range.upper_energy,
                 this_range.lower_energy, this_range.upper_energy );
        
        if( diff > 0.2*(std::min(prev_kev,this_kev)) )
          throw runtime_error( string(buffer) + " - and overlap is too large to fix." );
        
        input_warnings.push_back( string(buffer) + " - will remove overlap down the middle." );

        const double adjust = 0.51*diff;
        
        // TODO: should check that adjusting the ROIs to be not overlapping, doesnt result in invalid ROIs
        prev_range.upper_energy -= adjust;
        this_range.lower_energy += adjust;
      }//if( prev_range.upper_energy > this_range.lower_energy )
    }//for( size_t i = 1; i < energy_ranges.size(); ++i )
    
    
    if( energy_ranges.empty() )
      throw runtime_error( "No RoiRanges specified" );
    
    if( nuclides.empty() )
      throw runtime_error( "No nuclides specified" );
    
    
    // A helper function to print out the XML; right now just to stdout, but could be useful
    auto print_xml = [&](){
      try
      {
        using namespace rapidxml;
        
        rapidxml::xml_document<char> doc;
        
        xml_node<char> *base_node = doc.allocate_node( node_element, "RelActCalcAuto" );
        doc.append_node( base_node );
        
        xml_attribute<char> *attrib = doc.allocate_attribute( "version", "0" );
        base_node->append_attribute( attrib );
        
        const char *strvalue = doc.allocate_string( foreground_node->value(), foreground_node->value_size() );
        xml_node<char> *node = doc.allocate_node( node_element, "ForegroundFileName", strvalue );
        base_node->append_node( node );
        
        if( background_file.size() )
        {
          strvalue = doc.allocate_string( background_file.c_str() );
          node = doc.allocate_node( node_element, "BackgroundFileName", strvalue );
          base_node->append_node( node );
        }
        
        if( !output_html_name.empty() )
        {
          strvalue = doc.allocate_string( output_html_name.c_str() );
          node = doc.allocate_node( node_element, "OutputHtmlFileName", strvalue );
          base_node->append_node( node );
          append_attrib( node, "remark", "Can be either an absolute path, or relative to this files path." );
        }
        
        options.toXml( base_node );
        
        // Note that when we write <RoiAndNucsFromFile> as "true" it actually makes the XML so it
        //  wont be read back in
        
        if( extract_info_from_n42 )
        {
          node = doc.allocate_node( node_comment, "", "You must change RoiAndNucsFromFile to false,"
                                   " or delete RoiRangeList, NucInputInfoList, and FloatingPeakList"
                                   " nodes, or else you cant use this XML as input" );
          base_node->append_node( node );
        }//if( extract_info_from_n42 )
        
        node = doc.allocate_node( node_comment, "", "Optionally, we can extract ROI and Nuclide"
                                 " information from N42 files saved by InterSpec" );
        base_node->append_node( node );
        
        node = doc.allocate_node( node_element, "RoiAndNucsFromFile", (extract_info_from_n42 ? "true" : "false") );
        base_node->append_node( node );
        
        node = doc.allocate_node( node_element, "RoiRangeList" );
        base_node->append_node( node );
        for( const auto &range : energy_ranges )
          range.toXml( node );
        
        node = doc.allocate_node( node_element, "NucInputInfoList" );
        base_node->append_node( node );
        for( const auto &nuc : nuclides )
          nuc.toXml( node );
        
        if( !extra_peaks.empty() )
        {
          node = doc.allocate_node( node_element, "FloatingPeakList" );
          base_node->append_node( node );
          for( const auto &peak : extra_peaks )
            peak.toXml( node );
        }
       
        cout << "RelActCalcAuto input XML is:\n\n";
        rapidxml::print( std::cout, doc, 0 );
      }catch( std::exception &e )
      {
        cerr << "Failed to to write XML for RelActCalcAuto input: " << e.what() << endl;
      }
    };//print_xml(...)
    
    if( extract_info_from_n42 )
      print_xml();
    
    assert( foreground && (foreground->num_measurements() == 1) );
    
    shared_ptr<const DetectorPeakResponse> drf = foreground->detector();
    shared_ptr<const SpecUtils::Measurement> fore_meas = foreground->measurements()[0];
    shared_ptr<const SpecUtils::Measurement> back_meas;
    if( background )
      back_meas = background->measurements()[0];

    vector<shared_ptr<const PeakDef>> all_peaks;
    RelActAutoSolution solution = solve( options, energy_ranges, nuclides, extra_peaks, fore_meas, back_meas, drf, all_peaks );
    
    std::reverse( begin(input_warnings), end(input_warnings) );
    for( const string &w : input_warnings )
      solution.m_warnings.insert( std::begin(solution.m_warnings), "Preparing input: " + w );
      
    
    //Print out summary, etc.
    
    cout << "\n\n\n\n-----------------------------------------------------------\n\n"
    << "Done fitting:\n";
    
    solution.print_summary( std::cout );
    
    if( !output_html_name.empty() )
    {
      string out_file = output_html_name;
      if( !SpecUtils::is_absolute_path(out_file) )
        out_file = SpecUtils::append_path(SpecUtils::parent_path(xml_file_path), out_file);
      
      if( !SpecUtils::iequals_ascii( SpecUtils::file_extension(out_file), ".html") )
        out_file += ".html";
      
      try
      {
        ofstream output( out_file.c_str(), ios::binary | ios::out );
        
        if( output.is_open() )
          solution.print_html_report( output );
        else
          cerr << "\n\nFailed to open output html file '" << output_html_name << "' (as '"
          << out_file << "')" << endl;
      }catch( std::exception &e )
      {
        cerr << "Failed to write '" << out_file << "': " << e.what() << endl;
      }
    }//if( !output_file_name.empty() )
    
    
    cout << "\n\n\n";
    
  }catch( std::exception &e )
  {
    cerr << "RelAct test failed: " << e.what() << endl;
    return EXIT_FAILURE;
  }// try / catch
  
  return EXIT_SUCCESS;
}//void run_test()



const char *to_str( const FwhmForm form )
{
  switch( form )
  {
    case FwhmForm::Gadras:       return "Gadras";
    case FwhmForm::Polynomial_2: return "Polynomial_2";
    case FwhmForm::Polynomial_3: return "Polynomial_3";
    case FwhmForm::Polynomial_4: return "Polynomial_4";
    case FwhmForm::Polynomial_5: return "Polynomial_5";
    case FwhmForm::Polynomial_6: return "Polynomial_6";
  }//switch( form )
  
  assert( 0 );
  return "";
}//to_str( const FwhmForm form )


FwhmForm fwhm_form_from_str( const char *str )
{
  const size_t str_len = strlen(str);
  
  for( int iform = 0; iform <= static_cast<int>(FwhmForm::Polynomial_6); iform += 1 )
  {
    const FwhmForm x = FwhmForm(iform);
    const char *form_str = to_str( x );
    const bool case_sensitive = false;
    const size_t form_str_len = strlen(form_str);
    if( rapidxml::internal::compare(str, str_len, form_str, form_str_len, case_sensitive) )
       return x;
  }
  
  throw runtime_error( "fwhm_form_from_str(...): invalid input string '" + std::string(str) + "'" );
  return FwhmForm::Gadras;
}//FwhmForm fwhm_form_from_str(str)



void RoiRange::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "RoiRange::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "RoiRange" );
  parent->append_node( base_node );
  
  append_version_attrib( base_node, RoiRange::sm_xmlSerializationVersion );
  
  append_float_node( base_node, "LowerEnergy", lower_energy );
  append_float_node( base_node, "UpperEnergy", upper_energy );
  
  const char *cont_type_str = PeakContinuum::offset_type_str( continuum_type );
  append_string_node( base_node, "ContinuumType", cont_type_str );
  
  append_bool_node( base_node, "ForceFullRange", force_full_range );
  append_bool_node( base_node, "AllowExpandForPeakWidth", allow_expand_for_peak_width );
}//RoiRange::toXml(...)





void RoiRange::fromXml( const rapidxml::xml_node<char> *range_node )
{
  try
  {
    if( !range_node )
      throw runtime_error( "invalid input" );
    
    if( !rapidxml::internal::compare( range_node->name(), range_node->name_size(), "RoiRange", 8, false ) )
      throw std::logic_error( "invalid input node name" );
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( RoiRange::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    check_xml_version( range_node, RoiRange::sm_xmlSerializationVersion );
    
    lower_energy = get_float_node_value( range_node, "LowerEnergy" );
    upper_energy = get_float_node_value( range_node, "UpperEnergy" );
    
    const rapidxml::xml_node<char> *cont_type_node = XML_FIRST_NODE( range_node, "ContinuumType" );
    const string cont_type_str = SpecUtils::xml_value_str( cont_type_node );
    continuum_type = PeakContinuum::str_to_offset_type_str( cont_type_str.c_str(), cont_type_str.size() );
    
    force_full_range = get_bool_node_value( range_node, "ForceFullRange" );
    allow_expand_for_peak_width = get_bool_node_value( range_node, "AllowExpandForPeakWidth" );
  }catch( std::exception &e )
  {
    throw runtime_error( "RoiRange::fromXml(): " + string(e.what()) );
  }
}//void fromXml( const ::rapidxml::xml_node<char> *parent )



void NucInputInfo::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( nuclide );
  if( !nuclide )
    throw runtime_error( "NucInputInfo::toXml: null nuclide." );
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "NucInputInfo::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "NucInputInfo" );
  parent->append_node( base_node );
  
  append_version_attrib( base_node, NucInputInfo::sm_xmlSerializationVersion );
  
  append_string_node( base_node, "Nuclide", nuclide->symbol.c_str() );
  const string age_str = PhysicalUnits::printToBestTimeUnits(age, 10);
  append_string_node( base_node, "Age", age_str);
  append_bool_node( base_node, "FitAge", fit_age );
  
  if( !gammas_to_exclude.empty() )
  {
    xml_node<char> *exclude_node = doc->allocate_node( node_element, "GammasToExclude" );
    base_node->append_node( exclude_node );
    for( const double energy : gammas_to_exclude )
      append_float_node( exclude_node, "Energy", energy );
  }//if( !gammas_to_exclude.empty() )
  
  if( !peak_color_css.empty() )
    append_string_node( base_node, "PeakColorCss", peak_color_css);
}//void NucInputInfo::toXml(...)


void NucInputInfo::fromXml( const ::rapidxml::xml_node<char> *nuc_info_node )
{
  try
  {
    if( !nuc_info_node )
      throw runtime_error( "invalid input" );
    
    if( !rapidxml::internal::compare( nuc_info_node->name(), nuc_info_node->name_size(), "NucInputInfo", 12, false ) )
      throw std::logic_error( "invalid input node name" );
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( NucInputInfo::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    check_xml_version( nuc_info_node, NucInputInfo::sm_xmlSerializationVersion );
    
    const rapidxml::xml_node<char> *nuc_node = XML_FIRST_NODE( nuc_info_node, "Nuclide" );
    if( !nuc_node )
      throw runtime_error( "Missing 'Nuclide' node." );
    
    const string nuc_str = SpecUtils::xml_value_str( nuc_node );
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    nuclide = db->nuclide( nuc_str );
    if( !nuclide )
      throw runtime_error( "Invalid nuclide '" + nuc_str + "'" );
    
    const auto age_node = XML_FIRST_INODE( nuc_info_node, "Age");
    const string age_str = SpecUtils::xml_value_str( age_node );
    if( age_str.empty() || SpecUtils::iequals_ascii(age_str, "default") )
      age = PeakDef::defaultDecayTime(nuclide);
    else
      age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( age_str, nuclide->halfLife );
    
    if( XML_FIRST_NODE(nuc_info_node, "FitAge") )
      fit_age = get_bool_node_value( nuc_info_node, "FitAge" );
    else
      fit_age = false;
    
    gammas_to_exclude.clear();
    const rapidxml::xml_node<char> *exclude_node = XML_FIRST_NODE( nuc_info_node, "GammasToExclude" );
    XML_FOREACH_DAUGHTER( energy_node, exclude_node, "Energy" ) //ok if exclude_node is nullptr
    {
      const string strval = SpecUtils::xml_value_str(energy_node);
      double energy;
      if( !(stringstream(strval) >> energy) )
        throw runtime_error( "Invalid exclude energy '" + strval + "'" );
      gammas_to_exclude.push_back( energy );
    }//foreach( <Energy> node under <GammasToExclude> )
    
    const auto color_node = XML_FIRST_INODE( nuc_info_node, "PeakColorCss");
    peak_color_css = SpecUtils::xml_value_str( color_node );
  }catch( std::exception &e )
  {
    throw runtime_error( "NucInputInfo::fromXml(): " + string(e.what()) );
  }
}//void fromXml( const ::rapidxml::xml_node<char> *parent )



void FloatingPeak::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "FloatingPeak::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "FloatingPeak" );
  parent->append_node( base_node );
  
  append_version_attrib( base_node, FloatingPeak::sm_xmlSerializationVersion );
  append_float_node( base_node, "Energy", energy );
  append_bool_node( base_node, "ReleaseFwhm", release_fwhm );
}//void FloatingPeak::toXml(...)


void FloatingPeak::fromXml( const ::rapidxml::xml_node<char> *parent )
{
  try
  {
    if( !parent )
      throw runtime_error( "invalid input" );
    
    if( !rapidxml::internal::compare( parent->name(), parent->name_size(), "FloatingPeak", 12, false ) )
      throw std::logic_error( "invalid input node name" );
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( NucInputInfo::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    check_xml_version( parent, FloatingPeak::sm_xmlSerializationVersion );
    energy = get_float_node_value( parent, "Energy" );
    release_fwhm = get_bool_node_value( parent, "ReleaseFwhm" );
  }catch( std::exception &e )
  {
    throw runtime_error( "FloatingPeak::fromXml(): " + string(e.what()) );
  }
}//void fromXml(...)


void Options::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "Options::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "Options" );
  parent->append_node( base_node );
  
  append_version_attrib( base_node, Options::sm_xmlSerializationVersion );
  append_bool_node( base_node, "FitEnergyCal", fit_energy_cal );
  append_bool_node( base_node, "NucsOfElSameAge", nucs_of_el_same_age );
  
  const char *rell_eff_eqn_str = RelActCalc::to_str( rel_eff_eqn_type );
  auto rel_eff_node = append_string_node( base_node, "RelEffEqnType", rell_eff_eqn_str);
  append_attrib( rel_eff_node, "remark", "Possible values: Possible values: LnX, LnY, LnXLnY, \"FRAM Empirical\"" );
  
  append_int_node( base_node, "RelEffEqnOrder", rel_eff_eqn_order);
  
  const char *fwhm_form_str = to_str( fwhm_form );
  auto fwhm_node = append_string_node( base_node, "FwhmForm", fwhm_form_str );
  append_attrib( fwhm_node, "remark", "Possible values: Gadras, Polynomial_2, Polynomial_3, Polynomial_4, Polynomial_5, Polynomial_6" );
  
  append_string_node( base_node, "Title", spectrum_title );
}//void Options::toXml(...)


void Options::fromXml( const ::rapidxml::xml_node<char> *parent )
{
  try
  {
    if( !parent )
      throw runtime_error( "invalid input" );
    
    if( !rapidxml::internal::compare( parent->name(), parent->name_size(), "Options", 7, false ) )
      throw std::logic_error( "invalid input node name" );
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( Options::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    check_xml_version( parent, Options::sm_xmlSerializationVersion );
    
    fit_energy_cal = get_bool_node_value( parent, "FitEnergyCal" );
    nucs_of_el_same_age = get_bool_node_value( parent, "NucsOfElSameAge" );
    
    const rapidxml::xml_node<char> *rel_eff_eqn_node = XML_FIRST_NODE( parent, "RelEffEqnType" );
    const string rel_eff_str = SpecUtils::xml_value_str( rel_eff_eqn_node );
    
    rel_eff_eqn_type = RelActCalc::rel_eff_eqn_form_from_str( rel_eff_str.c_str() );
    
    const rapidxml::xml_node<char> *rel_eff_order_node = XML_FIRST_NODE( parent, "RelEffEqnOrder" );
    const string rel_eff_order_str = SpecUtils::xml_value_str( rel_eff_order_node );
    if( !(stringstream(rel_eff_order_str) >> rel_eff_eqn_order) )
      throw runtime_error( "Invalid 'RelEffEqnOrder' value '" + rel_eff_order_str + "'" );
    
    const rapidxml::xml_node<char> *fwhm_node = XML_FIRST_NODE( parent, "FwhmForm" );
    const string fwhm_str = SpecUtils::xml_value_str( fwhm_node );
    fwhm_form = fwhm_form_from_str( fwhm_str.c_str() );
    
    const rapidxml::xml_node<char> *title_node = XML_FIRST_NODE( parent, "Title" );
    spectrum_title = SpecUtils::xml_value_str( title_node );
  }catch( std::exception &e )
  {
    throw runtime_error( "Options::fromXml(): " + string(e.what()) );
  }
}//void Options::fromXml( const ::rapidxml::xml_node<char> *parent )


size_t num_parameters( const FwhmForm eqn_form )
{
  switch( eqn_form )
  {
    case FwhmForm::Gadras:       return 3;
    case FwhmForm::Polynomial_2: return 2;
    case FwhmForm::Polynomial_3: return 3;
    case FwhmForm::Polynomial_4: return 4;
    case FwhmForm::Polynomial_5: return 5;
    case FwhmForm::Polynomial_6: return 6;
  }//switch( m_options.fwhm_form )
  
  assert( 0 );
  throw runtime_error( "Invalid FwhmForm" );
  return 0;
}//size_t num_parameters( const FwhmForm eqn_form )


Options::Options()
: fit_energy_cal( false ),
  nucs_of_el_same_age( false ),
  rel_eff_eqn_type( RelActCalc::RelEffEqnForm::LnX ),
  rel_eff_eqn_order( 3 ),
  fwhm_form( FwhmForm::Polynomial_2 ),
  spectrum_title( "" )
{
}

RelActAutoSolution::RelActAutoSolution()
: m_status( RelActAutoSolution::Status::NotInitiated ),
  m_error_message( "" ),
  m_foreground{ nullptr },
  m_background{ nullptr },
  m_spectrum{ nullptr },
  m_rel_eff_form( RelActCalc::RelEffEqnForm::LnX ),
  m_rel_eff_coefficients{},
  m_rel_eff_covariance{},
  m_rel_activities{},
  m_rel_act_covariance{},
  m_fwhm_form( FwhmForm::Gadras ),
  m_fwhm_coefficients{},
  m_fwhm_covariance{},
  m_floating_peaks{},
  m_fit_peaks{},
  m_input_roi_ranges{},
  m_energy_cal_adjustments{ 0.0, 0.0 },
  m_drf{ nullptr },
  m_fit_energy_cal{ false, false },
  m_chi2( -1.0 ),
  m_dof( 0 ),
  m_num_function_eval_solution( 0 ),
  m_num_function_eval_total( 0 ),
  m_num_microseconds_eval( 0 )
{
  
}

std::ostream &RelActAutoSolution::print_summary( std::ostream &out ) const
{
  const double nsec_eval = 1.0E-6*m_num_microseconds_eval;
  
  out << "Computation took " << PhysicalUnits::printToBestTimeUnits(nsec_eval)
  << " with " << m_num_function_eval_solution << " function calls to solve, and "
  << (m_num_function_eval_total - m_num_function_eval_solution)
  << " more for covariance\n";
  
  switch( m_status )
  {
    case Status::Success:
      out << "Solution was successfully found.\n";
      break;
      
    case Status::NotInitiated:
      out << "Problem was not initialized.\n";
      break;
      
    case Status::FailedToSetupProblem:
      out << "Failed to setup problem.\n";
      break;
      
    case Status::FailToSolveProblem:
      out << "Failed to solve problem.\n";
      break;
  }//switch( m_status )
  
  if( m_error_message.size() )
    out << "Error: " << m_error_message << "\n";
  
  for( const string &warning: m_warnings )
    out << "Warning: " << warning << "\n";
  
  if( m_status != Status::Success )
    return;
  
  // Rake code from RelEff
  out << "Rel. Eff. Eqn.: y = "
  << RelActCalc::rel_eff_eqn_text(m_rel_eff_form,m_rel_eff_coefficients) << "\n"
  << "Activity ratios:\n";
  
  for( size_t i = 1; i < m_rel_activities.size(); ++i )
  {
    const NuclideRelAct &nuc_i = m_rel_activities[i];
      
    for( size_t j = 0; j < i; ++j )
    {
      const NuclideRelAct &nuc_j = m_rel_activities[j];
      
      const double act_ratio = activity_ratio(nuc_i.nuclide, nuc_j.nuclide);
      const double nuc_mass_ratio = mass_ratio(nuc_i.nuclide, nuc_j.nuclide);
      
      out << "\t" << std::setw(16) << (nuc_i.nuclide->symbol + " / " + nuc_j.nuclide->symbol)
      << "\tact: " << std::setw(10) << std::setprecision(4) << act_ratio
      << "\tmass: " << std::setw(10) << std::setprecision(4) << nuc_mass_ratio
      << "\n";
      out << "\t" << std::setw(16) << (nuc_j.nuclide->symbol + " / " + nuc_i.nuclide->symbol)
      << "\tact: " << std::setw(10) << std::setprecision(4) << 1.0/act_ratio
      << "\tmass: " << std::setw(10) << std::setprecision(4) << 1.0/nuc_mass_ratio
      << "\n";
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 0; i < used_isotopes.size(); ++i )
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  map<int,vector<const SandiaDecay::Nuclide *>> nucs_of_element;
  const SandiaDecay::Nuclide *u235 = nullptr, *pu240 = nullptr;
  for( const NuclideRelAct &nuc_rel_act : m_rel_activities )
  {
    const SandiaDecay::Nuclide * const nuc = nuc_rel_act.nuclide;
    
    nucs_of_element[nuc->atomicNumber].push_back( nuc );
    
    if( (nuc->atomicNumber == 92) && (nuc->massNumber == 235) )
      u235 = nuc;
    
    if( (nuc->atomicNumber == 94) && (nuc->massNumber == 240) )
      pu240 = nuc;
  }//for( const NuclideRelAct &nuc : m_rel_activities )
  
  if( u235 )
    out << "\nEnrichment " << 100.0*mass_enrichment_fraction(u235) << "% U235\n";
  
  if( pu240 )
    out << "\nEnrichment " << 100.0*mass_enrichment_fraction(pu240) << "% Pu240\n";
  
  for( auto &v : nucs_of_element )
  {
    vector<const SandiaDecay::Nuclide *> nucs = v.second;
    if( nucs.size() < 2 )
      continue;
    
    sort( begin(nucs), end(nucs), [](auto l, auto r){ return l->massNumber < r->massNumber; } );
    
    const SandiaDecay::Element * const el = db->element( nucs.front()->atomicNumber );
    assert( el );
    
    out << "For " << el->name << ":\n";
    for( const SandiaDecay::Nuclide *nuc : nucs )
      out << std::setw(5) << nuc->symbol << ": "
      << std::setw(10) << std::setprecision(4) << 100.0*mass_enrichment_fraction(nuc) << "%"
      << " of the " << el->name << ", by mass.\n";
  }//for( auto &v : nucs_of_element )
  
  return out;
}//RelActAutoSolution::print_summary(...)


void RelActAutoSolution::print_html_report( std::ostream &out ) const
{
  // TODO: we should probably use WTemplate to form this HTML - and then we could maybe share some models or whatever with the GUI.
  
  char buffer[1024*16] = { '\0' };
  
  const float live_time = m_spectrum ? m_spectrum->live_time() : 1.0f;
  
  stringstream results_html;
  
  results_html << "<div>&chi;<sup>2</sup>=" << m_chi2 << " and there were " << m_dof
               << " DOF (&chi;<sup>2</sup>/dof=" << m_chi2/m_dof << ")</div>\n";
  
  results_html << "<div class=\"releffeqn\">Rel. Eff. Eqn: y = "
               << RelActCalc::rel_eff_eqn_text( m_rel_eff_form, m_rel_eff_coefficients )
               << "</div>\n";
  
  double sum_rel_mass = 0.0;
  for( const auto &act : m_rel_activities )
    sum_rel_mass += act.rel_activity / act.nuclide->activityPerGram();
  
  results_html << "<table class=\"nuctable resulttable\">\n";
  results_html << "  <caption>Isotope relative activities and mass fractions</caption>\n";
  results_html << "  <thead><tr>"
                  " <th scope=\"col\">Nuclide</th>"
                  " <th scope=\"col\">Rel. Act</th>"
                  " <th scope=\"col\">Total Mas Frac.</th>"
                  " <th scope=\"col\">Enrichment</th>"
                  " </tr></thead>\n";
  results_html << "  <tbody>\n";
  for( const auto &act : m_rel_activities )
  {
    const double rel_mass = act.rel_activity / act.nuclide->activityPerGram();
    
    results_html << "  <tr><td>" << act.nuclide->symbol << "</td>"
    << "<td>" << act.rel_activity << " &pm; " << act.rel_activity_uncertainty << "</td>"
    << "<td>" << (100.0*rel_mass/sum_rel_mass) << "%</td>"
    << "<td>" << (100.0*mass_enrichment_fraction(act.nuclide)) << "%</td>"
    << "</tr>\n";
  }
  results_html << "  </tbody>\n"
  << "</table>\n\n";
  
  // Make the table of mass and activity ratios
  results_html << "<table class=\"massratiotable resulttable\">\n";
  results_html << "  <caption>Mass and Activity Ratios.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Nuclides</th>"
  "<th scope=\"col\">Mass Ratio</th>"
  "<th scope=\"col\">Activity Ratio</th>"
  "</tr></thead>\n";
  results_html << "  <tbody>\n";
  for( size_t i = 1; i < m_rel_activities.size(); ++i )
  {
    for( size_t j = 0; j < i; ++j )
    {
      const auto nuc_i = m_rel_activities[i];
      const auto nuc_j = m_rel_activities[j];
      assert( nuc_i.nuclide && nuc_j.nuclide );
      
      const double act_i = nuc_i.rel_activity;
      const double act_j = nuc_j.rel_activity;
      
      const double mass_i = act_i / nuc_i.nuclide->activityPerGram();
      const double mass_j = act_j / nuc_j.nuclide->activityPerGram();
      
      snprintf( buffer, sizeof(buffer), "<tr><td>%s</td><td>%1.6G</td><td>%1.6G</td></tr>\n",
               (nuc_i.nuclide->symbol + "/" + nuc_j.nuclide->symbol).c_str(),
               (mass_i / mass_j), (act_i / act_j) );
      results_html << buffer;
      
      snprintf( buffer, sizeof(buffer), "<tr><td>%s</td><td>%1.6G</td><td>%1.6G</td></tr>\n",
               (nuc_j.nuclide->symbol + "/" + nuc_i.nuclide->symbol).c_str(),
               (mass_j / mass_i), (act_j / act_i) );
      
      results_html << buffer;
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 0; i < used_isotopes.size(); ++i )
  results_html << "  </tbody>\n"
  << "</table>\n\n";
  
  // Make table giving info on each of the result peaks
  results_html << "<table class=\"peaktable resulttable\">\n";
  results_html << "  <caption>Peaks fit in analysis.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Energy (keV)</th>"
  "<th scope=\"col\">Nuclide</th>"
  "<th scope=\"col\">Yield</th>"
  "<th scope=\"col\">Amp</th>"
  "<th scope=\"col\">&Delta;Amp</th>"
  "<th scope=\"col\">CPS/Yield</th>"
  "<th scope=\"col\">Rel Eff</th>"
  "<th scope=\"col\">&Delta;Rel Eff</th>"
  "<th scope=\"col\">Cont. Type</th>"
  "<th scope=\"col\">Cont. Range</th>"
  "</tr></thead>\n";
  results_html << "  <tbody>\n";
  
  for( const PeakDef &info : m_fit_peaks )
  {
    const SandiaDecay::Nuclide *nuc = info.parentNuclide();
    
    const NuclideRelAct *nuc_info = nullptr;
    for( const auto &rel_act : m_rel_activities )
      nuc_info = (rel_act.nuclide == nuc) ? &rel_act : nuc_info;
    
    assert( !nuc || nuc_info );
    if( nuc && !nuc_info )
      throw logic_error( "Peak with nuc not in solution?" );
    
    if( !nuc || !nuc_info )
    {
      snprintf( buffer, sizeof(buffer),
              "  <tr><td>%.2f</td><td></td><td></td><td>%1.6G</td><td>%1.6G</td>"
               "<td></td><td></td><td></td><td></td><td></td><td></td></tr>\n",
               info.mean(), info.amplitude(), info.amplitudeUncert() );
    }else
    {
      double yield = 0;
      for( const pair<double,double> &energy_br : nuc_info->gamma_energy_br )
      {
        if( fabs(energy_br.first - info.gammaParticleEnergy()) < 1.0E-4 )
          yield += energy_br.second;
      }
      
      const double rel_act = nuc_info->rel_activity;
      const double cps_over_yield = info.amplitude() / (yield * live_time);
      const double meas_rel_eff = info.amplitude() / (yield * rel_act * live_time);
      const double fit_rel_eff = RelActCalc::eval_eqn( info.mean(), m_rel_eff_form, m_rel_eff_coefficients );
      
      //assert( fabs(meas_rel_eff - fit_rel_eff) < 0.1*std::max(meas_rel_eff,fit_rel_eff) );
      double fit_rel_eff_uncert = -1.0;
      
      if( m_rel_eff_covariance.size() )
      {
        try
        {
          fit_rel_eff_uncert = RelActCalc::eval_eqn_uncertainty( info.mean(), m_rel_eff_form, m_rel_eff_covariance );
          fit_rel_eff_uncert /= fit_rel_eff;
        }catch( std::exception &e )
        {
          cerr << "Error calling RelActCalc::eval_eqn_uncertainty: " << e.what() << endl;
        }
      }//if( m_rel_eff_covariance.size() )
      
      snprintf(buffer, sizeof(buffer),"  <tr>"
               "<td>%.2f</td><td>%s</td><td>%1.3G</td><td>%1.3G</td><td>%1.2G%%</td>"
               "<td>%1.3G</td><td>%1.6G</td><td>%1.2G%%</td>"
               "<td>%s</td><td>%.1f-%.1f</td>"
               "</tr>\n",
               info.mean(),
               nuc->symbol.c_str(),
               yield,
               info.amplitude(),
               info.amplitudeUncert() / info.amplitude(),
               cps_over_yield,
               fit_rel_eff,
               100*fit_rel_eff_uncert,
               PeakContinuum::offset_type_label( info.continuum()->type() ),
               info.continuum()->lowerEnergy(),
               info.continuum()->upperEnergy()
               );
    }//if( !nuc || !nuc_info ) / else
      
    results_html << buffer;
  }//for( const RelEff::PeakInfo &info : input_peaks )
  
  results_html << "  </tbody>\n"
  << "</table>\n\n";
  
  
  if( m_fit_energy_cal[0] || m_fit_energy_cal[1] )
  {
    results_html << "<div class=\"energycal\">\n";
    
    if( m_fit_energy_cal[0] && m_fit_energy_cal[1] )
      results_html << "Fit offset adjustment of " << m_energy_cal_adjustments[0]
                   << " keV and gain multiple of " << (1.0 + m_energy_cal_adjustments[1]) << ".\n";
    else if( m_fit_energy_cal[1] )
      results_html << "Fit gain multiple of 1 " << (m_energy_cal_adjustments[1] < 0.0 ? "- " : "+ ")
                   << fabs(m_energy_cal_adjustments[1]) << ".\n";
    else // if( m_fit_energy_cal[0] )
      results_html << "Fit offset adjustment of " << m_energy_cal_adjustments[0] << " keV.\n";
    
    results_html << "</div>\n";
  }//if( m_fit_energy_cal[0] || m_fit_energy_cal[1] )
  
  
  
  
  auto html_sanitize = []( string &val ){
    // We'll do some really stupid/simple sanitization
    SpecUtils::ireplace_all( val, "&", "&amp;" );
    SpecUtils::ireplace_all( val, "<", "&lt;" );
    SpecUtils::ireplace_all( val, ">", "&gt;" );
    SpecUtils::ireplace_all( val, "'", "&#39;" );
    SpecUtils::ireplace_all( val, "\"", "&quot;" );
  };
  
  
  results_html << "<div class=\"anacommand\">\n";
  results_html << "<table class=\"optionstable\">\n";
  results_html << "  <caption>Options used for analysis.</caption>\n";
  results_html << "  <tr><th scope=\"row\">fit_energy_cal</th><td><code>" << m_options.fit_energy_cal << "</code></td></tr>\n";
  results_html << "  <tr><th scope=\"row\">nucs_of_el_same_age</th><td><code>" << m_options.nucs_of_el_same_age << "</code></td></tr>\n";
  results_html << "  <tr><th scope=\"row\">rel_eff_eqn_type</th><td><code>"
               << RelActCalc::to_str(m_options.rel_eff_eqn_type) << "</code></td></tr>\n";
  results_html << "  <tr><th scope=\"row\">rel_eff_eqn_order</th><td><code>" << m_options.rel_eff_eqn_order << "</code></td></tr>\n";
  results_html << "  <tr><th scope=\"row\">fwhm_form</th><td><code>" << to_str(m_options.fwhm_form) << "</code></td></tr>\n";
  results_html << "</table>\n\n";
  results_html << "</div>\n";
  
  
 
  if( !m_warnings.empty() )
  {
    results_html << "<div class=\"warnings\">\n"
    << "<h3>Warnings</h3>\n";
    for( string warning : m_warnings )
    {
      html_sanitize( warning );
      
      results_html << "<div class=\"warningline\">" << warning << "</div>\n";
    }//for( string warning : warnings )
    
    results_html << "</div>\n";
  }//if( !warnings.empty() )
  

  const string current_time = Wt::WDateTime::currentDateTime().toString("yyyyMMdd hh:mm:ss").toUTF8();
  
  const double nsec_eval = 1.0E-6*m_num_microseconds_eval;
  results_html << "<div class=\"anacomputetime\">"
  << "Computation took " << PhysicalUnits::printToBestTimeUnits(nsec_eval)
  << " with " << m_num_function_eval_solution << " function calls to solve, and "
  << (m_num_function_eval_total - m_num_function_eval_solution)
  << " more for covariance\n"
  << "</div>\n";
  
  
  results_html << "<div class=\"anatime\">Analysis performed " << current_time
  << " with code compiled " __TIMESTAMP__
  << "</div>\n";
  
  
  // TODO: figure out how to reasonably plot RelEff values
  stringstream rel_eff_plot_values, add_rel_eff_plot_css;
  
  rel_eff_plot_values << "[";
  //rel_eff_plot_values << "{energy: 185, counts: 100, counts_uncert: 10, eff: 1.0, eff_uncert: 0.1, nuc_info: [{nuc: \"U235\", br: 0.2, rel_act: 300}]}"
  //", {energy: 1001, counts: 10, counts_uncert: 3, eff: 0.1, eff_uncert: 0.1, nuc_info: [{nuc: \"U238\", br: 0.1, rel_act: 100}]}";
  
  set<const SandiaDecay::Nuclide *> nuclides_with_colors, all_nuclides;
  
  size_t num_rel_eff_data_points = 0;
  for( const PeakDef &peak : m_fit_peaks )
  {
    // Skip peaks that are essentially zero amplitide
    if( peak.amplitude() < 0.1 )
      continue;
  
    const SandiaDecay::Nuclide *nuc = peak.parentNuclide();
    if( !nuc )
      continue;
    
    const NuclideRelAct *nuc_info = nullptr;
    for( const auto &rel_act : m_rel_activities )
      nuc_info = (rel_act.nuclide == nuc) ? &rel_act : nuc_info;
    
    assert( nuc_info );
    if( !nuc_info )
      continue;
    
    all_nuclides.insert( nuc );
    if( !nuclides_with_colors.count(nuc) && !peak.lineColor().isDefault() )
    {
      add_rel_eff_plot_css << "        .RelEffPlot circle." << nuc->symbol << "{ fill: " << peak.lineColor().cssText() << "; }\n";
      nuclides_with_colors.insert( nuc );
    }
    
    double yield = 0;
    for( const pair<double,double> &energy_br : nuc_info->gamma_energy_br )
    {
      if( fabs(energy_br.first - peak.gammaParticleEnergy()) < 1.0E-6 )
        yield += energy_br.second;
    }
      
    snprintf( buffer, sizeof(buffer), "{nuc: \"%s\", br: %1.6G, rel_act: %1.6G}",
                nuc_info->nuclide->symbol.c_str(), yield, nuc_info->rel_activity );

    const string isotopes_json = buffer;
    
    const double src_counts = nuc_info->rel_activity * yield * live_time;
    //if( src_counts < 1.0E-6 )
    //  continue;
    
    //if( peak.amplitude() < 1.0E-6 )
    //  continue;
    
    const double eff = peak.amplitude() / src_counts;
    
    //TODO: need to use proper uncertainty
    const double eff_uncert = 0; //peak.amplitudeUncert() / src_counts;
    
    snprintf( buffer, sizeof(buffer),
             "%s{energy: %.2f, counts: %1.7g, counts_uncert: %1.7g,"
             " eff: %1.6g, eff_uncert: %1.6g, nuc_info: [%s]}",
             (num_rel_eff_data_points ? ", " : ""),
             peak.mean(), peak.amplitude(), peak.amplitudeUncert(),
             eff, eff_uncert, isotopes_json.c_str() );
    
    //cout << "peak rel_eff: " << buffer << endl;
    
    if( IsNan(eff) || IsInf(eff)
       || IsNan(eff_uncert) || IsInf(eff_uncert)
       || IsNan(peak.amplitude()) || IsInf(peak.amplitude())
       || IsNan(peak.amplitudeUncert()) || IsInf(peak.amplitudeUncert())
       )
    {
      cerr << "Skipping Rel Eff plot point: " << buffer << endl;
    }else
    {
      rel_eff_plot_values << buffer;
      num_rel_eff_data_points += 1;
    }
  }//for( const PeakDef &peak : m_fit_peaks )
  
  rel_eff_plot_values << "]";
  
  
  // Incase user didnt associate color with a nuclide, assign some (pretty much randomly chosen)
  //  default colors.
  size_t unseen_nuc_index = 0;
  const vector<string> default_nuc_colors{ "#003f5c", "#ffa600", "#7a5195", "#ff764a", "#ef5675", "#374c80" };
  for( const auto *nuc : all_nuclides )
  {
    if( nuclides_with_colors.count(nuc) )
      continue;
    
    add_rel_eff_plot_css << "        .RelEffPlot circle." << nuc->symbol << "{ fill: "
                     << default_nuc_colors[unseen_nuc_index % default_nuc_colors.size()]
                     << "; }\n";
    
    unseen_nuc_index += 1;
  }//for( const auto *nuc : all_nuclides )
  

  auto load_file_contents = []( string filename ) -> string {
    string filepath = wApp ? wApp->docRoot() : "InterSpec_resources";
    filepath = SpecUtils::append_path(filepath, filename );
    
    vector<char> file_data;
    try
    {
      SpecUtils::load_file_data( filepath.c_str(), file_data );
    }catch( std::exception & )
    {
      throw std::runtime_error( "Failed to read " + filename );
    }
    
    return string( begin( file_data ), end( file_data ) );
  };//load_file_contents(...)
  
  
  string html = load_file_contents( "static_text/auto_rel_act_report.tmplt.html" );
  const string d3_js = load_file_contents( "d3.v3.min.js" );
  const string spectrum_chart_d3_js = load_file_contents( "SpectrumChartD3.js" );
  const string spectrum_chart_d3_css = load_file_contents( "SpectrumChartD3.css" );
  //const string spectrum_chart_d3_standalone_css = load_file_contents( "SpectrumChartD3StandAlone.css" );
  
  
  SpecUtils::ireplace_all( html, "\\;", ";" );
  
  SpecUtils::ireplace_all( html, "${D3_SCRIPT}", d3_js.c_str() );
  SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_JS}", spectrum_chart_d3_js.c_str() );
  SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_CSS}", spectrum_chart_d3_css.c_str() );
  
  string html_title = m_options.spectrum_title;
  SpecUtils::ireplace_all( html, "${TITLE}", html_title.c_str() );
  
  const string rel_eff_plot_js = load_file_contents( "RelEffPlot.js" );
  const string rel_eff_plot_css = load_file_contents( "RelEffPlot.css" );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_JS}", rel_eff_plot_js.c_str() );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_CSS}", rel_eff_plot_css.c_str() );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_ADDITIONAL_CSS}", add_rel_eff_plot_css.str().c_str() );
  
  SpecUtils::ireplace_all( html, "${REL_EFF_DATA_VALS}", rel_eff_plot_values.str().c_str() );
  
  const string rel_eff_eqn_js = RelActCalc::rel_eff_eqn_js_function( m_rel_eff_form, m_rel_eff_coefficients );
  SpecUtils::ireplace_all( html, "${FIT_REL_EFF_EQUATION}", rel_eff_eqn_js.c_str() );
  
  SpecUtils::ireplace_all( html, "${RESULTS_TXT}", results_html.str().c_str() );
  
  
  if( m_spectrum )
  {
    stringstream set_js_str;
        
    vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    D3SpectrumExport::write_js_for_chart( set_js_str, "specchart", "", "Energy (keV)", "Counts"  );
    
    set_js_str <<
    "  let spec_observer = new ResizeObserver(entries => {\n"
    "    for (let entry of entries) {\n"
    "      if (entry.target && (entry.target.id === \"specchart\")) {\n"
    "        spec_chart_specchart.handleResize(false);\n"
    "      }\n"
    "    }\n"
    "  });\n"
    "  spec_observer.observe( document.getElementById(\"specchart\") );\n"
    ;
    
    
    D3SpectrumExport::D3SpectrumChartOptions chart_options;
    chart_options.m_useLogYAxis = true;
    chart_options.m_legendEnabled = false;
    chart_options.m_compactXAxis = true;
    chart_options.m_allowDragRoiExtent = false;
    write_set_options_for_chart( set_js_str, "specchart", chart_options );
    set_js_str << "  spec_chart_specchart.setShowLegend(false);\n";
    
    // We'll set the initial display energy to be just the parts of the spectrum analyzed, plus a
    //  little padding
    float lower_energy = m_spectrum->gamma_energy_max();
    float upper_energy = m_spectrum->gamma_energy_min();
    for( const RoiRange &rr : m_input_roi_ranges )
    {
      lower_energy = std::min( lower_energy, static_cast<float>(rr.lower_energy) );
      upper_energy = std::max( upper_energy, static_cast<float>(rr.upper_energy) );
    }
    
    const float range = upper_energy - lower_energy;
    
    lower_energy -= 0.1f * range;
    upper_energy += 0.1f * range;
    
    lower_energy = std::max( m_spectrum->gamma_energy_min(), lower_energy );
    upper_energy = std::min( m_spectrum->gamma_energy_max(), upper_energy );
    
    if( upper_energy <= lower_energy )
    {
      lower_energy = m_spectrum->gamma_energy_min();
      upper_energy = m_spectrum->gamma_energy_max();
    }
    
    set_js_str << "  spec_chart_specchart.setXAxisRange("
               << lower_energy << ", " << upper_energy << ", false);\n";
        
    D3SpectrumExport::D3SpectrumOptions spec_options;
    spec_options.spectrum_type = SpecUtils::SpectrumType::Foreground;
    
    vector<shared_ptr<const PeakDef>> peaks;
    for( const auto &p : m_fit_peaks )
      peaks.push_back( make_shared<PeakDef>(p) );
    spec_options.peaks_json = PeakDef::peak_json( peaks, m_spectrum );
    
    const SpecUtils::Measurement * const meas_ptr = m_spectrum.get();
    D3SpectrumExport::write_and_set_data_for_chart( set_js_str, "specchart", { std::make_pair(meas_ptr,spec_options) } );
        
    
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_DIV}",
                            "<div id=\"specchart\" style=\"height: 30vw; flex: 1 2; overflow: hidden;\" class=\"SpecChart\"></div>" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_INIT_JS}", set_js_str.str().c_str() );
    
    
    SpecUtils::ireplace_all( html, "${CHART_SPACER_LEFT}", "" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_RIGHT}", "" );
  }else
  {
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_DIV}", "" );
    //SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_JS}", "" );
    //SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_CSS}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_INIT_JS}", "" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_LEFT}", "<div style=\"width: 10%\"> </div>" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_RIGHT}", "<div style=\"width: 15%\"> </div>" );
  }//if( m_spectrum ) / else
  
  out << html;

}//void print_html_report( std::ostream &out ) const



double RelActAutoSolution::mass_enrichment_fraction( const SandiaDecay::Nuclide *nuclide ) const
{
  const size_t nuc_index = nuclide_index( nuclide );
  const double rel_mass = m_rel_activities[nuc_index].rel_activity / nuclide->activityPerGram();
  
  double el_total_mass = 0.0;
  for( const NuclideRelAct &nuc : m_rel_activities )
  {
    if( nuc.nuclide->atomicNumber == nuclide->atomicNumber )
      el_total_mass += nuc.rel_activity / nuc.nuclide->activityPerGram();
  }
  
  // TODO: - For Pu or U, need to account for Pu242 and/or U236 by using something like by correlation
  //         Some potential papers to use to do this are:
  //         - Determination of 242Pu by correlation with 239Pu only, M.T.Swinhoe, T.Iwamoto, T.Tamura
  //           https://www.sciencedirect.com/science/article/pii/S0168900210000045
  //           (uses only U239 and has R.M.S deviation of 1.2% between 55 and 64% Pu239)
  //           Pu242 = 9.66E−3 * pow(Pu239,−3.83)  //I think - need to check if its mass fractions, or percentages, or relative activities, or what
  //           This paper also gives the values FRAM uses
  //        - Isotopic correlation for 242Pu composition prediction: Multivariate regression approach, Arnab Sarkar, et al
  //          https://www.sciencedirect.com/science/article/pii/S0969804314003820?via%3Dihub
  //          65.9 to 76% Pu239 - although actual values of coefficients dont seem to be given (???), and maybe only applicable to similar things as the samples used
  //        - Gunnink, R., 1982 (info in FRAM manual I think)
  //        - Bignan, G., et all 1998, ESARDA Bull. 28, 1-9
  //          https://esarda.jrc.ec.europa.eu/esarda-bulletin-n28_en#files
  //        - LA14018FRAMApplicationGuide.pdf
  
  
  return rel_mass / el_total_mass;
}//mass_enrichment_fraction


double RelActAutoSolution::mass_ratio( const SandiaDecay::Nuclide *numerator,
                                      const SandiaDecay::Nuclide *denominator ) const
{
  const double ratio = activity_ratio(numerator,denominator);
  
  return ratio * denominator->activityPerGram() / numerator->activityPerGram();
}//double mass_ratio(...)


double RelActAutoSolution::activity_ratio( const SandiaDecay::Nuclide *numerator,
                                          const SandiaDecay::Nuclide *denominator ) const
{
  const size_t i = nuclide_index(numerator);
  const size_t j = nuclide_index(denominator);
  
  const NuclideRelAct &nuc_i = m_rel_activities[i];
  const NuclideRelAct &nuc_j = m_rel_activities[j];
  
  const double act_i = nuc_i.rel_activity;
  const double act_j = nuc_j.rel_activity;
  
  return act_i / act_j;
}//double activity_ratio(...)


size_t RelActAutoSolution::nuclide_index( const SandiaDecay::Nuclide *nuclide ) const
{
  for( size_t i = 0; i < m_rel_activities.size(); ++i )
    if( nuclide == m_rel_activities[i].nuclide )
      return i;
  
  assert( 0 );
  throw runtime_error( "RelActAutoSolution: " + (nuclide ? nuclide->symbol : string()) + " is an invalid nuclide." );
  return m_rel_activities.size();
}//nuclide_index(...)

RelActAutoSolution solve( Options options,
                         std::vector<RoiRange> energy_ranges,
                         std::vector<NucInputInfo> nuclides,
                         std::vector<FloatingPeak> extra_peaks,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         std::shared_ptr<const DetectorPeakResponse> input_drf,
                         std::vector<std::shared_ptr<const PeakDef>> all_peaks
                         )
{
  return RelActAutoCostFcn::solve_ceres(
                     options,
                     energy_ranges,
                     nuclides,
                     extra_peaks,
                     foreground,
                     background,
                     input_drf,
                     all_peaks );
}//RelActAutoSolution


}//namespace RelActCalcAuto
