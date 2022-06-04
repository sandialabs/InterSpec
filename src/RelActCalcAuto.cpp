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
#include <limits>
#include <sstream>
#include <functional>
#include <type_traits>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "ceres/ceres.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{

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


void append_string_node( rapidxml::xml_node<char> *base_node, const char *node_name, const string &value )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  const char *strvalue = doc->allocate_string( value.c_str(), value.size() + 1 );
  xml_node<char> *node = doc->allocate_node( node_element, node_name, strvalue );
  base_node->append_node( node );
}//append_string_node(...)


void append_version_attrib( rapidxml::xml_node<char> *base_node, const int version )
{
  using namespace rapidxml;
  
  assert( base_node && base_node->document() );
  xml_document<char> *doc = base_node->document();
  
  char buffer[128];
  snprintf( buffer, sizeof(buffer), "%i", version );
  const char *value = doc->allocate_string( buffer );
  xml_attribute<char> *attrib = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attrib );
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
  vector<SandiaDecay::EnergyRatePair> nominal_gammas;
  
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
    mix.addAgedNuclideByActivity( nuclide, SandiaDecay::Bq, nominal_age );
    nominal_gammas = mix.gammas( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
    
    for( double energy : gammas_to_exclude )
    {
      const size_t did_remove = remove_gamma( energy, nominal_gammas );
      assert( did_remove );
    }
  }//NucInputGamma constructor
};//struct NucInputGamma



struct RelActAutoCostFcn /* : ROOT::Minuit2::FCNBase() */
{
  RelActCalcAuto::Options m_options;
  std::vector<NucInputGamma> m_nuclides;
  std::vector<RoiRangeChannels> m_energy_ranges;
  std::vector<RelActCalcAuto::FloatingPeak> m_extra_peaks;
  
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;
  std::shared_ptr<const SpecUtils::Measurement> m_background;
  std::shared_ptr<SpecUtils::Measurement> m_spectrum;
  
  float m_live_time;
  vector<float> m_channel_counts; //background subtracted channel counts
  vector<float> m_channel_count_uncerts; //e.g. sqrt( m_channel_counts[i] )
  std::shared_ptr<const SpecUtils::EnergyCalibration> m_energy_cal;
  
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
                     shared_ptr<const SpecUtils::Measurement> foreground,
                     shared_ptr<const SpecUtils::Measurement> background,
                     std::shared_ptr<const DetectorPeakResponse> drf
                           )
  : m_options( options ),
  m_nuclides{},
  m_energy_ranges{},
  m_extra_peaks( extra_peaks ),
  m_foreground( foreground ),
  m_background( background ),
  m_live_time( foreground ? foreground->live_time() : 0.0f ),
  m_channel_counts{},
  m_channel_count_uncerts{},
  m_energy_cal{},
  m_rel_eff_anchor_enhancement( 1000.0 ),
  m_drf( nullptr ),
  m_ncalls( 0 )
  {
    if( !foreground
       || (foreground->num_gamma_channels() < 128)
       || !foreground->energy_calibration()
       || !foreground->energy_calibration()->valid() )
      throw runtime_error( "RelActAutoCostFcn: invalid foreground spectrum." );
    
    if( background && (background->num_gamma_channels() != foreground->num_gamma_channels()) )
      throw runtime_error( "RelActAutoCostFcn: Diff num background/foreground channels." );
    
    m_energy_cal = foreground->energy_calibration();
    
    assert( foreground->gamma_counts() );
    m_channel_counts = *foreground->gamma_counts();
    m_channel_count_uncerts.resize( m_channel_counts.size(), 0.0 );
    
    if( !background )
    {
      for( size_t i = 0; i < m_channel_counts.size(); ++i )
      {
        const double counts = m_channel_counts[i];
        m_channel_count_uncerts[i] = (counts <= 1.0) ? 1.0 : sqrt( counts );
      }
      
      m_spectrum = make_shared<SpecUtils::Measurement>( *foreground );
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
      
      assert( background_counts.size() == m_channel_counts.size() );
      for( size_t i = 0; i < m_channel_counts.size(); ++i )
      {
        const double fore_counts = (m_channel_counts[i] < 0.0f) ? 0.0 : m_channel_counts[i];
        const double back_counts = (background_counts[i] < 0.0f) ? 0.0f : background_counts[i];
        const double uncert_2 = fore_counts*fore_counts + lt_sf*lt_sf*back_counts*back_counts;
        const double sub_val = fore_counts - lt_sf*back_counts;
        
        m_channel_counts[i] = static_cast<float>( std::max(sub_val,0.0) );
        m_channel_count_uncerts[i] = static_cast<float>( std::max( 1.0, sqrt(uncert_2) ) );
      }//for( loop over and set channel counts and uncertainties )
      
      m_spectrum = make_shared<SpecUtils::Measurement>( *foreground );
      m_spectrum->set_gamma_counts( make_shared<vector<float>>(m_channel_counts), foreground->live_time(), foreground->real_time() );
    }//if( !background ) / else
    
    
    //Need to initialize m_nuclides
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
                                               m_foreground->gamma_energy_min(),
                                               m_foreground->gamma_energy_max() );
        
        vector<shared_ptr<const PeakDef> > peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, nullptr, {}, false );
          
        auto all_peaks = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks), end(peaks) );
        
        new_drf->fitResolution( all_peaks, foreground, DetectorPeakResponse::kGadrasResolutionFcn );
        
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
    
    const bool highres = PeakFitUtils::is_high_res( foreground );
    
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
      const double min_br = numeric_limits<double>::epsilon();  //arbitrary
      const double num_fwhm_roi = 2.5; // arbitrary...
      
      vector<pair<double,double>> gammas_in_range;
      for( const auto &n : m_nuclides )
      {
        for( const auto &g : n.nominal_gammas )
        {
          if( (g.numPerSecond > min_br)
             && (g.energy >= roi_range.lower_energy)
             && (g.energy <= roi_range.upper_energy) )
          {
            double energy_sigma;
            float min_sigma, max_sigma;
            expected_peak_width_limits( g.energy, highres, min_sigma, max_sigma );
            
            if( m_drf && m_drf->hasResolutionInfo() )
            {
              energy_sigma = m_drf->peakResolutionSigma(g.energy);
              
              // A sanity check... maybe we dont want this?
              if( energy_sigma < min_sigma )
                energy_sigma = min_sigma;
              if( energy_sigma > max_sigma )
                energy_sigma = max_sigma;
            }else
            {
              energy_sigma = max_sigma;
            }
            
            double gamma_row_lower = g.energy - num_fwhm_roi*energy_sigma;
            double gamma_row_upper = g.energy + num_fwhm_roi*energy_sigma;
            
            if( !roi_range.allow_expand_for_peak_width )
            {
              gamma_row_lower = std::max( gamma_row_lower, roi_range.lower_energy );
              gamma_row_upper = std::min( gamma_row_upper, roi_range.upper_energy );
            }
            
            gammas_in_range.push_back( {gamma_row_lower,gamma_row_upper} );
          }
        }//for( const auto &g : n.nominal_gammas )
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
    
    // TODO: Figure out a proper value to set m_rel_eff_anchor_enhancement to, like maybe largest peak area divided m_live_time, or something like that
    m_rel_eff_anchor_enhancement = 0.0;
    for( const auto &r : m_energy_ranges )
      m_rel_eff_anchor_enhancement += m_spectrum->gamma_integral( r.lower_energy, r.upper_energy );
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
                                                        const std::shared_ptr<const DetectorPeakResponse> input_drf
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
    solution.m_rel_eff_form     = options.rel_eff_eqn_type;
    solution.m_fwhm_form        = options.fwhm_form;
    solution.m_input_roi_ranges = energy_ranges;
    solution.m_fit_energy_cal_adjustments = options.fit_energy_cal;
    if( input_drf && input_drf->isValid() && input_drf->hasResolutionInfo() )
      solution.m_drf = input_drf;
    
    RelActAutoCostFcn *cost_functor = nullptr;
    try
    {
      cost_functor = new RelActAutoCostFcn( options, energy_ranges, nuclides,
                                       extra_peaks, foreground, background, input_drf );
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
    solution.m_spectrum = cost_functor->m_spectrum;
    
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
        problem.SetParameterLowerBound(pars + 0, 0, -5.0 );
        problem.SetParameterUpperBound(pars + 0, 0, +5.0 );
      }else
      {
        //We'll only fit gain
        problem.SetParameterBlockConstant( pars + 0 );
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
    const size_t acts_start = rel_eff_start + options.rel_eff_eqn_order + 1;
    const size_t free_peak_start = acts_start + 2*nuclides.size();
    assert( (free_peak_start + 2*extra_peaks.size()) == cost_functor->number_parameters() );
    
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
        assert( !haveGadras || (drfpars.size() == 3) );
        const auto formToFit = haveGadras
                                ? DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial
                                : DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
        
        vector<float> sigma_coefs, sigma_coef_uncerts;
        MakeDrfFit::performResolutionFit( fake_peaks, formToFit,
                                       highres, num_fwhm_pars, sigma_coefs, sigma_coef_uncerts );
      
        assert( sigma_coefs.size() == num_fwhm_pars );
        assert( haveGadras || (sigma_coefs.size() == 3) );
        if( sigma_coefs.size() != num_fwhm_pars )
        
        drfpars = sigma_coefs;
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
    //  (e.g., we'll start with rel eff line == 1.0 for all energies)
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
              throw runtime_error( "When fitting the age of an element, all nuclides of the element"
                                   " must have the same age, and same value of to be fit or not." );
            
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
      mix.addAgedNuclideByActivity( nuc.nuclide, 1.0, nuc.age );
      const vector<SandiaDecay::EnergyRatePair> gammas = mix.photons(0.0);
      for( const auto &erange : energy_ranges )
      {
        vector<tuple<double,double,double>> energy_to_rel_act; //{energy, br, rel. act.}
        for( const SandiaDecay::EnergyRatePair &er : gammas )
        {
          if( er.numPerSecond <= std::numeric_limits<float>::min() )
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
          
          // The FWHM area covers about 76% of gaussian area
          const double rel_act = data_count / er.numPerSecond / foreground->live_time() / 0.76;
          
          energy_to_rel_act.emplace_back( er.energy, er.numPerSecond, rel_act );
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
    ceres_options.num_threads = 4; // TODO: SpecUtils::num_logical_cpu_cores() or num_physical_cpu_cores()
    
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_options, &problem, &summary);
    //std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
    
    solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
    
    
    
    
    solution.m_num_function_eval = cost_functor->m_ncalls;
    
    // blah blah blah - need to fill out all the answer info - also need to update energy calibration of m_spectrum, as well as shift the
    
    /*
    std::vector<double> m_rel_eff_coefficients;
    std::vector<std::vector<double>> m_rel_eff_covariance;
    std::vector<NuclideRelAct> m_rel_activities;
    std::vector<std::vector<double>> m_rel_act_covariance;
    std::vector<double> m_fwhm_coefficients;
    std::vector<std::vector<double>> m_fwhm_covariance;
    std::vector<FloatingPeakResult> m_floating_peaks;
    std::vector<PeakDef> m_fit_peaks;
    double m_energy_cal_adjustments[2];
    bool m_fit_energy_cal_adjustments;
    double m_chi2;
    size_t m_dof;
     // Update solution.m_spectrum with updated energy calibration
    */
    
    
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
      
    // TODO: is this the best way to apply corrections?
    return (1.0 + x[1]) * (energy + x[0]);
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
  
  
  PeaksForEnergyRange peaks_for_energy_range( const RoiRangeChannels &range, const std::vector<double> &x ) const
  {
    const size_t num_channels = range.num_channels;
    
    // We will use "adjusted" to refer to energies that have been mapped into the spectrums original
    //  energy calibrations
    const double adjusted_lower_energy = apply_energy_cal_adjustment( range.lower_energy, x );
    const double adjusted_upper_energy = apply_energy_cal_adjustment( range.upper_energy, x );
    
    pair<size_t,size_t> channel_range = range.channel_range( adjusted_lower_energy, adjusted_upper_energy, num_channels, m_energy_cal );
    const size_t first_channel = channel_range.first;
    const size_t last_channel = channel_range.first;
    
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
      vector<SandiaDecay::EnergyRatePair> gammas;
      
      const double rel_act = relative_activity( nucinfo.nuclide, x );
      
      cout << "peaks_for_energy_range: Relative activity of " << nucinfo.nuclide->symbol << " is " << PhysicalUnits::printToBestActivityUnits(rel_act) << endl;
      
      if( is_fixed_age(nucinfo.nuclide) )
      {
        gammas = nucinfo.nominal_gammas;
      }else
      {
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nucinfo.nuclide, SandiaDecay::Bq, age(nucinfo.nuclide,x) );
        gammas = mix.gammas( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
        
        for( double energy : nucinfo.gammas_to_exclude )
        {
          const size_t did_remove = nucinfo.remove_gamma( energy, gammas );
          assert( did_remove );
        }
      }//if( age is fixed ) / else( age may vary )
      
      for( const SandiaDecay::EnergyRatePair &gamma : gammas )
      {
        // Filter out zero-amplitude gammas
        // TODO: - come up with more intelligent lower bound of gamma rate to bother with
        if( gamma.numPerSecond < std::numeric_limits<float>::epsilon() )
          continue;
        
        // Filter the energy range based on true energies (i.e., not adjusted energies)
        if( (gamma.energy < range.lower_energy) || (gamma.energy > range.upper_energy) )
          continue;
        
        nuclides_used.insert( nucinfo.nuclide );
        
        // We compute the relative efficiency and FWHM based off of "true" energy
        const double rel_eff = relative_eff( gamma.energy, x );
        const double peak_fwhm = fwhm( gamma.energy, x );
        
        const double peak_amplitude = rel_act * m_live_time * rel_eff;
        const double peak_mean = apply_energy_cal_adjustment( gamma.energy, x );
        
        cout << "peaks_for_energy_range: Peak at " << gamma.energy << " keV for " << nucinfo.nuclide->symbol
        << " has a mean of " << peak_mean << " keV, FWHM=" << peak_fwhm << ", RelEff=" << rel_eff
        << ", and AMP=" << peak_amplitude << endl;
        
        peaks.emplace_back( peak_mean, peak_fwhm/2.35482, peak_amplitude );
        
        PeakDef &new_peak = peaks.back();

        // TODO: (low priority) Instead of calling #PeakDef::findNearestPhotopeak, should instead keep gamma source information from decay... and should actually implement this in SandiaDecay...
        const bool xray_only = false;
        size_t transition_index = 0;
        const SandiaDecay::Transition *transition = nullptr;
        PeakDef::SourceGammaType gamma_type;
        PeakDef::findNearestPhotopeak( nucinfo.nuclide, gamma.energy, 0.0, xray_only,
                                       transition, transition_index, gamma_type );
        assert( transition );
        if( transition )
          new_peak.setNuclearTransition( nucinfo.nuclide, transition,
                                           static_cast<int>(transition_index), gamma_type );
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
      
      cout << "peaks_for_energy_range: free peak at " << peak.energy << " has a FWHM=" << peak_fwhm << " and AMP=" << peak_amp << endl;
      
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
    
    
    cerr << "peaks_for_energy_range: returning " << peaks.size() << " in range [" << range.lower_energy << ", "
    << range.upper_energy << "] keV." << endl;
    
    return answer;
  }//PeaksForEnergyRange peaks_for_energy_range( const double lower_energy, const double upper_energy ) const
  
  
  
  void eval( const std::vector<double> &x, double *residuals ) const
  {
    m_ncalls += 1;
    
    assert( x.size() == number_parameters() );
    assert( residuals );
    assert( !m_energy_ranges.empty() );
    
    size_t residual_index = 0;
    for( const RoiRangeChannels &energy_range : m_energy_ranges )
    {
      const PeaksForEnergyRange info = peaks_for_energy_range( energy_range, x );
      
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
        
        cout << "eval: For channel " << channel << " (" << lower_energy << " to " << upper_energy
        << " keV) there are " << data_counts << " data counts, and currently fitting for "
        << gaussian_area << " gaussian counts, and " << continuum_counts
        << " continuum counts (total fit: " << (gaussian_area + continuum_counts) << ")" << endl;
        
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
    const char *xml_file_path = "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/simple_u_test.xml";
    
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
    
    const auto options_node = get_required_node(base_node, "Options");
    Options options;
    options.fromXml( options_node );

    std::vector<RoiRange> energy_ranges;
    const auto roi_ranges_node = get_required_node(base_node, "RoiRangeList");
    XML_FOREACH_DAUGHTER(roi_range_node, roi_ranges_node, "RoiRange")
    {
      RoiRange range;
      range.fromXml( roi_range_node );
      energy_ranges.push_back( range );
    }
    
    if( energy_ranges.empty() )
      throw runtime_error( "No RoiRanges specified" );
    
    std::vector<NucInputInfo> nuclides;
    const auto nucs_node = get_required_node(base_node, "NucInputInfoList");
    XML_FOREACH_DAUGHTER(nuc_node, nucs_node, "NucInputInfo")
    {
      NucInputInfo nuc;
      nuc.fromXml( nuc_node );
      nuclides.push_back( nuc );
    }
    
    if( nuclides.empty() )
      throw runtime_error( "No nuclides specified" );
    
    
    std::vector<FloatingPeak> extra_peaks;
    const auto float_peak_node = get_required_node(base_node, "FloatingPeakList");
    XML_FOREACH_DAUGHTER(float_peak_node, nucs_node, "FloatingPeak")
    {
      FloatingPeak peak;
      peak.fromXml( float_peak_node );
      extra_peaks.push_back( peak );
    }
    
    
    assert( foreground && (foreground->num_measurements() == 1) );
    
    shared_ptr<const DetectorPeakResponse> drf = foreground->detector();
    shared_ptr<const SpecUtils::Measurement> fore_meas = foreground->measurements()[0];
    shared_ptr<const SpecUtils::Measurement> back_meas;
    if( background )
      back_meas = background->measurements()[0];

    RelActAutoSolution solution = solve( options, energy_ranges, nuclides, extra_peaks, fore_meas, back_meas, drf );
    
    //Print out summarry, etc.
  }catch( std::exception &e )
  {
    cerr << "RelAct test failed: " << e.what() << endl;
    return EXIT_FAILURE;
  }// try / catch
  
  return EXIT_SUCCESS;
  
  /*
   Create input from XML:
   

   Then create XML to represent Options, FloatingPeak, NucInputInfo, RoiRange...
   
   Then create something to printout summary of RelActAutoSolution ....  Then wait on making sure things actually work before going much further
   

   
   
   
   
   <Nuclides><Nuclide>U235</Nuclide><Nuclide>U238</Nuclide>...</Nuclides>
   <EnergyRanges><Range><Lower
   */
  
  
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
    age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( age_str, nuclide->halfLife );
    
    fit_age = get_bool_node_value( nuc_info_node, "FitAge" );
    
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
  append_string_node( base_node, "RelEffEqnType", rell_eff_eqn_str);
  
  append_int_node( base_node, "RelEffEqnOrder", rel_eff_eqn_order);
  
  const char *fwhm_form_str = to_str( fwhm_form );
  append_string_node( base_node, "FwhmForm", fwhm_form_str );
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
  fwhm_form( FwhmForm::Polynomial_2 )
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
  m_fit_energy_cal_adjustments( false ),
  m_chi2( -1.0 ),
  m_dof( 0 ),
  m_num_function_eval( 0 ),
  m_num_microseconds_eval( 0 )
{
  
}



RelActAutoSolution solve( Options options,
                         std::vector<RoiRange> energy_ranges,
                         std::vector<NucInputInfo> nuclides,
                         std::vector<FloatingPeak> extra_peaks,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         std::shared_ptr<const DetectorPeakResponse> input_drf
                         )
{
  return RelActAutoCostFcn::solve_ceres(
                     options,
                     energy_ranges,
                     nuclides,
                     extra_peaks,
                     foreground,
                     background,
                     input_drf );
}//RelActAutoSolution


}//namespace RelActCalcAuto
