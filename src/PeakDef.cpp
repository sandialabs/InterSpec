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

#include <regex>
#include <memory>
#include <iostream>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/poisson.hpp>

#include <Wt/Json/Value>
#include <Wt/Json/Object>
#include <Wt/WApplication>
#include <Wt/Json/Serializer>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace std;
using SpecUtils::Measurement;

/** Peak XML minor version  changes:
 20231031: Minor version incremented from 1 to 2 to account for new peak skew models, and removal of Landau skew model
 */
const int PeakDef::sm_xmlSerializationMajorVersion = 0;
const int PeakDef::sm_xmlSerializationMinorVersion = 2;

const bool PeakDef::sm_defaultUseForDrfIntrinsicEffFit = true;
const bool PeakDef::sm_defaultUseForDrfFwhmFit = true;
const bool PeakDef::sm_defaultUseForDrfDepthOfInteractionFit = false;

/** Version 1 adds "FlatStep", "LinearStep", and "BiLinearStep" continuum types
 
 */
const int PeakContinuum::sm_xmlSerializationVersion = 1;

namespace
{
  //clones 'source' into the document that 'result' is a part of.
  //  'result' is cleared and set lexically equal to 'source'.
  void clone_node_deep( const ::rapidxml::xml_node<char> *source,
                        ::rapidxml::xml_node<char> *result )
  {
    using namespace ::rapidxml;
    
    xml_document<char> *doc = result->document();
    if( !doc )
      throw runtime_error( "clone_node_deep: insert result into document before calling" );
    
    result->remove_all_attributes();
    result->remove_all_nodes();
    result->type(source->type());
    
    
    // Clone name and value
    char *str = doc->allocate_string( source->name(), source->name_size() );
    result->name( str, source->name_size());
    
    if( source->value() )
    {
      str = doc->allocate_string( source->value(), source->value_size() );
      result->value( str, source->value_size() );
    }
    
    // Clone child nodes and attributes
    for( xml_node<char> *child = source->first_node(); child; child = child->next_sibling() )
    {
      xml_node<char> *clone = doc->allocate_node( child->type() );
      result->append_node( clone );
      clone_node_deep( child, clone );
    }
    
    for( xml_attribute<char> *attr = source->first_attribute(); attr; attr = attr->next_attribute())
    {
      const char *name = doc->allocate_string( attr->name(), attr->name_size() );
      const char *value = doc->allocate_string( attr->value(), attr->value_size() );
      xml_attribute<char> *clone = doc->allocate_attribute(name, value, attr->name_size(), attr->value_size());
      result->append_attribute( clone );
    }
  }//void clone_node_deep(...)
  
  
  bool gives_off_gammas( const SandiaDecay::Nuclide *nuc )
  {
    if( !nuc )
      return false;
      
    for( const SandiaDecay::Transition *t : nuc->decaysToChildren )
    {
      for( const SandiaDecay::RadParticle &p : t->products )
      {
        if( p.type == SandiaDecay::GammaParticle )
          return true;
      }
        
      if( gives_off_gammas( t->child ) )
        return true;
    }//for( size_t decN = 0; decN < decaysToChildren.size(); ++decN )
      
    return false;
  }//bool gives_off_gammas( const SandiaDecay::Nuclide *nuc )
}//namespace



void findROIEnergyLimits( double &lowerEnengy, double &upperEnergy,
                         const PeakDef &peak, const std::shared_ptr<const Measurement> &data )
{
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  if( continuum->energyRangeDefined() )
  {
    lowerEnengy  = continuum->lowerEnergy();
    upperEnergy = continuum->upperEnergy();
    return;
  }//if( continuum->energyRangeDefined() )
  
  if( !data || (data->num_gamma_channels() < 2) )
  {
    lowerEnengy = peak.lowerX();
    upperEnergy = peak.upperX();
    return;
  }//if( !data )
  
  const size_t lowbin = findROILimit( peak, data, false );
  const size_t upbin  = findROILimit( peak, data, true );
  if( lowbin == 0 )
    lowerEnengy = data->gamma_channel_center( lowbin );
  else
    lowerEnengy = data->gamma_channel_lower( lowbin );
  
  if( (upbin+1) >= data->num_gamma_channels() )
    upperEnergy = data->gamma_channel_center( std::min(upbin,data->num_gamma_channels()-1) );
  else
    upperEnergy = data->gamma_channel_upper( upbin );
  

}//void findROIEnergyLimits(...)


#define PRINT_ROI_DEBUG_INFO 0

#if( PRINT_ROI_DEBUG_INFO )
namespace
{
  class DebugLog
  {
    //Make it so cout/cerr statments always end up non-interleaved when multiple
    // threads are calling cout/cerr
  public:
    explicit DebugLog( std::ostream &os ) : os(os) {}
    ~DebugLog() { os << ss.rdbuf() << std::flush; }
    template <typename T>
    DebugLog& operator<<(T const &t){ ss << t; return *this;}
  private:
    std::ostream &os;
    std::stringstream ss;
  };
}//namespace
#endif



size_t findROILimitHighRes( const PeakDef &peak, const std::shared_ptr<const Measurement> &dataH, bool high )
{
  if( !dataH || !dataH->energy_calibration() || !dataH->energy_calibration()->valid() )
    return 0;
  
  const float mean = static_cast<float>( peak.mean() );
  const float fwhm = static_cast<float>( peak.fwhm() );
  
  // The plan is to iterate over the V&V test files and adjust these next quantities to most closely
  //  match the human selected ROI widths.
  //  Also, need to investigate highres_shrink_roi(...) to work better, for example the Cs137 peak of a shielded specturm
  const double min_width_mult = 0.75;
  const float nominal_fwhm_mult = 2.0f;
  const float steep_fwhm_mult = 3.5f;
  const float steep_continuum_limit = 0.35f;
  const float feature_nsigma_limit = 2.25f;
  
  
  const int direction = high ? 1 : -1;
  float nominal_low_edge = mean - (nominal_fwhm_mult * fwhm);
  size_t nominal_low_channel = dataH->find_gamma_channel( nominal_low_edge );
  
  float nominal_up_edge = mean + (nominal_fwhm_mult * fwhm);
  size_t nominal_up_channel = dataH->find_gamma_channel( nominal_up_edge );
  
  double nominal_data_area = dataH->gamma_channels_sum( nominal_low_channel, nominal_up_channel );
  
  double coefficients[2];
  const size_t nSideChanel = 3;
  PeakContinuum::eqn_from_offsets( nominal_low_channel, nominal_up_channel,
                                   mean, dataH, nSideChanel, nSideChanel, coefficients[1], coefficients[0] );
  
  // cout << "coefficients[0]=" << coefficients[0] << ", coefficients[1]=" << coefficients[1]
  //  << "\n\tdataH->gamma_channel_lower(nominal_low_channel)=" << dataH->gamma_channel_lower(nominal_low_channel)
  //  << "\n\tdataH->gamma_channel_upper(nominal_up_channel)=" << dataH->gamma_channel_upper(nominal_up_channel)
  //  << endl << endl;
  
  double nominal_cont_area = PeakContinuum::offset_eqn_integral( coefficients,
                                                 PeakContinuum::OffsetType::Linear,
                                                 dataH->gamma_channel_lower(nominal_low_channel),
                                                 dataH->gamma_channel_upper(nominal_up_channel),
                                                 mean );
  
  const double nominal_peak_area = (nominal_data_area - nominal_cont_area);
  const double nominal_uncert_ratio = nominal_peak_area / sqrt(nominal_data_area);
  
  //cout << "\n\n\nmean=" << mean << ", nominal_data_area=" << nominal_data_area << endl;
  //cout << "nominal_peak_area=" << nominal_peak_area << ", nominal_cont_area=" << nominal_cont_area
  //     << ", nominal_uncert_ratio=" << nominal_uncert_ratio << endl;
  
  const double width = dataH->gamma_channel_upper(nominal_up_channel) - dataH->gamma_channel_lower(nominal_low_channel);
  const double offset_contrib = coefficients[0]*width;
  const double slope_contrib = 0.5*coefficients[1]*width*width;
  
  //cout << "For " << mean << " kev, offset contributes " << offset_contrib << ", while slope contributes "
  //     << slope_contrib << "--->" << slope_contrib/offset_contrib << endl;
  
  if( (offset_contrib <= 0.0) || ((direction*slope_contrib/offset_contrib) > steep_continuum_limit) ) // 0.35 is arbitrary
  {
    if( high )
    {
      nominal_up_edge = mean + (steep_fwhm_mult * fwhm);
      nominal_up_channel = dataH->find_gamma_channel( nominal_up_edge );
    }else
    {
      nominal_low_edge = mean - (steep_fwhm_mult * fwhm);
      nominal_low_channel = dataH->find_gamma_channel( nominal_low_edge );
    }
  }//if( coefficients[1] < -30 )
  
  
  const size_t feature_detect_start = dataH->find_gamma_channel( mean + direction*min_width_mult*fwhm );
  const size_t feature_detect_stop = high ? nominal_up_channel : nominal_low_channel;

  
  // Feature detect start FWHM mult
  const size_t nchannel = dataH->num_gamma_channels();
  const size_t nprev_avrg = 4;
  for( size_t channel = feature_detect_start; channel != feature_detect_stop; channel += direction )
  {
    if( (channel <= nprev_avrg) || ((channel + nprev_avrg + 1) >= nchannel) )
      continue;
    
    const size_t prev_sum_start = channel - direction*nprev_avrg;
    const size_t prev_sum_end = channel - direction;
    
    const float prev_sum = dataH->gamma_channels_sum( prev_sum_start, prev_sum_end );
    const float val = dataH->gamma_channel_content( channel );
    const float val_next = dataH->gamma_channel_content( channel + direction );

    const float prev_avrg = prev_sum / nprev_avrg;
    const float prev_avrg_uncert = std::max( 1.0f*nprev_avrg, std::sqrt(prev_sum) ) / nprev_avrg;
    
    const float max_allowable = std::ceil(prev_avrg + feature_nsigma_limit*prev_avrg_uncert) +  0.001f;
    
    //cout << "mean=" << mean << ", channel=" << channel << ", prev_avrg(" << prev_sum_start
    //     <<  "," << prev_sum_end<< ")=" << prev_avrg
    //     << ", val=" << val << ", nextval=" << val_next << ", max_allowable=" << max_allowable << ""
    //     << endl;
    
    // We'll require two bins to be outside of tolerance
    if( (val > max_allowable && val_next > max_allowable) )
      return channel - direction;
  }//for( int bin = minBin; bin > lastbin; --bin )
  

  return high ? nominal_up_channel : nominal_low_channel;
}//findROILimitHighRes(...)
 

size_t findROILimit( const PeakDef &peak, const std::shared_ptr<const Measurement> &dataH, bool high )
{
  if( !peak.gausPeak() )
    return dataH->find_gamma_channel( (high ? peak.upperX() : peak.lowerX()) );
  

  //This implemntation is an adaptation of how PCGAP defines a region of
  //  interest.
  //  The basic idea is to include a maximum of 11.75 sigma away from the mean
  //  of the peak, but then start at ~1.5 sigma from mean, and try to detect if
  //  a new feature is occuring, and if so, stop the region of interest there.
  //  A feature is "detected" if the value of the bin contents exceeds 2.5 sigma
  //  from the "expected" value, where the expected value starts off being the
  //  smallest bin value so far (well, this bin averaged with the bins on either
  //  side of it), and then each preceeding bin is added to this background
  //  value.
  //
  //References for PCGAP are at:
  //  http://www.inl.gov/technicalpublications/Documents/3318133.pdf
  //  http://www.osti.gov/bridge/servlets/purl/800710-Zc4iYJ/native/800710.pdf
  
  typedef int indexing_t;
  //  typedef size_t indexing_t;
  
  const indexing_t nchannel = (!dataH ? indexing_t(0): static_cast<indexing_t>(dataH->num_gamma_channels()));
  
  if( nchannel<128 )
    throw runtime_error( "findROILimit(...): Invalid input" );
  
  const bool highres = PeakFitUtils::is_high_res(dataH);
  
  
  const vector<float> &contents = *dataH->gamma_channel_contents();
  
  
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  double lowxrange  = continuum->lowerEnergy();
  double highxrange = continuum->upperEnergy();
  
  const bool definedRange = (lowxrange != highxrange);
  
  if( definedRange )
    return dataH->find_gamma_channel( (high ? (highxrange-0.00001) : lowxrange) );
  
  if( PeakFitUtils::is_high_res(dataH) )
    return findROILimitHighRes( peak, dataH, high );
  
  
  const double mean = peak.mean();
  const double sigma = peak.sigma();
  
  const int direction = high ? 1 : -1;
  lowxrange = mean + direction * 7.5*sigma;  //2.3 FWHM
  
  const indexing_t nSideChannel = 1;
  
  indexing_t startchannel = static_cast<indexing_t>( dataH->find_gamma_channel( mean + direction*1.5*sigma ) );
  if( !high && startchannel < nSideChannel )
    startchannel = nSideChannel;
  
  indexing_t minChannel = startchannel;
  double minVal = contents[startchannel];
  indexing_t nbackbin = 1 + 2*nSideChannel;
  float backval = dataH->gamma_channels_sum( startchannel - nSideChannel,
                                             startchannel + nSideChannel );
  backval = std::max( backval, static_cast<float>(nbackbin) );
  
  const indexing_t meanchannel = static_cast<indexing_t>( dataH->find_gamma_channel( mean ) );
  indexing_t lastchannel  = static_cast<indexing_t>( dataH->find_gamma_channel( lowxrange + direction*0.0001 ) );
  
  //Make sure were not looking to far and that loop will terminate
  if( high ) lastchannel = std::min( lastchannel, nchannel-2 );
  else       lastchannel = std::max( lastchannel, indexing_t(1) );
  if( (high && (startchannel > lastchannel)) || (!high && (startchannel < lastchannel)) )
    startchannel = lastchannel;
  if( high || lastchannel>=direction )
    lastchannel += direction;
  
  
#if( PRINT_ROI_DEBUG_INFO )
  const bool debug_this_peak = (fabs(peak.mean() - 12.15) < 1.0); // && direction>0;
  const vector<float> &energies = *dataH->gamma_channel_energies();
  
  if( debug_this_peak )
    DebugLog(cerr) << "\n\n\n\nTo start with, lowxrange=" << lowxrange << " lastbin="
    << lastchannel << ", x(lastbin)=" << energies[lastchannel]
    << " for peak at mean=" << mean << ", sigma=" << sigma << "\n";
#endif
  
  //Lets find the bin with the smallest contents, and
  for( indexing_t channel = startchannel + direction*nSideChannel;
       channel != lastchannel && channel>nSideChannel; channel += direction )
  {
    assert( channel < (dataH->num_gamma_channels() + 100) );
    const float val = contents[channel];
    
    if( val <= minVal && (!highres || (contents[channel+direction] <= minVal)) )
    {
      minVal = val;
      minChannel = channel;
      nbackbin = 1 + 2*nSideChannel;
      backval = dataH->gamma_channels_sum( channel-nSideChannel,
                                           channel+nSideChannel );
      backval = std::max( backval, static_cast<float>(nbackbin) );
      channel += direction;
      if( channel == lastchannel )
        break;
#if( PRINT_ROI_DEBUG_INFO )
      if( debug_this_peak )
        DebugLog(cerr) << "FindMin: New background val at " << energies[channel]
        << ", backval/nbackbin=" << backval/nbackbin << "\n";
#endif
    }else
    {
      //If val is greater than backval, we _may_ be hitting a new feature, like
      //  a new peak, if we are, lets stop going any further away from the mean.
      //  Not detecting this may cause us to extend _past_ the new feature and
      //  find a global minimum we clearly dount want
      
      float background = backval / nbackbin;
      float background_sigma = sqrt(backval) / nbackbin;
      
      //For low resolution spectra, lets estimate the slope of the continuum, and
      //  use this to correct maximum_allowable number of counts; without doing this
      //  the lower ROI range for peaks on a falling continuum can be much to short.
      //  (not tested for HPGe)

      if( (nbackbin > 2) && !highres && channel > nbackbin )
      {
        try
        {
          const float *x = &(*dataH->channel_energies())[0];
          const float *y = &(*dataH->gamma_counts())[0];
          x = x + channel - direction - ((direction < 0) ? 0 : nbackbin);
          y = y + channel - direction - ((direction < 0) ? 0 : nbackbin);
          
          vector<double> coeffs, uncerts;
          fit_to_polynomial( x, y, nbackbin, 1, coeffs, uncerts );
          
          const float thisx = ((direction < 0) ? x[direction] : x[nbackbin+1]);
          background = std::max( coeffs[0] + thisx*coeffs[1], 1.0*nbackbin );
        }catch(...)
        {
        }
      }//if( nbackbin > 2 )

      const float sigma = sqrt( background_sigma*background_sigma + background );
      //      double max_allowable = std::ceil(background + 2.575829*sigma) + 0.001;
      float max_allowable = std::ceil(background + 3.0f*sigma) + 0.001f;
      
      if( background < 20.0f )
      {
        //boost::math::quantile shows up pretty significantly on the profiler,
        //  so we'll reduce the accuracy a bit.  I've arbitrarily selected a
        //  precision of 6 decimal places, even though this is probably more
        //  than we need (numerically weird stuff happens in places I dont
        //  understand, so erroring on the side of caution).
        using boost::math::policies::digits10;
        using boost::math::poisson_distribution;
        typedef boost::math::policies::policy<digits10<6> > my_pol_6;
        
        const poisson_distribution<float, my_pol_6 > pois( background );
        max_allowable = boost::math::quantile( pois, 0.99f );
      }//if( background < 20 )
      
      
      
      if( val > max_allowable && (!highres || contents[channel+direction] > max_allowable) )
      {
        //XXX - the below 3 is purely empircal, and meant to help avoid
        //      contamination due to the new feature
        if( channel >= 3*direction )
          lastchannel = channel - 3*direction;
        
#if( PRINT_ROI_DEBUG_INFO )
        if( debug_this_peak )
          DebugLog(cerr) << "FindMin: Found last bin at " << dataH->gamma_channel_center(lastchannel)
          << ", val=" << val << ", max_allowable=" << max_allowable << "\n";
#endif
        break;
      }else
      {
        ++nbackbin;
        backval += val;
#if( PRINT_ROI_DEBUG_INFO )
        if( debug_this_peak )
          DebugLog(cerr) << "FindMin: bin " << channel << ", x(bin)=" << energies[channel]
          << ", val=" << val << ", max_allowable=" << max_allowable
          << ", background=" << background << ", sigma=" << sigma << "\n";
#endif
      }
    }//if( val <= minVal ) / else
  }//for( ; bin > lastbin; --bin )
  
  nbackbin = 1 + 2*nSideChannel;
  backval = dataH->gamma_channels_sum( minChannel - nSideChannel,
                                       minChannel + nSideChannel );
  backval = std::max( backval, static_cast<float>(nbackbin) );
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "1) Bin with smallest contends at "
    << dataH->gamma_channel_center(minChannel)
    << ", backval=" << (backval/3.0) << ", lastbin=" << lastchannel
    << " x(lastbin)=" << dataH->gamma_channel_center(lastchannel) << "\n";
#endif
  
  //Make sure were not looking to far and that loop will terminate
  if( high ) lastchannel = std::min( lastchannel, nchannel-2 );
  else       lastchannel = std::max( lastchannel, indexing_t(1) );
  if( (high && (minChannel > lastchannel)) || (!high && (minChannel < lastchannel)) )
    minChannel = lastchannel;
  
  if( lastchannel || direction>0 )
    lastchannel += direction;
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "2) Bin with smallest contends at " << energies[minChannel]
    << ", backval=" << (backval/3.0) << ", lastbin=" << lastchannel
    << " x(lastbin)=" << energies[lastchannel] << "\n";
#endif
  
  
  if( direction < 0 && ((float(lastchannel)/float(nchannel)) < 0.04) )
  {
    size_t lower_channel = 0, upper_channel = 0;
    ExperimentalPeakSearch::find_spectroscopic_extent( dataH, lower_channel, upper_channel );
    if( static_cast<int>(lower_channel) >= lastchannel )
    {
      lastchannel = lower_channel ? static_cast<indexing_t>(lower_channel - 1) : 0;
      minChannel = std::max( minChannel, lastchannel );
    }
  }//if( direction < 0 && dataH->GetBinCenter(lastbin) < 100.0 )
  
  
  for( indexing_t channel = minChannel + direction*nSideChannel;
      channel != lastchannel && channel != meanchannel; channel += direction )
  {
    const float val = contents[channel];
    const float nextval = (channel>1 && (nchannel-channel)>0)  //probably is fine, but we'll check JIC
                          ? contents[channel+direction]
                          : contents[channel];
    
    float back = backval / nbackbin;
    float back_uncert = std::sqrt( backval ) / nbackbin;
    const float sigma = std::sqrt( back + back_uncert*back_uncert );
    //    float max_allowable = std::ceil(back + 2.575829*sigma) +  0.001;
    //    float min_allowable = std::floor(back - 2.575829*sigma) - 0.001;
    float max_allowable = std::ceil(back + 2.8f*sigma) +  0.001f;
    float min_allowable = std::floor(back - 2.8f*sigma) - 0.001f;
    
    if( back < 20.0f )
    {
      //Will use a reduced precision poisson_distribution to speed things up,
      //  since we really dont care past about 3 decimal places (but erring on
      //  the side of caution since I dont fully understand implications of
      //  reducing the accuracy).
      using boost::math::policies::digits10;
      using boost::math::poisson_distribution;
      typedef boost::math::policies::policy<digits10<6> > my_pol_6;
      
      const poisson_distribution<float, my_pol_6 > pois( back );
      max_allowable = boost::math::quantile( pois, 0.99f );
      min_allowable = boost::math::quantile( pois, 0.01f );
    }//if( val < 20.0 )
    
    //For high resolution spectra we'll require two bins to be outside of
    //  tollerance, since this will preserve the intent, but allow single bin
    //  spikes (which I swear are more common than poisson!) to not mess up the
    //  ROI.  It could probably be done for low resolution spectra, but I havent
    //  tested if (since single lowres peaks ROIs typically get calculated by
    //  find_roi_for_2nd_deriv_candidate(...) anyway
    if( (val>max_allowable && (!highres ||nextval>max_allowable))
        || (val<min_allowable && (!highres || nextval<min_allowable)) )
    {
#if( PRINT_ROI_DEBUG_INFO )
      if( debug_this_peak )
        DebugLog(cerr) << "Setting bin to " << (channel - direction) << " from expected "
        << lastchannel << " at energy=" << energies[channel-direction] << " kev"
        << ", this bring limit to " << (mean-energies[channel-direction])/sigma
        << " sigma from mean, val=" << val << ", min_allowable="
        << min_allowable << ", max_allowable=" << max_allowable << "\n";
#endif
      lastchannel = channel - (channel>0 ? direction : 0);
      break;
    }//if( val>max_allowable || val<min_allowable )
#if( PRINT_ROI_DEBUG_INFO )
    else
    {
      if( debug_this_peak )
        DebugLog(cerr) << "FindLimit: bin " << channel << ", x(bin)=" << energies[channel]
        << ", val=" << val
        << ", min_allowable=" << min_allowable
        << ", max_allowable=" << max_allowable
        << ", back=" << back << ", sigma=" << sigma << "\n";
    }
#endif
    
    nbackbin++;
    backval += val;
  }//for( int bin = minBin; bin > lastbin; --bin )
  
  
  //In principle, lastbin is the furthest from the mean we can end up, with
  //  the maximum being 11.75*sigma, or whereever a new feature was detected
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "minBin was " << minChannel << ", lastbin=" << lastchannel
    << " x(lastbin)=" << energies[lastchannel] << "\n";
#endif
  
  
  //Try to detect if there is a signficant skew on the peak by comparing
  //  4 to 7 sigma, to 7 to 11.75 sigma (or wherever is lastbin) to see if they
  //  are statistically compatible; if they are, just have ROI go to 7 sigma
  const int mean_channel = dataH->find_gamma_channel( mean );
  const int good_cont_channel = dataH->find_gamma_channel( mean + direction*7.05*sigma );
  if( (std::abs(int(lastchannel)-mean_channel) > std::abs(good_cont_channel-mean_channel)) )
  {
    const indexing_t nearest_channel = dataH->find_gamma_channel( mean + direction*3.5*sigma );
    const bool isNotDecreasing = isStatisticallyGreaterOrEqual( nearest_channel, good_cont_channel,
                                            good_cont_channel, lastchannel, dataH, 3.0 );

    
    
    if( high || isNotDecreasing )
    {
#if( PRINT_ROI_DEBUG_INFO )
      if( debug_this_peak )
        DebugLog(cerr) << "Setting lastchannel to " << lastchannel << "\n";
#endif
      lastchannel = good_cont_channel;
    }
    //    else
    //    {
    //      //now check from 11.75 to to 16.5 sigma, to see if we should include down
    //      //  to there.
    //      const int start = good_cont_bin;
    //      const size_t end_channel = dataH->find_gamma_channel( mean + direction*16.5*sigma );
    //      isNotDecreasing = isStatisticallyGreaterOrEqual( good_cont_bin, lastbin,
    //                                                      start, end_channel, dataH, 2.0 );
    //      if( !isNotDecreasing )
    //        lastbin = end;
    //    }
  }//if( abs(lastbin-mean_bin) > abs(good_cont_bin-mean_bin) )
  
  
  float val = dataH->gamma_channel_center(lastchannel);
  if( direction < 0 )
  {
    if( ((mean-val)/sigma) < 1.75 )
      lastchannel = dataH->find_gamma_channel( mean - 1.75*sigma );
  }else
  {
    if( ((val-mean)/sigma) < 1.75 )
      lastchannel = dataH->find_gamma_channel( mean + 1.75*sigma );
  }
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "Returning bin " << lastchannel
    << ", x=" << dataH->gamma_channel_center(lastchannel) << "\n";
#endif
  
  return lastchannel;
}//int findROILimit(...)




bool isStatisticallyGreaterOrEqual( const size_t start1, const size_t end1,
                                    const size_t start2, const size_t end2,
                                    const std::shared_ptr<const Measurement> &dataH,
                                    const double nsigma )
{
  size_t lowerbin = std::min( start1, end1 );
  size_t upperbin = std::max( start1, end1 );
  const double lower_area = dataH->gamma_channels_sum( lowerbin, upperbin );
  const int num_near_mean_bins = upperbin - lowerbin + 1;
  const double avrg_near_mean_area = lower_area / num_near_mean_bins;
  const double avrg_near_mean_uncert = sqrt(lower_area) / num_near_mean_bins;
  
  lowerbin = std::min( start2, end2 );
  upperbin = std::max( start2, end2 );
  const double upper_area = dataH->gamma_channels_sum( lowerbin, upperbin );
  const int num_tail_bins = upperbin - lowerbin + 1;
  const double tail_area = upper_area / num_tail_bins;
  const double avrg_tail_uncert = sqrt(upper_area) / num_tail_bins;
  
  const double uncert = sqrt(avrg_near_mean_uncert*avrg_near_mean_uncert
                             + avrg_tail_uncert*avrg_tail_uncert);
  
/*
#if( PRINT_ROI_DEBUG_INFO )
    DebugLog(cerr) << "Found avrg_near_mean_area=" << avrg_near_mean_area << "+-" << avrg_near_mean_uncert
    << ", tail_area=" << tail_area << "+-" << avrg_near_mean_uncert
    << ", uncert=" << uncert
    << " ---> (avrg_near_mean_area-2.0*uncert)=" << (avrg_near_mean_area-nsigma*uncert)
    << ", nsigma=" << ((avrg_near_mean_area-tail_area)/uncert)
    << "\n";
#endif
*/
  
  return (tail_area > (avrg_near_mean_area-nsigma*uncert));
}//isStatisticallyGreaterOrEqual( ... )


void estimatePeakFitRange( const PeakDef &peak, const std::shared_ptr<const Measurement> &dataH,
                           size_t &lower_channel, size_t &upper_channel )
{
  const size_t nchannel = dataH ? dataH->num_gamma_channels() : size_t(0);
  if( !nchannel )
    return;
  
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  double lowxrange  = continuum->lowerEnergy();
  double highxrange = continuum->upperEnergy();
  
  const bool definedRange = (lowxrange != highxrange);
  if( definedRange )
  {
    lower_channel  = max( dataH->find_gamma_channel(lowxrange), size_t(0) );
    upper_channel = min( dataH->find_gamma_channel(highxrange-0.00001), nchannel-1 );
    return;
  }//if( definedRange )
  
  const double mean = peak.mean();
  const double sigma = peak.gausPeak() ? peak.sigma() : 0.5*0.25*peak.roiWidth();
  
  
  if( continuum->type() == PeakContinuum::External )
  {
    lowxrange  = mean - 4.0*sigma;
    highxrange = mean + 4.0*sigma;
    
    lower_channel = dataH->find_gamma_channel(lowxrange);
    upper_channel = dataH->find_gamma_channel(highxrange);
    return;
  }//if( peak.m_offsetType == PeakDef::External )
  
  
  const bool polyContinuum = continuum->isPolynomial();
  if( polyContinuum )
  {
    lower_channel = findROILimit( peak, dataH, false );
    upper_channel = findROILimit( peak, dataH, true );
  }else
  {
    lower_channel = dataH->find_gamma_channel( mean - 4.0*sigma );
    upper_channel = dataH->find_gamma_channel( mean + 4.0*sigma );
  }
  
  if( lower_channel > upper_channel )
      std::swap( lower_channel, upper_channel );
      
  //Lets avoid some wierd going to too small of peak widths
  const size_t numfitbin = upper_channel - lower_channel;
  if( numfitbin <= 9 )  //9 chosen arbitrarily
  {
    lower_channel -= (10-numfitbin)/2;
    upper_channel += (10-numfitbin)/2;
  }//if( numfitbin <= 9 )
}//void setPeakXLimitsFromData( PeakDef &peak, const std::shared_ptr<const Measurement> &dataH )


ostream &operator<<( std::ostream &stream, const PeakContinuum &cont )
{
  switch( cont.type() )
  {
    case PeakContinuum::NoOffset:
      stream << "Underfined continuum";
    break;
      
    case PeakContinuum::External:
      stream << "Globally defined continuum";
    break;
      
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    {
      const char * const names[] = {"Flat", "Linear", "Bi-linear"};
      stream << names[cont.type() - PeakContinuum::FlatStep] << " step with coefficients {";
      for( size_t i = 0; i < cont.m_values.size(); ++i )
        stream << (i?", ":"") << cont.m_values[i];
      stream << "} relative to " << cont.m_referenceEnergy << " keV";
      break;
    }
      
    case PeakContinuum::Constant:   case PeakContinuum::Linear:
    case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
      stream << "Polynomial continuum with values {";
      for( size_t i = 0; i < cont.m_values.size(); ++i )
        stream << (i?", ":"") << cont.m_values[i];
      stream << "} relative to " << cont.m_referenceEnergy << " keV";
    break;
  }//switch( m_type )

  stream << ", valid from " << cont.lowerEnergy()
         << " to " << cont.upperEnergy() << " keV";
  
  return stream;
}//operator<<( std::ostream &stream, const PeakContinuum &cont )


std::ostream &operator<<( std::ostream &stream, const PeakDef &peak )
{
  stream << "mean=" << peak.m_coefficients[PeakDef::Mean];
  if( peak.m_uncertainties[PeakDef::Mean] > 0.0 )
    stream << "+-" << peak.m_uncertainties[PeakDef::Mean];
  stream << ", sigma=" << peak.m_coefficients[PeakDef::Sigma];
  if( peak.m_uncertainties[PeakDef::Sigma] > 0.0 )
    stream << "+-" << peak.m_uncertainties[PeakDef::Sigma];
  stream << ", amplitude=" << peak.m_coefficients[PeakDef::GaussAmplitude];
  if( peak.m_uncertainties[PeakDef::GaussAmplitude] > 0.0 )
    stream << "+-" << peak.m_uncertainties[PeakDef::GaussAmplitude];

  if( peak.m_transition )
  {
    const SandiaDecay::RadParticle *particle = NULL;
    const SandiaDecay::Transition *transition = peak.m_transition;
    const int index = peak.m_radparticleIndex;
    const int nproducts = static_cast<int>( transition->products.size() );
    if( (index>=0)  && (index<nproducts) )
      particle = &(transition->products[index]);
    stream << ", decay="
           << (transition->parent ? transition->parent->symbol : string("N/A") )
           << "->"
           << (transition->child ? transition->child->symbol : string("N/A") )
           << " " << (particle ? particle->energy : -1.0) << " keV";
  }else if( peak.m_sourceGammaType == PeakDef::AnnihilationGamma )
  {
    stream << ", Annihilation Gamma";
  }
  
  stream << ", " << *peak.m_continuum
         << ", chi2=" << peak.m_coefficients[PeakDef::Chi2DOF]
         << ", Skew0=" << peak.m_coefficients[PeakDef::SkewPar0]
         << ", Skew1=" << peak.m_coefficients[PeakDef::SkewPar1]
         << ", Skew2=" << peak.m_coefficients[PeakDef::SkewPar2]
         << ", Skew3=" << peak.m_coefficients[PeakDef::SkewPar3]
         << std::flush;
  return stream;
}//std::ostream &operator<<( std::ostream &stream, const PeakDef &peak )



#if( PERFORM_DEVELOPER_CHECKS )
void PeakDef::equalEnough( const PeakDef &lhs, const PeakDef &rhs )
{
  char buffer[512];
  
  if( lhs.m_userLabel != rhs.m_userLabel )
    throw runtime_error( "PeakDef user label for LHS ('"
                        + lhs.m_userLabel + "') doesnt match RHS ('"
                        + rhs.m_userLabel + "')" );
  
  if( lhs.m_type != rhs.m_type )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef peak type doesnt match, %i vs %i",
             int(lhs.m_type), int(rhs.m_type) );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_skewType != rhs.m_skewType )
  {
    snprintf(buffer, sizeof(buffer), "PeakDef skew type doesnt match, %i vs %i",
             int(lhs.m_skewType), int(rhs.m_skewType) );
    throw runtime_error( buffer );
  }
  
  
  for( CoefficientType t = CoefficientType(0); t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const double a = lhs.m_coefficients[t];
    const double b = rhs.m_coefficients[t];
    const double diff = fabs( a - b );
    
    if( diff > 1.0E-6*max(fabs(a),fabs(b)) )
    {
      snprintf(buffer, sizeof(buffer),
               "PeakDef coefficient %s of LHS (%1.8E) vs RHS (%1.8E) is out of tolerance.",
               to_string(t), lhs.m_coefficients[t], rhs.m_coefficients[t] );
      throw runtime_error( buffer );
    }
  }
  
  for( CoefficientType t = CoefficientType(0); t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const double a = lhs.m_uncertainties[t];
    const double b = rhs.m_uncertainties[t];
    const double diff = fabs( a - b );

    if( ((a > 0.0) || (b > 0.0)) && (diff > 1.0E-6*max(fabs(a),fabs(b))) )
    {
      snprintf(buffer, sizeof(buffer),
               "PeakDef uncertainty %s of LHS (%1.8E) vs RHS (%1.8E) is out of tolerance.",
               to_string(t), lhs.m_uncertainties[t], rhs.m_uncertainties[t] );
      throw runtime_error( buffer );
    }
  }
  
  for( CoefficientType t = CoefficientType(0); t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    // The "Fit For" values are not recorded for skew parameters that arent applicable, so we
    //  will get an erroneous fail, if we dont skip them, because we dont reset them when we change
    //  types
    bool skipSkew = false;
    switch( t )
    {
      case PeakDef::SkewPar0: case PeakDef::SkewPar1:
      case PeakDef::SkewPar2: case PeakDef::SkewPar3:
      {
        const int skew_par_num = static_cast<int>(t) - static_cast<int>(PeakDef::SkewPar0);
        const int nskew = static_cast<int>( PeakDef::num_skew_parameters( lhs.m_skewType ));
        skipSkew = (skew_par_num >= nskew);
        break;
      }
        
      default:
        break;
    }//switch( t )
    
    if( !skipSkew && (lhs.m_fitFor[t] != rhs.m_fitFor[t]) )
    {
      snprintf(buffer, sizeof(buffer),
               "PeakDef fit for %s of LHS (%i) vs RHS (%i) doesnt match.",
               to_string(t), int(lhs.m_fitFor[t]), int(rhs.m_fitFor[t]) );
      throw runtime_error( buffer );
    }
  }//for( check "Fit For" bools of all the parameters )
  
  
  if( !!lhs.m_continuum != !!rhs.m_continuum )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef availablity of continuum of LHS (%i) vs RHS (%i) continuums doesnt match.",
             int(!!lhs.m_continuum), int(!!rhs.m_continuum) );
    throw runtime_error( buffer );
  }

  if( !!lhs.m_continuum )
    PeakContinuum::equalEnough( *lhs.m_continuum, *rhs.m_continuum );
 
  if( lhs.m_parentNuclide != rhs.m_parentNuclide )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef parent nuclide of LHS (%s) vs RHS (%s) doesnt match.",
        (lhs.m_parentNuclide ? lhs.m_parentNuclide->symbol.c_str() : "none"),
        (rhs.m_parentNuclide ? rhs.m_parentNuclide->symbol.c_str() : "none") );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_parentNuclide && (lhs.m_transition != rhs.m_transition) )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef nuclide transition of LHS (%s -> %s) vs RHS (%s -> %s) doesnt match.",
             (lhs.m_transition->parent ? lhs.m_transition->parent->symbol.c_str() : "none"),
             (lhs.m_transition->child ? lhs.m_transition->child->symbol.c_str() : "none"),
             (rhs.m_transition->parent ? rhs.m_transition->parent->symbol.c_str() : "none"),
             (rhs.m_transition->child ? rhs.m_transition->child->symbol.c_str() : "none") );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_parentNuclide && (lhs.m_radparticleIndex != rhs.m_radparticleIndex) )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef particle index of LHS (%i) vs RHS (%i) doesnt match.",
             lhs.m_radparticleIndex, rhs.m_radparticleIndex );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_sourceGammaType != rhs.m_sourceGammaType )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef is annihilation of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_sourceGammaType), int(rhs.m_sourceGammaType) );
    throw runtime_error( buffer );
  }
  
//  std::vector< CandidateNuclide > m_candidateNuclides;
  
  if( lhs.m_xrayElement != rhs.m_xrayElement )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef xray of LHS (%s) vs RHS (%s) doesnt match.",
             (lhs.m_xrayElement ? lhs.m_xrayElement->symbol.c_str() : "none"),
             (rhs.m_xrayElement ? rhs.m_xrayElement->symbol.c_str() : "none") );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.m_xrayEnergy - rhs.m_xrayEnergy) > 0.001 )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef xray energy of LHS (%1.8E keV) vs RHS (%1.8E keV) doesnt match.",
             lhs.m_xrayEnergy, rhs.m_xrayEnergy );
    throw runtime_error( buffer );
  }

  if( lhs.m_reaction != rhs.m_reaction )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef reaction of LHS (%s) vs RHS (%s) doenst match.",
             (lhs.m_reaction ? lhs.m_reaction->name().c_str() : "none"),
             (rhs.m_reaction ? rhs.m_reaction->name().c_str() : "none") );
    throw runtime_error( buffer );
  }

  if( fabs(lhs.m_reactionEnergy - rhs.m_reactionEnergy) > 0.001 )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef reaction energy energy of LHS (%1.8E keV) vs RHS (%1.8E keV) doesnt match.",
             lhs.m_reactionEnergy, rhs.m_reactionEnergy );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_useForEnergyCal != rhs.m_useForEnergyCal )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef use for calibration of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_useForEnergyCal), int(rhs.m_useForEnergyCal) );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_useForShieldingSourceFit != rhs.m_useForShieldingSourceFit )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef use for shielding source fit of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_useForShieldingSourceFit), int(rhs.m_useForShieldingSourceFit) );
    throw runtime_error( buffer );
  }
  
  
  if( lhs.m_useForManualRelEff != rhs.m_useForManualRelEff )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef use for rel act from peaks of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_useForManualRelEff), int(rhs.m_useForManualRelEff) );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_useForDrfIntrinsicEffFit != rhs.m_useForDrfIntrinsicEffFit )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef use for DRF Abs Eff fit of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_useForDrfIntrinsicEffFit), int(rhs.m_useForDrfIntrinsicEffFit) );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_useForDrfFwhmFit != rhs.m_useForDrfFwhmFit )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef use for DRF FWHM fit of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_useForDrfFwhmFit), int(rhs.m_useForDrfFwhmFit) );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_useForDrfDepthOfInteractionFit != rhs.m_useForDrfDepthOfInteractionFit )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakDef use for DRF Depth of interaction fit of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_useForDrfDepthOfInteractionFit), int(rhs.m_useForDrfDepthOfInteractionFit) );
    throw runtime_error( buffer );
  }
}//void equalEnough( const PeakDef &lhs, const PeakDef &rhs )


void PeakContinuum::equalEnough( const PeakContinuum &lhs, const PeakContinuum &rhs )
{
  
  char buffer[512];
  
  if( lhs.m_type != rhs.m_type )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum type of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_type), int(rhs.m_type) );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.m_lowerEnergy - rhs.m_lowerEnergy) > 0.0001 )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum ROI lower energy of LHS (%1.8E keV) vs RHS (%1.8E keV) doesnt match.",
             lhs.m_lowerEnergy, rhs.m_lowerEnergy );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.m_upperEnergy - rhs.m_upperEnergy) > 0.0001 )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum ROI upper energy of LHS (%1.8E keV) vs RHS (%1.8E keV) doesnt match.",
             lhs.m_upperEnergy, rhs.m_upperEnergy );
    throw runtime_error( buffer );
  }

  if( fabs(lhs.m_referenceEnergy - rhs.m_referenceEnergy) > 0.0001 )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum reference energy of LHS (%1.8E keV) vs RHS (%1.8E keV) doesnt match.",
             lhs.m_referenceEnergy, rhs.m_referenceEnergy );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_values.size() != rhs.m_values.size() )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum number of coefficients of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_values.size()), int(rhs.m_values.size()) );
    throw runtime_error( buffer );
  }

  if( lhs.m_uncertainties.size() != rhs.m_uncertainties.size() )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum number of uncertainties of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_uncertainties.size()), int(rhs.m_uncertainties.size()) );
    throw runtime_error( buffer );
  }


  if( lhs.m_fitForValue.size() != rhs.m_fitForValue.size() )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum number of fit for variables of LHS (%i) vs RHS (%i) doesnt match.",
             int(lhs.m_fitForValue.size()), int(rhs.m_fitForValue.size()) );
    throw runtime_error( buffer );
  }
  
  if( lhs.m_values.size() != lhs.m_uncertainties.size()
     || lhs.m_values.size() != lhs.m_fitForValue.size() )
    throw runtime_error( "PeakContinuum something totally whack with number of coefficents somewhere!" );
  
  
  for( size_t i = 0; i < lhs.m_values.size(); ++i )
  {
    if( fabs(lhs.m_values[i]-rhs.m_values[i]) > (1.0E-5 * std::max(fabs(lhs.m_values[i]),fabs(rhs.m_values[i])))  )
    {
      snprintf(buffer, sizeof(buffer),
               "PeakContinuum value of %ith variables of LHS (%1.8E) vs RHS (%1.8E) doesnt match within tolerance.",
               int(i), lhs.m_values[i], rhs.m_values[i] );
      throw runtime_error( buffer );
    }

    if( fabs(lhs.m_uncertainties[i]-rhs.m_uncertainties[i]) > (1.0E-4 * std::max(fabs(lhs.m_uncertainties[i]),fabs(rhs.m_uncertainties[i])))  )
    {
      snprintf(buffer, sizeof(buffer),
               "PeakContinuum value of %ith uncertainty of LHS (%1.8E) vs RHS (%1.8E) doesnt match within tolerance.",
               int(i), lhs.m_uncertainties[i], rhs.m_uncertainties[i] );
      throw runtime_error( buffer );
    }

    if( lhs.m_fitForValue[i] != rhs.m_fitForValue[i] )
    {
      snprintf(buffer, sizeof(buffer),
               "PeakContinuum value of %ith fit for of LHS (%i) vs RHS (%i) doesnt match.",
               int(i), int(lhs.m_fitForValue[i]), int(rhs.m_fitForValue[i]) );
      throw runtime_error( buffer );
    }
  }
  
  if( !!lhs.m_externalContinuum != !!rhs.m_externalContinuum )
  {
    snprintf(buffer, sizeof(buffer),
             "PeakContinuum availablity of external continuum LHS (%i) vs RHS (%i) continuums doent match.",
             int(!!lhs.m_externalContinuum), int(!!rhs.m_externalContinuum) );
    throw runtime_error( buffer );
  }
  
  if( !!lhs.m_externalContinuum )
  {
    try
    {
      Measurement::equal_enough( *lhs.m_externalContinuum, *rhs.m_externalContinuum );
    }catch( std::exception &e )
    {
      snprintf( buffer, sizeof(buffer), "PeakContinuum caught testing external continuum: %s", e.what() );
      throw runtime_error( buffer );
    }//try / catch
  }
}//void equalEnough( const PeakContinuum &lhs, const PeakContinuum &rhs )
#endif //PERFORM_DEVELOPER_CHECKS


PeakDef::PeakDef()
{
  reset();
}


void PeakDef::reset()
{
  m_userLabel                 = "";
  m_type                      = GaussianDefined;
  m_skewType                  = PeakDef::NoSkew;

  m_parentNuclide             = NULL;
  m_transition                = NULL;
  m_radparticleIndex          = -1;
  m_sourceGammaType           = NormalGamma;
  m_useForEnergyCal           = true;
  m_useForShieldingSourceFit  = false;
  m_useForManualRelEff        = true;
  
  m_useForDrfIntrinsicEffFit       = PeakDef::sm_defaultUseForDrfIntrinsicEffFit;
  m_useForDrfFwhmFit               = PeakDef::sm_defaultUseForDrfFwhmFit;
  m_useForDrfDepthOfInteractionFit = PeakDef::sm_defaultUseForDrfDepthOfInteractionFit;
  
  m_xrayElement               = NULL;
  m_xrayEnergy                = 0.0;
  m_reaction                  = NULL;
  m_reactionEnergy            = 0.0;

  m_lineColor                 = Wt::WColor();
  
  std::shared_ptr<PeakContinuum> newcont = std::make_shared<PeakContinuum>();
  m_continuum = newcont;
  
  for( CoefficientType t = CoefficientType(0);
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    m_coefficients[t] = 0.0;
    m_uncertainties[t] = -1.0;
    
    switch( t )
    {
      case PeakDef::Mean:
      case PeakDef::Sigma:
      case PeakDef::GaussAmplitude:
      case PeakDef::SkewPar0:
      case PeakDef::SkewPar1:
      case PeakDef::SkewPar2:
      case PeakDef::SkewPar3:
        m_fitFor[t] = true;
      break;
        
      case PeakDef::Chi2DOF:
      case PeakDef::NumCoefficientTypes:
        m_fitFor[t] = false;
      break;
    }//switch( type )
  }//for( loop over coefficients )
}//void PeakDef::reset()


PeakDef::PeakDef( double m, double s, double a )
{
  reset();
  m_coefficients[PeakDef::Mean] = m;
  m_coefficients[PeakDef::Sigma] = s;
  m_coefficients[PeakDef::GaussAmplitude] = a;
}



PeakDef::PeakDef( double xlow, double xhigh, double mean,
                    std::shared_ptr<const Measurement> data, std::shared_ptr<const Measurement> background )
{
  reset();
  m_type = PeakDef::DataDefined;
  m_coefficients[PeakDef::Mean] = mean;
  m_continuum->setRange( xlow, xhigh );
  
  if( !data )
    return;

  m_continuum->setType( PeakContinuum::External );
  m_continuum->setExternalContinuum( background );
  m_continuum->setRange( xlow, xhigh );
  
  m_coefficients[PeakDef::GaussAmplitude] = gamma_integral( data, xlow, xhigh );
  if( background )
    m_coefficients[PeakDef::GaussAmplitude] -= gamma_integral( background, xlow, xhigh );
}//PeakDef( constructor )



const char *PeakDef::to_string( const CoefficientType type )
{
  switch( type )
  {
    case PeakDef::Mean:                return "Centroid";
    case PeakDef::Sigma:               return "Width";
    case PeakDef::GaussAmplitude:      return "Amplitude";
    case PeakDef::SkewPar0:            return "Skew0";
    case PeakDef::SkewPar1:            return "Skew1";
    case PeakDef::SkewPar2:            return "Skew2";
    case PeakDef::SkewPar3:            return "Skew3";
    case PeakDef::Chi2DOF:             return "Chi2";
    case PeakDef::NumCoefficientTypes: return "";
  }//switch( type )

  return "";
}//const char *PeakDef::to_string( const CoefficientType type )


const char *PeakDef::to_string( const SkewType type )
{
  switch( type )
  {
    case PeakDef::NoSkew:                 return "NoSkew";
    case PeakDef::Bortel:                 return "ExGauss"; //Could instead use "Bortel"
    case PeakDef::GaussExp:               return "GaussExp";
    case PeakDef::CrystalBall:            return "CrystalBall"; //Might change to "CB"
    case PeakDef::ExpGaussExp:            return "ExpGaussExp";
    case PeakDef::DoubleSidedCrystalBall: return "DoubleSidedCrystalBall"; //Might change to "DSCB"
  }//switch( skew_type )
  
  assert( 0 );
  throw runtime_error( "PeakDef::to_string(SkewType): invalid SkewType" );
  return "";
}//const char *to_string( const SkewType type )


const char *PeakDef::to_label( const SkewType type )
{
  switch( type )
  {
    case PeakDef::NoSkew:                 return "None";
    case PeakDef::Bortel:                 return "Exp*Gauss";
    case PeakDef::GaussExp:               return "GaussExp";
    case PeakDef::CrystalBall:            return "Crystal Ball";
    case PeakDef::ExpGaussExp:            return "ExpGaussExp";
    case PeakDef::DoubleSidedCrystalBall: return "Double Crystal Ball";
  }//switch( skew_type )
  
  assert( 0 );
  throw runtime_error( "PeakDef::to_string(SkewType): invalid SkewType" );
  return "";
}//const char *to_label( const SkewType type );


PeakDef::SkewType PeakDef::skew_from_string( const string &skew_type_str )
{
  if( SpecUtils::iequals_ascii(skew_type_str,"NoSkew") )
    return PeakDef::SkewType::NoSkew;
  
  if( SpecUtils::iequals_ascii(skew_type_str,"Bortel")
          || SpecUtils::iequals_ascii(skew_type_str,"ExGauss")  )
    return PeakDef::SkewType::Bortel;
  
  if( SpecUtils::iequals_ascii(skew_type_str,"GaussExp") )
    return PeakDef::SkewType::GaussExp;
  
  if( SpecUtils::iequals_ascii(skew_type_str,"CrystalBall")
          || SpecUtils::iequals_ascii(skew_type_str,"CB") )
    return PeakDef::SkewType::CrystalBall;
  
  if( SpecUtils::iequals_ascii(skew_type_str,"ExpGaussExp") )
    return PeakDef::SkewType::ExpGaussExp;
  
  if( SpecUtils::iequals_ascii(skew_type_str,"DoubleSidedCrystalBall")
          || SpecUtils::iequals_ascii(skew_type_str,"DSCB") )
    return PeakDef::SkewType::DoubleSidedCrystalBall;
  
  
  throw runtime_error( "Invalid peak skew type: " + string(skew_type_str) );
  
  return PeakDef::SkewType::NoSkew;
}//SkewType skew_from_string( const char *skew_str )


bool PeakDef::skew_parameter_range( const SkewType skew_type, const CoefficientType coef,
                                   double &lower_value, double &upper_value,
                                   double &starting_value, double &step_size )
{
  lower_value = upper_value = starting_value = step_size = 0.0;
  
  switch( skew_type )
  {
    case NoSkew:
      return false;
      
    case SkewType::Bortel:
    {
      if( coef != CoefficientType::SkewPar0 )
        return false;
      
      starting_value = 2.0;
      step_size = 1.0;
      lower_value = 0.0; //Below 0.005 would be numerically bad, but the Bortel function should protect against it.
      upper_value = 15;
      
      break;
    }//case SkewType::Bortel:
      
    case SkewType::CrystalBall:
    case SkewType::DoubleSidedCrystalBall:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar2: //alpha (right)
          if( skew_type == SkewType::CrystalBall )
            return false;
          // fall-though intentional
        case CoefficientType::SkewPar0: //alpha (left)
          starting_value = 2; // Saying skew becomes significant after 2 sigma, is maybe reasonable
          step_size = 0.5;
          lower_value = 0.5;  // You should at least be gaussian for half a sigma
          upper_value = 4.0;  // If you are gaussian all the way out to 4 sigma, you dont need skew
          break;
        
          
        case CoefficientType::SkewPar3: //n (right)
          if( skew_type == SkewType::CrystalBall )
            return false;
          // fall-though intentional
        case CoefficientType::SkewPar1: //n (left)
          // The valid values of `n` is probably a bit more complicated than just the range
          starting_value = 15;
          step_size = 2;
          lower_value = 1.05; //1.0 would be divide by zero
          upper_value = 100;  //much higher than this and we run into numerical issues
          break;
          
        default:
          return false;
      }//switch( coef )
      
      break;
    }//case SkewType::CrystalBall, SkewType::DoubleSidedCrystalBall:
      
    case SkewType::GaussExp:
    case SkewType::ExpGaussExp:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar1:
          if( skew_type == SkewType::GaussExp )
            return false;
          // fall-though intentional
        case CoefficientType::SkewPar0:
          starting_value = 1;  //A pretty good amount of skew
          step_size = 0.2;
          lower_value = 0.15;  //this is a huge amount of skew
          upper_value = 3.25;  //you really cant see any skew above ~2.75
          break;
          
        default:
          return false;
      }//switch( coef )
      
      break;
    }//case SkewType::ExpGaussExp: case SkewType::GaussExp:
  }//switch( skew_type )
  
  return true;
}//void skew_parameter_range(...)


size_t PeakDef::num_skew_parameters( const SkewType skew_type )
{
  switch ( skew_type )
  {
    case NoSkew:                 return 0;
    case Bortel:                 return 1;
    case CrystalBall:            return 2;
    case DoubleSidedCrystalBall: return 4;
    case GaussExp:               return 1;
    case ExpGaussExp:            return 2;
  }//switch ( skew_type )
  
  assert( 0 );
  throw std::logic_error( "Invalid peak skew" );
  return 0;
}//size_t num_skew_parameters( const SkewType skew_type );


bool PeakDef::is_energy_dependent( const SkewType skew_type, const CoefficientType coef )
{
  switch( skew_type )
  {
    case NoSkew:
      return false;
      
    case SkewType::Bortel:
      assert( coef == CoefficientType::SkewPar0 );
      return (coef == CoefficientType::SkewPar0);
      
    case SkewType::CrystalBall:
      assert( (coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1) );
      return (coef == CoefficientType::SkewPar0);
      
    case SkewType::DoubleSidedCrystalBall:
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar2));
      
    case SkewType::GaussExp:
      assert( coef == CoefficientType::SkewPar0 );
      return (coef == CoefficientType::SkewPar0);
      
    case SkewType::ExpGaussExp:
      assert( (coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1) );
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1));
  }//switch( skew_type )
  
  assert( 0 );
  return false;
}//bool is_energy_dependent( const SkewType skew_type, const CoefficientType coefficient )


double PeakDef::extract_energy_from_peak_source_string( std::string &str )
{  
  std::smatch energy_match;
  const std::regex energy_regexp("(?:^|\\s)((((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?)\\s*(kev|mev|ev|$))",
                                 regex::ECMAScript | regex::icase );
  
  if( !std::regex_search(str, energy_match, energy_regexp) )
    return -1.0;
  
  
  const string &val = energy_match[2];
  const string &units = energy_match[7];
  const string &total_match = energy_match[0];
  
  //cout << "Match for '"  << str << "': ";
  //for( const auto i : energy_match )
  //  cout << "'" << i << "', ";
  //cout << ", val='" << val << ", units='" << units << "'" << endl;
  
  double energy = -1.0;
  if( !(stringstream(val) >> energy) )
  {
    assert( 0 ); //should ever get here if regex is well formed
    return -1.0;
  }
  
  if( SpecUtils::iequals_ascii(units, "kev") )
  {
    // Nothing to do here
  }else if( SpecUtils::iequals_ascii(units, "mev") )
  {
    energy *= 1000.0;
  }else if( SpecUtils::iequals_ascii(units, "ev") )
  {
    energy /= 1000.0;
  }else if( units.empty() )
  {
    if( (energy < 295.0) && (str.find('.') == string::npos) )
    {
      // If value is less than 295, and there is no decimal, then its possible we've picked up on
      //  isotope number (e.g., the 235 from U235), so for the moment, we'll reject this match,
      //  unless we unambiguously know there were numbers or reaction before what we think is the
      //  energy
      const auto matchpos = str.find(total_match);
      auto first_num_pos = str.find_first_of( "0123456789)" );
      //if( first_num_pos > matchpos )
      //  first_num_pos = str.find_first_of( "x-ray" );
      //if( first_num_pos > matchpos )
      //  first_num_pos = str.find_first_of( "xray" );
      //if( first_num_pos > matchpos )
      //  first_num_pos = str.find_first_of( "x ray" );
      
      if( first_num_pos < matchpos )
      {
        // We have an x-ray, or there were numbers, or closing parenthesis before our match
      }else
      {
        return -1.0;
      }
    }//if( (energy < 295.0) && (str.find('.') == string::npos) )
  }
  
  SpecUtils::ireplace_all( str, total_match.c_str(), " " );
  SpecUtils::trim( str );
  
  return energy;
};//extract_energy_from_peak_source_string


void PeakDef::gammaTypeFromUserInput( std::string &txt,
                                      PeakDef::SourceGammaType &type )
{
  
  type = PeakDef::NormalGamma;
  
  if( SpecUtils::icontains( txt, "s.e." ) )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, "s.e.", "" );
  }else if( SpecUtils::icontains( txt, "single escape" ) )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, "single escape", "" );
  }else if( SpecUtils::iequals_ascii( txt, "s.e." )
     || SpecUtils::iequals_ascii( txt, "se" )
     || SpecUtils::iequals_ascii( txt, "escape" ) )
  {
    type = PeakDef::SingleEscapeGamma;
    txt = "";
  }else if( SpecUtils::icontains( txt, "se " ) && txt.size() > 5 )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, "se ", "" );
  }else if( SpecUtils::icontains( txt, " se" ) && txt.size() > 5 )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, " se", "" );
  }else if( SpecUtils::icontains( txt, "d.e." ) )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, "d.e.", "" );
  }else if( SpecUtils::icontains( txt, "double escape" ) )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, "double escape", "" );
  }else if( SpecUtils::icontains( txt, "de " ) && txt.size() > 5 )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, "de ", "" );
  }else if( SpecUtils::icontains( txt, " de" ) && txt.size() > 5 )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, " de", "" );
  }else if( SpecUtils::iequals_ascii( txt, "d.e." )
           || SpecUtils::iequals_ascii( txt, "de" )  )
  {
    type = PeakDef::DoubleEscapeGamma;
    txt = "";
  }else if( SpecUtils::icontains( txt, "x-ray" )
      || SpecUtils::icontains( txt, "xray" )
     || SpecUtils::icontains( txt, "x ray" ) )
  {
    type = PeakDef::XrayGamma;
    SpecUtils::ireplace_all( txt, "xray", "" );
    SpecUtils::ireplace_all( txt, "x-ray", "" );
    SpecUtils::ireplace_all( txt, "x ray", "" );
  }
}//PeakDef::SourceGammaType gammaType( std::string txt )


const Wt::WColor &PeakDef::lineColor() const
{
  return m_lineColor;
}

void PeakDef::setLineColor( const Wt::WColor &color )
{
  m_lineColor = color;
}


void PeakContinuum::toXml( rapidxml::xml_node<char> *parent, const int contId ) const
{
  using namespace rapidxml;
  
  xml_document<char> *doc = parent ? parent->document() : (xml_document<char> *)0;
  
  if( !doc )
    throw runtime_error( "PeakContinuum::toXml(...): invalid input" );
  
  char buffer[128];
  xml_node<char> *node = 0;
  xml_node<char> *cont_node = doc->allocate_node( node_element, "PeakContinuum" );
  
  // A reminder double check these logics when changing PeakContinuum::sm_xmlSerializationVersion
  static_assert( PeakContinuum::sm_xmlSerializationVersion == 1,
                "PeakContinuum::toXml needs to be updated for new serialization version." );
  
  // For version 1.0.8 and newer InterSpec, we will attempt to let InterSpec v1.0.7 and older be
  //  able to read the peaks in N42 files, as long as no stepped-continuums are used.
  int version = PeakContinuum::sm_xmlSerializationVersion;
  switch( m_type )
  {
    case NoOffset: case External: case Constant: case Linear: case Quadratic: case Cubic:
      // Nothing changed for these continuum types between version 0 and version 1.
      version = 0;
      break;
      
    case FlatStep: case LinearStep: case BiLinearStep:
      // These continuums were added for serialization version 1, starting with InterSpec v1.0.8.
      version = 1;
      break;
  }//switch( m_type )
  
  snprintf( buffer, sizeof(buffer), "%i", version );
  const char *val = doc->allocate_string( buffer );
  xml_attribute<char> *att = doc->allocate_attribute( "version", val );
  cont_node->append_attribute( att );

  snprintf( buffer, sizeof(buffer), "%i", contId );
  val = doc->allocate_string( buffer );
  att = doc->allocate_attribute( "id", val );
  cont_node->append_attribute( att );
    
  parent->append_node( cont_node );
  
  const char *type = offset_type_str(m_type);
  node = doc->allocate_node( node_element, "Type", type );
  cont_node->append_node( node );
  
  // Note: it could be energy range is not defined, in which case we could just
  //       not write "LowerEnergy" and "UpperEnergy"
  snprintf( buffer, sizeof(buffer), "%1.8e", m_lowerEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "LowerEnergy", val );
  cont_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%1.8e", m_upperEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "UpperEnergy", val );
  cont_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%1.8e", m_referenceEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "ReferenceEnergy", val );
  cont_node->append_node( node );
  
  if( m_type != NoOffset && m_type != External )
  {
    stringstream valsstrm, uncertstrm, fitstrm;
    for( size_t i = 0; i < m_values.size(); ++i )
    {
      const char *spacer = (i ? " " : "");
      snprintf( buffer, sizeof(buffer), "%1.8e", m_values[i] );  
      valsstrm << spacer << buffer;
    
      snprintf( buffer, sizeof(buffer), "%1.8e", m_uncertainties[i] );  
      uncertstrm << spacer << buffer;
    
      fitstrm << spacer << (m_fitForValue[i] ? '1': '0');
    }//for( size_t i = 0; i < m_values.size(); ++i )
    
    xml_node<char> *coeffs_node = doc->allocate_node( node_element, "Coefficients" );
    cont_node->append_node( coeffs_node );
        
    val = doc->allocate_string( valsstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Values", val );
    coeffs_node->append_node( node );
    
    val = doc->allocate_string( uncertstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Uncertainties", val );
    coeffs_node->append_node( node );
    
    val = doc->allocate_string( fitstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Fittable", val );
    coeffs_node->append_node( node );
  }//if( m_type != NoOffset && m_type != External )
  
  if( !!m_externalContinuum )
  {
    stringstream contXml;
    m_externalContinuum->write_2006_N42_xml( contXml );
    //We actually need to parse the XML here, and then insert it into the hierarchy
    
    const string datastr = contXml.str();
    std::unique_ptr<char []> data( new char [datastr.size()+1] );
    memcpy( data.get(), datastr.c_str(), datastr.size()+1 );
    
    xml_document<char> contdoc;
    const int flags = rapidxml::parse_normalize_whitespace
                     | rapidxml::parse_trim_whitespace;
    contdoc.parse<flags>( data.get() );
    
    node = doc->allocate_node( node_element, "ExternalContinuum", val );
    cont_node->append_node( node );
    
    xml_node<char> *spec_node = contdoc.first_node( "Measurement", 11 );
    if( !spec_node )
      throw runtime_error( "Didnt get expected Measurement node" );
    spec_node = spec_node->first_node( "Spectrum", 8 );
    if( !spec_node )
      throw runtime_error( "Didnt get expected Spectrum node" );
    
    xml_node<char> *new_spec_node = doc->allocate_node( node_element );
    node->append_node( new_spec_node );
    
    clone_node_deep( spec_node, new_spec_node );
  }//if( !!m_externalContinuum )
}//void PeakContinuum::toXml(...)





void PeakContinuum::fromXml( const rapidxml::xml_node<char> *cont_node, int &contId )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !cont_node )
    throw runtime_error( "PeakContinuum::fromXml(...): invalid input" );
  
  if( !compare( cont_node->name(), cont_node->name_size(), "PeakContinuum", 13, false ) )
    throw std::logic_error( "PeakContinuum::fromXml(...): invalid input node name" );
  
  xml_attribute<char> *att = cont_node->first_attribute( "version", 7 );
  
  int version;
  if( !att || !att->value() || (sscanf(att->value(), "%i", &version)!=1) )
    throw runtime_error( "PeakContinuum invalid version" );
  
  // A reminder double check these logics when changing PeakContinuum::sm_xmlSerializationVersion
  static_assert( PeakContinuum::sm_xmlSerializationVersion == 1,
                "PeakContinuum::toXml needs to be updated for new serialization version." );
  
  // Serialization version 1 is backwards compatible with version 0 for de-serialization, so no
  //  changes to this code is needed.
  if( (version < 0) || (version > PeakContinuum::sm_xmlSerializationVersion) )
    throw runtime_error( "Invalid PeakContinuum version: " + std::to_string(version) + ".  "
                    + "Only up to version " + to_string(PeakContinuum::sm_xmlSerializationVersion)
                    + " supported." );
  
  att = cont_node->first_attribute( "id", 2 );
  if( !att || !att->value() || (sscanf(att->value(), "%i", &contId)!=1) )
    throw runtime_error( "PeakContinuum invalid ID" );

  xml_node<char> *node = cont_node->first_node( "Type", 4 );

  if( !node || !node->value() )
    throw runtime_error( "PeakContinuum not Type node" );
  
  m_type = str_to_offset_type_str( node->value(), node->value_size() );
  
  float dummyval;
  node = cont_node->first_node( "LowerEnergy", 11 );
  if( node )
  {
    if( !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
      throw runtime_error( "Continuum didnt have valid LowerEnergy" );
    
    m_lowerEnergy = dummyval;
    
    node = cont_node->first_node( "UpperEnergy", 11 );
    if( !node || !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
      throw runtime_error( "Continuum didnt have UpperEnergy" );
    m_upperEnergy = dummyval;
  }else
  {
    m_lowerEnergy = m_upperEnergy = 0.0;
    node = cont_node->first_node( "UpperEnergy", 11 );
    if( node )
      throw runtime_error( "Continuum didnt have LowerEnergy, but did have UpperEnergy" );
  }//if( have <LowerEnergy> ) / else
  
  
  node = cont_node->first_node( "ReferenceEnergy", 15 );
  if( node )
  {
    if( !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
      throw runtime_error( "Continuum didnt have ReferenceEnergy" );
    m_referenceEnergy = dummyval;
  }else
  {
    if( m_lowerEnergy != m_upperEnergy )
      throw runtime_error( "Continuum didnt have ReferenceEnergy, but did have energy range defined" );
    m_referenceEnergy = 0.0;
  }//if( have <ReferenceEnergy> ) / else
  
  if( m_type != NoOffset && m_type != External )
  {
    xml_node<char> *coeffs_node = cont_node->first_node("Coefficients",12);
    if( !coeffs_node )
      throw runtime_error( "Continuum didnt have Coefficients node" );
    
    std::vector<float> contents;
    node = coeffs_node->first_node( "Values", 6 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Values" );
    
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_values.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_values[i] = contents[i]; 
    
    
    node = coeffs_node->first_node( "Uncertainties", 13 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Uncertainties" );  
    
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_uncertainties.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_uncertainties[i] = contents[i]; 
    
    
    node = coeffs_node->first_node( "Fittable", 8 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Fittable" );  
    
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_fitForValue.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_fitForValue[i] = (contents[i] > 0.5f); 
    
    if( m_values.size() != m_uncertainties.size() 
        || m_fitForValue.size() != m_values.size() )
      throw runtime_error( "Continuum coefficients not consistent" );
  }else
  {
    m_values.clear();
    m_uncertainties.clear();
    m_fitForValue.clear();
  }//if( m_type != NoOffset && m_type != External ) / else
  
  
  node = cont_node->first_node( "ExternalContinuum", 17 );
  if( node )
  {
    node = node->first_node( "Spectrum", 8 );
    if( !node )
      throw runtime_error( "Spectrum node expected under ExternalContinuum" );
    std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
    m_externalContinuum = meas;
      
    meas->set_info_from_2006_N42_spectrum_node( node );
  }//if( node )
}//void PeakContinuum::fromXml(...)



rapidxml::xml_node<char> *PeakDef::toXml( rapidxml::xml_node<char> *parent,
                     rapidxml::xml_node<char> *continuum_parent,
           std::map<std::shared_ptr<PeakContinuum>,int> &continuums ) const
{
  using namespace rapidxml;
  
  xml_document<char> *doc = parent ? parent->document() : (xml_document<char> *)0;
  
  if( !doc )
    throw runtime_error( "PeakDef::toXml(...): invalid input" );
  
  if( !m_continuum )
    throw logic_error( "PeakDef::toXml(...): continuum should be valid" );
  
  if( !continuums.count(m_continuum) )
  {
    const int index = static_cast<int>( continuums.size() + 1 );
    m_continuum->toXml( continuum_parent, index );
    continuums[m_continuum] = index;
  }//if( !continuums.count(m_continuum) )
  
  char buffer[128];
  const int contID = continuums[m_continuum];
  
  xml_node<char> *node = 0;
  xml_node<char> *peak_node = doc->allocate_node( node_element, "Peak" );
  
  snprintf( buffer, sizeof(buffer), "%i.%i",
            PeakDef::sm_xmlSerializationMajorVersion,
            PeakDef::sm_xmlSerializationMinorVersion );
  const char *val = doc->allocate_string( buffer );
  xml_attribute<char> *att = doc->allocate_attribute( "version", val );
  peak_node->append_attribute( att );
  
  snprintf( buffer, sizeof(buffer), "%i", contID );
  val = doc->allocate_string( buffer );
  att = doc->allocate_attribute( "continuumID", val );
  peak_node->append_attribute( att );
  
  parent->append_node( peak_node );
  
  if( m_userLabel.size() )
  {
    val = doc->allocate_string( m_userLabel.c_str() );
    node = doc->allocate_node( node_element, "UserLabel", val );
    peak_node->append_node( node );
  }//if( m_userLabel.size() )
  
  
  switch( m_type )
  {
    case GaussianDefined: val = "GaussianDefined"; break;
    case DataDefined:     val = "DataDefined";     break;
  }//switch( m_type )
  
  node = doc->allocate_node( node_element, "Type", val );
  peak_node->append_node( node );
  
  
  switch( m_skewType )
  {
    case NoSkew:                 val = "NoSkew";                 break;
    case Bortel:                 val = "ExGauss";                break;
    //case Doniach:                val = "Doniach";                break;
    case CrystalBall:            val = "CrystalBall";            break;
    case DoubleSidedCrystalBall: val = "DoubleSidedCrystalBall"; break;
    case GaussExp:               val = "GaussExp";               break;
    case ExpGaussExp:            val = "ExpGaussExp";            break;
  }//switch( m_skewType )
  
  node = doc->allocate_node( node_element, "Skew", val );
  peak_node->append_node( node );
  
  
  for( CoefficientType t = CoefficientType(0); 
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const char *label = to_string( t );
    
    // In going from PeakDef XML version 0.1 to 0.2, we need to have "LandauAmplitude",
    //  "LandauMode", and "LandauSigma" elements written, even if there is no skew, since
    //  InterSpec v1.0.11 and before always expect these elements.
    //  If we have a non-Landau skew, then older versions of InterSpec will fail anyway
    const char *extra_label = nullptr;
    switch( t )
    {
      case Mean: case Sigma: case GaussAmplitude:
      case Chi2DOF: case NumCoefficientTypes:
        break;
        
      case SkewPar0:
        if( (m_skewType == NoSkew) )
          extra_label = "LandauAmplitude";
        if( m_skewType == NoSkew )
          label = nullptr;
        break;
        
      case SkewPar1:
        if( (m_skewType == NoSkew) )
          extra_label = "LandauMode";
        if( (m_skewType != CrystalBall)
           && (m_skewType != DoubleSidedCrystalBall)
           && (m_skewType != ExpGaussExp) )
        {
          label = nullptr;
        }
        break;
        
      case SkewPar2:
        if( (m_skewType == NoSkew) )
          extra_label = "LandauSigma";
        if( m_skewType != DoubleSidedCrystalBall )
          label = nullptr;
        break;
        
      case SkewPar3:
        if( m_skewType != DoubleSidedCrystalBall )
          label = nullptr;
        break;
    }//switch( t )
    
    auto add_node = [this, &buffer, doc, peak_node, t]( const char *label, const bool is_deprecated ){
      double value = m_coefficients[t];
      double uncert = m_uncertainties[t];
      
      snprintf( buffer, sizeof(buffer), "%1.8e %1.8e", value, uncert );
      const char *val = doc->allocate_string( buffer );
      xml_node<char> *node = doc->allocate_node( node_element, label, val );
      xml_attribute<char> *att = doc->allocate_attribute( "fit", (m_fitFor[t] ? "true" : "false") );
      node->append_attribute( att );
      
      if( is_deprecated )
      {
        att = doc->allocate_attribute( "remark", "Deprecated element" );
        node->append_attribute( att );
      }
      
      peak_node->append_node( node );
    };//add_node
    
    if( label )
      add_node( label, false );
    if( extra_label )
      add_node( extra_label, true );
  }//for(...)
  
  /// TODO: Need to deprecate 'forCalibration' in favor of 'useForEnergyCalibration' the next
  ///       increment of PeakDef::sm_xmlSerializationMajorVersion
  att = doc->allocate_attribute( "forCalibration", (m_useForEnergyCal ? "true" : "false") );
  peak_node->append_attribute( att );
  
  att = doc->allocate_attribute( "useForEnergyCalibration", (m_useForEnergyCal ? "true" : "false") );
  peak_node->append_attribute( att );
  
  att = doc->allocate_attribute( "source", (m_useForShieldingSourceFit ? "true" : "false") );
  peak_node->append_attribute( att );
  
  att = doc->allocate_attribute( "useForManualRelEff", (m_useForManualRelEff ? "true" : "false") );
  peak_node->append_attribute( att );
  
  // Dont bother writing useForDrfIntrinsicEffFit, useForDrfFwhmFit, useForDrfDepthOfInteractionFit,
  //  unless their values have been set to true (when de-serializing them we will set to false
  //  if the attributes arent found)
  if( m_useForDrfIntrinsicEffFit != PeakDef::sm_defaultUseForDrfIntrinsicEffFit )
  {
    att = doc->allocate_attribute( "useForDrfIntrinsicEffFit", (m_useForDrfIntrinsicEffFit ? "true" : "false") );
    peak_node->append_attribute( att );
  }
  
  if( m_useForDrfFwhmFit != PeakDef::sm_defaultUseForDrfFwhmFit )
  {
    att = doc->allocate_attribute( "useForDrfFwhmFit", (m_useForDrfFwhmFit ? "true" : "false") );
    peak_node->append_attribute( att );
  }
  
  if( m_useForDrfDepthOfInteractionFit != PeakDef::sm_defaultUseForDrfDepthOfInteractionFit )
  {
    att = doc->allocate_attribute( "useForDrfDepthOfInteractionFit", (m_useForDrfDepthOfInteractionFit ? "true" : "false") );
    peak_node->append_attribute( att );
  }
  
  if( !m_lineColor.isDefault() )
  {
    //Added 20181027 without incrementing XML version since we're making it optional
    val = doc->allocate_string( m_lineColor.cssText(false).c_str() );  //Note: not including alpha because of Wt bug
    node = doc->allocate_node( node_element, "LineColor", val );
    peak_node->append_node( node );
  }//
  
  
  const char *gammaTypeVal = 0;
  switch( m_sourceGammaType )
  {
    case PeakDef::NormalGamma:       gammaTypeVal = "NormalGamma";       break;
    case PeakDef::AnnihilationGamma: gammaTypeVal = "AnnihilationGamma"; break;
    case PeakDef::SingleEscapeGamma: gammaTypeVal = "SingleEscapeGamma"; break;
    case PeakDef::DoubleEscapeGamma: gammaTypeVal = "DoubleEscapeGamma"; break;
    case PeakDef::XrayGamma:         gammaTypeVal = "XrayGamma";         break;
  }//switch( m_sourceGammaType )

  
  if( m_parentNuclide )
  {    
    xml_node<char> *nuc_node = doc->allocate_node( node_element, "Nuclide" );
    peak_node->append_node( nuc_node );
    
    val = doc->allocate_string( m_parentNuclide->symbol.c_str() );
    node = doc->allocate_node( node_element, "Name", val );
    nuc_node->append_node( node );
    
    if( m_transition )
    {
      string transistion_parent, decay_child;

      const SandiaDecay::Nuclide *trans_parent = m_transition->parent;
      transistion_parent = trans_parent->symbol;
      if( m_transition->child )
        decay_child = m_transition->child->symbol;
      const double energy = m_transition->products[m_radparticleIndex].energy;
      
      val = doc->allocate_string( transistion_parent.c_str() );
      node = doc->allocate_node( node_element, "DecayParent", val );
      nuc_node->append_node( node );
      
      val = doc->allocate_string( decay_child.c_str() );
      node = doc->allocate_node( node_element, "DecayChild", val );
      nuc_node->append_node( node );
      
      snprintf( buffer, sizeof(buffer), "%1.8e", energy );
      val = doc->allocate_string( buffer );
      node = doc->allocate_node( node_element, "DecayGammaEnergy", val );
      nuc_node->append_node( node );
    }//if( m_transition )
    
    node = doc->allocate_node( node_element, "DecayGammaType", gammaTypeVal );
    nuc_node->append_node( node );
  }//if( m_parentNuclide )
  
  if( m_xrayElement )
  {
    xml_node<char> *xray_node = doc->allocate_node( node_element, "XRay" );
    peak_node->append_node( xray_node );
    
    val = doc->allocate_string( m_xrayElement->symbol.c_str() );
    node = doc->allocate_node( node_element, "Element", val );
    xray_node->append_node( node );
    
    snprintf( buffer, sizeof(buffer), "%1.8e", m_xrayEnergy );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "Energy", val );
    xray_node->append_node( node );
  }//if( m_xrayElement )
  
  if( m_reaction )
  {
    xml_node<char> *rctn_node = doc->allocate_node( node_element, "Reaction" );
    peak_node->append_node( rctn_node );
    
    val = doc->allocate_string( m_reaction->name().c_str() );
    node = doc->allocate_node( node_element, "Name", val );
    rctn_node->append_node( node );
    
    snprintf( buffer, sizeof(buffer), "%1.8e", m_reactionEnergy );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "Energy", val );
    rctn_node->append_node( node );
    
    node = doc->allocate_node( node_element, "Type", gammaTypeVal );
    rctn_node->append_node( node );
  }//if( m_reaction )
  
  return peak_node;
}//rapidxml::xml_node<char> *toXml(...)



void PeakDef::fromXml( const rapidxml::xml_node<char> *peak_node,
             const std::map<int,std::shared_ptr<PeakContinuum> > &continuums )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !peak_node )
    throw logic_error( "PeakDef::fromXml(...): invalid input node" );
  
  if( !compare( peak_node->name(), peak_node->name_size(), "Peak", 4, false ) )
    throw std::logic_error( "PeakDef::fromXml(...): invalid input node name" );
  
  reset();
  
  int contID;
  xml_attribute<char> *att = peak_node->first_attribute( "continuumID", 11 );
  if( !att )
    throw runtime_error( "No continuum ID" );
  if( sscanf( att->value(), "%i", &contID ) != 1 )
    throw runtime_error( "Non integer continuum ID" );
  
  std::map<int,std::shared_ptr<PeakContinuum> >::const_iterator contpos;
  contpos = continuums.find( contID );
  if( contpos == continuums.end() )
    throw runtime_error( "Couldnt find valid continuum for peak" );
  
  m_continuum = contpos->second;
  
  /// TODO: Need to deprecate 'forCalibration' in favor of 'useForEnergyCalibration' next
  ///       PeakDef::sm_xmlSerializationMajorVersion version increment.
  att = peak_node->first_attribute( "forCalibration", 14 );
  //if( !att )
  //  att = peak_node->first_attribute( "useForEnergyCalibration", 23 );
  
  if( !att )
    throw runtime_error( "missing forCalibration attribute" );
  
  m_useForEnergyCal = compare(att->value(),att->value_size(),"true",4,false);
  if( !m_useForEnergyCal && !compare(att->value(),att->value_size(),"false",5,false) )
    throw runtime_error( "invalid forCalibration value" );
  
  att = peak_node->first_attribute( "source", 6 );
  if( !att )
    throw runtime_error( "missing source attribute" );
  m_useForShieldingSourceFit = compare(att->value(),att->value_size(),"true",4,false);
  if( !m_useForShieldingSourceFit && !compare(att->value(),att->value_size(),"false",5,false) )
    throw runtime_error( "invalid source value" );
  
  // useForManualRelEff was added June 2022, so for backward compatibility, dont require it
  att = peak_node->first_attribute( "useForManualRelEff", 18 );
  if( att )
  {
    m_useForManualRelEff = compare(att->value(),att->value_size(),"true",4,false);
    if( !m_useForManualRelEff && !compare(att->value(),att->value_size(),"false",5,false) )
      throw runtime_error( "invalid useForManualRelEff value" );
  }else
  {
    m_useForManualRelEff = true;
  }
  
  m_useForDrfIntrinsicEffFit = PeakDef::sm_defaultUseForDrfIntrinsicEffFit;
  att = peak_node->first_attribute( "useForDrfIntrinsicEffFit", 24 );
  if( !att )
    att = peak_node->first_attribute( "useForDrfFit", 12 );
  if( att )
    m_useForDrfIntrinsicEffFit = compare(att->value(),att->value_size(),"true",4,false);
  
  m_useForDrfFwhmFit = PeakDef::sm_defaultUseForDrfFwhmFit;
  att = peak_node->first_attribute( "useForDrfFwhmFit", 16 );
  if( att )
    m_useForDrfFwhmFit = compare(att->value(),att->value_size(),"true",4,false);
  
  m_useForDrfDepthOfInteractionFit = PeakDef::sm_defaultUseForDrfDepthOfInteractionFit;
  att = peak_node->first_attribute( "useForDrfDepthOfInteractionFit", 30 );
  if( att )
    m_useForDrfDepthOfInteractionFit = compare(att->value(),att->value_size(),"true",4,false);
  
  att = peak_node->first_attribute( "version", 7 );
  if( !att || !att->value_size() )
    throw runtime_error( "missing version attribute" );
  
  int majorVersion = 0, minorVersion = 0;
  const char * const version_begin = att->value();
  const char * const version_end = version_begin + att->value_size();
  
  if( sscanf( version_begin, "%i", &majorVersion ) != 1 )
    throw runtime_error( "Non integer version number" );
  
  if( majorVersion != PeakDef::sm_xmlSerializationMajorVersion )
    throw runtime_error( "Invalid peak version" );
  
  const char *version_period = std::find( version_begin, version_end, '.' );
  if( (version_period != version_end) && ((version_period+1) != version_end) )
  {
    if( sscanf( (version_period + 1), "%i", &minorVersion ) != 1 )
      throw runtime_error( "Non integer minor version number" );
  }
  
  // In the future we could use majorVersion and minorVersion to adjust parsing behaviour
  
  xml_node<char> *node = peak_node->first_node("UserLabel",9);
  if( node && node->value() )
    m_userLabel = node->value();
  
  node = peak_node->first_node("Type",4);
  if( !node || !node->value() )
    throw runtime_error( "No peak type" );
  
  if( compare(node->value(),node->value_size(),"GaussianDefined",15,false) )
    m_type = GaussianDefined;
  else if( compare(node->value(),node->value_size(),"DataDefined",11,false) )
    m_type = DataDefined;
  else
    throw runtime_error( "Invalid peak type" );
  
  node = peak_node->first_node("Skew",4);
  //if( !node || !node->value() )
  //  throw runtime_error( "No peak skew type" );
  
  if( !node || !node->value_size() || compare(node->value(),node->value_size(),"NoSkew",6,false) )
    m_skewType = NoSkew;
  else if( compare(node->value(),node->value_size(),"LandauSkew",10,false) )
    m_skewType = NoSkew;  //Pre 20231101, LandauSkew was defined, but I dont think ever use, so just throw it awa
  else if( compare(node->value(),node->value_size(),"Bortel",6,false)
          || compare(node->value(),node->value_size(),"ExGauss",7,false) )
    m_skewType = Bortel;
  //else if( compare(node->value(),node->value_size(),"Doniach",7,false) )
  //  m_skewType = Doniach;
  else if( compare(node->value(),node->value_size(),"CrystalBall",11,false)
          || compare(node->value(),node->value_size(),"CB",2,false) )
    m_skewType = CrystalBall;
  else if( compare(node->value(),node->value_size(),"DoubleSidedCrystalBall",22,false)
          || compare(node->value(),node->value_size(),"DSCB",4,false))
    m_skewType = DoubleSidedCrystalBall;
  else if( compare(node->value(),node->value_size(),"GaussExp",8,false) )
    m_skewType = GaussExp;
  else if( compare(node->value(),node->value_size(),"ExpGaussExp",11,false) )
    m_skewType = ExpGaussExp;
  else
    throw runtime_error( "Invalid peak skew type" );
  
  for( CoefficientType t = CoefficientType(0); 
      t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    // We dont expect all the skew coefficients
    bool want_par = true;
    switch( t )
    {
      case Mean: case Sigma: case GaussAmplitude:
      case Chi2DOF: case NumCoefficientTypes:
        want_par = true;
        break;
        
      case SkewPar0:
        want_par = (m_skewType != NoSkew);
        break;
        
      case SkewPar1:
        want_par = ((m_skewType == CrystalBall) || (m_skewType == DoubleSidedCrystalBall));
        break;
        
      case SkewPar2:
      case SkewPar3:
        want_par = (m_skewType == DoubleSidedCrystalBall);
        break;
    }//switch( t )
    
    if( !want_par )
    {
      m_coefficients[t] = 0.0;
      m_uncertainties[t] = 0.0;
      m_fitFor[t] = false;
      continue;
    }//if( !want_par )
    
    const char *label = to_string( t );
    node = peak_node->first_node(label);
    if( !node || !node->value() )
      throw runtime_error( "No coefficient " + string(label) );
    
    float dblval = 0.0f, dbluncrt = 0.0f;
    const int nread = sscanf(node->value(), "%g %g", &dblval, &dbluncrt);
    if( (nread != 1) && (nread != 2) )
      throw runtime_error( "unable to read value or uncert for " + string(label) );
    
    m_coefficients[t] = dblval;
    m_uncertainties[t] = dbluncrt;
    
    att = node->first_attribute("fit",3);
    if( !att || !att->value() )
    {
      switch( t )
      {
        case Mean: case Sigma: case GaussAmplitude:
        case SkewPar0: case SkewPar1: case SkewPar2: case SkewPar3:
          throw runtime_error( "No fit attribute for " + string(label) );
          break;
          
        case Chi2DOF: case NumCoefficientTypes:
          m_fitFor[t] = false;
          break;
      }//switch( t )
    }else
    {
      m_fitFor[t] = compare(att->value(),att->value_size(),"true",4,false);
      if( !m_fitFor[t] && !compare(att->value(),att->value_size(),"false",5,false) )
        throw runtime_error( "invalid fit value" );
    }
  }//for(...)

  
  xml_node<char> *line_color_node = peak_node->first_node("LineColor",9);
  if( line_color_node && (line_color_node->value_size() >= 7) )
  {
    //Added 20181027 without incrementing XML version since we're making it optional
    const string color = string( line_color_node->value(), line_color_node->value_size() );
    try
    {
      m_lineColor = Wt::WColor(color);
    }catch(...)
    {
      m_lineColor = Wt::WColor();
    }
  }else
  {
    m_lineColor = Wt::WColor();
  }
  
  xml_node<char> *nuc_node = peak_node->first_node("Nuclide",7);
  xml_node<char> *xray_node = peak_node->first_node("XRay",4);
  xml_node<char> *rctn_node = peak_node->first_node("Reaction",8);
  
  try
  {
  
    if( nuc_node )
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      xml_node<char> *name_node = nuc_node->first_node("Name",4);
      xml_node<char> *p_node = nuc_node->first_node("DecayParent",11);
      xml_node<char> *c_node = nuc_node->first_node("DecayChild",10);
      xml_node<char> *e_node = nuc_node->first_node("DecayGammaEnergy",16);
      xml_node<char> *type_node = nuc_node->first_node("DecayGammaType",14);
    
      const bool isNormalNucTrans = (p_node && c_node && e_node && name_node->value()
                                    && p_node->value() && c_node->value() && e_node->value());
    
      bool gotGammaType = false;
      const char *typeval = type_node->value();
      const size_t typelen = type_node->value_size();
    
      if( compare( typeval, typelen, "NormalGamma", 11, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = PeakDef::NormalGamma;
      }else if( compare( typeval, typelen, "AnnihilationGamma", 17, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = PeakDef::AnnihilationGamma;
      }else if( compare( typeval, typelen, "SingleEscapeGamma", 17, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = PeakDef::SingleEscapeGamma;
      }else if( compare( typeval, typelen, "DoubleEscapeGamma", 17, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = PeakDef::DoubleEscapeGamma;
      }else if( compare( typeval, typelen, "XrayGamma", 9, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = PeakDef::XrayGamma;
      }

      if( !name_node || !gotGammaType )
        throw runtime_error( "Invalidly specified nuclide" );
    
      m_parentNuclide = db->nuclide( name_node->value() );
    
      if( isNormalNucTrans )
      {
        const SandiaDecay::Nuclide *parent = db->nuclide( p_node->value() );
        const SandiaDecay::Nuclide *child = db->nuclide( c_node->value() );
    
        if( !m_parentNuclide )
          throw runtime_error( "Invalid nuclide name " + string(name_node->value()) );
    
        float decay_gamma_energy;
        if( sscanf( e_node->value(), "%g", &decay_gamma_energy ) != 1 )
          throw runtime_error( "Invalid nuclide gamma energy" );
    
        for( size_t i = 0; i < parent->decaysToChildren.size(); ++i )
        {
          const SandiaDecay::Transition *trans = parent->decaysToChildren[i];
          if( trans->parent==parent && trans->child==child )
          {
            size_t nearest = 0;
            double delta_eneregy = 999.9;
            for( size_t j = 0; j < trans->products.size(); ++j )
            {
              const SandiaDecay::RadParticle &particle = trans->products[j];
              if( particle.type == SandiaDecay::GammaParticle
                  || (m_sourceGammaType == PeakDef::AnnihilationGamma
                      && particle.type == SandiaDecay::PositronParticle)
                  || (m_sourceGammaType == PeakDef::XrayGamma
                     && particle.type == SandiaDecay::XrayParticle)
                 )
              {
                const double de = fabs( decay_gamma_energy - particle.energy );
                if( de < delta_eneregy )
                {
                  delta_eneregy = de;
                  nearest = j;
                }//if( de < delta_eneregy )
              }//if( particle.type == SandiaDecay::GammaParticle )
            }//for( size_tj = 0; j < trans->products.size(); ++j )
        
            if( delta_eneregy > 1.0 )
              throw std::runtime_error( "Couldnt find gamma near in energy to "
                                     + string(e_node->value()) + " keV" );
        
            m_radparticleIndex = static_cast<int>(nearest);
            m_transition = trans;
            i = parent->decaysToChildren.size();
          }//if( this transition matches )
        }//for( size_t i = 0; i < parent->decaysToChildren.size(); ++i )
      
        if( !m_transition && (m_sourceGammaType != PeakDef::AnnihilationGamma) )
        {
          if( parent && child )
            throw std::runtime_error( "Couldnt find specified transition for "
                                     + parent->symbol + " to " + child->symbol );
          else
            throw std::runtime_error( "Couldnt find specified transition" );
        }//if( !m_transition )
      }//if( isNormalNucTrans )
    }//if( nuc_node )
  
  
    if( xray_node )
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
      xml_node<char> *el_node = xray_node->first_node("Element",7);
      xml_node<char> *energy_node = xray_node->first_node("Energy",6);
  
      if( !el_node || !el_node->value() || !energy_node || !energy_node->value() )
        throw runtime_error( "Ill specified xray" );
    
      m_xrayElement = db->element( el_node->value() );
      if( !m_xrayElement )
        throw std::runtime_error( "Couldnt retrieve x-ray element" );
    
      float dummyval;
      if( sscanf( energy_node->value(), "%g", &dummyval ) != 1 )
        throw runtime_error( "non numeric xray energy" );
      m_xrayEnergy = dummyval;
    }//if( xray_node )
  
    if( rctn_node )
    {
      const ReactionGamma *rctns = ReactionGammaServer::database();
    
      xml_node<char> *name_node   = rctn_node->first_node("Name",4);
      xml_node<char> *energy_node = rctn_node->first_node("Energy",6);
      xml_node<char> *type_node   = rctn_node->first_node("Type",4);
    
      if( !name_node || !name_node->value() || !energy_node || !energy_node->value() )
        throw runtime_error( "Ill specified reaction" );
    
      float dummyval;
      if( sscanf( energy_node->value(), "%g", &dummyval ) != 1 )
        throw runtime_error( "non numeric reaction energy" );
      m_reactionEnergy = dummyval;
  
      //We will default do NormalGamma since early versions of serializtion didnt
      //  write the type
      m_sourceGammaType = PeakDef::NormalGamma;
    
      if( type_node )
      {
        const char *typeval = type_node->value();
        const size_t typelen = type_node->value_size();
  
        if( compare( typeval, typelen, "NormalGamma", 11, false ) )
          m_sourceGammaType = PeakDef::NormalGamma;
        else if( compare( typeval, typelen, "AnnihilationGamma", 17, false ) )
          m_sourceGammaType = PeakDef::AnnihilationGamma;
        else if( compare( typeval, typelen, "SingleEscapeGamma", 17, false ) )
          m_sourceGammaType = PeakDef::SingleEscapeGamma;
        else if( compare( typeval, typelen, "DoubleEscapeGamma", 17, false ) )
          m_sourceGammaType = PeakDef::DoubleEscapeGamma;
        else if( compare( typeval, typelen, "XrayGamma", 9, false ) )
          m_sourceGammaType = PeakDef::XrayGamma;
      }//if( type_node )

    
      const string reaction = name_node->value();
      std::vector<const ReactionGamma::Reaction *> candidates;
      rctns->reactions( float(m_reactionEnergy-1.0),
                        float(m_reactionEnergy+1.0), candidates );
      for( size_t i = 0; i < candidates.size(); ++i )
        if( candidates[i]->name() == reaction )
          m_reaction = candidates[i];
      //Could have called ReactionGamma::gammas( const string &name,...) as well
      //  rather than doing the manual search above
    
      if( !m_reaction )
        throw std::runtime_error( "Couldnt find reaction '" + reaction + "'" );
    }//if( rctn_node )
  }catch( std::exception &e )
  {
    m_radparticleIndex = -1;
    m_transition = NULL;
    m_sourceGammaType = NormalGamma;
    m_parentNuclide = NULL;
    
    stringstream msg;
    msg << "Failed to assign peak at " << mean() << " keV to nuclide/xray/reaction: " << e.what();
    
    if( wApp )
    {
      passMessage( msg.str(), WarningWidget::WarningMsgHigh );
    }else
    {
      cerr << msg.str() << endl;
    }
  }
}//void fromXml(...)

#if( SpecUtils_ENABLE_D3_CHART )
std::string PeakDef::gaus_peaks_to_json(const std::vector<std::shared_ptr<const PeakDef> > &peaks,
                                  const std::shared_ptr<const SpecUtils::Measurement> &foreground )
{
  //Need to check all numbers to make sure not inf or nan
  
  stringstream answer;
  
  if (peaks.empty())
    return answer.str();

  std::shared_ptr<const PeakContinuum> continuum = peaks[0]->continuum();
  if (!continuum)
    throw runtime_error("gaus_peaks_to_json: invalid continuum");
  
  
  if( IsInf(continuum->lowerEnergy()) || IsNan(continuum->lowerEnergy()) )
    throw runtime_error( "Continuum lower energy is invalid" );
  
  if( IsInf(continuum->upperEnergy()) || IsNan(continuum->upperEnergy()) )
    throw runtime_error( "Continuum upper energy is invalid" );
  
  const char *q = "\""; // for creating valid json format

  answer << "{" << q << "type" << q << ":" << q;
  switch( continuum->type() )
  {
    case PeakContinuum::NoOffset:     answer << "NoOffset";     break;
    case PeakContinuum::Constant:     answer << "Constant";     break;
    case PeakContinuum::Linear:       answer << "Linear";       break;
    case PeakContinuum::Quadratic:    answer << "Quadratic";    break;
    case PeakContinuum::Cubic:        answer << "Cubic";        break;
    case PeakContinuum::FlatStep:     answer << "FlatStep";     break;
    case PeakContinuum::LinearStep:   answer << "LinearStep";   break;
    case PeakContinuum::BiLinearStep: answer << "BiLinearStep"; break;
    case PeakContinuum::External:     answer << "External";     break;
  }//switch( continuum->type() )
  
  
  // We use the peaks defined range, and not the continuum, as the continuum may not
  //  have the range defined (but normally should).
  //  This next statement assumes peaks are sorted in increasing mean (but this only
  //  matters if ROI range is not defined) - we could check this, but maybe not
  //  worth the overhead for the edge case, that we might never encounter.
  answer << q << "," << q << "lowerEnergy" << q << ":" << peaks.front()->lowerX()
         << "," << q << "upperEnergy" << q << ":" << peaks.back()->upperX();
  
  
  if( foreground && foreground->channel_energies() && foreground->channel_energies()->size() > 2 )
  {
    // If we want integer counts (for simple spectra only), should use gamma_channels_sum, otherwise
    //   gamma_integral(...) will almost always return fractional.
    //const size_t lower_channel = foreground->find_gamma_channel( continuum->lowerEnergy() );
    //const size_t upper_channel = foreground->find_gamma_channel( continuum->upperEnergy() );
    //const double sum = foreground->gamma_channels_sum( lower_channel, upper_channel );
    const double sum = foreground->gamma_integral( continuum->lowerEnergy(), continuum->upperEnergy() );
    answer << "," << q << "roiCounts" << q << ":" << sum;
  }//if( foreground )
  
  
  switch( continuum->type() )
  {
    case PeakContinuum::NoOffset:
      break;
      
    case PeakContinuum::Constant:
    case PeakContinuum::Linear:
    case PeakContinuum::Quadratic:
    case PeakContinuum::Cubic:
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    {
      if( IsInf(continuum->referenceEnergy()) || IsNan(continuum->referenceEnergy()) )
        throw runtime_error( "Continuum reference energy is invalid" );
      
      answer << "," << q << "referenceEnergy" << q << ":" << continuum->referenceEnergy();
      const vector<double> &values = continuum->parameters();
      const vector<double> &uncerts = continuum->uncertainties();
      answer << "," << q << "coeffs" << q << ":[";
      for (size_t i = 0; i < values.size(); ++i)
      {
        if( IsInf(values[i]) || IsNan(values[i]) )
          throw runtime_error( "Continuum coef is invalid" );
        
        answer << (i ? "," : "") << values[i];
      }
      answer << "]," << q << "coeffUncerts" << q << ":[";
      for (size_t i = 0; i < uncerts.size(); ++i)
        answer << (i ? "," : "") << ((IsInf(uncerts[i]) || IsNan(uncerts[i])) ? -1.0 : uncerts[i]); //we'll let uncertainties slide since we dont use them
      answer << "]";
      
      answer << "," << q << "fitForCoeff" << q << ":[";
      for (size_t i = 0; i < continuum->fitForParameter().size(); ++i)
        answer << (i ? "," : "") << (continuum->fitForParameter()[i] ? "true" : "false");
      answer << "]";
      
      if( ((continuum->type() == PeakContinuum::FlatStep)
           || (continuum->type() == PeakContinuum::LinearStep)
           || (continuum->type() == PeakContinuum::BiLinearStep) )
         && foreground && foreground->num_gamma_channels() )
      {
        const size_t nchannel = foreground->num_gamma_channels();
        
        //We'll put in the coefficients, but also the values to make things easy in JS
        size_t firstbin = foreground->find_gamma_channel( continuum->lowerEnergy() );
        size_t lastbin = foreground->find_gamma_channel( continuum->upperEnergy() );
        
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : nchannel;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : nchannel;
        
        
        // When the JSON of the spectrum chart is defined, D3SpectrumExport.cpp/write_spectrum_data_js(...)
        //  will send the energy calibration as coefficients sometimes, and lower channel energies other
        //  times.  If coefficients are sent, then the JS computes channel bounds using doubles.
        //  And somewhat surprisingly, rounding causes visual artifacts of the continuum when
        //  the 'continuumEnergies' and 'continuumCounts' arrays are used, which are only accurate
        //  to float levels.  So we will increase accuracy of the 'answer' stream here, which appears
        //  to be enough to avoid these artifacts, but also commented out is how we could compute
        //  to double precision to match what happens in the JS.
        const auto oldprecision = answer.precision();  //probably 6 always
        answer << std::setprecision(std::numeric_limits<float>::digits10 + 1);
        
        answer << "," << q << "continuumEnergies" << q << ":[";
        for (size_t i = firstbin; i <= lastbin; ++i)
          answer << ((i!=firstbin) ? "," : "") << foreground->gamma_channel_lower(i);
        answer << "]," << q << "continuumCounts" << q << ":[";
        
        /*
         //Implementation to compute energies to double precision - doesnt appear to be necessary.
        const SpecUtils::EnergyCalType caltype = foreground->energy_calibration_model();
        if(  (caltype == SpecUtils::EnergyCalType::Polynomial
             || caltype == SpecUtils::EnergyCalType::FullRangeFraction
             || caltype == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
           && foreground->deviation_pairs().empty() )
        {
          const size_t nchannel = foreground->num_gamma_channels();
          const vector<float> &coefs = foreground->calibration_coeffs();
          const vector<pair<float,float>> &dev_pairs = foreground->deviation_pairs();
         
          answer << "," << q << "continuumEnergies" << q << ":[";
          for (size_t i = firstbin; i <= lastbin; ++i)
          {
            double energy;
            if( caltype == SpecUtils::EnergyCalType::FullRangeFraction )
              energy = SpecUtils::fullrangefraction_energy( i, coefs, nchannel, dev_pairs );
            else
              energy = SpecUtils::polynomial_energy( i, coefs, dev_pairs );
            answer << ((i!=firstbin) ? "," : "") << energy;
          }
          answer << "]," << q << "continuumCounts" << q << ":[";
        }else
        {
          answer << "," << q << "continuumEnergies" << q << ":[";
          for (size_t i = firstbin; i <= lastbin; ++i)
            answer << ((i!=firstbin) ? "," : "") << foreground->gamma_channel_lower(i);
          answer << "]," << q << "continuumCounts" << q << ":[";
        }
         */
#ifndef _WIN32        
#warning "Can make computing continuumCounts for FlatStep/LinearStep/BiLinearStep so much more efficient"
#endif        
        for (size_t i = firstbin; i <= lastbin; ++i)
        {
          const float lower_x = foreground->gamma_channel_lower( i );
          const float upper_x = foreground->gamma_channel_upper( i );
          const float cont_counts = continuum->offset_integral( lower_x, upper_x, foreground );
          answer << ((i!=firstbin) ? "," : "") << cont_counts;
        }
        answer << "]";
        
        answer << std::setprecision(9);
      }//if( continuum->type() == FlatStep/LinearStep/BiLinearStep )
      
      break;
    }//polynomial continuum
      
    case PeakContinuum::External:
    {
      if( continuum->externalContinuum()
          && continuum->externalContinuum()->num_gamma_channels() )
      {
        std::shared_ptr<const Measurement> hist = continuum->externalContinuum();
        size_t firstbin = hist->find_gamma_channel(continuum->lowerEnergy());
        size_t lastbin = hist->find_gamma_channel(continuum->upperEnergy());
        
        const size_t nchannel = hist->num_gamma_channels();
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : nchannel;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : nchannel;
        
        //see comments above in FlatStep section
        const auto oldprecision = answer.precision();
        answer << std::setprecision(std::numeric_limits<float>::digits10 + 1);

        answer << "," << q << "continuumEnergies" << q << ":[";
        for (size_t i = firstbin; i <= lastbin; ++i)
          answer << ((i!=firstbin) ? "," : "") << hist->gamma_channel_lower(i);
        answer << "]," << q << "continuumCounts" << q << ":[";
        for (size_t i = firstbin; i <= lastbin; ++i)
          answer << ((i!=firstbin) ? "," : "") << hist->gamma_channel_content(i);
        answer << "]";
        
        answer << std::setprecision(9);
      }
    }//case PeakContinuum::External:
  }//switch( continuum->type() )


  answer << "," << q << "peaks" << q << ":[";
  for (size_t i = 0; i < peaks.size(); ++i)
  {
    const PeakDef &p = *peaks[i];
    if (continuum != p.continuum())
      throw runtime_error("gaus_peaks_to_json: peaks all must share same continuum");
    answer << (i ? "," : "") << "{";

    if( !p.userLabel().empty() )
    {
      string label = p.userLabel();
      
      // TODO: Implement better JSON escaping with a custom function, or move to using nlohmann JSON (which is now in SpecUtils anyway).  See also https://stackoverflow.com/questions/7724448/simple-json-string-escape-for-c
      if( label.find_first_of( "\"\\\b\f\n\r\t" ) != string::npos ) //not exact, but close enough
      {
        Wt::Json::Object val( {{"userLabel", Wt::Json::Value( Wt::WString::fromUTF8(label) )} } );
        label = Wt::Json::serialize( val );
        SpecUtils::trim(label);
        if( label.size() && label.front() == '{' )
          label = label.substr(1);
        if( label.size() && label.back() == '}' )
          label = label.substr(0,label.size()-1);
        SpecUtils::trim(label);
        
        answer << label << ",";
      }else
      {
        answer << q << "userLabel" << q << ":" << q << label << q << ",";
      }
    }

    if (!p.lineColor().isDefault())
      answer << q << "lineColor" << q << ":" << q << p.lineColor().cssText(false) << q << ",";

    answer << q << "type" << q << ":";
    switch( p.type() )
    {
      case PeakDef::GaussianDefined:
        answer << q << "GaussianDefined" << q << ",";
      break;
        
      case PeakDef::DataDefined:
        answer << q << "DataDefined" << q << ",";
      break;
    }//switch( p.type() )

    
    double dist_norm = 0.0;
    answer << q << "skewType" << q << ":";
    
    switch( p.skewType() )
    {
      case PeakDef::NoSkew:
        answer << q << "NoSkew" << q << ",";
        break;
      
      case Bortel:
        answer << q << "ExGauss" << q << ",";
        break;
        
      case CrystalBall:
        answer << q << "CB" << q << ",";
        dist_norm = PeakDists::crystal_ball_norm( p.coefficient(CoefficientType::Sigma),
                                          p.coefficient(CoefficientType::SkewPar0),
                                      p.coefficient(CoefficientType::SkewPar1) );
        break;
        
      case DoubleSidedCrystalBall:
        answer << q << "DSCB" << q << ",";
        dist_norm = PeakDists::DSCB_norm( p.coefficient(CoefficientType::SkewPar0),
                                                   p.coefficient(CoefficientType::SkewPar1),
                                                   p.coefficient(CoefficientType::SkewPar2),
                                                   p.coefficient(CoefficientType::SkewPar3) );
        
        break;
        
      case GaussExp:
        answer << q << "GaussExp" << q << ",";
        dist_norm = PeakDists::gauss_exp_norm( p.coefficient(CoefficientType::Sigma),
                                   p.coefficient(CoefficientType::SkewPar0) );
        break;
        
      case ExpGaussExp:
        answer << q << "ExpGaussExp" << q << ",";
        dist_norm = PeakDists::exp_gauss_exp_norm( p.coefficient(CoefficientType::Sigma),
                                       p.coefficient(CoefficientType::SkewPar0),
                                       p.coefficient(CoefficientType::SkewPar1) );
        break;
    }//switch( p.type() )

    if( p.skewType() != PeakDef::NoSkew )
      answer << q << "DistNorm" << q << ":" << dist_norm << ",";
    
    if( (p.type() == PeakDef::GaussianDefined) && (p.skewType() != PeakDef::NoSkew) )
    {
      double hidden_frac = 1.0E-6; //
      
      try
      {
        pair<double,double> vis_limits;
        
        switch( p.skewType() )
        {
          case NoSkew:
            vis_limits.first = p.mean() - 5.0*p.sigma();
            vis_limits.second = p.mean() + 5.0*p.sigma();
            break;
            
          case Bortel:
            vis_limits = PeakDists::bortel_coverage_limits( p.mean(), p.sigma(),
                                                           p.coefficient(CoefficientType::SkewPar0),
                                                           hidden_frac );
            break;
            
          case GaussExp:
            vis_limits = PeakDists::gauss_exp_coverage_limits( p.mean(), p.sigma(),
                                                              p.coefficient(CoefficientType::SkewPar0),
                                                              hidden_frac );
            break;
          case CrystalBall:
            try
            {
              vis_limits = PeakDists::crystal_ball_coverage_limits( p.mean(), p.sigma(),
                                                                 p.coefficient(CoefficientType::SkewPar0),
                                                                 p.coefficient(CoefficientType::SkewPar1),
                                                                 hidden_frac );
            }catch( std::exception & )
            {
              // CB dist can have really long tail, causing the coverage limits to fail, because
              //  of unreasonable values - in this case we'll use the entire ROI.
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }
            break;
          case ExpGaussExp:
            vis_limits = PeakDists::exp_gauss_exp_coverage_limits( p.mean(), p.sigma(),
                                                                  p.coefficient(CoefficientType::SkewPar0),
                                                                  p.coefficient(CoefficientType::SkewPar1),
                                                                  hidden_frac );
            break;
          case DoubleSidedCrystalBall:
            // For largely skewed peaks, going out t 1E-6 is unreasonable, so we could limit
            //  this to 1E-3, to improve chances of success.
            //if( p.coefficient(CoefficientType::SkewPar1) < 5.0
            //   || p.coefficient(CoefficientType::SkewPar3) < 5.0)
            //  hidden_frac = 1.0E-3; //DSCB doesnt behave well for power-law peaks...
            
            try
            {
              const double mean = p.mean();
              const double sigma = p.sigma();
              const double left_skew = p.coefficient(CoefficientType::SkewPar0);
              const double left_n = p.coefficient(CoefficientType::SkewPar1);
              const double right_skew = p.coefficient(CoefficientType::SkewPar2);
              const double right_n = p.coefficient(CoefficientType::SkewPar3);
              const double p = hidden_frac;
              
              vis_limits = PeakDists::double_sided_crystal_ball_coverage_limits( mean, sigma,
                                              left_skew, left_n, right_skew, right_n, hidden_frac );
            }catch( std::exception & )
            {
              // DSCB dist can have really long tail, causing the coverage limits to fail, because
              //  of unreasonable values - in this case we'll use the entire ROI.
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }//try / catch
            
            break;
        }//switch( p.skewType() )
        
        answer << q << "visRange" << q << ":[" << vis_limits.first << "," << vis_limits.second << "],";
      }catch( std::exception &e )
      {
        cerr << "Failed to get limits for peak type " << PeakDef::to_string(p.skewType())
            << ": " << e.what() << endl
            << "\tFor peak with skew=[" << p.coefficient(CoefficientType::SkewPar0)
            << ", " << p.coefficient(CoefficientType::SkewPar1)
            << ", " << p.coefficient(CoefficientType::SkewPar2)
            << ", " << p.coefficient(CoefficientType::SkewPar3)
            << "], mean=" << p.mean() << ", and sigma=" << p.sigma()
            << " with prob=" << hidden_frac << endl;
      }//try / catch
    }//if( we should give a hint to the JS about range to draw skewed peaks )
    
    
    
    for (PeakDef::CoefficientType t = PeakDef::CoefficientType(0);
      t < PeakDef::NumCoefficientTypes; t = PeakDef::CoefficientType(t + 1))
    {
      // Dont include non-useful skew coefficients
      switch( t )
      {
        case SkewPar0:
        case SkewPar1:
        case SkewPar2:
        case SkewPar3:
        {
          double lower, upper, starting, step;
          if( !PeakDef::skew_parameter_range( p.m_skewType, t, lower, upper, starting, step ) )
            continue;
          
          break;
        }//case A Skew Parameter
          
        default:
          break;
      }//switch( t )
      
      double coef = p.coefficient(t), uncert = p.uncertainty(t);
      assert( !IsInf(coef) && !IsNan(coef) || (t == PeakDef::Chi2DOF) ); // PeakDef::Chi2DOF may be +inf, if it wasnt evaluated
      if( IsInf(coef) || IsNan(coef) )
      {
        //throw runtime_error( "Peak ceoff is inf or nan" );
        coef = 0.0;
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Peak coefficient is Inf or NaN" );
#endif//#if( PERFORM_DEVELOPER_CHECKS )
      }
      
      if( IsInf(uncert) || IsNan(uncert) )
        uncert = 0.0;
      
      answer << q << PeakDef::to_string(t) << q << ":[" << coef
        << "," << uncert << "," << (p.fitFor(t) ? "true" : "false")
        << "],";
    }//for(...)
    
    // The peak amplitude CPS, is only used for the little info box when you hover-over/tap a peak,
    //   so we'll just format the text here; may change in the future.
    if( (p.type() == PeakDef::GaussianDefined) && (p.amplitudeUncert() > 0.0f)
       && foreground && (foreground->live_time() > 0.0f) )
    {
      const float lt = foreground->live_time();
      const float cps = p.amplitude() / lt;
      const float cpsUncert = p.amplitudeUncert() / lt;
      const string uncertstr = PhysicalUnits::printValueWithUncertainty( cps, cpsUncert, 4 );
      answer << q << "cpsTxt" << q << ":" << q << uncertstr << q << ",";
    }//if( we can give peak CPS )
    
    answer << q << "useForEnergyCalibration" << q << ":" << (p.useForEnergyCalibration() ? "true" : "false")
           << "," << q << "forSourceFit" << q << ":" << (p.useForShieldingSourceFit() ? "true" : "false");

    if( p.useForDrfIntrinsicEffFit() != PeakDef::sm_defaultUseForDrfIntrinsicEffFit )
      answer << "," << q << "useForDrfIntrinsicEffFit" << q << ":" << (p.useForDrfIntrinsicEffFit() ? "true" : "false");
    
    if( p.useForDrfFwhmFit() != PeakDef::sm_defaultUseForDrfFwhmFit )
      answer << "," << q << "useForDrfFwhmFit" << q << ":" << (p.useForDrfFwhmFit() ? "true" : "false");
    
    if( p.useForDrfDepthOfInteractionFit() != PeakDef::sm_defaultUseForDrfDepthOfInteractionFit)
      answer << "," << q << "useForDrfDepthOfInteractionFit" << q << ":" << (p.useForDrfDepthOfInteractionFit() ? "true" : "false");
    
    
    /*
    const char *gammaTypeVal = 0;
    switch( p.sourceGammaType() )
    {
      case PeakDef::NormalGamma:       gammaTypeVal = "NormalGamma";       break;
      case PeakDef::AnnihilationGamma: gammaTypeVal = "AnnihilationGamma"; break;
      case PeakDef::SingleEscapeGamma: gammaTypeVal = "SingleEscapeGamma"; break;
      case PeakDef::DoubleEscapeGamma: gammaTypeVal = "DoubleEscapeGamma"; break;
      case PeakDef::XrayGamma:         gammaTypeVal = "XrayGamma";         break;
    }//switch( p.sourceGammaType() )

    if( p.parentNuclide() || p.xrayElement() || p.reaction() )
      answer << "," << q << "sourceType" << q << ":" << q << gammaTypeVal << q;
    */
    
    const auto srcTypeJSON = [q]( const PeakDef::SourceGammaType gamtype ) -> std::string{
      const char *type = "";
      switch( gamtype )
      {
        case NormalGamma:       return ""; break;
        case AnnihilationGamma: type = "annih."; break;
        case XrayGamma:         type = "x-ray";  break;
        case DoubleEscapeGamma: type = "D.E.";   break;
        case SingleEscapeGamma: type = "S.E.";   break;
      }//switch( p.sourceGammaType() )
      
      return string(",") + q + "type" + q + ":" + q + type + q;
    };//srcTypeJSON lamda
    
    
    if( p.parentNuclide() )
    {
      const SandiaDecay::Transition *trans = p.nuclearTransition();
      const SandiaDecay::RadParticle *decayPart = p.decayParticle();

      answer << "," << q << "nuclide" << q << ": { " << q << "name" << q << ": " << q << p.parentNuclide()->symbol << q;
      if( trans && decayPart )
      {
        string transistion_parent, decay_child;
        const SandiaDecay::Nuclide *trans_parent = trans->parent;
        transistion_parent = trans_parent->symbol;
        if( trans->child )
          decay_child = trans->child->symbol;

        answer << "," << q << "decayParent" << q << ":" << q << transistion_parent << q;
        answer << "," << q << "decayChild" << q << ":" << q << decay_child << q;
        const float energy = (IsInf(decayPart->energy) || IsNan(decayPart->energy)) ? 0.0f : decayPart->energy;
        answer << "," << q << "energy" << q << ":" << energy << "";
        answer << srcTypeJSON( p.sourceGammaType() );
      }else if( p.sourceGammaType() == AnnihilationGamma )
      {
        answer << "," << q << "energy" << q << ":511.006";
        answer << srcTypeJSON( p.sourceGammaType() );
      }//if( m_transition )

      answer << "}";
    }//if( m_parentNuclide )

    if( p.xrayElement() )
    {
      const float energy = (IsInf(p.xrayEnergy()) || IsNan(p.xrayEnergy())) ? 0.0f : p.xrayEnergy();
        
      answer << "," << q << "xray" << q
        << ":{"
             << q << "name" << q << ":" << q << p.xrayElement()->name << q << ","
             << q << "energy" << q << ":" << energy
             << srcTypeJSON( p.sourceGammaType() )
        << "}";
    }//if( m_xrayElement )


    if( p.reaction() )
    {
      const float energy = (IsInf(p.reactionEnergy()) || IsNan(p.reactionEnergy())) ? 0.0f : p.reactionEnergy();
      
      answer << "," << q << "reaction" << q
        << ":{"
             << q << "name" << q << ":" << q << p.reaction()->name() << q << ","
             << q << "energy" << q << ":" << energy
             << srcTypeJSON( p.sourceGammaType() )
        << "}";
    }//if( m_reaction )

    answer << "}";
  }//for( size_t i = 0; i < peaks.size(); ++i )

  answer << "]}";

  return answer.str();
}//std::string PeakDef::gaus_peaks_to_json(...)


string PeakDef::peak_json(const vector<std::shared_ptr<const PeakDef> > &inpeaks,
                          const std::shared_ptr<const SpecUtils::Measurement> &foreground )
{
  if (inpeaks.empty())
    return "[]";

  typedef std::map< std::shared_ptr<const PeakContinuum>, vector<std::shared_ptr<const PeakDef> > > ContinuumToPeakMap_t;

  ContinuumToPeakMap_t continuumToPeaks;
  for (size_t i = 0; i < inpeaks.size(); ++i)
    continuumToPeaks[inpeaks[i]->continuum()].push_back(inpeaks[i]);

  string json = "[";
  for (const ContinuumToPeakMap_t::value_type &vt : continuumToPeaks)
    json += ((json.size()>2) ? "," : "") + gaus_peaks_to_json(vt.second,foreground);

  json += "]";
  return json;
}//string peak_json( inpeaks )
#endif //#if( SpecUtils_ENABLE_D3_CHART )

std::shared_ptr<PeakContinuum> PeakDef::continuum()
{
  return m_continuum;
}

std::shared_ptr<const PeakContinuum> PeakDef::continuum() const
{
  return m_continuum;
}

std::shared_ptr<PeakContinuum> PeakDef::getContinuum()
{
  return m_continuum;
}

void PeakDef::setContinuum( std::shared_ptr<PeakContinuum> continuum )
{
  if( !continuum )
    throw runtime_error( "PeakDef::setContinuum(...): invalid input" );
  m_continuum = continuum;
}//void setContinuum(...)


void PeakDef::makeUniqueNewContinuum()
{
  m_continuum = std::make_shared<PeakContinuum>(*m_continuum);
}


//The below should in principle take care of gaussian area and the skew area
double PeakDef::peakArea() const
{
  return m_coefficients[PeakDef::GaussAmplitude];
}//double peakArea() const


double PeakDef::peakAreaUncert() const
{
  return m_uncertainties[PeakDef::GaussAmplitude];
}//double peakAreaUncert() const


void PeakDef::setPeakArea( const double a )
{
  m_coefficients[PeakDef::GaussAmplitude] = a;
}//void PeakDef::setPeakArea( const double a )


void PeakDef::setPeakAreaUncert( const double uncert )
{
  m_uncertainties[PeakDef::GaussAmplitude] = uncert;
}//void setPeakAreaUncert( const double a )


double PeakDef::areaFromData( std::shared_ptr<const Measurement> data ) const
{
  double sumval = 0.0;
  if( !data || !m_continuum || !m_continuum->energyRangeDefined() )
    return sumval;
  
  try
  {
    const float energyStart = (float)m_continuum->lowerEnergy();
    const float energyEnd = (float)m_continuum->upperEnergy();
    
    const size_t lower_channel = data->find_gamma_channel( energyStart );
    const size_t upper_channel = data->find_gamma_channel(energyEnd);
    
    for( size_t i = lower_channel; i <= upper_channel; ++i )
    {
      const float e0 = std::max( energyStart, data->gamma_channel_lower(i) );
      const float e1 = std::min( energyEnd, data->gamma_channel_upper(i) );
      const double data_area_i = data->gamma_integral(e0, e1);
      const double cont_area_1 = m_continuum->offset_integral(e0, e1, data);
      if( data_area_i > cont_area_1 )
        sumval += (data_area_i - cont_area_1);
    }
  }catch( std::exception &e )
  {
    cerr << "Caught exception in PeakDef::areaFromData(): " << e.what() << endl;
    return 0.0;
  }
  
  return sumval;
}//double areaFromData( std::shared_ptr<const Measurement> data ) const;

bool PeakDef::lessThanByMean( const PeakDef &lhs, const PeakDef &rhs )
{
  return (lhs.m_coefficients[PeakDef::Mean] < rhs.m_coefficients[PeakDef::Mean]);
}//lessThanByMean(...)

bool PeakDef::lessThanByMeanShrdPtr( const std::shared_ptr<const PeakDef> &lhs,
                                  const std::shared_ptr<const PeakDef> &rhs )
{
  if( !lhs || !rhs )
    return (lhs < rhs);
  return lessThanByMean( *lhs, *rhs );
}


bool PeakDef::causilyDisconnected( const PeakDef &lower_peak,
                                    const PeakDef &upper_peak,
                                    const double ncausality,
                                   const bool useRoiAsWell )
{
  if( lower_peak.continuum() == upper_peak.continuum() )
    return false;
  
  double lower_upper( 0.0 ), upper_lower( 0.0 );

  if( lower_peak.mean() < upper_peak.mean() )
  {
    lower_upper = lower_peak.gausPeak() ? lower_peak.mean() + ncausality*lower_peak.sigma() : lower_peak.upperX();
    upper_lower = upper_peak.gausPeak() ? upper_peak.mean() - ncausality*upper_peak.sigma() : upper_peak.lowerX();
    if( useRoiAsWell )
    {
      lower_upper = std::max( lower_upper, lower_peak.upperX() );
      upper_lower = std::min( upper_lower, upper_peak.lowerX() );
    }
  }else
  {
    lower_upper = upper_peak.gausPeak() ? upper_peak.mean() + ncausality*upper_peak.sigma() : upper_peak.upperX();
    upper_lower = lower_peak.gausPeak() ? lower_peak.mean() - ncausality*lower_peak.sigma() : lower_peak.lowerX();
    if( useRoiAsWell )
    {
      lower_upper = std::max( lower_upper, upper_peak.upperX() );
      upper_lower = std::min( upper_lower, lower_peak.lowerX() );
    }
  }//if( lower_peak.mean() < upper_peak.mean() ) / else

  return (upper_lower > lower_upper);
}//bool PeakDef::causilyDisconnected


bool PeakDef::causilyConnected( const PeakDef &lower_peak,
                                   const PeakDef &upper_peak,
                                   const double ncausality,
                                   const bool useRoiAsWell )
{
  return !causilyDisconnected( lower_peak, upper_peak, ncausality, useRoiAsWell );
}


bool PeakDef::ageFitNotAllowed( const SandiaDecay::Nuclide *nuc )
{
  if( !nuc || nuc->decaysToStableChildren() )
    return true;
  
  //now check for cases like Cs137 where the isotope reaches prompt and
  //  secular equilibrium very quickly (half life for these less than a day)
  //  and these time spans are less than the parents half life
  const double hl = nuc->halfLife;
  const double prompt = nuc->promptEquilibriumHalfLife();
  const double secular = nuc->secularEquilibriumHalfLife();
  
  //simpleFast: will catch cases like Cs137 where the element decays to
  //  very short lived descendants
  const bool simpleFast = ( (prompt > DBL_MIN)    //can obtain prompt equilibrium
                           && (prompt < 0.01*hl) //prompt half life less than 1% parents half life
                           && (secular < hl)     //can obtain secular equilibrium
                           && (prompt < 86400.0*SandiaDecay::second) //half lives for both of these is less than a day (an arbitrary time period)
                           && (secular < 86400.0*SandiaDecay::second) );
  
  if( simpleFast )
    return true;
  
  
  //catch cases like W187 who none of its descendants give off gammas
  for( const SandiaDecay::Transition *t : nuc->decaysToChildren )
    if( gives_off_gammas( t->child ) )
      return false;
  
  //if( prompt && (decsendnats of promp nuclide dont give off gammas) )
  //  return true;
  //return false;
  
  return true;
  /*
   //Check for cases like W187, whose none of its prodigeny give off gammas
   const vector<const SandiaDecay::Transition *> &decays = nuc->decaysToChildren;
   
   
   for( size_t decN = 0; decN < decays.size(); ++decN )
   {
   const SandiaDecay::Transition *trans = decays[decN];
   if( trans->child )
   {
   }
   }//for( size_t decN = 0; decN < decays.size(); ++decN )
   
   
   return false;
   */
}//bool ageFitNotAllowed( const SandiaDecay::Nuclide *nuc )


double PeakDef::defaultDecayTime( const SandiaDecay::Nuclide *nuclide, string *stranswer )
{
  const char *decayTimeStr = nullptr;
  double decaytime = 0;
  if( nuclide->canObtainSecularEquilibrium() )
  {
    decaytime = 10.0*nuclide->secularEquilibriumHalfLife();
  }else
  {
    decaytime = 7.0 * nuclide->promptEquilibriumHalfLife();
  }
  
  if( decaytime <= 0.0 || decaytime < 0.5*nuclide->halfLife )
  {
    decaytime = 2.5 * nuclide->halfLife;
    decayTimeStr = "2.5 HL";
  }
  
  if( nuclide->decaysToStableChildren() )
  {
    decaytime = 0.0;
    decayTimeStr = "0 s";
  }
  
  // If half-life is great than 100 years, or if U or Pu isotope, and half-life is greater than 2
  //  years, set default decay time to 20 years.
  if( (nuclide->halfLife > 100.0*SandiaDecay::year)
     || ( (nuclide->halfLife > 2.0*SandiaDecay::year)
         && ((nuclide->atomicNumber == 92) || (nuclide->atomicNumber == 94))) )
  {
    decaytime = 20.0 * SandiaDecay::year;
    decayTimeStr = "20 y";
  }
  
  //I *think* promptEquilibriumHalfLife() can maybe give a large value for an
  //  isotope (although I dont know of any examples of this) with a small half
  //  life, so we'll protect against it.
  if( decaytime > 100.0*nuclide->halfLife )
  {
    decaytime = 7.0 * nuclide->halfLife;
    decayTimeStr = "7 HL";
  }
  
  if( stranswer )
  {
    if( decayTimeStr )
      *stranswer = decayTimeStr;
    else
      *stranswer = PhysicalUnitsLocalized::printToBestTimeUnits( decaytime, 2, SandiaDecay::second );
  }//if( stranswer )
  
  return decaytime;
}//string defaultDecayTime(..)


void PeakDef::findNearestPhotopeak( const SandiaDecay::Nuclide *nuclide,
                                     const double energy,
                                     const double windowHalfWidth,
                                     const bool xraysOnly,
                                     const SandiaDecay::Transition *&transition,
                                     size_t &transition_index,
                                     SourceGammaType &sourceGammaType )
{
  transition = NULL;
  transition_index = 0;
  sourceGammaType = PeakDef::NormalGamma;
  
  if( !nuclide )
    return;
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuclide,1.0) );
  
  const double decaytime = defaultDecayTime( nuclide );
  
  vector<SandiaDecay::EnergyRatePair> gammas
  = mixture.gammas( decaytime,
                   SandiaDecay::NuclideMixture::OrderByAbundance, true );
  
  if( xraysOnly )
    gammas = mixture.xrays( decaytime, SandiaDecay::NuclideMixture::OrderByAbundance );
  
  if( gammas.empty() )
    return;
  
  map<const SandiaDecay::Transition *, vector<size_t> > ec_trans;
  
  
  double best_delta_e = 99999.9;
  SandiaDecay::EnergyRatePair nearest_gamma = gammas[0];
  for( const SandiaDecay::EnergyRatePair &gamma : gammas )
  {
    const double delta_e = fabs( gamma.energy - energy );
    if( delta_e < best_delta_e )
    {
      best_delta_e = delta_e;
      nearest_gamma = gamma;
    }//if( delta_e < best_delta_e )
  }//for( const SandiaDecay::EnergyRatePair &gamma : gammas )
  
  //loop over the decays and find the gamma nearest 'energy'
  best_delta_e = 99999.9;
  
  const vector<SandiaDecay::NuclideActivityPair> activities
                                                = mixture.activity( decaytime );
  double minEnergy = energy - windowHalfWidth;
  double maxEnergy = energy + windowHalfWidth;
  if( windowHalfWidth <= 0.0 )
  {
    minEnergy = -std::numeric_limits<double>::max();
    maxEnergy = std::numeric_limits<double>::max();
  }//if( windowHalfWidth <= 0.0 )
  
  
  double max_intensity = 0.0;
  for( const SandiaDecay::NuclideActivityPair &activity : activities )
  {
    for( const SandiaDecay::Transition *trans : activity.nuclide->decaysToChildren )
    {
      const vector<SandiaDecay::RadParticle> &products = trans->products;
      for( size_t index = 0; index < products.size(); ++index )
      {
        const SandiaDecay::RadParticle &product = products[index];
        
        if( xraysOnly && (product.type != SandiaDecay::XrayParticle) )
          continue;
        
        float energy;
        if( product.type == SandiaDecay::GammaParticle || product.type == SandiaDecay::XrayParticle )
          energy = product.energy;
        else if( product.type == SandiaDecay::PositronParticle )
          energy = 510.998910f;
        else
          continue;
        
        if( energy > maxEnergy || energy < minEnergy )
          continue;
        
        double intensity = activity.activity * trans->branchRatio * product.intensity;
        max_intensity = max( max_intensity, intensity );
      }//for( size_t index = 0; index < products.size(); ++index )
    }//for( const SandiaDecay::Transition *trans : activity.nuclide->decaysToChildren )
  }//for( SandiaDecay::NuclideActivityPair activity : activities )
  
  double positronfrac = 0.0;
  
  for( const SandiaDecay::NuclideActivityPair &activity : activities )
  {
    for( const SandiaDecay::Transition *trans : activity.nuclide->decaysToChildren )
    {
      const vector<SandiaDecay::RadParticle> &products = trans->products;
      for( size_t index = 0; index < products.size(); ++index )
      {
        const SandiaDecay::RadParticle &product = products[index];
        
        if( product.type == SandiaDecay::PositronParticle )
        {
          const double energy = 510.998910;
          if( energy > maxEnergy || energy < minEnergy )
            continue;
          
          double intensity = activity.activity * trans->branchRatio * product.intensity;
          positronfrac += intensity / max_intensity;
          ec_trans[trans].push_back( index );
        }else if( product.type == SandiaDecay::GammaParticle || product.type == SandiaDecay::XrayParticle )
        {
          if( product.energy > maxEnergy || product.energy < minEnergy )
            continue;
          
          if( xraysOnly && (product.type != SandiaDecay::XrayParticle) )
            continue;
          
          double intensity = activity.activity * trans->branchRatio * product.intensity;
          const double fracIntensity = intensity / max_intensity;
          
          // Previous to 20220617, a value of 1E-10 was used for the minimum BR, however, this fails
          //  for some non-user defined peaks (specifically in the Relative Activity calculations);
          //  so when search window is not being used (i.e., probably an automated process), the
          //  value is now small enough it seems to work well, but I also decreased the value when
          //  there is a window being used
          //  However, we could probably just eliminate this check, as all it really is for is to
          //  reduce the chances of using a gamma with a small BR, that just happens to be really
          //  close to our wanted energy.
          const double minRelativeBr = (windowHalfWidth <= 0.0) ? 1.0E-20 : 1.0E-12;
          
          if( fracIntensity < minRelativeBr )
            continue;
          
          if( fracIntensity <= 0.0 )
            continue;
          
          const double delta_e = fabs( product.energy - energy );
          double scaleDeltaE = (0.1*windowHalfWidth + delta_e) / fracIntensity;
          if( windowHalfWidth <= 0.0 )
            scaleDeltaE = delta_e;
          
          if( scaleDeltaE <= best_delta_e )
          {
            best_delta_e       = scaleDeltaE;
            transition         = trans;
            transition_index   = index;
            sourceGammaType = ((product.type==SandiaDecay::GammaParticle)
                                ? PeakDef::NormalGamma : PeakDef::XrayGamma);
          }//if( delta_e < best_delta_e )
        }//if( positron ) / else if( gamma )
      }//for( const SandiaDecay::RadParticle &particle : products )
    }//for( const SandiaDecay::Transition *trans : decays)
  }//for( SandiaDecay::NuclideActivityPair activity : activities )
  
  if( positronfrac > 0.0 )
  {
    const double delta_e = fabs( 510.998910 - energy );
    double scaleDeltaE = (0.1*windowHalfWidth + delta_e) / positronfrac;
    if( windowHalfWidth <= 0.0 )
      scaleDeltaE = delta_e;
    
    if( scaleDeltaE <= best_delta_e )
    {
      sourceGammaType = PeakDef::AnnihilationGamma;
      
      if( ec_trans.size()==1 && ec_trans.begin()->second.size()==1 )
      {
        transition         = ec_trans.begin()->first;
        transition_index   = ec_trans.begin()->second[0];
      }else
      {
        transition = NULL;
        transition_index = 0;
      }//if( there is only one EC decay ) / else
    }//if( delta_e < best_delta_e )
  }//if( positronfrac > 0.0 )
  
}//void findNearestPhotopeak(...)



const SandiaDecay::EnergyIntensityPair *PeakDef::findNearestXray(
                                const SandiaDecay::Element *el, double energy )
{
  if( !el )
    return NULL;
  
  const size_t nxray = el->xrays.size();
  size_t nearest = 0;
  double nearestDE = 999999.9;
  
  for( size_t i = 0; i < nxray; ++i )
  {
    const double thisDE = fabs( el->xrays[i].energy - energy );
    if( thisDE < nearestDE )
    {
      nearest = i;
      nearestDE = thisDE;
    }//if( thisDE < nearestE )
  }//for( size_t i = 0; i < nxray; ++i )
  
  if( nearestDE > 999999.0 )
    return NULL;
  
  return &(el->xrays[nearest]);
}//double findNearestXray( const SandiaDecay::Element *el, double energy )


/*
PeakDef::PeakDef( const PeakDef &rhs )
{
  *this = rhs;
}

const PeakDef &PeakDef::operator=( const PeakDef &rhs )
{
  if( (&rhs) == this )
    return *this;

  m_type       = rhs.m_type;
  m_offsetType = rhs.m_offsetType;

  m_parentNuclide            = rhs.m_parentNuclide;
  m_transition               = rhs.m_transition;
  m_radparticleIndex         = rhs.m_radparticleIndex;
  m_useForEnergyCal        = rhs.m_useForEnergyCal;
  m_useForShieldingSourceFit = rhs.m_useForShieldingSourceFit;
  m_useForManualRelEff = rhs.m_useForManualRelEff;
 
  m_useForDrfIntrinsicEffFit = rhs.m_useForDrfIntrinsicEffFit;
  m_useForDrfFwhmFit = rhs.m_useForDrfFwhmFit;
  m_useForDrfDepthOfInteractionFit = rhs.m_useForDrfDepthOfInteractionFit;
 
  const size_t nbytes = sizeof(m_coefficients[0])*NumCoefficientTypes;
  memcpy( &(m_coefficients[0]), &(rhs.m_coefficients[0]), nbytes );
  memcpy( &(m_uncertainties[0]), &(rhs.m_uncertainties[0]), nbytes );
}//const PeakDef &operator=( const PeakDef &rhs )
*/

bool PeakDef::operator==( const PeakDef &rhs ) const
{
  for( CoefficientType t = CoefficientType(0);
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    if( m_coefficients[t] != rhs.m_coefficients[t] )
      return false;
  }

  return m_type==rhs.m_type
         && (*m_continuum == *rhs.m_continuum)
      && m_parentNuclide==rhs.m_parentNuclide
      && m_transition==rhs.m_transition
      && m_sourceGammaType==rhs.m_sourceGammaType
      && m_radparticleIndex==rhs.m_radparticleIndex
      && m_useForEnergyCal==rhs.m_useForEnergyCal
      && m_useForShieldingSourceFit==rhs.m_useForShieldingSourceFit
      && m_useForManualRelEff==rhs.m_useForManualRelEff
      && m_useForDrfIntrinsicEffFit == rhs.m_useForDrfIntrinsicEffFit
      && m_useForDrfFwhmFit == rhs.m_useForDrfFwhmFit
      && m_useForDrfDepthOfInteractionFit == rhs.m_useForDrfDepthOfInteractionFit
      && m_lineColor==rhs.m_lineColor
      ;
}//PeakDef::operator==


void PeakDef::clearSources()
{
  m_radparticleIndex = -1;
  m_parentNuclide = nullptr;
  m_transition = nullptr;
  m_xrayElement = nullptr;
  m_reaction = nullptr;
  m_sourceGammaType = PeakDef::NormalGamma;
  m_xrayEnergy = m_reactionEnergy = 0.0f;
}//void PeakDef::clearSources()


bool PeakDef::hasSourceGammaAssigned() const
{
  return ((m_parentNuclide && m_transition && m_radparticleIndex>=0) 
    || (m_parentNuclide && (m_sourceGammaType == PeakDef::SourceGammaType::AnnihilationGamma))
    || m_xrayElement || m_reaction);
}

void PeakDef::setNuclearTransition( const SandiaDecay::Nuclide *parentNuclide,
                                    const SandiaDecay::Transition *transition,
                                    const int index,
                                    const SourceGammaType sourceType )
{
  m_transition = transition;
  m_sourceGammaType = sourceType;
  
  if( m_transition && (index >= 0) && (index < static_cast<int>(m_transition->products.size())) )
    m_radparticleIndex = index;
  else
    m_radparticleIndex = -1;
  m_parentNuclide = parentNuclide;

  if( m_transition || sourceType==PeakDef::AnnihilationGamma )
  {
    m_xrayElement = NULL;
    m_xrayEnergy = 0.0;
    m_reaction = NULL;
    m_reactionEnergy = 0.0;
  }//if( m_transition )
}//void setNuclearTransition( SandiaDecay::Transition *transition, int radParticle )


void PeakDef::setXray( const SandiaDecay::Element *el, const float energy )
{
  m_xrayElement = el;
  m_xrayEnergy = energy;
  
  if( el )
  {
    m_radparticleIndex = -1;
    m_parentNuclide = NULL;
    m_transition = NULL;
    m_reaction = NULL;
    m_reactionEnergy = 0.0f;
    m_sourceGammaType = PeakDef::NormalGamma;
  }//if( el )
}//void setXray( const SandiaDecay::Element *el, const float energy )


void PeakDef::setReaction( const ReactionGamma::Reaction *rctn,
                           const float energy,
                           const SourceGammaType sourceType )
{
  m_reaction = rctn;
  m_reactionEnergy = energy;
  
  if( rctn )
  {
    m_radparticleIndex = -1;
    m_parentNuclide = NULL;
    m_transition = NULL;
    m_xrayElement = NULL;
    m_xrayEnergy = 0.0f;
    if( sourceType == PeakDef::AnnihilationGamma )
      throw runtime_error( "PeakDef::setReaction can not be called with source"
                           " type PeakDef::AnnihilationGamma" );
    m_sourceGammaType = sourceType;
  }//if( rctn )
}//void setReaction( const ReactionGamma::Reaction *rctn, double energy )



const SandiaDecay::Transition *PeakDef::nuclearTransition() const
{
  return m_transition;
}//const Transition *nuclearTransition() const


const SandiaDecay::Nuclide *PeakDef::parentNuclide() const
{
  return m_parentNuclide;
}//const SandiaDecay::Nuclide *parentNuclide() const;

int PeakDef::decayParticleIndex() const
{
  return m_radparticleIndex;
}

const SandiaDecay::RadParticle *PeakDef::decayParticle() const
{
  if( !m_transition || (m_radparticleIndex < 0)
      || (m_radparticleIndex>=static_cast<int>(m_transition->products.size())) )
    return NULL;

  return &(m_transition->products[m_radparticleIndex]);
}//const RadParticle *decayParticle() const


PeakDef::SourceGammaType PeakDef::sourceGammaType() const
{
  return m_sourceGammaType;
}

const SandiaDecay::Element *PeakDef::xrayElement() const
{
  return m_xrayElement;
}

float PeakDef::xrayEnergy() const
{
  return m_xrayEnergy;
}

const ReactionGamma::Reaction *PeakDef::reaction() const
{
  return m_reaction;
}

float PeakDef::reactionEnergy() const
{
  return m_reactionEnergy;
}


float PeakDef::gammaParticleEnergy() const
{
  if( m_parentNuclide )
  {
    switch( m_sourceGammaType )
    {
      case NormalGamma:
      case XrayGamma:
        if( m_transition && (m_radparticleIndex >= 0) )
          return m_transition->products.at(m_radparticleIndex).energy;
      break;
      case AnnihilationGamma:
        return static_cast<float>( 510.99891*SandiaDecay::keV );
      case SingleEscapeGamma:
        if( m_transition && (m_radparticleIndex >= 0) )
          return m_transition->products.at(m_radparticleIndex).energy - 510.9989f;
      break;
      case DoubleEscapeGamma:
        if( m_transition && (m_radparticleIndex >= 0) )
          return m_transition->products.at(m_radparticleIndex).energy - 2.0f*510.9989f;
      break;
    }
  }//if( m_parentNuclide )
  
  if( m_xrayElement )
    return m_xrayEnergy;
  
  if( m_reaction )
  {
    switch( m_sourceGammaType )
    {
      case NormalGamma:
      case XrayGamma:
        return m_reactionEnergy;
        
      case AnnihilationGamma:
        return static_cast<float>( 510.99891*SandiaDecay::keV );
        
      case SingleEscapeGamma:
        return m_reactionEnergy - 510.9989f;
        
      case DoubleEscapeGamma:
        return m_reactionEnergy - 2.0f*510.9989f;
    }//switch( m_sourceGammaType )
  }//if( m_reaction )
  
  throw runtime_error( "Peak doesnt have a gamma associated with it" );
  
  return 0.0f;
}//float gammaParticleEnergy() const


const std::vector< PeakDef::CandidateNuclide > &PeakDef::candidateNuclides() const
{
  return m_candidateNuclides;
}//const std::vector< CandidateNuclide > &candidateNuclides() const


void PeakDef::addCandidateNuclide( const PeakDef::CandidateNuclide &candidate )
{
  m_candidateNuclides.push_back( candidate );
}//void addCandidateNuclide( const CandidateNuclide &candidate )


void PeakDef::addCandidateNuclide( const SandiaDecay::Nuclide * const nuclide,
                                   const SandiaDecay::Transition * const transition,
                                   const int radparticleIndex,
                                   const SourceGammaType sourceGammaType,
                                   const float weight )
{
  CandidateNuclide candidate;
  candidate.nuclide = nuclide;
  candidate.transition = transition;
  candidate.radparticleIndex = radparticleIndex;
  candidate.weight = weight;
  candidate.sourceGammaType = sourceGammaType;
  m_candidateNuclides.push_back( candidate );
}//void addCandidateNuclide(...)

void PeakDef::setCandidateNuclides( const std::vector<PeakDef::CandidateNuclide> &cand )
{
  m_candidateNuclides = cand;
}//void setCandidateNuclides( std::vector<const CandidateNuclide> &candidates )



void PeakDef::inheritUserSelectedOptions( const PeakDef &parent,
                                          const bool inheritNonFitForValues )
{
  if( parent.m_parentNuclide )
  {
    setNuclearTransition( parent.m_parentNuclide, parent.m_transition,
                          parent.m_radparticleIndex, parent.m_sourceGammaType );
  }else if( parent.m_xrayElement )
  {
    setXray( parent.m_xrayElement, parent.m_xrayEnergy );
  }else if( parent.m_reaction )
  {
    setReaction( parent.m_reaction, parent.m_reactionEnergy,
                 parent.m_sourceGammaType );
  }else
  {
    m_reaction = NULL;
    m_transition = NULL;
    m_xrayElement = NULL;
    m_parentNuclide = NULL;
    m_radparticleIndex = -1;
    m_sourceGammaType = NormalGamma;
    m_xrayEnergy = m_reactionEnergy = 0.0;
  }//if( parent nuclide ) / else / else
  
  m_lineColor = parent.lineColor();
  
  m_userLabel = parent.m_userLabel;
  m_useForEnergyCal = parent.m_useForEnergyCal;
  m_useForShieldingSourceFit = parent.m_useForShieldingSourceFit;
  m_useForManualRelEff = parent.m_useForManualRelEff;
  m_useForDrfIntrinsicEffFit = parent.m_useForDrfIntrinsicEffFit;
  m_useForDrfFwhmFit = parent.m_useForDrfFwhmFit;
  m_useForDrfDepthOfInteractionFit = parent.m_useForDrfDepthOfInteractionFit;
  
  
  for( CoefficientType t = CoefficientType(0);
      t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    m_fitFor[t] = parent.m_fitFor[t];
    
    if( inheritNonFitForValues && !m_fitFor[t] )
    {
      switch( t )
      {
        case PeakDef::Mean:             case PeakDef::Sigma:
        case PeakDef::GaussAmplitude:
        case PeakDef::SkewPar0: case PeakDef::SkewPar1:
        case PeakDef::SkewPar2: case PeakDef::SkewPar3:
          m_coefficients[t]  = parent.m_coefficients[t];
          m_uncertainties[t] = parent.m_uncertainties[t];
        break;
        
        case PeakDef::Chi2DOF: case PeakDef::NumCoefficientTypes:
          break;
      }//switch( t )
    }//if( inheritNonFitForValues )
  }//for( loop over PeakDef::CoefficientType )
  
  
  const auto &rhs_cont = parent.m_continuum;
  if( m_continuum && rhs_cont )
  {
    //Currently not copying the following values of continuum:
    //double m_lowerEnergy, m_upperEnergy;
    //double m_referenceEnergy;
    //std::vector<double> m_values, m_uncertainties;
    
    const vector<bool> rhs_fit_fors = rhs_cont->fitForParameter();
    const vector<bool> orig_fit_fors = m_continuum->fitForParameter();
    for( size_t i = 0; (i < rhs_fit_fors.size()) && (i < orig_fit_fors.size()); ++i )
    {
      m_continuum->setPolynomialCoefFitFor( i, rhs_fit_fors[i] );
    }
    
    // If we wanted to copy values of coefficients, we could use:
    /*
    if( (m_continuum->type() == rhs_cont->type())
       && (m_continuum->type() != PeakContinuum::OffsetType::External ) )
    {
      switch( m_continuum->type() )
      {
        case PeakContinuum::NoOffset:
          break;
        case PeakContinuum::External:
          //if( inheritNonFitForValues )
          //  m_continuum->setExternalContinuum( rhs_cont->externalContinuum() );
          break;
          
        case PeakContinuum::Constant:
        case PeakContinuum::Linear:
        case PeakContinuum::Quadratic:
        case PeakContinuum::Cubic:
        case PeakContinuum::FlatStep:
        case PeakContinuum::LinearStep:
        case PeakContinuum::BiLinearStep:
        {
          const vector<bool> rhs_fit_fors = rhs_cont->fitForParameter();
          assert( rhs_fit_fors.size() == m_continuum->fitForParameter().size() );
          
          for( size_t i = 0; i < rhs_fit_fors.size(); ++i )
          {
            m_continuum->setPolynomialCoefFitFor( i, rhs_fit_fors[i] );
          
            if( inheritNonFitForValues && !rhs_fit_fors[i] )
            {
              const vector<double> &rhs_pars = rhs_cont->parameters();
              const vector<double> &rhs_uncerts = rhs_cont->uncertainties();
              assert( rhs_pars.size() == rhs_fit_fors.size() );
              assert( rhs_pars.size() == rhs_uncerts.size() );
              
              m_continuum->setPolynomialCoef( i, rhs_pars[i] );
              m_continuum->setPolynomialUncert( i, rhs_uncerts[i] );
            }
          }//for( size_t i = 0; i < rhs_fit_fors.size(); ++i )
          
          break;
        }
      }//switch( m_continuum->type() )
    }
     */
  }//if( m_continuum )
  
}//void inheritUserSelectedOptions(...)


bool PeakContinuum::parametersProbablySet() const
{
  switch( m_type )
  {
    case NoOffset:
      return true;
    
    case Constant:     case Linear:
    case Quadratic:    case Cubic:
    case FlatStep:     case LinearStep:
    case BiLinearStep:
    {
      for( const auto &v : m_values )
      {
        if( v != 0.0 )
          return true;
      }
      break;
    }// polynomial or step continuum
    
    case External:
      return !!m_externalContinuum;
    break;
  }//switch( m_type )
  
  return false;
}//bool defined() const



double PeakDef::lowerX() const
{
  if( m_continuum->PeakContinuum::energyRangeDefined() )
    return m_continuum->lowerEnergy();

  const double mean = m_coefficients[PeakDef::Mean];
  const double sigma = m_coefficients[PeakDef::Sigma];
  
  switch( m_skewType )
  {
    case NoSkew:
      return mean - 4.0*sigma;
      
    case Bortel:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation
      const double skew = m_coefficients[PeakDef::SkewPar0];
      double lowx = mean - 8.0*sigma;
      while( lowx < (mean - 4.0*sigma) )
      {
        const double cumulative = PeakDists::bortel_indefinite_integral(lowx, mean, sigma, skew);
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case Bortel
      
    case CrystalBall:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation
      const double alpha = m_coefficients[PeakDef::SkewPar0];
      const double n = m_coefficients[PeakDef::SkewPar1];
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -alpha) )
      {
        const double cumulative = PeakDists::crystal_ball_tail_indefinite_t( sigma, alpha, n, (lowx - mean)/sigma );
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case CrystalBall
      
    case DoubleSidedCrystalBall:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double alpha_left = m_coefficients[PeakDef::SkewPar0];
      const double n_left = m_coefficients[PeakDef::SkewPar1];
      const double alpha_right = m_coefficients[PeakDef::SkewPar2];
      const double n_right = m_coefficients[PeakDef::SkewPar3];
      const double norm = PeakDists::DSCB_norm( alpha_left, n_left, alpha_right, n_right );
      
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -alpha_left) )
      {
        const double t = (lowx - mean) / sigma;
        const double cumulative = norm * PeakDists::DSCB_left_tail_indefinite_non_norm_t( alpha_left, n_left, t );
        
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case CrystalBall
      
    case GaussExp:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double skew = m_coefficients[PeakDef::SkewPar0];
      
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -skew) )
      {
        const double cumulative = PeakDists::gauss_exp_indefinite(mean, sigma, skew, lowx);
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case GaussExp:
      
      
    case ExpGaussExp:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double skew_left = m_coefficients[PeakDef::SkewPar0];
      const double skew_right = m_coefficients[PeakDef::SkewPar1];
      
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -skew_left) )
      {
        const double cumulative = PeakDists::exp_gauss_exp_indefinite( mean, sigma, skew_left, skew_right, lowx);
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case ExpGaussExp:
  }//switch( m_skewType )

  assert( 0 );
  return mean - 4.0*sigma;
}//double lowerX() const


double PeakDef::upperX() const
{
  if( m_continuum->PeakContinuum::energyRangeDefined() )
    return m_continuum->upperEnergy();
  
  const double mean = m_coefficients[PeakDef::Mean];
  const double sigma = m_coefficients[PeakDef::Sigma];
  
  switch( m_skewType )
  {
    case NoSkew:
    case GaussExp:
    case Bortel:
    case CrystalBall:
      return mean + 4.0*sigma;
      break;
      
    case DoubleSidedCrystalBall:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double alpha_left = m_coefficients[PeakDef::SkewPar0];
      const double n_left = m_coefficients[PeakDef::SkewPar1];
      const double alpha_right = m_coefficients[PeakDef::SkewPar2];
      const double n_right = m_coefficients[PeakDef::SkewPar3];
      
      const double oneOverSqrt2 = boost::math::constants::one_div_root_two<double>(); //0.70710678118654752440
      const double sqrtPiOver2 = boost::math::constants::root_half_pi<double>();      //1.2533141373155002512078826424
      
      
      double start_cumulative = 0.0;
      start_cumulative += PeakDists::DSCB_left_tail_indefinite_non_norm_t( alpha_left,
                                                      n_left, -alpha_left );
      start_cumulative += PeakDists::DSCB_gauss_indefinite_non_norm_t( alpha_right );
      start_cumulative -= PeakDists::DSCB_gauss_indefinite_non_norm_t( -alpha_left );
      start_cumulative -= PeakDists::DSCB_right_tail_indefinite_non_norm_t( alpha_right, n_right, alpha_right );

      const double norm = PeakDists::DSCB_norm( alpha_left, n_left, alpha_right, n_right );
      start_cumulative *= norm;
      
      double upper_t = alpha_right + 1;
      while( upper_t < 8 )
      {
        const double cumulative = start_cumulative
                + norm*PeakDists::DSCB_right_tail_indefinite_non_norm_t( alpha_right, n_right, upper_t );
        if( cumulative >= 0.975 )
          return (upper_t*sigma + mean);
        upper_t += 0.1;
      }
      return mean + 8.0*sigma;
    }//case DoubleSidedCrystalBall:
      
    case ExpGaussExp:
    {
      // Find `x` so that only 2.5% of the distribution is to the right, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double skew_left = m_coefficients[PeakDef::SkewPar0];
      const double skew_right = m_coefficients[PeakDef::SkewPar1];
      
      double upperx = mean + skew_right*sigma;
      while( upperx < (mean + 8.0*sigma) )
      {
        const double cumulative = PeakDists::exp_gauss_exp_indefinite( mean, sigma, skew_left, skew_right, upperx);
        if( cumulative >= 0.975 )
          return upperx;
        upperx += 0.1*sigma;
      }
      return mean + 8.0*sigma;
      break;
    }//case DoubleSidedCrystalBall:, case ExpGaussExp:
  }//switch( m_skewType )
  
  assert( 0 );
  return mean + 4.0*sigma;
}//double upperX() const


double PeakDef::gauss_integral( const double x0, const double x1 ) const
{
  const double &mean = m_coefficients[CoefficientType::Mean];
  const double &sigma = m_coefficients[CoefficientType::Sigma];
  const double &amp = m_coefficients[CoefficientType::GaussAmplitude];
  
  switch( m_skewType )
  {
    case SkewType::NoSkew:
      return amp*PeakDists::gaussian_integral( mean, sigma, x0, x1 );
          
    case SkewType::Bortel:
      return amp*PeakDists::bortel_integral( mean, sigma, m_coefficients[CoefficientType::SkewPar0], x0, x1 );
      
    //case SkewType::Doniach:
    //  return amp*doniach_integral( x0, x1, mean, sigma, m_coefficients[CoefficientType::SkewPar0] );
      
    case SkewType::CrystalBall:
      return amp*PeakDists::crystal_ball_integral( mean, sigma,
                            m_coefficients[CoefficientType::SkewPar0],
                            m_coefficients[CoefficientType::SkewPar1],
                            x0, x1 );
      
    case SkewType::DoubleSidedCrystalBall:
      return amp*PeakDists::double_sided_crystal_ball_integral( mean, sigma,
                                       m_coefficients[CoefficientType::SkewPar0],
                                       m_coefficients[CoefficientType::SkewPar1],
                                       m_coefficients[CoefficientType::SkewPar2],
                                       m_coefficients[CoefficientType::SkewPar3],
                                                         x0, x1 );
      break;
      
    case SkewType::GaussExp:
      return amp*PeakDists::gauss_exp_integral( mean, sigma,
                                         m_coefficients[CoefficientType::SkewPar0], x0, x1 );
      break;
      
    case SkewType::ExpGaussExp:
      return amp*PeakDists::exp_gauss_exp_integral( mean, sigma,
                                             m_coefficients[CoefficientType::SkewPar0],
                                             m_coefficients[CoefficientType::SkewPar1], x0, x1 );
      break;
  };//enum SkewType

  assert( 0 );
  throw runtime_error( "Invalid skew type" );
  return 0.0;
}//double gauss_integral( const double x0, const double x1 ) const;


void PeakDef::gauss_integral( const float *energies, double *channels, const size_t nchannel ) const 
{
  const double mean = m_coefficients[PeakDef::Mean];
  const double sigma = m_coefficients[PeakDef::Sigma];
  const double amp = m_coefficients[PeakDef::GaussAmplitude];
  
  PeakDists::photopeak_function_integral( mean, sigma, amp, m_skewType,
                                       m_coefficients + PeakDef::SkewPar0,
                                       nchannel, energies, channels );
}//void gauss_integral( const float * const energies, double *channels, const size_t nchannel )


double PeakDef::offset_integral( const double x0, const double x1,
                                 const std::shared_ptr<const SpecUtils::Measurement> &data ) const
{
  return m_continuum->offset_integral( x0, x1, data );
}//double offset_integral( const double x0, const double x1 ) const


////////////////////////////////////////////////////////////////////////////////

/** Text appropriate for use as a label for the continuum type in the gui. */
const char *PeakContinuum::offset_type_label_tr( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case PeakContinuum::NoOffset:     return "pct-none";
    case PeakContinuum::Constant:     return "pct-constant";
    case PeakContinuum::Linear:       return "pct-linear";
    case PeakContinuum::Quadratic:    return "pct-quadratic";
    case PeakContinuum::Cubic:        return "pct-cubic";
    case PeakContinuum::FlatStep:     return "pct-flat-step";
    case PeakContinuum::LinearStep:   return "pct-linear-step";
    case PeakContinuum::BiLinearStep: return "pct-bilinear-step";
    case PeakContinuum::External:     return "pct-global";
  }//switch( type )
  
  assert( 0 );
  return "pct-invalid";
}


const char *PeakContinuum::offset_type_str( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case NoOffset:     return "NoOffset";
    case Constant:     return "Constant";
    case Linear:       return "Linear";
    case Quadratic:    return "Quardratic"; //Note mispelling of "Quardratic" is left for backwards compatibility (InterSpec v1.0.6 and before), but should eventually be fixed...
    case Cubic:        return "Cubic";
    case FlatStep:     return "FlatStep";
    case LinearStep:   return "LinearStep";
    case BiLinearStep: return "BiLinearStep";
    case External:     return "External";
  }//switch( m_type )
  
  return "InvalidOffsetType";
}


size_t PeakContinuum::num_parameters( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case OffsetType::NoOffset:
    case OffsetType::External:
      return 0;
      
    case OffsetType::Constant:
    case OffsetType::Linear:
    case OffsetType::Quadratic:
    case OffsetType::Cubic:
      return static_cast<size_t>(type);
      
    case OffsetType::FlatStep:
    case OffsetType::LinearStep:
    case OffsetType::BiLinearStep:
      return 2 + (type - FlatStep);
  }//switch( type )
  
  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum polynomial type." );
  
  return 0;
}//size_t num_parameters( const OffsetType type );


bool PeakContinuum::is_step_continuum( const OffsetType type )
{
  switch( type )
  {
    case PeakContinuum::NoOffset: case PeakContinuum::External:
    case PeakContinuum::Constant: case PeakContinuum::Linear:
    case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
      return false;
      
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
      return true;
  }//switch( cont->type() )
  
  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum polynomial type." );
  return false;
}//bool is_step_continuum( const OffsetType type );


PeakContinuum::OffsetType PeakContinuum::str_to_offset_type_str( const char * const str, const size_t len )
{
  using ::rapidxml::internal::compare;
  
  /*
   // Alternative implementation for this function that is untested
  for( OffsetType type = OffsetType(0); type <= External; type = OffsetType(type+1) )
  {
    const char * const teststr = offset_type_str(type);
    const size_t teststr_len = strlen(teststr);
    if( compare(str,len,teststr,teststr_len,false) )
      return type;
  }
  if( compare(str,len,"Quadratic",9,false) )
    return PeakContinuum::Quadratic;
  throw runtime_error( "Invalid continuum type" );
  */
  
  if( compare(str,len,"NoOffset",8,false) )
    return PeakContinuum::NoOffset;
  
  if( compare(str,len,"Constant",8,false) )
    return PeakContinuum::Constant;
  
  if( compare(str,len,"Linear",6,false) )
    return PeakContinuum::Linear;
  
  if( compare(str,len,"Quardratic",10,false) || compare(str,len,"Quadratic",9,false) )
    return PeakContinuum::Quadratic;
  
  if( compare(str,len,"Cubic",5,false) )
    return PeakContinuum::Cubic;
  
  if( compare(str,len,"FlatStep",8,false) )
    return PeakContinuum::FlatStep;
  
  if( compare(str,len,"LinearStep",10,false) )
    return PeakContinuum::LinearStep;
  
  if( compare(str,len,"BiLinearStep",12,false) )
    return PeakContinuum::BiLinearStep;
  
  if( compare(str,len,"External",8,false) )
    return PeakContinuum::External;
  
  throw runtime_error( "Invalid continuum type" );
}//str_to_offset_type_str(...)



PeakContinuum::PeakContinuum()
: m_type( PeakContinuum::NoOffset ),
  m_lowerEnergy( 0.0 ),
  m_upperEnergy( 0.0 ),
  m_referenceEnergy( 0.0 )
{
}//PeakContinuum constructor

bool PeakContinuum::operator==( const PeakContinuum &rhs ) const
{
  return m_type==rhs.m_type
         && m_lowerEnergy == rhs.m_lowerEnergy
         && m_lowerEnergy == rhs.m_lowerEnergy
         && m_upperEnergy == rhs.m_upperEnergy
         && m_referenceEnergy == rhs.m_referenceEnergy
         && m_values == rhs.m_values
         && m_uncertainties == rhs.m_uncertainties
         && m_fitForValue == rhs.m_fitForValue
         && m_externalContinuum == rhs.m_externalContinuum;
}

void PeakContinuum::setParameters( double referenceEnergy,
                                   const std::vector<double> &x,
                                   const std::vector<double> &uncertainties )
{
  // First check size of inputs are valid
  const size_t num_expected_pars = PeakContinuum::num_parameters( m_type );
  
  if( x.size() != num_expected_pars )
    throw runtime_error( "PeakContinuum::setParameters invalid parameter size" );
  
  if( !uncertainties.empty() && (uncertainties.size() != num_expected_pars) )
    throw runtime_error( "PeakContinuum::setParameters invalid uncert size" );
  
  m_values = x;
  m_referenceEnergy = referenceEnergy;
  m_fitForValue.resize( m_values.size(), true );
  m_uncertainties = uncertainties;
  m_uncertainties.resize( m_values.size(), 0.0 );
}//void setParameters(...)



bool PeakContinuum::setPolynomialCoefFitFor( size_t polyCoefNum, bool fit )
{
  if( polyCoefNum >= m_fitForValue.size() )
    return false;
  
  m_fitForValue[polyCoefNum] = fit;
  return true;
}//bool setPolynomialCoefFitFor( size_t polyCoefNum, bool fit )


bool PeakContinuum::setPolynomialCoef( size_t polyCoef, double val )
{
  if( polyCoef >= m_values.size() )
    return false;
  
  m_values[polyCoef] = val;
  return true;
}


bool PeakContinuum::setPolynomialUncert( size_t polyCoef, double val )
{
  if( polyCoef >= m_uncertainties.size() )
    return false;
  
  m_uncertainties[polyCoef] = val;
  return true;
}

void PeakContinuum::setParameters( double referenceEnergy,
                                   const double *parameters,
                                   const double *uncertainties )
{
  if( !parameters )
    throw runtime_error( "PeakContinuum::setParameters invalid parameters" );
  
  vector<double> uncerts, values;
  
  switch( m_type )
  {
    case NoOffset: case External:
      throw runtime_error( "PeakContinuum::setParameters invalid m_type" );
      
    case Constant:   case Linear:
    case Quadratic: case Cubic:
    case FlatStep:
    case LinearStep:
    case BiLinearStep:
    {
      const size_t npar = num_parameters(m_type);
      values.insert( end(values), parameters, parameters + npar );
      if( uncertainties )
        uncerts.insert( end(uncerts), uncertainties, uncertainties + npar );
      
      break;
    }
  };//switch( m_type )
  
  setParameters( referenceEnergy, values, uncerts );
}//setParameters


void PeakContinuum::setExternalContinuum( const std::shared_ptr<const Measurement> &data )
{
  if( m_type != External )
    throw runtime_error( "PeakContinuum::setExternalContinuum invalid m_type" );
 
  m_externalContinuum = data;
}//setExternalContinuum(...)


void PeakContinuum::setRange( const double lower, const double upper )
{
  m_lowerEnergy = lower;
  m_upperEnergy = upper;
  if( m_lowerEnergy > m_upperEnergy )
    std::swap( m_lowerEnergy, m_upperEnergy );
}//void setRange( const double lowerenergy, const double upperenergy )


bool PeakContinuum::energyRangeDefined() const
{
  return (m_lowerEnergy != m_upperEnergy);
}//bool energyRangeDefined() const


bool PeakContinuum::isPolynomial() const
{
  switch( m_type )
  {
    case NoOffset:   case External:
      return false;
      
    case Constant:   case Linear:
    case Quadratic: case Cubic:
    case FlatStep:  case LinearStep:
    case BiLinearStep:
      return true;
      
    break;
  }//switch( m_type )
  
  return false;
}//bool isPolynomial() const


void PeakContinuum::setType( PeakContinuum::OffsetType type )
{
  const PeakContinuum::OffsetType oldType = m_type;
  
  m_type = type;
  
  switch( m_type )
  {
    case NoOffset:
      m_values.clear();
      m_uncertainties.clear();
      m_fitForValue.clear();
      m_externalContinuum.reset();
      //m_lowerEnergy = m_upperEnergy = m_referenceEnergy = 0.0;
    break;
      
    case Constant:
      m_values.resize( 1, 0.0 );
      m_uncertainties.resize( 1, 0.0 );
      m_fitForValue.resize( 1, true );
      m_externalContinuum.reset();
    break;
    
    case Linear:
      m_values.resize( 2, 0.0 );
      m_uncertainties.resize( 2, 0.0 );
      m_fitForValue.resize( 2, true );
      m_externalContinuum.reset();
      
      if( oldType == FlatStep )
        m_values[1] = m_uncertainties[1] = 0.0;
      break;
      
    case Quadratic:
      m_values.resize( 3, 0.0 );
      m_uncertainties.resize( 3, 0.0 );
      m_fitForValue.resize( 3, true );
      m_externalContinuum.reset();
      if( (oldType == BiLinearStep) || (oldType == LinearStep) )
        m_values[2] = m_uncertainties[2] = 0.0;
      break;
      
    case Cubic:
      m_values.resize( 4, 0.0 );
      m_uncertainties.resize( 4, 0.0 );
      m_fitForValue.resize( 4, true );
      m_externalContinuum.reset();
      if( oldType == BiLinearStep )
        m_values[2] = m_values[3] = m_uncertainties[2] = m_uncertainties[3] = 0.0;
      break;
      
    case FlatStep:
      if( (oldType == LinearStep) && (m_values.size() > 2) ) //Preserve the amount of "step"
      {
        m_values[1] = m_values[2];
        m_uncertainties[1] = m_uncertainties[2];
      }else
      {
        if( (m_values.size() > 1) && (m_uncertainties.size() > 1) )
          m_values[1] = m_uncertainties[1] = 0.0;
      }
      
      m_values.resize( 2, 0.0 );
      m_uncertainties.resize( 2, 0.0 );
      m_fitForValue.resize( 2, true );
      m_externalContinuum.reset();
    break;
      
    case LinearStep:
      m_values.resize( 3, 0.0 );
      m_uncertainties.resize( 3, 0.0 );
      m_fitForValue.resize( 3, true );
      m_externalContinuum.reset();
      
      if( oldType == FlatStep ) //Preserve the amount of "step"
      {
        m_values[2] = m_values[1];
        m_uncertainties[2] = m_uncertainties[1];
        m_values[1] = m_uncertainties[1] = 0.0;
      }else
      {
        //If going from BiLinear/Quadratic/Cubic to LinearStep the last parameter makes no sense to keep around
        m_values[2] = m_uncertainties[2] = 0.0;
      }
      
      break;
      
    case BiLinearStep:
      if( (oldType == FlatStep) && (m_values.size() > 1) )
        m_values[1] = m_uncertainties[1] = 0.0;
         
      m_values.resize( 4, 0.0 );
      m_uncertainties.resize( 4, 0.0 );
      m_fitForValue.resize( 4, true );
      m_externalContinuum.reset();
      //If going from Cubic to BiLinearStep the last two parameters makes no sense to keep around
      //  so just the right line equal to the left line
      m_values[2] = m_values[0];
      m_values[3] = m_values[1];
      m_uncertainties[2] = m_uncertainties[0];
      m_uncertainties[3] = m_uncertainties[1];
      break;
      
    case External:
      m_values.clear();
      m_uncertainties.clear();
      m_fitForValue.clear();
      m_referenceEnergy = 0.0;
    break;
  };//switch( m_type )
  
}//void setType( PeakContinuum::OffsetType type )


void PeakContinuum::calc_linear_continuum_eqn( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                              const double reference_energy,
                                              const double roi_start, const double roi_end,
                               const size_t num_lower_channels,
                               const size_t num_upper_channels )
{
  assert( data );
  assert( reference_energy >= roi_start );
  assert( reference_energy <= roi_end );
  assert( roi_end >= roi_start );
  assert( num_lower_channels > 0 );
  assert( num_upper_channels > 0 );
  
  if( !data || !data->energy_calibration() || !data->energy_calibration()->valid() )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: invalid input data" );
  
  if( (reference_energy < roi_start) || (reference_energy > roi_end) )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: reference energy must be within ROI" );
  
  if( roi_end < roi_start )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: lower energy greater than upper energy" );
  
  if( !num_lower_channels || !num_upper_channels )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: number of above/below channels must not be zero" );
  
  m_referenceEnergy = reference_energy;
  m_lowerEnergy = roi_start;
  m_upperEnergy = roi_end;
  m_type = PeakContinuum::Linear;
  m_values.resize( 2 );
  m_fitForValue.resize( 2, true );
  
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  assert( cal && cal->valid() );
  
  // We'll try to make the side channels be completely independent from the channels in the ROI, but
  //  we'll let there be up 10% of a channel overlap.
  double lower_cont_bound = cal->channel_for_energy( roi_start );
  double upper_cont_bound = cal->channel_for_energy( roi_end );
  
  double frac_part = lower_cont_bound - std::floor(lower_cont_bound);
  if( frac_part > 0.9 )
    lower_cont_bound = std::round( lower_cont_bound );
  else
    lower_cont_bound = std::floor( lower_cont_bound );
  
  frac_part = upper_cont_bound - std::floor(upper_cont_bound);
  if( frac_part < 0.1 )
    upper_cont_bound = std::round( upper_cont_bound );
  else
    upper_cont_bound = std::floor( upper_cont_bound );
  upper_cont_bound = std::max( upper_cont_bound, 1.0 );
  
  const size_t roi_first_channel = static_cast<size_t>( lower_cont_bound );
  // Upper cont bound give the channel number whos lower edge defines the upper bound on ROI, so we
  //  will to subtract 1 from it to find the last channel _in_ the ROI, however, we dont want to
  //  make the upper ROI channel be lower than the low one
  const size_t roi_last_channel = std::max( roi_first_channel, static_cast<size_t>( upper_cont_bound - 1.0 ) );
  
  double &m = m_values[1];
  double &b = m_values[0];
  
  PeakContinuum::eqn_from_offsets( roi_first_channel, roi_last_channel, reference_energy, data,
                                   num_lower_channels, num_upper_channels, m, b );
}//PeakContinuum::calc_linear_continuum_eqn


double PeakContinuum::offset_integral( const double x0, const double x1,
                                       const std::shared_ptr<const SpecUtils::Measurement> &data  ) const
{
  
  // A lambda to integrate data for step function continuums
  auto integrate_for_step = [this,x0,x1,data]( double &roi_lower, double &roi_upper, double &cumulative_data, double &roi_data_sum ) {
    
    if( !data || !data->num_gamma_channels() )
      throw runtime_error( "PeakContinuum::offset_integral: invalid data spectrum passed in" );
    
    roi_lower = lowerEnergy();
    roi_upper = upperEnergy();
    
    // To be consistent with how fit_amp_and_offset(...) handles things, we will do our own
    //  summing here, rather than calling Measurement::gamma_integral.
    const vector<float> &counts = *data->gamma_counts();
    const size_t mid_channel = data->find_gamma_channel( 0.5*(x0 + x1) );
    const size_t lower_channel = data->find_gamma_channel( roi_lower );
    const size_t upper_channel = data->find_gamma_channel( roi_upper );
    
    // Adjusting roi_lower and roi_upper for rounding to channel edges to match fit_amp_and_offset(...)
    roi_lower = data->gamma_channel_lower(lower_channel);
    roi_upper = data->gamma_channel_upper(upper_channel);
    
    //assert( mid_channel >= lower_channel ); //PeakDef::gaus_peaks_to_json will call a few channels below and above ROI extents
    //assert( mid_channel <= upper_channel );
    assert( lower_channel <= upper_channel );
    assert( upper_channel < counts.size() );
    
    cumulative_data = 0.0;
    roi_data_sum = 0.0;
    for( size_t i = lower_channel; i <= upper_channel; ++i )
    {
      roi_data_sum += counts[i];
      cumulative_data += (i < mid_channel) ? counts[i] : 0.0f;
    }
    
    if( (mid_channel >= lower_channel) && (mid_channel <= upper_channel) )
      cumulative_data += 0.5*counts[mid_channel];
  };//integrate_for_step lamda
  
  
  switch( m_type )
  {
    case NoOffset:
      return 0.0;
    
    case Constant: case Linear: case Quadratic: case Cubic:
      return offset_eqn_integral( &(m_values[0]), m_type, x0, x1, m_referenceEnergy );
    
    case FlatStep:
    case LinearStep:
    {
      assert( m_values.size() == (2 + (m_type - FlatStep)) );
      
      if( !data || !data->num_gamma_channels() )
        throw runtime_error( "PeakContinuum::offset_integral: invalid data spectrum passed in" );
      
      // If you change any of this, make sure you update fit_amp_and_offset(...) as well
      
      // TODO: this is not efficient, tested, or correct if integration limits do not correspond to channel limits
#ifndef _WIN32  
      #warning "TODO: Flat/Linear/BiLinear-Step offset_integral is not efficient, tested, or correct"
#endif
      double roi_lower, roi_upper, cumulative_data, roi_data_sum;
      integrate_for_step( roi_lower, roi_upper, cumulative_data, roi_data_sum );
      
      const double x0_rel = x0 - m_referenceEnergy;
      const double x1_rel = x1 - m_referenceEnergy;
      
      const double frac_data = cumulative_data / roi_data_sum;
      
      const double offset = m_values[0]*(x1_rel - x0_rel);
      const double linear = ((m_type == FlatStep) ? 0.0 :  0.5*m_values[1]*(x1_rel*x1_rel - x0_rel*x0_rel));
      const size_t step_index = ((m_type == FlatStep) ? 1 : 2);
      const double step_contribution = m_values[step_index] * frac_data * (x1_rel - x0_rel);
       
      const double answer = offset + linear + step_contribution;
      
      return std::max( answer, 0.0 );
    }//case FlatStep and LinearStep

      
    case BiLinearStep:
    {
      assert( m_values.size() == 4 );
      
      double roi_lower, roi_upper, cumulative_data, roi_data_sum;
      integrate_for_step( roi_lower, roi_upper, cumulative_data, roi_data_sum );
      
      const double x0_rel = x0 - m_referenceEnergy;
      const double x1_rel = x1 - m_referenceEnergy;
      
      const double frac_data = cumulative_data / roi_data_sum;
      
      const double left_poly = m_values[0]*(x1_rel - x0_rel) + 0.5*m_values[1]*(x1_rel*x1_rel - x0_rel*x0_rel);
      const double right_poly = m_values[2]*(x1_rel - x0_rel) + 0.5*m_values[3]*(x1_rel*x1_rel - x0_rel*x0_rel);
      
      const double contrib = ((1.0 - frac_data) * left_poly) + (frac_data * right_poly);
      
      return contrib;
    }//case BiLinearStep:
      
    case External:
      if( !m_externalContinuum )
        return 0.0;
      return gamma_integral( m_externalContinuum, x0, x1 );
    break;
  };//enum OffsetType
  
  return 0.0;
}//double offset_integral( const double x0, const double x1 ) const



void PeakContinuum::eqn_from_offsets( size_t lowchannel,
                             size_t highchannel,
                             const double reference_energy,
                             const std::shared_ptr<const SpecUtils::Measurement> &data,
                             const size_t num_lower_channels,
                             const size_t num_upper_channels,
                             double &m, double &b )
{
  // y = m*x + b, where y is density of continuum counts, per keV, at energy x
  
  
  if( !data || !data->energy_calibration() || !data->energy_calibration()->valid() )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: invalid input data" );
  
  if( lowchannel > highchannel )
    throw runtime_error( "PeakContinuum::eqn_from_offsets: lower channel greater than upper channel" );
  
  if( !num_lower_channels || !num_upper_channels )
    throw runtime_error( "PeakContinuum::eqn_from_offsets: number of above/below channels must not be zero" );
  
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  assert( cal && cal->valid() );
  const size_t nchannel = cal->num_channels();
  const size_t last_channel = nchannel ? nchannel - 1 : size_t(0);
  
  if( (num_upper_channels >= nchannel) || (num_lower_channels >= nchannel) )
    throw runtime_error( "PeakContinuum::eqn_from_offsets: number of above/below channels is to large" );
  
  lowchannel = std::max( lowchannel, num_lower_channels );
  if( (highchannel + num_upper_channels) >= nchannel )
    highchannel = nchannel - num_upper_channels - 1;
    
  
  // Now get the first and last channels for regions above/below ROI - note that these will be
  //  inclusive; e.g., if num_lower_channels==1, then first and last channel will be equal.
  //  We will also make sure to not go below zero, but if we go above last channel in spectrum,
  //  we'll just be lazy and assume they have zero values anyway.
  
  const size_t lower_cont_first_channel = lowchannel - num_lower_channels;
  const size_t lower_cont_last_channel =  lowchannel - 1;
  
  const size_t upper_cont_first_channel = highchannel + 1;
  const size_t upper_cont_last_channel = upper_cont_first_channel + num_upper_channels - 1;
  
  const double lower_count_sum = data->gamma_channels_sum( lower_cont_first_channel, lower_cont_last_channel );
  const double upper_count_sum = data->gamma_channels_sum( upper_cont_first_channel, upper_cont_last_channel );
  
  const double lower_low_energy = cal->energy_for_channel( std::min(lower_cont_first_channel, last_channel) );
  const double lower_up_energy = cal->energy_for_channel( std::min(lower_cont_last_channel + 1, last_channel) );
  
  const double upper_low_energy = cal->energy_for_channel( std::min(upper_cont_first_channel, last_channel) );
  const double upper_up_energy = cal->energy_for_channel( std::min(upper_cont_last_channel + 1, last_channel) );
  
  const double lower_dx = lower_up_energy - lower_low_energy;
  const double upper_dx = upper_up_energy - upper_low_energy;
  
  const double y1 = lower_count_sum / (1.0 + lower_cont_last_channel - lower_cont_first_channel);
  const double y2 = upper_count_sum / (1.0 + upper_cont_last_channel - upper_cont_first_channel);
  
  const double lower_x1_rel   = lower_low_energy - reference_energy;
  const double lower_x2_rel   = lower_up_energy - reference_energy;
  const double upper_x1_rel   = upper_low_energy - reference_energy;
  const double upper_x2_rel   = upper_up_energy - reference_energy;
  const double lower_sqr_diff = lower_x2_rel*lower_x2_rel - lower_x1_rel*lower_x1_rel;
  const double upper_sqr_diff = upper_x2_rel*upper_x2_rel - upper_x1_rel*upper_x1_rel;
  
  // The equation "mx + b" gives us the density of continuum counts at energy = x, so we will
  //  integrate this equation over the region below the ROI, and above the ROI, and do a little
  //  algebra so that we solve for m and b, such that we get the exact amount of counts from the
  //  equation, as is in data
  //
  //lower_count_sum = lower_dx * b + 0.5*m*lower_sqr_diff
  //upper_count_sum = upper_dx * b + 0.5*m*upper_sqr_diff
  //  ==> b = (lower_count_sum - 0.5*m*lower_sqr_diff)/lower_dx
  //  ==> upper_count_sum = upper_dx * ((lower_count_sum - 0.5*m*lower_sqr_diff)/lower_dx) + 0.5*m*upper_sqr_diff
  //      upper_count_sum = upper_dx*lower_count_sum/lower_dx - upper_dx*0.5*m*lower_sqr_diff/lower_dx + 0.5*m*upper_sqr_diff
  //      upper_count_sum = upper_dx*lower_count_sum/lower_dx + 0.5*m*(upper_sqr_diff - upper_dx*lower_sqr_diff/lower_dx)
  //      upper_count_sum - upper_dx*lower_count_sum/lower_dx = 0.5*m*(upper_sqr_diff - upper_dx*lower_sqr_diff/lower_dx)
  //  ==> m = 2*(upper_count_sum - upper_dx*lower_count_sum/lower_dx) / (upper_sqr_diff - upper_dx*lower_sqr_diff/lower_dx)
  
  if( fabs( (upper_sqr_diff - (upper_dx*lower_sqr_diff/lower_dx)) ) < FLT_EPSILON )
  {
    // Extraordinarily unlikely to happen; just assign the continuum to be flat, and average density
    //  of below and above ROI
    m = 0.0;
    b = (0.5 * lower_count_sum / lower_dx) + (0.5 * upper_count_sum / upper_dx);
  }else
  {
    const double ratio_dx = upper_dx / lower_dx;
    const double mult_m = 2.0 / (upper_sqr_diff - ratio_dx*lower_sqr_diff);
    
    // Calc m first, since b is dependent on m.
    m = mult_m * (upper_count_sum - (upper_dx*lower_count_sum/lower_dx));
    b = (lower_count_sum - 0.5*m*lower_sqr_diff) / lower_dx;
    
    // TODO: also give uncertainty of parameters, dont think its to hard, e.g. (but needs to be
    //       double checked, I did the differentiation/summing correctly in my head):
    //  const double uncert_m = mult_m * sqrt( upper_count_sum + ratio_dx*ratio_dx*lower_count_sum );
    //  const double uncert_b = (1.0 / lower_dx) * sqrt( lower_count_sum + 0.5*0.5*lower_sqr_diff*lower_sqr_diff*uncert_m*uncert_m )
  }
  
#if( PERFORM_DEVELOPER_CHECKS )
  const double coefs[2] = { b, m };
  const double eqn_lower_counts = PeakContinuum::offset_eqn_integral( coefs,
                                                                     PeakContinuum::OffsetType::Linear,
                                                                     lower_low_energy, lower_up_energy,
                                                                     reference_energy );
  
  const double eqn_upper_counts = PeakContinuum::offset_eqn_integral( coefs,
                                                                     PeakContinuum::OffsetType::Linear,
                                                                     upper_low_energy, upper_up_energy,
                                                                     reference_energy );
  
  if( fabs(eqn_lower_counts - lower_count_sum) > 1.0E-4*std::max(eqn_lower_counts,lower_count_sum)
     && (lower_count_sum > 1.0E-3)  )
  {
    char buffer[512] = { '\0' };
    snprintf( buffer, sizeof(buffer), "Region below ROI data integral not equal to equation integral."
             "\n\tReference energy: %f"
             "\n\tlower_low_energy: %f"
             "\n\tlower_up_energy: %f"
             "\n\teqn_lower_counts: %f, data lower_count_sum: %f"
             "\n\tcoefficients: %f, %f\n",
             reference_energy, lower_low_energy, lower_up_energy,
             eqn_lower_counts, lower_count_sum, b, m );
    log_developer_error( __func__, buffer );
    assert( 0 );
  }
  
  if( (fabs(eqn_upper_counts - upper_count_sum) > 1.0E-4*std::max(eqn_upper_counts,upper_count_sum))
     && (upper_count_sum > 1.0E-3) )
  {
    char buffer[512] = { '\0' };
    snprintf( buffer, sizeof(buffer), "Region Above ROI data integral not equal to equation integral."
             "\n\tReference energy: %f"
             "\n\tupper_low_energy: %f"
             "\n\tupper_up_energy: %f"
             "\n\teqn_upper_counts: %f, data upper_count_sum: %f"
             "\n\tcoefficients: %f, %f\n",
             reference_energy, upper_low_energy, upper_up_energy,
             eqn_upper_counts, upper_count_sum, b, m );
    log_developer_error( __func__, buffer );
    assert( 0 );
  }
#endif
  
  if( IsNan(m) || IsInf(m) || IsNan(b) || IsInf(b) )
  {
    // Shouldnt ever get here.
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, "PeakContinuum::eqn_from_offsets(...): Invalid results" );
#else
    cerr << "PeakContinuum::eqn_from_offsets(...): Invalid results" << endl;
#endif
    m = b = 0.0;
  }//if( an invalid value of m or b )
  
  //cout << "For reference energy " << reference_energy << " fit ceofs {" << b << ", " << m << "}" << endl;
}//void PeakContinuum::eqn_from_offsets(...)


double PeakContinuum::offset_eqn_integral( const double *coefs,
                                          PeakContinuum::OffsetType type,
                                          double x0, double x1,
                                          const double peak_mean )
{
  double answer = 0.0;
  x0 -= peak_mean;
  x1 -= peak_mean;
  
  //Explicitly evaluating the ppolynomial speeds up peak fitting by about a factor of two - suprising!
  //const int maxorder = static_cast<int>( type );
  //for( int order = 0; order < maxorder; ++order )
  //{
  //  const double exp = order + 1.0;
  //  answer += (coefs[order]/exp) * (std::pow(x1,exp) - std::pow(x0,exp));
  //}//for( int order = 0; order < maxorder; ++order )
  
  switch( type )
  {
    case NoOffset: case External: case FlatStep: case LinearStep: case BiLinearStep:
      throw runtime_error( "PeakContinuum::offset_eqn_integral(...) may only be"
                          " called for polynomial backgrounds" );
      
    case Cubic:
      answer += 0.25*coefs[3]*(x1*x1*x1*x1 - x0*x0*x0*x0);
      //fallthrough intentional
      
    case Quadratic:
      answer += 0.333333333333333*coefs[2]*(x1*x1*x1 - x0*x0*x0);
      //fallthrough intentional
      
    case Linear:
      answer += 0.5*coefs[1]*(x1*x1 - x0*x0);
      //fallthrough intentional
      
    case Constant:
      answer += coefs[0]*(x1 - x0);
      break;
  };//enum OffsetType
  
  return std::max( answer, 0.0 );
}//offset_eqn_integral(...)




void PeakContinuum::translate_offset_polynomial( double *new_coefs,
                                                 const double *old_coefs,
                                                 PeakContinuum::OffsetType type,
                                                 const double new_center,
                                                 const double old_center )
{
  switch( type )
  {
    case NoOffset:
      return;
    
    case Constant:
      new_coefs[0] = old_coefs[0];
      break;
      
    case Linear:
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      break;
      
    case Quadratic:
    case Cubic:
      throw runtime_error( "translate_offset_polynomial does not yet support "
                           "quadratic or cubic polynomials" );
    
    case LinearStep:
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      new_coefs[2] = old_coefs[2];
      break;
      
    case BiLinearStep:
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      new_coefs[2] = old_coefs[2] + old_coefs[3] * (new_center - old_center);
      new_coefs[3] = old_coefs[3];
      break;
    
    case FlatStep:
      new_coefs[0] = old_coefs[0];
      new_coefs[1] = old_coefs[1];
      break;
      
    case External:
      throw runtime_error( "translate_offset_polynomial does not support external continuum" );
  }//switch( type )
}//void translate_offset_polynomial(...)

