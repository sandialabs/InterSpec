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
#include <cmath>
#include <array>
#include <deque>
#include <tuple>
#include <chrono>
#include <thread>
#include <limits>
#include <cstdlib>
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
#include "Wt/WLocalDateTime"

#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#endif

#include "ceres/ceres.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "InterSpec/PeakDists_imp.hpp"
#include "InterSpec/RelActCalc_imp.hpp"
#include "InterSpec/RelActCalcAuto_imp.hpp"
#include "InterSpec/RelActCalc_CeresJetTraits.hpp"

using namespace std;
using namespace XmlUtils;

/** Right now if the user specifies a peak energy, we will just multiple the peak sigma (i.e.
 FWHM/2.35482) by the following value to define the ROI on either side of the mean.
 This value is arbitrary, and this is a very niave way to do things, but good enough for
 development purposes, at the moment.
 */
#define DEFAULT_PEAK_HALF_WIDTH_SIGMA 5.0


namespace
{
const double ns_decay_act_mult = SandiaDecay::MBq;
  
struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct


void sort_rois_by_energy( vector<RelActCalcAuto::RoiRange> &rois )
{
  std::sort( begin(rois), end(rois), []( const RelActCalcAuto::RoiRange &lhs,
                                        const RelActCalcAuto::RoiRange &rhs ) -> bool {
    return lhs.lower_energy < rhs.lower_energy;
  });
}//void sort_rois( vector<RoiRange> &rois )
 
  
  /*
   /// Function used to convert parameters of one FWHM type, to another.
   ///  Currently it goes from GADRAS to all other types, for a set of
   ///  nominal high and low resolution parameters.
  void fit_nominal_gadras_pars()
  {
    auto fitFromGadras = []( const vector<float> &gadras_coefs, const bool highres, const RelActCalcAuto::FwhmForm form_to_fit ){
      assert( gadras_coefs.size() == 3 );
      
      const float lower_energy = 50;
      const float upper_energy = 3000;
      const float delta_energy = 0.05*(upper_energy - lower_energy);
      
      const auto fwhm_fcn = [gadras_coefs]( float energy ) -> float {
        return DetectorPeakResponse::peakResolutionFWHM( energy,
                                                        DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn,
                                                        gadras_coefs ) / 2.35482f;
      };//fwhm_fcn lamda
      
      
      auto fake_peaks = make_shared<std::deque< std::shared_ptr<const PeakDef> > >();
      for( float ene = lower_energy; ene <=(1.001*upper_energy); ene += delta_energy )
      {
        const auto sigma = fwhm_fcn(ene);
        auto p = make_shared<PeakDef>( ene, sigma, 1000.0 );
        p->setSigmaUncert( 0.1*sigma );
        fake_peaks->push_back( p );
      }
      
      DetectorPeakResponse::ResolutionFnctForm drf_form_to_fit;
      int fit_order;
      switch( form_to_fit )
      {
        case RelActCalcAuto::FwhmForm::Gadras:
          fit_order = 3;
          drf_form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
          break;
          
        case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
          fit_order = 3;
          drf_form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse;
          break;
   
        case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
          fit_order = 2;
          drf_form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
          break;
   
        case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
          fit_order = 2;
          drf_form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
        break;
          
        case RelActCalcAuto::FwhmForm::Polynomial_2:
        case RelActCalcAuto::FwhmForm::Polynomial_3:
        case RelActCalcAuto::FwhmForm::Polynomial_4:
        case RelActCalcAuto::FwhmForm::Polynomial_5:
        case RelActCalcAuto::FwhmForm::Polynomial_6:
          drf_form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
          fit_order = num_parameters(form_to_fit);
          break;
      }//switch( rel_act_fwhm_form )
      
      vector<float> answer, uncerts;
      MakeDrfFit::performResolutionFit( fake_peaks, drf_form_to_fit, highres, fit_order, answer, uncerts );
      
      for( size_t i = 0; i < answer.size(); ++i )
        cout << "          parameters[fwhm_start + " << i << "] = " << answer[i] << endl;
    };//fitFromGadras lamda
    
    const vector<float> highres_pars{1.54f,0.264f,0.33f};
    const vector<float> lowres_pars{-6.5f,7.5f,0.55f};
    
    cout << "        if( highres )" << endl;
    cout << "        {" << endl;
    cout << "          switch( cost_functor->m_options.fwhm_form )" << endl;
    cout << "          {" << endl;
    cout << "            case RelActCalcAuto::FwhmForm::Gadras:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::Gadras );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse );
    cout << "            break;" << endl;
   
    cout << "            case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_2:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::Polynomial_2 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_3:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::Polynomial_3 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_4:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::Polynomial_4 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_5:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::Polynomial_5 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_6:" << endl;
    fitFromGadras( highres_pars, true, RelActCalcAuto::FwhmForm::Polynomial_6 );
    cout << "            break;" << endl;
    cout << "          }" << endl;
    cout << "        }else" << endl;
    cout << "        {" << endl;
    cout << "          switch( cost_functor->m_options.fwhm_form )" << endl;
    cout << "          {" << endl;
    cout << "            case RelActCalcAuto::FwhmForm::Gadras:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::Gadras );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse );
    cout << "            break;" << endl;
   
    cout << "            case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_2:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::Polynomial_2 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_3:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::Polynomial_3 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_4:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::Polynomial_4 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_5:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::Polynomial_5 );
    cout << "            break;" << endl;
    
    cout << "            case RelActCalcAuto::FwhmForm::Polynomial_6:" << endl;
    fitFromGadras( lowres_pars, false, RelActCalcAuto::FwhmForm::Polynomial_6 );
    cout << "            break;" << endl;
    cout << "          }" << endl;
    cout << "        }" << endl;
    cout << "Search for '----fit_nominal_gadras_pars----' and put these values there" << endl;
    cout << endl;
    assert( 0 );
  }//void fit_nominal_gadras_pars()
*/
  
void fill_in_default_start_fwhm_pars( std::vector<double> &parameters, size_t fwhm_start, bool highres, RelActCalcAuto::FwhmForm fwhm_form )
{
  if( highres )
  {
    // The following parameters fit from the GADRAS parameters {1.54f, 0.264f, 0.33f}, using the
    // fit_nominal_gadras_pars() commented out above.
    // ----fit_nominal_gadras_pars----
    switch( fwhm_form )
    {
      case RelActCalcAuto::FwhmForm::Gadras:
        parameters[fwhm_start + 0] = 1.54;
        parameters[fwhm_start + 1] = 0.264;
        parameters[fwhm_start + 2] = 0.33;
        break;
        
      case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
        parameters[fwhm_start + 0] = 1.86745;
        parameters[fwhm_start + 1] = 0.00216761;
        parameters[fwhm_start + 2] = 27.9835;
        break;
        
      case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
        parameters[fwhm_start + 0] = 1.0;
        parameters[fwhm_start + 1] = 0.035;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_2:
        parameters[fwhm_start + 0] = 2.10029;
        parameters[fwhm_start + 1] = 2.03657;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_3:
        parameters[fwhm_start + 0] = 2.26918;
        parameters[fwhm_start + 1] = 1.54837;
        parameters[fwhm_start + 2] = 0.192;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_4:
        parameters[fwhm_start + 0] = 2.49021;
        parameters[fwhm_start + 1] = 0.346357;
        parameters[fwhm_start + 2] = 1.3902;
        parameters[fwhm_start + 3] = -0.294974;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_5:
        parameters[fwhm_start + 0] = 2.5667;
        parameters[fwhm_start + 1] = -0.333729;
        parameters[fwhm_start + 2] = 2.59812;
        parameters[fwhm_start + 3] = -0.991013;
        parameters[fwhm_start + 4] = 0.124209;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_6:
        parameters[fwhm_start + 0] = 2.51611;
        parameters[fwhm_start + 1] = 0.318145;
        parameters[fwhm_start + 2] = 0.838887;
        parameters[fwhm_start + 3] = 0.734524;
        parameters[fwhm_start + 4] = -0.568632;
        parameters[fwhm_start + 5] = 0.0969217;
        break;
    }//switch( options.fwhm_form )
  }else
  {
    // The following parameters fit from the GADRAS parameters {-6.5f, 7.5f, 0.55f}, using the
    // fit_nominal_gadras_pars() commented out above.
    switch( fwhm_form )
    {
      case RelActCalcAuto::FwhmForm::Gadras:
        parameters[fwhm_start + 0] = -6.5;
        parameters[fwhm_start + 1] = 7.5;
        parameters[fwhm_start + 2] = 0.55;
        break;
        
      case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
        parameters[fwhm_start + 0] = -592.865;
        parameters[fwhm_start + 1] = 4.44776;
        parameters[fwhm_start + 2] = 21173.6;
        break;
      
      case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
        parameters[fwhm_start + 0] = -7.0;
        parameters[fwhm_start + 1] = 2.0;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_2:
        parameters[fwhm_start + 0] = -146.632;
        parameters[fwhm_start + 1] = 3928.7;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_3:
        parameters[fwhm_start + 0] = -101.518;
        parameters[fwhm_start + 1] = 3037.91;
        parameters[fwhm_start + 2] = 555.973;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_4:
        parameters[fwhm_start + 0] = -68.9708;
        parameters[fwhm_start + 1] = 2334.29;
        parameters[fwhm_start + 2] = 1873.59;
        parameters[fwhm_start + 3] = -428.651;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_5:
        parameters[fwhm_start + 0] = -37.7014;
        parameters[fwhm_start + 1] = 1605.68;
        parameters[fwhm_start + 2] = 4201.01;
        parameters[fwhm_start + 3] = -2230.26;
        parameters[fwhm_start + 4] = 383.314;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_6:
        parameters[fwhm_start + 0] = -11.4349;
        parameters[fwhm_start + 1] = 947.497;
        parameters[fwhm_start + 2] = 7088.34;
        parameters[fwhm_start + 3] = -5991.02;
        parameters[fwhm_start + 4] = 2192.23;
        parameters[fwhm_start + 5] = -286.53;
        break;
    }//switch( options.fwhm_form )
  }//if( highres ) / else
}//fill_in_default_start_fwhm_pars()


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
        // We need to subtract `diff` channels to our channel range; diff will usually be 1, but maybe it could be 2
        //
        //  Example case
        //    first_channel_f = 4910.48
        //    last_channel_f  = 5366.53
        //    wanted_nchan    = 457
        //    -->
        //      first_channel_i=4910, first_fractional=0.48, (i.e., we'll be off by 1 -  0.48 =0.52 if we add 1 to first_channel_i)
        //      last_channel_i=5367,  last_fractional=-0.47  (i.e., we'll be off by 1 + -0.47 =1.53 if we subtract 1 from last_channel_f)
        //      In this case niave_nchan=458, so we want to subtract 1 channel from the range
        //
        //    Second case
        //    first_channel_f = 4910.45
        //    last_channel_f  = 5366.48
        //    wanted_nchan    = 457
        //    -->
        //      first_channel_i=4910, first_fractional=0.45,  (i.e., we'll be off by 1 - 0.45 =0.55 if we add 1 from first_channel_i)
        //      last_channel_i=5367, last_fractional=-0.52    (i.e., we'll be off by 1 + -0.52=0.48 if we subtract 1 from last_channel_f)
        //      In this case niave_nchan=458, so we want to subtract 1 channel from the range
        
        if( (diff - first_fractional) < (diff + last_fractional) )
        {
          first_channel += diff;
        }else if( last_channel >= diff )
        {
          last_channel -= diff;
        }else
        {
          first_channel = 0;
          last_channel = wanted_nchan - 1;
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
  /** The activity multiple to use for this nuclide.
  
  These are the initially estimated activities of the nuclide, before the proper fit, and
  are used to keep the Ceres paramters in reasonable ranges.
  */
  double activity_multiple = 1.0;

  /** The age multiple to use for this nuclide.
  
  This is the initial age estimate, before the proper fit, and is used to keep the Ceres paramters
  in reasonable ranges.
  
  Will be set to -1.0 if age is controlled by another nuclide.
  */
  double age_multiple = 1.0;

  struct EnergyYield
  {
    double energy = 0.0;
    double yield = 0.0;
      
    size_t transition_index = 0;
    const SandiaDecay::Transition *transition = nullptr;
    PeakDef::SourceGammaType gamma_type;
  };
  
  vector<EnergyYield> nominal_gammas;
  
  /** Returns the decay gammas along with their rate, and the transition that led to them (i.e.,
   where in the decay chain they come from).
   
   Note that the returned result may have multiple entries with the exact same energy (i.e., from
   different decays).
   */
  static vector<EnergyYield> decay_gammas( const SandiaDecay::Nuclide * const parent,
                                          const double age,
                                          const std::vector<double> &gammas_to_exclude )
  {
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( parent, ns_decay_act_mult, age );
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mix.activity( 0.0 );
    
    // all_gamma_transitions may contain duplicate energies - we will combine these below
    vector<EnergyYield> all_gamma_transitions;
    all_gamma_transitions.reserve( 1536 ); //avoid a few resizes; for Pu241, at nominal age, the actual number of entries is 1229
    
    EnergyYield annihilationInfo;
    annihilationInfo.energy = 510.998910;
    annihilationInfo.yield = 0.0;
    annihilationInfo.transition = nullptr;
    annihilationInfo.transition_index = 0;
    annihilationInfo.gamma_type = PeakDef::SourceGammaType::AnnihilationGamma;
    size_t num_annih_trans = 0;
    
    for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
    {
      const SandiaDecay::Nuclide * const nuclide = activities[nucIndex].nuclide;
      const double activity = activities[nucIndex].activity / ns_decay_act_mult;
      
      if( activity <= 0.0 )
        continue;
      
      const size_t n_decaysToChildren = nuclide->decaysToChildren.size();
      
      for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
      {
        const SandiaDecay::Transition * const transition = nuclide->decaysToChildren[decayIndex];
        const size_t n_products = transition->products.size();
        
        for( size_t productNum = 0; productNum < n_products; ++productNum )
        {
          const SandiaDecay::RadParticle &particle = transition->products[productNum];
          if( particle.type == SandiaDecay::ProductType::GammaParticle )
          {
            bool exclude = false;
            for( const double exclude_energy : gammas_to_exclude )
              exclude = (exclude || (fabs(particle.energy - exclude_energy) < 0.001));
            
            if( exclude )
              continue;
            
            all_gamma_transitions.push_back( {} );
            EnergyYield &info = all_gamma_transitions.back();
            info.energy = particle.energy;
            info.yield = activity * particle.intensity * transition->branchRatio;
            
            info.transition_index = productNum;
            info.transition = transition;
            info.gamma_type = PeakDef::SourceGammaType::NormalGamma;
          }else if( particle.type == SandiaDecay::ProductType::PositronParticle )
          {
            annihilationInfo.yield += 2.0 * activity * particle.intensity * transition->branchRatio;
           
            num_annih_trans += 1;
            annihilationInfo.transition_index = productNum;
            annihilationInfo.transition = transition;
          }//if( particle.type is gamma ) / else if( position )
        }//for( size_t productNum = 0; productNum < n_products; ++productNum )
      }//for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
    }//for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
    
    
    if( annihilationInfo.yield > 0.0 )
    {
      bool exclude = false;
      for( const double exclude_energy : gammas_to_exclude )
        exclude = (exclude || (fabs(annihilationInfo.energy - exclude_energy) < 0.001));
      
      if( num_annih_trans != 1 )
      {
        annihilationInfo.transition = nullptr;
        annihilationInfo.transition_index = 0;
      }
      
      if( !exclude )
        all_gamma_transitions.push_back( annihilationInfo );
    }//if( annihilationInfo.yield > 0.0 )
    
    
    std::sort( begin(all_gamma_transitions), end(all_gamma_transitions),
               []( const auto &lhs, const auto &rhs ) { return lhs.energy < rhs.energy; });
    
    if( all_gamma_transitions.empty() )
      return all_gamma_transitions;
    
    // Now go throw and combine duplicate energies
    vector<EnergyYield> answer;
    answer.reserve( all_gamma_transitions.size() );
    answer.push_back( all_gamma_transitions.front() );
    size_t energy_start_index = 0;
    for( size_t current_index = 1; current_index < all_gamma_transitions.size(); ++current_index )
    {
      const EnergyYield &current = all_gamma_transitions[current_index];
      EnergyYield &prev = answer.back();
      if( current.energy == prev.energy )
      {
        prev.yield += current.yield;
        
        // We will assign the transition to the largest yield of this energy (not that it matters
        //  that much, I guess, given peaks only get a single source transition, but it seems like
        //  the right thing to do anyway).
        double max_yield = -1;
        for( size_t i = energy_start_index; i <= current_index; ++i )
        {
          if( all_gamma_transitions[i].yield > max_yield )
          {
            max_yield = all_gamma_transitions[i].yield;
            prev.gamma_type = all_gamma_transitions[i].gamma_type;
            prev.transition = all_gamma_transitions[i].transition;
            prev.transition_index = all_gamma_transitions[i].transition_index;
          }//
        }//for( size_t i = energy_start_index; i < current_index; ++i )
      }else
      {
        energy_start_index = current_index;
        answer.push_back( current );
      }
    }//for( size_t current_index = 1; current_index < all_gamma_transitions.size(); ++current_index )
    
    return answer;
  }//static vector<EnergyYield> decay_gammas( const SandiaDecay::Nuclide * const parent, ... )

  
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
    
    nominal_gammas = decay_gammas( nuclide, nominal_age, gammas_to_exclude );
  }//NucInputGamma constructor
};//struct NucInputGamma

  /** Setup the physical model shield parameters for the Ceres problem - used by both auto and manual.
   * Assumes each parameter is its own parameter block.
   *
   * @param lower_bounds The lower bounds the parameters can take - its optional to fill these out
   * @param upper_bounds The upper bounds the parameters can take - its optional to fill these out
   * @param constant_parameters The parameters that should be held constant in the problem
   * @param pars The array of parameters.
   * @param start_ind The starting index, of `pars`, for this shielding.  This is the index
   *                  Atomic Number of the shielding (if material isnt defined), and the next position
   *                  in the parameters is the Areal Density.
   * @param opt The physical model shield input to use.
   */
void setup_physical_model_shield_par( vector<optional<double>> &lower_bounds,
                                     vector<optional<double>> &upper_bounds,
                                     vector<int> &constant_parameters,
                                        double * const pars,
                                        const size_t start_ind,
                                        const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt )
{
  const size_t an_index = start_ind + 0;
  const size_t ad_index = start_ind + 1;
  
  if( !opt || (!opt->fit_atomic_number && !opt->material && ((opt->atomic_number < 1.0) || (opt->atomic_number > 98.0))) )
  {
    pars[an_index] = 0.0;
    pars[ad_index] = 0.0;
    constant_parameters.push_back( static_cast<int>(an_index) );
    constant_parameters.push_back( static_cast<int>(ad_index) );
    return;
  }//if( !opt )
  
  opt->check_valid();
  
  if( opt->material )
  {
    if( opt->fit_atomic_number )
      throw runtime_error( "You can not fit AN when defining a material" );
    
    pars[an_index] = 0.0;
    constant_parameters.push_back( static_cast<int>(an_index) );
  }else
  {
    double lower_an = opt->lower_fit_atomic_number;
    double upper_an = opt->upper_fit_atomic_number;
    if( (lower_an == upper_an) && (lower_an == 0.0) )
    {
      lower_an = 1.0;
      upper_an = 98.0;
    }
    
    double an = opt->atomic_number;
    if( (an == 0.0) && opt->fit_atomic_number )
      an = 0.5*(opt->lower_fit_atomic_number + opt->upper_fit_atomic_number);
    
    if( (an < 1) || (an > 98) || (an < lower_an) || (an > upper_an) )
      throw runtime_error( "Self-atten AN is invalid" );
    
    pars[an_index] = an / RelActCalc::ns_an_ceres_mult;
    
    if( opt->fit_atomic_number )
    {
      if( (lower_an < 1) || (upper_an > 98) || (upper_an <= lower_an) )
        throw runtime_error( "Self-atten AN limits is invalid" );
      
      lower_bounds[an_index] = lower_an / RelActCalc::ns_an_ceres_mult;
      upper_bounds[an_index] = upper_an / RelActCalc::ns_an_ceres_mult;
          
      // We could Quantize AN... but I'm not sure the below does it
      //std::vector<int> quantized_params(1, 0);  // 0 means continuous, 1 means quantized
      //quantized_params[an_index] = 1;  // Assuming an_index is the index of the AN parameter
      //problem->SetParameterization( pars + an_index,
      //     new ceres::SubsetParameterization(1, quantized_params));
    }else
    {
      constant_parameters.push_back( static_cast<int>(an_index) );
    }
  }//if( opt.material ) / else
  
  const double max_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
  double ad = opt->areal_density / PhysicalUnits::g_per_cm2;
  double lower_ad = opt->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
  double upper_ad = opt->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
  
  if( (lower_ad == upper_ad) && (lower_ad == 0.0) )
  {
    lower_ad = 0.0;
    upper_ad = max_ad;
  }
  
  if( (ad == 0.0) && opt->fit_areal_density )
  {
    ad = 2.5; // We want something away from zero, because Ceres doesnt like zero values much - 2.5 is arbitrary
    //  ad = 0.5*(lower_ad + upper_ad); //Something like 250 would be way too much
  }
  if( (ad < 0.0) || (ad > max_ad) )
    throw runtime_error( "Self-atten AD is invalid" );
  
  pars[ad_index] = ad;
  
  if( opt->fit_areal_density )
  {
    // Check for limits
    if( (lower_ad < 0.0) || (upper_ad > max_ad) || (lower_ad >= upper_ad) )
      throw runtime_error( "Self-atten AD limits is invalid" );
    
    lower_bounds[ad_index] = lower_ad;
    upper_bounds[ad_index] = upper_ad;
  }else
  {
    constant_parameters.push_back( static_cast<int>(ad_index) );
  }
}//void setup_physical_model_shield_par( ceres::Problem... )


struct RelActAutoCostFcn /* : ROOT::Minuit2::FCNBase() */
{
  RelActCalcAuto::Options m_options;
  std::vector<std::vector<NucInputGamma>> m_nuclides; //has same number of entries as `m_options.rel_eff_curves`
  std::vector<RoiRangeChannels> m_energy_ranges;
  
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
  
  /** Indexes of where information begin in the parameters vector - currently used as confirmation of index, but
   it probably makes sense to start using at some point.
   */
  size_t m_energy_cal_par_start_index;
  size_t m_fwhm_par_start_index;
  size_t m_rel_eff_par_start_index;
  size_t m_acts_par_start_index;
  size_t m_free_peak_par_start_index;
  size_t m_skew_par_start_index;
  size_t m_add_br_uncert_start_index; //Will only be valid if `m_peak_ranges_with_uncert` is not empty
  
  bool m_skew_has_energy_dependance;
  
  std::vector<std::pair<double,double>> m_peak_ranges_with_uncert;
  static constexpr double sm_peak_range_uncert_par_offset = 5.0;
  
  std::shared_ptr<std::atomic_bool> m_cancel_calc;
  
  /** just for debug purposes, we'll keep track of how many times the eval function gets called. */
  mutable std::atomic<size_t> m_ncalls;
  
  /** For tracking where we spend time. */
  mutable std::atomic<size_t> m_nanoseconds_spent_in_eval;
  
  // If class to implement cancelling a calculation
  class CheckCeresTerminateCallback : public ceres::IterationCallback
  {
    shared_ptr<atomic_bool> m_cancel;
    
  public:
    CheckCeresTerminateCallback( shared_ptr<atomic_bool> &cancel_calc ) : m_cancel( cancel_calc ) {}
    
    virtual ceres::CallbackReturnType operator()(const ceres::IterationSummary & )
    {
      return (m_cancel && m_cancel->load()) ? ceres::CallbackReturnType::SOLVER_ABORT
                                            : ceres::CallbackReturnType::SOLVER_CONTINUE;
    }
  };//class CheckCeresTerminateCallback
  
  
  RelActAutoCostFcn( RelActCalcAuto::Options options,
                     shared_ptr<const SpecUtils::Measurement> spectrum,
                     const vector<float> &channel_counts_uncert,
                     std::shared_ptr<const DetectorPeakResponse> drf,
                     std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                     std::shared_ptr<std::atomic_bool> cancel_calc
                           )
  : m_options( options ),
  m_nuclides{},
  m_energy_ranges{},
  m_spectrum( spectrum ),
  m_live_time( spectrum ? spectrum->live_time() : 0.0f ),
  m_channel_counts( (spectrum && spectrum->gamma_counts()) ? *spectrum->gamma_counts() : vector<float>() ),
  m_channel_count_uncerts( channel_counts_uncert ),
  m_energy_cal( spectrum ? spectrum->energy_calibration() : nullptr ),
  m_rel_eff_anchor_enhancement( 1000.0 ),
  m_drf( nullptr ),
  m_energy_cal_par_start_index( std::numeric_limits<size_t>::max() ),
  m_fwhm_par_start_index( std::numeric_limits<size_t>::max() ),
  m_rel_eff_par_start_index( std::numeric_limits<size_t>::max() ),
  m_acts_par_start_index( std::numeric_limits<size_t>::max() ),
  m_free_peak_par_start_index( std::numeric_limits<size_t>::max() ),
  m_skew_par_start_index( std::numeric_limits<size_t>::max() ),
  m_add_br_uncert_start_index( std::numeric_limits<size_t>::max() ),
  m_skew_has_energy_dependance( false ),
  m_cancel_calc( cancel_calc ),
  m_ncalls( 0 ),
  m_nanoseconds_spent_in_eval( size_t(0) )
  {
    if( !spectrum || (spectrum->num_gamma_channels() < 128) )
      throw runtime_error( "RelActAutoCostFcn: invalid spectrum." );
    
    if( !m_energy_cal || !m_energy_cal->valid() )
      throw runtime_error( "RelActAutoCostFcn: invalid energy calibration." );
    
    if( spectrum->num_gamma_channels() != channel_counts_uncert.size() )
      throw runtime_error( "RelActAutoCostFcn: Diff number of spectrum channels and channel counts uncert." );
    
    const bool highres = PeakFitUtils::is_high_res( spectrum );
    vector<RelActCalcAuto::RoiRange> energy_ranges = options.rois;
    vector<RelActCalcAuto::FloatingPeak> extra_peaks = options.floating_peaks;
    
    if( energy_ranges.empty() )
      throw runtime_error( "RelActAutoCostFcn: no energy ranges specified." );
    
    m_nuclides.resize( options.rel_eff_curves.size() );
    
    for( size_t rel_eff_index = 0; rel_eff_index < options.rel_eff_curves.size(); ++rel_eff_index )
    {
      const auto &rel_eff_input = options.rel_eff_curves[rel_eff_index];
      const vector<RelActCalcAuto::NucInputInfo> &input_nuclides = rel_eff_input.nuclides;
      
      if( input_nuclides.empty() )
        throw runtime_error( "RelActAutoCostFcn: no nuclides specified." );
    
      vector<NucInputGamma> &these_nuclides = m_nuclides[rel_eff_index];

      for( const auto &n : input_nuclides )
      {
        if( !n.nuclide )
          throw runtime_error( "RelActAutoCostFcn: null Nuclide." );
      
        for( const auto &pn : these_nuclides )
        {
          if( n.nuclide == pn.nuclide )
            throw runtime_error( "RelActAutoCostFcn: duplicate nuclide (" + n.nuclide->symbol + ")." );
        }
      
        these_nuclides.emplace_back( n );
      }//for( const auto &n : input_nuclides )

      if( rel_eff_input.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        if( rel_eff_input.phys_model_self_atten )
          rel_eff_input.phys_model_self_atten->check_valid();
        for( const auto &shield : rel_eff_input.phys_model_external_atten )
          shield->check_valid();
      }else
      {
        if( rel_eff_input.phys_model_self_atten
            || !rel_eff_input.phys_model_external_atten.empty() )
          throw runtime_error( "RelActAutoCostFcn: self attenuating, or external shielding can not"
                            " be defined unless using RelEffEqnForm::FramPhysicalModel." );
      }//if( RelEffEqnForm::FramPhysicalModel ) / else

      if( !drf && (rel_eff_input.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) )
      {
        const string drfpaths = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GenericGadrasDetectors" );
        const string drf_dir = SpecUtils::append_path( drfpaths, highres ? "HPGe 40%" : "LaBr 10%" );
        auto new_drf = make_shared<DetectorPeakResponse>();
        try
        {
          new_drf->fromGadrasDirectory( drf_dir );
          new_drf->setFwhmCoefficients( {}, DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm );
        }catch( std::exception &e )
        {
          throw runtime_error( "RelActAutoCostFcn: failed to open default DRF." );
        }
      
        drf = new_drf;
      }//if( no DRF, and were doing physical model )
    }//for( const auto &rel_eff_input : options.rel_eff_curves )
    
    
    if( cancel_calc && cancel_calc->load() )
      throw runtime_error( "User cancelled calculation." );
    
    
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
                                               spectrum->gamma_energy_max(),
                                               DetectorPeakResponse::EffGeometryType::FarField );
        
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
    for( size_t roi_index = 0; roi_index < m_options.rois.size(); ++roi_index )
    {
      const RelActCalcAuto::RoiRange &roi_range = m_options.rois[roi_index];
      
      if( (roi_range.lower_energy >= roi_range.upper_energy)
         || (roi_range.lower_energy < 0.0) )
        throw runtime_error( "RelActAutoCostFcn: mal-formed energy range" );
      
      // We'll check the ranges dont overlap (but they can touch)
      for( size_t other_roi = roi_index + 1; other_roi < m_options.rois.size(); ++other_roi )
      {
        const RelActCalcAuto::RoiRange &other_roi_range = m_options.rois[other_roi];
        
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
      
      if( cancel_calc && cancel_calc->load() )
        throw runtime_error( "User cancelled calculation." );
      
      //We'll try to limit/break-up energy ranges
      const double min_br = numeric_limits<double>::min();  //arbitrary
      const double num_sigma_half_roi = DEFAULT_PEAK_HALF_WIDTH_SIGMA;
      
      vector<pair<double,double>> gammas_in_range;
      
      // Define a helper function to add a gamma at an energy into \c gammas_in_range
      auto add_peak_to_range = [&gammas_in_range, this, &spectrum, highres, &roi_range, num_sigma_half_roi]( const double energy ){
        double energy_sigma;
        float min_sigma, max_sigma;
        expected_peak_width_limits( energy, highres, spectrum, min_sigma, max_sigma );
        
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
        
        double gamma_row_lower = energy - num_sigma_half_roi*energy_sigma;
        double gamma_row_upper = energy + num_sigma_half_roi*energy_sigma;
        
        if( !roi_range.allow_expand_for_peak_width )
        {
          gamma_row_lower = std::max( gamma_row_lower, roi_range.lower_energy );
          gamma_row_upper = std::min( gamma_row_upper, roi_range.upper_energy );
        }
        
        gammas_in_range.push_back( {gamma_row_lower,gamma_row_upper} );
      };//add_peak_to_range(...) lamda
      
      for( const vector<NucInputGamma> &nucs : m_nuclides )
      {
        for( const auto &n : nucs )
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
        }//for( const auto &n : nucs )
      }//for( const vector<NucInputGamma> &nucs : m_nuclides )
        
      for( const auto &peak : m_options.floating_peaks )
      {
        if( (peak.energy >= roi_range.lower_energy) && (peak.energy <= roi_range.upper_energy) )
          add_peak_to_range( peak.energy );
      }//for( const auto &peak : m_options.floating_peaks )
      
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
        
        cout << "Adding energy range [" << this_range.lower_energy << ", " << this_range.upper_energy << "]\n";
        
        m_energy_ranges.emplace_back( m_energy_cal, this_range );
      }//for( loop over gammas_in_range )
    }//for( const RelActCalcAuto::RoiRange &input : energy_ranges )
    
    
    std::sort( begin(m_energy_ranges), end(m_energy_ranges),
              []( const RoiRangeChannels &lhs, const RoiRangeChannels &rhs ) -> bool {
      return lhs.lower_energy < rhs.lower_energy;
    } );
    
    
    if( m_energy_ranges.empty() )
      throw runtime_error( "RelActAutoCostFcn: no gammas in the defined energy ranges." );
    
    if( cancel_calc && cancel_calc->load() )
      throw runtime_error( "User cancelled calculation." );
    
    // Check to make sure the floating peak is in a ROI range
    for( const RelActCalcAuto::FloatingPeak &peak : m_options.floating_peaks )
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
    }//for( const auto &peak : m_options.floating_peaks )
    
    
    // TODO: Figure out a proper value to set m_rel_eff_anchor_enhancement to, like maybe largest peak area divided m_live_time, or something like that
    m_rel_eff_anchor_enhancement = 0.0;
    for( const auto &r : m_energy_ranges )
      m_rel_eff_anchor_enhancement += spectrum->gamma_integral( r.lower_energy, r.upper_energy );
    
    
    // Start assigning where we expect to see each 'type' of parameter in Ceres par vectr
    m_energy_cal_par_start_index = 0;
    m_fwhm_par_start_index = RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars;
    
    const size_t num_fwhm_pars = num_parameters( options.fwhm_form );
    m_rel_eff_par_start_index = m_fwhm_par_start_index + num_fwhm_pars;
    
    size_t num_rel_eff_par = 0;
    for( const auto &rel_eff_curve : options.rel_eff_curves )
    {
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        num_rel_eff_par += 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2;
      else
        num_rel_eff_par += rel_eff_curve.rel_eff_eqn_order + 1;
    }//for( const auto &rel_eff_curve : options.rel_eff_curves )
    m_acts_par_start_index = m_rel_eff_par_start_index + num_rel_eff_par;
    
    size_t num_acts_par = 0;
    for( const auto &rel_eff_curve : options.rel_eff_curves )
      num_acts_par += 2*rel_eff_curve.nuclides.size();
    m_free_peak_par_start_index = m_acts_par_start_index + num_acts_par;
    
    const size_t num_free_peak_par = 2*extra_peaks.size();
    m_skew_par_start_index = m_free_peak_par_start_index + num_free_peak_par;
  }//RelActAutoCostFcn constructor.
  
  
  
  size_t number_parameters() const
  {
    size_t num_pars = 0;
    assert( (m_energy_cal_par_start_index == std::numeric_limits<size_t>::max())
           || m_energy_cal_par_start_index == num_pars );
    
    // Energy calibration; we will always have these, even if fixed values
    num_pars += RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars;
    
    assert( (m_fwhm_par_start_index == std::numeric_limits<size_t>::max())
           || m_fwhm_par_start_index == num_pars );
    
    // The FWHM equation
    num_pars += num_parameters( m_options.fwhm_form );
    
    assert( (m_rel_eff_par_start_index == std::numeric_limits<size_t>::max())
           || m_rel_eff_par_start_index == num_pars );
    
    // The Relative Eff coefficients
    num_pars += num_total_rel_eff_coefs();
    
    assert( (m_acts_par_start_index == std::numeric_limits<size_t>::max())
           || m_acts_par_start_index == num_pars );
    
    // The Activities; one parameter for activity, one for age (which may be fixed)
    for( const auto &n : m_nuclides )
      num_pars += 2*n.size();
    
    assert( (m_free_peak_par_start_index == std::numeric_limits<size_t>::max())
           || m_free_peak_par_start_index == num_pars );
    
    // Floating peaks; one parameter for amplitude, one for FWHM (which will usually be unused)
    num_pars += 2*m_options.floating_peaks.size();
    
    assert( (m_skew_par_start_index == std::numeric_limits<size_t>::max())
           || m_skew_par_start_index == num_pars );
    
    // Peak skew parameters; two sets of these, with some coefficients in the upper set
    //  maybe not being used
    const size_t num_skew = PeakDef::num_skew_parameters( m_options.skew_type );
    num_pars += 2*num_skew;
    
    assert( (m_add_br_uncert_start_index == numeric_limits<size_t>::max()) == m_peak_ranges_with_uncert.empty() );
    assert( (m_add_br_uncert_start_index == std::numeric_limits<size_t>::max())
           || (num_pars == m_add_br_uncert_start_index) );
    
    num_pars += m_peak_ranges_with_uncert.size();
    
    // Anything else?
    
    return num_pars;
  }//number_parameters()
  
  size_t number_residuals() const
  {
    // Number of gamma channels, plus one to anchor relative eff (if not Physical Rel Eff)
    size_t num_resids = 0;
    
    for( const auto &rel_eff : m_options.rel_eff_curves )
      num_resids += (rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel);
    
    for( const auto &r : m_energy_ranges )
      num_resids += r.num_channels;
    
    num_resids += m_peak_ranges_with_uncert.size();
    
    return num_resids;
  }//size_t number_residuals() const
  
  /** Solve the problem, using the Ceres optimizer. */
  static RelActCalcAuto::RelActAutoSolution solve_ceres( RelActCalcAuto::Options options,
                                                        std::shared_ptr<const SpecUtils::Measurement> foreground,
                                                        std::shared_ptr<const SpecUtils::Measurement> background,
                                                        const std::shared_ptr<const DetectorPeakResponse> input_drf,
                                                        std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                                                        std::shared_ptr<std::atomic_bool> cancel_calc
                                                        )
  {
    const vector<RelActCalcAuto::RoiRange> &energy_ranges = options.rois;
    const vector<RelActCalcAuto::FloatingPeak> &extra_peaks = options.floating_peaks;

    const auto start_time = std::chrono::high_resolution_clock::now();
    
    const bool highres = PeakFitUtils::is_high_res( foreground );
    
    RelActCalcAuto::RelActAutoSolution solution;
    
    DoWorkOnDestruct setFinalTime( [&solution,start_time](){
      const auto end_time = std::chrono::high_resolution_clock::now();
      solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
    });
    
    // We will make a (potentially background subtracted) spectrum - we also need to track the
    //  uncertainties in each channel
    shared_ptr<SpecUtils::Measurement> spectrum;
    vector<float> channel_counts, channel_count_uncerts;

    try
    {
      solution.m_foreground       = foreground;
      solution.m_background       = background;
      solution.m_options          = options;
      assert( !options.rel_eff_curves.empty() );
      if( options.rel_eff_curves.empty() )
        throw runtime_error( "Need at least ine relative efficiency curve" );
      
      solution.m_fwhm_form = options.fwhm_form;
      for( const auto &rel_eff_curve : options.rel_eff_curves )
        solution.m_rel_eff_forms.push_back( rel_eff_curve.rel_eff_eqn_type );
    
      if( input_drf && input_drf->isValid() && input_drf->hasResolutionInfo() )
        solution.m_drf = input_drf;
    
      const auto check_rel_eff_form = [&]( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve ){
        switch( rel_eff_curve.rel_eff_eqn_type )
        {
          case RelActCalc::RelEffEqnForm::LnX:
          case RelActCalc::RelEffEqnForm::LnY:
          case RelActCalc::RelEffEqnForm::LnXLnY:
          case RelActCalc::RelEffEqnForm::FramEmpirical:
          {
            if( rel_eff_curve.phys_model_self_atten || !rel_eff_curve.phys_model_external_atten.empty() )
            {
              throw runtime_error( "solve_ceres: You cannot provide self attenuating, or external"
                                " attenuating options unless using RelEffEqnForm::FramPhysicalModel." );
            }
        
            //if( rel_eff_curve.rel_eff_eqn_order == 0 )
            //  throw runtime_error( "Relative efficiency order must be at least 1 (but should probably be at least 2)" );
            break;
          }
        
          case RelActCalc::RelEffEqnForm::FramPhysicalModel:
          {
            if( rel_eff_curve.rel_eff_eqn_order != 0 )
              throw runtime_error( "solve_ceres: relative eff. eqn order must be zero for"
                                   " RelEffEqnForm::FramPhysicalModel." );
            break;
          }
        }//switch( solution.m_rel_eff_form )
      };//check_rel_eff_form
    
      for( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve : options.rel_eff_curves )
        check_rel_eff_form( rel_eff_curve );
    
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
    
      // Manual "relative_activity" assumes a measurement of 1-second (or rather peaks are in CPS),
      //  but this "auto" relative activity takes into account live_time
      const double live_time = spectrum->live_time();
    
      if( (!cancel_calc || !cancel_calc->load()) && all_peaks.empty() )
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
    
    
      // Check that the specified Pu242 by correlation method is valid.
      auto check_pu_corr_method = [&]( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve ){
        const vector<RelActCalcAuto::NucInputInfo> &nuclides = rel_eff_curve.nuclides;
      
        // Make a lamda that returns the mass-number of Plutonium nuclides; or if Am241 is present,
        //  will insert 241
        auto pu_iso_present = [&]() -> set<short> {
          set<short> answer;
          for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
          {
            if( nuc.nuclide
               && ((nuc.nuclide->atomicNumber == 94)
                   || ((nuc.nuclide->atomicNumber == 95) && (nuc.nuclide->massNumber == 241) ) ) )
            {
              answer.insert(nuc.nuclide->massNumber);
            }
          }//for( loop over input nucldies )
        
          return answer;
        };//pu_iso_present
      
        switch( rel_eff_curve.pu242_correlation_method )
        {
          case RelActCalc::PuCorrMethod::Bignan95_PWR:
          case RelActCalc::PuCorrMethod::Bignan95_BWR:
          {
            // Need Pu238, Pu239, and Pu240
            const set<short> punucs = pu_iso_present();
            
            if( punucs.count(242) )
              throw runtime_error( "You can not specify a Pu242 correction method AND include Pu242"
                                  " in the nuclides you are fitting in the spectrum." );
            
            if( !punucs.count(238) || !punucs.count(239) || !punucs.count(240) )
            {
              string msg = "Pu242 correction method of Bignan95 was specified, but problem did not"
              " contain";
              if( !punucs.count(238) )
                msg += " Pu238";
              if( !punucs.count(239) )
                msg += string(punucs.count(238) ? "" : ",") + " Pu239";
              if( !punucs.count(240) )
                msg += string((punucs.count(238) && punucs.count(239)) ? "" : ",") + " Pu240";
              msg += ", as required.";
              
              throw runtime_error( msg );
            }//if( !have_pu238 || !have_pu239 || !have_pu240 )
            
            break;
          }// Bignan95_PWR or Bignan95_BWR
            
          case RelActCalc::PuCorrMethod::ByPu239Only:
          {
            // Need Pu239, plus one other Pu isotope, or I guess we'll accept Am241
            const set<short> punucs = pu_iso_present();
            
            if( punucs.count(242) )
              throw runtime_error( "You can not specify a Pu242 correction method AND include Pu242"
                                  " in the nuclides you are fitting in the spectrum." );
            
            if( !punucs.count(239) || (punucs.size() < 2) )
            {
              string msg = "Pu242 correction method  using Pu239-only was specified, but problem did"
              " not contain";
              
              if( punucs.count(239) )
                msg += " any other Pu nuclides or Am241.";
              else
                msg += " Pu239.";
              
              throw runtime_error( msg );
            }//if( !have_pu238 || !have_pu239 || !have_pu240 )
            break;
          }//case RelActCalc::PuCorrMethod::ByPu239Only:
            
          case RelActCalc::PuCorrMethod::NotApplicable:
            break;
        }//switch( options.pu242_correlation_method )
      };//check_pu_corr_method


      for( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve : options.rel_eff_curves )
        check_pu_corr_method( rel_eff_curve );
    }catch( std::exception &e )
    {
      solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      solution.m_error_message = "Error initializing problem: " + string(e.what());
      
      return solution;
    }// try / catch to
    
    assert( !!spectrum );
    assert( !channel_counts.empty() );
    assert( !channel_count_uncerts.empty() );
    
    // Now use peaks_in_rois and manual stuff to try and estimate initial RelEff and activities
    
    shared_ptr<RelActAutoCostFcn> cost_functor;
    try
    {
      if( cancel_calc && cancel_calc->load() )
        throw runtime_error( "User cancelled calculation." );
      
      auto functor = new RelActAutoCostFcn( options, spectrum, channel_count_uncerts,
                                           input_drf, all_peaks, cancel_calc );
      cost_functor.reset( functor );
    }catch( std::exception &e )
    {
      if( cancel_calc && cancel_calc->load() )
        solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::UserCanceled;
      else
        solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      
      solution.m_error_message = "Error initializing problem: " + string(e.what());
      
      return solution;
    }//try / catch
  
    
    solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    solution.m_drf = cost_functor->m_drf;
    solution.m_spectrum = make_shared<SpecUtils::Measurement>( *spectrum );
    
    solution.m_final_roi_ranges.clear();
    for( const RoiRangeChannels &roi : cost_functor->m_energy_ranges )
      solution.m_final_roi_ranges.push_back( roi );

    // `cost_functor` hasnt had `m_peak_ranges_with_uncert` initialized yet, which we cant do
    //  until after we get our initial activity/rel-eff estimates setup
    const size_t num_pars_initial = cost_functor->number_parameters();
    
    vector<double> parameters( num_pars_initial, 0.0 );
    double *pars = &parameters[0];
    
    vector<int> constant_parameters;
    vector<std::optional<double>> lower_bounds( num_pars_initial ), upper_bounds( num_pars_initial );
    
    
    assert( cost_functor->m_energy_cal && cost_functor->m_energy_cal->valid() );
    for( size_t i = 0; i < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars; ++i )
      parameters[i] = RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset;
    
    if( options.fit_energy_cal )
    {
      static_assert( ((RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars == 2)
                      || (RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars == 3)),
                    "RelActAutoSolution::sm_num_energy_cal_pars only implemented for 2 and 3" );
      
      assert( !cost_functor->m_energy_ranges.empty() );
     
      // For each order we are fitting, we will fit the number of keV the adjustment makes to
      //  the far right-hand side of the spectrum.
      // i.e. To apply the correction to the original energy calibration:
      //
      // Offset:
      //  - For all energy cal types, par[0] will move the whole spectrum left or right that amount
      // Gain:
      //  - Polynomial: divide `parameters[1]` by number of channels, and add to gain.
      //  - FRF: add parameters[1] to gain.
      //  - EnergyCalType::LowerChannelEdge: let i = energy/max_energy, add `parameters[1]*i`
      // Cubic:
      //  - Polynomial: divide `parameters[2]` by number of channels squared, and add to cubic term.
      //  - FRF: add parameters[2] to cubic term.
      //  - EnergyCalType::LowerChannelEdge: not applicable - would need to do math
      //
      // To un-apply:
      
      
      const double lowest_energy = cost_functor->m_energy_ranges.front().lower_energy;
      const double highest_energy = cost_functor->m_energy_ranges.back().upper_energy;
      
      // Completely arbitrary
      if( ((lowest_energy < 200) && (highest_energy > 600))
         || ((lowest_energy < 120) && (highest_energy > 250)) )
      {
        solution.m_fit_energy_cal[0] = true;
        solution.m_fit_energy_cal[1] = true; 
      }else if( highest_energy < 200 )
      {
        //We'll only fit offset
        solution.m_fit_energy_cal[0] = true;
      }else
      {
        // We'll only fit gain
        solution.m_fit_energy_cal[1] = true;
      }
      
      if( solution.m_fit_energy_cal[0] )
      {
        constexpr double offset = RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset;
        constexpr double limit = RelActCalcAuto::RelActAutoSolution::sm_energy_offset_range_keV;
        constexpr double cal_mult = RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
        
        lower_bounds[0] = offset - (limit/cal_mult);
        upper_bounds[0] = offset + (limit/cal_mult);
      }else
      {
        constant_parameters.push_back( 0 );
      }//if( solution.m_fit_energy_cal[0] ) / else
      
      
      if( solution.m_fit_energy_cal[1] )
      {
        constexpr double offset = RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset;
        constexpr double gain_limit = RelActCalcAuto::RelActAutoSolution::sm_energy_gain_range_keV;
        constexpr double cal_mult = RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
        
        lower_bounds[1] = offset - (gain_limit/cal_mult);
        upper_bounds[1] = offset + (gain_limit/cal_mult);
      }else
      {
        constant_parameters.push_back( 1 );
      }//if( solution.m_fit_energy_cal[1] ) / else
      
      if( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars > 2 )
      {
        static_assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 3, "" );
        
        if( solution.m_fit_energy_cal[1]
           && (cost_functor->m_energy_ranges.size() > 4)
           && PeakFitUtils::is_high_res( cost_functor->m_spectrum )
           && (cost_functor->m_energy_cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
        {
          constexpr double offset = RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset;
          constexpr double cubic_limit = RelActCalcAuto::RelActAutoSolution::sm_energy_gain_range_keV;
          constexpr double cal_mult = RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
          
          lower_bounds[2] = offset - (cubic_limit/cal_mult);
          upper_bounds[2] = offset + (cubic_limit/cal_mult);
          solution.m_fit_energy_cal[2] = true;
        }else
        {
          for( size_t i = 2; i < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars; ++i )
          {
            constant_parameters.push_back( static_cast<int>(i) );
            solution.m_fit_energy_cal[i] = false;
          }
        }
      }//if( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars > 2 )
    }else
    {
      for( size_t i = 0; i < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars; ++i )
        constant_parameters.push_back( static_cast<int>(i) );
    }
    
    shared_ptr<const DetectorPeakResponse> res_drf = cost_functor->m_drf;

    
    
#ifndef NDEBUG
    {// Begin some dev checks
      size_t num_rel_eff_par = 0;
      for( const auto &rel_eff_curve : options.rel_eff_curves )
      {
        if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          num_rel_eff_par += 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2;
        else
          num_rel_eff_par += rel_eff_curve.rel_eff_eqn_order + 1;
      }//for( const auto &rel_eff_curve : options.rel_eff_curves )
      
      size_t num_acts_par = 0;
      for( const auto &rel_eff_curve : options.rel_eff_curves )
        num_acts_par += 2*rel_eff_curve.nuclides.size();
      
      assert( cost_functor->m_energy_cal_par_start_index == 0 );
      assert( cost_functor->m_fwhm_par_start_index == RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars );
      assert( cost_functor->m_rel_eff_par_start_index == (cost_functor->m_fwhm_par_start_index + num_parameters(options.fwhm_form)) );
      assert( cost_functor->m_acts_par_start_index == (cost_functor->m_rel_eff_par_start_index + num_rel_eff_par) );
      assert( cost_functor->m_free_peak_par_start_index == (cost_functor->m_acts_par_start_index + num_acts_par) );
      assert( cost_functor->m_skew_par_start_index == (cost_functor->m_free_peak_par_start_index + 2*extra_peaks.size()) );
      assert( (cost_functor->m_skew_par_start_index + 2*PeakDef::num_skew_parameters(options.skew_type)) == cost_functor->number_parameters() );
    }//End some dev checks
#endif //NDEBUG
    
    try
    {
      // We'll convert from one FWHM type to another here; throwing an exception if we fail
      if( !res_drf || !res_drf->hasResolutionInfo() )
        throw runtime_error( "No starting FWHM info was passed in, or derived from the spectrum;"
                             " please select a DRF with FWHM information." );
    
      const RelActCalcAuto::FwhmForm rel_act_fwhm_form = cost_functor->m_options.fwhm_form;
      const DetectorPeakResponse::ResolutionFnctForm drf_fwhm_type = res_drf->resolutionFcnType();
      vector<float> drfpars = res_drf->resolutionFcnCoefficients();
      
      assert( drf_fwhm_type != DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm );
      
      
      bool needToFitOtherType = false;
      
      switch( rel_act_fwhm_form )
      {
        case RelActCalcAuto::FwhmForm::Gadras:
          // We are fitting to the GADRAS functional form
          assert( num_parameters(options.fwhm_form) == 3 );
          needToFitOtherType = (drf_fwhm_type != DetectorPeakResponse::kGadrasResolutionFcn);
        break;
          
        case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
          needToFitOtherType = (drf_fwhm_type != DetectorPeakResponse::kSqrtEnergyPlusInverse);
        break;
          
        case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
          needToFitOtherType = (drf_fwhm_type != DetectorPeakResponse::kConstantPlusSqrtEnergy);
        break;
          
        case RelActCalcAuto::FwhmForm::Polynomial_2:
        case RelActCalcAuto::FwhmForm::Polynomial_3:
        case RelActCalcAuto::FwhmForm::Polynomial_4:
        case RelActCalcAuto::FwhmForm::Polynomial_5:
        case RelActCalcAuto::FwhmForm::Polynomial_6:
          assert( num_parameters(rel_act_fwhm_form) == (static_cast<size_t>(rel_act_fwhm_form)-1) );
          
          needToFitOtherType = ((drf_fwhm_type != DetectorPeakResponse::kSqrtPolynomial)
                                || (drfpars.size() != num_parameters(rel_act_fwhm_form)) );
        break;
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
      
        DetectorPeakResponse::ResolutionFnctForm formToFit;
        const size_t num_fwhm_pars = num_parameters(options.fwhm_form);
        switch( rel_act_fwhm_form )
        {
          case RelActCalcAuto::FwhmForm::Gadras:
            assert( num_fwhm_pars == 3 );
            formToFit = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
            break;
            
          case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
            assert( num_fwhm_pars == 3 );
            formToFit = DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse;
            break;
            
          case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
            assert( num_fwhm_pars == 2 );
            formToFit = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
            break;
         
          case RelActCalcAuto::FwhmForm::Polynomial_2:
          case RelActCalcAuto::FwhmForm::Polynomial_3:
          case RelActCalcAuto::FwhmForm::Polynomial_4:
          case RelActCalcAuto::FwhmForm::Polynomial_5:
          case RelActCalcAuto::FwhmForm::Polynomial_6:
            formToFit = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
            break;
        }//switch( rel_act_fwhm_form )
        
        vector<float> new_sigma_coefs, sigma_coef_uncerts;
        MakeDrfFit::performResolutionFit( fake_peaks, formToFit, highres, static_cast<int>(num_fwhm_pars),
                                         new_sigma_coefs, sigma_coef_uncerts );
      
        assert( new_sigma_coefs.size() == num_fwhm_pars );
        
        drfpars = new_sigma_coefs;
      }//if( needToFitOtherType )
      
      if( drfpars.size() != num_parameters(options.fwhm_form) )
      {
        assert( 0 );
        throw logic_error( "Unexpected num parameters from fit to FWHM function (logic error)." );
      }
      
      for( size_t i = 0; i < drfpars.size(); ++i )
        parameters[cost_functor->m_fwhm_par_start_index + i] = drfpars[i];
    }catch( std::exception &e )
    {
      // We failed to convert from one FWHM type to another - we'll just use some default values
      //  (I dont expect this to happen very often at all)
      solution.m_warnings.push_back( "Failed to create initial FWHM estimation, but will continue anyway: "
                                    + string(e.what()) );
    
      fill_in_default_start_fwhm_pars( parameters, cost_functor->m_fwhm_par_start_index, highres, options.fwhm_form );
    }//try / catch
  
          
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
          peak.m_mean = peak.m_energy;
          peak.m_fwhm = p->gausPeak() ? p->fwhm() : (2.35482 * 0.25 * p->roiWidth());
          peak.m_counts = p->amplitude();
          peak.m_counts_uncert = p->amplitudeUncert();
          peak.m_base_rel_eff_uncert = 0.1; //TODO: do we want this?
          peaks_in_range.push_back( peak );
          
          debug_manual_display_peaks.push_back( p );
        }//if( use_peak )
      }//for( const shared_ptr<const PeakDef> &p : all_peaks )
      
      const double real_time = (foreground && (foreground->real_time() > 0))
                                 ? foreground->real_time() : -1.0f;

      assert( cost_functor->m_nuclides.size() == options.rel_eff_curves.size() );
      if( cost_functor->m_nuclides.size() != options.rel_eff_curves.size() )
        throw logic_error( "Number of nuclides in cost functor does not match number of relative efficiency curves." );

      // Loop over each relative efficiency curve and estimate the starting parameters for the "auto" fit
      //  by doing a rough/course/poor manual fit.
      size_t num_re_nucs_seen = 0, num_re_curve_pars_seen = 0; // To help us keep track of rel eff and act/age index in parameters vector
      for( size_t re_eff_index = 0; re_eff_index < options.rel_eff_curves.size(); ++re_eff_index )
      {
        const double base_rel_eff_uncert = 1.0;
        
        const auto &rel_eff_curve = options.rel_eff_curves[re_eff_index];
        assert( (rel_eff_curve.rel_eff_eqn_order != 0)
               || (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) );
        
        // TODO: add something like `options.correct_for_decay_during_meas` or maybe this option should be a per-nuclide option.
        //const bool correct_for_decay = (real_time > 0.0) ? options.correct_for_decay_during_meas : false;
        const bool correct_for_decay = false;
        
        const size_t this_rel_eff_start = cost_functor->rel_eff_eqn_start_parameter(re_eff_index);
        assert( this_rel_eff_start == (cost_functor->m_rel_eff_par_start_index + num_re_curve_pars_seen) );
        
        const std::vector<NucInputGamma> &nucs = cost_functor->m_nuclides[re_eff_index];
        
        assert( nucs.size() == rel_eff_curve.nuclides.size() );
        if( nucs.size() != rel_eff_curve.nuclides.size() )
          throw logic_error( "Number of nuclides in cost functor does not match number of relative efficiency curve." );
        
        
        // Estimate the initial rel-eff equation
        //  Some things to consider are that we have a really poor estimate of activities right now
        //  and that we want the Rel Eff curve to be 1.0 at the lower energy of the lowest ROI...
        //  Maybe the thing to do is to use the auto-fit peaks, implement the manual fitting back-end
        //   and assigning of gammas to peaks, then use that to first solve for relative activity and
        //   relative efficiency
        bool succesfully_estimated_re_and_ra = false;
        try
        {
          // Lets do an initial setup of parameters, but we'll (hopefully) override this using estimates from
          // the manual rel eff solution.
          switch( rel_eff_curve.rel_eff_eqn_type )
          {
            case RelActCalc::RelEffEqnForm::LnX:
            case RelActCalc::RelEffEqnForm::LnY:
            case RelActCalc::RelEffEqnForm::LnXLnY:
            case RelActCalc::RelEffEqnForm::FramEmpirical:
            {
              // Note that \p rel_eff_order is the number of energy-dependent terms we will fit for in the
              //  relative efficiency equation (i.e., we will also have 1 non-energy dependent term, so we
              //  will fit for (1 + rel_eff_order) parameters).
              //
              //  We'll start with all values for rel eff equation at zero, except maybe the first one
              //  (e.g., we'll start with rel eff line == 1.0 for all energies), and then after we estimate
              //  starting activities, we'll do a little better job estimating things.
              for( size_t rel_eff_index = 0; rel_eff_index <= rel_eff_curve.rel_eff_eqn_order; ++rel_eff_index )
              {
                parameters[this_rel_eff_start + rel_eff_index] = 0.0;
              }
              
              if( (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::LnX)
                 || (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::LnY) )
              {
                parameters[this_rel_eff_start + 0] = 1.0;
              }
              
              break;
            }//case all models except FramPhysicalModel
              
            case RelActCalc::RelEffEqnForm::FramPhysicalModel:
            {
              const auto &self_atten_opt = rel_eff_curve.phys_model_self_atten;
              setup_physical_model_shield_par( lower_bounds, upper_bounds, constant_parameters, pars, this_rel_eff_start, self_atten_opt );
              
              for( size_t ext_ind = 0; ext_ind < rel_eff_curve.phys_model_external_atten.size(); ++ext_ind )
              {
                const size_t start_index = this_rel_eff_start + 2 + 2*ext_ind;
                const auto &opt = rel_eff_curve.phys_model_external_atten[ext_ind];
                setup_physical_model_shield_par( lower_bounds, upper_bounds, constant_parameters, pars, start_index, opt );
              }//for( size_t ext_ind = 0; ext_ind < options.phys_model_external_atten.size(); ++ext_ind )
              
              // Not sure what to do with the b and c values of the Modified Hoerl function ( E^b * c^(1/E) ).
              const size_t b_index = this_rel_eff_start + 2 + 2*rel_eff_curve.phys_model_external_atten.size();
              const size_t c_index = b_index + 1;
              
              pars[b_index] = 0.0;  //(energy/1000)^b - start b at 0, so term is 1.0
              pars[c_index] = 1.0;  //c^(1000/energy) - start c at 1, so term is 1
              if( !rel_eff_curve.phys_model_use_hoerl )
              {
                constant_parameters.push_back( static_cast<int>(b_index) );
                constant_parameters.push_back( static_cast<int>(c_index) );
              }else
              {
                lower_bounds[b_index] = 0.0;
                upper_bounds[b_index] = 2.0;
                lower_bounds[c_index] = 1.0E-6;  //e.x, pow(-0.1889,1000/124.8) is NaN
                upper_bounds[c_index] = 3.0;
              }
              break;
            }//case RelActCalc::RelEffEqnForm::FramPhysicalModel:
          }//switch( solution.m_rel_eff_form )
         
          
          // Now we'll put together the "manual" approximation for R.E. and try to solve it
          vector<RelActCalcManual::SandiaDecayNuc> nuc_sources;
          for( const RelActCalcAuto::NucInputInfo &info : nucs )
          {
            RelActCalcManual::SandiaDecayNuc nucinfo;
            nucinfo.nuclide = info.nuclide;
            nucinfo.age = info.age;
            nucinfo.correct_for_decay_during_meas = correct_for_decay;
            nuc_sources.push_back( nucinfo );
          }//for( const RelActCalcAuto::NucInputInfo &info : rel_eff_curve.nuclides )
          
          
          // TODO: For cbnm9375, fill_in_nuclide_info() and add_nuclides_to_peaks() produce nearly identical results (like tiny rounding errors on yields), but this causes a notable difference in the final "auto" solution, although this manual solution apears the same - really should figure this out - and then get rid of add_nuclides_to_peaks(...) - maybe this is all a testimate to how brittle somethign else is... (20250211: should rechech this, since a number of instabilities have been corrected)
          const double cluster_num_sigma = 1.5;
          const auto peaks_with_nucs = add_nuclides_to_peaks( peaks_in_range, nuc_sources, real_time, cluster_num_sigma );
          
          
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
              if( fabs(floater.energy - p.m_energy) < 1.0*(p.m_fwhm/2.35482) )
                is_floater_peak = true;
            }//for( const RelActCalcAuto::FloatingPeak &floater : extra_peaks )
            
            
            if( !is_floater_peak && !p.m_source_gammas.empty() )
              peaks_with_sources.push_back( p );
          }//for( const auto &p : peaks_with_nucs )
          
          
          size_t manual_rel_eff_order = rel_eff_curve.rel_eff_eqn_order;
          set<string> manual_nucs;
          for( const auto &p : peaks_with_sources )
          {
            for( const auto &l : p.m_source_gammas )
              manual_nucs.insert( l.m_isotope );
          }
          
          int manual_num_peaks = static_cast<int>( peaks_with_sources.size() );
          int manual_num_isos = static_cast<int>( manual_nucs.size() );
          int manual_num_rel_eff = static_cast<int>( manual_rel_eff_order + 1 );
          int num_free_pars = manual_num_peaks - (manual_num_rel_eff + manual_num_isos - 1);
          
          if( (manual_num_peaks - manual_num_isos) < 1 )
            throw runtime_error( "Not enough fit-peaks to perform initial manual rel. eff. estimation of parameters." );
          
          if( num_free_pars < 0 )
            manual_rel_eff_order = manual_num_peaks - manual_num_isos;
          
          RelActCalcManual::RelEffInput manual_input;
          manual_input.peaks = peaks_with_sources;
          manual_input.phys_model_detector = res_drf;
          manual_input.eqn_form = rel_eff_curve.rel_eff_eqn_type;
          manual_input.use_ceres_to_fit_eqn = (manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel);
          manual_input.eqn_order = manual_rel_eff_order;
          manual_input.phys_model_self_atten = rel_eff_curve.phys_model_self_atten;
          manual_input.phys_model_external_attens = rel_eff_curve.phys_model_external_atten;
          manual_input.phys_model_use_hoerl = false;
          
          RelActCalcManual::RelEffSolution manual_solution
                               = RelActCalcManual::solve_relative_efficiency( manual_input );
          
          if( manual_rel_eff_order < rel_eff_curve.rel_eff_eqn_order )
            solution.m_warnings.push_back( "Due to a lack of manually fit peaks, the relative"
                                          " efficiency equation order had to be reduced for initial"
                                          " estimate of relative efficiencies and activities." );
          
          const string re_id_str =  (options.rel_eff_curves.size() > 1) ? ("RE_" + std::to_string(re_eff_index)) : string("");
          cout << "Initial manual " << re_id_str << " estimate:" << endl;
          manual_solution.print_summary( cout );
          
          //ofstream debug_manual_html( "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/initial_manual" + re_id_str + "_estimate.html" );
          //manual_solution.print_html_report( debug_manual_html, options.spectrum_title, spectrum, debug_manual_display_peaks );
          
          //Need to fill out rel eff starting values and rel activities starting values
          
          if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
            throw runtime_error( manual_solution.m_error_message );
          
          assert( (manual_solution.m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel)
                 || (manual_solution.m_rel_eff_eqn_coefficients.size() == (manual_rel_eff_order + 1)) );
          
          if( manual_rel_eff_order != rel_eff_curve.rel_eff_eqn_order )
          {
            size_t num_rel_eff_coefs = rel_eff_curve.rel_eff_eqn_order + 1;
            if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
              num_rel_eff_coefs = (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2);
            
            manual_solution.m_rel_eff_eqn_coefficients.resize( num_rel_eff_coefs, 0.0 );
            manual_solution.m_rel_eff_eqn_covariance.resize( num_rel_eff_coefs );
            for( auto &v : manual_solution.m_rel_eff_eqn_covariance )
              v.resize( num_rel_eff_coefs, 0.0 );
          }//
          
          string rel_eff_eqn_str;
          if( manual_solution.m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            rel_eff_eqn_str = RelActCalc::rel_eff_eqn_text( manual_solution.m_input.eqn_form, manual_solution.m_rel_eff_eqn_coefficients );
            for( size_t i = 0; i < manual_solution.m_rel_eff_eqn_coefficients.size(); ++i )
              parameters[this_rel_eff_start + i] = manual_solution.m_rel_eff_eqn_coefficients[i];
          }else
          {
            const bool html_format = false;
            rel_eff_eqn_str = manual_solution.rel_eff_eqn_txt( html_format );
            
            parameters[this_rel_eff_start + 0] = 0.0;
            parameters[this_rel_eff_start + 1] = 0.0;
            if( rel_eff_curve.phys_model_self_atten )
            {
              size_t manual_index = 0;
              if( rel_eff_curve.phys_model_self_atten->fit_atomic_number )
              {
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                parameters[this_rel_eff_start + 0] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                manual_index += 1;
              }
              
              assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
              parameters[this_rel_eff_start + 1] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
              manual_index += 1;
              
              for( size_t i = 0; i < rel_eff_curve.phys_model_external_atten.size(); ++i )
              {
                const auto &ext_att = rel_eff_curve.phys_model_external_atten[i];
                assert( ext_att );
                if( !ext_att )
                  continue;
                
                if( ext_att->fit_atomic_number )
                {
                  assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                  parameters[this_rel_eff_start + 2 + 2*i + 0] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                  manual_index += 1;
                }
                
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                parameters[this_rel_eff_start + 2 + 2*i + 1] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                manual_index += 1;
              }//for( loop over options.phys_model_external_atten )
              
              const size_t b_index = this_rel_eff_start + 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 0;
              const size_t c_index = this_rel_eff_start + 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 1;
              if( manual_input.phys_model_use_hoerl )
              {
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                parameters[b_index] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                manual_index += 1;
                
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                parameters[c_index] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                manual_index += 1;
              }else
              {
                pars[b_index] = 0.0;  //(energy/1000)^b
                pars[c_index] = 1.0;  //c^(1000/energy)
              }//if( manual_input.phys_model_use_hoerl ) / else
            }//if( options.phys_model_self_atten )
          }//if( no Physical Model ) / else
          cout << "Starting with initial rel. eff. eqn = " << rel_eff_eqn_str << endl;
          
          const double live_time = (spectrum && (spectrum->live_time() > 0)) ? solution.m_spectrum->live_time() : 1.0f;
          
          for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
          {
            const RelActCalcAuto::NucInputInfo &nuc = rel_eff_curve.nuclides[nuc_num];
            const size_t act_index = cost_functor->nuclide_parameter_index(nuc.nuclide, re_eff_index);
            double rel_act = manual_solution.relative_activity( nuc.nuclide->symbol ) / live_time;
            
            /*
#warning "doing terrible dev hack to check seperating rel eff curves"
            if( re_eff_index == 0 )
            {
              if( nuc.nuclide->symbol != "U238" )
              {
                rel_act *= 0.1;
              }else
              {
                rel_act *= 10;
              }
            }else
            {
              if( nuc.nuclide->symbol == "U235" )
              {
                rel_act *= 0.1;
              }else
              {
                rel_act *= 10;
              }
            }
            */

            cout << "Updating initial activity estimate for " << nuc.nuclide->symbol << " from "
            << parameters[act_index] << " to " << rel_act << endl;
            
            assert( cost_functor->m_nuclides[re_eff_index][nuc_num].nuclide == nuc.nuclide );
            
            if( (rel_act > 1.0E-12) && !std::isinf(rel_act) && !std::isnan(rel_act) )
              cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = rel_act;
            else
              cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = 1.0;
            parameters[act_index] = 1.0;

            lower_bounds[act_index] = 0.0;
          }//for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
          
          
          succesfully_estimated_re_and_ra = true;
        }catch( std::exception &e )
        {
          cerr << "Failed to do initial estimate ov RelEff curve: " << e.what() << endl;
          
          solution.m_warnings.push_back( "Initial estimate of relative efficiency curve failed ('"
                                        + string(e.what())
                                        + "'), using a flat line as a starting point" );
        }//try / catch
          
          
        for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
        {
          const RelActCalcAuto::NucInputInfo &nuc = rel_eff_curve.nuclides[nuc_num];
          const size_t act_index = cost_functor->nuclide_parameter_index( nuc.nuclide, re_eff_index );
          const size_t age_index = act_index + 1;
          
          assert( nuc.nuclide );
          assert( nuc.age >= 0.0 );
          if( !nuc.nuclide || nuc.nuclide->isStable() )
            throw runtime_error( "Invalid nuclide" );
          
          if( nuc.age < 0.0 )
            throw runtime_error( "Invalid nuclide (" + nuc.nuclide->symbol + ") age" );
          
          assert( re_eff_index < cost_functor->m_nuclides.size() );
          assert( nuc_num < cost_functor->m_nuclides[re_eff_index].size() );
          assert( cost_functor->m_nuclides[re_eff_index][nuc_num].nuclide == nuc.nuclide );
          
          NucInputGamma &nuc_info = cost_functor->m_nuclides[re_eff_index][nuc_num];
          nuc_info.age_multiple = nuc.age;
          parameters[age_index] = 1.0;

          if( nuc_info.age_multiple <= 1.0*PhysicalUnits::second )
            nuc_info.age_multiple = PeakDef::defaultDecayTime( nuc.nuclide, nullptr );
          
          if( nuc_info.age_multiple <= 1.0*PhysicalUnits::second )
          {
            nuc_info.age_multiple = 1.0;
            parameters[age_index] = nuc.age;
          }

          bool fit_age = nuc.fit_age;
          
          // If ages of nuclides of an element are set to all be the same, we will use the age slot
          //  of the first nuclide to control the age of all nuclides for that element; we'll check
          //  for this here, and if we have a nuclide of the same element preceding \c nuc_num, we'll
          //  fix the age here
          if( (nuc_num > 0) && rel_eff_curve.nucs_of_el_same_age )
          {
            bool found_control_age = false;
            for( size_t i = 0; i < nuc_num; ++i )
            {
              if( rel_eff_curve.nuclides[i].nuclide->atomicNumber == nuc.nuclide->atomicNumber )
              {
                if( (rel_eff_curve.nuclides[i].age != nuc.age) || (rel_eff_curve.nuclides[i].fit_age != nuc.fit_age) )
                  throw runtime_error( "When its specified that all nuclides of the element"
                                      " must have the same age, and same (initial) age value and"
                                      " wether or not to be fit must be specified." );
                
                fit_age = false;
                parameters[age_index] = -1.0; //so we will catch logic error as an assert later.
                nuc_info.age_multiple = -1.0;

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
            if( rel_eff_curve.nucs_of_el_same_age )
            {
              for( size_t i = 0; i < rel_eff_curve.nuclides.size(); ++i )
              {
                if( rel_eff_curve.nuclides[i].nuclide->atomicNumber == nuc.nuclide->atomicNumber )
                  half_life = std::max(half_life, rel_eff_curve.nuclides[i].nuclide->halfLife);
              }
            }//if( options.nucs_of_el_same_age )
            
            double max_age = std::max( 5.0*nuc.age, 15.0*half_life );
            
            // We'll clamp the upper age, to the range where humans may have done the seperation
            // TODO: is there any problem with a nuclide >100 years, that isnt effectively infinite?
            max_age = std::min( max_age, 120*PhysicalUnits::year );

            
            lower_bounds[age_index] = 0.0;
            upper_bounds[age_index] = max_age / nuc_info.age_multiple;
          }else
          {
            constant_parameters.push_back( static_cast<int>(age_index) );
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
            
            assert( cost_functor->m_nuclides[re_eff_index][nuc_num].nuclide == nuc.nuclide );
            if( !top_energy_to_rel_act.empty() )
            {
              parameters[act_index] = 1.0;
              cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = std::get<2>( top_energy_to_rel_act[top_energy_to_rel_act.size()/2] );
              cout << "Setting initial relative activity for " << nuc.nuclide->symbol << " to " << cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple << " bq" << endl;
            }else
            {
              parameters[act_index] = 1.0;
              cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = 100*PhysicalUnits::bq;
              solution.m_warnings.push_back( "Could not estimate a starting activity for "
                                            + nuc.nuclide->symbol + ".  This may be because there are no"
                                            " significant gammas for the nuclide in the selected energy"
                                            " ranges." );
            }
          }//if( !succesfully_estimated_re_and_ra )
          
          lower_bounds[act_index] = 0.0;
        }//for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
        
        num_re_nucs_seen += rel_eff_curve.nuclides.size();
        if( rel_eff_curve.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
          num_re_curve_pars_seen += rel_eff_curve.rel_eff_eqn_order + 1;
        else
          num_re_curve_pars_seen += 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2;
      }//for( const auto &rel_eff_curve : options.rel_eff_curves )
        
      // The parameters will have entries for two sets of peak-skew parameters; one for
      //  the lowest energy of the problem, and one for the highest; in-between we will
      //  scale the energy-dependent skew parameters.  If there is no energy dependence for
      //  a parameter, or the problem doesnt span a large enough energy range to have a
      //  dependence, parameters for the second skew will be fixed garbage values (i.e.,
      //  unused)
      const size_t num_skew_coefs = PeakDef::num_skew_parameters( options.skew_type );
      for( size_t i = 0; i < num_skew_coefs; ++i )
      {
        const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
        double lower, upper, starting, step;
        const bool use = PeakDef::skew_parameter_range( options.skew_type, ct,
                                                        lower, upper, starting, step );
        assert( use );
        if( !use )
          throw logic_error( "Inconsistent skew parameter thing" );
          
        bool fit_energy_dep = PeakDef::is_energy_dependent( options.skew_type, ct );
          
        if( fit_energy_dep && (energy_ranges.size() == 1) )
        {
          // TODO: if statistically significant peaks across less than ~100 keV, then set fit_energy_dep to false
#ifdef _MSC_VER
#pragma message( "if statistically significant peaks across less than ~100 keV, then set fit_energy_dep to false" )
#else
#warning "if statistically significant peaks across less than ~100 keV, then set fit_energy_dep to false"
#endif
          const double dx = (energy_ranges.front().upper_energy - energy_ranges.front().lower_energy);
          fit_energy_dep = (dx > 100.0);
        }
          
        // If we will be fitting an energy dependence, make sure the cost functor knows this
        cost_functor->m_skew_has_energy_dependance |= fit_energy_dep;
          
        const size_t skew_start = cost_functor->m_skew_par_start_index;
        const size_t num_skew_coefs = PeakDef::num_skew_parameters( options.skew_type );
        
        parameters[skew_start + i] = starting;
        lower_bounds[skew_start + i] = lower;
        upper_bounds[skew_start + i] = upper;
          
        // Specify ranges for second set of skew parameters
        if( !fit_energy_dep )
        {
          parameters[skew_start + i + num_skew_coefs] = -999.9;
          constant_parameters.push_back( static_cast<int>(skew_start + i + num_skew_coefs) );
        }else
        {
          parameters[skew_start + i + num_skew_coefs] = starting;
          lower_bounds[skew_start + i + num_skew_coefs] = lower;
          upper_bounds[skew_start + i + num_skew_coefs] = upper;
        }
      }//for( size_t i = 0; i < (num_skew_par/2); ++i )
      
    
      // Floating peaks; one parameter for amplitude, one for FWHM (which will usually be unused)
    for( size_t extra_peak_index = 0; extra_peak_index < extra_peaks.size(); ++extra_peak_index )
    {
      const RelActCalcAuto::FloatingPeak &peak = extra_peaks[extra_peak_index];
      assert( peak.energy > 0.0 );
      if( peak.energy <= 0.0 )
        throw runtime_error( "Invalid floating peak energy." );
      
      const size_t free_peak_start = cost_functor->m_free_peak_par_start_index;
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
      lower_bounds[amp_index] = 0.0;
      
      parameters[fwhm_index] = 1.0;
      
      if( peak.release_fwhm )
      {
        // We'll say the width of this peak can not be less than 25% of what is expected
        //  TODO: we probably also want to put a lower bound of 3 or 4 channels as well
        lower_bounds[fwhm_index] = 0.25;
        
        // TODO: set an upper bound on peak width - this needs to be a convolution of Q value and detector resolution - maybe an input to this function
        //  Right now we'll just say a factor of 4 - which is plenty, unless its caused by a large Q reaction.
        upper_bounds[fwhm_index] = 4.0;
      }else
      {
        constant_parameters.push_back( static_cast<int>(fwhm_index) );
      }
    }//for( size_t extra_peak_index = 0; extra_peak_index < extra_peaks.size(); ++extra_peak_index )
    
    
    // If we are adding an additional uncertainty onto "branching ratios", what we will actually
    //  do is cluster things together that will effectively appear as one peak.  We will then
    //  add each clustered peak as its own residual, to allow it to vary
    if( options.additional_br_uncert > 0.0 )
    {
      const double cluster_num_sigma = 1.5;
      cost_functor->m_peak_ranges_with_uncert = cost_functor->cluster_photopeaks( cluster_num_sigma, parameters );
      
      const size_t num_pars = parameters.size() + cost_functor->m_peak_ranges_with_uncert.size();
      parameters.resize( num_pars, sm_peak_range_uncert_par_offset );
      pars = &parameters[0];
      lower_bounds.resize( num_pars, optional<double>(0.0) );
      upper_bounds.resize( num_pars );
      
      const size_t num_skew_coefs = PeakDef::num_skew_parameters( options.skew_type );
      cost_functor->m_add_br_uncert_start_index = cost_functor->m_skew_par_start_index + 2*num_skew_coefs;
      
      assert( num_pars == cost_functor->number_parameters() );
      assert( num_pars == (cost_functor->m_add_br_uncert_start_index + cost_functor->m_peak_ranges_with_uncert.size()) );
    }//if( options.additional_br_uncert > 0.0 )
    
    const size_t num_pars = cost_functor->number_parameters();
    
    ceres::CostFunction *cost_function = nullptr;
    
    // We cant currently use auto-diff if we are fitting energy calibration, or any nuclide ages, so lets check for that.
    //  TODO: write our own CostFunction class that does numeric differentiation for the energy cal and ages, and auto-diff for the rest
    bool use_auto_diff = !options.fit_energy_cal;
    for( const auto &rel_eff_curve : cost_functor->m_options.rel_eff_curves )
    {
      for( const auto &nuc : rel_eff_curve.nuclides )
        use_auto_diff = (use_auto_diff && !nuc.fit_age);
    }
    
#if( PERFORM_DEVELOPER_CHECKS )
    //Test auto diff vs numerical diff
    auto test_gradients = [constant_parameters,&cost_functor]( RelActAutoCostFcn *fcn, const vector<double> &x ){
      ceres::DynamicAutoDiffCostFunction<RelActAutoCostFcn,32> auto_diff( fcn, ceres::DO_NOT_TAKE_OWNERSHIP );
      auto_diff.SetNumResiduals( static_cast<int>(fcn->number_residuals()) );
      auto_diff.AddParameterBlock( static_cast<int>(fcn->number_parameters()) );
      
      ceres::NumericDiffOptions num_diff_options;
      num_diff_options.relative_step_size = 1E-4;
      auto num_diff = ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn>( fcn, ceres::DO_NOT_TAKE_OWNERSHIP, num_diff_options );
      num_diff.SetNumResiduals( static_cast<int>(fcn->number_residuals()) );
      num_diff.AddParameterBlock( static_cast<int>(fcn->number_parameters()) );
      
      const size_t num_res = fcn->number_residuals();
      const size_t num_pars = fcn->number_parameters();
      //assert( num_res == 1 );
      assert( num_pars == x.size() );
      
      // Compute the residuals and gradients using automatic differentiation
      vector<double> pars = x;
      std::vector<double*> parameter_blocks = { &(pars[0]) };
      std::vector<double> residuals_auto( num_res, 0.0 ), jacobians_auto( num_pars*num_res, 0.0 );
      std::vector<double*> jacobian_ptrs_auto = { &jacobians_auto[0] };
      auto_diff.Evaluate( parameter_blocks.data(), residuals_auto.data(), jacobian_ptrs_auto.data());
      
      
      std::vector<double> residuals_numeric( num_res, 0.0 ), jacobians_numeric( num_pars*num_res, 0.0 );
      std::vector<double*> jacobian_ptrs_numeric = { &(jacobians_numeric[0]) };
      num_diff.Evaluate( parameter_blocks.data(), residuals_numeric.data(), jacobian_ptrs_numeric.data() );
      
      // TODO: should make function to give parameter name!
      
      cout << "Non-equal of numeric and auto-diff Jacobians\n";
      cout << setw(7) << "Index" << setw(7) << "ParName" << setw(10) << "ResNum"
      << setw(12) << "ResVal" << setw(12) << "ParVal"
      << setw(12) << "Auto" << setw(12) << "Numeric"
      << setw(12) << "Diff" << setw(12) << "FracDiff"
      << endl;
      for( size_t i = 0; i < num_pars*num_res; ++i )
      {
        // i = residual_index * num_pars + parameter_index.
        const int par_num = static_cast<int>(i % num_pars);
        const string par_name = cost_functor->parameter_name(par_num);
        
        // Const parameters may not be used, so we'll ignore them.
        const auto const_pos = std::find(begin(constant_parameters), end(constant_parameters), par_num );
        if( const_pos != end(constant_parameters) )
          continue;
        
        const size_t residual_num = (i - par_num) / num_pars;
        const double diff = fabs(jacobians_auto[i] - jacobians_numeric[i]);
        const double frac_diff = diff / std::max( fabs(jacobians_auto[i]), fabs(jacobians_numeric[i]) );
        if( (diff > 1.0E-18) && (frac_diff > 1.0E-4) )  //This is totally arbitrary
        {
          cout << setprecision(4)
          << setw(7) << i << setw(10) << par_name << setw(7) << residual_num
          << setw(12) << residuals_numeric[residual_num] << setw(12) << x[par_num]
          << setw(12) << jacobians_auto[i] << setw(12) << jacobians_numeric[i]
          << setw(12) << diff << setw(12) << frac_diff
          << endl;
        }
      }
      cout << "---" << endl;
    };//test_gradients(...)
    
    //test_gradients( cost_functor.get(), parameters );
#endif
    
    if( use_auto_diff )
    {
      // It looks like using Jets<double,32> is fastest for an example problem, however it really
      //  doesnt seem to matter much (Wall times of {0.1728, 0.1789, 0.1665, 0.1608}, for strides of
      //  {4, 8, 16, 32}, respectively
      ceres::DynamicAutoDiffCostFunction<RelActAutoCostFcn,32> *dyn_auto_diff_cost_function
            = new ceres::DynamicAutoDiffCostFunction<RelActAutoCostFcn,32>( cost_functor.get(),
                                                                      ceres::DO_NOT_TAKE_OWNERSHIP );
      dyn_auto_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
      dyn_auto_diff_cost_function->AddParameterBlock( static_cast<int>(cost_functor->number_parameters()) );

      cost_function = dyn_auto_diff_cost_function;
    }else
    {
      ceres::NumericDiffOptions num_diff_options;
      num_diff_options.relative_step_size = 1E-2;
      ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn> *dyn_num_diff_cost_function
            = new ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn>( cost_functor.get(),
                                                ceres::DO_NOT_TAKE_OWNERSHIP, num_diff_options );
      dyn_num_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
      dyn_num_diff_cost_function->AddParameterBlock( static_cast<int>(cost_functor->number_parameters()) );

      cost_function = dyn_num_diff_cost_function;
    }//if( use_auto_diff ) / else
    
    
    //auto cost_function = new ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn>( cost_functor.get(),
    //                                                                ceres::DO_NOT_TAKE_OWNERSHIP );
    //cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
    //cost_function->AddParameterBlock( static_cast<int>(cost_functor->number_parameters()) );
    
    
    ceres::Problem problem;
    ceres::LossFunction *lossfcn = nullptr;
    // For an example problem, the loss function didnt seem to have a impact, or if they did, it wasnt good.
    //lossfcn = new ceres::HuberLoss( 5.0);  //The Huber loss function is quadratic for small residuals and linear for large residuals - probably what we would want to use
    //lossfcn = new ceres::CauchyLoss(15.0); //The Cauchy loss function is less sensitive to large residuals than the Huber loss.
    //lossfcn = new ceres::SoftLOneLoss(5.0);
    //lossfcn = new ceres::TukeyLoss(10.0); //Quadratic for small residuals and zero for large residuals - not good if initial guess isnt great
    
    problem.AddResidualBlock( cost_function, lossfcn, pars );
    problem.AddParameterBlock( pars, static_cast<int>(num_pars) );
    
    if( !constant_parameters.empty() )
    {
      ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_pars), constant_parameters );
      problem.SetManifold( pars, subset_manifold );
    }
    
    
    for( size_t i = 0; i < num_pars; ++i )
    {
      if( lower_bounds[i].has_value() )
        problem.SetParameterLowerBound(pars, static_cast<int>(i), *lower_bounds[i] );
      if( upper_bounds[i].has_value() )
        problem.SetParameterUpperBound(pars, static_cast<int>(i), *upper_bounds[i] );
    }//for( size_t i = 0; i < num_pars; ++i )
    
    
    // Okay - we've set our problem up
    ceres::Solver::Options ceres_options;
    ceres_options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH
    ceres_options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG
    
    //use_nonmonotonic_steps: set to true to allow to "jump over boulders".
    // Using function_tolerance=1e-7, initial_trust_region_radius=1E5, using use_nonmonotonic_steps
    //  false: 52.88% enrich, 261 function calls, 64 iterations, Chi2=0.664167 (terminate due to function tol reached)
    //  true:  52.85% enrich, 126 function calls, 23 iterations, Chi2=0.664162 (terminate due to function tol reached)
    ceres_options.use_nonmonotonic_steps = true;
    ceres_options.max_consecutive_nonmonotonic_steps = 5;
    
    
    // Trust region minimizer settings.
    //
    // Initial trust region has notable impact on getting the right answer; it also has a notable
    //  impact on computation time.
    //  Default value is 1E4.
    //  Notes for an example U problem (50% enrich, two objects), using auto-differentiation, and function_tolerance=1E-9, use_nonmonotonic_steps=false
    //  - 1e4: 42.19% enrich,  87 function calls,  16 iterations, Chi2=0.78602  (terminate due to function tol reached)
    //  - 1e5: 52.87% enrich, 479 function calls, 127 iterations, Chi2=0.664166 (terminate due to function tol reached)
    //  - 1e6: 52.89% enrich, 856 function calls, 232 iterations, Chi2=0.664172 (terminate due to function tol reached)
    //
    //  According to a LLM, the trust region is an area in the parameter space.
    //   So we will need to revisit this value if we scale the various parameters to be closer to reasonable.
    ceres_options.initial_trust_region_radius = 1e5;
    
    ceres_options.max_trust_region_radius = 1e16;
    
    // Minimizer terminates when the trust region radius becomes smaller than this value.
    ceres_options.min_trust_region_radius = 1e-32;
    // Lower bound for the relative decrease before a step is accepted.
    ceres_options.min_relative_decrease = 1e-3;
    
    
    ceres_options.linear_solver_type = ceres::DENSE_QR; //ceres::DENSE_SCHUR, ceres::DENSE_NORMAL_CHOLESKY, ceres::ITERATIVE_SCHUR
    
    ceres_options.minimizer_progress_to_stdout = true; //true;
    ceres_options.logging_type = ceres::PER_MINIMIZER_ITERATION;
    ceres_options.max_num_iterations = 1000;
    //ceres_options.max_solver_time_in_seconds = 120.0;
    
    // Changing function_tolerance from 1e-9 to 1e-7, and using initial_trust_region_radius=1E5, and use_nonmonotonic_steps=false
    //  - 52.88% enrich, 261 function calls, 64 iterations, Chi2=0.664167 (terminate due to function tol reached)
    ceres_options.function_tolerance = 1e-7; //1e-9;
    
    ceres_options.gradient_tolerance = 1.0E-4*ceres_options.function_tolerance; // Per documentation of `gradient_tolerance`
    // parameter_tolerance seems to be what terminates the minimization on some example problems.
    //  It looks to be:
    //    |step|_2 <= parameter_tolerance * ( |x|_2 +  parameter_tolerance)
    //  Where `|step|_2` and `|x|_2` indicate the norm of the parameters - e.g.,
    ceres_options.parameter_tolerance = 1e-11; //Default value is 1e-8.  Using 1e-11 gets the cost change down to 0.1 in Chi2, for an example Physical Model
    
    // TODO: there are a ton of ceres::Solver::Options that might be useful for us to set
    
    std::unique_ptr<CheckCeresTerminateCallback> terminate_callback;
    if( cancel_calc )
    {
      terminate_callback.reset( new CheckCeresTerminateCallback(cancel_calc) );
      ceres_options.callbacks.push_back( terminate_callback.get() );
    }
    
    
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
        solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
        solution.m_error_message += "The L-M solving failed - no convergence.";
        break;
        
      case ceres::FAILURE:
        solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
        solution.m_error_message += "The L-M solving failed.";
        break;
        
      case ceres::USER_FAILURE:
        if( cancel_calc && cancel_calc->load() )
        {
          solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::UserCanceled;
          solution.m_error_message += "Calculation was cancelled.";
        }else
        {
          // I dont think we should get here.
          assert( 0 );
          solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
          solution.m_error_message += "The L-M solving failed.";
        }
        break;
    }//switch( summary.termination_type )
    
    const bool success = (solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success);
    
    solution.m_num_function_eval_solution = static_cast<int>( cost_functor->m_ncalls );
    solution.m_num_microseconds_in_eval = static_cast<int>( cost_functor->m_nanoseconds_spent_in_eval / 1000 );  //convert from nanoseconds to micro
    
    {
      const auto now_time = std::chrono::high_resolution_clock::now();
      const auto dt = 1.0*std::chrono::duration_cast<std::chrono::nanoseconds>(now_time - start_time).count();
      const double frac = cost_functor->m_nanoseconds_spent_in_eval / dt;
      
      cout << "Spent " << 0.001*cost_functor->m_nanoseconds_spent_in_eval
      << " us in eval, and a total time of " << 0.001*dt << " us in fcnt - this is "
      << frac << " fraction of the time" << endl;
    }
    
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //
    cov_options.num_threads = ceres_options.num_threads;
    
    vector<double> uncertainties( num_pars, 0.0 ), uncerts_squared( num_pars, 0.0 );
    
    if( success )
    {
      ceres::Covariance covariance(cov_options);
      vector<pair<const double*, const double*> > covariance_blocks;
      covariance_blocks.emplace_back( pars, pars );
      
      //auto add_cov_block = [&covariance_blocks, &pars]( size_t start, size_t num ){
      //  for( size_t i = start; i < (start + num); ++i )
      //    for( size_t j = start; j < i; ++j )
      //      covariance_blocks.push_back( make_pair( pars + i, pars + j ) );
      //};
      
      //add_cov_block( rel_eff_start, num_rel_eff_par );
      //add_cov_block( acts_start, num_acts_par );
      //add_cov_block( fwhm_start, num_fwhm_pars );
      solution.m_covariance.clear();
      solution.m_final_uncertainties.clear();
      
      if( !covariance.Compute(covariance_blocks, &problem) )
      {
        cerr << "Failed to compute final covariances!" << endl;
        solution.m_warnings.push_back( "Failed to compute final covariances." );
      }else
      {
        // row-major order: the elements of the first row are consecutively given, followed by second
        //                  row contents, etc
        vector<double> row_major_covariance( num_pars * num_pars );
        const vector<const double *> const_par_blocks( 1, pars );
        
        const bool success = covariance.GetCovarianceMatrix( const_par_blocks, row_major_covariance.data() );
        assert( success );
        if( !success )
          throw runtime_error( "Failed to get covariance matrix - maybe didnt add all covariance blocks?" );
        
        solution.m_covariance.resize( num_pars, vector<double>(num_pars, 0.0) );
        
        for( size_t row = 0; row < num_pars; ++row )
        {
          for( size_t col = 0; col < num_pars; ++col )
            solution.m_covariance[row][col] = row_major_covariance[row*num_pars + col];
        }//for( size_t row = 0; row < num_nuclides; ++row )
        
        for( size_t i = 0; i < num_pars; ++i )
        {
          uncerts_squared[i] = solution.m_covariance[i][i];
          if( uncerts_squared[i] >= 0.0 )
            uncertainties[i] = sqrt( uncerts_squared[i] );
          else
            uncertainties[i] = std::numeric_limits<double>::quiet_NaN();
        }
        
        solution.m_final_uncertainties = uncertainties;
        
        // TODO: check the covariance is actually defined right (like not swapping row/col, calling the right places, etc).
        auto get_cov_block = [&solution]( const size_t start, const size_t num, vector<vector<double>> &cov ){
          cov.clear();
          cov.resize( num, vector<double>(num,0.0) );
          for( size_t i = start; i < (start + num); ++i )
          {
            for( size_t j = start; j < (start + num); ++j )
              cov[j-start][i-start] = solution.m_covariance[i][j];
          }
        };
        
        solution.m_rel_eff_covariance.resize( cost_functor->m_options.rel_eff_curves.size() );
        solution.m_rel_act_covariance.resize( cost_functor->m_options.rel_eff_curves.size() );
        for( size_t rel_eff_index = 0; rel_eff_index < cost_functor->m_options.rel_eff_curves.size(); ++rel_eff_index )
        {
          const size_t rel_eff_start = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
          const size_t num_rel_eff_par = cost_functor->rel_eff_eqn_num_parameters(rel_eff_index);
          get_cov_block( rel_eff_start, num_rel_eff_par, solution.m_rel_eff_covariance[rel_eff_index] );
        
          const size_t acts_start = cost_functor->rel_act_start_parameter( rel_eff_index );
          const size_t num_acts_par = cost_functor->rel_act_num_parameters( rel_eff_index );
          get_cov_block( acts_start, num_acts_par, solution.m_rel_act_covariance );
        }
        
        const size_t fwhm_start = cost_functor->m_fwhm_par_start_index;
        const size_t num_fwhm_pars = num_parameters(options.fwhm_form);
        get_cov_block( fwhm_start, num_fwhm_pars, solution.m_fwhm_covariance );
      }//if( we failed to get covariance ) / else
    }//if( solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success )
    
    solution.m_num_function_eval_total = static_cast<int>( cost_functor->m_ncalls );
  
    solution.m_final_parameters = parameters;
    
    solution.m_parameter_names.resize( parameters.size() );
    solution.m_parameter_were_fit.resize( parameters.size(), true );
    for( size_t i = 0; i < parameters.size(); ++i )
    {
      solution.m_parameter_names[i] = cost_functor->parameter_name( i );
      solution.m_parameter_were_fit[i] = std::find( begin(constant_parameters), end(constant_parameters), static_cast<int>(i) ) == end(constant_parameters);
    }


    for( size_t i = 0; i < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars; ++i )
      solution.m_energy_cal_adjustments[i] = parameters[i];
    
    // We will use `new_cal` to move peaks from their true energy, to the spectrums original
    //  energy calibration - so we will un-apply the energy calibration correction
    shared_ptr<const SpecUtils::EnergyCalibration> new_cal = solution.m_foreground->energy_calibration();
    if( success && options.fit_energy_cal )
    {
      shared_ptr<SpecUtils::EnergyCalibration> modified_new_cal = solution.get_adjusted_energy_cal();
      solution.m_spectrum->set_energy_calibration( modified_new_cal );
      new_cal = modified_new_cal;
    }//if( options.fit_energy_cal )
    
    if( options.additional_br_uncert > 0.0 )
    {
      solution.m_add_br_uncert_start_index = cost_functor->m_add_br_uncert_start_index;
      solution.m_peak_ranges_with_uncert = cost_functor->m_peak_ranges_with_uncert;
    }//if( options.additional_br_uncert > 0.0 )
  
    vector<PeakDef> fit_peaks;
    for( const auto &range : cost_functor->m_energy_ranges )
    {
      PeaksForEnergyRange these_peaks = cost_functor->peaks_for_energy_range( range, parameters );
      fit_peaks.insert( end(fit_peaks), begin(these_peaks.peaks), end(these_peaks.peaks) );
    }
    
    std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMean );
    
    solution.m_fit_peaks_in_spectrums_cal = fit_peaks;
    
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

    assert( cost_functor->m_nuclides.size() == cost_functor->m_options.rel_eff_curves.size() );

    solution.m_rel_eff_forms.clear();
    solution.m_rel_eff_forms.resize( cost_functor->m_options.rel_eff_curves.size(), RelActCalc::RelEffEqnForm::LnX );

    solution.m_rel_eff_coefficients.clear();
    solution.m_rel_eff_coefficients.resize( cost_functor->m_options.rel_eff_curves.size() );

    solution.m_rel_activities.clear();
    solution.m_rel_activities.resize( cost_functor->m_options.rel_eff_curves.size() );
    
    for( size_t rel_eff_index = 0; rel_eff_index < cost_functor->m_options.rel_eff_curves.size(); ++rel_eff_index )
    {
      const vector<NucInputGamma> &input_nuclides = cost_functor->m_nuclides[rel_eff_index];
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      
      const size_t rel_eff_start = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
      const size_t num_rel_eff_par = cost_functor->rel_eff_eqn_num_parameters(rel_eff_index);
      const auto rel_eff_iter = begin(parameters) + rel_eff_start;
      solution.m_rel_eff_coefficients[rel_eff_index].insert( end(solution.m_rel_eff_coefficients[rel_eff_index]),
                                      rel_eff_iter, rel_eff_iter + num_rel_eff_par );
      solution.m_rel_eff_forms[rel_eff_index] = rel_eff_curve.rel_eff_eqn_type;


      for( size_t act_index = 0; act_index < input_nuclides.size(); ++act_index )
      {
        const NucInputGamma &nuc_input = input_nuclides[act_index];
      
        RelActCalcAuto::NuclideRelAct nuc_output;
        nuc_output.nuclide = nuc_input.nuclide;
        nuc_output.age = cost_functor->age( nuc_input.nuclide, rel_eff_index, parameters );
        nuc_output.age_was_fit = nuc_input.fit_age;
        nuc_output.rel_activity = cost_functor->relative_activity( nuc_input.nuclide, rel_eff_index, parameters );
      
        nuc_output.age_uncertainty = cost_functor->age( nuc_input.nuclide, rel_eff_index, uncertainties );
        nuc_output.rel_activity_uncertainty = cost_functor->relative_activity( nuc_input.nuclide, rel_eff_index, uncertainties );
      
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
          const size_t par_act_index = cost_functor->nuclide_parameter_index( nuc_input.nuclide, rel_eff_index );
          const size_t par_age_index = par_act_index + 1;
          solution.m_warnings.push_back( "Variance for nuclide " + nuc_input.nuclide->symbol
                                      + " age, is invalid ("
                                      + std::to_string(uncerts_squared[par_age_index]) + ")" );
        }
      
        if( IsNan(nuc_output.rel_activity_uncertainty) )
        {
          const size_t par_act_index = cost_functor->nuclide_parameter_index( nuc_input.nuclide, rel_eff_index );
          solution.m_warnings.push_back( "Variance for nuclide " + nuc_input.nuclide->symbol
                                      + " activity, is invalid ("
                                      + std::to_string(uncerts_squared[par_act_index]) + ")" );
        }

        solution.m_rel_activities[rel_eff_index].push_back( nuc_output );
      }//for( size_t act_index = 0; act_index < cost_functor->m_nuclides.size(); ++ act_index )
    }//for( size_t rel_eff_index = 0; rel_eff_index < cost_functor->m_options.rel_eff_curves.size(); ++rel_eff_index )
    
    
    assert( cost_functor->m_options.rel_eff_curves.size() == solution.m_rel_activities.size() );
    if( cost_functor->m_options.rel_eff_curves.size() != solution.m_rel_activities.size() )
      throw logic_error( "size(rel_eff_curves) != size(solution.m_rel_activities)" );
    
    // If we want to correct for Pu242, we wont alter solution.m_rel_activities, but place the
    //  corrected Pu mass fractions in solution.m_corrected_pu
    solution.m_corrected_pu.clear();
    solution.m_corrected_pu.resize( cost_functor->m_options.rel_eff_curves.size() );
    for( size_t rel_eff_index = 0; rel_eff_index < cost_functor->m_options.rel_eff_curves.size(); ++rel_eff_index )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      const vector<RelActCalcAuto::NuclideRelAct> &rel_acts = solution.m_rel_activities[rel_eff_index];
      
      if( rel_eff_curve.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
      {
        try
        {
          RelActCalc::Pu242ByCorrelationInput raw_rel_masses;
          
          double pu_total_mass = 0.0, raw_rel_mass = 0.0;
          for( const RelActCalcAuto::NuclideRelAct &nuc : rel_acts )
          {
            assert( nuc.nuclide );
            if( !nuc.nuclide || (nuc.nuclide->atomicNumber < 94) )
              continue;
            
            const double rel_mass = nuc.rel_activity / nuc.nuclide->activityPerGram();
            
            if( nuc.nuclide->atomicNumber == 94 )
            {
              pu_total_mass += rel_mass;
              switch( nuc.nuclide->massNumber )
              {
                case 238: raw_rel_masses.pu238_rel_mass = rel_mass; break;
                case 239: raw_rel_masses.pu239_rel_mass = rel_mass; break;
                case 240: raw_rel_masses.pu240_rel_mass = rel_mass; break;
                case 241: raw_rel_masses.pu241_rel_mass = rel_mass; break;
                case 242:
                  assert( 0 );
                  throw std::logic_error( "Pu242 shouldnt be in the input nuclides if a Pu242"
                                         " correlation correction method was specified." );
                  break;
                default:  raw_rel_masses.other_pu_mass = rel_mass;  break;
              }//switch( nuc.nuclide->massNumber )
            }else if( (nuc.nuclide->atomicNumber == 95) && (nuc.nuclide->massNumber == 241) )
            {
              // TODO: need to account for half-life of Am-241 and back decay the equivalent Pu241 mass
              pu_total_mass += (242.0/241.0) * rel_mass;
              raw_rel_masses.am241_rel_mass = rel_mass;
            }
          }//for( const NuclideRelAct &nuc : m_rel_activities )
          
          // We dont have to divide by `pu_total_mass`, but we will, just for debuging.
          raw_rel_mass /= pu_total_mass;
          raw_rel_masses.pu238_rel_mass /= pu_total_mass;
          raw_rel_masses.pu239_rel_mass /= pu_total_mass;
          raw_rel_masses.pu240_rel_mass /= pu_total_mass;
          raw_rel_masses.pu241_rel_mass /= pu_total_mass;
          raw_rel_masses.am241_rel_mass /= pu_total_mass;
          raw_rel_masses.other_pu_mass  /= pu_total_mass;
          
          const RelActCalc::Pu242ByCorrelationOutput corr_output
          = RelActCalc::correct_pu_mass_fractions_for_pu242( raw_rel_masses,
                                                            rel_eff_curve.pu242_correlation_method );
          
          if( !corr_output.is_within_range )
            solution.m_warnings.push_back( "The fit Pu enrichment is outside range validated in the"
                                          " literature for the Pu242 correction by correlation." );
          
          solution.m_corrected_pu[rel_eff_index].reset( new RelActCalc::Pu242ByCorrelationOutput(corr_output) );
        }catch( std::exception &e )
        {
          solution.m_warnings.push_back( "Correcting for Pu242 content failed: " + string(e.what()) );
        }//try / catch
      }//if( options.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
    }//for( size_t rel_eff_index = 0; rel_eff_index < cost_functor->m_options.rel_eff_curves.size(); ++rel_eff_index )
    
    
    solution.m_fwhm_form = options.fwhm_form;
    solution.m_fwhm_coefficients.clear();
    
    for( size_t i = 0; i < num_parameters(options.fwhm_form); ++i )
      solution.m_fwhm_coefficients.push_back( parameters[cost_functor->m_fwhm_par_start_index + i] );
    
    // If we are fitting the energy cal, we need to adjust the ROI ranges to the spectrum cal
    solution.m_final_roi_ranges_in_spectrum_cal.clear();
    for( const RoiRangeChannels &roi : cost_functor->m_energy_ranges )
    {
      RelActCalcAuto::RoiRange roi_range = roi;
      if( options.fit_energy_cal )
      {
        roi_range.lower_energy = cost_functor->apply_energy_cal_adjustment( roi_range.lower_energy, parameters );
        roi_range.upper_energy = cost_functor->apply_energy_cal_adjustment( roi_range.upper_energy, parameters ); 
      } 

      solution.m_final_roi_ranges_in_spectrum_cal.push_back( roi_range );
    }//for( const RoiRangeChannels &roi : cost_functor->m_energy_ranges )


    for( size_t i = 0; i < cost_functor->m_options.floating_peaks.size(); ++i )
    {
      const size_t amp_index = cost_functor->m_free_peak_par_start_index + 2*i + 0;
      const size_t fwhm_index = amp_index + 1;
      
      RelActCalcAuto::FloatingPeakResult peak;
      peak.energy = cost_functor->m_options.floating_peaks[i].energy;
      peak.amplitude = parameters[amp_index];
      peak.fwhm = parameters[fwhm_index];
      
      if( !cost_functor->m_options.floating_peaks[i].release_fwhm )
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
    }//for( size_t i = 0; i < cost_functor->m_options.floating_peaks.size(); ++i )
    
    
    // Fill out PhysicalModel stuff, if we are fitting that
    solution.m_phys_model_results.clear();
    solution.m_phys_model_results.resize( cost_functor->m_options.rel_eff_curves.size() );
    const size_t num_rel_eff_curves = cost_functor->m_options.rel_eff_curves.size();
    for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      const size_t this_rel_eff_start_index = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
      const vector<double> &rel_eff_coefficients = solution.m_rel_eff_coefficients[rel_eff_index];
      
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        assert( rel_eff_coefficients.size() == (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2) );
        
        RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo phys_model_result;
        
        auto get_shield_info = [&]( const RelActCalc::PhysicalModelShieldInput &input,
                                   const size_t start_index )
          -> RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo {
          
          RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo answer;
          
          answer.material = input.material;
          answer.atomic_number_was_fit = input.fit_atomic_number;
          if( input.fit_atomic_number )
          {
            assert( !input.material );
            answer.atomic_number = rel_eff_coefficients[start_index + 0] * RelActCalc::ns_an_ceres_mult;
            assert( (answer.atomic_number >= 1.0) && (answer.atomic_number <= 98.0) );
            
            if( uncertainties.size() > (this_rel_eff_start_index + start_index) )
              answer.atomic_number_uncert = uncertainties[this_rel_eff_start_index + start_index + 0] * RelActCalc::ns_an_ceres_mult;
          }else if( !input.material )
          {
            assert( (rel_eff_coefficients[start_index + 0] * RelActCalc::ns_an_ceres_mult)
                   == input.atomic_number );
            answer.atomic_number = input.atomic_number;
            assert( (answer.atomic_number >= 1.0) && (answer.atomic_number <= 98.0) );
          }
          
          if( !answer.material && ((answer.atomic_number < 1.0) || (answer.atomic_number > 98.0)) )
            throw std::logic_error( "AN out of range" );
          
          answer.areal_density_was_fit = input.fit_areal_density;
          if( input.fit_areal_density )
          {
            answer.areal_density = rel_eff_coefficients[start_index + 1] * PhysicalUnits::g_per_cm2;
            assert( answer.areal_density >= 0.0 );
            
            if( uncertainties.size() > (this_rel_eff_start_index + start_index + 1) )
              answer.areal_density_uncert = uncertainties[this_rel_eff_start_index + start_index + 1] * PhysicalUnits::g_per_cm2; //or sqrt( m_rel_eff_covariance[start_index + 1][start_index + 1] )
          }else
          {
            assert( rel_eff_coefficients[start_index + 1] == (input.areal_density / PhysicalUnits::g_per_cm2) );
            answer.areal_density = input.areal_density;
          }
          
          return answer;
        };//auto get_shield_info lamda to fill out
        
        
        if( rel_eff_curve.phys_model_self_atten
           && (rel_eff_curve.phys_model_self_atten->material
               || ((rel_eff_curve.phys_model_self_atten->atomic_number >= 1.0) && (rel_eff_curve.phys_model_self_atten->atomic_number <= 98.0))
               || rel_eff_curve.phys_model_self_atten->fit_atomic_number ) )
        {
          phys_model_result.self_atten = get_shield_info( *rel_eff_curve.phys_model_self_atten, 0 );
        }//if( options.phys_model_self_atten )
        
        for( size_t ext_ind = 0; ext_ind < rel_eff_curve.phys_model_external_atten.size(); ++ext_ind )
        {
          auto res = get_shield_info( *rel_eff_curve.phys_model_external_atten[ext_ind], 2 + 2*ext_ind );
          phys_model_result.ext_shields.push_back( std::move(res) );
        }//for( loop over external attens )
        
        if( rel_eff_curve.phys_model_use_hoerl )
        {
          // Modified Hoerl corrections only present if fitting the Hoerl function was selected
          const size_t b_index = 2 + 2*rel_eff_curve.phys_model_external_atten.size();
          const size_t c_index = b_index + 1;
          assert( c_index < rel_eff_coefficients.size() );
          
          phys_model_result.hoerl_b = rel_eff_coefficients[b_index];
          phys_model_result.hoerl_c = rel_eff_coefficients[c_index];
          
          if( uncertainties.size() > (this_rel_eff_start_index + b_index) )
            phys_model_result.hoerl_b_uncert = uncertainties[this_rel_eff_start_index + b_index];
          if( uncertainties.size() > (this_rel_eff_start_index + c_index) )
            phys_model_result.hoerl_b_uncert = uncertainties[this_rel_eff_start_index + c_index];
        }//if( options.phys_model_use_hoerl )
        
        solution.m_phys_model_results[rel_eff_index] = std::move(phys_model_result);
      }//if( FramPhysicalModel )
    }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    
    
    // We will fit the "best" peak amplitudes.  To do this we will cluster the peaks into regions
    //  and fit them on a ROI by ROI basis; we will only fit the amplitude multiple in each region.
    //
    // blah blah blah
    //  For plotting the Rel. Eff. points
    {
     // const double cluster_num_sigma = 1.5;
     // cost_functor->m_peak_ranges_with_uncert = cost_functor->cluster_photopeaks( cluster_num_sigma, parameters );
     
    }
    
    /*
    if( options.additional_br_uncert > 0.0 )
    {
      cout << "ROIS were adjusted as follows:" << endl;
      for( size_t roi_index = 0; roi_index < cost_functor->m_peak_ranges_with_uncert.size(); ++roi_index )
      {
        const auto range = cost_functor->m_peak_ranges_with_uncert[roi_index];
        const size_t par_index = cost_functor->m_add_br_uncert_start_index + roi_index;
        const double multiple = solution.m_final_parameters[par_index]/sm_peak_range_uncert_par_offset;
        cout << "\t[" << setprecision(4) << setw(6) << range.first << "," << setprecision(4) << setw(5) << range.second << "]: "
        << multiple << endl;
      }
      cout << endl;
    }//if( options.additional_br_uncert > 0.0 )
     */
    
    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    cost_functor->eval( parameters, residuals.data() );
    solution.m_chi2 = 0.0;
    for( const double v : residuals )
      solution.m_chi2 += v*v;
    
    // TODO: need to setup the DOF
    solution.m_dof = 0;
    solution.m_warnings.push_back( "Not currently calculating DOF - need to implement" );
    //solution.m_dof = residuals.size() - ;
    
    //solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
    
    return solution;
  }//RelActCalcAuto::RelActAutoSolution solve_ceres( ceres::Problem &problem )
  
  
  vector<pair<double,double>> cluster_photopeaks( const double cluster_num_sigma, const std::vector<double> &x ) const
  {
    // TODO: this function is by no means optimized - it does stupid O(N^2) things.
    assert( m_spectrum );
    if( !m_spectrum || (m_spectrum->num_gamma_channels() < 7) )
      throw std::runtime_error( "Valid spectrum is needed for clustering gammas." );
    
    const double live_time = (m_spectrum && (m_spectrum->live_time() > 0.0))
                              ? m_spectrum->live_time() : 1.0;
    
    // First we will get all gammas in our rois
    vector<pair<double,double>> gammas_by_counts;  //pair<energy,counts>
    
    const auto lessThanByEnergy = []( const pair<double,double> &lhs, const pair<double,double> &rhs ){
      return lhs.first < rhs.first;
    };
    
    const auto moreThanByNumPerSecond = []( const pair<double,double> &lhs, const pair<double,double> &rhs ){
      return lhs.second > rhs.second;
    };
    
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    
    for( size_t rel_eff_index = 0; rel_eff_index < m_nuclides.size(); ++rel_eff_index )
    {
      const vector<NucInputGamma> &rel_rff_nucs = m_nuclides[rel_eff_index];
      
      for( const NucInputGamma &nuc_input : rel_rff_nucs )
      {
        for( const NucInputGamma::EnergyYield &line_info : nuc_input.nominal_gammas )
        {
          const double energy = line_info.energy;
          
          bool in_a_roi = false;
          for( size_t i = 0; !in_a_roi && (i < m_energy_ranges.size()); ++i )
          {
            in_a_roi = ((energy >= m_energy_ranges[i].lower_energy)
                        && (energy <= m_energy_ranges[i].upper_energy));
          }
          
          if( !in_a_roi )
            continue;
          
          const double rel_act = relative_activity( nuc_input.nuclide, rel_eff_index, x );
          const double rel_eff = relative_eff( energy, rel_eff_index, x );
          
          const double counts = rel_act * rel_eff * line_info.yield * live_time;
          gammas_by_counts.emplace_back( energy, counts );
        }//for( const NucInputGamma::EnergyYield &line_info : nuc_input.nominal_gammas )
      }//for( const NucInputGamma &nuc_input : m_nuclides )
    }//
    
    vector<pair<double,double>> gammas_by_energy = gammas_by_counts;
    
    std::sort( begin(gammas_by_energy), end(gammas_by_energy), lessThanByEnergy );
    std::sort( begin(gammas_by_counts), end(gammas_by_counts), moreThanByNumPerSecond );
    
    vector<pair<double,double>> answer;
    
    for( const pair<double,double> &energy_counts : gammas_by_counts )
    {
      auto ene_pos = std::lower_bound( begin(gammas_by_energy), end(gammas_by_energy),
                                      energy_counts, lessThanByEnergy );
      if( ene_pos == end(gammas_by_energy) )
        continue;
      
      if( ene_pos->first != energy_counts.first )
        continue; //We've already removed erp from gammas_by_energy
      
      const double energy = energy_counts.first;
      const double counts = energy_counts.second;
      
      const double sigma = fwhm(energy, x) / 2.35482;
      
      double lower = energy - cluster_num_sigma * sigma;
      double upper = energy + cluster_num_sigma * sigma;
      
      // Now remove all gammas from `gammas_by_energy` that are between lower and upper.
      //  Note that if a gamma is equal to the lower or upper bound, it will be removed
      //  (I'm pretty sure)
      const auto start_remove = std::lower_bound(begin(gammas_by_energy), end(gammas_by_energy),
                                                 make_pair(lower,0.0), lessThanByEnergy );
      const auto end_remove = std::upper_bound(begin(gammas_by_energy), end(gammas_by_energy),
                                               make_pair(upper,0.0), lessThanByEnergy );
      
      const double counts_in_region = std::accumulate(start_remove, end_remove, 0.0,
        []( const double &sum, const pair<double,double> &el ){
        return sum + el.second;
      });
      
      const double data_area = m_spectrum->gamma_integral( static_cast<float>(lower),
                                                          static_cast<float>(upper) );
      
#ifndef NDEBUG
      //cout << "For [" << lower << "," << upper << "), will remove :{";
      for( auto iter = start_remove; iter != end_remove; ++iter )
      {
        assert( (iter->first >= lower) && (iter->first <= upper) );
        //cout << iter->first << ",";
      }
      //cout << "}" << endl;
      
      for( const auto &roi : answer )
      {
        assert( (energy <= roi.first) || (energy >= roi.second) );
      }
#endif //#ifndef NDEBUG
      
      gammas_by_energy.erase( start_remove, end_remove );
      
      const double signif = counts_in_region / sqrt(data_area);
      //cout << "For [" << lower << "," << upper << "), there are data_area="
      //<< data_area << ", counts_in_region=" << counts_in_region << ", signif=" << signif
      //<< endl;
      
      if( (data_area > 10) && (counts_in_region > 10) && (signif > 0.01) )
      {
        // Our new ROI could overlap with existing ROI, so we'll fix this up before inserting
        // TODO: right now we are just shrinking the new (smaller) region to not invade the previous region - we should do some sort of weighted balancing
        for( const auto &prev : answer )
        {
          if( (lower < prev.second) && (lower > prev.first) )
            lower = prev.second;
          if( (upper > prev.first) && (upper < prev.second) )
            upper = prev.first;
        }
        
        assert( lower <= upper );
        if( lower <= upper )
         answer.emplace_back( lower, upper );
      }
      
    }//for( const SandiaDecay::EnergyRatePair &erp : gammas_by_counts )
    
    std::sort( begin(answer), end(answer), []( const pair<double,double> &lhs, const pair<double,double> &rhs ){
      return lhs.first < rhs.first;
    });
    
    assert( std::is_sorted( begin(answer), end(answer) ) );
    
    for( size_t i = 1; i < answer.size(); ++i )
    {
      assert( answer[i-1].second <= answer[i].first );
    }
    
    return answer;
  }//vector<pair<double,double>> cluster_photopeaks( const std::vector<double> &x )

  
  template<typename T>
  T fwhm( const T energy, const std::vector<T> &x ) const
  {
    const auto drf_start = x.data() + RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars;
    assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars == m_fwhm_par_start_index );
    const size_t num_drf_par = num_parameters(m_options.fwhm_form);
    
    return eval_fwhm( energy, m_options.fwhm_form, drf_start, num_drf_par );
  }//float fwhm(...)
  
  
  template<typename P, typename T>
  void set_peak_skew( P &peak, const std::vector<T> &x ) const
  {
    const size_t skew_start = RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars
                              + num_parameters(m_options.fwhm_form)
                              + num_total_rel_eff_coefs()
                              + num_total_act_coefs()
                              + 2*m_options.floating_peaks.size();
    assert( skew_start == m_skew_par_start_index );
    
    if( m_options.skew_type == PeakDef::SkewType::NoSkew )
      return;
    
    peak.setSkewType( m_options.skew_type );
    
    const size_t num_skew = PeakDef::num_skew_parameters( m_options.skew_type );
    assert( x.size() >= (skew_start + 2*num_skew) );
    vector<T> skew_pars( begin(x)+skew_start, begin(x)+skew_start+num_skew );
    assert( skew_pars.size() == num_skew );
    
    if( m_skew_has_energy_dependance )
    {
      const double lower_energy = m_spectrum->gamma_channel_lower(0);
      const double upper_energy = m_spectrum->gamma_channel_upper( m_spectrum->num_gamma_channels() - 1 );
      const T mean_frac = (peak.mean() - lower_energy) / (upper_energy - lower_energy);
      
      for( size_t i = 0; i < skew_pars.size(); ++i )
      {
        const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + i);
        if( PeakDef::is_energy_dependent( m_options.skew_type, ct ) )
        {
          const T lower_en_val = x[skew_start + i];
          const T upper_en_val = x[skew_start + i + num_skew];
          assert( upper_en_val > -999.0 );//should NOT have value -999.9
          
          const T val = lower_en_val + mean_frac*(upper_en_val - lower_en_val);
          peak.set_coefficient( val, ct );
        }else
        {
          assert( x[skew_start + num_skew + i] < -999.0 ); //should have value -999.9
          const T val = x[skew_start + i];
          peak.set_coefficient( val, ct );
        }
      }
    }else
    {
      for( size_t i = 0; i < skew_pars.size(); ++i )
      {
        assert( x[skew_start + num_skew + i] < -999.0 ); //should have value -999.9
        const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + i);
        const T val = x[skew_start + i];
        peak.set_coefficient( val, ct );
      }
    }//if( m_skew_has_energy_dependance )
  }//set_peak_skew(...)
  
  /// Returns the index for this nuclide within the relative efficiency curve
  size_t nuclide_index( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index ) const
  {
    assert( nuc );
    if( !nuc )
      throw runtime_error( "Null nuclide" );
    
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." );
    
    const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    size_t nuc_index = rel_eff_curve.nuclides.size();
    for( size_t i = 0; i < nuc_index; ++i )
    {
      if( rel_eff_curve.nuclides[i].nuclide == nuc )
        nuc_index = i;
    }
    
    assert( nuc_index < rel_eff_curve.nuclides.size() );
    if( nuc_index >= rel_eff_curve.nuclides.size() )
      throw std::logic_error( "Invalid nuclide." );
    
    return nuc_index;
  }//size_t nuclide_index( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index )
  

  /** Returns the index for this nuclide's activity (the age is +1) within the Ceres parameter vector */
  size_t nuclide_parameter_index( const SandiaDecay::Nuclide *nuclide, const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." );
    
    size_t index = m_acts_par_start_index;
    for( size_t i = 0; i < rel_eff_index; ++i )
      index += 2*m_options.rel_eff_curves[i].nuclides.size();
  
    return index + 2*nuclide_index(nuclide, rel_eff_index);
  }//size_t nuclide_parameter_index( const SandiaDecay::Nuclide *nuclide, const size_t rel_eff_index )


  size_t rel_act_start_parameter( const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    
    size_t index = m_acts_par_start_index;
    for( size_t i = 0; i < rel_eff_index; ++i )
      index += 2*m_options.rel_eff_curves[i].nuclides.size();
    
    assert( index < number_parameters() );

    return index;
  }//size_t rel_act_start_parameter( const size_t rel_eff_index ) const

  size_t rel_act_num_parameters( const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." );

    return 2*m_options.rel_eff_curves[rel_eff_index].nuclides.size();
  }//size_t rel_act_num_parameters( const size_t rel_eff_index ) const
  
  size_t rel_eff_eqn_start_parameter( const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." );
    
    size_t rel_eff_start_index = RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars
                                       + num_parameters(m_options.fwhm_form);
    for( size_t i = 0; i < rel_eff_index; ++i )
    {
      const auto &rel_eff_curve = m_options.rel_eff_curves[i];
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        rel_eff_start_index += (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2);
      else
        rel_eff_start_index += (rel_eff_curve.rel_eff_eqn_order + 1);
    }
    
    return rel_eff_start_index;
  }//size_t rel_eff_eqn_start_parameter( const size_t rel_eff_index ) const

  
  size_t rel_eff_eqn_num_parameters( const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." );
    
    const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      return (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2);
    
    return (rel_eff_curve.rel_eff_eqn_order + 1);
  }//size_t rel_eff_eqn_num_parameters( const size_t rel_eff_index ) const
  
  
  size_t num_total_rel_eff_coefs() const
  {
    size_t num_coefs = 0;
    for( const auto &rel_eff_curve : m_options.rel_eff_curves )
    {
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        num_coefs += (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2);
      else
        num_coefs += (rel_eff_curve.rel_eff_eqn_order + 1);
    }
    return num_coefs;
  }//size_t num_total_rel_eff_coefs() const
  
  
  size_t num_total_act_coefs() const
  {
    size_t num_coefs = 0;
    for( const auto &nucs : m_nuclides )
      num_coefs += 2*nucs.size();
    return num_coefs;
  }//size_t num_total_act_coefs() const
  
  
  std::string parameter_name( const size_t index ) const
  {
    // We'll try to return 12 character, or less names.
    assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars == m_fwhm_par_start_index );
    
    if( index < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars )
    {
      if( index == 0 )
        return "EneOffset";
      if( index == 1 )
        return "EneGain";
      if( index == 2 )
        return "EneQuad";
      assert( 0 );
      static_assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 3, "Need to update sm_num_energy_cal_pars" );
      return "Unknown";
    }//if( index < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars )
    
    
    if( index < m_rel_eff_par_start_index )
    {
      // Its a FWHM parameter
      return "FWHM_" + std::to_string(index - m_fwhm_par_start_index);
    }//if( index < m_rel_eff_par_start_index )
    
    
    if( index < m_acts_par_start_index )
    {
      size_t coef_num = m_rel_eff_par_start_index;

      const size_t num_rel_eff_curves = m_options.rel_eff_curves.size();
      for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
      {
        const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
        size_t num_coefs_this_rel_eff = rel_eff_curve.rel_eff_eqn_order + 1;
        if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          num_coefs_this_rel_eff += 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2;

        const size_t sub_ind = index - coef_num;
        if( sub_ind >= num_coefs_this_rel_eff )
        {
          coef_num += num_coefs_this_rel_eff;
          continue;
        }

        const string re_ind = num_rel_eff_curves > 1 ? std::to_string(rel_eff_index) : "";
        if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        {
          if( sub_ind == 0 )
            return "SAtt" + re_ind + "(AN)";
          if( sub_ind == 1 )
            return "SAtt" + re_ind + "(AD)";
        
          if( sub_ind < (2 + 2*rel_eff_curve.phys_model_external_atten.size()) )
          {
            const size_t ext_atten_num = (sub_ind - 2) / 2;
            if( (sub_ind % 2) == 0 )
              return "EAtt" + re_ind + std::to_string(ext_atten_num) + "(AN)";
            return "EAtt" + re_ind + std::to_string(ext_atten_num) + "(AD)";
          }
        
          assert( sub_ind >= (2 + 2*rel_eff_curve.phys_model_external_atten.size()) );
          const size_t hoerl_num = sub_ind - 2 - 2*rel_eff_curve.phys_model_external_atten.size();
          if( hoerl_num == 0 )
            return "Hoerl" + re_ind + "(b)";
        
          assert( hoerl_num == 1 );
          if( hoerl_num == 1 )
            return "Hoerl" + re_ind + "(c)";
        }else
        {
          return "RE_" + re_ind + std::to_string(index - sub_ind);
        }//if( physical model ) / else
      }//for( const auto &rel_eff_curve : m_options.rel_eff_curves )

      assert( 0 );
      throw std::logic_error( "Logic for determining Physical Model coefficient name is bad." );
    }//if( index < m_acts_par_start_index )
    
    if( index < m_free_peak_par_start_index )
    {
      // Each index gets two parameters, activity, and age
      size_t coef_num = m_acts_par_start_index;
      const size_t num_rel_eff_curves = m_options.rel_eff_curves.size();
      for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
      {
        const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
        const size_t act_index = index - coef_num;
        const size_t num_coefs_this_rel_eff = 2*rel_eff_curve.nuclides.size();
        if( act_index >= num_coefs_this_rel_eff )
        {
          coef_num += num_coefs_this_rel_eff;
          continue;
        }

        const string re_ind = num_rel_eff_curves > 1 ? std::to_string(rel_eff_index) : "";
        const size_t act_num = act_index / 2;
        assert( act_num < rel_eff_curve.nuclides.size() );
        if( act_num >= rel_eff_curve.nuclides.size() )
          throw runtime_error( "Logic for determining nuclide index is total whack, yo" );
      
        assert( rel_eff_curve.nuclides[act_num].nuclide );
        if( (act_index % 2) == 0 )
          return "Act" + re_ind + "(" + rel_eff_curve.nuclides[act_num].nuclide->symbol + ")";
        return "Age" + re_ind + "(" + rel_eff_curve.nuclides[act_num].nuclide->symbol + ")";
      }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    }//if( index < m_free_peak_par_start_index )
    
    if( index < m_skew_par_start_index )
    {
      return "FPeak_" + std::to_string(index - m_free_peak_par_start_index);
    }//if( index < m_skew_par_start_index )
    
    if( index < m_add_br_uncert_start_index )
    {
      // TODO: we could give a better name for the skew parameter here
      return "Skew_" + std::to_string(index - m_skew_par_start_index);
    }//if( index < m_add_br_uncert_start_index )
    
    assert( m_add_br_uncert_start_index != std::numeric_limits<size_t>::max() );
    assert( (m_add_br_uncert_start_index + m_peak_ranges_with_uncert.size()) == number_parameters() );
    
    if( m_add_br_uncert_start_index == std::numeric_limits<size_t>::max() )
      throw runtime_error( "Whack computing of parameter name - things not legit." );
    
    if( index >= (m_add_br_uncert_start_index + m_peak_ranges_with_uncert.size()) )
      throw runtime_error( "Bad computing of parameter name - too large of index" );
    
    const size_t range_ind = index - m_add_br_uncert_start_index;
    const int mid_energy = static_cast<int>( std::round(0.5*(m_peak_ranges_with_uncert[range_ind].second
                                      + m_peak_ranges_with_uncert[range_ind].second)) );
    
    return "dBr" + std::to_string(mid_energy);
  }//std::string parameter_name( const size_t index ) const
  
  
  template<typename T>
  T relative_activity( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( x.size() == number_parameters() );
    assert( rel_eff_index < m_nuclides.size() );
    if( rel_eff_index >= m_nuclides.size() )
      throw std::logic_error( "relative_activity: invalid relative efficiency curve index." );
    
    const size_t parent_nuc_par_index = nuclide_parameter_index( nuc, rel_eff_index );
    const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];

    const size_t nuc_index = nuclide_index( nuc, rel_eff_index );
    const NucInputGamma &nuc_info = nuclides[nuc_index];
    
    assert( nuc_info.activity_multiple > 0.0 );

    return nuc_info.activity_multiple * x[parent_nuc_par_index];
  }//double relative_activity(...)
  
  
  /** Returns the nuclide that is responsible for setting the passed in nuclides age.
   Will return the input nuclide if it controls its own age.
   */
  const SandiaDecay::Nuclide *age_controlling_nuc( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." ); 
    
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    if( m_nuclides.size() != m_options.rel_eff_curves.size() )
      throw std::logic_error( "Invalid number of nuclides." );

    const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];
    const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    if( !rel_eff_curve.nucs_of_el_same_age )
      return nuc;
    
    const size_t nuc_index = nuclide_index( nuc, rel_eff_index );
    assert( nuc_index < nuclides.size() );
    if( nuc_index >= nuclides.size() )
      throw std::logic_error( "Invalid nuclides index." );
    
    for( size_t i = 0; i < nuc_index; ++i )
    {
      if( nuclides[i].nuclide->atomicNumber == nuclides[nuc_index].nuclide->atomicNumber )
        return nuclides[i].nuclide;
    }
    
    return nuc;
  }//age_controlling_nuc(...)
  
  
  template<typename T>
  T age( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( rel_eff_index < m_nuclides.size() );
    if( rel_eff_index >= m_nuclides.size() )
      throw std::logic_error( "age: invalid releff passed in." );
    
    const size_t act_start_index = RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars
                                   + num_parameters(m_options.fwhm_form)
                                   + num_total_rel_eff_coefs();
    assert( act_start_index == m_acts_par_start_index );
    
    // The `parent_nuc` will often times be `nuc`.
    const SandiaDecay::Nuclide *parent_nuc = age_controlling_nuc( nuc, rel_eff_index );
    const size_t parent_par_index = nuclide_parameter_index( parent_nuc, rel_eff_index );
    const size_t par_nuc_index = nuclide_index( parent_nuc, rel_eff_index );

    const T age_scaler = x[parent_par_index + 1];
    const T age = age_scaler * m_nuclides[rel_eff_index][par_nuc_index].age_multiple;

    // We'll allow a little bit of slack on the age...
    assert( age >= static_cast<double>(-std::numeric_limits<float>::epsilon()) );
    if( age < static_cast<double>(-std::numeric_limits<float>::epsilon()) )
      throw runtime_error( "Negative age for " + m_nuclides[rel_eff_index][par_nuc_index].nuclide->symbol + " found." );
    
    return (max)( age, T(0.0) );
  }//age( nuclide )
  
  
  bool is_fixed_age( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index ) const
  {
    const SandiaDecay::Nuclide *parent_nuc = age_controlling_nuc( nuc, rel_eff_index );
    const size_t nuc_index = nuclide_index( parent_nuc, rel_eff_index );
    return !m_nuclides[rel_eff_index][nuc_index].fit_age;
  }//
  
  
  
  template<typename T>
  struct PhysModelRelEqnDef
  {
    shared_ptr<const DetectorPeakResponse> det;
    
    std::optional<RelActCalc::PhysModelShield<T>> self_atten;
    std::vector<RelActCalc::PhysModelShield<T>> external_attens;
    
    std::optional<T> hoerl_b, hoerl_c;
  };//struct PhysModelRelEqnDef
  
  
  template<typename T>
  static PhysModelRelEqnDef<T> make_phys_eqn_input( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve,
                                                std::shared_ptr<const DetectorPeakResponse> drf,
                                                const std::vector<T> &coeffs,
                                                const size_t rel_eff_start )
  {
    assert( (rel_eff_start + 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2) <= coeffs.size() );
    if( (rel_eff_start + 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2) > coeffs.size() )
      throw std::logic_error( "make_phys_eqn_input: whack number of inputs." );
    
    PhysModelRelEqnDef<T> answer;
    answer.det = drf;
    
    if( rel_eff_curve.phys_model_self_atten )
    {
      RelActCalc::PhysModelShield<T> atten;
      atten.material = rel_eff_curve.phys_model_self_atten->material;
      if( !atten.material )
      {
        atten.atomic_number = coeffs[rel_eff_start + 0] * RelActCalc::ns_an_ceres_mult;
        assert( (atten.atomic_number >= 1.0) && (atten.atomic_number <= 98.0) );
        if( (atten.atomic_number < 1.0) || (atten.atomic_number > 98.0) )
          throw std::logic_error( "make_phys_eqn_input: whack self-atten AN" );
      }
      atten.areal_density = coeffs[rel_eff_start + 1] * PhysicalUnits::g_per_cm2;
      
      answer.self_atten = std::move(atten);
    }//if( rel_eff_curve.phys_model_self_atten )
    
    for( size_t ext_ind = 0; ext_ind < rel_eff_curve.phys_model_external_atten.size(); ++ext_ind )
    {
      assert( rel_eff_curve.phys_model_external_atten[ext_ind] );
      const RelActCalc::PhysicalModelShieldInput &ext_opt = *rel_eff_curve.phys_model_external_atten[ext_ind];
      
      RelActCalc::PhysModelShield<T> atten;
      atten.material = ext_opt.material;
      if( !atten.material )
      {
        atten.atomic_number = coeffs[rel_eff_start + 2 + 2*ext_ind + 0] * RelActCalc::ns_an_ceres_mult;
        assert( (atten.atomic_number >= 1.0) && (atten.atomic_number <= 98.0) );
        if( (atten.atomic_number < 1.0) || (atten.atomic_number > 98.0) )
          throw std::logic_error( "make_phys_eqn_input: whack external AN" );
      }
      atten.areal_density = coeffs[rel_eff_start + 2 + 2*ext_ind + 1] * PhysicalUnits::g_per_cm2;
     
      answer.external_attens.push_back( std::move(atten) );
    }//for( loop over external attenuators )
    
    if( rel_eff_curve.phys_model_use_hoerl )
    {
      const size_t b_index = rel_eff_start + 2 + 2*rel_eff_curve.phys_model_external_atten.size();
      const size_t c_index = b_index + 1;
      
      //TODO: maybe should consult: phys_model_use_hoerl?
      const T b = coeffs[b_index]; //(energy/1000)^b
      const T c = coeffs[c_index]; //c^(1000/energy)
      
      answer.hoerl_b = coeffs[b_index];
      answer.hoerl_c = coeffs[c_index];
    }//if( (b != 0.0) || (c != 1.0) )
    
    return answer;
  }//PhysModelRelEqnDef make_phys_eqn_input(...)
  
  
  template<typename T>
  T relative_eff( const double energy, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "relative_eff(): invalid rel eff curve index" );
    
    const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    const size_t rel_eff_start_index = rel_eff_eqn_start_parameter(rel_eff_index);
    const size_t num_rel_eff_par = rel_eff_eqn_num_parameters(rel_eff_index);
    
    assert( (rel_eff_start_index + num_rel_eff_par) < x.size() );
    assert( (rel_eff_index != 0) || (rel_eff_start_index == m_rel_eff_par_start_index) );
    
    const T * const rel_ef_pars = &(x[rel_eff_start_index]);
    
    if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      assert( (rel_eff_start_index + num_rel_eff_par) < x.size() );
      assert( m_drf );
      
      const PhysModelRelEqnDef<T> re_input = make_phys_eqn_input( rel_eff_curve,
                                                             m_drf, x, rel_eff_start_index );
      
      return RelActCalc::eval_physical_model_eqn_imp( energy, re_input.self_atten,
                                                 re_input.external_attens, re_input.det.get(),
                                                 re_input.hoerl_b, re_input.hoerl_c );
    }//if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    
    assert( rel_eff_curve.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel );
    
    return RelActCalc::eval_eqn_imp( energy, rel_eff_curve.rel_eff_eqn_type, rel_ef_pars, num_rel_eff_par );
  }//
  
  
  /** Translates from a "true" energy (e.g., that of a gamma, or ROI bounds), to the energy
   of m_energy_cal.
   */
  template<typename T>
  T apply_energy_cal_adjustment( double energy, const std::vector<T> &x ) const
  {
    // TODO: Auto differentiation doesnt _seem_ to capture effects of energy calibration
    //assert( std::is_same_v<T, double> );
    
    assert( x.size() > RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars );
    assert( 0 == m_energy_cal_par_start_index );
    
    if( !m_options.fit_energy_cal )
    {
      if constexpr ( !std::is_same_v<T, double> )
      {
        assert( abs(x[0] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) < 1.0E-5 );
        assert( abs(x[1] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) < 1.0E-5 );
      }
      
      return T(energy);
    }//if( we arent fitting energy cal )
    
    if( (x[0] == RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset)
       && (x[1] == RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) )
      return T(energy);
    
    // Check adjustments are near the limits we placed (which was [-5,5], and [0.985,1.015])
    assert( abs(x[0] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) <= RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset );
    assert( abs(x[1] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) <= RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset );
    assert( (RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 2)
           || (abs(x[2] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) <= RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) );
    
    //if( x[0] != RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset || x[1] != RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset
    //   || x[2] != RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset )
    //cout << "offset_delta=" << (x[0]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0) << ", "
    // << "gain_delta=" << (x[1]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0) << ", "
    //<< "cubic_delta=" << (x[2]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0) << endl;
    
    const T offest_adj = (x[0]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                            * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
    const T gain_adj = (x[1]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                            * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
    const T quad_adj = (RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars > 2)
                              ? ((x[2]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                                 * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple)
                              : T(0.0);
    
    
    
    // TODO: implement and try out optimization so that if there is no adjustment to be made, skip doing the energy cal adjustment.
    //       Note: we do need to be careful we dont make the cuts so large that the steps of the
    //       numerical differentiation around zero will fail (I did have trouble with this for
    //       relative activities and using FLT_EPSILON as a cutoff).
    //if( (fabs(x[0]) < 1.0E-9 ) && (fabs(x[1]) < 1.0E-12) )
    //  return energy;
    
    switch( m_energy_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      {
        const double num_channel = static_cast<double>( m_energy_cal->num_channels() );
        const vector<float> &orig_coeffs = m_energy_cal->coefficients();
        vector<T> coefs( orig_coeffs.size() );
        for( size_t i = 0; i < orig_coeffs.size(); ++i )
          coefs[i] = T( orig_coeffs[i] );
        
        assert( coefs.size() >= 2 );
        coefs[0] += offest_adj;
        coefs[1] += (gain_adj / num_channel);
        if( quad_adj != 0.0 )
        {
          if( coefs.size() > 2 )
            coefs[2] += (quad_adj / (num_channel*num_channel));
          else
            coefs.push_back( quad_adj / (num_channel*num_channel) );
        }//if( quad_adj != 0.0 )
        
        const auto &dev_pairs = m_energy_cal->deviation_pairs();
        const double channel = m_energy_cal->channel_for_energy( energy );
        
        //return SpecUtils::polynomial_energy( channel, coefs, dev_pairs );
        
        T val = coefs[0];
        for( size_t i = 1; i < coefs.size(); ++i )
          val += coefs[i] * pow( channel, static_cast<double>(i) );
        
        if( !dev_pairs.empty() )
        {
          if constexpr ( !std::is_same_v<T, double> )
            val += T(SpecUtils::deviation_pair_correction( static_cast<float>(val.a), dev_pairs ));
          else
            val += SpecUtils::deviation_pair_correction( static_cast<float>(val), dev_pairs );
        }//if( !dev_pairs.empty() )
        
        return val;
      }//case polynomial or FRF
      
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
      {
        const double num_channel = static_cast<double>( m_energy_cal->num_channels() );
        const vector<float> &orig_coeffs = m_energy_cal->coefficients();
        const size_t ncoeffs = std::min( orig_coeffs.size(), size_t(4) );
        assert( ncoeffs >= 2 );
        vector<T> coefs( ncoeffs );
        for( size_t i = 0; i < ncoeffs; ++i )
          coefs[i] = T( orig_coeffs[i] );
        
        coefs[0] += offest_adj;
        coefs[1] += gain_adj;
        
        if( quad_adj != 0.0 )
        {
          if( coefs.size() > 2 )
            coefs[2] += quad_adj;
          else
            coefs.push_back( quad_adj );
        }//if( quad_adj != 0.0 )
        
        const auto &dev_pairs = m_energy_cal->deviation_pairs();
        const double channel = m_energy_cal->channel_for_energy( energy );
        
        //return SpecUtils::fullrangefraction_energy( channel, coefs, num_channel, dev_pairs );
        
        const double frac_chan = channel / num_channel;
          

        T val = coefs[0];
        for( size_t c = 1; c < ncoeffs; ++c )
          val += coefs[c] * pow( frac_chan, static_cast<double>(c) );
        
        if( coefs.size() > 4 )
          val += coefs[4] / (1.0 + 60.0*frac_chan);
          
        if( !dev_pairs.empty() )
        {
          if constexpr ( !std::is_same_v<T, double> )
            val += T(SpecUtils::deviation_pair_correction( static_cast<float>(val.a), dev_pairs ));
          else
            val += SpecUtils::deviation_pair_correction( static_cast<float>(val), dev_pairs );
        }//if( !dev_pairs.empty() )
        
        return val;
      }//case polynomial or FRF
      
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      {
        const double lower_energy = m_energy_cal->lower_energy();
        const double range = m_energy_cal->upper_energy() - lower_energy;
        const double range_frac = (energy - lower_energy) / range;
        
        T new_energy = energy + offest_adj + (range_frac*gain_adj);
        assert( quad_adj == 0.0 );
        //if( quad_adj != 0.0 )
        //  new_energy += quad_adj*range_frac*range_frac;
        
        return new_energy;
      }//case LowerChannelEdge:
        
      case SpecUtils::EnergyCalType::InvalidEquationType:
        break;
    }//switch( m_energy_cal->type() )
    
    assert( 0 );
    throw runtime_error( "Energy cal must be valid" );
    
    return T(energy);
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
      assert( fabs(1.0 - x[1]) < std::numeric_limits<float>::epsilon() );
      
      return adjusted_energy;
    }//if( we arent fitting energy cal )
    
    if( (x[0] == 0.0) && (x[1] == 1.0) )
      return adjusted_energy;
    
    // Check adjustments are near the limits we placed (which was [-5,5], and [0.985,1.015])
    assert( fabs(x[0]) <= 5.1 );
    assert( (x[1] >= 0.984) && (x[1] <= 1.016) );
    
    const double offset_adj = (x[0]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                              * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
    const double gain_adj = (x[1]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                              * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
    const double quad_adj = (RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars > 2)
                      ? ((x[2]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                         * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple)
                      : 0.0;
    
    switch( m_energy_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      {
        const size_t num_channel = m_energy_cal->num_channels();
        vector<float> coefs = m_energy_cal->coefficients();
        assert( coefs.size() >= 2 );
        coefs[0] += offset_adj;
        coefs[1] += (gain_adj / num_channel);
        
        if( quad_adj != 0.0 )
        {
          if( coefs.size() > 2 )
            coefs[2] += (quad_adj / (num_channel*num_channel));
          else
            coefs.push_back( quad_adj / (num_channel*num_channel) );
        }//if( quad_adj != 0.0 )
        
        const auto &dev_pairs = m_energy_cal->deviation_pairs();
        
        // TODO: is there a cheaper way to do this?
        const double channel = SpecUtils::find_polynomial_channel( adjusted_energy, coefs, num_channel, dev_pairs );
        
        return m_energy_cal->energy_for_channel( channel );
      }//case polynomial
        
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
      {
        const size_t num_channel = m_energy_cal->num_channels();
        vector<float> coefs = m_energy_cal->coefficients();
        assert( coefs.size() >= 2 );
        coefs[0] += offset_adj;
        coefs[1] += gain_adj;
        
        if( quad_adj != 0.0 )
        {
          if( coefs.size() > 2 )
            coefs[2] += quad_adj;
          else
            coefs.push_back( quad_adj );
        }//if( quad_adj != 0.0 )
        
        const auto &dev_pairs = m_energy_cal->deviation_pairs();
        
        // TODO: is there a cheaper way to do this?
        const double channel = SpecUtils::find_fullrangefraction_channel( adjusted_energy, coefs, num_channel, dev_pairs );
        
        return m_energy_cal->energy_for_channel( channel );
      }//case FRF
        
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      {
        const double lower_energy = m_energy_cal->lower_energy();
        const double range = m_energy_cal->upper_energy() - lower_energy;
        
        assert( quad_adj == 0.0 );
        return (adjusted_energy - offset_adj + gain_adj*lower_energy/range)/(1 + gain_adj/range);
      }//case LowerChannelEdge:
        
      case SpecUtils::EnergyCalType::InvalidEquationType:
        break;
    }//switch( m_energy_cal->type() )
    
    assert( 0 );
    throw runtime_error( "Energy cal must be valid" );
    
    return adjusted_energy;
  }//double un_apply_energy_cal_adjustment( double energy, const std::vector<double> &x ) const
  
  
  /** A stand-in for the `PeakDef` class to allow auto-differentiation, and also simplify things  */
  template<typename T>
  struct PeakDefImp
  {
    T m_mean = T(0.0);
    T m_sigma = T(0.0);
    T m_amplitude = T(0.0);
    T m_skew_pars[4] = { T(0.0), T(0.0), T(0.0), T(0.0) };
    
    const SandiaDecay::Nuclide *m_parent_nuclide = nullptr;
    const SandiaDecay::Transition *m_transition = nullptr;
    size_t m_rad_particle_index = 0;
    
    //We will probably eventually implement x-rays and reactions...
    //const SandiaDecay::Element *m_xray_element = nullptr;
    //double m_xray_energy = 0.0l
    //const ReactionGamma::Reaction *m_reaction = nullptr;
    //double m_reaction_energy = 0.0l
    
    PeakDef::SkewType m_skew_type = PeakDef::SkewType::NoSkew;
    PeakDef::SourceGammaType m_gamma_type = PeakDef::SourceGammaType::NormalGamma;
    
    size_t m_rel_eff_index = std::numeric_limits<size_t>::max();
    
    const T &mean() const { return m_mean; }
    
    inline void setSkewType( const PeakDef::SkewType &type )
    {
      m_skew_type = type;
    }
    
    inline void set_coefficient( T val, const PeakDef::CoefficientType &coef )
    {
      const int index = static_cast<int>( coef - PeakDef::CoefficientType::SkewPar0 );
      assert( (index >= 0) && (index < 4) );
      m_skew_pars[index] = val;
    }
    
    void gauss_integral( const float *energies, T *channels, const size_t nchannel ) const
    {
      switch( m_skew_type )
      {
        case PeakDef::SkewType::NoSkew:
          PeakDists::gaussian_integral( m_mean, m_sigma, m_amplitude, energies, channels, nchannel );
          break;
          
        case PeakDef::SkewType::Bortel:
          PeakDists::bortel_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], energies, channels, nchannel );
          break;
          
        case PeakDef::SkewType::CrystalBall:
          PeakDists::crystal_ball_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
          break;
          
        case PeakDef::SkewType::DoubleSidedCrystalBall:
          PeakDists::double_sided_crystal_ball_integral( m_mean, m_sigma, m_amplitude,
                                             m_skew_pars[0], m_skew_pars[1],
                                             m_skew_pars[2], m_skew_pars[3],
                                             energies, channels, nchannel );
          break;
          
        case PeakDef::SkewType::GaussExp:
          PeakDists::gauss_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], energies, channels, nchannel );
          break;
          
        case PeakDef::SkewType::ExpGaussExp:
          PeakDists::exp_gauss_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
          break;
      }//switch( skew_type )
    }//void gauss_integral( const float *energies, double *channels, const size_t nchannel ) const
  };//struct PeakDefImp
  
  template<typename T>
  struct PeakContinuumImp
  {
    PeakContinuum::OffsetType m_type = PeakContinuum::OffsetType::NoOffset;
    
    T m_lower_energy = T(0.0);
    T m_upper_energy = T(0.0);
    T m_reference_energy = T(0.0);
    std::array<T,4> m_values = { T(0.0), T(0.0), T(0.0), T(0.0) };
  };//struct PeakContinuumImp
  
  
  template<typename T>
  struct PeaksForEnergyRangeImp
  {
    std::vector<PeakDefImp<T>> peaks;
    
    PeakContinuumImp<T> continuum;
    
    size_t first_channel;
    size_t last_channel;
    bool no_gammas_in_range;
    bool forced_full_range;
    
    /** Peak plus continuum counts for [first_channel, last_channel] */
    std::vector<T> peak_counts;
  };//struct PeaksForEnergyRangeImp
  
  template<typename T>
  PeaksForEnergyRangeImp<T> peaks_for_energy_range_imp( const RoiRangeChannels &range, const std::vector<T> &x ) const
  {
    const size_t num_channels = range.num_channels;
    
    // We will use "adjusted" to refer to energies that have been mapped into the spectrums original
    //  energy calibrations
    const T adjusted_lower_energy = apply_energy_cal_adjustment( range.lower_energy, x );
    const T adjusted_upper_energy = apply_energy_cal_adjustment( range.upper_energy, x );
    
    // TODO: Check this conversion from `Jet<>` to double doesnt mess anything up - I *think* this is _fine_...
    pair<size_t,size_t> channel_range;
    if constexpr ( !std::is_same_v<T, double> )
      channel_range = range.channel_range( adjusted_lower_energy.a, adjusted_upper_energy.a, num_channels, m_energy_cal );
    else
      channel_range = range.channel_range( adjusted_lower_energy, adjusted_upper_energy, num_channels, m_energy_cal );
    
    const size_t first_channel = channel_range.first;
    const size_t last_channel = channel_range.second;
    
    //const float adjusted_first_channel_lower = m_energy_cal->energy_for_channel(first_channel);
    //const float adjusted_last_channel_upper = m_energy_cal->energy_for_channel(last_channel + 0.999 );
    
    const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters( range.continuum_type ) );
    const bool is_step_continuum = PeakContinuum::is_step_continuum( range.continuum_type );
    
    size_t num_free_peak_pars = 0;
    set<const SandiaDecay::Nuclide *> nuclides_used;
    
    PeaksForEnergyRangeImp<T> answer;
    answer.first_channel = first_channel;
    answer.last_channel = last_channel;
    answer.no_gammas_in_range = false;
    answer.forced_full_range = range.force_full_range;
    
    vector<PeakDefImp<T>> &peaks = answer.peaks;
    
    // Go through and create peaks based on rel act, eff, etc
    for( size_t rel_eff_index = 0; rel_eff_index < m_nuclides.size(); ++rel_eff_index )
    {
      const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];
      for( const NucInputGamma &nucinfo : nuclides )
      {
        const vector<NucInputGamma::EnergyYield> *gammas = nullptr;
        std::unique_ptr<vector<NucInputGamma::EnergyYield>> aged_gammas_cache;
        
        const T rel_act = relative_activity( nucinfo.nuclide, rel_eff_index, x );
        //cout << "peaks_for_energy_range_imp: Relative activity of " << nucinfo.nuclide->symbol
        //     << " is " << PhysicalUnits::printToBestActivityUnits(rel_act) << endl;
        
        if( is_fixed_age(nucinfo.nuclide, rel_eff_index) )
        {
          gammas = &(nucinfo.nominal_gammas);
        }else
        {
          const T nuc_age = age(nucinfo.nuclide, rel_eff_index, x);
          
          double nuc_age_val;
          if constexpr ( !std::is_same_v<T, double> )
          {
            cerr << "Fitting age is not implemented for using Jets" << endl;
            assert(0);
            throw logic_error( "Fitting age is not implemented for using Jets" );
            
            nuc_age_val = nuc_age.a;
          }else
          {
            nuc_age_val = nuc_age;
          }
          
          aged_gammas_cache.reset( new vector<NucInputGamma::EnergyYield>() );
          *aged_gammas_cache = NucInputGamma::decay_gammas( nucinfo.nuclide, nuc_age_val, nucinfo.gammas_to_exclude );
          gammas = aged_gammas_cache.get();
        }//if( age is fixed ) / else( age may vary )
        
        assert( gammas );
        
        for( const NucInputGamma::EnergyYield &gamma : *gammas )
        {
          const double energy = gamma.energy;
          const double yield = gamma.yield;
          const size_t transition_index = gamma.transition_index;
          const SandiaDecay::Transition * const transition = gamma.transition;
          const PeakDef::SourceGammaType gamma_type = gamma.gamma_type;
          
          assert( transition || (gamma_type == PeakDef::SourceGammaType::AnnihilationGamma) );
          assert( !transition
                 || (gamma_type == PeakDef::SourceGammaType::AnnihilationGamma)
                 || (fabs(transition->products[transition_index].energy - energy) < 0.0001) );
          
          // Filter out zero-amplitude gammas
          // TODO: - come up with more intelligent lower bound of gamma rate to bother wit
          if( yield < std::numeric_limits<float>::min() )
            continue;
          
          // Filter the energy range based on true energies (i.e., not adjusted energies)
          if( (gamma.energy < range.lower_energy) || (gamma.energy > range.upper_energy) )
            continue;
          
          nuclides_used.insert( nucinfo.nuclide );
          
          // We compute the relative efficiency and FWHM based off of "true" energy
          const T rel_eff = relative_eff( gamma.energy, rel_eff_index, x );
          if( isinf(rel_eff) || isnan(rel_eff) )
            throw runtime_error( "peaks_for_energy_range_imp: inf or NaN rel. eff for "
                                + std::to_string(gamma.energy) + " keV."  );
          
          const T peak_fwhm = fwhm( T(gamma.energy), x );
          double fwhm;
          if constexpr ( !std::is_same_v<T, double> )
            fwhm = peak_fwhm.a;
          else
            fwhm = peak_fwhm;
          
          if( isinf(peak_fwhm) || isnan(peak_fwhm) )
          {
            stringstream msg;
            msg << "peaks_for_energy_range_imp: " << fwhm << " FWHM for "
            << std::setprecision(2) << gamma.energy << " keV, from pars={";
            const size_t num_drf_par = num_parameters(m_options.fwhm_form);
            for( size_t i = 0; i < num_drf_par; ++i )
            {
              double par_val;
              if constexpr ( !std::is_same_v<T, double> )
                par_val = x[2 + i].a;
              else
                par_val = x[2 + i];
              
              msg << (i ? ", " : "") << par_val;
            }
            msg << "}";
            
            throw runtime_error( msg.str() );
          }//if( IsInf(peak_fwhm) || IsNan(peak_fwhm) )
          
          // Do a sanity check to make sure peak isnt getting too narrow
          const double nchannel = m_energy_cal->channel_for_energy(gamma.energy + 0.5*fwhm)
                                           - m_energy_cal->channel_for_energy(gamma.energy - 0.5*fwhm);
          if( nchannel < 1.5 )
            throw runtime_error( "peaks_for_energy_range_imp: for peak at " + std::to_string(gamma.energy)
                                + " keV, FWHM=" + std::to_string(fwhm) + " which is only "
                                + std::to_string(nchannel) + "channels - too small." );
          
          if( peak_fwhm < 0.001 )
            throw runtime_error( "peaks_for_energy_range_imp: for peak at " + std::to_string(gamma.energy)
                                + " keV, FWHM=" + std::to_string(fwhm) + " which is too small." );
          
          T br_uncert_adj(1.0);
          if( (m_options.additional_br_uncert > 0.0) && !m_peak_ranges_with_uncert.empty() )
          {
            assert( m_add_br_uncert_start_index != std::numeric_limits<size_t>::max() );
            assert( (m_add_br_uncert_start_index + m_peak_ranges_with_uncert.size()) == number_parameters() );
            
            // TODO: we could use std::lower_bound to find the applicable `m_peak_ranges_with_uncert` element
            for( size_t range_index = 0; range_index < m_peak_ranges_with_uncert.size(); ++range_index )
            {
              const pair<double,double> &range = m_peak_ranges_with_uncert[range_index];
              
              if( (gamma.energy >= range.first) && (gamma.energy <= range.second) )
              {
                // If `x[m_add_br_uncert_start_index + range_index]` is zero, then our amplitude is zero; if it is two, then our amplitude is doubled.
                const T adjvalue = x[m_add_br_uncert_start_index + range_index]/sm_peak_range_uncert_par_offset - 1.0;
                br_uncert_adj = max(T(0.0), 1.0 + adjvalue);
                
                break;
              }//
            }//for( find range this gamma belongs to, if any )
          }//if( m_options.additional_br_uncert > 0.0 && !m_peak_ranges_with_uncert.empty() )
          
          
          const T peak_amplitude = rel_act * static_cast<double>(m_live_time) * rel_eff * yield * br_uncert_adj;
          if( isinf(peak_amplitude) || isnan(peak_amplitude) )
            throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak amplitude for "
                                + std::to_string(gamma.energy) + " keV.");
          
          const T peak_mean = apply_energy_cal_adjustment( gamma.energy, x );
          if( isinf(peak_mean) || isnan(peak_mean) )
            throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak mean for "
                                + std::to_string(gamma.energy) + " keV.");
          
          if( peak_amplitude < static_cast<double>( std::numeric_limits<float>::min() ) )
          {
            //cout << "peaks_for_energy_range_imp: Peak at " << gamma.energy << " keV for " << nucinfo.nuclide->symbol
            //<< " has a mean of " << peak_mean << " keV, FWHM=" << peak_fwhm << ", RelEff=" << rel_eff
            //<< ", and AMP=" << peak_amplitude << ", rel_act=" << rel_act << ", yield="
            //<< yield << ", m_live_time=" << m_live_time << endl;
            continue;
          }
          
          PeakDefImp<T> peak;
          peak.m_mean = peak_mean;
          peak.m_sigma = peak_fwhm/2.35482;
          peak.m_amplitude = peak_amplitude;
          if( transition || (gamma_type == PeakDef::SourceGammaType::AnnihilationGamma) )
          {
            peak.m_parent_nuclide = nucinfo.nuclide;
            peak.m_transition = transition;
            peak.m_rad_particle_index = transition_index;
            peak.m_gamma_type = gamma_type;
          }
          
          set_peak_skew( peak, x );
          peak.m_rel_eff_index = rel_eff_index;
          
          
          peaks.push_back( std::move(peak) );
        }//for( const SandiaDecay::EnergyRatePair &gamma : gammas )
      }//for( const NucInputGamma &nucinfo : m_nuclides )
    }//for( size_t rel_eff_index = 0; rel_eff_index < m_nuclides.size(); ++rel_eff_index )
    
    const size_t free_peaks_start_index = m_free_peak_par_start_index;
    
    for( size_t index = 0; index < m_options.floating_peaks.size(); ++index )
    {
      const RelActCalcAuto::FloatingPeak &peak = m_options.floating_peaks[index];
      
      if( (peak.energy < range.lower_energy) || (peak.energy > range.upper_energy) )
        continue;
      
      const size_t amp_index = free_peaks_start_index + 2*index + 0;
      const size_t fwhm_index = free_peaks_start_index + 2*index + 1;
      assert( fwhm_index < x.size() );
      
      num_free_peak_pars += 1;
      const T peak_amp = x[amp_index];
      if( isinf(peak_amp) || isnan(peak_amp) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak amplitude for "
                            + std::to_string(peak.energy) + " keV extra peak.");
      
      T peak_mean( peak.energy );
      if( peak.apply_energy_cal_correction )
        peak_mean = apply_energy_cal_adjustment( peak.energy, x );
      
      if( isinf(peak_mean) || isnan(peak_mean) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak mean for "
                            + std::to_string(peak.energy) + " keV extra peak.");
      
      const T peak_fwhm = x[fwhm_index] * fwhm(peak_mean, x);
      
      if( isinf(peak_fwhm) || isnan(peak_fwhm) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak FWHM for "
                            + std::to_string(peak.energy) + " keV extra peak.");
      
      assert( peak.release_fwhm || (abs(x[fwhm_index] - 1.0) < 1.0E-5) );
      
      num_free_peak_pars += peak.release_fwhm;
      
      //cout << "peaks_for_energy_range_imp: free peak at " << peak.energy << " has a FWHM=" << peak_fwhm << " and AMP=" << peak_amp << endl;
      PeakDefImp<T> imp_peak;
      imp_peak.m_mean = peak_mean;
      imp_peak.m_sigma = peak_fwhm/2.35482;
      imp_peak.m_amplitude = peak_amp;
      
      set_peak_skew( imp_peak, x );
      
      peaks.push_back( std::move(imp_peak) );
    }//for( const RelActCalcAuto::FloatingPeak &peak : m_options.floating_peaks )

    if( peaks.empty() )
    {
      // We will add a zero-amplitude peak, and fit the continuum, so this way we account for this
      //  region, even if age or something has drove all the gammas out of this regions
      cerr << "peaks_for_energy_range_imp: no peaks in range [" << range.lower_energy << ", "
           << range.upper_energy << "] keV." << endl;
      
      answer.no_gammas_in_range = true;
      
      const T middle_energy = 0.5*(adjusted_lower_energy + adjusted_upper_energy);
      T middle_fwhm = fwhm( middle_energy, x );
      if( isinf(middle_fwhm) || isnan(middle_fwhm) )
        middle_fwhm = T(1.0); //arbitrary
      
      
      PeakDefImp<T> peak;
      peak.m_mean = T(middle_energy);
      peak.m_sigma = T(middle_fwhm/2.35482);
      peak.m_amplitude = T(0.0);
      
      peaks.push_back( std::move(peak) );
    }//if( peaks.empty() )
    
    
    answer.continuum.m_type = range.continuum_type;
    answer.continuum.m_lower_energy = adjusted_lower_energy;
    answer.continuum.m_upper_energy = adjusted_upper_energy;
    answer.continuum.m_reference_energy = adjusted_lower_energy;
    
    assert( !peaks.empty() );
    assert( num_channels == ((1 + last_channel) - first_channel) );
    assert( first_channel < m_channel_counts.size() );
    assert( last_channel < m_channel_counts.size() );
    assert( m_energy_cal && m_energy_cal->channel_energies() );
    assert( m_energy_cal->channel_energies()->size() >= m_channel_counts.size() );
    
    assert( static_cast<size_t>(num_polynomial_terms) < answer.continuum.m_values.size() );
    
    const T ref_energy = answer.continuum.m_reference_energy;
    const float * const data = &(m_channel_counts[first_channel]);
    const float * const energies = &((*m_energy_cal->channel_energies())[first_channel]);
    
    T *continuum_coeffs = answer.continuum.m_values.data();
    
    std::vector<T> &peak_counts = answer.peak_counts;
    peak_counts.resize( num_channels );
    vector<PeakDefImp<T>> dummy;
    RelActCalcAuto::fit_continuum( energies, data, num_channels, num_polynomial_terms,
                                  is_step_continuum, ref_energy, peaks,
                                  continuum_coeffs,
                                  peak_counts.data() );
    
    /*
#ifndef NDEBUG
     // We use different linear algebra packages and methodologies to fit the coefficients
     //  in the below two functions, so we wont check they are the exact same, just that they
     //  are close
    {//begin test against fit_amp_and_offset
      vector<PeakDef> fixedAmpPeaks;
      for( size_t i = 0; i < peaks.size(); ++i )
      {
        PeakDef p( peaks[i].m_mean, peaks[i].m_sigma, peaks[i].m_amplitude );
        p.setSkewType( peaks[i].m_skew_type );
        const size_t num_skew_par = PeakDef::num_skew_parameters(peaks[i].m_skew_type);
        for( size_t i = 0; i < num_skew_par; ++i )
          p.set_coefficient( peaks[i].m_skew_pars[i], static_cast<PeakDef::CoefficientType>(PeakDef::CoefficientType::SkewPar0 + i) );
        fixedAmpPeaks.push_back( p );
      }
      
      vector<double> dummy_amps, continuum_coeffs_old, dummy_amp_uncert, continuum_uncerts;
      fit_amp_and_offset( energies, data, num_channels, num_polynomial_terms,
                         is_step_continuum, ref_energy, {}, {}, fixedAmpPeaks,
                                             PeakDef::SkewType::NoSkew, nullptr, dummy_amps,
                         continuum_coeffs_old, dummy_amp_uncert, continuum_uncerts );
      
      assert( static_cast<int>(continuum_coeffs_old.size()) == num_polynomial_terms );
      for( int i = 0; i < num_polynomial_terms; ++i )
      {
        const double new_val = continuum_coeffs[i];
        const double old_val = continuum_coeffs_old[i];
        const double max_val = max( fabs(new_val), fabs(old_val) );
        const double diff = fabs( new_val - old_val );
        if( (diff > 1.0E-12) && (diff > 1.0E-1*max_val) )
          cout << "Coefficient " << i << " has val=" << new_val << " where expected " << old_val << endl;
      }
    }//end test against fit_amp_and_offset
#endif // NDEBUG
     */
    
    
    for( int i = 0; i < num_polynomial_terms; ++i )
    {
      const T &val = continuum_coeffs[i];
      if( isinf(val) || isnan(val) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN continuum coefficient for range "
                            + std::to_string(range.lower_energy) + " to "
                            + std::to_string(range.upper_energy) + " keV" );
    }//for( const double &val : continuum_coeffs )
    
    // TODO: - currently not defining degrees of freedom well - not using number of relative efficiency terms, or FWHM terms at all, and just blindly using all activity and free peak terms.
    const double approx_dof = 1.0*range.num_channels - nuclides_used.size() - num_polynomial_terms - num_free_peak_pars;
    
    std::sort( begin(peaks), end(peaks), []( const auto &lhs, const auto &rhs ){
      return lhs.mean() < rhs.mean();
    } );
    
    return answer;
  }//PeaksForEnergyRange peaks_for_energy_range_imp( const double lower_energy, const double upper_energy ) const
  
  
  
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
    bool forced_full_range;
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
    const PeaksForEnergyRangeImp<double> computed_peaks = peaks_for_energy_range_imp( range, x );
    
    
    PeaksForEnergyRange answer;
    answer.first_channel = computed_peaks.first_channel;
    answer.last_channel = computed_peaks.last_channel;
    answer.no_gammas_in_range = computed_peaks.no_gammas_in_range;
    answer.forced_full_range = computed_peaks.forced_full_range;
    
    for( size_t i = 0; i < computed_peaks.peaks.size(); ++i )
    {
      const PeakDefImp<double> &comp_peak = computed_peaks.peaks[i];
      PeakDef peak( comp_peak.m_mean, comp_peak.m_sigma, comp_peak.m_amplitude );
      peak.setSkewType( comp_peak.m_skew_type );
      const size_t num_skew = PeakDef::num_skew_parameters( comp_peak.m_skew_type );
      for( size_t skew_index = 0; skew_index < num_skew; ++skew_index )
      {
        const PeakDef::CoefficientType coef
             = static_cast<PeakDef::CoefficientType>( static_cast<int>(PeakDef::CoefficientType::SkewPar0) + num_skew );
        peak.set_coefficient( comp_peak.m_skew_pars[skew_index], coef );
      }//for( size_t skew_index = 0; skew_index < num_skew; ++skew_index )
      
      peak.setNuclearTransition( comp_peak.m_parent_nuclide, comp_peak.m_transition,
                                static_cast<int>(comp_peak.m_rad_particle_index), comp_peak.m_gamma_type );
      
      if( i != 0 )
      {
        peak.setContinuum( answer.peaks.front().continuum() );
      }else
      {
        shared_ptr<PeakContinuum> cont = peak.continuum();
        assert( cont );
        const PeakContinuumImp<double> &comp_cont = computed_peaks.continuum;
        cont->setRange( comp_cont.m_lower_energy, comp_cont.m_upper_energy );
        cont->setType( comp_cont.m_type );
        const size_t npar = PeakContinuum::num_parameters( comp_cont.m_type );
        assert( npar <= comp_cont.m_values.size() );
        const vector<double> pars( begin(comp_cont.m_values), begin(comp_cont.m_values) + npar );
        cont->setParameters( comp_cont.m_reference_energy, pars, {} );
      }//if( i != 0 )
      
      if( comp_peak.m_rel_eff_index < m_nuclides.size() )
      {
        const NucInputGamma *nuc_info = nullptr;
        const vector<NucInputGamma> &rel_eff_nucs = m_nuclides[comp_peak.m_rel_eff_index];
        for( size_t nuc_index = 0; !nuc_info && (nuc_index < rel_eff_nucs.size()); ++nuc_index )
        {
          if( rel_eff_nucs[nuc_index].nuclide == comp_peak.m_parent_nuclide )
            nuc_info = &(rel_eff_nucs[nuc_index]);
        }//
        
        assert( nuc_info );
        if( nuc_info && !nuc_info->peak_color_css.empty() )
          peak.setLineColor( Wt::WColor( Wt::WString::fromUTF8(nuc_info->peak_color_css) ) );
      }//if( comp_peak.m_rel_eff_index < m_nuclides.size() )
      
      answer.peaks.push_back( peak );
    }//for( size_t i = 0; i < computed_peaks.peaks.size(); ++i )
    
    std::sort( begin(answer.peaks), end(answer.peaks), &PeakDef::lessThanByMean );
    
    return answer;
  }//PeaksForEnergyRange peaks_for_energy_range( const double lower_energy, const double upper_energy ) const
  
  
  template<typename T>
  void eval( const std::vector<T> &x, T *residuals ) const
  {
    m_ncalls += 1;
    
    const auto start_time = std::chrono::high_resolution_clock::now();
    DoWorkOnDestruct incrementTimeInFncn( [this,start_time](){
      const auto end_time = std::chrono::high_resolution_clock::now();
      m_nanoseconds_spent_in_eval += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    });
    
    
    assert( x.size() == number_parameters() );
    assert( residuals );
    assert( !m_energy_ranges.empty() );
    
    const auto do_cancel = [this,residuals](){
      const size_t n = number_residuals();
      for( size_t i = 0; i < n; ++i )
        residuals[i] = T(0.0);
    };
    
    if( m_cancel_calc && m_cancel_calc->load() )
    {
      do_cancel();
      return;
    }//if( m_cancel_calc && m_cancel_calc->load() )
    
    
    // Zero out residuals, although I'm not sure if we actually need to do this.
    const size_t num_residuals = number_residuals();
    if constexpr ( !std::is_same_v<T, double> )
    {
      for( size_t i = 0; i < num_residuals; ++i )
        residuals[i] = T(0.0);
    }else
    {
      memset( residuals, 0, sizeof(double) * num_residuals );
    }
    
    // We'll do the simplest parallelization we can by computing peaks in multiple threads; this
    //  actually seems to be reasonably effective in filling up the CPU cores during fitting.
    //  (setting ceres_options.num_threads >1 doesnt seem to do much (any?) good)
    SpecUtilsAsync::ThreadPool pool;
    vector<PeaksForEnergyRangeImp<T>> peaks_in_ranges_imp( m_energy_ranges.size() );
    
    for( size_t i = 0; i < m_energy_ranges.size(); ++i )
    {
      pool.post( [i,&peaks_in_ranges_imp,this,&x](){
        peaks_in_ranges_imp[i] = peaks_for_energy_range_imp( m_energy_ranges[i], x );
      } );
    }//
    pool.join();
    
    
    if( m_cancel_calc && m_cancel_calc->load() )
    {
      do_cancel();
      return;
    }//if( m_cancel_calc && m_cancel_calc->load() )
    
    assert( SpecUtilsAsync::num_logical_cpu_cores() > 0 );

    
    assert( peaks_in_ranges_imp.size() == m_energy_ranges.size() );
    
    size_t residual_index = 0;
    
    const bool try_better_par = true;
    vector<mutex> roi_mutexes( m_energy_ranges.size() );
    
    for( size_t roi_index = 0; roi_index < m_energy_ranges.size(); ++roi_index )
    {
      const RoiRangeChannels &energy_range = m_energy_ranges[roi_index];
      const PeaksForEnergyRangeImp<T> &info = peaks_in_ranges_imp[roi_index];
      
      assert( info.peaks.size() );
      assert( info.first_channel < m_channel_counts.size() );
      assert( info.last_channel < m_channel_counts.size() );
      assert( info.last_channel >= info.first_channel );
      assert( m_channel_counts.size() == m_channel_count_uncerts.size() );
      
      const size_t residual_start_index = residual_index;
      
      const size_t nchannels_roi = (info.last_channel - info.first_channel) + 1;
      residual_index += nchannels_roi;
      
      assert( info.peak_counts.size() == nchannels_roi );
      
      T * const this_residual = residuals + residual_start_index;
      for( size_t index = 0; index < nchannels_roi; ++index )
      {
        const size_t data_index = info.first_channel + index;
        
        const double data_counts = m_channel_counts[data_index];
        const double data_uncert = m_channel_count_uncerts[data_index];
        const T peak_area = info.peak_counts[index];
        
        this_residual[index] = (data_counts - peak_area) / data_uncert;
      }
    }//for( loop over m_energy_ranges )
    
    
    // See TODO above about calculations methods giving slightly different end-results, unless be
    //  truncate the accuracy of the residuals to be floats.
    //cerr << "\n\nRounding residuals to floats for debug\n\n" << endl;
    //const size_t nresid = number_residuals();
    //for( size_t i = 0; (i+1) < nresid; ++i )
    //  residuals[i] = static_cast<float>(residuals[i]);
    
    for( size_t rel_eff_index = 0; rel_eff_index < m_options.rel_eff_curves.size(); ++rel_eff_index )
    {
      const auto &rel_eff = m_options.rel_eff_curves[rel_eff_index];
      
      if( rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        assert( rel_eff_index
               || ((residual_index + m_options.rel_eff_curves.size() + m_peak_ranges_with_uncert.size()) == number_residuals()) );
        
        // Now make sure the relative efficiency curve is anchored to 1.0 (this removes the degeneracy
        //  between the relative efficiency amplitude, and the relative activity amplitudes).
        const double lowest_energy = m_energy_ranges.front().lower_energy;
        const T lowest_energy_rel_eff = relative_eff( lowest_energy, rel_eff_index, x );
        residuals[residual_index] = m_rel_eff_anchor_enhancement * (1.0 - lowest_energy_rel_eff);
        ++residual_index;
      }//if( m_options.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    }//for( size_t rel_eff_index = 0; rel_eff_index < m_options.rel_eff_curves.size(); ++rel_eff_index )
    
    if( !m_peak_ranges_with_uncert.empty() )
    {
      assert( m_options.additional_br_uncert > 0.0 );
      
      for( size_t range_index = 0; range_index < m_peak_ranges_with_uncert.size(); ++range_index )
      {
        const size_t par_index = m_add_br_uncert_start_index + range_index;
        assert( par_index < x.size() );
        
        // Normalize the parameter so nominal is 0, and a 100% change is +-1;
        const T norm_val = x[par_index]/sm_peak_range_uncert_par_offset - 1.0;
        residuals[residual_index] = norm_val / m_options.additional_br_uncert;
        ++residual_index;
      }
    }//if( !m_peak_ranges_with_uncert.empty() )
    
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
  template<typename T>
  bool operator()( T const *const *parameters, T *residuals ) const
  {
    try
    {
      const size_t num_pars = number_parameters();
      vector<T> pars( parameters[0], parameters[0] + num_pars );
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
    MaterialDB matdb;
    const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
    try
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      matdb.parseGadrasMaterialFile( materialfile, db, false );
    }catch( std::exception &e )
    {
      throw runtime_error( "Couldnt initialize material database." );
    }
    
    
    
    const char *xml_file_path = "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/simple_pu_test.xml";
    //const char *xml_file_path = "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/thor_core_614_668_kev_test.xml";
    //const char *xml_file_path = "/Users/wcjohns/rad_ana/InterSpec_RelAct/RelActTest/LaBr_pu_test.xml";
    
    rapidxml::file<char> input_file( xml_file_path );
    
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( input_file.data() );
    
    const auto base_node = get_required_node(&doc, "RelActCalcAuto");
  
    RelActCalcAuto::RelActAutoGuiState state;
    state.deSerialize( base_node, &matdb );
    
    
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
    
    bool extract_info_from_n42 = false;
    try
    {
      extract_info_from_n42 = get_bool_node_value(base_node, "RoiAndNucsFromFile");
    }catch(std::exception &)
    {
      //<RoiAndNucsFromFile> is optional, so will get here if it doesnt exist (or invalid value in it)
    }
    
    vector<string> input_warnings;
    Options options = state.options;
    
    if( state.options.rel_eff_curves.size() != 1 )
      throw runtime_error( "run_tests: onle supports exactly one rel-eff curve" );
    
    if( extract_info_from_n42 )
    {
      vector<RoiRange> energy_ranges = options.rois;
      vector<FloatingPeak> extra_peaks = options.floating_peaks;
      
      const RelActCalcAuto::RelEffCurveInput &rel_eff = state.options.rel_eff_curves[0];
      vector<NucInputInfo> nuclides = rel_eff.nuclides;
      
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
      if( rel_eff.nucs_of_el_same_age )
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
      
      
      
      options.rois = energy_ranges;
      options.floating_peaks = extra_peaks;
      
      assert( options.rel_eff_curves.size() == 1 );
      if( options.rel_eff_curves.size() != 1 )
        throw runtime_error( "options.rel_eff_curves.size() != 1" );
      auto &updated_rel_eff = options.rel_eff_curves[0];
      updated_rel_eff.nuclides = nuclides;
    }//if( extract_info_from_n42 )
    
    sort_rois_by_energy( options.rois );
    
    std::sort( begin(options.floating_peaks), end(options.floating_peaks),
      []( const FloatingPeak &lhs, const FloatingPeak &rhs ) -> bool {
        return lhs.energy < rhs.energy;
    });
    
    
    // Make sure ranges dont overlap - if they do, split the difference.
    //  This only catches simple small
    for( size_t i = 1; i < options.rois.size(); ++i )
    {
      RoiRange &prev_range = options.rois[i-1];
      RoiRange &this_range = options.rois[i];
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
    
    
    if( options.rois.empty() )
      throw runtime_error( "No RoiRanges specified" );
    
    if( options.rel_eff_curves.empty() )
      throw runtime_error( "No RelEffCurveInput specified" );
    
    for( const auto &rel_eff : options.rel_eff_curves )
    {
      if( rel_eff.nuclides.empty() )
        throw runtime_error( "No nuclides specified" );
    }//for( const auto &rel_eff : options.rel_eff_curves )
    
    
    // A helper function to print out the XML; right now just to stdout, but could be useful
    auto print_xml = [&](){
      try
      {
        using namespace rapidxml;
        
        rapidxml::xml_document<char> doc;
        
        RelActCalcAuto::RelActAutoGuiState updated_state;
        
        updated_state.options = options;
        xml_node<char> *base_node = updated_state.serialize( &doc );
        
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
    RelActAutoSolution solution = solve( options, fore_meas, back_meas, drf, all_peaks );
    
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
    case FwhmForm::SqrtEnergyPlusInverse:  return "SqrtEnergyPlusInverse";
    case FwhmForm::ConstantPlusSqrtEnergy: return "ConstantPlusSqrtEnergy";
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

  
/** This function is the same as `DetectorPeakResponse::peakResolutionFWHM(...)`, but templated to allow Jets
 TODO: refactor this function and the equivalent `DetectorPeakResponse` function into a single imlpementation.
 */
template<typename T>
T peakResolutionFWHM( T energy, DetectorPeakResponse::ResolutionFnctForm fcnFrm,
                     const T * const pars, const size_t num_pars )
{
  switch( fcnFrm )
  {
    case DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn:
    {
      if( num_pars != 3 )
        throw std::runtime_error( "RelActCalcAuto::peakResolutionSigma():"
                                 " pars not defined" );
      const T &a = pars[0];
      const T &b = pars[1];
      const T &c = pars[2];
      
      if( (energy >= 661.0) || (abs(a) < T(1.0E-6)) )
        return 6.61 * b * pow(energy/661.0, c);
      
      if( a < 0.0 )
      {
        const T p = pow( c, T(1.0/log(1.0-a)) );
        return 6.61 * b * pow(energy/661.0, p);
      }//if( a < 0.0 )
      
      if( a > 6.61*b )
        return a;
      
      const T A7 = sqrt( pow(6.61*b, 2.0) - a*a )/6.61;
      return sqrt(a*a + pow(6.61 * A7 * pow(energy/661.0, c), 2.0));
    }//case kGadrasResolutionFcn:
      
    case DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse:
    {
      if( num_pars != 3 )
        throw std::runtime_error( "RelActCalcAuto::peakResolutionSigma():"
                                 " pars not defined" );
      energy /= PhysicalUnits::keV;
      
      return sqrt(pars[0] + pars[1]*energy + pars[2]/energy);
    }//case kSqrtEnergyPlusInverse:
      
    case DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy:
    {
      if( num_pars != 2 )
        throw std::runtime_error( "RelActCalcAuto::peakResolutionSigma():"
                                 " pars not defined" );
      energy /= PhysicalUnits::keV;
      
      return pars[0] + pars[1]*sqrt(energy);
    }//case kConstantPlusSqrtEnergy:
      
    case DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial:
    {
      if( num_pars < 1 )
        throw runtime_error( "RelActCalcAuto::peakResolutionSigma():"
                            " pars not defined" );
      
      energy /= PhysicalUnits::MeV;
      //return  A1 + A2*std::pow( energy + A3*energy*energy, A4 );
      
      T val(0.0);
      for( size_t i = 0; i < num_pars; ++i )
        val += pars[i] * pow(energy, static_cast<double>(i) );
      return sqrt( val );
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm:
      throw std::runtime_error( "RelActCalcAuto::peakResolutionSigma():"
                               " Resolution not defined" );
      break;
  }//switch( m_resolutionForm )
  
  //Lets keep MSVS happy
  assert(0);
  return T(0.0);
}//static float peakResolutionFwhmGadras(...)

  
template<typename T>
T eval_fwhm( const T energy, const FwhmForm form, const T * const pars, const size_t num_pars )
{
  DetectorPeakResponse::ResolutionFnctForm fctntype = DetectorPeakResponse::kNumResolutionFnctForm;
  switch( form )
  {
    case RelActCalcAuto::FwhmForm::Gadras:
      fctntype = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
      break;
      
    case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
      fctntype = DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse;
      break;
      
    case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
      fctntype = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
      break;
      
    case RelActCalcAuto::FwhmForm::Polynomial_2:
    case RelActCalcAuto::FwhmForm::Polynomial_3:
    case RelActCalcAuto::FwhmForm::Polynomial_4:
    case RelActCalcAuto::FwhmForm::Polynomial_5:
    case RelActCalcAuto::FwhmForm::Polynomial_6:
      fctntype = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      break;
  }//switch( m_options.fwhm_form )
  
  assert( fctntype != DetectorPeakResponse::kNumResolutionFnctForm );
  if( fctntype == DetectorPeakResponse::kNumResolutionFnctForm )
    throw runtime_error( "eval_fwhm: invalid FwhmForm" );
  
  //return DetectorPeakResponse::peakResolutionFWHM( energy, fctntype, drfx );
  return peakResolutionFWHM( energy, fctntype, pars, num_pars );
}//float eval_fwhm( const float energy, const FwhmForm form, const vector<float> &drfx )


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
    XML_FOREACH_CHILD( energy_node, exclude_node, "Energy" ) //ok if exclude_node is nullptr
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
  append_bool_node( base_node, "ApplyEnergyCalCorrection", apply_energy_cal_correction );
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
    
    try
    {
      // Dont require there to be a <ApplyEnergyCalCorrection> node
      apply_energy_cal_correction = get_bool_node_value( parent, "ApplyEnergyCalCorrection" );
    }catch( ... )
    {
      apply_energy_cal_correction = true;
    }
  }catch( std::exception &e )
  {
    throw runtime_error( "FloatingPeak::fromXml(): " + string(e.what()) );
  }
}//void fromXml(...)


rapidxml::xml_node<char> *Options::toXml( rapidxml::xml_node<char> *parent ) const
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
  
  const char *fwhm_form_str = to_str( fwhm_form );
  auto fwhm_node = append_string_node( base_node, "FwhmForm", fwhm_form_str );
  append_attrib( fwhm_node, "remark", "Possible values: Gadras, Polynomial_2, Polynomial_3, Polynomial_4, Polynomial_5, Polynomial_6" );
  
  append_string_node( base_node, "Title", spectrum_title );
  
  append_string_node( base_node, "SkewType", PeakDef::to_string(skew_type) );
  
  append_float_node( base_node, "AddUncert", additional_br_uncert );
  
  
  xml_node<char> *rel_eff_node = doc->allocate_node( node_element, "RelEffCurveInputs" );
  base_node->append_node( rel_eff_node );
  for( const auto &curve : rel_eff_curves )
    curve.toXml( rel_eff_node );

  if( !rois.empty() )
  {
    xml_node<char> *node = doc->allocate_node( node_element, "RoiRangeList" );
    base_node->append_node( node );
    for( const auto &range : rois )
      range.toXml( node );
  }
  
  if( !floating_peaks.empty() )
  {
    xml_node<char> *node = doc->allocate_node( node_element, "FloatingPeakList" );
    base_node->append_node( node );
    for( const auto &peak : floating_peaks )
      peak.toXml( node );
  }//if( !floating_peaks.empty() )
  
  
  return base_node;
}//rapidxml::xml_node<char> *Options::toXml(...)


void Options::fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB )
{
  try
  {
    if( !parent )
      throw runtime_error( "invalid input" );
    
    if( !rapidxml::internal::compare( parent->name(), parent->name_size(), "Options", 7, false ) )
      throw std::logic_error( "invalid input node name" );
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( Options::sm_xmlSerializationVersion == 2,
                  "needs to be updated for new serialization version." );
    
    check_xml_version( parent, Options::sm_xmlSerializationVersion );
    
    fit_energy_cal = get_bool_node_value( parent, "FitEnergyCal" );
    
    const rapidxml::xml_node<char> *fwhm_node = XML_FIRST_NODE( parent, "FwhmForm" );
    const string fwhm_str = SpecUtils::xml_value_str( fwhm_node );
    fwhm_form = fwhm_form_from_str( fwhm_str.c_str() );
    
    const rapidxml::xml_node<char> *title_node = XML_FIRST_NODE( parent, "Title" );
    spectrum_title = SpecUtils::xml_value_str( title_node );
    
    // skew_type added 20231111; we wont require it.
    skew_type = PeakDef::SkewType::NoSkew;
    const rapidxml::xml_node<char> *skew_node = XML_FIRST_NODE( parent, "SkewType" );
    const string skew_str = SpecUtils::xml_value_str( skew_node );
    if( !skew_str.empty() )
      skew_type = PeakDef::skew_from_string( skew_str );
    
    // Additional uncertainty added 202250125, so we'll allow it to be missing.
    try
    {
      additional_br_uncert = get_float_node_value( parent, "AddUncert" );
      additional_br_uncert = std::max( additional_br_uncert, 0.0 );
    }catch( ... )
    {
      additional_br_uncert = 0.0;
    }
    
    rel_eff_curves.clear();
    const int version = get_int_attribute( parent, "version" );
    if( version >= 2 )
    {
      const rapidxml::xml_node<char> *rel_eff_node = XML_FIRST_NODE( parent, "RelEffCurveInputs" );
      if( rel_eff_node )
      {
        XML_FOREACH_CHILD( curve_node, rel_eff_node, "RelEffCurveInput" )
        {
          RelEffCurveInput curve;
          curve.fromXml( curve_node, materialDB );
          rel_eff_curves.push_back( curve );
        }
      }//if( rel_eff_node )
    }else
    {
      RelEffCurveInput rel_eff_input;
      rel_eff_input.fromXml( parent, materialDB );
      rel_eff_curves.push_back( rel_eff_input );
    }//if( version >= 2 ) / else

    const rapidxml::xml_node<char> *node = XML_FIRST_NODE(parent, "RoiRangeList");
    if( !node && parent->parent() )
    {
      assert( get_int_attribute(parent, "version") <= 2 );
      node = XML_FIRST_NODE( parent->parent(), "RoiRangeList" );
    }
    
    if( node )
    {
      XML_FOREACH_CHILD( roi_node, node, "RoiRange" )
      {
        RelActCalcAuto::RoiRange roi;
        roi.fromXml( roi_node );
        rois.push_back( roi );
      }
    }//if( <RoiRangeList> )
  
  
    node = XML_FIRST_NODE(parent, "FloatingPeakList");
    if( !node && parent->parent() )
    {
      assert( get_int_attribute(parent, "version") <= 2 );
      node = XML_FIRST_NODE( parent->parent(), "FloatingPeakList" );
    }
    
    if( node )
    {
       XML_FOREACH_CHILD( peak_node, node, "FloatingPeak" )
      {
        RelActCalcAuto::FloatingPeak peak;
        peak.fromXml( peak_node );
        floating_peaks.push_back( peak );
      }
    }//if( <FloatingPeakList> )
  }catch( std::exception &e )
  {
    throw runtime_error( "Options::fromXml(): " + string(e.what()) );
  }
}//void Options::fromXml( const ::rapidxml::xml_node<char> *parent )


size_t num_parameters( const FwhmForm eqn_form )
{
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Polynomial_2) == 3, "FwhmForm enum needs updating" );
  
  switch( eqn_form )
  {
    case FwhmForm::Gadras:       return 3;
    case FwhmForm::SqrtEnergyPlusInverse:  return 3;
    case FwhmForm::ConstantPlusSqrtEnergy: return 2;
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

RelEffCurveInput::RelEffCurveInput()
: nuclides{},
  nucs_of_el_same_age( true ),
  rel_eff_eqn_type( RelActCalc::RelEffEqnForm::LnX ),
  rel_eff_eqn_order( 0 ),
  phys_model_self_atten( nullptr ),
  phys_model_external_atten{},
  phys_model_use_hoerl( true ),
  pu242_correlation_method( RelActCalc::PuCorrMethod::NotApplicable )
{
}

rapidxml::xml_node<char> *RelEffCurveInput::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "RoiRange::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "RelEffCurveInput" );
  parent->append_node( base_node );
  append_version_attrib( base_node, RelEffCurveInput::sm_xmlSerializationVersion );

  if( !nuclides.empty() )
  {
    xml_node<char> *nuc_node = doc->allocate_node( node_element, "NucInputInfoList" );
    base_node->append_node( nuc_node );
    for( const auto &nuc : nuclides )
      nuc.toXml( nuc_node );
  }

  append_bool_node( base_node, "NucsOfElSameAge", nucs_of_el_same_age );
  
  const char *rell_eff_eqn_str = RelActCalc::to_str( rel_eff_eqn_type );
  auto rel_eff_node = append_string_node( base_node, "RelEffEqnType", rell_eff_eqn_str);
  append_attrib( rel_eff_node, "remark", "Possible values: Possible values: LnX, LnY, LnXLnY, \"FRAM Empirical\", \"FRAM Physical\"" );
  
  if( (rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel) || (rel_eff_eqn_order > 0) )
    append_int_node( base_node, "RelEffEqnOrder", rel_eff_eqn_order);
  

  // Even if not RelActCalc::RelEffEqnForm::FramPhysicalModel, we'll still write the state of
  //  the shielding widgets to the XML.
  if( phys_model_self_atten || !phys_model_external_atten.empty() )
  {
    xml_node<char> *phys_node = doc->allocate_node( node_element, "PhysicalModelOptions" );
    base_node->append_node( phys_node );
    
    if( phys_model_self_atten )
    {
      xml_node<char> *self_atten_node = doc->allocate_node( node_element, "SelfAtten" );
      phys_node->append_node( self_atten_node );
      phys_model_self_atten->toXml( self_atten_node );
    }//if( phys_model_self_atten )
    
    if( !phys_model_external_atten.empty() )
    {
      xml_node<char> *ext_atten_node = doc->allocate_node( node_element, "ExtAtten" );
      phys_node->append_node( ext_atten_node );
      for( const auto &ext : phys_model_external_atten )
        ext->toXml( ext_atten_node );
    }//if( !phys_model_external_atten.empty() )
    
    XmlUtils::append_bool_node( phys_node, "PhysModelUseHoerl", phys_model_use_hoerl );
  }//if( phys_model_self_atten || !phys_model_external_atten.empty() )

  if( pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
  {
    const string &method_str = RelActCalc::to_str( pu242_correlation_method );
    append_string_node( base_node, "Pu242CorrMethod", method_str );
  }

  return base_node;
}//RelEffCurveInput::toXml(...)


void RelEffCurveInput::fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB )
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "RelEffCurveInput::fromXml: invalid parent." );

  // We'll allow the node name to be either "RelEffCurveInput" or "Options" nodes as input,
  //  but if its "Options", then it must be version 0 or 1.
  const bool is_rel_eff_curve_input = rapidxml::internal::compare( parent->name(), parent->name_size(), "RelEffCurveInput", 16, false );  
  const bool is_options = !is_rel_eff_curve_input && rapidxml::internal::compare( parent->name(), parent->name_size(), "Options", 7, false );
  if( !is_rel_eff_curve_input && !is_options )
    throw std::logic_error( "invalid input node name" );
  
  check_xml_version( parent, is_rel_eff_curve_input ? RelEffCurveInput::sm_xmlSerializationVersion : 1 );
  

  const rapidxml::xml_node<char> *node = XML_FIRST_NODE(parent, "NucInputInfoList");
  if( !node && parent->parent() && is_options )  //Compatibility with Options::sm_xmlSerializationVersion versions 0 and 1
    node = XML_FIRST_NODE(parent->parent(), "NucInputInfoList");
  
  if( node )
  {
    XML_FOREACH_CHILD( nuc_node, node, "NucInputInfo" )
    {
      RelActCalcAuto::NucInputInfo nuc;
      nuc.fromXml( nuc_node );
      nuclides.push_back( nuc );
    }
  }//if( <NucInputInfoList> )

  nucs_of_el_same_age = get_bool_node_value( parent, "NucsOfElSameAge" );
    
  const rapidxml::xml_node<char> *rel_eff_eqn_node = XML_FIRST_NODE( parent, "RelEffEqnType" );
  const string rel_eff_str = SpecUtils::xml_value_str( rel_eff_eqn_node );
    
  rel_eff_eqn_type = RelActCalc::rel_eff_eqn_form_from_str( rel_eff_str.c_str() );
    
  const rapidxml::xml_node<char> *rel_eff_order_node = XML_FIRST_NODE( parent, "RelEffEqnOrder" );

  if( (rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel) || rel_eff_order_node )
  {
    const string rel_eff_order_str = SpecUtils::xml_value_str( rel_eff_order_node );
    if( !(stringstream(rel_eff_order_str) >> rel_eff_eqn_order) )
      throw runtime_error( "Invalid 'RelEffEqnOrder' value '" + rel_eff_order_str + "'" );
  }else
  {
    rel_eff_eqn_order = 0;
  }
  

  // Even if not RelActCalc::RelEffEqnForm::FramPhysicalModel, we'll still try to grab the state
  phys_model_use_hoerl = true;
  phys_model_self_atten.reset();
  phys_model_external_atten.clear();
  const rapidxml::xml_node<char> *phys_model_node = XML_FIRST_NODE( parent, "PhysicalModelOptions" );
    
  if( (rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) && !phys_model_node )
    throw runtime_error( "Missing PhysicalModelOptions node" );
    
  if( phys_model_node )
  {
    const rapidxml::xml_node<char> *self_atten_par_node = XML_FIRST_NODE( phys_model_node, "SelfAtten" );
    const rapidxml::xml_node<char> *self_atten_node = SpecUtils::xml_first_node(self_atten_par_node, "PhysicalModelShield");
      
    if( self_atten_node )
    {
      auto att = make_shared<RelActCalc::PhysicalModelShieldInput>();
      att->fromXml( self_atten_node, materialDB );
      phys_model_self_atten = att;
    }//if( self_atten_node )
      
      
    const rapidxml::xml_node<char> *ext_atten_node = XML_FIRST_NODE( phys_model_node, "ExtAtten" );
    if( ext_atten_node )
    {
      XML_FOREACH_CHILD( att_node, ext_atten_node, "PhysicalModelShield" )
      {
        auto att = make_shared<RelActCalc::PhysicalModelShieldInput>();
        att->fromXml( att_node, materialDB );
        phys_model_external_atten.push_back( att );
      }
    }//if( ext_atten_node )
      
    phys_model_use_hoerl = XmlUtils::get_bool_node_value(phys_model_node, "PhysModelUseHoerl" );
  }//if( phys_model_node )

  pu242_correlation_method = RelActCalc::PuCorrMethod::NotApplicable;
  const rapidxml::xml_node<char> *pu242_corr_node = XML_FIRST_NODE( parent, "Pu242CorrMethod" );
  const string pu242_corr_str = SpecUtils::xml_value_str( pu242_corr_node );
  if( !pu242_corr_str.empty() )
  {
    bool found = false;
    for( int i = 0; i <= static_cast<int>(RelActCalc::PuCorrMethod::NotApplicable); ++i )
    {
      const auto method = RelActCalc::PuCorrMethod(i);
      const std::string &method_str = RelActCalc::to_str( method );
      if( SpecUtils::iequals_ascii(pu242_corr_str, method_str) )
      {
        pu242_correlation_method = method;
        found = true;
        break;
      }
    }//for( loop over RelActCalc::PuCorrMethod types )
      
    assert( found );
  }//if( !pu242_corr_str.empty() )
}//void RelEffCurveInput::fromXml(...)
  
  

Options::Options()
: fit_energy_cal( false ),
  fwhm_form( FwhmForm::Polynomial_2 ),
  spectrum_title( "" ),
  skew_type( PeakDef::SkewType::NoSkew ),
  additional_br_uncert( 0.0 ),
  rel_eff_curves{}
{
}
  
  
RelActAutoGuiState::RelActAutoGuiState()
  : options{},
  background_subtract( false ),
  show_ref_lines( false ),
  lower_display_energy( 0.0 ),
  upper_display_energy( 0.0 )
{
}
  

::rapidxml::xml_node<char> *RelActAutoGuiState::serialize( ::rapidxml::xml_node<char> *parent_node ) const
{
  assert( parent_node );
  
  using namespace rapidxml;
  xml_document<char> * const doc = parent_node->document();
  
  assert( doc );
  if( !doc )
    throw runtime_error( "RelActAutoGuiState::serialize: no xml_document" );
  
  xml_node<char> *base_node = doc->allocate_node( node_element, "RelActCalcAuto" );
  parent_node->append_node( base_node );
  
  xml_attribute<char> *attrib = doc->allocate_attribute( "version", "0" );
  base_node->append_attribute( attrib );
  
  // Elements in the offline-xml we dont deal with here
  //  base_node->append_node( <"ForegroundFileName"> );
  //  base_node->append_node( <"BackgroundFileName"> );
  //  base_node->append_node( <"OutputHtmlFileName"> );
  //  base_node->append_node( <"RoiAndNucsFromFile"> );
  
  rapidxml::xml_node<char> *options_node = options.toXml( base_node );
  
  {// Begin put extra options
    if( background_subtract )
    {
      const char *val = background_subtract ? "true" : "false";
      xml_node<char> *node = doc->allocate_node( node_element, "BackgroundSubtract", val );
      options_node->append_node( node );
    }//if( we have a background to subtract )
  }// End put extra options
  
  
  {//Begin show_ref_lines
    const char *val = show_ref_lines ? "true" : "false";
    xml_node<char> *node = doc->allocate_node( node_element, "ShowRefGammaLines", val );
    options_node->append_node( node );
  }//End show_ref_lines
  
  if( lower_display_energy < upper_display_energy )
  {
    const string lower_energy_str = std::to_string( lower_display_energy );
    const string upper_energy_str = std::to_string( upper_display_energy );
    
    const char *val = doc->allocate_string( lower_energy_str.c_str() );
    xml_node<char> *node = doc->allocate_node( node_element, "LowerDisplayEnergy", val );
    options_node->append_node( node );
    
    val = doc->allocate_string( upper_energy_str.c_str() );
    node = doc->allocate_node( node_element, "UpperDisplayEnergy", val );
    options_node->append_node( node );
  }//if( lower_display_energy < upper_display_energy )
  
  return base_node;
}//::rapidxml::xml_node<char> *serialize( ::rapidxml::xml_node<char> *parent ) const
  

void RelActAutoGuiState::deSerialize( const rapidxml::xml_node<char> *base_node, MaterialDB *materialDb )
{
  using namespace rapidxml;
  
  assert( base_node );
  if( !base_node )
    throw runtime_error( "RelActAutoGuiState::deSerialize: nullptr passed in." );
  
  const string base_node_name = SpecUtils::xml_name_str(base_node);
  if( base_node_name != "RelActCalcAuto" )
    throw runtime_error( "RelActAutoGuiState::deSerialize: invalid node passed in named '"
                         + base_node_name + "'" );
  
  const xml_attribute<char> *attrib = XML_FIRST_ATTRIB(base_node, "version");
  const string base_version = SpecUtils::xml_value_str(attrib);
  if( !SpecUtils::istarts_with(base_version, "0") )
    throw runtime_error( "RelActAutoGuiState::deSerialize: invalid xml version='" + base_version + "'" );
  
  const xml_node<char> *node = XML_FIRST_NODE(base_node, "Options");
  if( !node )
    throw runtime_error( "RelActAutoGui::deSerialize: No <Options> node." );
  
  options.fromXml( node, materialDb );
  
  
  {// Begin get extra options
    const xml_node<char> *opt_node = XML_FIRST_NODE(node, "BackgroundSubtract");
    string val = SpecUtils::xml_value_str( opt_node );
    background_subtract = SpecUtils::iequals_ascii(val, "true");
    
    opt_node = XML_FIRST_NODE(node, "ShowRefGammaLines");
    val = SpecUtils::xml_value_str( opt_node );
    show_ref_lines = SpecUtils::iequals_ascii(val, "true");
    
    opt_node = XML_FIRST_NODE(node, "LowerDisplayEnergy");
    const string lower_display_energy_str = SpecUtils::xml_value_str( opt_node );
    opt_node = XML_FIRST_NODE(node, "UpperDisplayEnergy");
    const string upper_display_energy_str = SpecUtils::xml_value_str( opt_node );
    
    lower_display_energy = upper_display_energy = 0.0;
    if( !lower_display_energy_str.empty()
       && !(stringstream(lower_display_energy_str) >> lower_display_energy) )
    {
      // Maybe give a warning here?
    }
    
    if( !upper_display_energy_str.empty()
       && !(stringstream(upper_display_energy_str) >> upper_display_energy) )
    {
      // Maybe give a warning here?
    }
    
    
    /*
    opt_node = XML_FIRST_NODE(node, "UPuDataSrc");
    val = SpecUtils::xml_value_str( opt_node );
    if( !val.empty() )
    {
      using RelActCalcManual::PeakCsvInput::NucDataSrc;
      
      bool set_src = false;
      for( int i = 0; i < static_cast<int>(NucDataSrc::Undefined); ++i )
      {
        const NucDataSrc src = NucDataSrc(i);
        const char *src_str = RelActCalcManual::PeakCsvInput::to_str(src);
        if( val == src_str )
        {
          set_src = true;
          m_u_pu_data_source->setCurrentIndex( i );
          break;
        }
      }//for( int i = 0; i < static_cast<int>(NucDataSrc::Undefined); ++i )
      
      if( !set_src )
        cerr << "Failed to convert '" << val << "' into a NucDataSrc" << endl;
    }//if( UPuDataSrc not empty )
     */
  }// End get extra options
}//void deSerialize( const rapidxml::xml_node<char> *base_node )
  

RelActAutoSolution::RelActAutoSolution()
: m_status( RelActAutoSolution::Status::NotInitiated ),
  m_error_message( "" ),
  m_foreground{ nullptr },
  m_background{ nullptr },
  m_spectrum{ nullptr },
  m_rel_eff_forms{},
  m_rel_eff_coefficients{},
  m_rel_eff_covariance{},
  m_rel_activities{},
  m_rel_act_covariance{},
  m_fwhm_form( FwhmForm::SqrtEnergyPlusInverse ),
  m_fwhm_coefficients{},
  m_fwhm_covariance{},
  m_floating_peaks{},
  m_fit_peaks{},
  m_drf{ nullptr },
  m_energy_cal_adjustments{},
  m_fit_energy_cal{},
  m_chi2( -1.0 ),
  m_dof( 0 ),
  m_num_function_eval_solution( 0 ),
  m_num_function_eval_total( 0 ),
  m_num_microseconds_eval( 0 ),
  m_num_microseconds_in_eval( 0 )
{
  for( auto &i : m_energy_cal_adjustments )
    i = 0.0;
  for( auto &i : m_fit_energy_cal )
    i = false;
}
  
std::string RelActAutoSolution::rel_eff_txt( const bool html_format, const size_t rel_eff_index ) const
{
  assert( m_rel_eff_forms.size() == m_rel_eff_coefficients.size() );
  assert( m_rel_eff_forms.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_options.rel_eff_curves.size() );
  
  if( rel_eff_index >= m_options.rel_eff_curves.size() )
    throw logic_error( "RelActAutoSolution::rel_eff_txt: invalis rel eff index" );
  
  const RelActCalc::RelEffEqnForm eqn_form = m_rel_eff_forms[rel_eff_index];
  const vector<double> &coeffs = m_rel_eff_coefficients[rel_eff_index];
  
  if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::rel_eff_eqn_text( eqn_form, coeffs);
  
  const RelActAutoCostFcn::PhysModelRelEqnDef phys_in
           = RelActAutoCostFcn::make_phys_eqn_input( m_options.rel_eff_curves[rel_eff_index], m_drf, coeffs, 0 );
  
  return RelActCalc::physical_model_rel_eff_eqn_text( phys_in.self_atten, phys_in.external_attens,
                                phys_in.det, phys_in.hoerl_b, phys_in.hoerl_c, html_format );
}//std::string RelActAutoSolution::rel_eff_txt()
  

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
      
    case Status::UserCanceled:
      out << "User cancelled calculation.\n";
      break;
  }//switch( m_status )
  
  if( m_error_message.size() )
    out << "Error: " << m_error_message << "\n";
  
  for( const string &warning: m_warnings )
    out << "Warning: " << warning << "\n";
  
  if( m_status != Status::Success )
    return out;
  
  assert( m_final_parameters.size() == m_parameter_names.size() );
  out << "Raw Ceres par values: [";
  for( size_t i = 0; i < m_final_parameters.size(); ++i )
    out << ((i > 0) ? ", " : "") << "{" << m_parameter_names[i] << "=" << m_final_parameters[i] 
        << "," << (m_parameter_were_fit[i] ? "Fit" : "NotFit") << "}";
  out << "]\n\n";

  // Rake code from RelEff
  const size_t num_rel_eff = m_options.rel_eff_curves.size();
  for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff; ++rel_eff_index )
  {
    out << "Rel. Eff. Eqn." << ((num_rel_eff > 1) ? (" " + std::to_string(rel_eff_index)) : string()) << ": y = ";
    out << rel_eff_txt(false, rel_eff_index) << "\n";
  }
  
  bool ene_cal_fit = false;
  for( const bool &fit : m_fit_energy_cal )
    ene_cal_fit = (ene_cal_fit || fit);
  
  if( ene_cal_fit )
  {
    static_assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 3 );
    const char * const label[3] = { "offset", "gain", "cubic" };
    const char * const post_script[3] = { "keV", "keV/chnl", "keV/chnl2" };
    
    const auto cal = m_spectrum ? m_spectrum->energy_calibration() : nullptr;
    const size_t num_chan = (cal ? cal->num_channels() : 1);
    
    out << "Energy calibration was fit for: ";
    bool printed_out = false;
    for( size_t i = 0; i < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars; ++i )
    {
      if( m_fit_energy_cal[i] )
      {
        double value = (m_energy_cal_adjustments[i]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0);
        if( i == 0 )
          value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
        else if( i == 1 )
          value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / num_chan;
        else if( i == 2 )
          value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / (num_chan * num_chan);
        
        char buffer[32] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%1.6G", value );
        out << (printed_out ? ", " : "") << label[i] << "=" << buffer << " " << post_script[i];
        printed_out = true;
      }
    }//for( loop over energy cal ish )
    
    out << "\n";
  }//if( m_fit_energy_cal[0] || m_fit_energy_cal[1] )
  
  
  for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff; ++rel_eff_index )
  {
    const vector<NuclideRelAct> &rel_acts = m_rel_activities[rel_eff_index];
    if( rel_eff_index )
      out << "\n";
    out << "Activity ratios" << ((num_rel_eff > 1) ? (" " + std::to_string(rel_eff_index)) : string()) << ":\n";
    for( size_t i = 1; i < rel_acts.size(); ++i )
    {
      const NuclideRelAct &nuc_i = rel_acts[i];
      
      for( size_t j = 0; j < i; ++j )
      {
        const NuclideRelAct &nuc_j = rel_acts[j];
        
        const double act_ratio = activity_ratio(nuc_i.nuclide, nuc_j.nuclide, rel_eff_index);
        const double nuc_mass_ratio = mass_ratio(nuc_i.nuclide, nuc_j.nuclide, rel_eff_index);
        
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
  }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff; ++rel_eff_index )
  
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  for( size_t rel_eff_index = 0; rel_eff_index < m_rel_activities.size(); ++rel_eff_index )
  {
    const vector<NuclideRelAct> &rel_acts = m_rel_activities[rel_eff_index];
    
    map<int,vector<const SandiaDecay::Nuclide *>> nucs_of_element;
    const SandiaDecay::Nuclide *u235 = nullptr, *pu240 = nullptr;
    for( const NuclideRelAct &nuc_rel_act : rel_acts )
    {
      const SandiaDecay::Nuclide * const nuc = nuc_rel_act.nuclide;
      
      nucs_of_element[nuc->atomicNumber].push_back( nuc );
      
      if( (nuc->atomicNumber == 92) && (nuc->massNumber == 235) )
        u235 = nuc;
      
      if( (nuc->atomicNumber == 94) && (nuc->massNumber == 240) )
        pu240 = nuc;
    }//for( const NuclideRelAct &nuc : m_rel_activities )
    
    if( u235 )
      out << "\nEnrichment " << 100.0*mass_enrichment_fraction(u235,rel_eff_index) << "% U235\n";
    
    if( pu240 )
      out << "\nEnrichment " << 100.0*mass_enrichment_fraction(pu240,rel_eff_index) << "% Pu240\n";
    
    assert( m_corrected_pu.size() <= m_rel_activities.size() );
    
    // For Pu, print a corrected enrichment table
    const shared_ptr<const RelActCalc::Pu242ByCorrelationOutput> pu_corr = (m_corrected_pu.size() > rel_eff_index)
    ? m_corrected_pu[rel_eff_index]
    : nullptr;
    if( pu_corr )
    {
      out << "Plutonium Enrichment" << ((num_rel_eff > 1) ? (" " + std::to_string(rel_eff_index)) : string()) << ":\n";
      
      set<string> nuclides;
      for( const auto &act : rel_acts )
      {
        if( act.nuclide
           && ((act.nuclide->atomicNumber == 94)
               || ((act.nuclide->atomicNumber == 95) && (act.nuclide->massNumber == 241))) )
        {
          nuclides.insert( act.nuclide->symbol );
        }
      }//for( const auto &act : rel_acts )
      
      if( nuclides.count( "Pu238" ) )
        out << "\tPu238: " << SpecUtils::printCompact(100.0*pu_corr->pu238_mass_frac, 4) << "\n";
      if( nuclides.count( "Pu239" ) )
        out << "\tPu239: " << SpecUtils::printCompact(100.0*pu_corr->pu239_mass_frac, 4) << "\n";
      if( nuclides.count( "Pu240" ) )
        out << "\tPu40: " << SpecUtils::printCompact(100.0*pu_corr->pu240_mass_frac, 4) << "\n";
      if( nuclides.count( "Pu241" ) )
        out << "\tPu241: " << SpecUtils::printCompact(100.0*pu_corr->pu241_mass_frac, 4) << "\n";
      out << "\tPu242 (by corr): " << SpecUtils::printCompact(100.0*pu_corr->pu242_mass_frac, 4) << "\n";
      if( nuclides.count( "Am241" ) )
        out << "\tAm241: " << SpecUtils::printCompact(100.0*pu_corr->am241_mass_frac, 4) << "\n";
      out << "\n";
    }//if( pu_corr )
    
    
    for( auto &v : nucs_of_element )
    {
      vector<const SandiaDecay::Nuclide *> nucs = v.second;
      if( nucs.size() < 2 )
        continue;
      
      std::sort( begin(nucs), end(nucs), [](auto l, auto r){ return l->massNumber < r->massNumber; } );
      
      const SandiaDecay::Element * const el = db->element( nucs.front()->atomicNumber );
      assert( el );
      
      if( pu_corr && (el->atomicNumber == 94) )
        continue;
      
      out << "For " << el->name << ":\n";
      for( const SandiaDecay::Nuclide *nuc : nucs )
        out << std::setw(5) << nuc->symbol << ": "
        << std::setw(10) << std::setprecision(4) << 100.0*mass_enrichment_fraction(nuc,rel_eff_index) << "%"
        << " of the " << el->name << ", by mass.\n";
    }//for( auto &v : nucs_of_element )
  }//for( size_t rel_eff_index = 0; rel_eff_index < m_rel_activities.size(); ++rel_eff_index )
  
  
  // Now print out ratios of activities for Rel Eff cyrves that have same nuclides
  for( size_t outer_re_index = 0; outer_re_index < m_rel_activities.size(); ++outer_re_index )
  {
    for( const NuclideRelAct &outer_act : m_rel_activities[outer_re_index] )
    {
      for( size_t inner_re_index = outer_re_index + 1; inner_re_index < m_rel_activities.size(); ++inner_re_index )
      {
        for( const NuclideRelAct &inner_act : m_rel_activities[inner_re_index] )
        {
          if( outer_act.nuclide != inner_act.nuclide )
            continue;
            
          const double inner_activity = inner_act.rel_activity;
          const double outer_activity = outer_act.rel_activity;
          const double inner_counts = nuclide_counts(inner_act.nuclide, inner_re_index);
          const double outer_counts = nuclide_counts(outer_act.nuclide, outer_re_index);
          
          out << "Ratio of " << outer_act.nuclide->symbol << " between RE curves "
          << outer_re_index << " and " << inner_re_index << ":\n"
          << "\tActivity: " << outer_activity << " / " << inner_activity << " = " << (outer_activity/inner_activity) << "\n"
          << "\tCounts  : " << outer_counts << " / " << inner_counts << " = " << (outer_counts/inner_counts) << "\n"; 
        }//for( const NuclideRelAct &inner_act : m_rel_activities[inner_re_index] )
      }//for( loop over further RE curves )
    }//for( const NuclideRelAct &outer_act : m_rel_activities[outer_re_index] )
  }//for( loop over Rel Eff curves print out ratios of activities between curves )
  
  return out;
}//RelActAutoSolution::print_summary(...)




void RelActAutoSolution::print_html_report( std::ostream &out ) const
{
  // TODO: we should probably use WTemplate to form this HTML - and then we could maybe share some models or whatever with the GUI.
  
  assert( m_rel_eff_forms.size() == m_rel_eff_coefficients.size() );
  assert( m_rel_eff_forms.size() == m_options.rel_eff_curves.size() );
  
  char buffer[1024*16] = { '\0' };
  
  // First, load in template html, and do the problem-independent substitions
  auto load_file_contents = []( string filename ) -> string {
    Wt::WApplication *app = Wt::WApplication::instance();
    const string docroot = app ? app->docRoot() : ".";
    const string resource_path = SpecUtils::append_path( docroot, "InterSpec_resources" );
    const string filepath = SpecUtils::append_path( resource_path, filename );
    
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
  
  
  const float live_time = m_spectrum ? m_spectrum->live_time() : 1.0f;
  
  if( m_options.rel_eff_curves.size() > 1 )
    cerr << "\n\n\nWARNING - only printing results for first rel eff cyrve!!\n\n\n";
#warning "only printing HTML results for first rel eff curve!!"
  
  for( size_t rel_eff_index = 0; (rel_eff_index < 1) && (rel_eff_index < m_rel_eff_forms.size()); ++rel_eff_index )
  {
    const RelEffCurveInput &rel_eff = m_options.rel_eff_curves[rel_eff_index];
    const vector<NuclideRelAct> &rel_activities = m_rel_activities[rel_eff_index];
    
    shared_ptr<const RelActCalc::Pu242ByCorrelationOutput> pu_corr = (rel_eff_index >= m_corrected_pu.size())
                                                                   ? nullptr : m_corrected_pu[rel_eff_index];
    
    stringstream results_html;
    
    if( m_options.rel_eff_curves.size() > 1 )
      results_html << "<p><div style=\"color: red\">Only including results for first Rel. Eff. Curve!!</div></p>\n";
    
    results_html << "<div>&chi;<sup>2</sup>=" << m_chi2 << " and there were " << m_dof
    << " DOF (&chi;<sup>2</sup>/dof=" << m_chi2/m_dof << ")</div>\n";
    
    results_html << "<div class=\"releffeqn\">Rel. Eff. Eqn: y = "
    << rel_eff_txt(true,rel_eff_index)
    << "</div>\n";
    
    // For Pu, print a corrected enrichment table
    if( pu_corr )
    {
      results_html << "<table class=\"nuctable resulttable\">\n";
      results_html << "  <caption>Plutonium mass fractions</caption>\n";
      results_html << "  <thead><tr>"
      " <th scope=\"col\">Nuclide</th>"
      " <th scope=\"col\">% Pu Mass</th>"
      " </tr></thead>\n";
      results_html << "  <tbody>\n";
      
      set<string> nuclides;
      for( const auto &act : rel_activities )
      {
        if( act.nuclide
           && ((act.nuclide->atomicNumber == 94)
               || ((act.nuclide->atomicNumber == 95) && (act.nuclide->massNumber == 241))) )
        {
          nuclides.insert( act.nuclide->symbol );
        }
      }//for( const auto &act : rel_activities )
      
      if( nuclides.count( "Pu238" ) )
        results_html << "  <tr><td>Pu238</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu238_mass_frac, 4)
        << "</td></tr>\n";
      if( nuclides.count( "Pu239" ) )
        results_html << "  <tr><td>Pu239</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu239_mass_frac, 4)
        << "</td></tr>\n";
      if( nuclides.count( "Pu240" ) )
        results_html << "  <tr><td>Pu40</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu240_mass_frac, 4)
        << "</td></tr>\n";
      if( nuclides.count( "Pu241" ) )
        results_html << "  <tr><td>Pu241</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu241_mass_frac, 4)
        << "</td></tr>\n";
      results_html << "  <tr><td>Pu242 (by corr)</td><td>"
      << SpecUtils::printCompact(100.0*pu_corr->pu242_mass_frac, 4)
      << "</td></tr>\n";
      if( nuclides.count( "Am241" ) )
        results_html << "  <tr><td>Am241</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->am241_mass_frac, 4)
        << "</td></tr>\n";
      results_html << "  </tbody>\n"
      << "</table>\n\n";
    }//if( pu_corr )
    
    
    double sum_rel_mass = 0.0;
    for( const auto &act : rel_activities )
      sum_rel_mass += act.rel_activity / act.nuclide->activityPerGram();
    
    
    results_html << "<table class=\"nuctable resulttable\">\n";
    results_html << "  <caption>Relative activities and mass fractions</caption>\n";
    results_html << "  <thead><tr>"
    " <th scope=\"col\">Nuclide</th>"
    " <th scope=\"col\">Rel. Act.</th>"
    " <th scope=\"col\">Total Mas Frac.</th>"
    " <th scope=\"col\">Enrichment</th>"
    " </tr></thead>\n";
    results_html << "  <tbody>\n";
    
    
    for( const auto &act : rel_activities )
    {
      const double rel_mass = act.rel_activity / act.nuclide->activityPerGram();
      
      results_html << "  <tr><td>" << act.nuclide->symbol;
      
      if( act.nuclide
         && (act.nuclide->atomicNumber == 94)
         && (act.nuclide->massNumber == 242)
         && (rel_eff.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable) )
      {
        results_html << " (by corr)";
      }
      
      results_html << "</td>"
      << "<td>" << act.rel_activity << " &plusmn; " << act.rel_activity_uncertainty << "</td>"
      << "<td>" << (100.0*rel_mass/sum_rel_mass) << "%</td>"
      << "<td>" << (100.0*mass_enrichment_fraction(act.nuclide,rel_eff_index)) << "%</td>"
      << "</tr>\n";
    }//for( const auto &act : rel_activities )
    
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
    for( size_t i = 1; i < rel_activities.size(); ++i )
    {
      for( size_t j = 0; j < i; ++j )
      {
        const auto nuc_i = rel_activities[i];
        const auto nuc_j = rel_activities[j];
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
    
    shared_ptr<const Material> self_atten;
    if( rel_eff.phys_model_self_atten )
      self_atten = rel_eff.phys_model_self_atten->material;
    vector<shared_ptr<const Material>> external_attens;
    for( const auto &ext : rel_eff.phys_model_external_atten )
      external_attens.push_back( ext->material );
    
    assert( !self_atten || (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) );
    assert( external_attens.empty() || (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) );
    
    for( const PeakDef &info : m_fit_peaks )
    {
      const double energy = info.mean();
      const SandiaDecay::Nuclide *nuc = info.parentNuclide();
      
      const NuclideRelAct *nuc_info = nullptr;
      for( const auto &rel_act : rel_activities )
        nuc_info = (rel_act.nuclide == nuc) ? &rel_act : nuc_info;
      
      assert( !nuc || nuc_info );
      if( nuc && !nuc_info )
        throw logic_error( "Peak with nuc not in solution?" );
      
      if( !nuc || !nuc_info )
      {
        snprintf( buffer, sizeof(buffer),
                 "  <tr><td>%.2f</td><td></td><td></td><td>%1.6G</td><td>%1.6G</td>"
                 "<td></td><td></td><td></td><td></td><td></td><td></td></tr>\n",
                 energy, info.amplitude(), info.amplitudeUncert() );
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
        double fit_rel_eff = -1.0;
        
        if( rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
        {
          fit_rel_eff = RelActCalc::eval_eqn( energy, rel_eff.rel_eff_eqn_type, m_rel_eff_coefficients[rel_eff_index] );
        }else
        {
          assert( m_drf );
          
          const RelActAutoCostFcn::PhysModelRelEqnDef input
          = RelActAutoCostFcn::make_phys_eqn_input( rel_eff, m_drf, m_rel_eff_coefficients[rel_eff_index], 0 );
          
          fit_rel_eff = RelActCalc::eval_physical_model_eqn( energy, input.self_atten,
                                                            input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
        }//if( m_options.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        
        //assert( fabs(meas_rel_eff - fit_rel_eff) < 0.1*std::max(meas_rel_eff,fit_rel_eff) );
        double fit_rel_eff_uncert = -1.0;
        
        if( m_rel_eff_covariance.size() )
        {
          try
          {
            //TODO: should implement `rel_eff_eqn_uncert(energy)` member function, and get rid of manually calling functions
            //fit_rel_eff_uncert = rel_eff_eqn_uncert( energy );
            //fit_rel_eff_uncert /= fit_rel_eff;
            
            vector<vector<double>> cov;
            if( m_rel_eff_covariance.size() > rel_eff_index )
              cov = m_rel_eff_covariance[rel_eff_index];
            const vector<double> &rel_eff_pars = m_rel_eff_coefficients[rel_eff_index];
            if( rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
            {
              fit_rel_eff_uncert = RelActCalc::eval_eqn_uncertainty( energy, rel_eff.rel_eff_eqn_type,
                                                                    rel_eff_pars, cov );
              
            }else
            {
              const RelActAutoCostFcn::PhysModelRelEqnDef input
              = RelActAutoCostFcn::make_phys_eqn_input( rel_eff, m_drf, rel_eff_pars, 0 );
              
              fit_rel_eff_uncert = RelActCalc::eval_physical_model_eqn_uncertainty( energy,
                                                                                   input.self_atten, input.external_attens, input.det.get(),
                                                                                   input.hoerl_b, input.hoerl_c, cov );
              
              fit_rel_eff_uncert /= fit_rel_eff;
            }
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
                 energy,
                 nuc->symbol.c_str(),
                 yield,
                 info.amplitude(),
                 info.amplitudeUncert() / info.amplitude(),
                 cps_over_yield,
                 fit_rel_eff,
                 100*fit_rel_eff_uncert,
                 Wt::WString::tr(PeakContinuum::offset_type_label_tr( info.continuum()->type() )).toUTF8().c_str(),
                 info.continuum()->lowerEnergy(),
                 info.continuum()->upperEnergy()
                 );
      }//if( !nuc || !nuc_info ) / else
      
      results_html << buffer;
    }//for( const RelEff::PeakInfo &info : input_peaks )
    
    results_html << "  </tbody>\n"
    << "</table>\n\n";
    
    bool ene_cal_fit = false;
    for( const bool &fit : m_fit_energy_cal )
      ene_cal_fit = (ene_cal_fit || fit);
    
    if( ene_cal_fit )
    {
      static_assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 3 );
      const char * const label[3] = { "offset adjustment", "gain adjustment", "cubic adjust" };
      const char * const post_script[3] = { "keV", "keV/chnl", "keV/chnl2" };
      
      const auto cal = m_spectrum ? m_spectrum->energy_calibration() : nullptr;
      const size_t num_chan = (cal ? cal->num_channels() : 1);
      
      results_html << "<div class=\"energycal\">\n";
      results_html << "Fit ";
      
      bool printed = false;
      for( size_t i = 0; i < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars; ++i )
      {
        if( !m_fit_energy_cal[i] )
          continue;
        
        double value = (m_energy_cal_adjustments[i]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0);
        if( i == 0 )
          value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
        else if( i == 1 )
          value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / num_chan;
        else if( i == 2 )
          value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / (num_chan * num_chan);
        
        results_html << (printed ? " and " : "") << label[i] << " of " << value
        << " " << post_script[i];
        
        printed = true;
      }//for( loop over RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars )
      
      results_html << ".\n";
      
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
    results_html << "  <tr><th scope=\"row\">nucs_of_el_same_age</th><td><code>" << rel_eff.nucs_of_el_same_age << "</code></td></tr>\n";
    results_html << "  <tr><th scope=\"row\">rel_eff_eqn_type</th><td><code>"
    << RelActCalc::to_str(rel_eff.rel_eff_eqn_type) << "</code></td></tr>\n";
    results_html << "  <tr><th scope=\"row\">rel_eff_eqn_order</th><td><code>" << rel_eff.rel_eff_eqn_order << "</code></td></tr>\n";
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
    
    string utc_time, local_time;
    if( wApp )
    {
      utc_time = Wt::WDateTime::currentDateTime().toString("yyyyMMdd hh:mm:ss").toUTF8();
      local_time = Wt::WLocalDateTime::currentDateTime().toString("yyyyMMdd hh:mm:ss").toUTF8();
    }else
    {
      const auto utc_ts = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
      utc_time = SpecUtils::to_common_string( utc_ts, true );
      
      
      std::time_t current_time = std::chrono::system_clock::to_time_t(utc_ts);
      struct tm current_local_time;
#if( defined( WIN32 ) )
      localtime_s( &current_local_time, &current_time );
#else
      localtime_r( &current_time, &current_local_time );
#endif
      char buffer[64] = { '\0' };
      
      //"9-Sep-2014 03:02:15 PM"
      std::strftime( buffer, sizeof(buffer), "%e-%b-%Y %r", &current_local_time );
      local_time = buffer;
      SpecUtils::trim(local_time); //"%e" can be precedded by a space
      
      //const auto local_ts = to_time_point( boost::posix_time::second_clock::local_time() );
      //local_time = SpecUtils::to_common_string( local_ts, true );
    }
    
    const double nsec_eval = 1.0E-6*m_num_microseconds_eval;
    results_html << "<div class=\"anacomputetime\">"
    << "Computation took " << PhysicalUnits::printToBestTimeUnits(nsec_eval)
    << " with " << m_num_function_eval_solution << " function calls to solve, and "
    << (m_num_function_eval_total - m_num_function_eval_solution)
    << " more for covariance\n"
    << "</div>\n";
    
    
    results_html << "<div class=\"anatime\">Analysis performed " << local_time << " (" << utc_time
    << " UTC), with code compiled " __TIMESTAMP__ "</div>\n";
    
 
    stringstream rel_eff_plot_values, add_rel_eff_plot_css;
    rel_eff_json_data( rel_eff_plot_values, add_rel_eff_plot_css, rel_eff_index );

    SpecUtils::ireplace_all( html, "${REL_EFF_DATA_VALS}", rel_eff_plot_values.str().c_str() );
    SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_ADDITIONAL_CSS}", add_rel_eff_plot_css.str().c_str() );
  
    string rel_eff_eqn_js;
    
    if( m_rel_eff_forms[rel_eff_index] != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      rel_eff_eqn_js = RelActCalc::rel_eff_eqn_js_function( m_rel_eff_forms[rel_eff_index], m_rel_eff_coefficients[rel_eff_index] );
    }else
    {
      const RelActAutoCostFcn::PhysModelRelEqnDef input
                   = RelActAutoCostFcn::make_phys_eqn_input( rel_eff, m_drf, m_rel_eff_coefficients[rel_eff_index], 0 );
      
      rel_eff_eqn_js = RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                                                                          input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
    }//if( not FramPhysicalModel ) / else
    
    SpecUtils::ireplace_all( html, "${FIT_REL_EFF_EQUATION}", rel_eff_eqn_js.c_str() );
    
    SpecUtils::ireplace_all( html, "${RESULTS_TXT}", results_html.str().c_str() );
  }
  
  
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
    for( const RoiRange &rr : m_options.rois )
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


void RelActAutoSolution::rel_eff_json_data( std::ostream &rel_eff_plot_values,
                                            std::ostream &add_rel_eff_plot_css,
                                           const size_t rel_eff_index ) const
{
  // TODO: refactor RelEffChart::setData to take care of all this functionality, and just call it their
  
  assert( rel_eff_index < m_rel_activities.size() );
  
  rel_eff_plot_values << "[";
  //rel_eff_plot_values << "{energy: 185, counts: 100, counts_uncert: 10, eff: 1.0, eff_uncert: 0.1, nuc_info: [{nuc: \"U235\", br: 0.2, rel_act: 300}]}"
  //", {energy: 1001, counts: 10, counts_uncert: 3, eff: 0.1, eff_uncert: 0.1, nuc_info: [{nuc: \"U238\", br: 0.1, rel_act: 100}]}";
  char buffer[1024] = { '\0' };
  const float live_time = m_spectrum ? m_spectrum->live_time() : 1.0f;
  
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
    for( const auto &rel_act : m_rel_activities[rel_eff_index] )
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
  
}//std::string rel_eff_json_data() const



double RelActAutoSolution::mass_enrichment_fraction( const SandiaDecay::Nuclide *nuclide, const size_t rel_eff_index ) const
{
  if( !nuclide )
    throw runtime_error( "RelActAutoSolution::mass_enrichment_fraction: null nuclide." );

  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "mass_enrichment_fraction: invalid rel eff index" );
  
  const vector<NuclideRelAct> &rel_acts = m_rel_activities[rel_eff_index];
  
  // We will consider this to by Pu if either Pu or Am241
  const bool is_pu = ((nuclide->atomicNumber == 94)
                      || ((nuclide->atomicNumber == 95) && (nuclide->massNumber == 241)));
  
  auto is_wanted_element = [is_pu,nuclide]( const SandiaDecay::Nuclide *test ) -> bool {
    if( test->atomicNumber == nuclide->atomicNumber )
      return true;
    
    return is_pu
            && ((test->atomicNumber == 94)
                || ((test->atomicNumber == 95) && (test->massNumber == 241)));
  };//is_wanted_element lamda
  
  const size_t nuc_index = nuclide_index( nuclide, rel_eff_index );
  const double rel_mass = rel_acts[nuc_index].rel_activity / nuclide->activityPerGram();
    
  double el_total_mass = 0.0;
  for( const NuclideRelAct &nuc : rel_acts )
  {
    if( is_wanted_element(nuc.nuclide) )
      el_total_mass += nuc.rel_activity / nuc.nuclide->activityPerGram();
  }
    
  if( is_pu && (m_corrected_pu.size() > rel_eff_index) && m_corrected_pu[rel_eff_index] )
  {
    const auto &pu_corr = m_corrected_pu[rel_eff_index];
    if( (nuclide->atomicNumber == 95) && (nuclide->massNumber == 241) )
      return pu_corr->am241_mass_frac;
    
    switch( nuclide->massNumber )
    {
      case 238: return pu_corr->pu238_mass_frac; break;
      case 239: return pu_corr->pu239_mass_frac; break;
      case 240: return pu_corr->pu240_mass_frac; break;
      case 241: return pu_corr->pu241_mass_frac; break;
      case 242: return pu_corr->pu242_mass_frac; break;
      default: return (1.0 - pu_corr->pu242_mass_frac) * rel_mass / el_total_mass;
    }//switch( nuclide->massNumber )
    
    assert( 0 );
  }//if( is_pu && pu_corr )
  
  return rel_mass / el_total_mass;
}//mass_enrichment_fraction


double RelActAutoSolution::mass_ratio( const SandiaDecay::Nuclide *numerator,
                                      const SandiaDecay::Nuclide *denominator,
                                      const size_t rel_eff_index ) const
{
  const double ratio = activity_ratio(numerator, denominator, rel_eff_index);
  
  return ratio * denominator->activityPerGram() / numerator->activityPerGram();
}//double mass_ratio(...)


double RelActAutoSolution::activity_ratio( const SandiaDecay::Nuclide *numerator,
                                          const SandiaDecay::Nuclide *denominator,
                                          const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "activity_ratio: invalid rel eff index" );
  
  return rel_activity(numerator, rel_eff_index) / rel_activity(denominator, rel_eff_index);
}//double activity_ratio(...)


double RelActAutoSolution::rel_activity( const SandiaDecay::Nuclide *nuc, const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "rel_activity: invalid rel eff index" );

  if( !nuc )
    throw runtime_error( "RelActAutoSolution::activity_ratio null input" );
    
  // If Pu242, we will use the corrected mass ratio to Pu239 to get activity relative
  //  to Pu239, and return the multiple of this.
  //  The other Pu isotopes dont need a correction, I dont think
  const bool is_pu = ((nuc->atomicNumber == 94)
                        || ((nuc->atomicNumber == 95) && (nuc->massNumber == 241)));
  if( (nuc->atomicNumber == 94) && (nuc->massNumber == 242)
     && (m_corrected_pu.size() > rel_eff_index) && m_corrected_pu[rel_eff_index] )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    assert( db );
    const SandiaDecay::Nuclide * const pu239 = db->nuclide("Pu239");
    const SandiaDecay::Nuclide * const pu242 = db->nuclide("Pu242");
    assert( pu239 && pu242 );
      
    const double corr_rel_pu239_act = pu239->activityPerGram() * m_corrected_pu[rel_eff_index]->pu239_mass_frac;
    const double corr_rel_pu242_act = pu242->activityPerGram() * m_corrected_pu[rel_eff_index]->pu242_mass_frac;
      
    const double raw_pu239_activity = m_rel_activities[rel_eff_index][nuclide_index(pu239,rel_eff_index)].rel_activity;
    return raw_pu239_activity * corr_rel_pu242_act / corr_rel_pu242_act;
  }//if( Pu242, and make correction )
    
  const size_t nuc_index = nuclide_index( nuc, rel_eff_index );

  return m_rel_activities[rel_eff_index][nuc_index].rel_activity;
}//double rel_activity(...)
  

double RelActAutoSolution::nuclide_counts( const SandiaDecay::Nuclide *nuclide, const size_t rel_eff_index ) const
{
  double total_counts = 0.0;
  
  const double live_time = m_spectrum ? static_cast<double>(m_spectrum->live_time()) : 0.0;

  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "nuclide_counts: invalid rel eff index" ); 
    

  const vector<NuclideRelAct> &nuclides = m_rel_activities[rel_eff_index];
  vector<unique_ptr<vector<NucInputGamma::EnergyYield>>> aged_gammas_cache( nuclides.size() );  

  for( const RoiRange &range : m_final_roi_ranges )
  {
      
      for( size_t nuc_index = 0; nuc_index < nuclides.size(); ++nuc_index )
      {
        const NuclideRelAct &nucinfo = nuclides[nuc_index];
        if( nucinfo.nuclide != nuclide )
          continue;

        const double rel_act = rel_activity( nucinfo.nuclide, rel_eff_index ); //almost same as `nucinfo.rel_activity`, except with Pu corrections
        
        const vector<pair<double,double>> &gamma_energy_br = nucinfo.gamma_energy_br;

        
        for( const pair<double,double> &gamma : gamma_energy_br )
        {
          const double energy = gamma.first;
          const double yield = gamma.second;
          
          // Filter the energy range based on true energies (i.e., not adjusted energies), and skip zero-yield gammas
          if( (yield < std::numeric_limits<float>::min())
            || (energy < range.lower_energy) 
            || (energy > range.upper_energy) )
          {
            continue;
          }
          
          const double rel_eff = relative_efficiency( energy, rel_eff_index );
          if( isinf(rel_eff) || isnan(rel_eff) )
            throw runtime_error( "nuclide_counts: inf or NaN rel. eff for "
                                + std::to_string(energy) + " keV."  );
        
          double br_uncert_adj = 1.0;
          if( (m_options.additional_br_uncert > 0.0) && !m_peak_ranges_with_uncert.empty() )
          {
            assert( m_add_br_uncert_start_index != std::numeric_limits<size_t>::max() );
            assert( m_add_br_uncert_start_index < m_final_parameters.size() );
            assert( (m_add_br_uncert_start_index + m_peak_ranges_with_uncert.size()) == m_final_parameters.size() );
            
            // TODO: we could use std::lower_bound to find the applicable `m_peak_ranges_with_uncert` element
            for( size_t range_index = 0; range_index < m_peak_ranges_with_uncert.size(); ++range_index )
            {
              const pair<double,double> &range = m_peak_ranges_with_uncert[range_index];
              
              if( (energy >= range.first) && (energy <= range.second) )
              {
                const double adjvalue = m_final_parameters[m_add_br_uncert_start_index + range_index]/RelActAutoCostFcn::sm_peak_range_uncert_par_offset - 1.0;
                br_uncert_adj = max(0.0, 1.0 + adjvalue);
                
                break;
              }//
            }//for( find range this gamma belongs to, if any )
          }//if( m_options.additional_br_uncert > 0.0 && !m_peak_ranges_with_uncert.empty() )
          
          
          const double peak_amplitude = rel_act * live_time * rel_eff * yield * br_uncert_adj;
          if( isinf(peak_amplitude) || isnan(peak_amplitude) )
            throw runtime_error( "nuclide_counts: inf or NaN peak amplitude for "
                                + std::to_string(energy) + " keV.");
          
          if( peak_amplitude > static_cast<double>( std::numeric_limits<float>::min() ) )
            total_counts += peak_amplitude;
        }//for( const SandiaDecay::EnergyRatePair &gamma : gammas )
      }//for( const NucInputGamma &nucinfo : m_nuclides )
  }//for( const RoiRange &range : m_final_roi_ranges )
  
   return total_counts;
}//double nuclide_counts(...)


double RelActAutoSolution::relative_efficiency( const double energy, const size_t rel_eff_index ) const
{
  assert( m_rel_eff_forms.size() == m_rel_eff_coefficients.size() );
  
  assert( rel_eff_index < m_rel_eff_forms.size() );
  if( rel_eff_index >= m_rel_eff_forms.size() )
    throw std::logic_error( "relative_efficiency: invalid rel eff index" ); 
    
  const RelActCalc::RelEffEqnForm eqn_form = m_rel_eff_forms[rel_eff_index];
  const vector<double> &coeffs = m_rel_eff_coefficients[rel_eff_index];
  
  if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::eval_eqn( energy, eqn_form, coeffs );
    
  if( rel_eff_index >= m_phys_model_results.size() )
    throw std::logic_error( "relative_efficiency: invalid rel eff index" );
  
  const optional<PhysicalModelFitInfo> &phys_model_result = m_phys_model_results[rel_eff_index];
  
  if( !phys_model_result.has_value() )
    throw std::logic_error( "relative_efficiency: no Physical Model results" );
    
  optional<RelActCalc::PhysModelShield<double>> self_atten;
  vector<RelActCalc::PhysModelShield<double>> external_attens;

  if( phys_model_result->self_atten.has_value() )
  {
    RelActCalc::PhysModelShield<double> self_shield;
    self_shield.material = phys_model_result->self_atten->material;
    if( phys_model_result->self_atten->atomic_number.has_value() )
      self_shield.atomic_number = *phys_model_result->self_atten->atomic_number;
    self_shield.areal_density = phys_model_result->self_atten->areal_density;
    self_atten = std::move(self_shield);
  }
  
  for( const PhysicalModelFitInfo::ShieldInfo &ext_shield : phys_model_result->ext_shields )
  {
    RelActCalc::PhysModelShield<double> ext_shield_mdl;
    ext_shield_mdl.material = ext_shield.material;
    if( ext_shield.atomic_number.has_value() ) 
      ext_shield_mdl.atomic_number = *ext_shield.atomic_number;
    ext_shield_mdl.areal_density = ext_shield.areal_density;
    external_attens.push_back( std::move(ext_shield_mdl) );
  }

  return RelActCalc::eval_physical_model_eqn( energy, self_atten, external_attens, m_drf.get(), 
                                              phys_model_result->hoerl_b, phys_model_result->hoerl_c );
}//double relative_efficiency( const double energy, const size_t rel_eff_index ) const;


size_t RelActAutoSolution::nuclide_index( const SandiaDecay::Nuclide *nuclide, const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "nuclide_index: invalid rel eff index" );
  
  for( size_t i = 0; i < m_rel_activities[rel_eff_index].size(); ++i )
    if( nuclide == m_rel_activities[rel_eff_index][i].nuclide )
      return i;
  
  assert( 0 );
  throw runtime_error( "RelActAutoSolution: " + (nuclide ? nuclide->symbol : string()) + " is an invalid nuclide." );
  return m_rel_activities.size();
}//nuclide_index(...)


string RelActAutoSolution::rel_eff_eqn_js_function( const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_rel_eff_forms.size() );
  assert( m_rel_eff_coefficients.size() == m_rel_eff_forms.size() );
  
  assert( rel_eff_index < m_rel_eff_forms.size() );
  if( rel_eff_index >= m_rel_eff_forms.size() )
    throw std::logic_error( "rel_eff_eqn_js_function: invalid rel eff index" );
  
  const RelActCalc::RelEffEqnForm eqn_form = m_rel_eff_forms[rel_eff_index];
  const vector<double> &coeffs = m_rel_eff_coefficients[rel_eff_index];
  
  if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::rel_eff_eqn_js_function( eqn_form, coeffs );
  
  const RelActAutoCostFcn::PhysModelRelEqnDef input = RelActAutoCostFcn::make_phys_eqn_input(
                            m_options.rel_eff_curves[rel_eff_index], m_drf, coeffs, 0 );
  
    
  return RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                            input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
}//string rel_eff_eqn_js_function() const
  
  
std::shared_ptr<SpecUtils::EnergyCalibration> RelActAutoSolution::get_adjusted_energy_cal() const
{
  shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = m_foreground ? m_foreground->energy_calibration() : nullptr;
  if( !orig_cal || !orig_cal->valid() )
    return nullptr;
  
  auto new_cal = make_shared<SpecUtils::EnergyCalibration>( *orig_cal );
  
  if( !m_options.fit_energy_cal )
    return new_cal;
  
  assert( m_energy_cal_adjustments.size() == RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars );
  if( m_energy_cal_adjustments.size() != RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars )
    throw runtime_error( "RelActAutoSolution::get_adjusted_energy_cal: m_energy_cal_adjustments empty" );
  
  const size_t num_channel = orig_cal->num_channels();
  const SpecUtils::EnergyCalType energy_cal_type = orig_cal->type();
  
    
  const double offset_adj = (m_energy_cal_adjustments[0]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                              * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
  const double gain_adj = (m_energy_cal_adjustments[1]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                              * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
  const double quad_adj = (RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars > 2)
                      ? ((m_energy_cal_adjustments[2]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                              * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple)
                      : 0.0;
    
  switch( energy_cal_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    {
      vector<float> coefs = orig_cal->coefficients();
      assert( coefs.size() >= 2 );
      coefs[0] -= offset_adj;
      coefs[1] -= (gain_adj / num_channel);
        
      if( quad_adj != 0.0 )
      {
        if( coefs.size() > 2 )
          coefs[2] -= (quad_adj / (num_channel*num_channel));
        else
          coefs.push_back( -quad_adj / (num_channel*num_channel) );
      }//if( quad_adj != 0.0 )
        
      const auto &dev_pairs = orig_cal->deviation_pairs();
      new_cal->set_polynomial( num_channel, coefs, dev_pairs );
      break;
    }//case polynomial
        
    case SpecUtils::EnergyCalType::FullRangeFraction:
    {
      vector<float> coefs = orig_cal->coefficients();
      assert( coefs.size() >= 2 );
      coefs[0] -= offset_adj;
      coefs[1] -= gain_adj;
        
      if( quad_adj != 0.0 )
      {
        if( coefs.size() > 2 )
          coefs[2] -= quad_adj;
        else
          coefs.push_back( -quad_adj );
      }//if( quad_adj != 0.0 )
      
      const auto &dev_pairs = orig_cal->deviation_pairs();
      new_cal->set_full_range_fraction( num_channel, coefs, dev_pairs );
      break;
    }//case FullRangeFraction:
        
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    {
      assert( orig_cal->channel_energies() );
      if( !orig_cal->channel_energies() || orig_cal->channel_energies()->empty() )
        throw runtime_error( "Invalid lower channel energies???" );
        
      assert( quad_adj == 0.0 );
      const double lower_energy = orig_cal->lower_energy();
      const double range = orig_cal->upper_energy() - lower_energy;
        
      vector<float> lower_energies = *orig_cal->channel_energies();
      for( size_t i = 0; i < lower_energies.size(); ++i )
      {
        const float low_e = lower_energies[i];
        lower_energies[i] = (low_e - offset_adj + gain_adj*lower_energy/range)/(1 + gain_adj/range);
      }
        
      new_cal->set_lower_channel_energy( num_channel, std::move(lower_energies) );
      break;
    }//case LowerChannelEdge
        
    case SpecUtils::EnergyCalType::InvalidEquationType:
      break;
  }//switch( m_energy_cal->type() )
  
  return new_cal;
}//std::shared_ptr<SpecUtils::EnergyCalibration> get_adjusted_energy_cal() const
  
  
  
RelActAutoSolution solve( const Options options,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         std::shared_ptr<const DetectorPeakResponse> input_drf,
                         std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                         std::shared_ptr<std::atomic_bool> cancel_calc
                         )
{ 
  const std::vector<RoiRange> &energy_ranges = options.rois;
  const std::vector<FloatingPeak> &extra_peaks = options.floating_peaks;
  

  const RelActAutoSolution orig_sol = RelActAutoCostFcn::solve_ceres(
                     options,
                     foreground,
                     background,
                     input_drf,
                     all_peaks,
                     cancel_calc );
  
  bool all_roi_full_range = true;
  for( const auto &roi : energy_ranges )
    all_roi_full_range = (all_roi_full_range && roi.force_full_range && !roi.allow_expand_for_peak_width);
  
  if( all_roi_full_range
     || (orig_sol.m_status != RelActAutoSolution::Status::Success)
     || !orig_sol.m_spectrum )
  {
    return orig_sol;
  }
  
  // If we are here there was at least one ROI that didnt have force_full_range set, or had
  //  allow_expand_for_peak_width set.
  // So we will go through and adjust these ROIs based on peaks that are statistically significant,
  //  based on initial solution, and then re-fit.
    
  RelActAutoSolution current_sol = orig_sol;
  
  int num_function_eval_solution = orig_sol.m_num_function_eval_solution;
  int num_function_eval_total = orig_sol.m_num_function_eval_total;
  int num_microseconds_eval = orig_sol.m_num_microseconds_eval;
  int num_microseconds_in_eval = orig_sol.m_num_microseconds_in_eval;
  
  
  bool stop_iterating = false, errored_out_of_iterating = false;
  const size_t max_roi_adjust_iterations = 2;
  size_t num_roi_iters = 0;
  for( ; !stop_iterating && (num_roi_iters < max_roi_adjust_iterations); ++num_roi_iters )
  {
    assert( current_sol.m_spectrum
           && current_sol.m_spectrum->energy_calibration()
           && current_sol.m_spectrum->energy_calibration()->valid() );
    
    if( !current_sol.m_spectrum
       || !current_sol.m_spectrum->energy_calibration()
       || !current_sol.m_spectrum->energy_calibration()->valid() )
    {
      return orig_sol;
    }
    
    vector<RoiRange> fixed_energy_ranges;
    for( const RoiRange &roi : energy_ranges )
    {
      if( roi.force_full_range )
        fixed_energy_ranges.push_back( roi );
    }//for( const RoiRange &roi : energy_ranges )
    
    
    /** `orig_sol.m_spectrum` is background subtracted, and energy corrected, of those options were
     selected.
     */
    const auto spectrum = current_sol.m_spectrum;
    const float live_time = spectrum->live_time();
    const double num_sigma_half_roi = DEFAULT_PEAK_HALF_WIDTH_SIGMA;
    
    // We'll collect all the individual
    vector<RoiRange> significant_peak_ranges;
    
    /** Updates the passed in ROI to have limits for the peak mean energy passed in, and adds ROI
     to `significant_peak_ranges`.
     */
    auto add_updated_roi = [&significant_peak_ranges,
                             &current_sol,
                             &fixed_energy_ranges,
                             num_sigma_half_roi]( RoiRange roi, const double energy ){
      const double fwhm = eval_fwhm( energy, current_sol.m_fwhm_form, current_sol.m_fwhm_coefficients.data(), current_sol.m_fwhm_coefficients.size() );
      const double sigma = fwhm / 2.35482f;
      
      double roi_lower = energy - (num_sigma_half_roi * sigma);
      double roi_upper = energy + (num_sigma_half_roi * sigma);
      
      bool keep_roi = true;
      if( !roi.allow_expand_for_peak_width )
      {
        roi_lower = std::max( roi_lower, roi.lower_energy );
        roi_upper = std::min( roi_upper, roi.upper_energy );
      }else
      {
        // Make sure we havent expanded into the range of any ROIs with forced widths
        for( const RoiRange &fixed : fixed_energy_ranges )
        {
          const bool overlaps = ((roi.upper_energy > fixed.lower_energy)
                                 && (roi.lower_energy < fixed.upper_energy));
          if( overlaps )
          {
            // Do some development checks about bounds of ROI
            if( (roi_lower >= fixed.lower_energy) && (roi_upper <= fixed.upper_energy) )
            {
              // This ROI is completely within a fixed ROI - this shouldnt happen!
              assert( 0 );
              keep_roi = false;
              break;
            }//if( this peaks ROI is completely within a fixed ROI )
            
            if( (fixed.lower_energy >= roi_lower) && (fixed.upper_energy <= roi_upper) )
            {
              // This fixed ROI is completely within a the peaks ROI - this shouldnt happen!
              assert( 0 );
              keep_roi = false;
              break;
            }//if( the fixed ROI is completely within this peaks ROI )
            
            if( roi.upper_energy > fixed.upper_energy )
            {
              assert( roi_lower <= fixed.upper_energy );
              assert( roi_upper > fixed.upper_energy );
              roi_lower = std::max( roi_lower, fixed.upper_energy );
            }else
            {
              assert( roi_upper >= fixed.lower_energy );
              assert( roi_lower < fixed.lower_energy );
              roi_upper = std::min( roi_upper, fixed.lower_energy );
            }
            
            assert( roi_lower < roi_upper );
          }//if( overlaps )
        }//for( const RoiRange &fixed : fixed_energy_ranges )
      }//if( !roi.allow_expand_for_peak_width ) / else
      
      
      if( !keep_roi || (roi_lower >= roi_upper) )
      {
        // SHouldnt ever get here
        assert( 0 );
        return;
      }
      
      roi.force_full_range = true;
      roi.allow_expand_for_peak_width = false;
      roi.lower_energy = roi_lower;
      roi.upper_energy = roi_upper;
      
      significant_peak_ranges.push_back( roi );
    };//add_updated_roi
    
    
    vector<optional<RelActAutoCostFcn::PhysModelRelEqnDef<double>>> phys_model_inputs( current_sol.m_rel_eff_coefficients.size() );
    assert( options.rel_eff_curves.size() == current_sol.m_rel_activities.size() );

    
    // Note: we loop over original energy_ranges, not the energy ranges from the solution,
    //       (to avoid the ROIs from expanding continuously, and also we've marked the
    //        updated ROIs as force full-range)
    for( const RoiRange &roi : energy_ranges )
    {
      if( roi.force_full_range )
        continue;
      
      // TODO: should we group peaks together by nuclide?  I think so, but probably not a huge effect at first
      
      // Estimate peak amplitudes, so we can decide if they are significant enough.
      for( size_t rel_eff_index = 0; rel_eff_index < current_sol.m_rel_activities.size(); ++rel_eff_index )
      {
        const vector<NuclideRelAct> &rel_acts = current_sol.m_rel_activities[rel_eff_index];
        const RelActCalcAuto::RelEffCurveInput &rel_eff = options.rel_eff_curves[rel_eff_index];
        const vector<double> &rel_eff_coefs = current_sol.m_rel_eff_coefficients[rel_eff_index];
        const RelActCalc::RelEffEqnForm eqn_form = rel_eff.rel_eff_eqn_type;
        
        optional<RelActAutoCostFcn::PhysModelRelEqnDef<double>> &phys_model_input = phys_model_inputs[rel_eff_index];
        if( (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel) && !phys_model_input.has_value() )
          *phys_model_input = RelActAutoCostFcn::make_phys_eqn_input( rel_eff, current_sol.m_drf, rel_eff_coefs, 0 );
        
        
        for( const NuclideRelAct &rel_act : rel_acts )
        {
          assert( rel_act.nuclide );
          if( !rel_act.nuclide )
            continue;
          
          if( (rel_eff.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable)
             && (rel_act.nuclide->atomicNumber == 94)
             && (rel_act.nuclide->massNumber == 242) )
          {
            continue;
          }
          
          
          for( const pair<double,double> &energy_br : rel_act.gamma_energy_br )
          {
            const double &energy = energy_br.first;
            const double &br = energy_br.second;
            if( (energy < roi.lower_energy) || (energy > roi.upper_energy) )
              continue;
            
            double rel_eff_val = -1.0;
            
            if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
            {
              rel_eff_val = RelActCalc::eval_eqn( energy, eqn_form, rel_eff_coefs );
            }else
            {
              assert( phys_model_input.has_value() );
              rel_eff_val = RelActCalc::eval_physical_model_eqn( energy, phys_model_input->self_atten,
                                                                phys_model_input->external_attens, phys_model_input->det.get(),
                                                                phys_model_input->hoerl_b, phys_model_input->hoerl_c );
            }//if( options.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
            
            const double expected_counts = live_time * br * rel_eff_val * rel_act.rel_activity;
            const double fwhm = eval_fwhm( energy, current_sol.m_fwhm_form,
                                          current_sol.m_fwhm_coefficients.data(),
                                          current_sol.m_fwhm_coefficients.size() );
            const double sigma = fwhm / 2.35482f;
            
            const double peak_width_nsigma = 3.0;
            const float lower_energy = static_cast<float>(energy - (peak_width_nsigma * sigma));
            const float upper_energy = static_cast<float>(energy + (peak_width_nsigma * sigma));
            const double data_counts = spectrum->gamma_integral(lower_energy, upper_energy);
            
            // We'll use a very simple requirement that the expected peak area should be at least 3
            //  times the sqrt of the data area for the peak; this is of course just a coarse FOM
            //  TODO: the value of 3.0 was chosen "by eye" from only a couple example spectra - need to re-evaluate
            const double significance_limit = 3.0;
            const bool significant = (expected_counts > significance_limit*sqrt(data_counts));
            
            cout << "For " << energy << " keV " << rel_act.nuclide->symbol
            << " expect counts=" << expected_counts
            << ", and FWHM=" << fwhm
            << ", with data_counts=" << data_counts
            << ", is_significant=" << significant
            << endl;
            
            if( !significant )
              continue;
            
            add_updated_roi( roi, energy );
          }//for( const pair<double,double> &energy_br : rel_act.gamma_energy_br )
        }//for( const NuclideRelAct &rel_act : current_sol.m_rel_activities )
      }//for( const vector<NuclideRelAct> &rel_acts : current_sol.m_rel_activities )
    }//for( const RoiRange &roi : energy_ranges )
    
    // Now make sure all `extra_peaks` are in an energy range; if not add an energy range for them
    for( const FloatingPeak &peak : extra_peaks )
    {
      bool in_energy_range = false;
      for( const RoiRange &roi : significant_peak_ranges )
      {
        in_energy_range = ((peak.energy >= roi.lower_energy) && (peak.energy <= roi.upper_energy));
        if( in_energy_range )
        {
          add_updated_roi( roi, peak.energy );
          break;
        }
      }//for( const RoiRange &roi : significant_peak_ranges )
      
      if( !in_energy_range )
      {
        for( const RoiRange &roi : fixed_energy_ranges )
        {
          in_energy_range = ((peak.energy >= roi.lower_energy) && (peak.energy <= roi.upper_energy));
          if( in_energy_range )
            break;
        }//for( const RoiRange &roi : fixed_energy_ranges )
      }//if( we needed to check if the floating peak is in a fixed range )
      
      //If the fixed peak isnt in any of the ROI's from significant gammas, or forced-full-range
      //  ROIs, create a ROI for it
      if( !in_energy_range )
      {
        bool found_roi = false;
        for( const RoiRange &roi : energy_ranges )
        {
          found_roi = ((peak.energy >= roi.lower_energy) && (peak.energy <= roi.upper_energy));
          if( found_roi )
          {
            add_updated_roi( roi, peak.energy );
            break;
          }//
        }//for( const RoiRange &roi : energy_ranges )
        
        assert( found_roi );
      }//if( the fixed-peak isnt in any of the ROI's
    }//for( const FloatingPeak &peak : extra_peaks )
    
    //  Sort ROIs first so they are in increasing energy
    sort_rois_by_energy( significant_peak_ranges );
    
    // Now combine overlapping ranges in significant_peak_ranges
    vector<RoiRange> combined_sig_peak_ranges;
    for( size_t index = 0; index < significant_peak_ranges.size(); ++index )
    {
      RoiRange range = significant_peak_ranges[index];
      
      // We'll loop over remaining ROIs until there is not an overlap; note that
      //  `index` is incremented if we combine with a ROI.
      for( size_t j = index + 1; j < significant_peak_ranges.size(); ++j, ++index )
      {
        const RoiRange &next_range = significant_peak_ranges[j];
        
        if( next_range.lower_energy > range.upper_energy )
          break; // note: doesnt increment `index`
        
        range.upper_energy = next_range.upper_energy;
        const int cont_type_int = std::max( static_cast<int>(range.continuum_type),
                                           static_cast<int>(next_range.continuum_type) );
        range.continuum_type = PeakContinuum::OffsetType(cont_type_int);
      }//for( loop over the next peaks, until they shouldnt be combined )
      
      combined_sig_peak_ranges.push_back( range );
    }//for( size_t i = 0; i < significant_peak_ranges.size(); ++i )
    
    // Put all the ROIs into one vector, and sort them
    vector<RoiRange> updated_energy_ranges = fixed_energy_ranges;
    
    updated_energy_ranges.insert( end(updated_energy_ranges),
                                 begin(combined_sig_peak_ranges),
                                 end(combined_sig_peak_ranges) );
    
    sort_rois_by_energy( updated_energy_ranges );
    
    cout << "\nFor iteration " << num_roi_iters << ", the energy ranges being are going from the"
    " original:\n\t{";
    for( size_t i = 0; i < energy_ranges.size(); ++i )
      cout << (i ? ", " : "") << "{" << energy_ranges[i].lower_energy
           << ", " << energy_ranges[i].upper_energy << "}";
    cout << "}\nTo:\n\t{";
    for( size_t i = 0; i < updated_energy_ranges.size(); ++i )
      cout << (i ? ", " : "") << "{" << updated_energy_ranges[i].lower_energy
      << ", " << updated_energy_ranges[i].upper_energy << "}";
    cout << "}\n\n";
    
    
    try
    {
      auto updated_options = options;
      updated_options.rois = updated_energy_ranges;
      
      const RelActAutoSolution updated_sol
      = RelActAutoCostFcn::solve_ceres( updated_options, foreground, background, input_drf, all_peaks, cancel_calc );
      
      switch( updated_sol.m_status )
      {
        case RelActAutoSolution::Status::Success:
          break;
          
        case RelActAutoSolution::Status::NotInitiated:
          throw runtime_error( "After breaking up energy ranges, could not initialize finding solution." );
          
        case RelActAutoSolution::Status::FailedToSetupProblem:
          throw runtime_error( "After breaking up energy ranges, the setup for finding a solution became invalid." );
          
        case RelActAutoSolution::Status::FailToSolveProblem:
          throw runtime_error( "After breaking up energy ranges, failed to solve the problem." );
          
        case RelActAutoSolution::Status::UserCanceled:
          return updated_sol;
      }//switch( updated_sol.m_status )
      
      // If the number of ROIs didnt change for this iteration, we'll assume we're done - this is
      //  totally ignoring that the actual ranges of the ROIs could have significantly changed, and
      //  could change a little more
      //  TODO: evaluate if ROI ranges have changed enough to warrant another iteration of finding a solution
      const size_t num_prev_rois = current_sol.m_final_roi_ranges.size();
      stop_iterating = (updated_sol.m_final_roi_ranges.size() == num_prev_rois);
      
      current_sol = updated_sol;
      
      // Update tally of function calls and eval time (eval time will be slightly off, but oh well)
      num_function_eval_solution += current_sol.m_num_function_eval_solution;
      num_function_eval_total += current_sol.m_num_function_eval_total;
      num_microseconds_eval += current_sol.m_num_microseconds_eval;
      num_microseconds_in_eval += current_sol.m_num_microseconds_in_eval;
    }catch( std::exception &e )
    {
      stop_iterating = errored_out_of_iterating = true;
      current_sol.m_warnings.push_back( "Failed to break up energy ranges to ROIs with significant"
                                    " gamma counts: " + string(e.what())
                                    + " - the current solution is from before final ROI divisions" );
    }//try / catch
  }//for( size_t roi_iteration = 0; roi_iteration < max_roi_adjust_iterations; ++roi_iteration )
  
  
  current_sol.m_num_function_eval_solution = num_function_eval_solution;
  current_sol.m_num_function_eval_total = num_function_eval_total;
  current_sol.m_num_microseconds_eval = num_microseconds_eval;
  current_sol.m_num_microseconds_in_eval = num_microseconds_in_eval;
  
  
  if( !errored_out_of_iterating && !stop_iterating )
    current_sol.m_warnings.push_back( "Final ROIs based on gamma line significances may not have been found." );
  
  cout << "It took " << num_roi_iters << " iterations (of max " << max_roi_adjust_iterations
  << ") to find final ROIs to use;";
  if( errored_out_of_iterating )
    cout << " an error prevented finding final ROIs." << endl;
  else if( stop_iterating )
    cout << " final ROIs were found." << endl;
  else
    cout << " final ROIs were NOT found." << endl;
  
  return current_sol;
}//RelActAutoSolution

#if( PERFORM_DEVELOPER_CHECKS )

void RoiRange::equalEnough( const RoiRange &lhs, const RoiRange &rhs )
{
  if( fabs(lhs.lower_energy - rhs.lower_energy) > 1.0e-3 )
    throw std::runtime_error( "Lower energy in lhs and rhs are not the same" );
  
  if( fabs(lhs.upper_energy - rhs.upper_energy) > 1.0e-3 )
    throw std::runtime_error( "Upper energy in lhs and rhs are not the same" );
  
  if( lhs.continuum_type != rhs.continuum_type )
    throw std::runtime_error( "Continuum type in lhs and rhs are not the same" );
  
  if( lhs.force_full_range != rhs.force_full_range )
    throw std::runtime_error( "Force full range in lhs and rhs are not the same" );
  
  if( lhs.allow_expand_for_peak_width != rhs.allow_expand_for_peak_width )
    throw std::runtime_error( "Allow expand for peak width in lhs and rhs are not the same" );
}//RoiRange::equalEnough


void NucInputInfo::equalEnough( const NucInputInfo &lhs, const NucInputInfo &rhs )
{
  if( lhs.nuclide != rhs.nuclide )
    throw std::runtime_error( "Nuclide in lhs and rhs are not the same" );
  
  if( fabs(lhs.age - rhs.age) > 1.0e-5 * std::max(lhs.age, rhs.age) )
    throw std::runtime_error( "Age in lhs and rhs are not the same" );
  
  if( lhs.fit_age != rhs.fit_age )
    throw std::runtime_error( "Fit age in lhs and rhs are not the same" );
  
  if( lhs.gammas_to_exclude.size() != rhs.gammas_to_exclude.size() )
    throw std::runtime_error( "Number of gammas to exclude in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.gammas_to_exclude.size(); ++i )
    if( fabs(lhs.gammas_to_exclude[i] - rhs.gammas_to_exclude[i]) > 1.0e-3 )
      throw std::runtime_error( "Gamma to exclude in lhs and rhs are not the same" );
  
  if( lhs.peak_color_css != rhs.peak_color_css )
    throw std::runtime_error( "Peak color CSS in lhs and rhs are not the same" );
}//NucInputInfo::equalEnough


void FloatingPeak::equalEnough( const FloatingPeak &lhs, const FloatingPeak &rhs )
{
  if( fabs(lhs.energy - rhs.energy) > 1.0e-5 )
    throw std::runtime_error( "Energy in lhs and rhs are not the same" );
  
  if( lhs.release_fwhm != rhs.release_fwhm )
    throw std::runtime_error( "Release FWHM in lhs and rhs are not the same" );
  
  if( lhs.apply_energy_cal_correction != rhs.apply_energy_cal_correction )
    throw std::runtime_error( "Apply energy cal correction in lhs and rhs are not the same" );
}//FloatingPeak::equalEnough


void RelEffCurveInput::equalEnough( const RelEffCurveInput &lhs, const RelEffCurveInput &rhs )
{
  if( lhs.nuclides.size() != rhs.nuclides.size() )
    throw std::runtime_error( "Number of nuclides in lhs and rhs are not the same" );
  
  if( lhs.nucs_of_el_same_age != rhs.nucs_of_el_same_age )
    throw std::runtime_error( "Nuclides of same element same age in lhs and rhs are not the same" );
  
  if( lhs.rel_eff_eqn_type != rhs.rel_eff_eqn_type )
    throw std::runtime_error( "Relative efficiency equation type in lhs and rhs are not the same" );
  
  if( lhs.rel_eff_eqn_order != rhs.rel_eff_eqn_order )
    throw std::runtime_error( "Relative efficiency equation order in lhs and rhs are not the same" );
  
  if( !lhs.phys_model_self_atten != !rhs.phys_model_self_atten )
    throw std::runtime_error( "Physical model self attenuation in lhs and rhs are not the same" );
  
  if( lhs.phys_model_self_atten && rhs.phys_model_self_atten )
    RelActCalc::PhysicalModelShieldInput::equalEnough( *lhs.phys_model_self_atten, *rhs.phys_model_self_atten );
  
  if( lhs.phys_model_external_atten.size() != rhs.phys_model_external_atten.size() )
    throw std::runtime_error( "Number of physical model external attenuations in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.phys_model_external_atten.size(); ++i )
  {
    assert( lhs.phys_model_external_atten[i] );
    assert( rhs.phys_model_external_atten[i] );

    if( !lhs.phys_model_external_atten[i] || !rhs.phys_model_external_atten[i] )
      throw std::runtime_error( "Physical model external attenuation has nullptr" );
    
    RelActCalc::PhysicalModelShieldInput::equalEnough( *lhs.phys_model_external_atten[i], *rhs.phys_model_external_atten[i] );
  }
  
  if( lhs.phys_model_use_hoerl != rhs.phys_model_use_hoerl )
    throw std::runtime_error( "Physical model use Hoerl in lhs and rhs are not the same" );
  
  if( lhs.pu242_correlation_method != rhs.pu242_correlation_method )
    throw std::runtime_error( "Pu242 correlation method in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.nuclides.size(); ++i )
    RelActCalcAuto::NucInputInfo::equalEnough( lhs.nuclides[i], rhs.nuclides[i] );
}//RelEffCurveInput::equalEnough


void Options::equalEnough( const Options &lhs, const Options &rhs )
{
  if( lhs.fit_energy_cal != rhs.fit_energy_cal )
    throw std::runtime_error( "Fit energy calibration in lhs and rhs are not the same" );
  
  if( lhs.fwhm_form != rhs.fwhm_form )
    throw std::runtime_error( "FWHM form in lhs and rhs are not the same" );
  
  if( lhs.spectrum_title != rhs.spectrum_title )
    throw std::runtime_error( "Spectrum title in lhs and rhs are not the same" );
  
  if( lhs.skew_type != rhs.skew_type )
    throw std::runtime_error( "Skew type in lhs and rhs are not the same" );
  
  if( fabs(lhs.additional_br_uncert - rhs.additional_br_uncert) > 1.0e-5
    && ((lhs.additional_br_uncert > 0.0) || (rhs.additional_br_uncert > 0.0)) )
    throw std::runtime_error( "Additional BR uncertanty in lhs and rhs are not the same" );
  
  if( lhs.rel_eff_curves.size() != rhs.rel_eff_curves.size() )
    throw std::runtime_error( "Number of relative efficiency curves in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.rel_eff_curves.size(); ++i )
    RelActCalcAuto::RelEffCurveInput::equalEnough( lhs.rel_eff_curves[i], rhs.rel_eff_curves[i] );


  if( lhs.rois.size() != rhs.rois.size() )
    throw std::runtime_error( "Number of ROIs in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.rois.size(); ++i )
    RelActCalcAuto::RoiRange::equalEnough( lhs.rois[i], rhs.rois[i] );

  if( lhs.floating_peaks.size() != rhs.floating_peaks.size() )
    throw std::runtime_error( "Number of floating peaks in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.floating_peaks.size(); ++i )
    RelActCalcAuto::FloatingPeak::equalEnough( lhs.floating_peaks[i], rhs.floating_peaks[i] );
}//Options::equalEnough


void RelActAutoGuiState::equalEnough( const RelActAutoGuiState &lhs, const RelActAutoGuiState &rhs )
{
  RelActCalcAuto::Options::equalEnough( lhs.options, rhs.options );

  if( lhs.background_subtract != rhs.background_subtract )
    throw std::runtime_error( "Background subtract in lhs and rhs are not the same" );
  
  if( lhs.show_ref_lines != rhs.show_ref_lines )
    throw std::runtime_error( "Show ref lines in lhs and rhs are not the same" );
  
  if( fabs(lhs.lower_display_energy - rhs.lower_display_energy) > 1.0e-3 )
    throw std::runtime_error( "Lower display energy in lhs and rhs are not the same" );
  
  if( fabs(lhs.upper_display_energy - rhs.upper_display_energy) > 1.0e-3 )
    throw std::runtime_error( "Upper display energy in lhs and rhs are not the same" );
}//RelActAutoGuiState::equalEnough

#endif // PERFORM_DEVELOPER_CHECKS
}//namespace RelActCalcAuto
