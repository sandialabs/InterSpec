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

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

#include "Wt/WDateTime"
#include "Wt/WApplication"
#include "Wt/WLocalDateTime"

#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#undef ERROR
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
#include "InterSpec/RelEffChart.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/RelActCalc_imp.hpp"
#include "InterSpec/RelActCalcAuto_imp.hpp"
#include "InterSpec/RelActCalc_CeresJetTraits.hpp"

#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#undef ERROR
#endif

using namespace std;
using namespace XmlUtils;

/** Right now if the user specifies a peak energy, we will just multiple the peak sigma (i.e.
 FWHM/2.35482) by the following value to define the ROI on either side of the mean.
 This value is arbitrary, and this is a very niave way to do things, but good enough for
 development purposes, at the moment.
 */
#define DEFAULT_PEAK_HALF_WIDTH_SIGMA 5.0


// 20250324 HACK to test fitting peak skew, or use a hard-coded preset value
//#define PEAK_SKEW_HACK 0 //No hacking
//#define PEAK_SKEW_HACK 1  //Fix the peak skew, and dont fit it
#define PEAK_SKEW_HACK 2  //First fit without peak skew, then fit with peak skew, starting with hard-coded value for DSCB

namespace
{
void sort_rois_by_energy( vector<RelActCalcAuto::RoiRange> &rois )
{
  std::sort( begin(rois), end(rois), []( const RelActCalcAuto::RoiRange &lhs,
                                        const RelActCalcAuto::RoiRange &rhs ) -> bool {
    return lhs.lower_energy < rhs.lower_energy;
  });
}//void sort_rois( vector<RoiRange> &rois )
}//namespace

namespace RelActCalcAutoImp
{
const double ns_decay_act_mult = SandiaDecay::MBq;


// Forward declaration
struct RelActAutoCostFcn;

struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct
 

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
          
        case RelActCalcAuto::FwhmForm::NotApplicable:
          assert( 0 );
          throw runtime_error( "NotApplicable should not be used here" );
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
    cout << "            case RelActCalcAuto::FwhmForm::NotApplicable:" << endl;
    cout << "            assert( 0 );" << endl;
    cout << "            throw runtime_error( \"NotApplicable should not be used here\" );" << endl;
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
    cout << "            case RelActCalcAuto::FwhmForm::NotApplicable:" << endl;
    cout << "            assert( 0 );" << endl;
    cout << "            throw runtime_error( \"NotApplicable should not be used here\" );" << endl;
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
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 3)) );
        if( parameters.size() < (fwhm_start + 3) )
          parameters.resize( fwhm_start + 3, 0.0 );
        parameters[fwhm_start + 0] = 1.54;
        parameters[fwhm_start + 1] = 0.264;
        parameters[fwhm_start + 2] = 0.33;
        break;
        
      case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 3)) );
        if( parameters.size() < (fwhm_start + 3) )
          parameters.resize( fwhm_start + 3, 0.0 );
        parameters[fwhm_start + 0] = 1.86745;
        parameters[fwhm_start + 1] = 0.00216761;
        parameters[fwhm_start + 2] = 27.9835;
        break;
        
      case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 2)) );
        if( parameters.size() < (fwhm_start + 2) )
          parameters.resize( fwhm_start + 2, 0.0 );
        parameters[fwhm_start + 0] = 1.0;
        parameters[fwhm_start + 1] = 0.035;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_2:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 2)) );
        if( parameters.size() < (fwhm_start + 2) )
          parameters.resize( fwhm_start + 2, 0.0 );
        parameters[fwhm_start + 0] = 2.10029;
        parameters[fwhm_start + 1] = 2.03657;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_3:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 3)) );
        if( parameters.size() < (fwhm_start + 3) )
          parameters.resize( fwhm_start + 3, 0.0 );
        parameters[fwhm_start + 0] = 2.26918;
        parameters[fwhm_start + 1] = 1.54837;
        parameters[fwhm_start + 2] = 0.192;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_4:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 4)) );
        if( parameters.size() < (fwhm_start + 4) )
          parameters.resize( fwhm_start + 4, 0.0 );
        parameters[fwhm_start + 0] = 2.49021;
        parameters[fwhm_start + 1] = 0.346357;
        parameters[fwhm_start + 2] = 1.3902;
        parameters[fwhm_start + 3] = -0.294974;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_5:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 5)) );
        if( parameters.size() < (fwhm_start + 5) )
          parameters.resize( fwhm_start + 5, 0.0 );
        parameters[fwhm_start + 0] = 2.5667;
        parameters[fwhm_start + 1] = -0.333729;
        parameters[fwhm_start + 2] = 2.59812;
        parameters[fwhm_start + 3] = -0.991013;
        parameters[fwhm_start + 4] = 0.124209;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_6:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 6)) );
        if( parameters.size() < (fwhm_start + 6) )
          parameters.resize( fwhm_start + 6, 0.0 );
        parameters[fwhm_start + 0] = 2.51611;
        parameters[fwhm_start + 1] = 0.318145;
        parameters[fwhm_start + 2] = 0.838887;
        parameters[fwhm_start + 3] = 0.734524;
        parameters[fwhm_start + 4] = -0.568632;
        parameters[fwhm_start + 5] = 0.0969217;
        break;

      case RelActCalcAuto::FwhmForm::NotApplicable:
        assert( 0 );
        throw runtime_error( "NotApplicable should not be used here" );
        break;
    }//switch( options.fwhm_form )
  }else
  {
    // The following parameters fit from the GADRAS parameters {-6.5f, 7.5f, 0.55f}, using the
    // fit_nominal_gadras_pars() commented out above.
    switch( fwhm_form )
    {
      case RelActCalcAuto::FwhmForm::Gadras:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 3)) );
        if( parameters.size() < (fwhm_start + 3) )
          parameters.resize( fwhm_start + 3, 0.0 );
        parameters[fwhm_start + 0] = -6.5;
        parameters[fwhm_start + 1] = 7.5;
        parameters[fwhm_start + 2] = 0.55;
        break;
        
      case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 3)) );
        if( parameters.size() < (fwhm_start + 3) )
          parameters.resize( fwhm_start + 3, 0.0 );
        parameters[fwhm_start + 0] = -592.865;
        parameters[fwhm_start + 1] = 4.44776;
        parameters[fwhm_start + 2] = 21173.6;
        break;
      
      case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 2)) );
        if( parameters.size() < (fwhm_start + 2) )
          parameters.resize( fwhm_start + 2, 0.0 );
        parameters[fwhm_start + 0] = -7.0;
        parameters[fwhm_start + 1] = 2.0;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_2:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 2)) );
        if( parameters.size() < (fwhm_start + 2) )
          parameters.resize( fwhm_start + 2, 0.0 );
        parameters[fwhm_start + 0] = -146.632;
        parameters[fwhm_start + 1] = 3928.7;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_3:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 3)) );
        if( parameters.size() < (fwhm_start + 3) )
          parameters.resize( fwhm_start + 3, 0.0 );
        parameters[fwhm_start + 0] = -101.518;
        parameters[fwhm_start + 1] = 3037.91;
        parameters[fwhm_start + 2] = 555.973;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_4:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 4)) );
        if( parameters.size() < (fwhm_start + 4) )
          parameters.resize( fwhm_start + 4, 0.0 );
        parameters[fwhm_start + 0] = -68.9708;
        parameters[fwhm_start + 1] = 2334.29;
        parameters[fwhm_start + 2] = 1873.59;
        parameters[fwhm_start + 3] = -428.651;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_5:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 5)) );
        if( parameters.size() < (fwhm_start + 5) )
          parameters.resize( fwhm_start + 5, 0.0 );
        parameters[fwhm_start + 0] = -37.7014;
        parameters[fwhm_start + 1] = 1605.68;
        parameters[fwhm_start + 2] = 4201.01;
        parameters[fwhm_start + 3] = -2230.26;
        parameters[fwhm_start + 4] = 383.314;
        break;
        
      case RelActCalcAuto::FwhmForm::Polynomial_6:
        assert( parameters.empty() || (parameters.size() >= (fwhm_start + 6)) );
        if( parameters.size() < (fwhm_start + 6) )
          parameters.resize( fwhm_start + 6, 0.0 );
        parameters[fwhm_start + 0] = -11.4349;
        parameters[fwhm_start + 1] = 947.497;
        parameters[fwhm_start + 2] = 7088.34;
        parameters[fwhm_start + 3] = -5991.02;
        parameters[fwhm_start + 4] = 2192.23;
        parameters[fwhm_start + 5] = -286.53;
        break;

      case RelActCalcAuto::FwhmForm::NotApplicable:
        assert( 0 );
        throw runtime_error( "NotApplicable should not be used here" );
        break;
    }//switch( options.fwhm_form )
  }//if( highres ) / else
}//fill_in_default_start_fwhm_pars()


/** @brief Get the fwhm coefficients to use in the non-linear fit to data.
 * 
 * @param fwhm_form 
 * @param fwhm_estimation_method 
 * @param all_peaks 
 * @param highres 
 * @param lowest_energy The lowest energy in the analysis energy range; used only when converting 
 *        DetectorPeakResponse FWHM info to a different form, and is not used when fitting FWHM from data.
 * @param highest_energy The highest energy in the analysis energy range; used only when converting 
 *        DetectorPeakResponse FWHM info to a different form, and is not used when fitting FWHM from data.
 * @param input_drf The input detector response function; for some FwhmEstimationMethod, this
 *                  will be used to get the FWHM coefficients, otherwise the all_peaks will be used.
 * @param [out] paramaters The FWHM coefficients to use in the non-linear fit to data.
 * @param [out] warnings Any warnings generated by this function.
 * @return std::shared_ptr<const DetectorPeakResponse> The detector response function with the
 *         FWHM coefficients set - this may be a new object, or the same as the input_drf.  
 *         If a new object, it will have the efficiency of the input_drf (if not-null, otherwise a flat 
 *         efficiency of 1.0 if input was null), but the FWHM coefficients will be set to the values 
 *         returned by this function.  If a FwhmEstimationMethod is specified as FixedToDetectorEfficiency,
 *         or StartingFromDetectorEfficiency, then the returned detector, will be same as input.  If
 *         the FwhmEstimationMethod is StartFromDetEffOrPeaksInSpectrum, then the returned detector
 *         will only be different from the input if the FWHM had to be estimated from the data (or input was null)
 * 
 * Throws exception if error is encountered.
 */
std::shared_ptr<const DetectorPeakResponse> get_fwhm_coefficients( const RelActCalcAuto::FwhmForm fwhm_form, 
            const RelActCalcAuto::FwhmEstimationMethod fwhm_estimation_method,
            const std::vector<std::shared_ptr<const PeakDef>> all_peaks, 
            const bool highres,
            const double lowest_energy,
            const double highest_energy,
            std::shared_ptr<const DetectorPeakResponse> input_drf,
            vector<double> &paramaters,
            vector<string> &warnings )  
{
  paramaters.clear();
  const size_t num_fwhm_pars = num_parameters(fwhm_form);

  if( (fwhm_estimation_method == RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency)
      != (fwhm_form == RelActCalcAuto::FwhmForm::NotApplicable) )
  {
    throw runtime_error( "When FwhmEstimationMethod is specified as FixedToDetectorEfficiency,"
                         " FWHM form must be set to NotApplicable, and vice-versa." );
  }//if( check FWHM type and method are compatible )

  
  bool use_drf_fwhm = false, estimate_fwhm_from_data = false;
  switch( fwhm_estimation_method )
  {
    case RelActCalcAuto::FwhmEstimationMethod::StartingFromAllPeaksInSpectrum:      
    case RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum:      
      use_drf_fwhm = false;
      estimate_fwhm_from_data = true;
    break;

    case RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum:
      use_drf_fwhm = (input_drf && input_drf->hasResolutionInfo());
      estimate_fwhm_from_data = !use_drf_fwhm;
    break;

    case RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency:
      use_drf_fwhm = true;             //dead store
      estimate_fwhm_from_data = false; //dead store
         
      if( fwhm_form != RelActCalcAuto::FwhmForm::NotApplicable )
      {  
        throw runtime_error( "When FwhmEstimationMethod is specified as FixedToDetectorEfficiency,"
                             " FWHM form must be set to NotApplicable." );
      }//if( fwhm_form != RelActCalcAuto::FwhmForm::NotApplicable )

      if( input_drf && input_drf->hasResolutionInfo() )
        return input_drf;

      throw runtime_error( "The detector efficiency function does not have FWHM information;"
                           " please select a DRF with FWHM information, or select to use"
                           " FWHM from peaks in spectrum." );
    break;
      

    case RelActCalcAuto::FwhmEstimationMethod::StartingFromDetectorEfficiency:      
      use_drf_fwhm = true;
      estimate_fwhm_from_data = false;
    break;
  }//switch( fwhm_estimation_method )


  if( use_drf_fwhm && (!input_drf || !input_drf->hasResolutionInfo()) )
  {
    throw runtime_error( "The detector efficiency function does not have FWHM information;"
                          " please select a DRF with FWHM information, or select to use"
                          " FWHM from peaks in spectrum." );
  }//if( use_drf_fwhm && (!input_drf || !input_drf->hasResolutionInfo()) )


  if( use_drf_fwhm )
  {
    const DetectorPeakResponse::ResolutionFnctForm drf_fwhm_type = input_drf->resolutionFcnType();
    const vector<float> drfpars = input_drf->resolutionFcnCoefficients();
      
    assert( drf_fwhm_type != DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm );
    assert( fwhm_estimation_method != RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency );
      
    bool needToFitOtherType = false;
      
    switch( fwhm_form )
    {
      case RelActCalcAuto::FwhmForm::Gadras:
      // We are fitting to the GADRAS functional form
        assert( num_parameters(fwhm_form) == 3 );
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
        assert( num_parameters(fwhm_form) == (static_cast<size_t>(fwhm_form)-1) );
          
        needToFitOtherType = ((drf_fwhm_type != DetectorPeakResponse::kSqrtPolynomial)
                                || (drfpars.size() != num_parameters(fwhm_form)) );
      break;

      case RelActCalcAuto::FwhmForm::NotApplicable:
      {
        needToFitOtherType = false;
        throw runtime_error( "The detector efficiency function does not have FWHM information;"
                              " please select a DRF with FWHM information, or select to use"
                              " FWHM from peaks in spectrum." );
        break;
      }//case RelActCalcAuto::FwhmForm::NotApplicable:
    }//switch( fwhm_form )
     
    if( !needToFitOtherType )
    {
      // We are done here, return the input DetectorPeakResponse, and update the FWHM parameters, 
      //  with the values from the input DetectorPeakResponse
      assert( num_parameters(fwhm_form) == drfpars.size() );

      paramaters.resize( drfpars.size() );
      for( size_t i = 0; i < drfpars.size(); ++i )
        paramaters[i] = static_cast<double>( drfpars[i] );

      return input_drf;
    }//if( !needToFitOtherType )

    // We need to convert from the DetectorPeakResponse FWHM type to the FWHM type we want to use.
    // Make a vector of ~20 equally spaced peaks, with uncert 5% that peaks width -
    //  fairly arbitrary
    const double delta_energy = 0.05*(highest_energy - lowest_energy);
        
    auto fake_peaks = make_shared<std::deque< std::shared_ptr<const PeakDef> > >();
    for( double ene = lowest_energy; ene <=(1.001*highest_energy); ene += delta_energy )
    {
      const float sigma = input_drf->peakResolutionSigma(ene);
      auto p = make_shared<PeakDef>( ene, sigma, 1000.0 );
      p->setSigmaUncert( 0.05*sigma );
      fake_peaks->push_back( p );
    }
      
    DetectorPeakResponse::ResolutionFnctForm formToFit;
    switch( fwhm_form )
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

      case RelActCalcAuto::FwhmForm::NotApplicable:
        assert( 0 );
        throw logic_error( "FwhmForm::NotApplicable should not have made it here" );
      break;
    }//switch( fwhm_form )
        
    try
    {
      vector<float> new_sigma_coefs, sigma_coef_uncerts;
      MakeDrfFit::performResolutionFit( fake_peaks, formToFit, static_cast<int>(num_fwhm_pars),
                                          new_sigma_coefs, sigma_coef_uncerts );
      
      assert( new_sigma_coefs.size() == num_fwhm_pars );
          
      paramaters.resize( new_sigma_coefs.size() );
      for( size_t i = 0; i < new_sigma_coefs.size(); ++i )
        paramaters[i] = static_cast<double>( new_sigma_coefs[i] );

      estimate_fwhm_from_data = false;

      // We are done here, return the input DetectorPeakResponse, with the updated FWHM parameters
      return input_drf;
    }catch( std::exception &e )
    {
      if( fwhm_estimation_method == RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum )
      {
        estimate_fwhm_from_data = true;
      }else
      {
        // We failed to convert from one FWHM type to another
        throw runtime_error( "Failed to convert detectors efficiency functions FWHM info to "
                             + std::string(to_str(fwhm_form)) + " initial FWHM estimation."
                            " Issue:" + string(e.what()) );
      }//if( fwhm_estimation_method == StartFromDetEffOrPeaksInSpectrum ) / else
    }//try / catch
  }//if( use_drf_fwhm )

  assert( estimate_fwhm_from_data );
  assert( fwhm_form != RelActCalcAuto::FwhmForm::NotApplicable );
  assert( (fwhm_estimation_method != RelActCalcAuto::FwhmEstimationMethod::StartingFromDetectorEfficiency)
          && (fwhm_estimation_method != RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency) );

  DetectorPeakResponse::ResolutionFnctForm form_to_fit;
  int fit_order;
  switch( fwhm_form )
  {
    case RelActCalcAuto::FwhmForm::Gadras:
      fit_order = 3;
      form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
    break;
          
    case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:
      fit_order = 3;
      form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse;
    break;
   
    case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy:
      fit_order = 2;
      form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
    break;
          
    case RelActCalcAuto::FwhmForm::Polynomial_2:
    case RelActCalcAuto::FwhmForm::Polynomial_3:
    case RelActCalcAuto::FwhmForm::Polynomial_4:
    case RelActCalcAuto::FwhmForm::Polynomial_6:
    case RelActCalcAuto::FwhmForm::Polynomial_5:
      form_to_fit = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      fit_order = static_cast<int>( num_parameters(fwhm_form) );
    break;
          
    case RelActCalcAuto::FwhmForm::NotApplicable:
      assert( 0 );
      throw runtime_error( "NotApplicable should not be used here" );
    break;
  }//switch( fwhm_form )

  try
  {
    //Instead of the below, we could potentially do
    //  new_drf->fitResolution( peaks_deque, spectrum, form_to_fit );
    
    vector<float> fwhm_paramatersf, uncerts;
    auto peaks_deque = make_shared<deque<shared_ptr<const PeakDef>>>(begin(all_peaks), end(all_peaks));
    MakeDrfFit::performResolutionFit( peaks_deque, form_to_fit, fit_order, fwhm_paramatersf, uncerts );
      
    vector<pair<double,shared_ptr<const PeakDef>>> distances;
    for( const auto &p : *peaks_deque )
    {
      const float mean = static_cast<float>(p->mean());
      const float pred_fwhm = DetectorPeakResponse::peakResolutionFWHM( mean, form_to_fit, fwhm_paramatersf );

      const double frac_diff = fabs( p->fwhm() - pred_fwhm ) / p->fwhm();
      if( !IsNan(frac_diff) && !IsInf(frac_diff) )
        distances.emplace_back( frac_diff, p );
    }//for( const auto &p : initial_fit_peaks )
      
    std::sort( begin(distances), end(distances), []( const auto &lhs, const auto &rhs) -> bool {
      return lhs.first > rhs.first;
    } );

    // Limit to un-selecting max of 20% of peaks (arbitrarily chosen), if they deviate
    // more than 17.5% (again, arbitrarily chosen) from the fit.
    const size_t max_remove = (distances.size() < 8) ? size_t(0) : static_cast<size_t>( std::ceil( 0.2*distances.size() ) );
    auto filtered_peaks = make_shared<deque<shared_ptr<const PeakDef>>>();
    for( size_t index = 0, num_removed = 0; index < distances.size(); ++index )
    {
      if( (distances[index].first > 0.175 ) && (num_removed < max_remove) ) //0.175 chosen arbitrarily
      {
        ++num_removed;
        continue;
      }
        
      filtered_peaks->push_back( distances[index].second );
    }//for( size_t index = 0, num_removed = 0; index < distances.size(); ++index )

    if( (max_remove + filtered_peaks->size()) < distances.size() )
      cerr << "RelActCalcAuto - failed logic: max_remove=" << max_remove 
      << ", filtered_peaks->size()=" << filtered_peaks->size() 
      << ", distances.size()=" << distances.size() << endl;

    assert( (max_remove + filtered_peaks->size()) >= distances.size() );

    if( filtered_peaks->size() != all_peaks.size() )
    {
      try
      {
        vector<float> new_result, new_result_uncerts;
        MakeDrfFit::performResolutionFit( filtered_peaks, form_to_fit,
                                          fit_order, new_result, new_result_uncerts );
        fwhm_paramatersf = new_result;
        uncerts.swap( new_result_uncerts );
      }catch( std::exception &e )
      {
        warnings.push_back( "Failed to refine FWHM fit from data: " + string(e.what()) + ".  Will use initial estimate." );
      }
    }//if( filtered_peaks->size() != all_peaks.size() )

    paramaters.resize( fwhm_paramatersf.size() );
    for( size_t i = 0; i < fwhm_paramatersf.size(); ++i )
      paramaters[i] = static_cast<double>( fwhm_paramatersf[i] );
  }catch( std::exception &e )
  {
    paramaters.clear();
    fill_in_default_start_fwhm_pars( paramaters, 0, highres, fwhm_form );
    warnings.push_back( "Failed to estimate FWHM from data: " + string(e.what()) + ".  Using default FWHM parameters." );
  }
  
  vector<float> fwhm_pars_float( paramaters.size() );
  for( size_t i = 0; i < paramaters.size(); ++i )
    fwhm_pars_float[i] = static_cast<float>( paramaters[i] );

  shared_ptr<DetectorPeakResponse> new_drf;
  if( input_drf )
    new_drf = make_shared<DetectorPeakResponse>(*input_drf);
  
  if( !new_drf )
  {
    const string drfpaths = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GenericGadrasDetectors" );
    const string drf_dir = SpecUtils::append_path( drfpaths, highres ? "HPGe 40%" : "LaBr 10%" );
    try
    {
      new_drf = std::make_shared<DetectorPeakResponse>();
      new_drf->fromGadrasDirectory( drf_dir );
      new_drf->setFwhmCoefficients( {}, DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm ); //just to be sure
    }catch( std::exception &e )
    {
      //throw runtime_error( "RelActAutoCostFcn: failed to open default DRF." );
      cerr << "RelActAutoCostFcn: failed to open default DRF." << endl;
      assert( 0 );
      new_drf.reset();
    }
  }//if( !new_drf )
  
  assert( !!new_drf );
  if( !new_drf )
  {
    new_drf = make_shared<DetectorPeakResponse>( "FLAT", "FLAT" );
    const vector<float> drf_coefs{ 0.0f, 0.0f, 0.0f, 0.0f }, uncerts;
    new_drf->fromExpOfLogPowerSeriesAbsEff( drf_coefs, uncerts,
                                           25*PhysicalUnits::cm,
                                           2*PhysicalUnits::cm,
                                           PhysicalUnits::keV,
                                           static_cast<float>(lowest_energy),
                                           static_cast<float>(highest_energy),
                                           DetectorPeakResponse::EffGeometryType::FarField );
  }//
  
  new_drf->setFwhmCoefficients( fwhm_pars_float, form_to_fit );

  return new_drf;
};//get_fwhm_coefficients(...)



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
  
  Will be set to -1.0 if age is controlled by another nuclide, an x-ray, or a reaction.
  */
  double age_multiple = 1.0;

  struct EnergyYield
  {
    double energy = 0.0;
    double yield = 0.0;
      
    size_t transition_index = 0;
    const SandiaDecay::Transition *transition = nullptr;
    PeakDef::SourceGammaType gamma_type;
    
    const SandiaDecay::Element *element = nullptr;
    const ReactionGamma::Reaction *reaction = nullptr;
  };//struct EnergyYield
  
  std::shared_ptr<const vector<EnergyYield>> nominal_gammas;
  
  // Implemeneted below RelActAutoCostFcn, so we can access it
  NucInputGamma( const RelActCalcAuto::NucInputInfo &info, const RelActAutoCostFcn * const cost_fcn );
  NucInputGamma() = delete;
  NucInputGamma( const NucInputGamma & ) = default;
  NucInputGamma( NucInputGamma && ) = default;
  NucInputGamma &operator=( const NucInputGamma & ) = default;
  NucInputGamma &operator=( NucInputGamma && ) = default;
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
  /** How to perform differentiation for jacobians; either "auto" or numeric.

   We should use auto-diff generally, since it is faster, and should be less-error-prone, but to check on
   things, you can choose to use numeric differentiation.
   */
  static const bool sm_use_auto_diff = true;

  /** Ceres uses a compile time constant for the number of parameters it performs for "auto" differentiation for; this
   variable constrols how many parameters will be differentiated at a time.

   It looks like using Jets<double,32> is fastest for an example problem, however it really
   doesnt seem to matter much (Wall times of {0.1728, 0.1789, 0.1665, 0.1608}, for strides of
   {4, 8, 16, 32}, respectively
  */
  static const size_t sm_auto_diff_stride_size = 32;


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
  
#if( PEAK_SKEW_HACK == 2 )
  bool m_skip_skew; // 20250324 HACK to test fitting peak skew
#endif

  /** When m_options.additional_br_uncert we will allow each peak range (see `cluster_photopeaks(...)`) to vary within this uncertainty;
  This variable defines the lower and upper energy of each peak range.
   */
  std::vector<std::pair<double,double>> m_peak_ranges_with_uncert;

  /** We will use a Ceres paramater to scale the yield of each peak-range.
   The Ceres paramater is defined like:

   ```
     const T par_value = x[range_index];
     const T num_sigma_from_nominal = (par_value - sm_peak_range_uncert_offset)/sm_peak_range_uncert_par_scale;
     const T yield_multiple = 1.0 + num_sigma_from_nominal*m_options.additional_br_uncert;
   ```
   */
  static constexpr double sm_peak_range_uncert_offset = 1.0;
  static constexpr double sm_peak_range_uncert_par_scale = 0.1;

  
  mutable std::mutex m_aged_gammas_cache_mutex;
  /** When we are fitting the age of nuclides, performing the decay can take a significant amount of the time, so we'll cache results.
   
   Note the map key is {nuclide,age}, where age is kept as a float to hopefully help minimize how many entries get put into the cache.
   
   For an example CZT Pu problem, fitting the Pu age, the cache hit about 99% of calls,
   and in a debug build, the average call time to `eval(...)` went from 17 ms, to 5.5 ms.
   The cache size grew to about 1100 entries, probably taking a bit over 30 MB of memory.
   
   Note that this is usefull when fitting ages.
   */
  mutable std::map<std::pair<const void *,float>,std::shared_ptr<const std::vector<NucInputGamma::EnergyYield>>> m_aged_gammas_cache;
  

  /** The initial estimate of the area of each free peak.

   This is used to multiple the free peak area paramater by to get the peak amplitude.
   */
  vector<double> m_free_peak_area_multiples;
  
  /** Tracks if a solution has been found yet or not (e.g., minimization is complete).
   
   Currently, only used to make sure we arent calling functions intended for use only after
   an an intitial solution has been found (e.g., depends on covariance matrix).
   */
  bool m_solution_finished;
  
  /** Thread pool used to calculate the each ROI.
   
   Each ROI is submitted to this pool, but computation of peaks, and counts from the peaks,
   may be computed using `SpecUtilsAsync::ThreadPool` (when it seems benfecial).
   The reason `SpecUtilsAsync::ThreadPool` is not used recursively is that there can
   be some delays when recursively using this thread pool, especially using Apple GCD.
   
   We wont call `boost::asio::thread_pool::join`, so according to the boost
   header, we should be good to go from a safety standpoint...
   */
  mutable boost::asio::thread_pool m_pool;
  
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
  
  class CeresStepSummaryCallback : public ceres::IterationCallback
  {
    const RelActAutoCostFcn * const m_cost_fcn;
    
  public:
    CeresStepSummaryCallback( const RelActAutoCostFcn *cost_fcn ) : m_cost_fcn( cost_fcn ) {}
    
    virtual ceres::CallbackReturnType operator()(const ceres::IterationSummary &sum )
    {
      cout << "Iteration " << sum.iteration << ", " << sum.iteration_time_in_seconds << " seconds: " << endl;
      
      return ceres::CallbackReturnType::SOLVER_CONTINUE;
    }
  };//class CeresStepSummaryCallback
  
  
  
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
#if( PEAK_SKEW_HACK == 2 )
  m_skip_skew( false ), // 20250324 HACK to test fitting peak skew
#endif
  m_peak_ranges_with_uncert{},
  m_aged_gammas_cache_mutex{},
  m_aged_gammas_cache{},
  m_free_peak_area_multiples{},
  m_solution_finished( false ),
  m_pool{ std::max(4u, std::thread::hardware_concurrency()) },
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
    
    const size_t num_rel_eff_curves = options.rel_eff_curves.size();
    m_nuclides.resize( num_rel_eff_curves );
    
    for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_input = options.rel_eff_curves[rel_eff_index];
      const vector<RelActCalcAuto::NucInputInfo> &input_nuclides = rel_eff_input.nuclides;
      
      if( input_nuclides.empty() )
        throw runtime_error( "RelActAutoCostFcn: no nuclides specified." );
    
      vector<NucInputGamma> &these_nuclides = m_nuclides[rel_eff_index];

      for( const RelActCalcAuto::NucInputInfo &n : input_nuclides )
      {
        if( RelActCalcAuto::is_null(n.source) )
          throw runtime_error( "RelActAutoCostFcn: null Nuclide." );
      
        for( const NucInputGamma &pn : these_nuclides )
        {
          if( (n.source == pn.source) )
            throw runtime_error( "RelActAutoCostFcn: duplicate nuclide (" + RelActCalcAuto::to_name(n.source) + ")." );
        }

        these_nuclides.emplace_back( n, this );
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
    }//for( const auto &rel_eff_input : options.rel_eff_curves )
    
    
    if( cancel_calc && cancel_calc->load() )
      throw runtime_error( "User cancelled calculation." );
    
    
    if( !drf || !drf->isValid() || !drf->hasResolutionInfo() )
      throw runtime_error( "RelActAutoCostFcn: input DRF must be valid and have resolution info." );
    
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
          assert( n.nominal_gammas );
          for( const auto &g : *n.nominal_gammas )
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

      assert( m_drf && m_drf->isValid() && m_drf->hasResolutionInfo() );
      const double det_fwhm = m_drf->peakResolutionFWHM( peak.energy );
      
      double data_count = m_spectrum->gamma_integral(peak.energy - 0.5*det_fwhm, peak.energy + 0.5*det_fwhm);

      // Enforce using a peak area estimate of least 1 count; the fit paramater will scale this, so its fine.
      //  and we are only really concerned with really large peak areas for this scaling.
      data_count = std::max( 0.76, data_count );
      
      // Assume data is entirely due to peak; mean +- 0.5*fwhm is ~76% gaussian are
      m_free_peak_area_multiples.push_back( data_count / 0.76 );
    }//for( const auto &peak : m_options.floating_peaks )
    
    assert( m_free_peak_area_multiples.size() == m_options.floating_peaks.size() );
    
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
    
    const size_t num_rel_eff_curves = options.rel_eff_curves.size();
    
    RelActCalcAuto::RelActAutoSolution solution;
    
    DoWorkOnDestruct setFinalTime( [&solution,start_time](){
      const auto end_time = std::chrono::high_resolution_clock::now();
      solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
    });
    
    // We will make a (potentially background subtracted) spectrum - we also need to track the
    //  uncertainties in each channel
    shared_ptr<SpecUtils::Measurement> spectrum;
    vector<float> channel_counts, channel_count_uncerts;
    vector<double> starting_fwhm_paramaters; //The initial FWHM parameters to use in the fit

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

      if( !foreground )
        throw runtime_error( "Not valid foreground provided." );

      if( (foreground->live_time() < 0.01) || (foreground->real_time() < 0.01) )
        throw runtime_error( "Foreground must have non-zero live and real times." );

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
          const double uncert_2 = fore_counts + lt_sf*lt_sf*back_counts;
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
        all_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( spectrum, nullptr, {}, false, highres );
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
    
      // Get the initial FWHM parameters to use in the fit
      try
      {
        double lowest_energy = 10000.0, highest_energy = 0.0;
        for( const RelActCalcAuto::RoiRange &r : options.rois )
        {
          lowest_energy = std::min( lowest_energy, r.lower_energy );
          highest_energy = std::max( highest_energy, r.upper_energy );
        }
          
        // This next call may return `input_drf` if its valid and has FWHM info,
        //  otherwise it will add FWHM estimate, or if it isnt valid at all, it will
        //  load a default efficiency function
        solution.m_drf = get_fwhm_coefficients( options.fwhm_form, options.fwhm_estimation_method, all_peaks,
                                              highres, lowest_energy, highest_energy, input_drf,
                                              starting_fwhm_paramaters, solution.m_warnings );
      }catch( std::exception &e )
      {
        throw runtime_error( "Failed to get initial FWHM: " + string(e.what()) + "." );
      }//try / catch (getting initial FWHM parameters)

      assert( solution.m_drf  && solution.m_drf->isValid() && solution.m_drf->hasResolutionInfo() );
    
      // Check that the specified Pu242 by correlation method is valid.
      auto check_pu_corr_method = [&]( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve ){
        const vector<RelActCalcAuto::NucInputInfo> &nuclides = rel_eff_curve.nuclides;
      
        // Make a lamda that returns the mass-number of Plu nuclides
        auto pu_iso_present = [&]() -> set<short> {
          set<short> answer;
          for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
          {
            const SandiaDecay::Nuclide * const nuc_nuclide = RelActCalcAuto::nuclide(nuc.source);
            if( nuc_nuclide && (nuc_nuclide->atomicNumber == 94) )
              answer.insert(nuc_nuclide->massNumber);
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
              string msg = "Pu242 correlation method of Bignan95 was specified, but problem did not"
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
              string msg = "Pu242 correlation method  using Pu239-only was specified, but problem did"
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
    
      // Check nucliode constraints are all valid
      for( const RelActCalcAuto::RelEffCurveInput &rel_eff_curve : options.rel_eff_curves ) 
        rel_eff_curve.check_nuclide_constraints();

      options.check_same_hoerl_and_external_shielding_specifications();
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
                                           solution.m_drf, all_peaks, cancel_calc );
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
  
    
    solution.m_cost_functor = cost_functor;
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
    
            
    // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 2 )
    vector<int> initial_const_parameters;
#endif
    
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
      if( ((lowest_energy < 250) && (highest_energy > 650))
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

      assert( parameters.size() >= (cost_functor->m_fwhm_par_start_index + starting_fwhm_paramaters.size()) );
      assert( (options.fwhm_estimation_method != RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency) 
              || (starting_fwhm_paramaters.size() == 0) );
              
      if( parameters.size() < (cost_functor->m_fwhm_par_start_index + starting_fwhm_paramaters.size()) )
        throw logic_error( "Parameters vector is not large enough to hold FWHM parameters??? - shouldnt happen" );

      // Copy the FWHM parameters to the parameters vector
      for( size_t i = 0; i < starting_fwhm_paramaters.size(); ++i )
      {
        const size_t par_index = cost_functor->m_fwhm_par_start_index + i;
        parameters[par_index] = starting_fwhm_paramaters[i];
        if( options.fwhm_estimation_method == RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum )
          constant_parameters.push_back( static_cast<int>(par_index) );
      }

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

      assert( cost_functor->m_nuclides.size() == num_rel_eff_curves );
      if( cost_functor->m_nuclides.size() != num_rel_eff_curves )
        throw logic_error( "Number of nuclides in cost functor does not match number of relative efficiency curves." );

      // Loop over each relative efficiency curve and estimate the starting parameters for the "auto" fit
      //  by doing a rough/course/poor manual fit.
      const RelActCalcAuto::RelEffCurveInput *first_phys_model_curve = nullptr; //We'll track this for sharing external shielding and Hoerl function, if the user has selected it
      size_t first_phys_model_par_start = std::numeric_limits<size_t>::max();
      size_t num_re_nucs_seen = 0, num_re_curve_pars_seen = 0; // To help us keep track of rel eff and act/age index in parameters vector
      for( size_t re_eff_index = 0; re_eff_index < num_rel_eff_curves; ++re_eff_index )
      {
        const double base_rel_eff_uncert = 1.0;
        
        const auto &rel_eff_curve = options.rel_eff_curves[re_eff_index];
        //assert( (rel_eff_curve.rel_eff_eqn_order != 0)
        //       || (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) );

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
              
              pars[b_index] = (0.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;  //(energy/1000)^b - start b at 0, so term is 1.0
              pars[c_index] = (1.0/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;  //c^(1000/energy) - start c at 1, so term is 1
              if( !rel_eff_curve.phys_model_use_hoerl || (options.same_hoerl_for_all_rel_eff_curves && first_phys_model_curve) )
              {
                constant_parameters.push_back( static_cast<int>(b_index) );
                constant_parameters.push_back( static_cast<int>(c_index) );
                
                if( options.same_hoerl_for_all_rel_eff_curves && first_phys_model_curve )
                {
                  // This is just so we can assert on the values later, as a double check
                  pars[b_index] = -1.0;
                  pars[c_index] = -1.0;
                }
              }else
              {
                lower_bounds[b_index] = (0.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;
                upper_bounds[b_index] = (2.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;
                lower_bounds[c_index] = (1.0E-6/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;  //e.x, pow(-0.1889,1000/124.8) is NaN
                upper_bounds[c_index] = (3.0/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;
              }

              if( !first_phys_model_curve )
              {
                first_phys_model_curve = &rel_eff_curve;
                first_phys_model_par_start = this_rel_eff_start;
              }

              break;
            }//case RelActCalc::RelEffEqnForm::FramPhysicalModel:
          }//switch( solution.m_rel_eff_form )
         
          
          // Now we'll put together the "manual" approximation for R.E. and try to solve it
          vector<RelActCalcManual::SandiaDecayNuc> nuc_sources;
          for( const RelActCalcAuto::NucInputInfo &info : nucs )
          {
            RelActCalcManual::SandiaDecayNuc nucinfo;
            nucinfo.nuclide = RelActCalcAuto::nuclide(info.source);
            nucinfo.element = RelActCalcAuto::element(info.source);
            nucinfo.reaction = RelActCalcAuto::reaction(info.source);
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
          manual_input.phys_model_detector = solution.m_drf;
          manual_input.eqn_form = rel_eff_curve.rel_eff_eqn_type;
          manual_input.use_ceres_to_fit_eqn = (manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel);
          manual_input.eqn_order = manual_rel_eff_order;
          manual_input.phys_model_self_atten = rel_eff_curve.phys_model_self_atten;
          manual_input.phys_model_external_attens = rel_eff_curve.phys_model_external_atten;
          manual_input.phys_model_use_hoerl = false;


          const bool use_first_ext_shield = ((rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) 
                                            && options.same_external_shielding_for_all_rel_eff_curves
                                            && (first_phys_model_curve != (&rel_eff_curve)));
          const bool use_first_hoerl = ((rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) 
                                        && options.same_hoerl_for_all_rel_eff_curves
                                        && manual_input.phys_model_use_hoerl
                                        && (first_phys_model_curve != (&rel_eff_curve)));


          if( use_first_ext_shield )
          {
            for( size_t i = 0; i < manual_input.phys_model_external_attens.size(); ++i )
            {
              //assert( 0 ); // this has not been checked as of 20250317

              const size_t ext_att_index = first_phys_model_par_start + 2 + 2*i;
              const shared_ptr<const RelActCalc::PhysicalModelShieldInput> orig_ext_att = manual_input.phys_model_external_attens[i];
              assert( orig_ext_att );
              shared_ptr<RelActCalc::PhysicalModelShieldInput> updated_ext_att = make_shared<RelActCalc::PhysicalModelShieldInput>( *orig_ext_att );
              updated_ext_att->fit_atomic_number = false;
              updated_ext_att->fit_areal_density = false;
              if( !updated_ext_att->material )
                updated_ext_att->atomic_number = parameters[ext_att_index + 0] * RelActCalc::ns_an_ceres_mult;
              updated_ext_att->areal_density = parameters[ext_att_index + 1] * PhysicalUnits::g_per_cm2;
              
              manual_input.phys_model_external_attens[i] = updated_ext_att;
            }
          }//if( use_first_ext_shield )
          
          for( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint : rel_eff_curve.act_ratio_constraints )
          {
            RelActCalcManual::ManualActRatioConstraint manual_constraint;
            assert( !RelActCalcAuto::is_null(constraint.constrained_source) && !RelActCalcAuto::is_null(constraint.controlling_source) );
            if( RelActCalcAuto::is_null(constraint.constrained_source) || RelActCalcAuto::is_null(constraint.controlling_source) )
              throw runtime_error( "Invalid nuclide in activity ratio constraint" );
            
            manual_constraint.m_constrained_nuclide = RelActCalcAuto::to_name(constraint.constrained_source);
            manual_constraint.m_controlling_nuclide = RelActCalcAuto::to_name(constraint.controlling_source);
            manual_constraint.m_constrained_to_controlled_activity_ratio = constraint.constrained_to_controlled_activity_ratio;
            manual_input.act_ratio_constraints.push_back( manual_constraint );
          }

          if( use_first_hoerl )
          {
            // TODO: we should fix the Hoerl parameters here to what we got for first Phys Model Rel Eff Curve;
            //       but actually we arent ever fitting the Hoerl parameters, so this is moot (and we havent implemented a 
            //       way to contrain the Hoerl parameters yet).
          }//if( use_first_hoerl )

          for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )
          {
            // Manual mass fraction constraints only accept a single fixed mass fraction for an isotope.
            RelActCalcManual::MassFractionConstraint manual_constraint;
            manual_constraint.m_nuclide = RelActCalcAuto::to_name(constraint.nuclide);
            manual_constraint.m_mass_fraction_lower = constraint.lower_mass_fraction;
            manual_constraint.m_mass_fraction_upper = constraint.upper_mass_fraction;

            for( const RelActCalcAuto::NucInputInfo &nuc : rel_eff_curve.nuclides )
            {
              const SandiaDecay::Nuclide * const nuc_nuclide = RelActCalcAuto::nuclide(nuc.source);
              if( nuc_nuclide && (nuc_nuclide->atomicNumber == constraint.nuclide->atomicNumber) )
                manual_constraint.m_specific_activities[nuc_nuclide->symbol] = nuc_nuclide->activityPerGram();
            }

            manual_input.mass_fraction_constraints.push_back( manual_constraint );
          }//for( const MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )

          RelActCalcManual::RelEffSolution manual_solution
                               = RelActCalcManual::solve_relative_efficiency( manual_input );
          
          if( manual_rel_eff_order < rel_eff_curve.rel_eff_eqn_order )
            solution.m_warnings.push_back( "Due to a lack of manually fit peaks, the relative"
                                          " efficiency equation order had to be reduced for initial"
                                          " estimate of relative efficiencies and activities." );
          
          const string re_id_str =  (num_rel_eff_curves > 1) ? ("RE_" + std::to_string(re_eff_index)) : string("");
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
            
            parameters[this_rel_eff_start + 0] = 0.0;  //Atomic number
            parameters[this_rel_eff_start + 1] = 0.0;  //Areal density

            size_t manual_index = 0;
            if( rel_eff_curve.phys_model_self_atten )
            {
              if( rel_eff_curve.phys_model_self_atten->fit_atomic_number )
              {
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                parameters[this_rel_eff_start + 0] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index); //Atomic number; note both manual and auto RelEff use RelActCalc::ns_an_ceres_mult
                manual_index += 1;
              }
              
              assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
              parameters[this_rel_eff_start + 1] = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index); //Areal density; both manual and auto RelEff use g/cm2
              
              manual_index += 1;
            }//if( options.phys_model_self_atten )

            for( size_t i = 0; i < rel_eff_curve.phys_model_external_atten.size(); ++i )
            {
              const auto &ext_att = rel_eff_curve.phys_model_external_atten[i];
              assert( ext_att );
              if( !ext_att )
                continue;

              if( use_first_ext_shield )
              {
                parameters[this_rel_eff_start + 2 + 2*i + 0] = -1.0;
                parameters[this_rel_eff_start + 2 + 2*i + 1] = -1.0;
                manual_index += 1; //AN is only a parameter if being fit, AD is always a paramater
                continue;
              }

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
              assert( first_phys_model_curve );

              if( options.same_hoerl_for_all_rel_eff_curves && (first_phys_model_curve != (&rel_eff_curve)) )
              {
                assert( parameters[b_index] == -1.0 );
                assert( parameters[c_index] == -1.0 );
                manual_index += 2;
              }else
              {
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                const double b_par_val = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                parameters[b_index] = b_par_val;
                manual_index += 1;
                
                assert( manual_solution.m_rel_eff_eqn_coefficients.size() > manual_index );
                const double c_par_val = manual_solution.m_rel_eff_eqn_coefficients.at(manual_index);
                parameters[c_index] = c_par_val;
                manual_index += 1;
              }
            }else
            {
              if( options.same_hoerl_for_all_rel_eff_curves && (first_phys_model_curve != (&rel_eff_curve)) )
              {
                assert( parameters[b_index] == -1.0 );
                assert( parameters[c_index] == -1.0 );
              }else
              {
                parameters[b_index] = (0.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;  //(energy/1000)^b
                parameters[c_index] = (1.0/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;  //c^(1000/energy)
              }
            }//if( manual_input.phys_model_use_hoerl ) / else
          }//if( no Physical Model ) / else
          cout << "Starting with initial rel. eff. eqn = " << rel_eff_eqn_str << endl;
          
          const double live_time = (spectrum && (spectrum->live_time() > 0)) ? solution.m_spectrum->live_time() : 1.0f;
          
          for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
          {
            const RelActCalcAuto::NucInputInfo &nuc = rel_eff_curve.nuclides[nuc_num];
            const size_t act_index = cost_functor->nuclide_parameter_index(nuc.source, re_eff_index);
            double rel_act = manual_solution.relative_activity( nuc.name() ) / live_time;

            //Count how many rel eff curves this nuclide is in, and divide by that... we can probably come up with a better solution...
            size_t nuc_count = 0;
            for( const auto &rel_eff_curve : options.rel_eff_curves )
            {
              for( const RelActCalcAuto::NucInputInfo &other_nuc : rel_eff_curve.nuclides )
                nuc_count += (nuc.source == other_nuc.source);
            }
            assert( nuc_count != 0 );
            rel_act /= nuc_count;
            
            assert( cost_functor->m_nuclides[re_eff_index][nuc_num].source == nuc.source );

            bool is_constrained = false;
            for( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &nuc_constraint : rel_eff_curve.act_ratio_constraints )
            {
              assert( !RelActCalcAuto::is_null(nuc_constraint.constrained_source) );
              assert( !RelActCalcAuto::is_null(nuc_constraint.controlling_source) );

              is_constrained = (nuc_constraint.constrained_source == nuc.source); 
              if( is_constrained )
              {
                parameters[act_index] = -1.0;
                cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = -1.0;
                assert( std::find( constant_parameters.begin(), constant_parameters.end(), static_cast<int>(act_index) ) == constant_parameters.end() );
                constant_parameters.push_back( static_cast<int>(act_index) );
                break;
              }
            }//for( const auto &nuc_constraint : rel_eff_curve.act_ratio_constraints )

            for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )
            {
              if( !RelActCalcAuto::nuclide(nuc.source) )
                break;
                
              if( is_constrained )
              {
                assert( constraint.nuclide != RelActCalcAuto::nuclide(nuc.source) ); //We shouldnt have this nuclide controlled by another nuclide, and constrain its mass fraction
                break;
              }

              assert( !RelActCalcAuto::is_null(constraint.nuclide) );
              is_constrained = (RelActCalcAuto::SrcVariant(constraint.nuclide) == nuc.source);
              if( is_constrained )
              {
                //Rel Act paramater is constrained within 0.5 and 1.5, to make mass fraction go between lower and upper values
                //  If there are multiple mass fraction constraints on the same nuclide, and a particular Ceres parameter
                //  solution gives the sum of all the mass fractions to be greater than 1.0, we will just return a say
                //  that particular set of parameters is invalid, even though we could maybe do something more intelligent.

                //TODO: To model mass-fraction constraints, should switch to a paramter that gives total RelAct of an element, and then use a ceres::Manifold to make all the nuclides of the element add up to 1.0 (eg, on a surface).

                parameters[act_index] = 1.0;
                if( constraint.lower_mass_fraction != constraint.upper_mass_fraction )
                {
                  try
                  {
                    const double mf = manual_solution.mass_fraction( constraint.nuclide->symbol );
                    double frac = (mf - constraint.lower_mass_fraction) / (constraint.upper_mass_fraction - constraint.lower_mass_fraction);
                    assert( frac > -0.0002 && frac < 1.0002 );
                    frac = std::min( 1.0, std::max( frac, 0.0 ) );
                    parameters[act_index] = 0.5 + 0.5*frac;
                    cerr << "\n\nStill having trouble navigating to global minima\n\nHacked mass fraction to be halfway in the range it should!\n\n." << endl;
#pragma message("Still having trouble navigating to global minima for RelActAuto. Hacked mass fraction to be halfway in the range it should! - should undo this (but seems to help for some problems).") 
                  }catch( std::exception & )
                  {
                    cerr << "Warning: failed to get initial manual solution rel-eff fraction" << endl;
                  }
                }//if( constraint.lower_mass_fraction != constraint.upper_mass_fraction )

                cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = -1.0;
                
                //If the lower and upper mass fractions are the same, we can just fix the rel act to 1.0
                if( fabs(constraint.lower_mass_fraction - constraint.upper_mass_fraction) 
                    <= 1.0E-6*std::max(constraint.lower_mass_fraction, constraint.upper_mass_fraction) )
                {
                  assert( std::find( constant_parameters.begin(), constant_parameters.end(), static_cast<int>(act_index) ) == constant_parameters.end() );
                  constant_parameters.push_back( static_cast<int>(act_index) );
                }else
                {
                  lower_bounds[act_index] = 0.5;
                  upper_bounds[act_index] = 1.5;
                }
                break;
              }//if( is_constrained )
            }//for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )


            if( is_constrained )
              continue;

            assert( !RelActCalcAuto::is_null(nuc.source) );

            // Check if nuclide has relative activity limits defined, or maybe even fixed.
            if( nuc.min_rel_act.has_value() 
               && nuc.max_rel_act.has_value() 
                && (nuc.min_rel_act.value() == nuc.max_rel_act.value()) )
            {
              parameters[act_index] = 1.0;
              cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = nuc.min_rel_act.value();
              assert( cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple > 0.0 );
              constant_parameters.push_back( static_cast<int>(act_index) );
              cout << "Fixing " << nuc.name() << " rel act to " << nuc.min_rel_act.value() << " for rel eff curve " << re_eff_index << endl;
              continue;
            }

            if( nuc.starting_rel_act.has_value() )
            {
              cout << "Will use starting rel act of " << nuc.starting_rel_act.value() << " for " << nuc.name() 
                   << " on rel eff curve " << re_eff_index << ", instead of initial estimate of " << rel_act << endl;
              rel_act = nuc.starting_rel_act.value();
            }
            if( (rel_act < 1.0E-16) || IsInf(rel_act) || IsNan(rel_act) )
              rel_act = 1.0;
            
            lower_bounds[act_index] = (nuc.min_rel_act.has_value() ? nuc.min_rel_act.value() : 0.0) / rel_act;
            if( nuc.max_rel_act.has_value() )
              upper_bounds[act_index] = nuc.max_rel_act.value() / rel_act;

            cout << "Updating initial activity estimate for " << nuc.name() << " from "
                   << parameters[act_index] << " to " << rel_act << endl;          
            
            cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = rel_act;
            assert( rel_act > 0.0 );

            parameters[act_index] = 1.0;
          }//for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
          
          
          succesfully_estimated_re_and_ra = true;
        }catch( std::exception &e )
        {
          cerr << "Failed to do initial estimate on RelEff curve: " << e.what() << endl;
          
          solution.m_warnings.push_back( "Initial estimate of relative efficiency curve failed ('"
                                        + string(e.what())
                                        + "'), using a flat line as a starting point" );
        }//try / catch
          
          
        for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )
        {
          const RelActCalcAuto::NucInputInfo &nuc = rel_eff_curve.nuclides[nuc_num];
          
          const size_t act_index = cost_functor->nuclide_parameter_index( nuc.source, re_eff_index );
          const size_t age_index = act_index + 1;
          
          assert( !RelActCalcAuto::is_null(nuc.source) );

          const SandiaDecay::Nuclide * const nuc_nuclide = RelActCalcAuto::nuclide(nuc.source);
          const SandiaDecay::Element * const nuc_element = RelActCalcAuto::element(nuc.source);
          const ReactionGamma::Reaction * const nuc_reaction = RelActCalcAuto::reaction(nuc.source);

          assert( nuc_nuclide || nuc_element || nuc_reaction );
          assert( !nuc_nuclide || (nuc.age >= 0.0) );
          assert( nuc_nuclide || (nuc.age <= 0.0) );

          if( nuc_nuclide && nuc_nuclide->isStable() )
            throw runtime_error( "Invalid stable nuclide (" + nuc.name() + ")" );
          
          if( nuc_nuclide && (nuc.age < 0.0) )
            throw runtime_error( "Invalid nuclide (" + nuc.name() + ") age" );
          
          if( !nuc_nuclide && (nuc.age > 0.0) )
            throw runtime_error( "Invalid age for non-nuclide source (" + nuc.name() + ")" );
          
          assert( re_eff_index < cost_functor->m_nuclides.size() );
          assert( nuc_num < cost_functor->m_nuclides[re_eff_index].size() );
          assert( cost_functor->m_nuclides[re_eff_index][nuc_num].source == nuc.source );
          
          NucInputGamma &nuc_info = cost_functor->m_nuclides[re_eff_index][nuc_num];

          bool fit_age = nuc.fit_age;
          if( nuc_nuclide )
          {
            nuc_info.age_multiple = nuc.age;
            parameters[age_index] = 1.0;

            if( nuc_info.age_multiple <= 1.0*PhysicalUnits::second )
              nuc_info.age_multiple = PeakDef::defaultDecayTime( nuc_nuclide, nullptr );
          
            if( nuc_info.age_multiple <= 1.0*PhysicalUnits::second )
            {
              nuc_info.age_multiple = 1.0;
              parameters[age_index] = nuc.age;
            }
          }else
          {
            assert( RelActCalcAuto::element(nuc.source) || RelActCalcAuto::reaction(nuc.source) );
            assert( nuc.age <= 0.0 );

            fit_age = false;  // This will make age a constant parameter below
            nuc_info.age_multiple = 1.0;
            parameters[age_index] = 0.0;
          }
          
          // If ages of nuclides of an element are set to all be the same, we will use the age slot
          //  of the first nuclide to control the age of all nuclides for that element; we'll check
          //  for this here, and if we have a nuclide of the same element preceding \c nuc_num, we'll
          //  fix the age here
          if( (nuc_num > 0) && nuc_nuclide && rel_eff_curve.nucs_of_el_same_age )
          {
            bool found_control_age = false;
            for( size_t i = 0; i < nuc_num; ++i )
            {
              const RelActCalcAuto::NucInputInfo &other_nuc = rel_eff_curve.nuclides[i];
              const SandiaDecay::Nuclide * const other_nuc_nuclide = RelActCalcAuto::nuclide(other_nuc.source);

              if( other_nuc_nuclide && (other_nuc_nuclide->atomicNumber == nuc_nuclide->atomicNumber) )
              {
                if( (other_nuc.age != nuc.age) || (other_nuc.fit_age != nuc.fit_age) )
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
            double half_life = nuc_nuclide->halfLife;
            if( rel_eff_curve.nucs_of_el_same_age )
            {
              for( size_t i = 0; i < rel_eff_curve.nuclides.size(); ++i )
              {
                const RelActCalcAuto::NucInputInfo &other_nuc = rel_eff_curve.nuclides[i];
                const SandiaDecay::Nuclide * const other_nuc_nuclide = RelActCalcAuto::nuclide(other_nuc.source);

                if( other_nuc_nuclide && (other_nuc_nuclide->atomicNumber == nuc_nuclide->atomicNumber) )
                  half_life = std::max(half_life, other_nuc_nuclide->halfLife);
              }
            }//if( options.nucs_of_el_same_age )
            
            double max_age = std::max( 5.0*nuc.age, 15.0*half_life );
            
            // We'll clamp the upper age, to the range where humans may have done the seperation
            // TODO: is there any problem with a nuclide >100 years, that isnt effectively infinite?
            max_age = std::min( max_age, 120*PhysicalUnits::year );

            if( nuc.fit_age_max.has_value() )
              max_age = nuc.fit_age_max.value();
            
            const double min_age = nuc.fit_age_min.has_value() ? nuc.fit_age_min.value() : 0.0;

            lower_bounds[age_index] = min_age / nuc_info.age_multiple;
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
            
            vector<SandiaDecay::EnergyRatePair> gammas;
            
            if( nuc_nuclide )
            {
              SandiaDecay::NuclideMixture mix;
              mix.addAgedNuclideByActivity( nuc_nuclide, ns_decay_act_mult, nuc.age );
              gammas = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
            }else if( nuc_element )
            {
              gammas.reserve( nuc_element->xrays.size() );
              for( const SandiaDecay::EnergyIntensityPair &xray : nuc_element->xrays )
                gammas.emplace_back( ns_decay_act_mult*xray.intensity, xray.energy );
            }else if( nuc_reaction )
            {
              gammas.reserve( nuc_reaction->gammas.size() );
              for( const ReactionGamma::Reaction::EnergyYield &rctn_gamma : nuc_reaction->gammas )
                gammas.emplace_back( ns_decay_act_mult*rctn_gamma.abundance, rctn_gamma.energy );
            }else
            {
              assert( 0 );
              throw runtime_error( "invalid src" );
            }
            
            
            for( const RelActCalcAuto::RoiRange &erange : energy_ranges )
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
            

            bool is_constrained = false;
            for( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &nuc_constraint : rel_eff_curve.act_ratio_constraints )
            {
              is_constrained = (nuc_constraint.constrained_source == nuc.source);
              if( is_constrained )
              {
                parameters[act_index] = -1.0;
                cout << "Setting par " << act_index << " to -1.0 for " << nuc.name() << endl;

                cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = -1.0;
                assert( std::find( constant_parameters.begin(), constant_parameters.end(), static_cast<int>(act_index) ) == constant_parameters.end() );
                constant_parameters.push_back( static_cast<int>(act_index) );
                break;
              }
            }//for( const auto &nuc_constraint : rel_eff_curve.act_ratio_constraints )

            for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints ) 
            {
              if( !nuc_nuclide )
                break;

              if( is_constrained )
              {
                assert( constraint.nuclide != nuc_nuclide );
                break;
              }

              assert( constraint.nuclide );
              is_constrained = (constraint.nuclide == nuc_nuclide);
              if( is_constrained )
              {
                parameters[act_index] = 1.0;
                cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = -1.0;
                cout << "Setting par " << act_index << " to 1.0 for " << nuc.name() << endl;

                if( constraint.lower_mass_fraction == constraint.upper_mass_fraction )
                {
                  cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = -1.0;
                  assert( std::find( constant_parameters.begin(), constant_parameters.end(), static_cast<int>(act_index) ) == constant_parameters.end() );
                  constant_parameters.push_back( static_cast<int>(act_index) );
                }else
                {
                  lower_bounds[act_index] = 0.5;
                  upper_bounds[act_index] = 1.5;
                }
                break;
              }//if( is_constrained ) 
            }//for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )

            if( is_constrained )
              continue;
              

            assert( cost_functor->m_nuclides[re_eff_index][nuc_num].source == nuc.source );

            double rel_act_mult = 1.0;
            if( nuc.starting_rel_act.has_value() )
            {
              rel_act_mult = nuc.starting_rel_act.value();
            }else if( !top_energy_to_rel_act.empty() )
            {
              // Use mid BR
              //const tuple<double,double,double> &representative_act = top_energy_to_rel_act[top_energy_to_rel_act.size()/2]; //{energy, br, rel. act.}
              //rel_act_mult = std::get<2>( representative_act );
              
              // Use highest BR
              const tuple<double,double,double> &representative_act = top_energy_to_rel_act.front(); //{energy, br, rel. act.}
              rel_act_mult = std::get<2>( representative_act );

              // For physical model, our efficiency may be way smaller than one, so we will correct
              //  for this; the other Rel Eff equations we will assume a starting flat line at 1.0.
              if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
              {
                //const size_t this_rel_eff_start = cost_functor->rel_eff_eqn_start_parameter(re_eff_index);
                //const RelActAutoCostFcn::PhysModelRelEqnDef<double> re_input
                //          = RelActAutoCostFcn::make_phys_eqn_input( rel_eff_curve, cost_functor->m_drf,
                //                                                   parameters, this_rel_eff_start );
                //const double rel_eff = RelActCalc::eval_physical_model_eqn_imp( energy, re_input.self_atten,
                //                                           re_input.external_attens, re_input.det.get(),
                //                                           re_input.hoerl_b, re_input.hoerl_c );
                const double energy = std::get<0>( representative_act );
                const double rel_eff = cost_functor->relative_eff( energy, re_eff_index, parameters );
                assert( !IsNan(rel_eff) && !IsInf(rel_eff) );
                if( !IsNan(rel_eff) && !IsInf(rel_eff) && (rel_eff > 0.0) )
                {
                  rel_act_mult /= rel_eff;
                }
              }//if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
            }else
            {
              rel_act_mult = 100*PhysicalUnits::bq;
              solution.m_warnings.push_back( "Could not estimate a starting activity for "
                                            + nuc.name() + ".  This may be because there are no"
                                            " significant gammas for the nuclide in the selected energy"
                                            " ranges." );
            }

                        
            if( nuc.min_rel_act.has_value() )
              rel_act_mult = std::max( rel_act_mult, nuc.min_rel_act.value() );

            if( nuc.max_rel_act.has_value() )
              rel_act_mult = std::min( rel_act_mult, nuc.max_rel_act.value() );

            parameters[act_index] = 1.0;
            cost_functor->m_nuclides[re_eff_index][nuc_num].activity_multiple = rel_act_mult;
            assert( rel_act_mult > 0.0 );

#ifndef NDEBUG
            {// begin debug input
              for( const auto &erange : energy_ranges )
              {
                double nuc_sum = 0.0;
                for( const SandiaDecay::EnergyRatePair &er : gammas )
                {
                  if( er.numPerSecond <= std::numeric_limits<float>::min() ) //1.17549e-38
                    continue;
                  
                  if( (er.energy < erange.lower_energy) || (er.energy > erange.upper_energy) )
                    continue;
                  
                  const double energy = er.energy;
                  const double rel_eff = cost_functor->relative_eff( energy, re_eff_index, parameters );
                  const double yield = er.numPerSecond / ns_decay_act_mult;
                  
                  nuc_sum += rel_eff * yield * foreground->live_time() * rel_act_mult;
                }
                
                cout << "For initial estimate: " << nuc.name() << " contributes " << nuc_sum << " gammas to datas "
                << foreground->gamma_integral(erange.lower_energy, erange.upper_energy)
                << " counts in range [" << erange.lower_energy << ", " << erange.upper_energy << "]" << endl;
              }//for( const auto &erange : energy_ranges )
              
              cout << "Setting initial relative activity for " << nuc.name()
                   << " to " << rel_act_mult << " bq" << endl;
            }// end debug input
#endif

            if( nuc.min_rel_act.has_value()
             && nuc.max_rel_act.has_value() 
             && (nuc.min_rel_act.value() == nuc.max_rel_act.value()) )
            {
              constant_parameters.push_back( static_cast<int>(act_index) );
            }else
            {
              lower_bounds[act_index] = (nuc.min_rel_act.has_value() ? nuc.min_rel_act.value() : 0.0) / rel_act_mult;
              if( nuc.max_rel_act.has_value() )
                upper_bounds[act_index] = nuc.max_rel_act.value() / rel_act_mult;
            }
          }//if( !succesfully_estimated_re_and_ra )
        }//for( size_t nuc_num = 0; nuc_num < rel_eff_curve.nuclides.size(); ++nuc_num )


        // If we have multiple mass-fraction constraints for this element, we need to make sure thier initial
        //  values sum to less than 1.0.
        // Sum starting mass fractions for all constraints in this Rel Eff curve
        map<short int,double> elements_starting_mass_fracs, elements_lower_mass_fracs, elements_fixed_mass_frac;
        for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )
        {
          const SandiaDecay::Nuclide *nuc = constraint.nuclide;
          const size_t act_index = cost_functor->nuclide_parameter_index( nuc, re_eff_index );
          const double starting_dist = parameters[act_index] - 0.5;

          const double lower_frac = constraint.lower_mass_fraction;
          const double upper_frac = constraint.upper_mass_fraction;
          const double starting_mass_frac = lower_frac + starting_dist*(upper_frac - lower_frac);

          elements_lower_mass_fracs[nuc->atomicNumber] += lower_frac;
          if( lower_frac == upper_frac )
            elements_fixed_mass_frac[nuc->atomicNumber] += lower_frac;

          elements_starting_mass_fracs[nuc->atomicNumber] += starting_mass_frac;
        }//for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )

        //Now go threw and check if any sum is greater than or equal to 1.0
        for( const map<short int,double>::value_type &an_sum : elements_starting_mass_fracs )
        {
          const short int atomic_number = an_sum.first;
          const double starting_mass_frac = an_sum.second;
          if( starting_mass_frac < 1.0 )
            continue;

          // If we're here, we have need to reduce the mass fractions for this element
          const double lower_allowed_frac = elements_lower_mass_fracs[atomic_number];

          assert( lower_allowed_frac < 1.0 ); //`RelEffCurveInput::check_nuclide_constraints()` should have already checked this
          if( lower_allowed_frac >= 1.0 )
            throw std::logic_error( "Mass fraction constraint sums to over 1 - please fix `RelEffCurveInput::check_nuclide_constraints()` to catch this." );

          const double fixed_fraction = elements_fixed_mass_frac[atomic_number];
          const double initial_variable_amount = starting_mass_frac - fixed_fraction;
          const double target_var_frac = 0.5*((1.0 - fixed_fraction) + lower_allowed_frac);  //halfway between the smallest and largest amount allowed
          const double amount_need_reduced = initial_variable_amount - target_var_frac;
          const double frac_variable_reduce = amount_need_reduced / (initial_variable_amount - lower_allowed_frac);
          const double updated_variable_frac = initial_variable_amount - amount_need_reduced;

          //Example starting mass fraction 1.25, with 0.5 fixed, and lowest variable amount of 0.1
          //  fixed_fraction = 0.5
          //  starting_mass_frac = 1.25
          //  lower_allowed_frac = 0.1
          //  initial_variable_amount = 1.25 - 0.5 = 0.75
          //  target_var_frac = 0.5*((1.0 - 0.5) + 0.1) = 0.3;
          //  amount_need_reduced = 0.75 - 0.3 = 0.45;
          //  frac_variable_reduce = 0.45 / (0.75 - 0.1) = 0.6923;
          //  updated_variable_frac = 0.75 - 0.45 = 0.3;
          //So then 0.5 fixed, 0.3 variable, and 0.2 "other"

          assert( initial_variable_amount > 0.0 );
          assert( amount_need_reduced >= 0.0 );
          assert( (frac_variable_reduce > 0.0) && (frac_variable_reduce < 1.0) );

          double check_var_frac = 0.0;
          for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &c : rel_eff_curve.mass_fraction_constraints )
          {
            if( (c.nuclide->atomicNumber != atomic_number) || (c.lower_mass_fraction == c.upper_mass_fraction) )
              continue;

            const size_t act_index = cost_functor->nuclide_parameter_index( c.nuclide, re_eff_index );
            parameters[act_index] = 0.5 + ((1.0 - frac_variable_reduce)*(parameters[act_index] - 0.5));
            assert( (parameters[act_index] >= 0.5) && (parameters[act_index] <= 1.5) );

            check_var_frac += (c.lower_mass_fraction + (parameters[act_index] - 0.5)*(c.upper_mass_fraction - c.lower_mass_fraction));
          }//for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff_curve.mass_fraction_constraints )

          assert( fabs(check_var_frac - updated_variable_frac) < 0.0001 );
        }//for( const map<short int,double>::value_type &an_sum : elements_starting_mass_fracs )

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
          
        bool fit_skew_energy_dep = PeakDef::is_energy_dependent( options.skew_type, ct );

        if( fit_skew_energy_dep )
        {
          //We'll check if stat sig peaks span at least 100 keV (arbitrarily chosen), and if not,
          // drop fitting the energy dependence for skew
          const double min_energy_range_for_skew_ene_dep = 100;
          shared_ptr<const PeakDef> lowest_peak, highest_peak;
          for( const shared_ptr<const PeakDef> &p : all_peaks )
          {
            for( const auto &r : energy_ranges )
            {
              if( (p->mean() >= r.lower_energy) && (p->mean() <= r.upper_energy) )
              {
                if( !lowest_peak || (p->mean() < lowest_peak->mean()) )
                  lowest_peak = p;

                if( !highest_peak || (p->mean() > highest_peak->mean()) )
                  highest_peak = p;

                break;
              }//if( we found the ROI for this peak )
            }//for( const auto &r : energy_ranges )
          }//for( const shared_ptr<const PeakDef> &p : all_peaks )

          if( !lowest_peak
             || !highest_peak
             || (highest_peak->mean() < (lowest_peak->mean() + min_energy_range_for_skew_ene_dep)))
          {
            fit_skew_energy_dep = false;
          }
        }//if( fit_skew_energy_dep )

        // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 1 ) //Fix the peak skew, and dont fit it
        assert( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall );
        if( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall )
        {
          fit_skew_energy_dep = false;
        }
#endif
        
        // If we will be fitting an energy dependence, make sure the cost functor knows this
        cost_functor->m_skew_has_energy_dependance |= fit_skew_energy_dep;

        const size_t skew_start = cost_functor->m_skew_par_start_index;
        const size_t num_skew_coefs = PeakDef::num_skew_parameters( options.skew_type );
        
        parameters[skew_start + i] = starting;
        lower_bounds[skew_start + i] = lower;
        upper_bounds[skew_start + i] = upper;
        
        // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 1 ) //Fix the peak skew, and dont fit it
#pragma message( "Doing peak skew hack - fixing peak skew to preset values!" )
        assert( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall );
        if( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall )
        {
          cerr << "\n\n\nDoing peak skew hack! - fixing peak skew to preset values\n\n\n" << endl;
          
          const double coefs[4] = { 2.7005, 16.0960, 2.1167, 4.1235 };
          parameters[skew_start + i] = starting = coefs[i];
          constant_parameters.push_back( static_cast<int>(skew_start + i) );
          lower_bounds[skew_start + i] = std::nullopt;
          upper_bounds[skew_start + i] = std::nullopt;
        }
#elif( PEAK_SKEW_HACK == 2 ) //First fit without peak skew, then fit with
#pragma message( "Doing peak skew hack! - first fitting without skew, then fitting skew." )
        if( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall )
        {
          cerr << "\n\n\nDoing peak skew hack! - first fitting without skew, then fitting skew.\n\n\n" << endl;
          const double coefs[4] = { 2.7005, 16.0960, 2.1167, 4.1235 };
          parameters[skew_start + i] = starting = coefs[i];
          constant_parameters.push_back( static_cast<int>(skew_start + i) );
          initial_const_parameters.push_back( static_cast<int>(skew_start + i) );
        }
#endif
        
          
        // Specify ranges for second set of skew parameters
        if( !fit_skew_energy_dep )
        {
          parameters[skew_start + i + num_skew_coefs] = -999.9;
          constant_parameters.push_back( static_cast<int>(skew_start + i + num_skew_coefs) );
        }else
        {
          parameters[skew_start + i + num_skew_coefs] = starting;
          lower_bounds[skew_start + i + num_skew_coefs] = lower;
          upper_bounds[skew_start + i + num_skew_coefs] = upper;
          
          // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 1 ) //Fix the peak skew, and dont fit it
          if( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall )
          {
            constant_parameters.push_back( static_cast<int>(skew_start + i + num_skew_coefs) );
            lower_bounds[skew_start + i + num_skew_coefs] = std::nullopt;
            upper_bounds[skew_start + i + num_skew_coefs] = std::nullopt;
          }
#elif( PEAK_SKEW_HACK == 2 ) //First fit without peak skew, then fit with
          if( options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall )
          {
            constant_parameters.push_back( static_cast<int>(skew_start + i + num_skew_coefs) );
            initial_const_parameters.push_back( static_cast<int>(skew_start + i + num_skew_coefs) );
          }
#endif
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
      
      parameters[amp_index] = 1.0;   //See cost_functor->m_free_peak_area_multiples[extra_peak_index] for the initial peak amplitude estimate, which is what this paramater will multiply to get the actual area
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
        parameters[fwhm_index] = -1.0; //so we can assert on
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
      parameters.resize( num_pars, sm_peak_range_uncert_offset );
      pars = &parameters[0];
                
      const double val_with_zero_yield = sm_peak_range_uncert_offset - sm_peak_range_uncert_par_scale/options.additional_br_uncert;

      lower_bounds.resize( num_pars, optional<double>(val_with_zero_yield) );
      upper_bounds.resize( num_pars );
      
      const size_t num_skew_coefs = PeakDef::num_skew_parameters( options.skew_type );
      cost_functor->m_add_br_uncert_start_index = cost_functor->m_skew_par_start_index + 2*num_skew_coefs;
      
      assert( num_pars == cost_functor->number_parameters() );
      assert( num_pars == (cost_functor->m_add_br_uncert_start_index + cost_functor->m_peak_ranges_with_uncert.size()) );
    }//if( options.additional_br_uncert > 0.0 )
    
    const size_t num_pars = cost_functor->number_parameters();
    
    ceres::CostFunction *cost_function = nullptr;
    
#if( PERFORM_DEVELOPER_CHECKS )
    //Test auto diff vs numerical diff
    auto test_gradients = [constant_parameters,&cost_functor]( RelActAutoCostFcn *fcn, const vector<double> &x ){
      ceres::DynamicAutoDiffCostFunction<RelActAutoCostFcn,32> auto_diff( fcn, ceres::DO_NOT_TAKE_OWNERSHIP );
      auto_diff.SetNumResiduals( static_cast<int>(fcn->number_residuals()) );
      auto_diff.AddParameterBlock( static_cast<int>(fcn->number_parameters()) );
      
      ceres::NumericDiffOptions num_diff_options;
      num_diff_options.relative_step_size = 1E-6; //Default is 1E-6.
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
        if( (diff > 1.0E-18) && (frac_diff > 1.0E-2) )  //This is totally arbitrary
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
    
    test_gradients( cost_functor.get(), parameters );
#endif

    if( sm_use_auto_diff )
    {
      ceres::DynamicAutoDiffCostFunction<RelActAutoCostFcn,sm_auto_diff_stride_size> *dyn_auto_diff_cost_function
            = new ceres::DynamicAutoDiffCostFunction<RelActAutoCostFcn,sm_auto_diff_stride_size>( cost_functor.get(),
                                                                      ceres::DO_NOT_TAKE_OWNERSHIP );
      dyn_auto_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
      dyn_auto_diff_cost_function->AddParameterBlock( static_cast<int>(cost_functor->number_parameters()) );

      cost_function = dyn_auto_diff_cost_function;
    }else
    {
      ceres::NumericDiffOptions num_diff_options;
      num_diff_options.relative_step_size = 1E-4; //Default is 1E-6
      ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn> *dyn_num_diff_cost_function
            = new ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn>( cost_functor.get(),
                                                ceres::DO_NOT_TAKE_OWNERSHIP, num_diff_options );
      dyn_num_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
      dyn_num_diff_cost_function->AddParameterBlock( static_cast<int>(cost_functor->number_parameters()) );

      cost_function = dyn_num_diff_cost_function;
    }//if( sm_use_auto_diff ) / else

    
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
      // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 2 )
      const auto pos = std::find( begin(constant_parameters), end(constant_parameters), static_cast<int>(i) );
      if( pos != end(constant_parameters) )
      {
        // Test that if the paramater is totally const (e.g., not just `initial_const_parameters`), then we
        //  shouldnt have put lower/upper bounds on the paramater
        const auto init_pos = std::find( begin(initial_const_parameters), end(initial_const_parameters), static_cast<int>(i) );
        assert( (init_pos != end(initial_const_parameters))
               || (!lower_bounds[i].has_value() && !upper_bounds[i].has_value()) );
        
        continue;
      }
#endif
      
      if( lower_bounds[i].has_value() )
        problem.SetParameterLowerBound(pars, static_cast<int>(i), *lower_bounds[i] );
      if( upper_bounds[i].has_value() )
        problem.SetParameterUpperBound(pars, static_cast<int>(i), *upper_bounds[i] );
    }//for( size_t i = 0; i < num_pars; ++i )
    
    
    // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 2 )
    if( !initial_const_parameters.empty() )
      cost_functor->m_skip_skew = true;
#endif
    
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
    size_t num_fit_par = 0;
    double par_area = 1.0;
    for( size_t i = 0; i < parameters.size(); ++i )
    {
      const auto const_pos = std::find( begin(constant_parameters), end(constant_parameters), static_cast<int>(i) );
      if( const_pos == end(constant_parameters) )
      {
        num_fit_par += 1;
        par_area *= (parameters[i] > 0.1) ? parameters[i] : 1.0;
        cout << "Starting value of " << cost_functor->parameter_name(i) << ": " << parameters[i] << endl;
      }
    }
    cout << "Starting with parameter volume of " << par_area << " from " << num_fit_par << " paramaters." << endl;
    // The below is pretty arbitrary - and only kinda sort optimized on one problem
    ceres_options.initial_trust_region_radius = 100.0*std::min( std::max( par_area, 10.0 ), 10.0*num_fit_par );

    ceres_options.max_trust_region_radius = 1e16;

    // Minimizer terminates when the trust region radius becomes smaller than this value.
    ceres_options.min_trust_region_radius = 1e-32;
    // Lower bound for the relative decrease before a step is accepted.
    ceres_options.min_relative_decrease = 1e-3;
    
    
    ceres_options.linear_solver_type = ceres::DENSE_QR; //ceres::DENSE_SCHUR, ceres::DENSE_NORMAL_CHOLESKY, ceres::ITERATIVE_SCHUR
    
    ceres_options.minimizer_progress_to_stdout = true; //true;
    ceres_options.logging_type = ceres::PER_MINIMIZER_ITERATION;
    ceres_options.max_num_iterations = 50000;
#ifndef NDEBUG
    ceres_options.max_solver_time_in_seconds = 600.0; //10 minutes!
#else
    ceres_options.max_solver_time_in_seconds = 120.0;
#endif

    // Changing function_tolerance from 1e-9 to 1e-7, and using initial_trust_region_radius=1E5, and use_nonmonotonic_steps=false
    //  - 52.88% enrich, 261 function calls, 64 iterations, Chi2=0.664167 (terminate due to function tol reached)
    // Decreasing from 1e-7 to 1e9 increases number of calls/time by about 30%, for one example problem, but solution is better.
    ceres_options.function_tolerance = 1e-9;

    ceres_options.gradient_tolerance = 1.0E-4*ceres_options.function_tolerance; // Per documentation of `gradient_tolerance`
    // parameter_tolerance seems to be what terminates the minimization on some example problems.
    //  It looks to be:
    //    |step|_2 <= parameter_tolerance * ( |x|_2 +  parameter_tolerance)
    //  Where `|step|_2` and `|x|_2` indicate the norm of the parameters - e.g.,
    ceres_options.parameter_tolerance = 1e-11; //Default value is 1e-8.  Using 1e-11 gets the cost change down to 0.1 in Chi2, for an example Physical Model
    
    // TODO: there are a ton of ceres::Solver::Options that might be useful for us to set

    std::unique_ptr<RelActAutoCostFcn::CheckCeresTerminateCallback> terminate_callback;
    if( cancel_calc )
    {
      terminate_callback.reset( new RelActAutoCostFcn::CheckCeresTerminateCallback(cancel_calc) );
      ceres_options.callbacks.push_back( terminate_callback.get() );
    }
    
    
    //RelActAutoCostFcn::CeresStepSummaryCallback step_summary_callback( cost_functor.get() );
    //ceres_options.callbacks.push_back( &step_summary_callback );
    
    // Setting ceres_options.num_threads >1 doesnt seem to do much (any?) good - so instead we do
    //  some multiple threaded computations in RelActAutoSolution::eval(...)
    ceres_options.num_threads = static_cast<int>( std::thread::hardware_concurrency(), 2 );

    cost_functor->m_solution_finished = false;

    ceres::Solver::Summary summary;
    ceres::Solve(ceres_options, &problem, &summary);
    //std::cout << summary.BriefReport() << "\n";
    
    cost_functor->m_solution_finished = true;
    
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
    
    
    // 20250324 HACK to test fitting peak skew
#if( PEAK_SKEW_HACK == 2 )
    if( (summary.termination_type == ceres::CONVERGENCE) || (summary.termination_type == ceres::USER_SUCCESS) )
    {
      // 20250324 HACK to test fitting peak skew
      if( !initial_const_parameters.empty() && cost_functor->m_skip_skew )
      {
        cout << "20250324 HACK to test fitting peak skew - replacing the manifold!" << endl;
        cost_functor->m_skip_skew = false;
        for( const int par : initial_const_parameters )
        {
          const auto pos = std::find( begin(constant_parameters), end(constant_parameters), par );
          assert( pos != end(constant_parameters) );
          if( pos != end(constant_parameters) )
            constant_parameters.erase(pos);
          
          if( lower_bounds[par].has_value() )
            problem.SetParameterLowerBound(pars, par, *lower_bounds[par] );
          if( upper_bounds[par].has_value() )
            problem.SetParameterUpperBound(pars, par, *upper_bounds[par] );
        }
        
        ceres::Manifold *subset_manifold = nullptr;
        if( !constant_parameters.empty() )
          subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_pars), constant_parameters );
        
        problem.SetManifold( pars, subset_manifold );
        
        cout << "20250324 HACK to test fitting peak skew - reset manifold - now refitting." << endl;
        
        ceres::Solve(ceres_options, &problem, &summary);
        
        cout << "20250324 HACK to test fitting peak skew - Done refitting." << endl;
        
        std::cout << summary.FullReport() << "\n";
      }//
    }//if( (summary.termination_type == ceres::CONVERGENCE) || (summary.termination_type == ceres::USER_SUCCESS) )
#endif //#if( PEAK_SKEW_HACK == 2 )
    
    
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
    // DENSE_SVD: accurate but slow (only to be used for small to moderate sized problems). Handles full-rank as well as rank deficient Jacobians.
    // SPARSE_QR: fast, but not capable of computing the covariance if the Jacobian is rank deficient
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD; //SPARSE_QR;
    cov_options.num_threads = ceres_options.num_threads;
    //cov_options.min_reciprocal_condition_number = 1e-14;
    //cov_options.column_pivot_threshold = -1;
    cov_options.null_space_rank = -1; //default: 0;
    
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
      solution.m_phys_units_cov.clear();
      solution.m_final_uncertainties.clear();
      
      if( !covariance.Compute(covariance_blocks, &problem) )
      {
        //Possible reason we're here: rank deficient Jacobian is encountered
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
          throw runtime_error( "Failed to get covariance matrix." );
        
        solution.m_covariance.resize( num_pars, vector<double>(num_pars, 0.0) );
        
        for( size_t row = 0; row < num_pars; ++row )
        {
          for( size_t col = 0; col < num_pars; ++col )
            solution.m_covariance[row][col] = row_major_covariance[row*num_pars + col];
        }//for( size_t row = 0; row < num_pars; ++row )
        
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
        
        solution.m_rel_eff_covariance.resize( num_rel_eff_curves );
        solution.m_rel_act_covariance.resize( num_rel_eff_curves );
        for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
        {
          const size_t rel_eff_start = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
          const size_t num_rel_eff_par = cost_functor->rel_eff_eqn_num_parameters(rel_eff_index);
          vector<vector<double>> &cov = solution.m_rel_eff_covariance[rel_eff_index];
          
          get_cov_block( rel_eff_start, num_rel_eff_par, cov );
          
          for( size_t par_index = 0; par_index < num_rel_eff_par; ++par_index )
          {
            const double scale = cost_functor->parameter_scale_factor(rel_eff_start + par_index);
            for( size_t row = 0; row < num_rel_eff_par; ++row )
              cov[row][par_index] *= scale;
            for( size_t col = 0; col < num_rel_eff_par; ++col )
              cov[par_index][col] *= scale;
          }
        
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
    solution.m_parameter_scale_factors.resize( parameters.size(), 1.0 );
    
    solution.m_parameter_names.resize( parameters.size() );
    solution.m_parameter_were_fit.resize( parameters.size(), true );
    for( size_t i = 0; i < parameters.size(); ++i )
    {
      solution.m_parameter_names[i] = cost_functor->parameter_name( i );
      solution.m_parameter_were_fit[i] = std::find( begin(constant_parameters), end(constant_parameters), static_cast<int>(i) ) == end(constant_parameters);

      solution.m_parameter_scale_factors[i] = cost_functor->parameter_scale_factor( i );
      
      assert( (solution.m_parameter_scale_factors[i] > 0.0)
             || !solution.m_parameter_were_fit[i]
             || cost_functor->mass_constraint_multiple(i,parameters).has_value() );
    }//for( size_t i = 0; i < parameters.size(); ++i )

    if( !solution.m_covariance.empty() )
    {
      solution.m_phys_units_cov = solution.m_covariance;

      for( size_t row = 0; row < num_pars; ++row )
      {
        const optional<double> row_mass_mult = cost_functor->mass_constraint_multiple(row,parameters);
        const double row_mult = row_mass_mult.has_value() ? *row_mass_mult : solution.m_parameter_scale_factors[row];

        for( size_t col = 0; col < num_pars; ++col )
        {
          const optional<double> col_mass_mult = cost_functor->mass_constraint_multiple(col,parameters);
          const double col_mult = col_mass_mult.has_value() ? *col_mass_mult : solution.m_parameter_scale_factors[row];

          solution.m_phys_units_cov[row][col] *= row_mult*col_mult;
        }//for( size_t col = 0; col < num_pars; ++col )
      }//for( size_t row = 0; row < num_pars; ++row )
    }//if( !solution.m_covariance.empty() )

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
    vector<vector<PeakDef>> fit_peaks_for_each_curve( (num_rel_eff_curves > 1) ? num_rel_eff_curves : size_t(0) );
    
    for( const RoiRangeChannels &range : cost_functor->m_energy_ranges )
    {
      const PeaksForEnergyRange these_peaks = cost_functor->peaks_for_energy_range( range, parameters, {} );

      fit_peaks.insert( end(fit_peaks), begin(these_peaks.peaks), end(these_peaks.peaks) );
      
      if( num_rel_eff_curves > 1 )
      {
        shared_ptr<PeakContinuum> shared_continuum;
        for( size_t i = 0; i < num_rel_eff_curves; ++i )
        {
          PeaksForEnergyRange this_re_peaks = cost_functor->peaks_for_energy_range( range, parameters, {i} );
          
          //Set the peaks for this ROI to share the same PeakContinuum as for this ROI and other R.E. curves
          for( PeakDef &p : this_re_peaks.peaks )
          {
            if( !shared_continuum )
              shared_continuum = p.continuum();
            p.setContinuum( shared_continuum );
          }//for( PeakDef &p : this_re_peaks.peaks )
          
          vector<PeakDef> &re_peaks = fit_peaks_for_each_curve[i];
          re_peaks.insert( end(re_peaks), begin(this_re_peaks.peaks), end(this_re_peaks.peaks) );
        }//for( size_t i = 0; i < num_rel_eff_curves; ++i )
      }//if( more than R.E. curve )
    }//for( const RoiRangeChannels &range : cost_functor->m_energy_ranges )
    
    std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMean );
    for( vector<PeakDef> &re_peaks : fit_peaks_for_each_curve )
      std::sort( begin(re_peaks), end(re_peaks), &PeakDef::lessThanByMean );
    
    solution.m_fit_peaks_in_spectrums_cal = fit_peaks;
    if( num_rel_eff_curves == 1 )
      solution.m_fit_peaks_in_spectrums_cal_for_each_curve = vector<vector<PeakDef>>{1, fit_peaks};
    else
      solution.m_fit_peaks_in_spectrums_cal_for_each_curve = fit_peaks_for_each_curve;
    
    // \c fit_peaks are in the original energy calibration of the spectrum, we may need to adjust
    //  them to match the new energy calibration
    if( new_cal != cost_functor->m_energy_cal )
    {
      auto adjust_peaks = [cost_functor, new_cal]( vector<PeakDef> &input_peaks ){
        deque<shared_ptr<const PeakDef>> tmp_peaks;
        for( const auto &p : input_peaks )
          tmp_peaks.push_back( make_shared<const PeakDef>( p ) );
        
        auto adjusted_peaks = EnergyCal::translatePeaksForCalibrationChange( tmp_peaks,
                                                                            cost_functor->m_energy_cal, new_cal );
        input_peaks.clear();
        for( const auto &p : adjusted_peaks )
          input_peaks.push_back( *p );
      };//adjust_peaks lambda
      
      adjust_peaks( fit_peaks );
      
      for( vector<PeakDef> &re_peaks : fit_peaks_for_each_curve )
        adjust_peaks( re_peaks );
    }//if( new_cal != cost_functor->m_energy_cal )
    
    solution.m_fit_peaks = fit_peaks;
    if( num_rel_eff_curves == 1 )
      solution.m_fit_peaks_for_each_curve = vector<vector<PeakDef>>{1, fit_peaks};
    else
      solution.m_fit_peaks_for_each_curve = fit_peaks_for_each_curve;
    
    assert( solution.m_fit_peaks_for_each_curve.size() == num_rel_eff_curves );
    assert( solution.m_fit_peaks_in_spectrums_cal_for_each_curve.size() == num_rel_eff_curves );
    
    assert( cost_functor->m_nuclides.size() == num_rel_eff_curves );
    assert( cost_functor->m_options.rel_eff_curves.size() == num_rel_eff_curves );

    solution.m_rel_eff_forms.clear();
    solution.m_rel_eff_forms.resize( num_rel_eff_curves, RelActCalc::RelEffEqnForm::LnX );

    solution.m_rel_eff_coefficients.clear();
    solution.m_rel_eff_coefficients.resize( num_rel_eff_curves );

    solution.m_rel_activities.clear();
    solution.m_rel_activities.resize( num_rel_eff_curves );
    
    for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    {
      const vector<NucInputGamma> &input_nuclides = cost_functor->m_nuclides[rel_eff_index];
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      
      const size_t rel_eff_start = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
      const size_t num_rel_eff_par = cost_functor->rel_eff_eqn_num_parameters(rel_eff_index);
      const auto rel_eff_iter = begin(parameters) + rel_eff_start;
      solution.m_rel_eff_coefficients[rel_eff_index] = cost_functor->pars_for_rel_eff_curve( rel_eff_index, parameters );
      solution.m_rel_eff_forms[rel_eff_index] = rel_eff_curve.rel_eff_eqn_type;

      for( size_t act_index = 0; act_index < input_nuclides.size(); ++act_index )
      {
        const NucInputGamma &nuc_input = input_nuclides[act_index];
        const RelActCalcAuto::NucInputInfo &src = rel_eff_curve.nuclides[act_index];
        assert( nuc_input.source == src.source );
      
        RelActCalcAuto::NuclideRelAct nuc_output;
        nuc_output.source = nuc_input.source;
        nuc_output.age = cost_functor->age( nuc_input, rel_eff_index, parameters );
        nuc_output.age_was_fit = nuc_input.fit_age;
        nuc_output.rel_activity = cost_functor->relative_activity( nuc_input.source, rel_eff_index, parameters );
      
        nuc_output.age_uncertainty = cost_functor->age( nuc_input, rel_eff_index, uncertainties );
        
        bool is_mass_constrained = false;
        
        if( const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(nuc_output.source) )
        {
          for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint mass_cons : rel_eff_curve.mass_fraction_constraints )
            is_mass_constrained |= (mass_cons.nuclide == nuc);
        }

        if( is_mass_constrained )
        {
          // We'll multiple the uncertainty of mass-fraction paramater, by the derivative of RelAct

          nuc_output.rel_activity_uncertainty = -1;
          if( sm_use_auto_diff )
          {
            try
            {
              const size_t rel_act_index = cost_functor->nuclide_parameter_index( nuc_output.source, rel_eff_index );
              vector<ceres::Jet<double,sm_auto_diff_stride_size>> input_jets( begin(parameters), end(parameters) );
              input_jets[rel_act_index].v[0] = 1.0; //It doesnt matter which element of `v` we use to get the derivative, so we'll just use the first one.
              ceres::Jet<double,sm_auto_diff_stride_size> rel_act_jet = cost_functor->relative_activity(nuc_output.source, rel_eff_index, input_jets);
              assert( (nuc_output.rel_activity < 1.0E-12) || (fabs(nuc_output.rel_activity - rel_act_jet.a) < 1.0E-6*nuc_output.rel_activity) );
              const double derivative = rel_act_jet.v[0];
              nuc_output.rel_activity_uncertainty = uncertainties[rel_act_index] * derivative;
            }catch( std::exception & )
            {

            }//try / catch
          }else
          {
            //We wont bother implementing this, since we never plan to use numerical differentiation
#pragma message( "Calculating rel-act uncert when not using auto-diff is not implementing")
          }

        }else
        {
          nuc_output.rel_activity_uncertainty = cost_functor->relative_activity( nuc_input.source, rel_eff_index, uncertainties );
        }

        vector<SandiaDecay::EnergyRatePair> gammas;
        if( const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(nuc_output.source) )
        {
          SandiaDecay::NuclideMixture mix;
          mix.addAgedNuclideByActivity( nuc, ns_decay_act_mult, nuc_output.age );
          gammas = mix.photons( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
        }else if( const SandiaDecay::Element *el = RelActCalcAuto::element(nuc_output.source) )
        {
          for( const SandiaDecay::EnergyIntensityPair &xray : el->xrays )
            gammas.emplace_back( ns_decay_act_mult*xray.intensity, xray.energy );
        }else if( const ReactionGamma::Reaction *rctn = RelActCalcAuto::reaction(nuc_output.source) )
        {
          for( const ReactionGamma::Reaction::EnergyYield &reaction : rctn->gammas )
            gammas.emplace_back( ns_decay_act_mult*reaction.abundance, reaction.energy );
        }else
        {
          assert( 0 );
        }
        
        nuc_output.gamma_energy_br.clear();
        for( const auto &gamma : gammas )
        {
          const double yield = gamma.numPerSecond / ns_decay_act_mult;
          nuc_output.gamma_energy_br.push_back( {gamma.energy, yield} );
        }
      
        if( IsNan(nuc_output.age_uncertainty) )
        {
          const size_t par_act_index = cost_functor->nuclide_parameter_index( src.source, rel_eff_index );
          const size_t par_age_index = par_act_index + 1;
          solution.m_warnings.push_back( "Variance for nuclide " + nuc_input.name()
                                      + " age, is invalid ("
                                      + std::to_string(uncerts_squared[par_age_index]) + ")" );
        }
      
        if( IsNan(nuc_output.rel_activity_uncertainty) )
        {
          const size_t par_act_index = cost_functor->nuclide_parameter_index( src.source, rel_eff_index );
          solution.m_warnings.push_back( "Variance for nuclide " + nuc_input.name()
                                      + " activity, is invalid ("
                                      + std::to_string(uncerts_squared[par_act_index]) + ")" );
        }

        solution.m_rel_activities[rel_eff_index].push_back( nuc_output );
      }//for( size_t act_index = 0; act_index < cost_functor->m_nuclides.size(); ++ act_index )
    }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    
    
    assert( solution.m_rel_activities.size() == num_rel_eff_curves );
    if( solution.m_rel_activities.size() != num_rel_eff_curves )
      throw logic_error( "size(rel_eff_curves) != size(solution.m_rel_activities)" );
    
    // If we want to correct for Pu242, we wont alter solution.m_rel_activities, but place the
    //  corrected Pu mass fractions in solution.m_corrected_pu
    solution.m_corrected_pu.clear();
    solution.m_corrected_pu.resize( num_rel_eff_curves );
    for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      const vector<RelActCalcAuto::NuclideRelAct> &rel_acts = solution.m_rel_activities[rel_eff_index];
      
      if( rel_eff_curve.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
      {
        try
        {
          RelActCalc::Pu242ByCorrelationInput raw_rel_masses;

          set<double> pu_ages;
          double pu_total_mass = 0.0, raw_rel_mass = 0.0;
          for( const RelActCalcAuto::NuclideRelAct &nuc : rel_acts )
          {
            assert( !RelActCalcAuto::is_null(nuc.source) );
            const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(nuc.source);
            if( !nuclide || (nuclide->atomicNumber != 94) )
              continue;
            
            const double rel_mass = nuc.rel_activity / nuclide->activityPerGram();
            
            pu_ages.insert( nuc.age );
            pu_total_mass += rel_mass;
            switch( nuclide->massNumber )
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
            }//switch( nuclide->massNumber )

            if( pu_ages.size() == 0 )
            {
              assert( pu_total_mass == 0.0 );
            }else if( pu_ages.size() == 1 )
            {
              raw_rel_masses.pu_age = *begin(pu_ages);
            }else
            {
              const vector<double> ages( begin(pu_ages), end(pu_ages) );
              raw_rel_masses.pu_age = ages[ages.size() / 2]; //just take the median age...
            }
          }//for( const NuclideRelAct &nuc : m_rel_activities )
          
          // We dont have to divide by `pu_total_mass`, but we will, just for debuging.
          raw_rel_mass /= pu_total_mass;
          raw_rel_masses.pu238_rel_mass /= pu_total_mass;
          raw_rel_masses.pu239_rel_mass /= pu_total_mass;
          raw_rel_masses.pu240_rel_mass /= pu_total_mass;
          raw_rel_masses.pu241_rel_mass /= pu_total_mass;
          raw_rel_masses.other_pu_mass  /= pu_total_mass;
          
          const RelActCalc::Pu242ByCorrelationOutput corr_output
          = RelActCalc::correct_pu_mass_fractions_for_pu242( raw_rel_masses,
                                                            rel_eff_curve.pu242_correlation_method );
          
          if( !corr_output.is_within_range )
            solution.m_warnings.push_back( "The fit Pu enrichment is outside range validated in the"
                                          " literature for the Pu242 correlation." );
          
          solution.m_corrected_pu[rel_eff_index].reset( new RelActCalc::Pu242ByCorrelationOutput(corr_output) );
        }catch( std::exception &e )
        {
          solution.m_warnings.push_back( "Correcting for Pu242 content failed: " + string(e.what()) );
        }//try / catch
      }//if( options.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
    }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    
    
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


    assert( cost_functor->m_free_peak_area_multiples.size() == cost_functor->m_options.floating_peaks.size() );
    
    for( size_t i = 0; i < cost_functor->m_options.floating_peaks.size(); ++i )
    {
      const size_t amp_index = cost_functor->m_free_peak_par_start_index + 2*i + 0;
      const size_t fwhm_index = amp_index + 1;

      const double area_multiple = cost_functor->m_free_peak_area_multiples[i];
      
      RelActCalcAuto::FloatingPeakResult peak;
      peak.energy = cost_functor->m_options.floating_peaks[i].energy;
      peak.amplitude = parameters[amp_index] * area_multiple;
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
      
      peak.amplitude_uncert = uncertainties[amp_index] * area_multiple;
      
      
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
    solution.m_phys_model_results.resize( num_rel_eff_curves );
    
    size_t first_phys_model_index = num_rel_eff_curves;
    const RelActCalcAuto::RelEffCurveInput *first_phys_model_rel_eff_curve = nullptr;

    for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      const size_t this_rel_eff_start_index = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
      const vector<double> &rel_eff_coefficients = solution.m_rel_eff_coefficients[rel_eff_index];
      
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        assert( rel_eff_coefficients.size() == (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2) );

        if( !first_phys_model_rel_eff_curve )
        {
          first_phys_model_index = rel_eff_index;
          first_phys_model_rel_eff_curve = &rel_eff_curve;
        }

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
          
          phys_model_result.hoerl_b = (rel_eff_coefficients[b_index] - RelActCalc::ns_decay_hoerl_b_offset) * RelActCalc::ns_decay_hoerl_b_multiple;
          phys_model_result.hoerl_c = (rel_eff_coefficients[c_index] - RelActCalc::ns_decay_hoerl_c_offset) * RelActCalc::ns_decay_hoerl_c_multiple;
          
          if( uncertainties.size() > (this_rel_eff_start_index + b_index) )
            phys_model_result.hoerl_b_uncert = uncertainties[this_rel_eff_start_index + b_index] * RelActCalc::ns_decay_hoerl_b_multiple;

          if( uncertainties.size() > (this_rel_eff_start_index + c_index) )
            phys_model_result.hoerl_c_uncert = uncertainties[this_rel_eff_start_index + c_index] * RelActCalc::ns_decay_hoerl_c_multiple; 

          // If this isnt the first physical model, and we are using the same Hoerl function for all
          //  then we will use the Hoerl function from the first physical model.
          if( first_phys_model_rel_eff_curve && options.same_hoerl_for_all_rel_eff_curves && (rel_eff_index != first_phys_model_index) )
          {
            assert( phys_model_result.hoerl_b == solution.m_phys_model_results[first_phys_model_index]->hoerl_b );
            assert( phys_model_result.hoerl_c == solution.m_phys_model_results[first_phys_model_index]->hoerl_c );
            
            phys_model_result.hoerl_b_uncert = solution.m_phys_model_results[first_phys_model_index]->hoerl_b_uncert;
            phys_model_result.hoerl_c_uncert = solution.m_phys_model_results[first_phys_model_index]->hoerl_c_uncert;
          }
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
        const double multiple = 1.0 + (solution.m_final_parameters[par_index] - sm_peak_range_uncert_offset)/sm_peak_range_uncert_par_scale;
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
    
    
#ifndef NDEBUG
    // Now that have the solution filled out, lets do a few sanity checks
    for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      const size_t this_rel_eff_start_index = cost_functor->rel_eff_eqn_start_parameter(rel_eff_index);
      const vector<double> &rel_eff_coefficients = solution.m_rel_eff_coefficients[rel_eff_index];
      
      const vector<RelActCalcAuto::NuclideRelAct> &solution_rel_acts = solution.m_rel_activities[rel_eff_index];
      const vector<NucInputGamma> &input_nuclides = cost_functor->m_nuclides[rel_eff_index];
      assert( input_nuclides.size() == solution_rel_acts.size() );
      
      for( size_t act_index = 0; act_index < input_nuclides.size(); ++act_index )
      {
        const NucInputGamma &input_src = input_nuclides[act_index];
        const RelActCalcAuto::NuclideRelAct &solution_src = solution_rel_acts[act_index];
        assert( input_src.source == solution_src.source );
        
        const size_t solution_index = solution.fit_parameters_index_for_source( solution_src.source, rel_eff_index );
        const size_t implementation_index = cost_functor->nuclide_parameter_index( input_src.source, rel_eff_index );
        if( solution_index != implementation_index )
          cout << "solution_index != implementation_index - " << solution_index << " != " << implementation_index << endl;
        assert( solution_index == implementation_index );
      }//for( size_t act_index = 0; act_index < input_nuclides.size(); ++act_index )
    }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
#endif
    
    //solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
    
    return solution;
  }//RelActCalcAuto::RelActAutoSolution solve_ceres( ceres::Problem &problem )
  

  /** Returns the decay gammas along with their rate, and the transition that led to them (i.e.,
   where in the decay chain they come from).
   
   Note that the returned result may have multiple entries with the exact same energy (i.e., from
   different decays).  Also, that gammas with zero amplitude will also be returned, so this way no
   matter the age, the same number of gammas will be returned, in the same order.
   */
  std::shared_ptr<const vector<NucInputGamma::EnergyYield>> decay_gammas( const RelActCalcAuto::NucInputInfo &nuc_info,
                                          double age,
                                          const std::vector<double> &gammas_to_exclude ) const
  {
    assert( RelActCalcAuto::nuclide(nuc_info.source) || (nuc_info.age < 0.0) );
    assert( !RelActCalcAuto::nuclide(nuc_info.source) || (nuc_info.age >= 0.0) );
    
    assert( !RelActCalcAuto::is_null(nuc_info.source) );
        
    const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(nuc_info.source);
    const SandiaDecay::Element * const element = RelActCalcAuto::element(nuc_info.source);
    const ReactionGamma::Reaction * const reaction = RelActCalcAuto::reaction(nuc_info.source);

    assert( nuclide || (age <= 0.0) );
    age = nuclide ? age : 0.0; //for non nuclides, make sure age is 0.0

    const void *key_ptr = static_cast<const void *>( nuclide );
    if( !key_ptr )
      key_ptr = static_cast<const void *>( element );
    if( !key_ptr )
      key_ptr = static_cast<const void *>( reaction );
    if( !key_ptr )
      throw std::logic_error( "Invalid NucInputInfo" );
    
    
    const pair<const void *,float> key{ key_ptr, static_cast<float>(age) };
    
    // First check if we have already computed this
    {//begin lock on m_aged_gammas_cache_mutex
      std::lock_guard<std::mutex> cache_lock( m_aged_gammas_cache_mutex );
          
      const auto pos = m_aged_gammas_cache.find( key );
      if( pos != end(m_aged_gammas_cache) )
      {
        assert( pos->second );
        return pos->second;
      }
    }//end lock on m_aged_gammas_cache_mutex

    // We havent cached what we are being asked for, so we will compute it and cache the result
    shared_ptr<vector<NucInputGamma::EnergyYield>> answer = make_shared<vector<NucInputGamma::EnergyYield>>();

    //Taking arguments by reference for this next cleanup block because the cleanup variable
    // variable (addToCacheOnExit) will be destructed before any of the ones already defined
    // at this point.
    DoWorkOnDestruct addToCacheOnExit( [this, key, answer]() { 
      std::lock_guard<std::mutex> cache_lock( m_aged_gammas_cache_mutex );
      assert( answer );
      m_aged_gammas_cache[key] = answer;
            
      // TODO: limit size of cache
    } );

    if( element )
    {
      const vector<SandiaDecay::EnergyIntensityPair> &xrays = element->xrays;
      
      for( const SandiaDecay::EnergyIntensityPair &xray : xrays )
      {
        NucInputGamma::EnergyYield info;
        info.energy = xray.energy;
        info.yield = xray.intensity;
        info.transition_index = 0;
        info.transition = nullptr;
        info.gamma_type = PeakDef::SourceGammaType::XrayGamma;
        info.element = element;
        answer->push_back( std::move(info) );
      }//for( const SandiaDecay::EnergyIntensityPair &xray : xrays )

      return answer;
    }else if( reaction )
    {
      const vector<ReactionGamma::Reaction::EnergyYield> &gammas = reaction->gammas;
      
      for( const ReactionGamma::Reaction::EnergyYield &reaction_yield : gammas )
      {
        NucInputGamma::EnergyYield info;
        info.energy = reaction_yield.energy;
        info.yield = reaction_yield.abundance;
        info.gamma_type = PeakDef::SourceGammaType::NormalGamma;
        info.reaction = reaction;
        answer->push_back( std::move(info) );
      }//for( const ReactionGamma::EnergyYield &reaction : reactions )

      return answer;
    }//if( nuc_info.element ) / else if( nuc_info.reaction )
    
    const SandiaDecay::Nuclide * const parent = nuclide;
    assert( parent );
    
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( parent, ns_decay_act_mult, age );
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mix.activity( 0.0 );
    
    // all_gamma_transitions may contain duplicate energies - we will combine these below
    vector<NucInputGamma::EnergyYield> all_gamma_transitions;
    all_gamma_transitions.reserve( 1536 ); //avoid a few resizes; for Pu241, at nominal age, the actual number of entries is 1229
    
    bool has_annihilation = false;
    NucInputGamma::EnergyYield annihilationInfo;
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
      
      const size_t n_decaysToChildren = nuclide->decaysToChildren.size();
      
      for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
      {
        const SandiaDecay::Transition * const transition = nuclide->decaysToChildren[decayIndex];
        const size_t n_products = transition->products.size();
        
        for( size_t productNum = 0; productNum < n_products; ++productNum )
        {
          const SandiaDecay::RadParticle &particle = transition->products[productNum];
          if( (particle.type == SandiaDecay::ProductType::GammaParticle)
              || (particle.type == SandiaDecay::ProductType::XrayParticle) )
          {
            bool exclude = false;
            for( const double exclude_energy : gammas_to_exclude )
              exclude = (exclude || (fabs(particle.energy - exclude_energy) < 0.001));
            
            if( exclude )
              continue;
            
            all_gamma_transitions.push_back( {} );
            NucInputGamma::EnergyYield &info = all_gamma_transitions.back();
            info.energy = particle.energy;
            info.yield = activity * particle.intensity * transition->branchRatio;
            
            info.transition_index = productNum;
            info.transition = transition;
            info.gamma_type = (particle.type == SandiaDecay::ProductType::GammaParticle) 
                              ? PeakDef::SourceGammaType::NormalGamma
                              : PeakDef::SourceGammaType::XrayGamma;
          }else if( particle.type == SandiaDecay::ProductType::PositronParticle )
          {
            has_annihilation = true;
            annihilationInfo.yield += 2.0 * activity * particle.intensity * transition->branchRatio;
           
            num_annih_trans += 1;
            annihilationInfo.transition_index = productNum;
            annihilationInfo.transition = transition;
          }//if( particle.type is gamma ) / else if( position )
        }//for( size_t productNum = 0; productNum < n_products; ++productNum )
      }//for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
    }//for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
    
    if( has_annihilation )
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
    
    
    // Now go throw and combine duplicate energies
    answer->reserve( all_gamma_transitions.size() );
    answer->push_back( all_gamma_transitions.front() );
    size_t energy_start_index = 0;
    for( size_t current_index = 1; current_index < all_gamma_transitions.size(); ++current_index )
    {
      const NucInputGamma::EnergyYield &current = all_gamma_transitions[current_index];
      NucInputGamma::EnergyYield &prev = answer->back();
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
        answer->push_back( current );
      }
    }//for( size_t current_index = 1; current_index < all_gamma_transitions.size(); ++current_index )

    assert( answer );
    return answer;
  }//std::shared_ptr<const vector<EnergyYield>> decay_gammas( const SandiaDecay::Nuclide * const parent, ... )

  
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
        assert( nuc_input.nominal_gammas );
        for( const NucInputGamma::EnergyYield &line_info : *nuc_input.nominal_gammas )
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
          
          const double rel_act = relative_activity( nuc_input.source, rel_eff_index, x );
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
    
#ifndef NDEBUG
    // Make sure the scale factors are 1.0 for all FWHM parameters that are not FWHM parameters
    for( size_t i = 0; i < num_drf_par; ++i )
    {
      const double scale_factor = this->parameter_scale_factor( i + RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars );
      assert( scale_factor == 1.0 );
    }
#endif

    return eval_fwhm( energy, m_options.fwhm_form, drf_start, num_drf_par, m_drf );
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
    
#if( PEAK_SKEW_HACK == 2 )
    if( (m_options.skew_type == PeakDef::SkewType::NoSkew)
       || m_skip_skew ) // 20250324 HACK to test fitting peak skew
    {
      peak.setSkewType( PeakDef::SkewType::NoSkew );
      return;
    }
#endif
    
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
          assert( upper_en_val > -500.0 );//should NOT have value -999.9
          
          const T val = lower_en_val + mean_frac*(upper_en_val - lower_en_val);
          peak.set_coefficient( val, ct );
        }else
        {
          assert( x[skew_start + num_skew + i] < -500.0 ); //should have value -999.9, unless being numerically varies
          const T val = x[skew_start + i];
          peak.set_coefficient( val, ct );
        }
      }
    }else
    {
      for( size_t i = 0; i < skew_pars.size(); ++i )
      {
        assert( x[skew_start + num_skew + i] < -500.0 ); //should have value -999.9
        const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + i);
        const T val = x[skew_start + i];
        peak.set_coefficient( val, ct );
      }
    }//if( m_skew_has_energy_dependance )
  }//set_peak_skew(...)
  
  /// Returns the index for this nuclide within the relative efficiency curve
  size_t nuclide_index( const RelActCalcAuto::SrcVariant &src, const size_t rel_eff_index ) const
  {
    assert( !RelActCalcAuto::is_null(src) );
    
    if( RelActCalcAuto::is_null(src) )
      throw runtime_error( "nuclide_index: null src" );
    
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "nuclide_index: Invalid relative efficiency curve index." );
    
    const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    size_t nuc_index = rel_eff_curve.nuclides.size();
    for( size_t i = 0; i < nuc_index; ++i )
    {
      if( rel_eff_curve.nuclides[i].source == src )
        nuc_index = i;
    }
    
    assert( nuc_index < rel_eff_curve.nuclides.size() );
    if( nuc_index >= rel_eff_curve.nuclides.size() )
      throw std::logic_error( "Invalid nuclide." );
    
    return nuc_index;
  }//size_t nuclide_index( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index )
  

  /** Returns the index for this nuclide's activity (the age is +1) within the Ceres parameter vector */
  size_t nuclide_parameter_index( const RelActCalcAuto::SrcVariant &src, const size_t rel_eff_index ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." );
    
    size_t index = m_acts_par_start_index;
    for( size_t i = 0; i < rel_eff_index; ++i )
      index += 2*m_options.rel_eff_curves[i].nuclides.size();
  
    return index + 2*nuclide_index(src, rel_eff_index);
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
          num_coefs_this_rel_eff = 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2;

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
      
        assert( !RelActCalcAuto::is_null(rel_eff_curve.nuclides[act_num].source) );
        if( (act_index % 2) == 0 )
        {
          //Check if mass fraction constraint for nuclide, and if so
          const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(rel_eff_curve.nuclides[act_num].source);
          const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint *constraint = nullptr;
          if( nuc )
          {
            const auto start_iter = begin(rel_eff_curve.mass_fraction_constraints);
            const auto end_iter = end(rel_eff_curve.mass_fraction_constraints);
            const auto pos = std::find_if( start_iter, end_iter, [nuc]( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &mfc ){
              return mfc.nuclide == nuc;
            } );
            if( pos != end_iter )
              constraint = &(*pos);
          }

          if( constraint )
            return "MFrac" + re_ind + "(" + rel_eff_curve.nuclides[act_num].name() + ")";

          return "Act" + re_ind + "(" + rel_eff_curve.nuclides[act_num].name() + ")";
        }
        return "Age" + re_ind + "(" + rel_eff_curve.nuclides[act_num].name() + ")";
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
  

  double parameter_scale_factor( const size_t index ) const
  {
    assert( RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars == m_fwhm_par_start_index );
    
    if( index < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars )
    {
      // The offset, gain, and quad adjustments are `(par/sm_energy_par_offset - 1.0) * sm_energy_cal_multiple`
      //  This value is then added to the base paramater value.
      return RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset;
    }//if( index < RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars )
    
    
    if( index < m_rel_eff_par_start_index )
    {
      // Its a FWHM parameter
      return 1.0;  // No scale factor for FWHM parameters currently
    }//if( index < m_rel_eff_par_start_index )
    
    
    if( index < m_acts_par_start_index )
    {
      size_t coef_num = m_rel_eff_par_start_index;

      const size_t num_rel_eff_curves = m_options.rel_eff_curves.size();
      size_t first_physical_model_index = num_rel_eff_curves;
      for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
      {
        const auto &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
        size_t num_coefs_this_rel_eff = rel_eff_curve.rel_eff_eqn_order + 1;
        if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          num_coefs_this_rel_eff = 2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2;

        const size_t sub_ind = index - coef_num;
        if( sub_ind >= num_coefs_this_rel_eff )
        {
          coef_num += num_coefs_this_rel_eff;
          continue;
        }

        const string re_ind = num_rel_eff_curves > 1 ? std::to_string(rel_eff_index) : "";
        if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        {
          if( first_physical_model_index == num_rel_eff_curves )
            first_physical_model_index = rel_eff_index;

          if( sub_ind == 0 )
            return RelActCalc::ns_an_ceres_mult; //Atomic number
          if( sub_ind == 1 )
            return 1.0; // We just use 1.0 for areal density right now
        
          if( sub_ind < (2 + 2*rel_eff_curve.phys_model_external_atten.size()) )
          {
            const size_t ext_atten_num = (sub_ind - 2) / 2;
            if( (sub_ind % 2) == 0 )
              return RelActCalc::ns_an_ceres_mult;
            return 1.0; //AD
          }

          assert( sub_ind >= (2 + 2*rel_eff_curve.phys_model_external_atten.size()) );

          if( !rel_eff_curve.phys_model_use_hoerl )
            return 0.0;  // We arent using the Hoerl function, so no scale factor

          if( m_options.same_hoerl_for_all_rel_eff_curves && (rel_eff_index != first_physical_model_index) )
            return 0.0;  // We arent using the Hoerl function, so no scale factor

          const size_t hoerl_num = sub_ind - 2 - 2*rel_eff_curve.phys_model_external_atten.size();
          if( hoerl_num == 0 )
            return rel_eff_curve.phys_model_use_hoerl ? RelActCalc::ns_decay_hoerl_b_multiple : 0.0;
        
          assert( hoerl_num == 1 );
          if( hoerl_num == 1 )
            return rel_eff_curve.phys_model_use_hoerl ? RelActCalc::ns_decay_hoerl_c_multiple : 0.0;
        }else
        {
          return 1.0; // The polynomial-esque Rel. Eff. eqn coefficients are not scaled
        }//if( physical model ) / else
      }//for( const auto &rel_eff_curve : m_options.rel_eff_curves )

      assert( 0 );
      throw std::logic_error( "Logic for determining Physical Model coefficient name is bad." );
    }//if( index < m_acts_par_start_index )
    
    if( index < m_free_peak_par_start_index )
    {
      // This is a nuclide activity or age parameter; each source gets two parameters, activity, and age
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
      
        const NucInputGamma &src = m_nuclides[rel_eff_index][act_num];
        assert( !RelActCalcAuto::is_null(rel_eff_curve.nuclides[act_num].source) );
        if( (act_index % 2) == 0 )
          return (src.activity_multiple > 0.0) ? src.activity_multiple : 0.0;
        return (src.age_multiple > 0.0) ? src.age_multiple : 0.0;
      }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
    }//if( index < m_free_peak_par_start_index )
    
    if( index < m_skew_par_start_index )
    {
      // This is a free peak parameter; each peak gets one parameter - its amplitude.
      assert( m_free_peak_area_multiples.size() == m_options.floating_peaks.size() );
      if( m_free_peak_area_multiples.size() != m_options.floating_peaks.size() )
        throw std::logic_error( "Bad computing of parameter scale factor - free peak area multiples size mismatch" );  //shouldnt ever happen

      const size_t index_within_free_peaks = index - m_free_peak_par_start_index;
      const size_t free_peak_num = index_within_free_peaks / 2;
      assert( free_peak_num < m_options.floating_peaks.size() );
      if( free_peak_num >= m_options.floating_peaks.size() )
        throw std::logic_error( "Bad computing of parameter scale factor - too large of index" );
      
      if( (index_within_free_peaks % 2) == 0 )  // We want the amplitude.
        return m_free_peak_area_multiples[free_peak_num];

      if( !m_options.floating_peaks[free_peak_num].release_fwhm ) //If not fitting FWHM, return 0.0.
        return 0.0;
      return 1.0;
    }//if( index < m_skew_par_start_index )
    
    if( index < m_add_br_uncert_start_index )
    {
      return 1.0; // No scale factor for skew parameters currently
    }//if( index < m_add_br_uncert_start_index )
    
    assert( m_add_br_uncert_start_index != std::numeric_limits<size_t>::max() );
    assert( (m_add_br_uncert_start_index + m_peak_ranges_with_uncert.size()) == number_parameters() );
    
    if( m_add_br_uncert_start_index == std::numeric_limits<size_t>::max() )
      throw runtime_error( "Whack computing of scale - things not legit." );
    
    if( index >= (m_add_br_uncert_start_index + m_peak_ranges_with_uncert.size()) )
      throw runtime_error( "Bad computing of parameter name - too large of index" );
    
    //(par_value - sm_peak_range_uncert_offset)/sm_peak_range_uncert_par_scale;
    return 1.0 / sm_peak_range_uncert_par_scale;
  }//double parameter_scale_factor( const size_t index ) const
  

  template<typename T>
  T relative_activity( const RelActCalcAuto::SrcVariant &src, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( x.size() == number_parameters() );
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    assert( (rel_eff_index < m_nuclides.size()) && (rel_eff_index < m_options.rel_eff_curves.size()) );
    if( rel_eff_index >= m_nuclides.size() )
      throw std::logic_error( "relative_activity: invalid relative efficiency curve index." );
    
    assert( !RelActCalcAuto::is_null(src) );
    if( RelActCalcAuto::is_null(src) )
      throw runtime_error( "relative_activity: invalid source." );

    // Check for a source constraint
    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];

#ifndef NDEBUG
    {
      const size_t nuc_x_index = nuclide_parameter_index( src, rel_eff_index );
      const size_t nuc_index = nuclide_index( src, rel_eff_index ); // Will throw exception if not a nuclide in the rel eff curve
      const double act_multiple = parameter_scale_factor( nuc_x_index );
      const double age_multiple = parameter_scale_factor( nuc_x_index + 1 );
      assert( ((m_nuclides[rel_eff_index][nuc_index].age_multiple < 0.0) && (age_multiple == 0.0)) 
              || (age_multiple == m_nuclides[rel_eff_index][nuc_index].age_multiple) );
      assert( ((m_nuclides[rel_eff_index][nuc_index].activity_multiple < 0.0) && (act_multiple == 0.0))
              || (act_multiple == m_nuclides[rel_eff_index][nuc_index].activity_multiple) );
    }
#endif

    for( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &nuc_constraint : rel_eff_curve.act_ratio_constraints )
    {
      if( nuc_constraint.constrained_source == src )
      {
#ifndef NDEBUG
        const size_t nuc_x_index = nuclide_parameter_index( src, rel_eff_index );
        const size_t nuc_index = nuclide_index( src, rel_eff_index );
        assert( m_nuclides[rel_eff_index][nuc_index].source == src );
        assert( m_nuclides[rel_eff_index][nuc_index].activity_multiple == -1.0 );
        const T rel_act_x_val = x[nuc_x_index];
        // rel_act_x_val should be fixed to -1.0, but numeric differentiator may still vary a little, 
        //  or if we pass in paramater uncertainties, then the value will be 0.0.
        assert( rel_act_x_val <= 0.0 ); 
#endif

        const double act_ratio = nuc_constraint.constrained_to_controlled_activity_ratio;
        assert( act_ratio > 0.0 );
        const RelActCalcAuto::SrcVariant &control_nuc = nuc_constraint.controlling_source;
        assert( !RelActCalcAuto::is_null(control_nuc) );
        const T controlled_act = relative_activity( control_nuc, rel_eff_index, x );
        return act_ratio * controlled_act;
      }
    }//for( const auto &nuc_constraint : rel_eff_curve.act_ratio_constraints )

    // Check for a mass fraction constraint (only applicable for nuclide sources)
    const SandiaDecay::Nuclide *src_nuc = RelActCalcAuto::nuclide(src);
    if( src_nuc )
    {
      for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &mass_frac_constraint : rel_eff_curve.mass_fraction_constraints )
      {
        if( mass_frac_constraint.nuclide != src_nuc )
          continue;
        
        // If we are here, `src.nuclide` is a mass-constrained nuclide.

        // Sum the relative masses of the other nuclides of this element
        // and sum the mass-constrained portion of this element.
        T sum_unconstrained_rel_mass_of_el( 0.0 ); //Note this is rel act divide by specific activity, and does not add up to one
        T sum_constrained_frac_rel_mass_of_el( 0.0 ); // This will include `src.nuclide`, and be less than 1.0
        for( const RelActCalcAuto::NucInputInfo &nuclide : rel_eff_curve.nuclides )
        {
          const SandiaDecay::Nuclide *nuclide_nuc = RelActCalcAuto::nuclide(nuclide.source);
          if( !nuclide_nuc || (nuclide_nuc->atomicNumber != src_nuc->atomicNumber) )
            continue;
          
          const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint *mass_constraint = nullptr;  
          for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &other_constraint : rel_eff_curve.mass_fraction_constraints )
          {
            if( other_constraint.nuclide == nuclide_nuc )
            {
              mass_constraint = &other_constraint;
              break;
            }
          }
              
          if( mass_constraint )
          {
            const size_t nuc_x_index = nuclide_parameter_index( nuclide.source, rel_eff_index );
            const T rel_dist = (x[nuc_x_index] - 0.5); //Rel Act paramater is constrained within 0.5 and 1.5, to make mass fraction go between lower and upper
            assert( rel_dist >= (0.0 - 1.0E-6) );
            assert( rel_dist <= (1.0 + 1.0E-6) );
            
            const T rel_mass = (1.0 - rel_dist)*mass_constraint->lower_mass_fraction + rel_dist*mass_constraint->upper_mass_fraction;
            sum_constrained_frac_rel_mass_of_el += rel_mass;
          }else
          {
            const T rel_act = relative_activity( nuclide.source, rel_eff_index, x );
            const double specific_activity = nuclide_nuc->activityPerGram();
            sum_unconstrained_rel_mass_of_el += (rel_act / specific_activity);
          }//if( is_mass_constrained ) / else
        }//for( const NucInputInfo &nuclide : rel_eff_curve.nuclides )

        // If there are multiple mass fraction constraints on the same nuclide, and a particular Ceres parameter
        //  solution gives the sum of all the mass fractions to be greater than 1.0, we throw an exception, causing
        //  this particular set of parameters to be rejected.
        //  Also, right now we are requiring at least one nuclide for the element to not have a mass-fraction constraint,
        //  which we may relax in the future (its just a little easier to implement this way, because otherwise we would
        //  have to scale the rel_act of them all, which we could do...)
        if( sum_constrained_frac_rel_mass_of_el > 1.0 )
          throw runtime_error( "Sum of constrained mass fractions of element is greater than 1.0." );

        if( sum_constrained_frac_rel_mass_of_el < 0.0 )
          throw runtime_error( "Sum of constrained mass fractions of element is less than 0.0." );

        const T unconstrained_rel_mass_frac_of_el = 1.0 - sum_constrained_frac_rel_mass_of_el;
        

        const size_t this_nuc_x_index = nuclide_parameter_index( src, rel_eff_index );
        const T rel_dist = x[this_nuc_x_index] - 0.5;
        assert( (rel_dist >= (0.0 - 1.0E-6)) || (x[this_nuc_x_index] == 0.0) );
        assert( rel_dist <= (1.0 + 1.0E-6) || (x[this_nuc_x_index] == 0.0) );
        const T this_rel_mass_frac = mass_frac_constraint.lower_mass_fraction
                                     + rel_dist*(mass_frac_constraint.upper_mass_fraction - mass_frac_constraint.lower_mass_fraction);
        assert( this_rel_mass_frac >= (0.0 - 1.0E-6) );
        assert( this_rel_mass_frac <= (1.0 + 1.0E-6) );

        const T total_rel_mass = sum_unconstrained_rel_mass_of_el / unconstrained_rel_mass_frac_of_el;
        const T this_rel_mass = total_rel_mass * this_rel_mass_frac;
        const T this_rel_act = this_rel_mass * mass_frac_constraint.nuclide->activityPerGram();

        return this_rel_act;
      }//for( const auto &mass_frac_constraint : rel_eff_curve.mass_fraction_constraints )
    }//if( src.nuclide )
    
    // No nuclide constraint, so we can just use the normal logic
    const size_t parent_nuc_par_index = nuclide_parameter_index( src, rel_eff_index );
    const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];

    const size_t nuc_index = nuclide_index( src, rel_eff_index );
    const NucInputGamma &nuc_info = nuclides[nuc_index];
    assert( nuc_info.source == src );
    
    assert( nuc_info.activity_multiple > 0.0 );

    return nuc_info.activity_multiple * x[parent_nuc_par_index];
  }//double relative_activity(...)
  

  template<typename T>
  T mass_enrichment_fraction( const RelActCalcAuto::SrcVariant &src, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( x.size() == number_parameters() );
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    assert( (rel_eff_index < m_nuclides.size()) && (rel_eff_index < m_options.rel_eff_curves.size()) );
    if( rel_eff_index >= m_nuclides.size() )
      throw std::logic_error( "mass_enrichment_fraction: invalid relative efficiency curve index." );
    
      assert( m_solution_finished );
    if( !m_solution_finished )
      throw std::logic_error( "mass_enrichment_fraction: should only be called after minimiztion has been completed." );
    
    assert( !RelActCalcAuto::is_null(src) );
    if( RelActCalcAuto::is_null(src) )
      throw runtime_error( "mass_enrichment_fraction: invalid source." );

    const SandiaDecay::Nuclide *nuclide = RelActCalcAuto::nuclide(src);
    if( !nuclide )
      throw runtime_error( "mass_enrichment_fraction: source is not a nuclide." );

    // Check for a source constraint
    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];

    // Check for mass fraction constraints
    for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &mass_frac_constraint : rel_eff_curve.mass_fraction_constraints )
    {
      if( mass_frac_constraint.nuclide == nuclide )
      {
        // For mass-constrained nuclides, return the mass fraction directly
        const size_t nuc_x_index = nuclide_parameter_index( src, rel_eff_index );
        const T rel_dist = x[nuc_x_index] - T(0.5);
        const T mass_frac = mass_frac_constraint.lower_mass_fraction
                           + rel_dist*(mass_frac_constraint.upper_mass_fraction - mass_frac_constraint.lower_mass_fraction);
        return mass_frac;
      }
    }//for( const auto &mass_frac_constraint : rel_eff_curve.mass_fraction_constraints )
    
    // No constraints, calculate normal mass enrichment fraction
    T sum_rel_mass = T(0.0);
    T nuc_rel_mass = T(-1.0);
    
    for( const RelActCalcAuto::NucInputInfo &nuclide_info : rel_eff_curve.nuclides )
    {
      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(nuclide_info.source);
      if( !nuc || (nuc->atomicNumber != nuclide->atomicNumber) )
        continue;
      
      const T rel_act = relative_activity( nuclide_info.source, rel_eff_index, x );
      const T rel_mass = rel_act / nuc->activityPerGram();
      
      if( nuclide_info.source == src )
        nuc_rel_mass = rel_mass;
        
      sum_rel_mass += (std::max)( rel_mass, T(0.0) );
    }
    
    if( nuc_rel_mass < T(0.0) )
      return T(0.0);
    
    
    
    
    return nuc_rel_mass / sum_rel_mass;
  }//T mass_enrichment_fraction(...)


  /** If the passed in index is the RelativeActivity index of a mass-constrained nuclide, will return the
   multiple from the parameter value, to the relative activity value; this is useful for convertering the uncertainties
   of the parameter into RelativeActivity uncertainty, or correlations between relative activiies.

   @param cost_functor The cost functor used to solve the problem.
   @param index The index, in `x` parametes, that you want the multiple of.
   @param parameters The solved, or current, parameters of the problem; the returned value will depend on the current parameters.
   @returns If `index` is Rel. Act. index of a nuclide that is mass-constrained, returns `d{RelAct}/d{parameter[index]}`, otherwise returns `nullopt`
  */
  std::optional<double> mass_constraint_multiple( const size_t index, const vector<double> &parameters ) const
  {
    const RelActAutoCostFcn *cost_functor = this;
    assert( index < parameters.size() );
    // First check if `index` refers to the relative activity of a source.
    if( (index < cost_functor->m_acts_par_start_index)
       || (index >= cost_functor->m_free_peak_par_start_index) )
    {
      return std::nullopt;
    }

    // Get the rel eff curve index and source `index` refers to
    for( size_t rel_eff_index = 0; rel_eff_index < cost_functor->m_options.rel_eff_curves.size(); ++rel_eff_index )
    {
      const size_t relact_begin = cost_functor->rel_act_start_parameter( rel_eff_index );
      const size_t relact_end = relact_begin + cost_functor->rel_act_num_parameters( rel_eff_index );

      if( (index < relact_begin) || (index >= relact_end) )
        continue;

      if( ((index - relact_begin) % 2) != 0 )
        return std::nullopt;

      const RelActCalcAuto::RelEffCurveInput &rel_effs = cost_functor->m_options.rel_eff_curves[rel_eff_index];
      if( rel_effs.mass_fraction_constraints.empty() )
        return std::nullopt; //common case; skip looping over sources below

      for( size_t src_index = 0; src_index < rel_effs.nuclides.size(); src_index += 1 )
      {
        const RelActCalcAuto::SrcVariant &src = rel_effs.nuclides[src_index].source;
        const size_t src_par_index = cost_functor->nuclide_parameter_index( src, rel_eff_index );
        if( index != src_par_index )
          continue;

        const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(src);
        if( !nuc )
          return std::nullopt; //non-nucs dont have mass-fraction constraints.

        for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_effs.mass_fraction_constraints )
        {
          if( constraint.nuclide != nuc )
            continue;

          // If we are here, we have a constraint for this source.
          //Rel Act paramater is constrained within 0.5 and 1.5, to make mass fraction go between lower and upper:
          // i.e.: rel_mass = (1.0 - rel_dist)*mass_constraint->lower_mass_fraction + rel_dist*mass_constraint->upper_mass_fraction;
          //  So we will multiple by d{RelAct}/d{Par[index]} to convert to the multiple for Rel. Act.
          double derivative = std::numeric_limits<double>::infinity();
          const double rel_act = cost_functor->relative_activity(src, rel_eff_index, parameters);
          if( sm_use_auto_diff )
          {
            vector<ceres::Jet<double,sm_auto_diff_stride_size>> input_jets( begin(parameters), end(parameters) );
            input_jets[index].v[0] = 1.0; //It doesnt matter which element of `v` we use to get the derivative, so we'll just use the first one.
            ceres::Jet<double,sm_auto_diff_stride_size> rel_act_jet = cost_functor->relative_activity(src, rel_eff_index, input_jets);
            assert( (rel_act < 1.0E-12) || (fabs(rel_act - rel_act_jet.a) < 1.0E-8*rel_act) );
            derivative = rel_act_jet.v[0];
          }else
          {
            // We dont really need this, as we always use auto diff - but will leave this in here to potentually use as a check in this future.
            const double step_size = 1.0E-3;
            vector<double> local_pars = parameters;
            if( ((parameters[index] - step_size) > 0.5) && ((parameters[index] + step_size) < 1.5) )
            {
              local_pars[index] = parameters[index] - step_size;
              const double lower_act = cost_functor->relative_activity(src, rel_eff_index, local_pars);
              local_pars[index] = parameters[index] + step_size;
              const double upper_act = cost_functor->relative_activity(src, rel_eff_index, local_pars);
              derivative = (upper_act - lower_act) / (2*step_size);
            }else if( (parameters[index] - step_size) > 0.5 )
            {
              local_pars[index] = parameters[index] - step_size;
              const double lower_act = cost_functor->relative_activity(src, rel_eff_index, local_pars);
              derivative = (rel_act - lower_act) / step_size;
            }else
            {
              local_pars[index] = parameters[index] + step_size;
              const double upper_act = cost_functor->relative_activity(src, rel_eff_index, local_pars);
              derivative = (upper_act - rel_act) / step_size;
            }
          }//if( we can use auto-diff ) / else numeric diff


          return derivative;
        }//for( const RelEffCurveInput::MassFractionConstraint &constraint : rel_effs.mass_fraction_constraints )

        return std::nullopt; //Didnt have the constraint
      }//for( size_t src_index = 0; src_index < rel_effs.nuclides.size(); src_index += 1 )
    }//for( loop over rel eff curves )

    assert( 0 );
    throw std::logic_error( "Failed to find expected source to check for constraint" );
    return std::nullopt;
  };//mass_constraint_multiple lambda

  //template<typename T>
  //T relative_activity( const RelActCalcAuto::SrcVariant &src, const size_t rel_eff_index, const std::vector<T> &x ) const
  //{
  //  return relative_activity( input_src_info(src, rel_eff_index), rel_eff_index, x );
  //}
  
  /** Converts between SandiaDecay::Nuclide, SandiaDecay::Element, or Reaction to its NucInputInfo. */
  const RelActCalcAuto::NucInputInfo &input_src_info( const RelActCalcAuto::SrcVariant &src, const size_t rel_eff_index ) const
  {
    assert( !RelActCalcAuto::is_null(src) );
    if( RelActCalcAuto::is_null(src) )
      throw runtime_error( "input_src_info: null src requested." );

    assert( rel_eff_index < m_nuclides.size() );
    if( rel_eff_index >= m_nuclides.size() )
      throw runtime_error( "Invalid relative efficiency curve index." ); 

    const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];
    for( size_t i = 0; i < nuclides.size(); ++i )
    {
      if( nuclides[i].source == src )
        return nuclides[i];
    }

    throw runtime_error( "input_src_info: src not found in nuclide list." );
  }//input_src_info(...)
  
  /** Returns the nuclide that is responsible for setting the passed in nuclides age.
   Will return the input nuclide if it controls its own age.
   */
  const SandiaDecay::Nuclide *age_controlling_nuc( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index ) const
  {
    assert( nuc );
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "Invalid relative efficiency curve index." ); 
    
    assert( m_nuclides.size() == m_options.rel_eff_curves.size() );
    if( m_nuclides.size() != m_options.rel_eff_curves.size() )
      throw std::logic_error( "Invalid number of nuclides." );

    const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];
    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    if( !rel_eff_curve.nucs_of_el_same_age )
      return nuc;
    
    const RelActCalcAuto::NucInputInfo &src = input_src_info( nuc, rel_eff_index );
    
    const size_t nuc_index = nuclide_index( src.source, rel_eff_index );
    assert( nuc_index < nuclides.size() );
    if( nuc_index >= nuclides.size() )
      throw std::logic_error( "Invalid nuclides index." );

    const NucInputGamma &input_nuc = nuclides[nuc_index];
    assert( RelActCalcAuto::nuclide(input_nuc.source) == nuc );
    
    for( size_t i = 0; i < nuc_index; ++i )
    {
      const NucInputGamma &nuc_info = nuclides[i];
      const SandiaDecay::Nuclide * const nuc_info_nuc = RelActCalcAuto::nuclide(nuc_info.source);
      if( nuc_info_nuc && (nuc_info_nuc->atomicNumber == nuc->atomicNumber) )
        return nuc_info_nuc;
    }
    
    return nuc;
  }//age_controlling_nuc(...)
  
  
  template<typename T>
  T age( const RelActCalcAuto::NucInputInfo &src, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( rel_eff_index < m_nuclides.size() );
    if( rel_eff_index >= m_nuclides.size() )
      throw std::logic_error( "age: invalid releff passed in." );
   
   assert( !RelActCalcAuto::is_null(src.source) );
   if( RelActCalcAuto::is_null(src.source) )
     throw runtime_error( "age: invalid source." );
   
    const SandiaDecay::Nuclide * const src_nuc = RelActCalcAuto::nuclide(src.source);
    if( !src_nuc )
    {
      if constexpr ( !std::is_same_v<T, double> )
        assert( x[nuclide_parameter_index(src.source, rel_eff_index) + 1] <= 0.0 );
      else
        assert( x[nuclide_parameter_index(src.source, rel_eff_index) + 1] <= 1.0E-6 );
      
      return T(0.0);
    }
    
    const size_t act_start_index = RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars
                                   + num_parameters(m_options.fwhm_form)
                                   + num_total_rel_eff_coefs();
    assert( act_start_index == m_acts_par_start_index );
    
    // The `parent_nuc` will often times be `nuc`.
    const SandiaDecay::Nuclide *parent_nuc = age_controlling_nuc( src_nuc, rel_eff_index );
    
    const size_t parent_par_index = nuclide_parameter_index( parent_nuc, rel_eff_index );
    const size_t par_nuc_index = nuclide_index( parent_nuc, rel_eff_index );

    const T age_scaler = x[parent_par_index + 1];
    const T age = age_scaler * m_nuclides[rel_eff_index][par_nuc_index].age_multiple;

    // If doing numerical differentiation, we may get negative ages, so we wont strickly check things
    //assert( age >= static_cast<double>(-std::numeric_limits<float>::epsilon()) );
    //if( age < static_cast<double>(-std::numeric_limits<float>::epsilon()) )
    //  throw runtime_error( "Negative age for " + m_nuclides[rel_eff_index][par_nuc_index].name() + " found." );
    
    return (max)( age, T(0.0) );
  }//age( nuclide )
  
  
  bool is_fixed_age( const SandiaDecay::Nuclide * const nuc, const size_t rel_eff_index ) const
  {
    assert( nuc );
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
      
      const T &b = coeffs[b_index];
      const T &c = coeffs[c_index];
      
      answer.hoerl_b = (b - RelActCalc::ns_decay_hoerl_b_offset) * RelActCalc::ns_decay_hoerl_b_multiple;  //(energy/1000)^b
      answer.hoerl_c = (c - RelActCalc::ns_decay_hoerl_c_offset) * RelActCalc::ns_decay_hoerl_c_multiple;  //c^(1000/energy)
    }//if( (b != 0.0) || (c != 1.0) )
    
    return answer;
  }//static PhysModelRelEqnDef make_phys_eqn_input(...)
  

  template<typename T>
  PhysModelRelEqnDef<T> make_phys_eqn_input( const size_t rel_eff_curve_index, const std::vector<T> &x ) const 
  {
    if( rel_eff_curve_index >= m_options.rel_eff_curves.size() )
      throw std::logic_error( "make_phys_eqn_input: rel eff curve index " + to_string(rel_eff_curve_index) + " is out of range" );
    
    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_curve_index];
    
    assert( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel );
    if( rel_eff_curve.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      throw std::logic_error( "make_phys_eqn_input: rel eff curve index " + to_string(rel_eff_curve_index) + " is not a physical model" );
    
    const size_t rel_eff_start = rel_eff_eqn_start_parameter(rel_eff_curve_index);

    if( (!m_options.same_hoerl_for_all_rel_eff_curves && !m_options.same_external_shielding_for_all_rel_eff_curves) )
      return make_phys_eqn_input( rel_eff_curve, m_drf, x, rel_eff_start );

    const vector<T> coefs = pars_for_rel_eff_curve( rel_eff_curve_index, x );

    return make_phys_eqn_input( rel_eff_curve, m_drf, coefs, 0 );
  }//make_phys_eqn_input(...)
  

  /** Returns the parameters for the specified relative efficiency curve.
   Note: scale factors have not been applied, as `make_phys_eqn_input(...)` will apply them, or for non-physical
   models, no scale factors are used.
   */
  template<typename T>
  vector<T> pars_for_rel_eff_curve( const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "pars_for_rel_eff_curve(): invalid rel eff curve index" );

    const size_t rel_eff_start_index = rel_eff_eqn_start_parameter(rel_eff_index);
    const size_t num_rel_eff_par = rel_eff_eqn_num_parameters(rel_eff_index);
    assert( (rel_eff_start_index + num_rel_eff_par) < x.size() );
    if( (rel_eff_start_index + num_rel_eff_par) >= x.size() )
      throw std::logic_error( "pars_for_rel_eff_curve(): whack number of rel eff parameters" );

    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];

    vector<T> coefs( begin(x) + rel_eff_start_index, begin(x) + rel_eff_start_index + num_rel_eff_par );

    if( (rel_eff_curve.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel)
      || (!m_options.same_hoerl_for_all_rel_eff_curves && !m_options.same_external_shielding_for_all_rel_eff_curves) )
    {
      // The polynomial-esque Rel. Eff. eqn coefficients are not scaled
      //  and `make_phys_eqn_input(...)` applies the scale factors for physical model.
      // For non-physical model, we'll test first/last paramater scale is 1.0.
      assert( (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel)
             || (parameter_scale_factor(rel_eff_start_index) == 1.0) );
      assert( (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel)
             || (parameter_scale_factor(rel_eff_start_index + num_rel_eff_par - 1) == 1.0) );
      
      return coefs;
    }

    // Find first Physical Model, and adjust Hoerl B and C paramaters.
    const RelActCalcAuto::RelEffCurveInput *first_phys_model = nullptr;
    size_t first_phys_model_par_start_index = std::numeric_limits<size_t>::max();
    
    for( size_t i = 0; i < m_options.rel_eff_curves.size(); ++i )
    {
      const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[i];
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        if( rel_eff_index == i )
          return coefs;

        first_phys_model = &rel_eff_curve;
        first_phys_model_par_start_index = rel_eff_eqn_start_parameter(i);
        break;
      }
    }//for( size_t i = 0; i < m_options.rel_eff_curves.size(); ++i )

    assert( first_phys_model );
    if( !first_phys_model )
      throw std::logic_error( "make_phys_eqn_input: no physical model found" );

    assert( (2 + 2*first_phys_model->phys_model_external_atten.size() + 2) <= coefs.size() );
    if( (2 + 2*first_phys_model->phys_model_external_atten.size() + 2) > coefs.size() )
      throw std::logic_error( "make_phys_eqn_input: inconsistent number of rel eff parameters - shouldnt have happened - programming logic error." ); //shouldnt happen

    if( m_options.same_hoerl_for_all_rel_eff_curves )
    {
      const size_t first_hoerl_b_index = first_phys_model_par_start_index + 2 + 2*first_phys_model->phys_model_external_atten.size();
      const size_t first_hoerl_c_index = first_hoerl_b_index + 1;
      assert( first_hoerl_c_index < x.size() );
      if( first_hoerl_c_index >= x.size() )
        throw std::logic_error( "make_phys_eqn_input: first hoerl indices out of range" ); //shouldnt happen

      const size_t this_hoerl_b_coefs_index = 2 + 2*rel_eff_curve.phys_model_external_atten.size();  
      const size_t this_hoerl_c_coefs_index = this_hoerl_b_coefs_index + 1;
      assert( this_hoerl_c_coefs_index < coefs.size() );
      if( this_hoerl_c_coefs_index >= coefs.size() )
        throw std::logic_error( "make_phys_eqn_input: this hoerl indices out of range" ); //shouldnt happen

      assert( coefs[this_hoerl_b_coefs_index] == -1.0 );
      assert( coefs[this_hoerl_c_coefs_index] == -1.0 );

      coefs[this_hoerl_b_coefs_index] = x[first_hoerl_b_index];
      coefs[this_hoerl_c_coefs_index] = x[first_hoerl_c_index];
    }//if( m_options.same_hoerl_for_all_rel_eff_curves )

    if( m_options.same_external_shielding_for_all_rel_eff_curves )
    {
      for( size_t i = 0; i < first_phys_model->phys_model_external_atten.size(); ++i )
      {
        const size_t atten_index = first_phys_model_par_start_index + 2 + 2*i;
        coefs[2 + 2*i + 0] = x[atten_index + 0];  //AN
        coefs[2 + 2*i + 1] = x[atten_index + 1];  //AD
      }
    }//if( m_options.same_external_shielding_for_all_rel_eff_curves )

    return coefs;
  }//pars_for_rel_eff_curve( const size_t rel_eff_index, const std::vector<T> &x )


  template<typename T>
  T relative_eff( const double energy, const size_t rel_eff_index, const std::vector<T> &x ) const
  {
    assert( rel_eff_index < m_options.rel_eff_curves.size() );
    if( rel_eff_index >= m_options.rel_eff_curves.size() )
      throw runtime_error( "relative_eff(): invalid rel eff curve index" );
    
    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    
    const size_t rel_eff_start_index = rel_eff_eqn_start_parameter(rel_eff_index);
    const size_t num_rel_eff_par = rel_eff_eqn_num_parameters(rel_eff_index);
    
    assert( (rel_eff_start_index + num_rel_eff_par) < x.size() );
    assert( (rel_eff_index != 0) || (rel_eff_start_index == m_rel_eff_par_start_index) );
    
    const T * const rel_ef_pars = &(x[rel_eff_start_index]);
    
    if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      assert( (rel_eff_start_index + num_rel_eff_par) < x.size() );
      assert( m_drf );
      
      const PhysModelRelEqnDef<T> re_input = make_phys_eqn_input( rel_eff_index, x );
      
      // TODO: this next call can be a real bottlenck - specifically the `GammaInteractionCalc::transmition_length_coefficient` function can be really slow - we should some-how maybe cache these results
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

    // If we are doing numeric differentiation, and the calibration paramaters are all
    //  default/nominal values, then lets skip the math.
    // However, if we are doing auto-differentiation, we need to do the math so we can
    //  carry the gradient information into the answer.
    if constexpr ( std::is_same_v<T, double> )
    {
      if( (x[0] == RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset)
         && (x[1] == RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset)
         && ((RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 2)
             || (x[2] == RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset)) )
      {
        return T(energy);
      }
    }//if constexpr ( std::is_same_v<T, double> )

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
    
    assert( parameter_scale_factor(0) == (RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) );
    assert( parameter_scale_factor(1) == (RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) );
    assert( (RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars <= 2)
         || (parameter_scale_factor(2) == (RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset)) );

    
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
      assert( fabs(x[0] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) < 1.0E-5 );
      assert( fabs(x[1] - RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset) < 1.0E-5 );
      
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

  
  
  /** Computes peaks for a ROI range, given current paramaters
   @param range The channel range to generate peaks for
   @param x The Ceres paramaters to use to form the peaks
   @param multithread Wether to use a single, or multiple threads to comput the peaks
   */
  template<typename T>
  RelActCalcAuto::PeaksForEnergyRangeImp<T> peaks_for_energy_range_imp( const RoiRangeChannels &range,
                                                       const std::vector<T> &x,
                                                       const bool multithread ) const
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


    // Throws lambda if any peak quantities are NaN or Inf
    auto check_peak_reasonable = [this,&x]( const RelActCalcAuto::PeakDefImp<T> &peak, const double gamma_energy ){
      const T &peak_sigma = peak.m_sigma;

      double fwhm;
      if constexpr ( !std::is_same_v<T, double> )
        fwhm = 2.35482*peak_sigma.a;
      else
        fwhm = 2.35482 * peak_sigma;

      if( isinf(peak_sigma) || isnan(peak_sigma) )
      {
        stringstream msg;
        msg << "peaks_for_energy_range_imp: " << fwhm << " FWHM for "
        << std::setprecision(2) << gamma_energy << " keV, from pars={";
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
      }//if( IsInf(peak_sigma) || IsNan(peak_sigma) )

      // Do a sanity check to make sure peak isnt getting too narrow
      const double nchannel = m_energy_cal->channel_for_energy(gamma_energy + 0.5*fwhm)
                                       - m_energy_cal->channel_for_energy(gamma_energy - 0.5*fwhm);
      if( nchannel < 1.5 )
        throw runtime_error( "peaks_for_energy_range_imp: for peak at " + std::to_string(gamma_energy)
                            + " keV, FWHM=" + std::to_string(fwhm) + " which is only "
                            + std::to_string(nchannel) + "channels - too small." );

      if( fwhm < 0.001 )
        throw runtime_error( "peaks_for_energy_range_imp: for peak at " + std::to_string(gamma_energy)
                            + " keV, FWHM=" + std::to_string(fwhm) + " which is too small." );


      const T &peak_amplitude = peak.m_amplitude;
      if( isinf(peak_amplitude) || isnan(peak_amplitude) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak amplitude for "
                            + std::to_string(gamma_energy) + " keV.");

      const T &peak_mean = peak.m_mean;
      if( isinf(peak_mean) || isnan(peak_mean) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak mean for "
                            + std::to_string(gamma_energy) + " keV.");

    };//auto check_peak_reasonable

    // We will accept peaks outdide the ROI range, if they will still effect the ROI, so here we will define
    //  a peak for the lower and upper edges of the ROI, and use thier extent to approximate peaks
    //  above and below the ROI (in principle we should search around, but its not too important).
    //
    //  TODO: We are including peaks such that 1E-4 of them make it into the ROI - this is totally arbitrary, and instead we should take into account the peak area.
    RelActCalcAuto::PeakDefImp<T> lower_range_peak, upper_range_peak;
    lower_range_peak.m_mean = T(range.lower_energy);
    lower_range_peak.m_sigma = fwhm( T(range.lower_energy), x ) / 2.35482;
    lower_range_peak.m_amplitude = T(1.0);
    set_peak_skew( lower_range_peak, x );

    upper_range_peak.m_mean = T(range.upper_energy);
    upper_range_peak.m_sigma = fwhm( T(range.upper_energy), x ) / 2.35482;
    upper_range_peak.m_amplitude = T(1.0);
    set_peak_skew( upper_range_peak, x );

    check_peak_reasonable( lower_range_peak, range.lower_energy );
    check_peak_reasonable( upper_range_peak, range.upper_energy );

    // The Cyrstal Ball distributions can have reallllly long tails, so we will use a smaller coverage function
    //  for these.
    const bool is_crystal_ball = ((m_options.skew_type == PeakDef::SkewType::CrystalBall)
                                  || (m_options.skew_type == PeakDef::SkewType::DoubleSidedCrystalBall));

    const double missing_frac = is_crystal_ball ? 1.0E-3 : 1.0E-4;
    const pair<double,double> lower_peak_limits = lower_range_peak.peak_coverage_limits( missing_frac, 20.0 );
    const pair<double,double> upper_peak_limits = upper_range_peak.peak_coverage_limits( missing_frac, 20.0 );

    assert( range.lower_energy >= lower_peak_limits.first );
    assert( range.upper_energy <= upper_peak_limits.second );

    // TODO: for Crystal Ball dists, they can have really far-reaching tails - perhaps we should limit the max extent of peaks to save CPU, or whatever
    const double lower_mean = lower_peak_limits.first;
    const double upper_mean = upper_peak_limits.second;


    size_t num_free_peak_pars = 0;
    set<const void *> nuclides_used;
    
    RelActCalcAuto::PeaksForEnergyRangeImp<T> answer;
    answer.first_channel = first_channel;
    answer.last_channel = last_channel;
    answer.no_gammas_in_range = false;
    answer.forced_full_range = range.force_full_range;
    
    vector<RelActCalcAuto::PeakDefImp<T>> &peaks = answer.peaks;

    // Go through and create peaks based on rel act, eff, etc
    for( size_t rel_eff_index = 0; rel_eff_index < m_nuclides.size(); ++rel_eff_index )
    {
      const vector<NucInputGamma> &nuclides = m_nuclides[rel_eff_index];
      for( const NucInputGamma &src_info : nuclides )
      {
        const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(src_info.source);
        const SandiaDecay::Element * const el = RelActCalcAuto::element(src_info.source);
        const ReactionGamma::Reaction * const rctn = RelActCalcAuto::reaction(src_info.source);

        const bool fixed_age = ((!nuc) || is_fixed_age(nuc, rel_eff_index));
        double nuc_age_val = -1.0, forward_diff_nuc_age_val = -1.0, backward_diff_nuc_age_val = -1.0;
        std::shared_ptr<const vector<NucInputGamma::EnergyYield>> gammas, forward_diff_gammas, backward_diff_gammas;
        
        const T rel_act = relative_activity( src_info.source, rel_eff_index, x );
        //cout << "peaks_for_energy_range_imp: Relative activity of " << src_info.name()
        //     << " is " << PhysicalUnits::printToBestActivityUnits(rel_act) << endl;
        
        T nuc_age;
        if( fixed_age )
        {
          gammas = src_info.nominal_gammas;
          assert( gammas );
        }else
        {
          // If we are here, we are fitting the age - however, automatic differentiation does
          //  not work with age, so if our template type is a `Jet<>` we will do numeric 
          //  differentiation, and insert this info into the Jet.

          assert( nuc );

          nuc_age = age(src_info, rel_eff_index, x);
          
          if constexpr ( !std::is_same_v<T, double> )
          {
            nuc_age_val = nuc_age.a;
            gammas = this->decay_gammas( src_info, nuc_age_val, src_info.gammas_to_exclude );
            assert( gammas );
            
            // TODO: be a bit smarter about selecting dt
            // At larger ages, we dont expect much of a change in gamma BRs, so we'll be a
            //  bit larger of a time delta, to hopefully avoid getting zero derivative; this
            //  was not tested, and I'm not even sure its needed!
            double forward_dt = 0.001*nuc_age_val, backward_dt = -1.0;
            if( nuc_age_val > 0.1*nuc->halfLife )
            {
              if( nuc_age_val > 5.0*nuc->halfLife )
              {
                forward_dt = 0.1*nuc_age_val;
                backward_dt = 0.1*nuc_age_val;
              }else if( nuc_age_val > 2.5*nuc->halfLife )
              {
                forward_dt = 0.01*nuc_age_val;
                backward_dt = 0.01*nuc_age_val;
              }else
              {
                forward_dt = 0.001*nuc_age_val;
                backward_dt = 0.001*nuc_age_val;
              }
            }else
            {
              forward_dt = std::min( 0.01*nuc_age_val, 0.001*nuc->halfLife );
              if( nuc_age_val > 0.0 )
                backward_dt = std::min(0.1*nuc_age_val, 0.001*nuc->halfLife);
            }//if( nuc_age_val > 0.001*nuc->halfLife ) / else
            
            forward_diff_nuc_age_val = nuc_age_val + forward_dt;
            forward_diff_gammas = this->decay_gammas( src_info, forward_diff_nuc_age_val, src_info.gammas_to_exclude );

            assert( forward_diff_gammas );
            assert( forward_diff_gammas->size() == gammas->size() );

            if( backward_dt > 0.0 )
            {
              backward_diff_nuc_age_val = nuc_age_val - backward_dt;
              backward_diff_gammas = this->decay_gammas( src_info, backward_diff_nuc_age_val, src_info.gammas_to_exclude );
              assert( backward_diff_gammas );
              assert( backward_diff_gammas->size() == gammas->size() );
            }//if( nuc_age_val > 0.0 )
          }else
          {
            nuc_age_val = nuc_age;
            gammas = this->decay_gammas( src_info, nuc_age_val, src_info.gammas_to_exclude );
            assert( gammas );
          }
        }//if( age is fixed ) / else( age may vary )
        
        assert( gammas );
        
        for( size_t gamma_index = 0; gamma_index < gammas->size(); ++gamma_index )
        {
          // TODO: gammas right next to eachother may have exactly, or really close energies that we could combine into a single peak to save a little time - should consider doing this
          const NucInputGamma::EnergyYield &gamma = (*gammas)[gamma_index];
          const double energy = gamma.energy;
          T yield = T(gamma.yield);
          const size_t transition_index = gamma.transition_index;
          const SandiaDecay::Transition * const transition = gamma.transition;
          const PeakDef::SourceGammaType gamma_type = gamma.gamma_type;
          
          assert( !transition || transition->products.empty() || (transition_index < transition->products.size()) );
          assert( !nuc || (transition || (gamma_type == PeakDef::SourceGammaType::AnnihilationGamma)) );
          assert( !nuc
                 || !transition
                 || (gamma_type == PeakDef::SourceGammaType::AnnihilationGamma)
                 || (fabs(transition->products[transition_index].energy - energy) < 0.0001) );
          
          // Filter out zero-amplitude gammas
          // TODO: - come up with more intelligent lower bound of gamma rate to bother with
          const double min_br_wanted = std::numeric_limits<float>::min();
          if constexpr ( std::is_same_v<T, double> )
          {
            if( yield < min_br_wanted )
              continue;
          }else
          {
            if( fixed_age )
            {
              if( yield < min_br_wanted )
                continue;
            }else if( yield < min_br_wanted )
            {
              // only skip this gamma if nominal, and forward, and backward points are zero.
              assert( forward_diff_gammas );
              assert( forward_diff_gammas->size() == gammas->size() );
              assert( !backward_diff_gammas || (backward_diff_gammas->size() == gammas->size()) );
              
              const bool forward_is_zero = ((*forward_diff_gammas)[gamma_index].yield  < min_br_wanted);
              const bool backward_is_zero = (!backward_diff_gammas || ((*backward_diff_gammas)[gamma_index].yield  < min_br_wanted));
              if( forward_is_zero && backward_is_zero )
                continue;
            }//if( fixed age ) / else
          }//if( double ) / else ( ceres::Jet )
          
          // Filter the energy range based on true energies (i.e., not adjusted energies)
          if( (energy < lower_mean) || (energy > upper_mean) )
            continue;
          
          if( nuc )
            nuclides_used.insert( nuc );
          else if( el )
            nuclides_used.insert( el );
          else if( rctn )
            nuclides_used.insert( rctn );
          else
          {
            assert( 0 );
          }
          
          // We compute the relative efficiency and FWHM based off of "true" energy
          const T rel_eff = relative_eff( gamma.energy, rel_eff_index, x );
          if( isinf(rel_eff) || isnan(rel_eff) )
            throw runtime_error( "peaks_for_energy_range_imp: inf or NaN rel. eff for "
                                + std::to_string(gamma.energy) + " keV."  );

          
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
                const T par_value = x[m_add_br_uncert_start_index + range_index];
                const T num_sigma_from_nominal = (par_value - sm_peak_range_uncert_offset)/sm_peak_range_uncert_par_scale;

                br_uncert_adj = max(T(0.0), 1.0 + num_sigma_from_nominal*m_options.additional_br_uncert);
                
                break;
              }//
            }//for( find range this gamma belongs to, if any )
          }//if( m_options.additional_br_uncert > 0.0 && !m_peak_ranges_with_uncert.empty() )

          if constexpr ( !std::is_same_v<T, double> )
          {
            if( nuc && !is_fixed_age(nuc, rel_eff_index) )
            {
              // Here we will do the numeric differentiation and put the results in the Jet
              assert( forward_diff_gammas );
              assert( gamma_index < forward_diff_gammas->size() );
              const NucInputGamma::EnergyYield &gamma_forward = (*forward_diff_gammas)[gamma_index];
              assert( gamma_forward.energy == gamma.energy );
              // The `transition_index` can actually change with age.
              //  `decay_gammas(...)` notes "We will assign the transition to the largest yield", so since the yield
              //  changes with age, the transition may change.  It does really matter to us, since the yield was
              //  summed in `decay_gammas(..)` anyway.
              //assert( gamma_forward.transition_index == gamma.transition_index );
              //assert( gamma_forward.transition == gamma.transition );
              assert( gamma_forward.gamma_type == gamma.gamma_type );
              const double forward_yield = gamma_forward.yield;
              const double forward_derivative = (forward_yield - yield.a) / (forward_diff_nuc_age_val - nuc_age_val);

              if( backward_diff_gammas )
              {
                assert( gamma_index < backward_diff_gammas->size() );
                const NucInputGamma::EnergyYield &gamma_backward = (*backward_diff_gammas)[gamma_index];
                assert( gamma_backward.energy == gamma.energy );
                //See note above about why this next assert is commented out
                //assert( gamma_backward.transition_index == gamma.transition_index );
                //assert( gamma_backward.transition == gamma.transition );
                assert( gamma_backward.gamma_type == gamma.gamma_type );
                
                const double f_x = yield.a;
                const double f_plus = gamma_forward.yield;
                const double f_minus = gamma_backward.yield;
                const double h_fwd = forward_diff_nuc_age_val - nuc_age_val;
                const double h_bwd = nuc_age_val - backward_diff_nuc_age_val;
                const double hf2 = h_fwd * h_fwd;
                const double hb2 = h_bwd * h_bwd;
                const double hf_x_hb = h_fwd * h_bwd;
                const double hf_p_hb = h_fwd + h_bwd;

                const double numerator = (hb2 * f_plus) - (hf2 * f_minus) + (hf2 - hb2) * f_x;
                const double denominator = hf_x_hb * hf_p_hb;
                const double derivative = numerator / denominator;

                yield.v = derivative * nuc_age.v;
              }
            }//if( is_fixed_age(src_info.nuclide, rel_eff_index) )
          }//if( !std::is_same_v<T, double> )          

          const T peak_mean = apply_energy_cal_adjustment( gamma.energy, x );
          const T peak_amplitude = rel_act * static_cast<double>(m_live_time) * rel_eff * yield * br_uncert_adj;
          const T peak_fwhm = fwhm( T(gamma.energy), x );


          RelActCalcAuto::PeakDefImp<T> peak;
          peak.m_mean = peak_mean;
          peak.m_sigma = peak_fwhm/2.35482;
          peak.m_amplitude = peak_amplitude;
          peak.m_src_energy = gamma.energy;
          
          if( transition || (gamma_type == PeakDef::SourceGammaType::AnnihilationGamma) )
          {
            peak.m_parent_nuclide = nuc;
            peak.m_transition = transition;
            peak.m_rad_particle_index = transition_index;
            peak.m_gamma_type = gamma_type;
          }else if( el )
          {
            peak.m_xray_element = el;
            assert( gamma_type == PeakDef::SourceGammaType::XrayGamma );
            peak.m_gamma_type = PeakDef::SourceGammaType::XrayGamma;
          }else if( rctn )
          {
            peak.m_reaction = rctn;
            peak.m_gamma_type = gamma_type;
          }
          
          set_peak_skew( peak, x );
          peak.m_rel_eff_index = rel_eff_index;

          check_peak_reasonable( peak, gamma.energy );


          if( peak_amplitude < static_cast<double>( std::numeric_limits<float>::min() ) )
            continue;

          peaks.push_back( std::move(peak) );
        }//for( const SandiaDecay::EnergyRatePair &gamma : gammas )
      }//for( const NucInputGamma &src_info : m_nuclides )
    }//for( size_t rel_eff_index = 0; rel_eff_index < m_nuclides.size(); ++rel_eff_index )
    
    const size_t free_peaks_start_index = m_free_peak_par_start_index;
    assert( m_free_peak_area_multiples.size() == m_options.floating_peaks.size() );

    for( size_t index = 0; index < m_options.floating_peaks.size(); ++index )
    {
      const RelActCalcAuto::FloatingPeak &peak = m_options.floating_peaks[index];
      if( (peak.energy < lower_mean) || (peak.energy > upper_mean) )
        continue;
      
      const size_t amp_index = free_peaks_start_index + 2*index + 0;
      const size_t fwhm_index = free_peaks_start_index + 2*index + 1;
      assert( fwhm_index < x.size() );

      const double area_multiple = m_free_peak_area_multiples[index];
      
      num_free_peak_pars += 1;
      // TODO: we should use some sort of a scale-factor for the peak amplitude
      const T peak_amp = x[amp_index] * area_multiple;
      if( isinf(peak_amp) || isnan(peak_amp) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak amplitude for "
                            + std::to_string(peak.energy) + " keV extra peak.");
      
      T peak_mean( peak.energy );
      if( peak.apply_energy_cal_correction )
        peak_mean = apply_energy_cal_adjustment( peak.energy, x );
      
      if( isinf(peak_mean) || isnan(peak_mean) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak mean for "
                            + std::to_string(peak.energy) + " keV extra peak.");
      
      const T peak_fwhm = (peak.release_fwhm ? x[fwhm_index] : T(1.0)) * fwhm(peak_mean, x);
      
      if( isinf(peak_fwhm) || isnan(peak_fwhm) )
        throw runtime_error( "peaks_for_energy_range_imp: inf or NaN peak FWHM for "
                            + std::to_string(peak.energy) + " keV extra peak.");
      
      assert( peak.release_fwhm || (abs(x[fwhm_index] + 1.0) < 1.0E-5) );
      
      num_free_peak_pars += peak.release_fwhm;
      
      //cout << "peaks_for_energy_range_imp: free peak at " << peak.energy << " has a FWHM=" << peak_fwhm << " and AMP=" << peak_amp << endl;
      RelActCalcAuto::PeakDefImp<T> imp_peak;
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
      //cerr << "peaks_for_energy_range_imp: no peaks in range [" << range.lower_energy << ", "
      //     << range.upper_energy << "] keV." << endl;
      
      answer.no_gammas_in_range = true;
      
      const T middle_energy = 0.5*(adjusted_lower_energy + adjusted_upper_energy);
      T middle_fwhm = fwhm( middle_energy, x );
      if( isinf(middle_fwhm) || isnan(middle_fwhm) )
        middle_fwhm = T(1.0); //arbitrary
      
      
      RelActCalcAuto::PeakDefImp<T> peak;
      peak.m_mean = T(middle_energy);
      peak.m_sigma = T(middle_fwhm/2.35482);
      peak.m_amplitude = T(0.0);
      peak.m_rel_eff_index = std::numeric_limits<size_t>::max(); //should already be this value
      
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
    assert( m_channel_count_uncerts.size() == m_channel_counts.size() );
    assert( m_energy_cal && m_energy_cal->channel_energies() );
    assert( m_energy_cal->channel_energies()->size() >= m_channel_counts.size() );
    
    assert( static_cast<size_t>(num_polynomial_terms) < answer.continuum.m_values.size() );
    
    const T ref_energy = answer.continuum.m_reference_energy;
    const float * const data = &(m_channel_counts[first_channel]);
    const float * const data_uncerts = &(m_channel_count_uncerts[first_channel]);
    const float * const energies = &((*m_energy_cal->channel_energies())[first_channel]);

    std::vector<T> &peak_counts = answer.peak_counts;
    peak_counts.resize( num_channels, T(0.0) );
    T *continuum_coeffs = answer.continuum.m_values.data();

    vector<RelActCalcAuto::PeakDefImp<T>> dummy;

    PeakFit::fit_continuum( energies, data, data_uncerts, num_channels, num_polynomial_terms,
                                  is_step_continuum, ref_energy, peaks, multithread,
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
   @param rel_eff_indices GIves which Rel. Eff. curves to return peaks for.
          If empty, will return for all curves plus free peaks.
          If has all indices, then will not return free-peaks.
          Note: this argument being non-empty is actually less efficient than it being empty, because all peaks need to be
          calculated in order to fit the continuum anyway (so the filtering happens after this).
   
   @returns A peak for each gamma, of each nuclide, as well as each free-floating peak, in the
            problem, that is within the specified energy range, and in the specified Rel. Eff. curve(s).
   
    Currently this function is only used to get the final peaks fit, so not hyper optimized.
   */
  PeaksForEnergyRange peaks_for_energy_range( const RoiRangeChannels &range,
                                             const std::vector<double> &x,
                                             const std::set<size_t> &rel_eff_indices ) const
  {
    RelActCalcAuto::PeaksForEnergyRangeImp<double> computed_peaks = peaks_for_energy_range_imp( range, x, true );

    
    if( !rel_eff_indices.empty() )
    {
      vector<RelActCalcAuto::PeakDefImp<double>> filtered_peaks;
      for( const RelActCalcAuto::PeakDefImp<double> &p : computed_peaks.peaks )
      {
        if( rel_eff_indices.count(p.m_rel_eff_index) )
          filtered_peaks.push_back( p );
      }
      computed_peaks.peaks.swap( filtered_peaks );
    }//if( !rel_eff_indices.empty() )
    
    PeaksForEnergyRange answer;
    answer.first_channel = computed_peaks.first_channel;
    answer.last_channel = computed_peaks.last_channel;
    answer.no_gammas_in_range = computed_peaks.no_gammas_in_range;
    answer.forced_full_range = computed_peaks.forced_full_range;
    
    for( size_t i = 0; i < computed_peaks.peaks.size(); ++i )
    {
      const RelActCalcAuto::PeakDefImp<double> &comp_peak = computed_peaks.peaks[i];
      PeakDef peak( comp_peak.m_mean, comp_peak.m_sigma, comp_peak.m_amplitude );
      peak.setSkewType( comp_peak.m_skew_type );
      const size_t num_skew = PeakDef::num_skew_parameters( comp_peak.m_skew_type );
      for( size_t skew_index = 0; skew_index < num_skew; ++skew_index )
      {
        const PeakDef::CoefficientType coef
             = static_cast<PeakDef::CoefficientType>( static_cast<int>(PeakDef::CoefficientType::SkewPar0) + skew_index );
        peak.set_coefficient( comp_peak.m_skew_pars[skew_index], coef );
      }//for( size_t skew_index = 0; skew_index < num_skew; ++skew_index )
      
      if( comp_peak.m_parent_nuclide )
      {
        peak.setNuclearTransition( comp_peak.m_parent_nuclide, comp_peak.m_transition,
                                  static_cast<int>(comp_peak.m_rad_particle_index), comp_peak.m_gamma_type );
      }else if( comp_peak.m_xray_element )
      {
        peak.setXray( comp_peak.m_xray_element, comp_peak.m_src_energy );
      }else if( comp_peak.m_reaction )
      {
        peak.setReaction( comp_peak.m_reaction, comp_peak.m_src_energy, comp_peak.m_gamma_type );
      }
      
      if( i != 0 )
      {
        peak.setContinuum( answer.peaks.front().continuum() );
      }else
      {
        shared_ptr<PeakContinuum> cont = peak.continuum();
        assert( cont );
        const RelActCalcAuto::PeakContinuumImp<double> &comp_cont = computed_peaks.continuum;
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
          const NucInputGamma &input_src = rel_eff_nucs[nuc_index];
          if( (RelActCalcAuto::nuclide(input_src.source) == comp_peak.m_parent_nuclide)
             && (RelActCalcAuto::element(input_src.source) == comp_peak.m_xray_element)
             && (RelActCalcAuto::reaction(input_src.source) == comp_peak.m_reaction) )
          {
            nuc_info = &(rel_eff_nucs[nuc_index]);
          }
        }//
        
        assert( nuc_info );
        if( nuc_info && !nuc_info->peak_color_css.empty() )
          peak.setLineColor( Wt::WColor( Wt::WString::fromUTF8(nuc_info->peak_color_css) ) );
      }//if( comp_peak.m_rel_eff_index < m_nuclides.size() )
      
      // Put in a way to label the peak with the relative efficiency index, so we can tell which
      //  relative efficiency curve the peak came from, if there are multiple.
      if( m_nuclides.size() > 1 )
        peak.setUserLabel( "RelEff " + std::to_string(comp_peak.m_rel_eff_index) );

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
    
    vector<RelActCalcAuto::PeaksForEnergyRangeImp<T>> peaks_in_ranges_imp( m_energy_ranges.size() );

    // TODO: multi-thread computation needs to be looked at more hollistically, both here and in `PeakFit::fit_continuum(...)`, and possibly in peaks_for_energy_range_imp
    const bool multhread_each_roi = (m_energy_ranges.size() < 6);
    
    // Calling `m_pool.join()` puts the threadpool in a state where it will no longer work,
    //  so instead we will manually manage waiting on jobs to finish.
    std::mutex cv_mutex;
    std::condition_variable cv;
    size_t tasks_completed = 0;
    std::string exception_msg;
    const size_t num_energy_ranges = m_energy_ranges.size();
    
    for( size_t i = 0; i < num_energy_ranges; ++i )
    {
      boost::asio::post( m_pool,
                        [i,&peaks_in_ranges_imp,this,&x,multhread_each_roi,&cv,&cv_mutex,&tasks_completed,&exception_msg](){
        try
        {
          peaks_in_ranges_imp[i] = peaks_for_energy_range_imp( m_energy_ranges[i], x, multhread_each_roi );
        }catch( std::exception &e )
        {
          std::lock_guard<std::mutex> lock(cv_mutex);
          exception_msg = e.what();
        }
        
        
        std::lock_guard<std::mutex> lock(cv_mutex); //lock_guard is simpler than unique_lock, so use it here
        tasks_completed += 1;
        cv.notify_one();
      } );
    }//for( size_t i = 0; i < num_energy_ranges; ++i )
    
    {//begin wait for things to finish
      std::unique_lock<std::mutex> lock(cv_mutex);
      cv.wait(lock, [num_energy_ranges,&tasks_completed]() -> bool {
          return tasks_completed == num_energy_ranges;
      } );
    }//end wait for things to finish
    
    if( !exception_msg.empty() )
      throw runtime_error( exception_msg );
    
    if( m_cancel_calc && m_cancel_calc->load() )
    {
      do_cancel();
      return;
    }//if( m_cancel_calc && m_cancel_calc->load() )
    
    assert( SpecUtilsAsync::num_logical_cpu_cores() > 0 );


    
    assert( peaks_in_ranges_imp.size() == m_energy_ranges.size() );
    
    size_t residual_index = 0;
    
    const bool try_better_par = true;
    
    for( size_t roi_index = 0; roi_index < m_energy_ranges.size(); ++roi_index )
    {
      const RoiRangeChannels &energy_range = m_energy_ranges[roi_index];
      const RelActCalcAuto::PeaksForEnergyRangeImp<T> &info = peaks_in_ranges_imp[roi_index];

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
        assert( !isnan(this_residual[index]) && !isinf(this_residual[index]) );
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
        assert( !isnan(residuals[residual_index]) && !isinf(residuals[residual_index]) );
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
        const T num_sigma_from_nominal = (x[par_index] - sm_peak_range_uncert_offset)/sm_peak_range_uncert_par_scale;
        residuals[residual_index] = num_sigma_from_nominal; // m_options.additional_br_uncert;
        assert( !isnan(residuals[residual_index]) && !isinf(residuals[residual_index]) );
        ++residual_index;
      }
    }//if( !m_peak_ranges_with_uncert.empty() )
    
    
    assert( residual_index == number_residuals() );
#ifndef NDEBUG
    for( size_t i = 0; i < residual_index; ++i )
    {
      assert( !isnan(residuals[i]) && !isinf(residuals[i]) );
    }
#endif
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


NucInputGamma::NucInputGamma( const RelActCalcAuto::NucInputInfo &info, const RelActAutoCostFcn * const cost_fcn )
   : RelActCalcAuto::NucInputInfo( info )
{
  if( RelActCalcAuto::is_null(source) )
    throw runtime_error( "NucInputGamma: null Nuclide, Element, and Reaction." );
  
  const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(source);
  if( nuc && (age < 0.0) )
    throw runtime_error( "NucInputGamma: age may not be negative for nuclide (" + info.name() + ": "
                           + PhysicalUnits::printToBestTimeUnits(age)  + ")" );
  if( !nuc && (age > 0.0) )
    throw runtime_error( "NucInputGamma: age must be negative for Element or Reaction." );
    
  if( !nuc && fit_age  )
    throw runtime_error( "NucInputGamma: fit_age must be false for Element or Reaction." );
    
  assert( cost_fcn );
  nominal_gammas = cost_fcn->decay_gammas( info, info.age, gammas_to_exclude );
}//NucInputGamma constructor

}//namespace RelActCalcAutoImp


namespace RelActCalcAuto
{
const SandiaDecay::Nuclide *nuclide( const SrcVariant &src )
{
  if( std::holds_alternative<const SandiaDecay::Nuclide *>(src) )
    return std::get<const SandiaDecay::Nuclide *>(src);
  return nullptr;
}

const SandiaDecay::Element *element( const SrcVariant &src )
{
  if( std::holds_alternative<const SandiaDecay::Element *>(src) )
    return std::get<const SandiaDecay::Element *>(src);
  return nullptr;
}

const ReactionGamma::Reaction *reaction( const SrcVariant &src )
{
  if( std::holds_alternative<const ReactionGamma::Reaction *>(src) )
    return std::get<const ReactionGamma::Reaction *>(src);
  return nullptr;
}

bool is_null( const SrcVariant &src )
{
  return (!nuclide(src) && !element(src) && !reaction(src));
}


std::string to_name( const SrcVariant &src )
{
  const SandiaDecay::Nuclide *nuc = nuclide(src);
  if( nuc )
    return nuc->symbol;
  const SandiaDecay::Element *el = element(src);
  if( el )
    return el->symbol;
  const ReactionGamma::Reaction *rxn = reaction(src);
  if( rxn )
    return rxn->name();
  throw runtime_error( "to_name(SrcVariant): unknown source type" );
  return "";
}


SrcVariant source_from_string( const std::string &name )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "source_from_string: no decay database" );
  
  const SandiaDecay::Nuclide *nuc = db->nuclide( name );
  if( nuc )
    return nuc;
  
  const SandiaDecay::Element *el = db->element( name );
  if( el )
    return el;
  
  const ReactionGamma * const reaction_db = ReactionGammaServer::database();
  assert( reaction_db );
  if( !reaction_db )
    throw runtime_error( "Couldnt load reaction DB" );
      
  vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
  reaction_db->gammas( name, possible_rctns );
      
  // TODO: we are currently taking the first reaction; however, in principle there could be multiple - however, `ReactionGamma` doesnt have an interface to just return a reaction by name, I guess because
  for( const ReactionGamma::ReactionPhotopeak &r : possible_rctns )
    return r.reaction;

  return std::monostate();
}//source_from_string()


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
            nucinfo.source = n;
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
            if( en_nuc.second.count( RelActCalcAuto::nuclide(n.source) ) )
              common_age = std::max( common_age, n.age );
          }
          
          for( auto &n : nuclides )
          {
            if( en_nuc.second.count( RelActCalcAuto::nuclide(n.source) ) )
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
    case FwhmForm::Gadras:        return "Gadras";
    case FwhmForm::SqrtEnergyPlusInverse:  return "SqrtEnergyPlusInverse";
    case FwhmForm::ConstantPlusSqrtEnergy: return "ConstantPlusSqrtEnergy";
    case FwhmForm::Polynomial_2:  return "Polynomial_2";
    case FwhmForm::Polynomial_3:  return "Polynomial_3";
    case FwhmForm::Polynomial_4:  return "Polynomial_4";
    case FwhmForm::Polynomial_5:  return "Polynomial_5";
    case FwhmForm::Polynomial_6:  return "Polynomial_6";
    case FwhmForm::NotApplicable: return "NotApplicable";
  }//switch( form )
  
  assert( 0 );
  return "";
}//to_str( const FwhmForm form )


FwhmForm fwhm_form_from_str( const char *str )
{
  const size_t str_len = strlen(str);
  
  for( int iform = 0; iform <= static_cast<int>(FwhmForm::NotApplicable); iform += 1 )
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


const char *to_str( const FwhmEstimationMethod form )
{
  switch( form )
  {
    case FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum: 
      return "StartFromDetEffOrPeaksInSpectrum";
    case FwhmEstimationMethod::StartingFromAllPeaksInSpectrum: 
      return "StartingFromAllPeaksInSpectrum";
    case FwhmEstimationMethod::FixedToAllPeaksInSpectrum: 
      return "FixedToAllPeaksInSpectrum";
    case FwhmEstimationMethod::StartingFromDetectorEfficiency: 
      return "StartingFromDetectorEfficiency";
    case FwhmEstimationMethod::FixedToDetectorEfficiency: 
      return "FixedToDetectorEfficiency";
  }
  
  assert( 0 );
  return "InvalidFwhmEstimationMethod";
}//to_str( const FwhmEstimationMethod form )


FwhmEstimationMethod fwhm_estimation_method_from_str( const char *str )
{
  const size_t str_len = strlen(str);
  
  for( int iform = 0; iform <= static_cast<int>(FwhmEstimationMethod::FixedToDetectorEfficiency); iform += 1 )
  {
    const FwhmEstimationMethod x = FwhmEstimationMethod(iform);
    const char *form_str = to_str( x );
    const bool case_sensitive = false;
    const size_t form_str_len = strlen(form_str);
    if( rapidxml::internal::compare(str, str_len, form_str, form_str_len, case_sensitive) )
      return x;
  }
  
  throw runtime_error( "fwhm_estimation_method_from_str(...): invalid input string '" + std::string(str) + "'" );
  return FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
}//FwhmEstimationMethod fwhm_estimation_method_from_str(str)

  
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
T eval_fwhm( const T energy, const FwhmForm form, const T * const pars, const size_t num_pars, const shared_ptr<const DetectorPeakResponse> &drf )
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

    case RelActCalcAuto::FwhmForm::NotApplicable:
      assert( drf && drf->isValid() && drf->hasResolutionInfo() );
      if( !drf || !drf->hasResolutionInfo() )
        throw runtime_error( "RelActCalcAuto::eval_fwhm(): invalid detector response function - how did we get here???" );
      
      if constexpr ( std::is_same_v<T, double> )
      {
        return drf->peakResolutionFWHM( static_cast<float>(energy) );
      }else
      {
        return T(drf->peakResolutionFWHM( static_cast<float>(energy.a) ));
      }
      break;
  }//switch( m_options.fwhm_form )

  assert( fctntype != DetectorPeakResponse::kNumResolutionFnctForm );
  if( fctntype == DetectorPeakResponse::kNumResolutionFnctForm )
    throw runtime_error( "eval_fwhm: invalid FwhmForm" );
  
  //return
  const T answer = peakResolutionFWHM( energy, fctntype, pars, num_pars );

  if( isnan(answer) || isinf(answer) )
    throw runtime_error( "eval_fwhm: inf/NaN result." );
  
  if constexpr ( !std::is_same_v<T, double> )
  {
    if( answer.a <= 0.0 )
      throw runtime_error( "eval_fwhm: negative result." );
  }else
  {
    if( answer <= 0.0 )
      throw runtime_error( "eval_fwhm: negative result." );
  }
  
#ifndef NDEBUG
  float energy_kev;
  vector<float> drf_pars;

  if constexpr ( !std::is_same_v<T, double> )
  {
    energy_kev = static_cast<float>( energy.a );
    for( size_t i = 0; i < num_pars; ++i )
      drf_pars.push_back( static_cast<float>(pars[i].a) );
  }else
  {
    energy_kev = static_cast<float>( energy );
    for( size_t i = 0; i < num_pars; ++i )
      drf_pars.push_back( static_cast<float>(pars[i]) );
  }

  const double drf_answer = DetectorPeakResponse::peakResolutionFWHM( energy_kev, fctntype, drf_pars );

  assert( (abs(answer - drf_answer) < 0.01) || (abs(answer - drf_answer) < 0.01*max(answer,T(drf_answer))) );
#endif

  return answer;
}//T eval_fwhm( const T energy, const FwhmForm form, ... )


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


const std::string NucInputInfo::name() const
{
  if( RelActCalcAuto::is_null(source) )
    return "";

  return RelActCalcAuto::to_name(source);
}//NucInputInfo::name()

bool NucInputInfo::operator==( const NucInputInfo &rhs ) const
{
  return ((this->source == rhs.source) &&
         (this->age == rhs.age) &&
         (this->fit_age == rhs.fit_age) &&
         (this->fit_age_min == rhs.fit_age_min) &&
         (this->fit_age_max == rhs.fit_age_max) &&
         (this->min_rel_act == rhs.min_rel_act) &&
         (this->max_rel_act == rhs.max_rel_act) &&
         (this->starting_rel_act == rhs.starting_rel_act) &&
         (this->gammas_to_exclude == rhs.gammas_to_exclude) &&
         (this->peak_color_css == rhs.peak_color_css));
}

void NucInputInfo::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( !RelActCalcAuto::is_null(source) );
  if( RelActCalcAuto::is_null(source) )
    throw runtime_error( "NucInputInfo::toXml: null nuclide, element, or reaction." );
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "NucInputInfo::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "NucInputInfo" );
  parent->append_node( base_node );
  
  append_version_attrib( base_node, NucInputInfo::sm_xmlSerializationVersion );
  
  const string src_name = RelActCalcAuto::to_name(source);

  if( RelActCalcAuto::nuclide(source) )
  {
    append_string_node( base_node, "Nuclide", src_name );  
    const string age_str = PhysicalUnits::printToBestTimeUnits(age, 10);
    append_string_node( base_node, "Age", age_str);
    append_bool_node( base_node, "FitAge", fit_age );
    if( fit_age_min.has_value() )
      append_string_node( base_node, "FitAgeMin", PhysicalUnits::printToBestTimeUnits(fit_age_min.value()) );
    if( fit_age_max.has_value() )
      append_string_node( base_node, "FitAgeMax", PhysicalUnits::printToBestTimeUnits(fit_age_max.value()) );
  }else if( RelActCalcAuto::element(source) )
  {
    append_string_node( base_node, "Element", src_name );
  }else if( RelActCalcAuto::reaction(source) )
  {
    append_string_node( base_node, "Reaction", src_name );
  }
  
  if( min_rel_act.has_value() )
    append_float_node( base_node, "MinRelAct", min_rel_act.value() );
  if( max_rel_act.has_value() )
    append_float_node( base_node, "MaxRelAct", max_rel_act.value() );
  
  if( starting_rel_act.has_value() )
    append_float_node( base_node, "StartingRelAct", starting_rel_act.value() );

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
    
    age = -1.0;
    source = std::monostate();
    fit_age = false;
    fit_age_min = std::nullopt;
    fit_age_max = std::nullopt;
    min_rel_act = std::nullopt;
    max_rel_act = std::nullopt;
    starting_rel_act = std::nullopt;
    gammas_to_exclude.clear();


    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

    const rapidxml::xml_node<char> *nuc_node = XML_FIRST_NODE( nuc_info_node, "Nuclide" );
    const rapidxml::xml_node<char> *element_node = XML_FIRST_NODE( nuc_info_node, "Element" );
    const rapidxml::xml_node<char> *reaction_node = XML_FIRST_NODE( nuc_info_node, "Reaction" );
    
    if( !nuc_node && !element_node && !reaction_node )
      throw runtime_error( "Missing 'Nuclide', 'Element', or 'Reaction' node." );
    
    
    if( nuc_node )
    {
      const string nuc_str = SpecUtils::xml_value_str( nuc_node );
      const SandiaDecay::Nuclide * const nuc = db->nuclide( nuc_str );
      if( !nuc )
        throw runtime_error( "Invalid nuclide '" + nuc_str + "'" );
      source = nuc;
      const auto age_node = XML_FIRST_INODE( nuc_info_node, "Age");
      const string age_str = SpecUtils::xml_value_str( age_node );
      if( age_str.empty() || SpecUtils::iequals_ascii(age_str, "default") )
        age = PeakDef::defaultDecayTime(nuc);
      else
        age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife );
    
      if( XML_FIRST_NODE(nuc_info_node, "FitAge") )
        fit_age = get_bool_node_value( nuc_info_node, "FitAge" );

      const rapidxml::xml_node<char> *fit_age_min_node = XML_FIRST_NODE(nuc_info_node, "FitAgeMin");
      if( fit_age_min_node )
      {
        const string strval = SpecUtils::xml_value_str( fit_age_min_node );
        fit_age_min = PhysicalUnits::stringToTimeDurationPossibleHalfLife( strval, nuc->halfLife );
      }

      const rapidxml::xml_node<char> *fit_age_max_node = XML_FIRST_NODE(nuc_info_node, "FitAgeMax");
      if( fit_age_max_node )
      {
        const string strval = SpecUtils::xml_value_str( fit_age_max_node );
        fit_age_max = PhysicalUnits::stringToTimeDurationPossibleHalfLife( strval, nuc->halfLife );
      }
    }else if( element_node )
    {
      const string element_str = SpecUtils::xml_value_str( element_node );
      const SandiaDecay::Element * const element = db->element( element_str );
      if( !element )
        throw runtime_error( "Invalid element '" + element_str + "'" );
      source = element;
    }else if( reaction_node )
    {
      const string reaction_str = SpecUtils::xml_value_str( reaction_node );
      const ReactionGamma * const reaction_db = ReactionGammaServer::database();
      assert( reaction_db );
      if( !reaction_db )
        throw runtime_error( "Couldnt load reaction DB" );
      
      vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
      reaction_db->gammas( reaction_str, possible_rctns );
      
      const ReactionGamma::Reaction *reaction = nullptr;
      // TODO: we are currently taking the first reaction; however, in principle there could be multiple - however, `ReactionGamma` doesnt have an interface to just return a reaction by name, I guess because
      for( size_t i = 0; !reaction && (i < possible_rctns.size()); ++i )
        reaction = possible_rctns[i].reaction;
      
      if( !reaction )
        throw runtime_error( "Invalid reaction '" + reaction_str + "'" );

      source = reaction;
    }else
    {
      assert( 0 );
      throw runtime_error( "Missing 'Nuclide', 'Element', or 'Reaction' node." );
    }
        
    if( XML_FIRST_NODE(nuc_info_node, "MinRelAct") )
      min_rel_act = get_float_node_value( nuc_info_node, "MinRelAct" );
    if( XML_FIRST_NODE(nuc_info_node, "MaxRelAct") )
      max_rel_act = get_float_node_value( nuc_info_node, "MaxRelAct" );
    if( XML_FIRST_NODE(nuc_info_node, "StartingRelAct") )
      starting_rel_act = get_float_node_value( nuc_info_node, "StartingRelAct" );

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
    static_assert( FloatingPeak::sm_xmlSerializationVersion == 0,
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


void Options::check_same_hoerl_and_external_shielding_specifications() const
{
  if( !same_hoerl_for_all_rel_eff_curves && !same_external_shielding_for_all_rel_eff_curves )
    return;

  if( same_hoerl_for_all_rel_eff_curves )
  {
    size_t num_hoerl_curves = 0;
    for( const auto &rel_eff_curve : rel_eff_curves )
    {
      if( (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) 
        && rel_eff_curve.phys_model_use_hoerl )
        ++num_hoerl_curves;
    }
    if( num_hoerl_curves < 2 )
      throw logic_error( "If using the same Hoerl function for all relative efficiency curves, "
                        "you must have at least two relative efficiency curves with the Hoerl function enabled." );
  }//if( same_hoerl_for_all_rel_eff_curves )

  if( same_external_shielding_for_all_rel_eff_curves )
  {
    // Check at least 
    size_t num_external_shielding_curves = 0;
    const RelActCalcAuto::RelEffCurveInput *first_phys_model = nullptr;
    for( const auto &rel_eff_curve : rel_eff_curves )
    {
      if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel ) 
      {
        const vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>> &external_shieldings = rel_eff_curve.phys_model_external_atten;
        if( external_shieldings.empty() )
          throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                            "each relative efficiency curve must have an external shielding specification." );
        
        ++num_external_shielding_curves;
        if( !first_phys_model )
        {
          first_phys_model = &rel_eff_curve;
        }else
        {
          const vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>> &first_ext_shields = first_phys_model->phys_model_external_atten;

          if( external_shieldings.size() != first_ext_shields.size() )
            throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                              "all relative efficiency curves must have the same number of external shieldings specified." );
          for( size_t i = 0; i < external_shieldings.size(); ++i )
          {
            const RelActCalc::PhysicalModelShieldInput &ext_shield = *external_shieldings[i];
            const RelActCalc::PhysicalModelShieldInput &first_ext_shield = *first_ext_shields[i];

            if( (!ext_shield.material) != (!first_ext_shield.material) )
              throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                "all relative efficiency curves must consistently specify shielding materials." );

            if( ext_shield.material && first_ext_shield.material 
                && (ext_shield.material->name != first_ext_shield.material->name) )
              throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                "all relative efficiency curves must have the same external shielding materials ('" 
                                + ext_shield.material->name + "' vs '" + first_ext_shield.material->name + ")." );
            
            if( !ext_shield.material && (fabs(ext_shield.atomic_number - first_ext_shield.atomic_number) > 1.0e-6) )
              throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                "and they are not defined materials, the atomic numbers must be the same." );
            
            if( ext_shield.fit_atomic_number != first_ext_shield.fit_atomic_number )
              throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                "fitting for atomic number must be consistent." );
            
            if( !ext_shield.material && ext_shield.fit_atomic_number )
            {
              if( fabs(ext_shield.atomic_number - first_ext_shield.atomic_number) > 1.0e-6 )
                throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                  "and atomic number is being fit, the atomic numbers must be the same." );

              if( fabs(ext_shield.lower_fit_atomic_number - first_ext_shield.lower_fit_atomic_number) > 1.0e-6 )
                throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                  "and atomic number is being fit, the lower limits for the atomic numbers must be the same." );

              if( fabs(ext_shield.upper_fit_atomic_number - first_ext_shield.upper_fit_atomic_number) > 1.0e-6 )
                throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                  "and atomic number is being fit, the upper limits for the atomic numbers must be the same." );
            }//if( !ext_shield.material && ext_shield.fit_atomic_number )
            
            if( ext_shield.fit_areal_density != first_ext_shield.fit_areal_density )
              throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                "fitting for areal density must be consistent." );
            
            if( fabs(ext_shield.areal_density - first_ext_shield.areal_density) > 1.0e-6 )
                  throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                    "the areal densities must be the same." );

            if( ext_shield.fit_areal_density )
            {
              if( fabs(ext_shield.lower_fit_areal_density - first_ext_shield.lower_fit_areal_density) > 1.0e-6 )
                throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                  "and areal density is being fit, the lower limits for the areal densities must be the same." );

              if( fabs(ext_shield.upper_fit_areal_density - first_ext_shield.upper_fit_areal_density) > 1.0e-6 )
                throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                                  "and areal density is being fit, the upper limits for the areal densities must be the same." );
            }//if( ext_shield.fit_areal_density )
          }//for( size_t i = 0; i < external_shieldings.size(); ++i )
        }//if( !first_phys_model ) / else
      }//if( (rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) 
        
    }//for( const auto &rel_eff_curve : rel_eff_curves )


    if( num_external_shielding_curves < 2 )
      throw logic_error( "If using the same external shielding for all relative efficiency curves, "
                        "you must have at least two relative efficiency curves with the external shielding enabled." );
  }//if( same_external_shielding_for_all_rel_eff_curves )
}//void Options::check_same_hoerl_and_external_shielding_specifications() const


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
  
  const char *fwhm_estimation_method_str = to_str( fwhm_estimation_method );
  auto fwhm_estimation_method_node = append_string_node( base_node, "FwhmEstimationMethod", fwhm_estimation_method_str );
  append_attrib( fwhm_estimation_method_node, "remark", "Possible values: StartFromDetEffOrPeaksInSpectrum, StartingFromAllPeaksInSpectrum, FixedToAllPeaksInSpectrum, StartingFromDetectorEfficiency, FixedToDetectorEfficiency" );

  append_string_node( base_node, "Title", spectrum_title );
  
  append_string_node( base_node, "SkewType", PeakDef::to_string(skew_type) );
  
  append_float_node( base_node, "AddUncert", additional_br_uncert );
  
  
  xml_node<char> *rel_eff_node = doc->allocate_node( node_element, "RelEffCurveInputs" );
  base_node->append_node( rel_eff_node );
  for( const auto &curve : rel_eff_curves )
    curve.toXml( rel_eff_node );

  append_bool_node( base_node, "SameHoerlForAllRelEffCurves", same_hoerl_for_all_rel_eff_curves );
  append_bool_node( base_node, "SameExternalShieldingForAllRelEffCurves", same_external_shielding_for_all_rel_eff_curves );

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
    
    const rapidxml::xml_node<char> *fwhm_estimation_method_node = XML_FIRST_NODE( parent, "FwhmEstimationMethod" );
    if( fwhm_estimation_method_node )
    {
      const string method_str = SpecUtils::xml_value_str( fwhm_estimation_method_node );
      fwhm_estimation_method = fwhm_estimation_method_from_str( method_str.c_str() );
    }else
    {
      fwhm_estimation_method = FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
    }

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

      try
      {
        same_hoerl_for_all_rel_eff_curves = get_bool_node_value( parent, "SameHoerlForAllRelEffCurves" );
      }catch( ... )
      {
        same_hoerl_for_all_rel_eff_curves = false;
      }

      try
      {
        same_external_shielding_for_all_rel_eff_curves = get_bool_node_value( parent, "SameExternalShieldingForAllRelEffCurves" );
      }catch( ... )
      {
        same_external_shielding_for_all_rel_eff_curves = false;
      }
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
    case FwhmForm::Gadras:        return 3;
    case FwhmForm::SqrtEnergyPlusInverse:  return 3;
    case FwhmForm::ConstantPlusSqrtEnergy: return 2;
    case FwhmForm::Polynomial_2:  return 2;
    case FwhmForm::Polynomial_3:  return 3;
    case FwhmForm::Polynomial_4:  return 4;
    case FwhmForm::Polynomial_5:  return 5;
    case FwhmForm::Polynomial_6:  return 6;
    case FwhmForm::NotApplicable: return 0;
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


void RelEffCurveInput::check_nuclide_constraints() const
{
  // Make sure nuclides in constraints are non-null, not the same nuclide, and are in the nuclides 
  //  list that we are fitting for.
  for( const RelEffCurveInput::ActRatioConstraint &nuc_constraint : act_ratio_constraints )
  {
    if( nuc_constraint.constrained_to_controlled_activity_ratio <= 0.0 )
      throw logic_error( "Constrained to controlled activity ratio is less than or equal to 0.0." );

    if( !nuclide(nuc_constraint.constrained_source)
    && !element(nuc_constraint.constrained_source)
    && !reaction(nuc_constraint.constrained_source) )
      throw logic_error( "Constrained nuclide is nullptr." );

    if( !nuclide(nuc_constraint.controlling_source)
        && !element(nuc_constraint.controlling_source)
       && !reaction(nuc_constraint.controlling_source) )
      throw logic_error( "Controlling nuclide is nullptr." );

    if( nuc_constraint.constrained_source == nuc_constraint.controlling_source )
      throw logic_error( "Constrained and controlling nuclides are the same." );

    // Check that the constrained nuclide is a nuclide in this RelEffCurve
    bool is_constrained_nuclide_in_curve = false, is_controlling_nuclide_in_curve = false;
    for( const NucInputInfo &nuclide : nuclides )
    {
      if( nuc_constraint.constrained_source == nuclide.source )
        is_constrained_nuclide_in_curve = true;

      if( nuc_constraint.controlling_source == nuclide.source )
        is_controlling_nuclide_in_curve = true;
    }//for( const auto &nuclide : nuclides )
        
    if( !is_constrained_nuclide_in_curve )
      throw logic_error( "Constrained nuclide is not in the relative efficiency curve." );

    if( !is_controlling_nuclide_in_curve )
      throw logic_error( "Controlling nuclide is not in the relative efficiency curve." );

    for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &mass_fraction_constraint : mass_fraction_constraints )
    {
      const SandiaDecay::Nuclide * const mass_frac_nuc = mass_fraction_constraint.nuclide;
      if( !mass_frac_nuc )
        throw logic_error( "Mass fraction constraint nuclide is nullptr." );

      if( mass_frac_nuc == nuclide(nuc_constraint.constrained_source) )
        throw logic_error( "Constrained nuclide " + mass_frac_nuc->symbol + " is also a mass fraction constraint." );

      if( mass_frac_nuc == nuclide(nuc_constraint.controlling_source) )
      {
        // We'll let a mass fraction constrained nuclide constrol the activity for a nuclide of another element
        // but not of the same element.
        const SandiaDecay::Nuclide *constrained_nuclide = nuclide(nuc_constraint.constrained_source);
        if( constrained_nuclide && (mass_frac_nuc->atomicNumber == constrained_nuclide->atomicNumber) )
          throw logic_error( "Constrained nuclide " + constrained_nuclide->symbol 
                             + " is controlled by a mass fraction constrained nuclide of the same element (not currently supported)." );
      }//if( mass_fraction_constraint.nuclide == nuc_constraint.controlling_nuclide )
    }//for( const RelEffCurveInput::MassFractionConstraint &mass_fraction_constraint : mass_fraction_constraints )
  }//for( const RelEffCurveInput::ActRatioConstraint &nuc_constraint : act_ratio_constraints )

  // Check that the constrained nuclide is only controlled by one nuclide
  for( const ActRatioConstraint &nuc_constraint : act_ratio_constraints )
  {
    size_t count = 0;
    for( const ActRatioConstraint &other_constraint : act_ratio_constraints )
    {
      if( other_constraint.constrained_source == nuc_constraint.constrained_source )
        ++count;
    }

    if( count > 1 )
      throw logic_error( "Constrained nuclide is controlled by more than one nuclide." );
  }//for( const auto &nuc_constraint : act_ratio_constraints )

  // Make sure no duplicate constraints (although we should have caught this in the above 
  //  check for duplicate constrained sources)
  for( size_t i = 1; i < act_ratio_constraints.size(); ++i )
  {
    const RelEffCurveInput::ActRatioConstraint &outer_constraint = act_ratio_constraints[i];
    for( size_t j = 0; j < i; ++j )
    {
      const RelEffCurveInput::ActRatioConstraint &inner_constraint = act_ratio_constraints[j];
      if( (outer_constraint.constrained_source == inner_constraint.constrained_source)
        && (outer_constraint.controlling_source == inner_constraint.controlling_source) )
        {
          throw logic_error( "Duplicate nuclide constraints." );
        }
    }
  }//for( size_t i = 0; i < act_ratio_constraints.size(); ++i )

  // Now we need to walk the chain of constraints to make sure we dont have a cycle
  // e.g. {constrained: U235, controlling: U238} -> {constrained: U238, controlling: U234} -> {constrained: U234, controlling: U235}
  for( size_t outer_index = 0; outer_index < act_ratio_constraints.size(); ++outer_index )
  { 
    const RelEffCurveInput::ActRatioConstraint &outer_constraint = act_ratio_constraints[outer_index];
    
    set<size_t> visited_constraints;
    visited_constraints.insert( outer_index );

    bool found_controller = true;
    SrcVariant current_controller = outer_constraint.controlling_source;  // e.g. U238

    while( found_controller )
    {
      found_controller = false;
      
      for( size_t inner_index = 0; inner_index < act_ratio_constraints.size(); ++inner_index )
      {
        const RelEffCurveInput::ActRatioConstraint &inner_constraint = act_ratio_constraints[inner_index];
        if( current_controller == inner_constraint.constrained_source )
        {
          if( visited_constraints.count( inner_index ) )
            throw logic_error( "Cycle in nuclide constraints." );

          found_controller = true;
          current_controller = inner_constraint.controlling_source; // e.g. U234
          visited_constraints.insert( inner_index );
          break;
        }
      }//for( size_t inner_index = 0; inner_index < rel_eff_curve.act_ratio_constraints.size(); ++inner_index )
    }//while( found_constroller )
  }//for( const RelEffCurveInput::ActRatioConstraint &nuc_constraint : rel_eff_curve.act_ratio_constraints )


  //Now check the mass fraction constraints (we have already checked that the nuclide is not constrained by an act ratio constraint)
  std::map<short int, size_t> num_constrained_nucs_of_el;
  std::map<short int, double> lower_mass_fraction_sums_of_el;
  for( size_t index = 0; index < mass_fraction_constraints.size(); ++index )
  {
    const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint = mass_fraction_constraints[index];
    if( !constraint.nuclide )
      throw logic_error( "Mass fraction constraint nuclide is nullptr." );

    if( constraint.lower_mass_fraction < 0.0 || constraint.lower_mass_fraction > 1.0 )
      throw logic_error( "Mass fraction constraint (" + constraint.nuclide->symbol 
                            + ") has lower mass fraction (" + std::to_string(constraint.lower_mass_fraction) 
                            + ") which is not in [0.0, 1.0." );

    if( constraint.lower_mass_fraction > constraint.upper_mass_fraction )
      throw logic_error( "Mass fraction constraint (" + constraint.nuclide->symbol 
                            + ") has lower mass fraction (" + std::to_string(constraint.lower_mass_fraction) 
                            + ") which is greater than the upper mass fraction (" + std::to_string(constraint.upper_mass_fraction) + ")." );

    if( constraint.upper_mass_fraction < 0.0 || constraint.upper_mass_fraction > 1.0 )
      throw logic_error( "Mass fraction constraint (" + constraint.nuclide->symbol 
                            + ") has upper mass fraction (" + std::to_string(constraint.upper_mass_fraction) 
                            + ") which is not in (0.0, 1.0)." );

    if( (constraint.upper_mass_fraction == constraint.lower_mass_fraction)
        && ((constraint.lower_mass_fraction <= 0.0) || (constraint.lower_mass_fraction >= 1.0)) )
    {
      throw logic_error( "Mass fraction constraint (" + constraint.nuclide->symbol 
                            + ") has upper and lower mass fractions set to the same value (" 
                            + std::to_string(constraint.upper_mass_fraction) + "), which must be in (0, 1)." );
    }

    lower_mass_fraction_sums_of_el[constraint.nuclide->atomicNumber] += constraint.lower_mass_fraction;

    // Check that the constrained nuclide is a nuclide in this RelEffCurve
    size_t num_src_nucs_for_element = 0;
    bool is_constrained_nuclide_in_curve = false;
    
    num_constrained_nucs_of_el[constraint.nuclide->atomicNumber] += 1;

    for( const NucInputInfo &src : nuclides )
    {
      const SandiaDecay::Nuclide *src_nuc = nuclide(src.source);
      if( constraint.nuclide == src_nuc )
        is_constrained_nuclide_in_curve = true;

      if( src_nuc && (src_nuc->atomicNumber == constraint.nuclide->atomicNumber) )
        ++num_src_nucs_for_element;
    }//for( const NucInputInfo &nuclide : nuclides )

    if( !is_constrained_nuclide_in_curve )
      throw logic_error( "Mass fraction constraint nuclide (" + constraint.nuclide->symbol 
                            + ") is not in the relative efficiency curve." );

    if( num_src_nucs_for_element < 2  )
      throw logic_error( "Constrained nuclide (" + constraint.nuclide->symbol 
                            + ") must have at least 2 nuclides of the same element in the relative efficiency curve." );

    for( size_t i = index + 1; i < mass_fraction_constraints.size(); ++i )
    {
      const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &other_constraint = mass_fraction_constraints[i];
      if( other_constraint.nuclide == constraint.nuclide )
        throw logic_error( "Multiple mass fraction constraints for the same nuclide (" + constraint.nuclide->symbol + ")." );
    }//for( size_t i = index + 1; i < mass_fraction_constraints.size(); ++i )
  }//for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &mass_fraction_constraint : mass_fraction_constraints )

  // Check that the mass fraction sum of each element is 1.0
  for( const auto &[el, lower_mass_fraction_sum] : lower_mass_fraction_sums_of_el )
  {
    if( lower_mass_fraction_sum >= 1.0 )
      throw logic_error( "The lower mass fraction sum of element (AN=" + std::to_string(el) 
                          + ") is above 1.0 (sum=" + std::to_string(lower_mass_fraction_sum) + ")." );
  }//for( const [const short int &el, const double &lower_mass_fraction_sum] : lower_mass_fraction_sums_of_el )

  // Check that there is at least one nuclide of each element that is not constrained; to do
  // this we need to count the number of nuclides of each element total, and compare that to the
  // number of constrained nuclides of that element (num_constrained_nucs_of_el).
  std::map<short int, size_t> num_nucs_of_el;
  for( const NucInputInfo &src : nuclides )
  {
    const SandiaDecay::Nuclide *src_nuc = nuclide(src.source);
    if( src_nuc )
      num_nucs_of_el[src_nuc->atomicNumber] += 1;
  }
  
  for( const auto &[el, num_constrained_nucs] : num_constrained_nucs_of_el )
  {
    assert( num_nucs_of_el[el] > 0 );

    if( num_nucs_of_el[el] <= num_constrained_nucs )
      throw logic_error( "For each element with a mass fraction constraint, you must have at least one nuclide which is not constrained." );
  }//for( const [const short int &el, const size_t &num_nucs] : num_constrained_nucs_of_el )


  // Check that if rel act min/max/starting values are set, then nuclide is not constrained to another nuclide
  for( const NucInputInfo &nuc : nuclides )
  {
    assert( !is_null(nuc.source) );                                                                         
    if( is_null(nuc.source) )
      throw logic_error( "Nuclide, element, or reaction is nullptr." );
      
    if( nuc.min_rel_act.has_value() || nuc.max_rel_act.has_value() || nuc.starting_rel_act.has_value() )
    {
      for( const ActRatioConstraint &constraint : act_ratio_constraints )
      {
        if( constraint.constrained_source == nuc.source )
          throw runtime_error( "Nuclide " + nuc.name() + " is constrained to another nuclide"
                               " - you can not specify its min, max, or starting activity." );
      }//for( const ActRatioConstraint &constraint : act_ratio_constraints )

      for( const MassFractionConstraint &constraint : mass_fraction_constraints )
      {
        assert( constraint.nuclide );
        if( !constraint.nuclide )
          continue;

        if( constraint.nuclide == nuclide(nuc.source) )
          throw runtime_error( "Nuclide " + nuc.name() + " is a mass fraction constraint."
                               " - you can not specify its min, max, or starting activity." );
      }//for( const MassFractionConstraint &constraint : mass_fraction_constraints )
    }//if( nuc.min_rel_act.has_value() || nuc.max_rel_act.has_value() || nuc.starting_rel_act.has_value() )
    
    if( nuc.min_rel_act.has_value() && (nuc.min_rel_act.value() < 0.0) )
      throw runtime_error( "Nuclide " + nuc.name() + " has min rel act set to "
                        + std::to_string(nuc.min_rel_act.value()) + " - min rel act must be greater or equal to 0.0." );

    if( nuc.max_rel_act.has_value() && (nuc.max_rel_act.value() < 0.0) )
      throw runtime_error( "Nuclide " + nuc.name() + " has max rel act set to "
                      + std::to_string(nuc.max_rel_act.value()) + " - max rel act must be greater or equal to 0.0." );

    if( nuc.starting_rel_act.has_value() && (nuc.starting_rel_act.value() < 0.0) )
      throw runtime_error( "Nuclide " + nuc.name() + " has starting rel act set to "
                          + std::to_string(nuc.starting_rel_act.value())
                          + " - starting rel act must be greater or equal to 0.0." );

    //Check that if min_rel_act is set, and max_rel_act is set, then max_rel_act is greater than min_rel_act
    if( (nuc.min_rel_act.has_value() && nuc.max_rel_act.has_value())
        && (nuc.min_rel_act.value() > nuc.max_rel_act.value()) )
    {
      const string msg = "Nuclide " + nuc.name() + " has min rel act set to " + std::to_string(nuc.min_rel_act.value())
                         + " and max rel act set to " + std::to_string(nuc.max_rel_act.value())
                         + " - max rel act must be greater or equal to min rel act.";
      throw runtime_error( msg );
    }
    
    if( (nuc.starting_rel_act.has_value() && nuc.min_rel_act.has_value()) 
        && (nuc.starting_rel_act.value() < nuc.min_rel_act.value()) )
    {
      const string msg = "Nuclide " + nuc.name() + " has starting rel act set to "
                        + std::to_string(nuc.starting_rel_act.value())
                        + " - starting rel act must be greater or equal to min rel act.";
      throw runtime_error( msg );
    }

    if( (nuc.starting_rel_act.has_value() && nuc.max_rel_act.has_value()) 
        && (nuc.starting_rel_act.value() > nuc.max_rel_act.value()) )
    {
      const string msg = "Nuclide " + nuc.name() + " has starting rel act set to "
                          + std::to_string(nuc.starting_rel_act.value())
                          + " - starting rel act must be less or equal to max rel act.";
      throw runtime_error( msg );
    }

    const SandiaDecay::Nuclide *src_nuc = nuclide(nuc.source);
    if( !src_nuc )
    {
      if( nuc.age > 0.0 )
        throw runtime_error( "non-nuclide " + nuc.name() + " has age set to " + std::to_string(nuc.age)
                             + " - age must be less than or equal to 0.0." );
      if( nuc.fit_age )
        throw runtime_error( "non-nuclide " + nuc.name() + " has fit age set to true - fit age must be false." );
      if( nuc.fit_age_min.has_value() || nuc.fit_age_max.has_value() )
        throw runtime_error( "non-nuclide " + nuc.name() + " has fit age min or max set, but fit age is not set." );
    }else
    {
      if( !nuc.fit_age )
      {
        if( nuc.fit_age_min.has_value() || nuc.fit_age_max.has_value() )
          throw runtime_error( "Nuclide " + nuc.name() + " has fit age min or max set, but fit age is not set." );
      }

      if( nuc.age < 0.0 )
        throw runtime_error( "Nuclide " + nuc.name() + " has age set to " + std::to_string(nuc.age)
                            + " - age must be greater or equal to 0.0." );

      if( nuc.fit_age_min.has_value() && nuc.fit_age_max.has_value() )
      {
        if( nuc.fit_age_min.value() >= nuc.fit_age_max.value() )
          throw runtime_error( "Nuclide " + nuc.name() + " has fit age min set to "
                              + std::to_string(nuc.fit_age_min.value()) + " and fit age max set to "
                              + std::to_string(nuc.fit_age_max.value())
                              + " - fit age max must be greater or equal to fit age min." );
      }

      if( nuc.fit_age_min.has_value() && nuc.age < nuc.fit_age_min.value() )
        throw runtime_error( "Nuclide " + nuc.name() + " has age set to " + std::to_string(nuc.age)
                            + " - age must be greater or equal to fit age min." );

      if( nuc.fit_age_max.has_value() && nuc.age > nuc.fit_age_max.value() )
        throw runtime_error( "Nuclide " + nuc.name() + " has age set to " + std::to_string(nuc.age)
                            + " - age must be less or equal to fit age max." );
    }//if( !nuclide ) / else
  }//for( const NucInputInfo &nuc : nuclides )

  if( nucs_of_el_same_age )
  {
    for( const NucInputInfo &inner_nuc : nuclides )
    {
      const SandiaDecay::Nuclide *inner_src_nuc = nuclide(inner_nuc.source);
      if( !inner_src_nuc )
        continue;
      
      for( const NucInputInfo &outer_nuc : nuclides )
      {
        const SandiaDecay::Nuclide *outer_src_nuc = nuclide(outer_nuc.source);
        if( !outer_src_nuc )
          continue;
        
        if( inner_src_nuc->atomicNumber != outer_src_nuc->atomicNumber )
          continue;

        if( inner_nuc.age != outer_nuc.age )
          throw runtime_error( "When nucs_of_el_same_age is true, all nuclides of same element must have the same age input." );

        if( inner_nuc.fit_age != outer_nuc.fit_age )
          throw runtime_error( "When nucs_of_el_same_age is true, all nuclides of same element must have the same fit_age input." );

        if( (inner_nuc.fit_age_min.has_value() != outer_nuc.fit_age_min.has_value())
            || (inner_nuc.fit_age_max.has_value() != outer_nuc.fit_age_max.has_value())
            || (inner_nuc.fit_age_min.has_value() && (inner_nuc.fit_age_min.value() != outer_nuc.fit_age_min.value()))
            || (inner_nuc.fit_age_max.has_value() && (inner_nuc.fit_age_max.value() != outer_nuc.fit_age_max.value())) )
        {
          const string msg = "When nucs_of_el_same_age is true, all nuclides of same element must have the same"
          " fit_age_min and fit_age_max input."
          " - nuclide " + inner_nuc.name() + " has fit_age_min set to "
          + std::to_string(inner_nuc.fit_age_min.value_or(0.0)) + " and fit_age_max set to "
          + std::to_string(inner_nuc.fit_age_max.value_or(0.0))
          + ", but nuclide " + outer_nuc.name() + " has fit_age_min set to "
          + std::to_string(outer_nuc.fit_age_min.value_or(0.0))
          + " and fit_age_max set to " + std::to_string(outer_nuc.fit_age_max.value_or(0.0)) + ".";
          throw runtime_error( msg );
        }//if( same age limits not specified for all nuclides of same element )
      }//for( const NucInputInfo &outer_nuc : nuclides )
    }//for( const NucInputInfo &inner_nuc : nuclides )
  }//if( nucs_of_el_same_age )
}//void RelEffCurveInput::check_nuclide_constraints() const

    
RelEffCurveInput::ActRatioConstraint RelEffCurveInput::ActRatioConstraint::from_mass_ratio( const SandiaDecay::Nuclide *constrained, 
                                              const SandiaDecay::Nuclide *controlling, 
                                            const double mass_ratio_constrained_to_controlling )
{
  assert( constrained );
  assert( controlling );
  assert( mass_ratio_constrained_to_controlling > 0.0 );
  
  ActRatioConstraint answer;
  answer.constrained_source = constrained;
  answer.controlling_source = controlling;
  answer.constrained_to_controlled_activity_ratio = RelActCalc::mass_ratio_to_act_ratio(constrained, controlling, mass_ratio_constrained_to_controlling);

  return answer;
}


void RelEffCurveInput::ActRatioConstraint::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  if( is_null(controlling_source) || is_null(constrained_source) )
    throw logic_error( "Controlling or constrained nuclide is nullptr." );

  if( constrained_to_controlled_activity_ratio <= 0.0 )
    throw logic_error( "Constrained to controlled activity ratio is less than or equal to 0.0." );

  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "ActRatioConstraint::toXml: invalid parent." );

  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "ActRatioConstraint", nullptr, 18, 0 );
  parent->append_node( base_node );
  append_version_attrib( base_node, RelEffCurveInput::ActRatioConstraint::sm_xmlSerializationVersion );
  
  append_string_node( base_node, "ControllingSource", to_name(controlling_source) );
  append_string_node( base_node, "ConstrainedSource", to_name(constrained_source) );
  append_float_node( base_node, "ActivityRatio", constrained_to_controlled_activity_ratio );
}//void RelEffCurveInput::ActRatioConstraint::toXml( ::rapidxml::xml_node<char> *parent ) const


void RelEffCurveInput::ActRatioConstraint::fromXml( const ::rapidxml::xml_node<char> *constraint_node )
{
  using namespace rapidxml;
  
  if( !constraint_node )
    throw runtime_error( "invalid input" );
    
  if( !rapidxml::internal::compare( constraint_node->name(), constraint_node->name_size(), "ActRatioConstraint", 18, false ) )
    throw std::logic_error( "invalid input node name" );
    
  // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
  static_assert( RelEffCurveInput::ActRatioConstraint::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
  XmlUtils::check_xml_version( constraint_node, RelEffCurveInput::ActRatioConstraint::sm_xmlSerializationVersion );

  const rapidxml::xml_node<char> *controlling_node = XML_FIRST_INODE( constraint_node, "ControllingSource" );
  if( !controlling_node )
    controlling_node = XML_FIRST_INODE( constraint_node, "ControllingNuclide" );
  if( !controlling_node )
    throw std::runtime_error( "No <ControllingSource> or <ControllingNuclide> node" );

  const rapidxml::xml_node<char> *constrained_node = XML_FIRST_INODE( constraint_node, "ConstrainedSource" );
  if( !constrained_node )
    constrained_node = XML_FIRST_INODE( constraint_node, "ConstrainedNuclide" );
  if( !constrained_node )
    throw std::runtime_error( "No <ConstrainedSource> or <ConstrainedNuclide> node" );
  
  const string controlling_nuclide_name = SpecUtils::xml_value_str( controlling_node );
  const string constrained_nuclide_name = SpecUtils::xml_value_str( constrained_node );
  const double activity_ratio = XmlUtils::get_float_node_value( constraint_node, "ActivityRatio" );
  if( activity_ratio <= 0.0 )
    throw runtime_error( "Activity ratio is less than or equal to 0.0." );

  controlling_source = source_from_string( controlling_nuclide_name );
  if( is_null(controlling_source) )
    throw runtime_error( "Controlling nuclide ('" + controlling_nuclide_name + ") not found in decay database." );

  constrained_source = source_from_string( constrained_nuclide_name );
  if( is_null(constrained_source) )
    throw runtime_error( "Constrained nuclide ('" + constrained_nuclide_name + ") not found in decay database." );

  constrained_to_controlled_activity_ratio = activity_ratio;
}//void RelEffCurveInput::ActRatioConstraint::fromXml( const ::rapidxml::xml_node<char> *parent )  


#if( PERFORM_DEVELOPER_CHECKS )
void RelEffCurveInput::ActRatioConstraint::equalEnough( const ActRatioConstraint &lhs, const ActRatioConstraint &rhs )
{
  if( fabs(lhs.constrained_to_controlled_activity_ratio - rhs.constrained_to_controlled_activity_ratio) > 1e-6 )
    throw logic_error( "Constrained to controlled activity ratio is not equal." );

  if( lhs.controlling_source != rhs.controlling_source )
    throw logic_error( "Controlling nuclide is not equal." );

  if( lhs.constrained_source != rhs.constrained_source )
    throw logic_error( "Constrained nuclide is not equal." );
}//void RelEffCurveInput::ActRatioConstraint::equalEnough(...)
#endif


void RelEffCurveInput::MassFractionConstraint::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "RoiRange::toXml: invalid parent." );

  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "MassFractionConstraint", nullptr, 22, 0 );
  parent->append_node( base_node );
  append_version_attrib( base_node, MassFractionConstraint::sm_xmlSerializationVersion );

  append_string_node( base_node, "Nuclide", nuclide->symbol );
  if( lower_mass_fraction == upper_mass_fraction )
  {
    append_float_node( base_node, "MassFraction", lower_mass_fraction );  
  }else
  {
    append_float_node( base_node, "LowerMassFraction", lower_mass_fraction );
    append_float_node( base_node, "UpperMassFraction", upper_mass_fraction );
  }
}//void MassFractionConstraint::toXml( ::rapidxml::xml_node<char> *parent ) const


void RelEffCurveInput::MassFractionConstraint::fromXml( const ::rapidxml::xml_node<char> *parent )
{
  using namespace rapidxml;
  
  if( !parent )
    throw runtime_error( "invalid input" ); 
    
  if( !rapidxml::internal::compare( parent->name(), parent->name_size(), "MassFractionConstraint", 22, false ) )
    throw std::logic_error( "invalid input node name" );
    
  // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
  static_assert( MassFractionConstraint::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
                  
  const rapidxml::xml_node<char> *nuclide_node = XmlUtils::get_required_node( parent, "Nuclide" );
  const string nuclide_name = SpecUtils::xml_value_str( nuclide_node );

  double lower_mass_fraction_value = 0.0, upper_mass_fraction_value = 0.0;
  const rapidxml::xml_node<char> *mass_fraction_node = XML_FIRST_NODE( parent, "MassFraction" );
  
  if( mass_fraction_node )
  {
    if( !SpecUtils::parse_double(mass_fraction_node->value(), mass_fraction_node->value_size(), lower_mass_fraction_value) )
      throw std::runtime_error( "Invalid MassFraction XML value ('" + SpecUtils::xml_value_str(mass_fraction_node) + "')." );
    upper_mass_fraction_value = lower_mass_fraction_value;
  }else
  {
    lower_mass_fraction_value = XmlUtils::get_float_node_value( parent, "LowerMassFraction" );
    upper_mass_fraction_value = XmlUtils::get_float_node_value( parent, "UpperMassFraction" );
  }
  
  if( lower_mass_fraction_value < 0.0 || lower_mass_fraction_value >= 1.0 )
    throw runtime_error( "Lower mass fraction must be at least 0.0 and less than 1.0." );

  if( upper_mass_fraction_value <= 0.0 || upper_mass_fraction_value > 1.0 )
    throw runtime_error( "Upper mass fraction must be greater than 0.0 and less than or equal to 1.0." );

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "No decay database found." );  

  const SandiaDecay::Nuclide *nuclide_ptr = db->nuclide( nuclide_name );
  if( !nuclide_ptr )
    throw runtime_error( "Nuclide '" + nuclide_name + "' not found in decay database." );

  nuclide = nuclide_ptr;
  lower_mass_fraction = lower_mass_fraction_value;
  upper_mass_fraction = upper_mass_fraction_value;
}//void MassFractionConstraint::fromXml( const ::rapidxml::xml_node<char> *parent )


#if( PERFORM_DEVELOPER_CHECKS )
void RelEffCurveInput::MassFractionConstraint::equalEnough( const MassFractionConstraint &lhs, const MassFractionConstraint &rhs )
{
  if( lhs.nuclide != rhs.nuclide )
    throw logic_error( "Nuclide is not equal." );

  if( fabs(lhs.lower_mass_fraction - rhs.lower_mass_fraction) > 1e-7 )
    throw logic_error( "Lower mass fraction is not equal." );

  if( fabs(lhs.upper_mass_fraction - rhs.upper_mass_fraction) > 1e-7 )
    throw logic_error( "Upper mass fraction is not equal." );
}//void MassFractionConstraint::equalEnough(...)
#endif


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

  append_string_node( base_node, "Name", name );

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
  if( (rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel)
     || phys_model_self_atten
     || !phys_model_external_atten.empty() )
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

  if( !act_ratio_constraints.empty() )
  {
    xml_node<char> *act_ratio_constraints_node = doc->allocate_node( node_element, "ActRatioConstraints" );
    base_node->append_node( act_ratio_constraints_node );
    for( const auto &constraint : act_ratio_constraints )
      constraint.toXml( act_ratio_constraints_node );
  }//if( !act_ratio_constraints.empty() )

  if( !mass_fraction_constraints.empty() )
  {
    xml_node<char> *mass_fraction_constraints_node = doc->allocate_node( node_element, "MassFractionConstraints" );
    base_node->append_node( mass_fraction_constraints_node );
    for( const auto &constraint : mass_fraction_constraints )
      constraint.toXml( mass_fraction_constraints_node );
  }//if( !mass_fraction_constraints.empty() ) 

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
  
  const rapidxml::xml_node<char> *name_node = XML_FIRST_NODE(parent, "Name");
  name = SpecUtils::xml_value_str( name_node );

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

  // Get the ActRatioConstraints
  const rapidxml::xml_node<char> *act_ratio_constraints_node = XML_FIRST_NODE( parent, "ActRatioConstraints" );
  if( act_ratio_constraints_node )
  {
    XML_FOREACH_CHILD( constraint_node, act_ratio_constraints_node, "ActRatioConstraint" )
    {
      ActRatioConstraint constraint;
      constraint.fromXml( constraint_node );

      // Check that the nuclides are valid, and in this Rel Eff Curve
      if( is_null(constraint.controlling_source) || is_null(constraint.constrained_source) )
        throw runtime_error( "Invalid nuclide in ActRatioConstraint" ); 

      if( constraint.controlling_source == constraint.constrained_source )
        throw runtime_error( "Controlling and constrained sources cannot be the same" );

      if( constraint.constrained_to_controlled_activity_ratio <= 0.0 )
        throw runtime_error( "Constrained to controlled activity ratio must be greater than 0.0" );

      bool has_constrained = false, has_controlling = false;
      for( const NucInputInfo &src : nuclides )
      {
        assert( !is_null(src.source) );

        has_constrained |= (constraint.constrained_source == src.source);
        has_controlling |= (constraint.controlling_source == src.source);
      }   

      if( !has_constrained )
        throw runtime_error( "Constrained nuclide not found in nuclides" );

      if( !has_controlling )
        throw runtime_error( "Controlling nuclide not found in nuclides" );

      act_ratio_constraints.push_back( constraint );
    }//for( each ActRatioConstraint )
  }//if( act_ratio_constraints_node )

  // Get the MassFractionConstraints
  const rapidxml::xml_node<char> *mass_fraction_constraints_node = XML_FIRST_NODE( parent, "MassFractionConstraints" );
  if( mass_fraction_constraints_node )
  {
    XML_FOREACH_CHILD( constraint_node, mass_fraction_constraints_node, "MassFractionConstraint" )
    {
      MassFractionConstraint constraint;
      constraint.fromXml( constraint_node );

      bool has_nuclide = false;
      size_t num_nucs_of_el = 0;
      for( const NucInputInfo &nuc : nuclides )
      {
        const SandiaDecay::Nuclide *src_nuc = nuclide(nuc.source);
        has_nuclide |= (src_nuc == constraint.nuclide);
        if( src_nuc && (src_nuc->atomicNumber == constraint.nuclide->atomicNumber) )
          ++num_nucs_of_el;
      }   

      if( !has_nuclide )
        throw runtime_error( "MassFractionConstraint nuclide '" + constraint.nuclide->symbol + "' not found in problems nuclides." );

      if( num_nucs_of_el < 2 )
        throw runtime_error( "MassFractionConstraint for " + constraint.nuclide->symbol 
                                + " can only be used when there is two or more nuclides for an element." );

      mass_fraction_constraints.push_back( constraint );
    }//for( each MassFractionConstraint )
  }//if( mass_fraction_constraints_node )
}//void RelEffCurveInput::fromXml(...)
  
  

Options::Options()
: fit_energy_cal( false ),
  fwhm_form( FwhmForm::Polynomial_2 ),
  fwhm_estimation_method( FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum ),
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
  m_final_parameters{},
  m_final_uncertainties{},
  m_covariance{},
  m_parameter_names{},
  m_parameter_were_fit{},
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
  
  const RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef phys_in
           = RelActCalcAutoImp::RelActAutoCostFcn::make_phys_eqn_input( m_options.rel_eff_curves[rel_eff_index], m_drf, coeffs, 0 );
  
  return RelActCalc::physical_model_rel_eff_eqn_text( phys_in.self_atten, phys_in.external_attens,
                                phys_in.det, phys_in.hoerl_b, phys_in.hoerl_c, html_format );
}//std::string RelActAutoSolution::rel_eff_txt()
  

std::ostream &RelActAutoSolution::print_summary( std::ostream &out ) const
{
  char buffer[128] = { '\0' };

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
  //out << "Raw Ceres par values: [";
  //for( size_t i = 0; i < m_final_parameters.size(); ++i )
  //  out << ((i > 0) ? ", " : "") << "{" << m_parameter_names[i] << "=" << m_final_parameters[i]
  //      << "," << (m_parameter_were_fit[i] ? "Fit" : "NotFit") << "}";
  //out << "]\n\n";

  // Rel Eff code from RelEff
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
      const SandiaDecay::Nuclide * const nuc_i_nuclide = RelActCalcAuto::nuclide(nuc_i.source);

      for( size_t j = 0; j < i; ++j )
      {
        const NuclideRelAct &nuc_j = rel_acts[j];
        const SandiaDecay::Nuclide * const nuc_j_nuclide = RelActCalcAuto::nuclide(nuc_j.source);

        
        const double act_ratio = activity_ratio(nuc_i.source, nuc_j.source, rel_eff_index);
        const double nuc_mass_ratio = (nuc_i_nuclide && nuc_j_nuclide) ? mass_ratio(nuc_i_nuclide, nuc_j_nuclide, rel_eff_index) : 0.0;
        
        out << "\t" << std::setw(16) << (nuc_i.name() + " / " + nuc_j.name())
        << "\tact: " << std::setw(10) << std::setprecision(4) << act_ratio;
        try
        {
          const double uncert = activity_ratio_uncertainty(nuc_i.source, rel_eff_index, nuc_j.source, rel_eff_index );
          out << " ± " << uncert;
        }catch( std::exception &e )
        {
          cerr << "\nFailed to get act ratio uncert for "
          << RelActCalcAuto::to_name(nuc_i.source) << "/" << RelActCalcAuto::to_name(nuc_j.source)
          << " - " << e.what() << endl;
        }
        out << ",\tmass: " << std::setw(10) << std::setprecision(4) << nuc_mass_ratio;
        out << "\n";
        
        if( nuc_i_nuclide && nuc_j_nuclide )
        {
          out << "\t" << std::setw(16) << (nuc_j.name() + " / " + nuc_i.name())
          << "\tact: " << std::setw(10) << std::setprecision(4) << 1.0/act_ratio;
          try
          {
            const double uncert = activity_ratio_uncertainty(nuc_j.source, rel_eff_index, nuc_i.source, rel_eff_index );
            out << " ± " << uncert;
          }catch( std::exception &e )
          {
            cerr << "\nFailed to get act ratio uncert for "
            << RelActCalcAuto::to_name(nuc_j.source) << "/" << RelActCalcAuto::to_name(nuc_i.source)
            << " - " << e.what() << endl;
          }
          
          out << ",\tmass: " << std::setw(10) << std::setprecision(4) << 1.0/nuc_mass_ratio;
          out << "\n";
        }
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
      const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(nuc_rel_act.source);
      if( !nuc )
        continue;
      
      nucs_of_element[nuc->atomicNumber].push_back( nuc );
      
      if( (nuc->atomicNumber == 92) && (nuc->massNumber == 235) )
        u235 = nuc;
      
      if( (nuc->atomicNumber == 94) && (nuc->massNumber == 240) )
        pu240 = nuc;
    }//for( const NuclideRelAct &nuc : m_rel_activities )

    const auto print_enrich = [&out,this,rel_eff_index]( const SandiaDecay::Nuclide *const nuc ){
      try
      {
        char buffer[128] = { '\0' };
        const pair<double,optional<double>> enrich_frac = mass_enrichment_fraction(nuc,rel_eff_index);
        const double enrich = 100.0*enrich_frac.first;
        snprintf( buffer, sizeof(buffer), "%.2f", enrich );
        out << "\nEnrichment " << buffer << "% " << nuc->symbol;
        
        if( enrich_frac.second.has_value() )
        {
          const double neg_2sigma = 100.0*(enrich_frac.first - 2.0*enrich_frac.second.value());
          const double pos_2sigma = 100.0*(enrich_frac.first + 2.0*enrich_frac.second.value());
          
          snprintf( buffer, sizeof(buffer), " (2σ: %.2f%%, %.2f%%)", neg_2sigma, pos_2sigma );
          out << buffer;
        }
      }catch( std::exception &e )
      {
        // If covariance matrix couldnt be computed
        cerr << "Failed to get enrichment 2σ for " << nuc->symbol << ": " << e.what() << endl;
      }//try / catch

      out << "\n";
    };//print_enrich lambda

    if( u235 )
      print_enrich( u235 );

    if( pu240 )
      print_enrich( pu240 );

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
        const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(act.source);
        if( nuc && (nuc->atomicNumber == 94) )
        {
          nuclides.insert( act.name() );
        }
      }//for( const auto &act : rel_acts )
      
      if( nuclides.count( "Pu238" ) )
        out << "\tPu238: " << SpecUtils::printCompact(100.0*pu_corr->pu238_mass_frac, 4) << "\n";
      if( nuclides.count( "Pu239" ) )
        out << "\tPu239: " << SpecUtils::printCompact(100.0*pu_corr->pu239_mass_frac, 4) << "\n";
      if( nuclides.count( "Pu240" ) )
        out << "\tPu240: " << SpecUtils::printCompact(100.0*pu_corr->pu240_mass_frac, 4) << "\n";
      if( nuclides.count( "Pu241" ) )
        out << "\tPu241: " << SpecUtils::printCompact(100.0*pu_corr->pu241_mass_frac, 4) << "\n";
      out << "\tPu242 (by corr): " << SpecUtils::printCompact(100.0*pu_corr->pu242_mass_frac, 4) << "\n";
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
      {
        try
        {
          const pair<double,optional<double>> enrich_frac = mass_enrichment_fraction(nuc,rel_eff_index);
          const double frac = enrich_frac.first;
          out << std::setw(5) << nuc->symbol << ": "
          << std::setw(10) << std::setprecision(4) << (100.0*frac) << "%"
          << " of the " << el->name << ", by mass.\n";
        }catch( std::exception & )
        {
          assert( 0 ); //shouldnt throw, I dont thing
        }
      }
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
          if( outer_act.source != inner_act.source )
            continue;
            
          const double inner_activity = inner_act.rel_activity;
          const double outer_activity = outer_act.rel_activity;
          const double inner_counts = nuclide_counts(inner_act.source, inner_re_index);
          const double outer_counts = nuclide_counts(outer_act.source, outer_re_index);
          
          out << "Ratio of " << outer_act.name() << " between RE curves "
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
  
  assert( m_rel_eff_coefficients.empty() || (m_rel_eff_forms.size() == m_rel_eff_coefficients.size()) );
  assert( m_rel_eff_forms.empty() || (m_rel_eff_forms.size() == m_options.rel_eff_curves.size()) );
  

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

  const bool have_multiple_rel_eff = (m_rel_eff_forms.size() > 1);
    
  stringstream results_html;

  if( m_status == Status::Success )
  {
    results_html << "<div>&chi;<sup>2</sup>=" << m_chi2 << " and there were " << m_dof
    << " DOF (&chi;<sup>2</sup>/dof=" << m_chi2/m_dof << ")</div>\n";
  }else
  {
    string fail_reason;
    switch( m_status )
    {
      case Status::NotInitiated:
        fail_reason = "Not initiated";
        break;
      case Status::FailedToSetupProblem:
        fail_reason = "Failed to setup problem";
        break;
      case Status::FailToSolveProblem:
        fail_reason = "Failed to solve problem";
        break;
      case Status::UserCanceled:
        fail_reason = "User canceled";
        break;
      case Status::Success:
        assert( false );
        break;
    }

    results_html << "<div style=\"color: red; font-size: 150%;\">Calculation failed, " << fail_reason << ": " << m_error_message << "</div>\n";
  }//if( m_status == Status::Success ) /  else

  for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_coefficients.size(); ++rel_eff_index )
  {
    results_html << "<div class=\"releffeqn\">Rel. Eff. Eqn" 
    << (m_rel_eff_forms.size() > 1 ? (" " + std::to_string(rel_eff_index)) : string()) 
    << ": y = " << rel_eff_txt(true,rel_eff_index)
    << "</div>\n";
  }//for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )


  for( size_t rel_eff_index = 0; rel_eff_index < m_options.rel_eff_curves.size(); ++rel_eff_index )
  {
    if( (rel_eff_index >= m_rel_eff_forms.size())
        || (rel_eff_index >= m_rel_eff_coefficients.size())
        || (rel_eff_index >= m_rel_activities.size()) )
    {
      assert( m_status != Status::Success );
      break;
    }

    const RelEffCurveInput &rel_eff = m_options.rel_eff_curves[rel_eff_index];
    const vector<NuclideRelAct> &rel_activities = m_rel_activities[rel_eff_index];
    
    shared_ptr<const RelActCalc::Pu242ByCorrelationOutput> pu_corr = (rel_eff_index >= m_corrected_pu.size())
                                                                   ? nullptr : m_corrected_pu[rel_eff_index];
    
    // For Pu, print a corrected enrichment table
    if( pu_corr )
    {
      results_html << "<table class=\"nuctable resulttable\">\n"
      "  <caption>Plutonium mass fractions"
      << (have_multiple_rel_eff ? (" Rel. Eff. " + std::to_string(rel_eff_index)) : string())
      << "</caption>\n"
      "  <thead><tr>"
      " <th scope=\"col\">Nuclide</th>"
      " <th scope=\"col\">% Pu Mass</th>"
      " </tr></thead>\n"
      "  <tbody>\n";

      double pu239_act = 0.0; //We will need this to convert Pu242 rel mass, to a matching rel act
      set<string> pu_nuclides;
      set<double> pu_nuclide_ages;
      vector<tuple<const SandiaDecay::Nuclide *,double>> pu_rel_acts;
      for( const NuclideRelAct &act : rel_activities )
      {
        const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(act.source);
        if( nuc && (nuc->atomicNumber == 94) )
        {
          if( nuc->massNumber == 239 )
            pu239_act = act.rel_activity;
          pu_nuclides.insert( act.name() );
          pu_nuclide_ages.insert( act.age );
          pu_rel_acts.emplace_back( nuc, act.rel_activity );
        }
      }//for( const auto &act : rel_activities )
      
      if( pu_nuclides.count( "Pu238" ) )
        results_html << "  <tr><td>Pu238</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu238_mass_frac, 4)
        << "</td></tr>\n";
      if( pu_nuclides.count( "Pu239" ) )
        results_html << "  <tr><td>Pu239</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu239_mass_frac, 4)
        << "</td></tr>\n";
      if( pu_nuclides.count( "Pu240" ) )
        results_html << "  <tr><td>Pu240</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu240_mass_frac, 4)
        << "</td></tr>\n";
      if( pu_nuclides.count( "Pu241" ) )
        results_html << "  <tr><td>Pu241</td><td>"
        << SpecUtils::printCompact(100.0*pu_corr->pu241_mass_frac, 4)
        << "</td></tr>\n";
      results_html << "  <tr><td>Pu242 (by corr)</td><td>"
      << SpecUtils::printCompact(100.0*pu_corr->pu242_mass_frac, 4)
      << "</td></tr>\n";
      results_html << "  </tbody>\n"
      << "</table>\n\n";

      {//Begin put pu242 into `pu_rel_acts`
        const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
        assert( db );
        const SandiaDecay::Nuclide * const pu239 = db->nuclide("Pu239");
        const SandiaDecay::Nuclide * const pu242 = db->nuclide("Pu242");
        assert( pu239 && pu242 );
        const double pu242_rel_act = pu239_act * (pu_corr->pu242_mass_frac * pu242->activityPerGram())
                                                     / (pu_corr->pu239_mass_frac * pu239->activityPerGram());

        pu_rel_acts.emplace_back( pu242, pu242_rel_act );
      }//End put pu242 into `pu_rel_acts`

      if( (pu_nuclide_ages.size() == 1) && (pu_nuclides.size() > 1) )
      {
        const double back_decay_time = *begin(pu_nuclide_ages);
        vector<tuple<const SandiaDecay::Nuclide *,double,double>> nuc_act_massfrac_at_t0
                                          = RelActCalc::back_decay_relative_activities( back_decay_time, pu_rel_acts );
        std::sort( begin(nuc_act_massfrac_at_t0), end(nuc_act_massfrac_at_t0), []( const auto &lhs, const auto &rhs ) -> bool {
          return (get<0>(lhs)->symbol < get<0>(rhs)->symbol);
        } );

        results_html << "<table class=\"nuctable resulttable\">\n"
        "  <caption>Plutonium mass fractions at T=0"
        << (have_multiple_rel_eff ? (" Rel. Eff. " + std::to_string(rel_eff_index)) : string())
        << "</caption>\n"
        "  <thead><tr>"
        " <th scope=\"col\">Nuclide</th>"
        " <th scope=\"col\">% Pu Mass</th>"
        " </tr></thead>\n"
        "  <tbody>\n";
        for( tuple<const SandiaDecay::Nuclide *,double,double> nuc_act_mass : nuc_act_massfrac_at_t0 )
        {
          results_html << "  <tr><td>" << get<0>(nuc_act_mass)->symbol << "</td><td>"
          << SpecUtils::printCompact(100.0*get<2>(nuc_act_mass), 4)
          << "</td></tr>\n";
        }
        results_html << "  </tbody>\n"
        << "</table>\n\n";
      }//if( (pu_nuclide_ages.size() == 1) && (pu_nuclides.size() > 1) )
    }//if( pu_corr )


    double sum_rel_mass = 0.0;
    for( const auto &act : rel_activities )
    {
      const SandiaDecay::Nuclide * const nuc = nuclide(act.source);
      if( nuc )
        sum_rel_mass += act.rel_activity / nuc->activityPerGram();
    }

    
    results_html << "<table class=\"nuctable resulttable\">\n"
    "  <caption>Relative activities and mass fractions"
    << (m_rel_eff_forms.size() > 1 ? (" Rel. Eff. " + std::to_string(rel_eff_index)) : string())
    << "</caption>\n"
    "  <thead><tr>"
    " <th scope=\"col\">Nuclide</th>"
    " <th scope=\"col\">Rel. Act.</th>"
    " <th scope=\"col\">Total Mas Frac.</th>"
    " <th scope=\"col\">Enrichment</th>"
    " <th scope=\"col\">Enrich 2&sigma;</th>"
    " <th scope=\"col\">Det. Counts</th>"
    " </tr></thead>\n"
    "  <tbody>\n";
    
    
    for( const NuclideRelAct &act : rel_activities )
    {
      results_html << "  <tr><td>" << act.name();
      
      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(act.source);
      if( nuc
         && (nuc->atomicNumber == 94)
         && (nuc->massNumber == 242)
         && (rel_eff.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable) )
      {
        results_html << " (by corr)";
      }
      
      results_html << "</td>"
      << "<td>" << act.rel_activity << " &plusmn; " << act.rel_activity_uncertainty << "</td>";

      if( nuc )
      {
        const double rel_mass = act.rel_activity / nuc->activityPerGram();
        results_html << "<td>" << (100.0*rel_mass/sum_rel_mass) << "%</td>";
        try
        {
          const pair<double,optional<double>> enrich_frac = mass_enrichment_fraction(nuc,rel_eff_index);
          results_html << "<td>" << (100.0*enrich_frac.first) << "%</td>";
          if( enrich_frac.second.has_value() )
          {
            const double minus_2sigma = 100.0*(enrich_frac.first - 2.0*enrich_frac.second.value());
            const double plus_2sigma = 100.0*(enrich_frac.first - 2.0*enrich_frac.second.value());
            results_html << "<td>" << minus_2sigma << "%, " << plus_2sigma << "%</td>";
          }else
          {
            results_html << "<td>--</td>";
          }
        }catch( std::exception &e )
        {
          results_html << "<td>--</td><td>--</td>";
        }
      }else
      {
        results_html << "<td></td>" << "<td></td>" << "<td></td>";
      }
      
      const double nuc_counts = nuclide_counts( act.source, rel_eff_index );
      results_html << "<td>" << nuc_counts << "</td>";
      results_html << "</tr>\n";
    }//for( const auto &act : rel_activities )
    
    results_html << "  </tbody>\n"
    << "</table>\n\n";
  }//for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )

  for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )
  {
    if( (rel_eff_index >= m_rel_eff_forms.size())
        || (rel_eff_index >= m_rel_eff_coefficients.size())
        || (rel_eff_index >= m_rel_activities.size()) )
    {
      assert( m_status != Status::Success );
      break;
    }

    const RelEffCurveInput &rel_eff = m_options.rel_eff_curves[rel_eff_index];
    const vector<NuclideRelAct> &rel_activities = m_rel_activities[rel_eff_index];
    
    // Make the table of mass and activity ratios
    results_html << "<table class=\"massratiotable resulttable\">\n";
    results_html << "  <caption>Mass and Activity Ratios."
    << (m_rel_eff_forms.size() > 1 ? (" Rel. Eff. " + std::to_string(rel_eff_index)) : string())
    << "</caption>\n"
    "  <thead><tr>"
    "<th scope=\"col\">Nuclides</th>"
    "<th scope=\"col\">Mass Ratio</th>"
    "<th scope=\"col\">Activity Ratio</th>"
    "<th scope=\"col\">Uncertainty</th>"
    "</tr></thead>\n"
    "  <tbody>\n";
    for( size_t i = 1; i < rel_activities.size(); ++i )
    {
      for( size_t j = 0; j < i; ++j )
      {
        const NuclideRelAct &nuc_i = rel_activities[i];
        const NuclideRelAct &nuc_j = rel_activities[j];

        assert( !is_null(nuc_i.source) && !is_null(nuc_j.source) );
        if( is_null(nuc_i.source) || is_null(nuc_j.source) )
          continue;

        const SandiaDecay::Nuclide * const nuc_i_nuc = nuclide(nuc_i.source);
        const SandiaDecay::Nuclide * const nuc_j_nuc = nuclide(nuc_j.source);

        const double act_i = nuc_i.rel_activity;
        const double act_j = nuc_j.rel_activity;
        const std::string name_i = to_name(nuc_i.source);
        const std::string name_j = to_name(nuc_j.source);

        string mass_i_j = "--", mass_j_i = "--";
        if( nuc_i_nuc && nuc_j_nuc )
        {
          const double mass_i = act_i / nuc_i_nuc->activityPerGram();
          const double mass_j = act_j / nuc_j_nuc->activityPerGram();

          snprintf( buffer, sizeof(buffer), "%1.6G", (mass_i / mass_j) );
          mass_i_j = buffer;

          snprintf( buffer, sizeof(buffer), "%1.6G", (mass_j / mass_i) );
          mass_j_i = buffer;
        }

        string uncert_i_j_str = "--", uncert_j_i_str = "--";
        try
        {
          const double uncert_i_j = activity_ratio_uncertainty( nuc_i.source, rel_eff_index, nuc_j.source, rel_eff_index );
          const double uncert_j_i = activity_ratio_uncertainty( nuc_j.source, rel_eff_index, nuc_i.source, rel_eff_index );

          snprintf( buffer, sizeof(buffer), "%1.6G%%", 100.0*uncert_i_j/(act_i/act_j) );
          uncert_i_j_str = buffer;

          snprintf( buffer, sizeof(buffer), "%1.6G%%", 100.0*uncert_j_i/(act_j/act_i) );
          uncert_j_i_str = buffer;
        }catch( std::exception &e )
        {
          // We get here if there was an issue computing the covatriance
        }//try / catch

        snprintf( buffer, sizeof(buffer), "<tr><td>%s</td><td>%s</td><td>%1.6G</td><td>%s</td></tr>\n",
                 (name_i + "/" + name_j).c_str(), mass_i_j.c_str(), (act_i / act_j), uncert_i_j_str.c_str() );
        results_html << buffer;
        
        snprintf( buffer, sizeof(buffer), "<tr><td>%s</td><td>%s</td><td>%1.6G</td><td>%s</td></tr>\n",
                 (name_j + "/" + name_i).c_str(), mass_j_i.c_str(), (act_j / act_i), uncert_j_i_str.c_str() );

        results_html << buffer;
      }//for( size_t j = 0; j < i; ++j )
    }//for( size_t i = 0; i < used_isotopes.size(); ++i )
    results_html << "  </tbody>\n"
    << "</table>\n\n";
  }//for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )



  // Make table giving info on each of the result peaks
  results_html << "<table class=\"peaktable resulttable\">\n";
  results_html << "  <caption>Peaks fit in analysis.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Energy (keV)</th>"
  "<th scope=\"col\">Nuclide</th>"
  "<th scope=\"col\">Yield</th>"
  << (have_multiple_rel_eff ? "<th scope=\"col\">Rel. Eff. #</th>" : "")
  << "<th scope=\"col\">Amp</th>"
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
    const double energy = info.mean();
    SrcVariant peak_src;
    const SandiaDecay::Nuclide * const nuc = info.parentNuclide();
    const SandiaDecay::Element * const el = info.xrayElement();
    const ReactionGamma::Reaction * const reaction = info.reaction();
    if( nuc )
      peak_src = nuc;
    else if( el )
      peak_src = el;
    else if( reaction )
      peak_src = reaction;
    assert( !is_null(peak_src) );
    if( is_null(peak_src) )
      throw logic_error( "Peak with null src not in solution?" );

    int rel_eff_index = 0;
    if( have_multiple_rel_eff )
    {
      string label = info.userLabel();
      const bool is_rel_eff_label = SpecUtils::istarts_with( label, "RelEff " );
      assert( is_rel_eff_label );
      if( is_rel_eff_label )
      {
        label = label.substr(7);
        const bool ok = SpecUtils::parse_int( label.c_str(), label.size(), rel_eff_index );
        assert( ok );
        assert( rel_eff_index >= 0 && rel_eff_index < static_cast<int>(m_rel_eff_forms.size()) );
        if( rel_eff_index < 0 || rel_eff_index >= static_cast<int>(m_rel_eff_forms.size()) )
          rel_eff_index = 0;
      }
    }//if( have_multiple_rel_eff )

    assert( (rel_eff_index >= 0) && (rel_eff_index < static_cast<int>(m_rel_eff_forms.size())) );
    const RelEffCurveInput &rel_eff = m_options.rel_eff_curves[rel_eff_index];
    const vector<NuclideRelAct> &rel_activities = m_rel_activities[rel_eff_index];
    
    const NuclideRelAct *nuc_info = nullptr;
    for( const NuclideRelAct &rel_act : rel_activities )
      nuc_info = (rel_act.source == peak_src) ? &rel_act : nuc_info;
    
    assert( nuc_info );
    if( !nuc_info )
      throw logic_error( "Peak with nuc not in solution?" );
      

    shared_ptr<const Material> self_atten;
    if( rel_eff.phys_model_self_atten )
      self_atten = rel_eff.phys_model_self_atten->material;
    vector<shared_ptr<const Material>> external_attens;
    for( const auto &ext : rel_eff.phys_model_external_atten )
      external_attens.push_back( ext->material );
    
    assert( !self_atten || (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) );
    assert( external_attens.empty() || (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) );

    if( !nuc || !nuc_info )
    {
      const char *rel_eff_label = (have_multiple_rel_eff ? "<td></td>" : "");
      
      snprintf( buffer, sizeof(buffer),
                 "  <tr><td>%.2f</td><td></td><td></td>%s<td>%1.6G</td><td>%1.6G</td>"
                 "<td></td><td></td><td></td><td></td><td></td><td></td></tr>\n",
                 energy, rel_eff_label, info.amplitude(), info.amplitudeUncert() );
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
          
        const RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef input
        = RelActCalcAutoImp::RelActAutoCostFcn::make_phys_eqn_input( rel_eff, m_drf, m_rel_eff_coefficients[rel_eff_index], 0 );
          
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
            const RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef input
              = RelActCalcAutoImp::RelActAutoCostFcn::make_phys_eqn_input( rel_eff, m_drf, rel_eff_pars, 0 );
              
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
        
      std::string format_str = "  <tr>"s
                 "<td>%.2f</td>"    // energy
                 "<td>%s</td>"      // nuclide
                 "<td>%1.3G</td>"s   // yield
                 + (have_multiple_rel_eff ? ("<td>" + std::to_string(rel_eff_index) + "</td>") : ""s)
                 + "<td>%1.3G</td>"   // amplitude
                 "<td>%1.2G%%</td>" // amplitude uncertainty
                 "<td>%1.3G</td>"   // cps/yield
                 "<td>%1.6G</td>"   // fit rel eff
                 "<td>%1.2G%%</td>" // fit rel eff uncertainty
                 "<td>%s</td>"      // continuum type
                 "<td>%.1f-%.1f</td>" // continuum range
                 "</tr>\n"s;

      snprintf(buffer, sizeof(buffer), format_str.c_str(),
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

    for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )
    {
      const RelEffCurveInput &rel_eff = m_options.rel_eff_curves[rel_eff_index];
      results_html << "  <tr><th scope=\"row\">nucs_of_el_same_age</th><td><code>" << rel_eff.nucs_of_el_same_age << "</code></td></tr>\n";
      results_html << "  <tr><th scope=\"row\">rel_eff_eqn_type</th><td><code>"
      << RelActCalc::to_str(rel_eff.rel_eff_eqn_type) << "</code></td></tr>\n";

      if( rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        // TODO: add physical model info here
      }else
      {
        results_html << "  <tr><th scope=\"row\">rel_eff_eqn_order</th><td><code>" << rel_eff.rel_eff_eqn_order << "</code></td></tr>\n";
      }
    }//for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )
    
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
    
  vector<RelEffChart::ReCurveInfo> rel_eff_info_sets;
  
  assert( m_fit_peaks_for_each_curve.size() == m_rel_eff_forms.size() );
  
  for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )
  {
    const RelEffCurveInput &rel_eff = m_options.rel_eff_curves[rel_eff_index];
  
    RelEffChart::ReCurveInfo info;

    info.live_time = m_spectrum ? m_spectrum->live_time() : 1.0;
    info.fit_peaks = m_fit_peaks_for_each_curve[rel_eff_index];
    info.rel_acts = m_rel_activities[rel_eff_index];
    info.re_curve_name = Wt::WString::fromUTF8(rel_eff.name);

    try
    {
      info.re_curve_eqn_txt = Wt::WString::fromUTF8( "y = " + rel_eff_txt( false, rel_eff_index ) ); 

      if( rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        info.js_rel_eff_eqn = RelActCalc::rel_eff_eqn_js_function( rel_eff.rel_eff_eqn_type, m_rel_eff_coefficients[rel_eff_index] );
      }else
      {
        const RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef input
                   = RelActCalcAutoImp::RelActAutoCostFcn::make_phys_eqn_input( rel_eff, m_drf, m_rel_eff_coefficients[rel_eff_index], 0 );
      
        info.js_rel_eff_eqn = RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                                                                          input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
      }//if( not FramPhysicalModel ) / else
    }catch( std::exception &e )
    {
      cerr << "Error creating RelEffChart::ReCurveInfo: " << e.what() << endl;
    }

    rel_eff_info_sets.push_back( info );
  }//for( size_t rel_eff_index = 0; rel_eff_index < m_rel_eff_forms.size(); ++rel_eff_index )

  
  const string rel_eff_json = RelEffChart::jsonForData( rel_eff_info_sets );

  SpecUtils::ireplace_all( html, "${REL_EFF_JSON}", rel_eff_json.c_str() );
  
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



pair<double,optional<double>> RelActAutoSolution::mass_enrichment_fraction( const SandiaDecay::Nuclide *nuclide, const size_t rel_eff_index ) const
{
#define USE_RelActAutoCostFcn_MASS_FRAC_IMP 1
  
#if( USE_RelActAutoCostFcn_MASS_FRAC_IMP )
  
  if( !m_cost_functor )
    throw std::runtime_error( "RelActAutoSolution::mass_enrichment_fraction(): cost_functor not set." );
  
  const vector<double> &x = m_final_parameters;
  const vector<vector<double>> &covariance = m_covariance;
  
  if( covariance.empty() )
  {
    const double frac = m_cost_functor->mass_enrichment_fraction( nuclide, rel_eff_index, x );
    return { frac, optional<double>{} };
  }
  
  const size_t num_pars = m_cost_functor->number_parameters();
  assert( covariance.size() == num_pars );
  if( covariance.size() != num_pars )
    throw std::logic_error( "mass_enrichment_fraction: covariance matrix has wrong number of rows." );
  for( const auto &row : covariance )
  {
    assert( row.size() == num_pars );
    if( row.size() != num_pars )
      throw std::logic_error( "mass_enrichment_fraction: covariance matrix has wrong number of columns." );
  }
  
  
  // We will use `J * cov * J^T` (where J is jacobian - e.g, d{MassFrac(nuclide)}/d{paramater}) to comput the
  //  uncertainty, using auto-differentiation to compute the Jacobian vector.  Only the elements corresponding
  //  to the activities of the other nuclides of this element will be non-zero in principle; however, we will
  //  compute all element, because there may be some activity ratio, or mass-ratio constraints that would
  //  complicate things, so we'll just be consistent.
  double enrichment = -1.0;
  vector<double> jacobian( num_pars, 0.0 );
  for( size_t i = 0; i < num_pars; i += RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size )
  {
    vector<ceres::Jet<double,RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size>> x_local( begin(x), end(x) );
    for( size_t j = 0; (j < RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size) && (i+j < num_pars); ++j )
      x_local[i+j].v[j] = 1.0;
    const ceres::Jet<double,RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size> enrichment_jet
                      = m_cost_functor->mass_enrichment_fraction( nuclide, rel_eff_index, x_local );
    for( size_t j = 0; (j < RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size) && ((i+j) < num_pars); ++j )
      jacobian[i+j] = enrichment_jet.v[j];
    enrichment = enrichment_jet.a;
  }
    
  double uncertainty = 0.0;
  for( size_t i = 0; i < num_pars; ++i )
  {
    for( size_t j = 0; j < num_pars; ++j )
      uncertainty += jacobian[i]*covariance[i][j]*jacobian[j];
  }
    
  uncertainty = sqrt( uncertainty );
    
  return {enrichment, uncertainty};
#else
  static_assert( 0, "RelActAutoSolution::mass_enrichment_fraction(...): this implementation needs a little work to return both nominal answer, as well as uncertainty." );
  // This is based upon `RelEffSolution::mass_fraction( const std::string &nuclide, const double num_sigma ) const`
  //  If any issues are found in this function, please also check that function.
  const double num_sigma = 1.0;
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  if( !nuclide ) //Will be nullptr for reactions and x-rays
    throw runtime_error( "RelActAutoSolution::mass_enrichment_fraction(nullptr, " + std::to_string(rel_eff_index)
                         + "): invalid nuclide" );

  if( m_phys_units_cov.empty() && (num_sigma != 0.0) )
    throw runtime_error( "RelActAutoSolution::mass_enrichment_fraction(" + nuclide->symbol + ", "
                      + std::to_string(rel_eff_index)  + ", " + std::to_string(num_sigma) + "): no valid covariance" );

  const NuclideRelAct &nuc_info = nucinfo( nuclide, rel_eff_index );

  assert( rel_eff_index < m_options.rel_eff_curves.size() );
  const RelEffCurveInput &re_curve = m_options.rel_eff_curves.at(rel_eff_index);

  const bool is_pu = (nuclide->atomicNumber == 94);
  const bool using_pu242_corr = (is_pu && (m_corrected_pu.size() > rel_eff_index) && m_corrected_pu[rel_eff_index]);
  if( using_pu242_corr && (num_sigma == 0.0) )
  {
    const vector<NuclideRelAct> &rel_acts = m_rel_activities[rel_eff_index];
    const size_t nuc_index = nuclide_index( nuclide, rel_eff_index );
    assert( RelActCalcAuto::nuclide(rel_acts[nuc_index].source) == nuclide );
    const double rel_mass = rel_acts[nuc_index].rel_activity / nuclide->activityPerGram();

    double el_total_mass = 0.0;
    for( const NuclideRelAct &nuc : rel_acts )
    {
      const SandiaDecay::Nuclide * const nuc_nuclide = RelActCalcAuto::nuclide(nuc.source);
      if( nuc_nuclide && (nuc_nuclide->atomicNumber == nuclide->atomicNumber) )
        el_total_mass += nuc.rel_activity / nuc_nuclide->activityPerGram();
    }

    const shared_ptr<const RelActCalc::Pu242ByCorrelationOutput> &pu_corr = m_corrected_pu[rel_eff_index];

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

    assert( (rel_mass >= 0.0) && (rel_mass <= el_total_mass) );

    return rel_mass / el_total_mass;
  }//if( is pu, and using Pu correction, and nominal answer is wanted )


 // If this nuclide was constrained, then we need to find the ultimate controlling nuclide.
  double nuc_mult = 1.0;
  SrcVariant src( nuclide );
  const bool nuc_was_constrained = walk_to_controlling_nuclide( src, rel_eff_index, nuc_mult );

 #pragma message( "RelActAutoSolution: Need to check if we are computing mass_fraction with uncertainties correctly when there are constraints (looks correct with a simple example)." )
  if( nuc_was_constrained )
    cerr << "RelActAutoSolution: Need to check if we are computing mass_fraction with uncertainties correctly when there are constraints." << endl;


  // The Covariance matrix is in terms of fit paramaters - so we need to multiple by the nuclides mutliepl to get to activity.
  const size_t nuc_act_par_index = fit_parameters_index_for_source( src, rel_eff_index );
  assert( nuc_act_par_index < m_parameter_scale_factors.size() );

  const double cov_nuc_nuc = m_phys_units_cov[nuc_act_par_index][nuc_act_par_index];
  const double sqrt_cov_nuc_nuc = sqrt(cov_nuc_nuc);

#ifndef NDEBUG
  // Check that relative activity uncertainties have been computed compatible with what we are
  //  assuming here (and no funny business has happened).
  if( fabs(sqrt_cov_nuc_nuc - nuc_info.rel_activity_uncertainty) > 1.0E-6*(std::max)(sqrt_cov_nuc_nuc, nuc_info.rel_activity_uncertainty) )
  {
    cout << "sqrt_cov_nuc_nuc = " << sqrt_cov_nuc_nuc << " (="<< sqrt_cov_nuc_nuc << ")" << endl;
    cout << "m_rel_activities[nuc_index].m_rel_activity_uncert = " << nuc_info.rel_activity_uncertainty << endl;
    cout << "m_rel_activities[nuc_index].m_rel_activity = " << nuc_info.rel_activity << endl;
  }
  assert( fabs(sqrt_cov_nuc_nuc - nuc_info.rel_activity_uncertainty) < 1.0E-6*(std::max)(nuc_info.rel_activity_uncertainty, sqrt_cov_nuc_nuc) );
#endif 

  double sum_rel_mass = 0.0, nuc_rel_mas = -1.0;

  set<double> pu_ages;//only used if using_pu242_corr
  RelActCalc::Pu242ByCorrelationInput raw_pu_masses; //only used if using_pu242_corr


#define TRY_ALT_ENRICH_UNCERT_METHOD 0

#if( TRY_ALT_ENRICH_UNCERT_METHOD )
  //I think the currently implemented way of getting mass fraction uncertainty is not forrect, instead, if we
  //  think of enrichment fraction as:
  //let y = (RelAct[nuclide]/c_nuc) / ( (RelAct[0]/c_0) + (RelAct[1]/c_1) + ... );
  //  where c_0 is el_nuc_0->activityPerGram(), c_1 is el_nuc_1->activityPerGram(), etc
  //Then we can get the Jacobian J=[dy/dRelAct0, dy/dRelAct1, ...]
  //  and then uncert_squared = J*[Cov]*J^T
  //I tried to do this in this section - but its not quite working out yet (may give negative uncert2), so leaving not working for now

  vector<double> rel_masses, controlled_rel_act_mults, rel_act_scale_factors, act_per_grams;
  vector<size_t> src_act_par_indexes;
  size_t src_act_par_index = m_rel_activities[rel_eff_index].size(), src_index = m_rel_activities[rel_eff_index].size();
  double sum_nominal_rel_masses = 0.0;
#endif

  //for( size_t rel_eff_loop_index = 0; rel_eff_loop_index < m_rel_activities.size(); ++rel_eff_loop_index )
  {
    const size_t rel_eff_loop_index = rel_eff_index;
    const vector<NuclideRelAct> &rel_acts = m_rel_activities[rel_eff_loop_index];
    for( size_t rel_act_loop_index = 0; rel_act_loop_index < rel_acts.size(); ++rel_act_loop_index )
    {
      const NuclideRelAct &act = rel_acts[rel_act_loop_index];
      SrcVariant act_src = act.source;
      const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(act_src);
      if( !nuc || (nuc->atomicNumber != nuclide->atomicNumber) )
        continue;

      double loop_nuc_mult = 1.0;
      walk_to_controlling_nuclide( act_src, rel_eff_loop_index, loop_nuc_mult );

      const size_t act_par_index = fit_parameters_index_for_source( act_src, rel_eff_loop_index );

      const double rel_act_scale_factor = m_parameter_scale_factors[act_par_index];
      double fit_act_for_index = loop_nuc_mult * rel_act_scale_factor * m_final_parameters[act_par_index];

      // Need to check if mass-fraction constrained nuc, and if so, correct `fit_act_for_index`.
      const RelActCalcAuto::RelEffCurveInput &input = m_options.rel_eff_curves[rel_eff_index];
      const vector<RelEffCurveInput::MassFractionConstraint> &mass_frac_consts = input.mass_fraction_constraints;
      for( const RelEffCurveInput::MassFractionConstraint &constraint : mass_frac_consts )
      {
        if( constraint.nuclide == nuc )
        {
          fit_act_for_index = rel_activity( act_src, rel_eff_index );
          break;
        }//if( constraint.nuclide == nuc )
      }//for( const MassFractionConstraint &constraint : mass_frac_consts )

      assert( fabs(fit_act_for_index - rel_activity(act_src, rel_eff_index)) < 1.0E-6*fit_act_for_index );

      const double cov_nuc_index = loop_nuc_mult*m_phys_units_cov[nuc_act_par_index][act_par_index]; //TODO: this may need another loop_nuc_mult type thing multiplied
      const double varied_fit_act_for_index = fit_act_for_index + (cov_nuc_index / cov_nuc_nuc) * num_sigma * sqrt_cov_nuc_nuc;

      const double rel_act = varied_fit_act_for_index;
      const double rel_mass = rel_act / nuc->activityPerGram();

#if( TRY_ALT_ENRICH_UNCERT_METHOD )
      if( nuc == nuclide )
      {
        assert( (src_act_par_index == rel_acts.size()) && (src_index == rel_acts.size()) );

        src_index = rel_masses.size();
        src_act_par_index = act_par_index;
      }

      const double nominal_rel_mass = fit_act_for_index / nuc->activityPerGram();
      rel_masses.push_back( nominal_rel_mass );
      src_act_par_indexes.push_back( act_par_index );
      controlled_rel_act_mults.push_back( loop_nuc_mult );
      rel_act_scale_factors.push_back( rel_act_scale_factor );
      act_per_grams.push_back( nuc->activityPerGram() );
      sum_nominal_rel_masses += nominal_rel_mass;
#endif

      if( using_pu242_corr && (nuc->atomicNumber == 94) )
      {
        assert( act.age >= 0.0 );
        if( act.age >= 0.0 )
          pu_ages.insert( act.age );

        switch( nuc->massNumber )
        {
          case 238: raw_pu_masses.pu238_rel_mass = rel_mass; break;
          case 239: raw_pu_masses.pu239_rel_mass = rel_mass; break;
          case 240: raw_pu_masses.pu240_rel_mass = rel_mass; break;
          case 241: raw_pu_masses.pu241_rel_mass = rel_mass; break;
          case 242:
            assert( 0 );
            throw std::logic_error( "Pu242 cant be fit nuclide Pu242 correlation correction method was specified." );
            break;
          default:  raw_pu_masses.other_pu_mass = rel_mass;  break;
        }//switch( nuclide->massNumber )
      }//if( using_pu242_corr && (nuc->atomicNumber == 94) )

      sum_rel_mass += (std::max)( rel_mass, 0.0 );
      //cout << "mass_frac(" << nuclide->symbol << "): for nuc " << nuc->symbol << " nom_act=" << fit_act_for_index
      //<< ", varied=" << varied_fit_act_for_index << ", cov_nuc_index=" << cov_nuc_index << ", cov_nuc_nuc=" << cov_nuc_nuc
      //<< ", num_sigma=" << num_sigma << ", sqrt_cov_nuc_nuc=" << sqrt_cov_nuc_nuc << endl;

      if( (rel_eff_loop_index == rel_eff_index) && (nuclide == nuc) )
        nuc_rel_mas = rel_mass;
    }//for( size_t rel_act_loop_index = 0; rel_act_loop_index < rel_acts.size(); ++rel_act_loop_index )
  }//for( size_t rel_eff_loop_index = 0; rel_eff_loop_index < m_rel_activities.size(); ++rel_eff_loop_index )

#if( TRY_ALT_ENRICH_UNCERT_METHOD )
  assert( src_index < m_rel_activities[rel_eff_index].size() );
  vector<double> jacobians( src_act_par_indexes.size(), 0.0 );
  for( size_t row = 0; row < src_act_par_indexes.size(); ++row )
  {
    const double rel_mass = rel_masses[row];
    const double rel_act_scale_factor = rel_act_scale_factors[row];
    const double loop_nuc_mult = controlled_rel_act_mults[row];

    if( src_index == row )
      jacobians[row] = (sum_nominal_rel_masses - rel_mass) / (sum_nominal_rel_masses*sum_nominal_rel_masses * act_per_grams[row]);
    else
      jacobians[row] = -rel_masses[src_index] / (sum_nominal_rel_masses*sum_nominal_rel_masses* act_per_grams[row] );
  }//for( size_t i = 0; i < src_act_par_indexes.size(); ++i )

  double uncert_sq = 0.0;
  for( size_t i = 0; i < src_act_par_indexes.size(); ++i )
  {
    for( size_t j = 0; j < src_act_par_indexes.size(); ++j )
      uncert_sq += jacobians[i]*jacobians[j]*m_phys_units_cov[src_act_par_indexes[i]][src_act_par_indexes[j]];
  }

  const double uncert = sqrt( uncert_sq );
  cout << "Using updated uncert method, we would get: " << rel_masses[src_index] << " + " << num_sigma << "*" << uncert << endl;
  return rel_masses[src_index] + num_sigma*uncert;
#endif

  // If we are doing Pu242 correction, we have to recompute the correction, based on the varied mass fractions
  if( using_pu242_corr )
  {
    const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
    assert( rel_eff_curve.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable );

    assert( num_sigma != 0.0 );  //We already dealt with the nominal case above.

    if( pu_ages.size() == 0 )
    {
      assert( sum_rel_mass == 0.0 );
      raw_pu_masses.pu_age = 0.0; //no correction for age - but we shouldnt get here anyway
    }else if( pu_ages.size() == 1 )
    {
      raw_pu_masses.pu_age = *begin(pu_ages);
    }else
    {
      const vector<double> ages( begin(pu_ages), end(pu_ages) );
      raw_pu_masses.pu_age = ages[ages.size() / 2]; //just take the median age...
    }

    // We dont have to divide `pu_total_mass`, but we will, just for debuging.
    raw_pu_masses.pu238_rel_mass /= sum_rel_mass;
    raw_pu_masses.pu239_rel_mass /= sum_rel_mass;
    raw_pu_masses.pu240_rel_mass /= sum_rel_mass;
    raw_pu_masses.pu241_rel_mass /= sum_rel_mass;
    raw_pu_masses.other_pu_mass  /= sum_rel_mass;

    const RelActCalc::Pu242ByCorrelationOutput corr_output
           = RelActCalc::correct_pu_mass_fractions_for_pu242( raw_pu_masses, rel_eff_curve.pu242_correlation_method );

    switch( nuclide->massNumber )
    {
      case 238: return corr_output.pu238_mass_frac;
      case 239: return corr_output.pu239_mass_frac;
      case 240: return corr_output.pu240_mass_frac;
      case 241: return corr_output.pu241_mass_frac;
      case 242: return corr_output.pu242_mass_frac;
      default:
        assert( 0 );
        throw runtime_error( "Unhandled Pu isotope for mass fraction amount" );
    }//switch( nuclide->massNumber )
  }//if( using_pu242_corr )


  if( nuc_rel_mas < 0.0 ) // This can happen when we go down a couple sigma
    return 0.0;

  //cout << "nuc_rel_mas = " << nuc_rel_mas << ", sum_rel_mass = " << sum_rel_mass << endl;
  //cout << "nuc_rel_mas / sum_rel_mass = " << nuc_rel_mas / sum_rel_mass << endl;
  return nuc_rel_mas / sum_rel_mass;
#endif //#if( USE_RelActAutoCostFcn_MASS_FRAC_IMP ) / else
}//mass_enrichment_fraction


double RelActAutoSolution::mass_ratio( const SandiaDecay::Nuclide *numerator,
                                      const SandiaDecay::Nuclide *denominator,
                                      const size_t rel_eff_index ) const
{
  const double ratio = activity_ratio(SrcVariant(numerator), SrcVariant(denominator), rel_eff_index);
  
  return ratio * denominator->activityPerGram() / numerator->activityPerGram();
}//double mass_ratio(...)


const NuclideRelAct &RelActAutoSolution::nucinfo( const SrcVariant src, const size_t rel_eff_index ) const
{
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "nucinfo: invalid rel eff index" );

  const vector<NuclideRelAct> &nuclides = m_rel_activities[rel_eff_index];

  for( size_t nuc_index = 0; nuc_index < nuclides.size(); ++nuc_index )
  {
    const NuclideRelAct &nucinfo = nuclides[nuc_index];
    if( nucinfo.source == src )
      return nucinfo;
  }//for( size_t nuc_index = 0; nuc_index < nuclides.size(); ++nuc_index )

  throw runtime_error( "nucinfo: nuclide not found in rel eff curve" );
}

double RelActAutoSolution::activity_ratio( const SrcVariant &numerator,
                                          const SrcVariant &denominator,
                                          const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "activity_ratio: invalid rel eff index" );
  
  return rel_activity(numerator, rel_eff_index) / rel_activity(denominator, rel_eff_index);
}//double activity_ratio(...)


size_t RelActAutoSolution::fit_parameters_index_for_source( const SrcVariant &src, const size_t rel_eff_index ) const
{
  // Calculate the starting index of activity parameters for this rel_eff_index
  // This follows the same structure as in the RelActAutoCostFcn class
  size_t acts_start_index = RelActAutoSolution::sm_num_energy_cal_pars; // Energy cal parameters
  
  // Add FWHM parameters
  acts_start_index += num_parameters( m_options.fwhm_form );
  
  // Add relative efficiency parameters for all curves up to this one
  for( size_t i = 0; i < m_options.rel_eff_curves.size(); ++i )
  {
    const auto &rel_eff_curve = m_options.rel_eff_curves[i];
    if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      acts_start_index += (2 + 2*rel_eff_curve.phys_model_external_atten.size() + 2); // AN + AD + external shields + hoerl b,c
    else
      acts_start_index += (rel_eff_curve.rel_eff_eqn_order + 1);
  }
  
  // Add activity parameters for all nuclides in curves before this one
  for( size_t i = 0; i < rel_eff_index; ++i )
  {
    acts_start_index += 2 * m_options.rel_eff_curves[i].nuclides.size();
  }

  const size_t nuc_index = nuclide_index( src, rel_eff_index );
  return acts_start_index + 2 * nuc_index;
}//size_t fit_parameters_index_for_source(...)


bool RelActAutoSolution::walk_to_controlling_nuclide( SrcVariant &src, const size_t rel_eff_index, double &multiple ) const
{
  assert( multiple == 1.0 );
  multiple = 1.0;
  assert( rel_eff_index < m_rel_activities.size() );
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "walk_to_controlling_nuclide: invalid rel eff index" );
  if( m_rel_activities.size() != m_options.rel_eff_curves.size() )
    throw std::logic_error( "walk_to_controlling_nuclide: size(m_rel_activities) != size(m_options.rel_eff_curves)" );

  const vector<NuclideRelAct> &rel_acts = m_rel_activities[rel_eff_index];
  
  // Check source is valid one in the rel eff curve
  bool found_src = false;
  for( size_t i = 0; !found_src && (i < m_rel_activities[rel_eff_index].size()); ++i )
    found_src = ( src == m_rel_activities[rel_eff_index][i].source );
  if( !found_src )
    throw std::logic_error( "walk_to_controlling_nuclide: src not found in rel eff curve" );

#ifndef NDEBUG
  const SrcVariant original_src = src;
#endif

  const RelActCalcAuto::RelEffCurveInput &rel_eff_curve = m_options.rel_eff_curves[rel_eff_index];
  
  if( rel_eff_curve.act_ratio_constraints.empty() )
    return false;

  SrcVariant controller_src = std::monostate();

  bool found_controller = false;
  size_t sentinel = 0; //Dont need, but just to check the logic for development

  while( controller_src != src )
  {
    sentinel += 1;
    assert( sentinel < 100 );
    if( sentinel > 1000 )
      throw logic_error( "RelActAutoSolution::activity_ratio_uncert: possible infinite loop - logic error" );

    controller_src = src;
    for( const RelEffCurveInput::ActRatioConstraint &constraint : rel_eff_curve.act_ratio_constraints )
    {
      if( constraint.constrained_source == src )
      {
        controller_src = constraint.controlling_source;
        multiple *= constraint.constrained_to_controlled_activity_ratio;
        found_controller = true;
        break;
      }
    }
  }//while( controller_index != src_index )

#ifndef NDEBUG
  const size_t orig_src_index = nuclide_index( original_src, rel_eff_index ); 
  const size_t final_src_index = nuclide_index( src, rel_eff_index ); 
  assert( fabs((multiple * m_rel_activities[rel_eff_index][final_src_index].rel_activity_uncertainty) - m_rel_activities[rel_eff_index][orig_src_index].rel_activity_uncertainty) < 1e-6 );
#endif

  return found_controller;
}//bool walk_to_controlling_nuclide( size_t &iso_index, double &multiple ) const;


double RelActAutoSolution::activity_ratio_uncertainty( SrcVariant numerator, size_t numerator_rel_eff_index,
                                                      SrcVariant denominator, size_t denominator_rel_eff_index ) const
{
  // Check if covariance matrix is available
  if( m_phys_units_cov.empty() || m_final_parameters.empty() )
    throw std::logic_error( "activity_ratio_uncertainty: covariance matrix or parameters not available" );
  
  // Validate inputs
  assert( !RelActCalcAuto::is_null(numerator) );
  assert( !RelActCalcAuto::is_null(denominator) );
  if( RelActCalcAuto::is_null(numerator) || RelActCalcAuto::is_null(denominator) )
    throw std::logic_error( "activity_ratio_uncertainty: null source provided" );
  
  assert( numerator_rel_eff_index < m_options.rel_eff_curves.size() );
  if( numerator_rel_eff_index >= m_options.rel_eff_curves.size() )
    throw std::logic_error( "activity_ratio_uncertainty: invalid numerator rel eff index" );

  assert( denominator_rel_eff_index < m_options.rel_eff_curves.size() );
  if( denominator_rel_eff_index >= m_options.rel_eff_curves.size() )
    throw std::logic_error( "activity_ratio_uncertainty: invalid denominator rel eff index" );
  
  const RelEffCurveInput &num_rel_eff_curve = m_options.rel_eff_curves[numerator_rel_eff_index];
  const RelEffCurveInput &denom_rel_eff_curve = m_options.rel_eff_curves[denominator_rel_eff_index];
  
  double num_mult = 1.0, den_mult = 1.0;
  if( !num_rel_eff_curve.act_ratio_constraints.empty() )
    walk_to_controlling_nuclide( numerator, numerator_rel_eff_index, num_mult );
  
  if( !denom_rel_eff_curve.act_ratio_constraints.empty() )
    walk_to_controlling_nuclide( denominator, denominator_rel_eff_index, den_mult );

  if( (numerator == denominator) && (numerator_rel_eff_index == denominator_rel_eff_index) )
    return 0.0;
  
  
  // TODO: implement support when mass fractions are constrained
  const auto mass_frac_constraint = []( const SrcVariant &src, const RelEffCurveInput &curve )
              -> const RelEffCurveInput::MassFractionConstraint *{
    const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(src);
    if( !nuc )
      return nullptr;
    const auto start_iter = begin(curve.mass_fraction_constraints);
    const auto end_iter = end(curve.mass_fraction_constraints);
    const auto pos = std::find_if( start_iter, end_iter, [&]( const RelEffCurveInput::MassFractionConstraint &mfc ){
      return mfc.nuclide == nuc;
    } );
    return (pos != end_iter) ? &(*pos) : nullptr;
  };


  const size_t num_activity_param_index = fit_parameters_index_for_source( numerator, numerator_rel_eff_index );
  const size_t den_activity_param_index = fit_parameters_index_for_source( denominator, denominator_rel_eff_index );


  // Check parameter indices are valid
  if( (num_activity_param_index >= m_final_parameters.size())
     || (den_activity_param_index >= m_final_parameters.size())
     || (num_activity_param_index >= m_phys_units_cov.size())
     || (den_activity_param_index >= m_phys_units_cov.size())
     || ((num_activity_param_index < m_phys_units_cov.size()) && (den_activity_param_index >= m_phys_units_cov[num_activity_param_index].size()))
     || ((den_activity_param_index < m_phys_units_cov.size()) && (num_activity_param_index >= m_phys_units_cov[den_activity_param_index].size()))
     || (num_activity_param_index >= m_parameter_scale_factors.size())
     || (den_activity_param_index >= m_parameter_scale_factors.size())
     )
  {
    throw std::logic_error( "activity_ratio_uncertainty: parameter index out of bounds" );
  }

  // Get the parameter values
  const double num_activity = rel_activity( numerator, numerator_rel_eff_index );
  const double den_activity = rel_activity( denominator, denominator_rel_eff_index );

  #ifndef NDEBUG
  {// Begin make sure we are handling `m_parameter_scale_factors` correctly
    const RelEffCurveInput::MassFractionConstraint *numerator_mfc = mass_frac_constraint( numerator, num_rel_eff_curve );
    const RelEffCurveInput::MassFractionConstraint *denominator_mfc = mass_frac_constraint( denominator, denom_rel_eff_curve );

    if( !numerator_mfc )
    {
      const double num_scale_factor = m_parameter_scale_factors[num_activity_param_index];
      const double num_activity_other_calc = m_final_parameters[num_activity_param_index] * num_scale_factor;
      const double num_actual_rel_act = m_rel_activities[numerator_rel_eff_index][nuclide_index(numerator, numerator_rel_eff_index)].rel_activity;
      assert( ((std::max)(num_activity_other_calc, num_actual_rel_act) < 1.0E-6)
              || (fabs(num_activity_other_calc - num_actual_rel_act) < 1e-6*(std::max)(num_activity_other_calc, num_actual_rel_act)) );
      assert( ((std::max)(num_activity_other_calc, num_activity) < 1.0E-6)
              || (fabs(num_activity_other_calc - num_activity) < 1e-6*(std::max)(num_activity_other_calc, num_activity)) );
    }

    if( !denominator_mfc )
    {
      const double den_scale_factor = m_parameter_scale_factors[den_activity_param_index];
      const double den_activity_other_calc = m_final_parameters[den_activity_param_index] * den_scale_factor;
      const double den_actual_rel_act = m_rel_activities[denominator_rel_eff_index][nuclide_index(denominator, denominator_rel_eff_index)].rel_activity;
      assert( ((std::max)(den_activity_other_calc, den_actual_rel_act) < 1.0E-6)
              || (fabs(den_activity_other_calc - den_actual_rel_act) < 1e-6*(std::max)(den_activity_other_calc, den_actual_rel_act)) );
      assert( ((std::max)(den_activity, den_actual_rel_act) < 1.0E-6)
              || (fabs(den_activity - den_actual_rel_act) < 1e-6*(std::max)(den_activity, den_actual_rel_act)) );
    }
  }// End make sure we are handling `m_parameter_scale_factors` correctly
  #endif
  
  // Check denominator is not zero
  if( fabs(den_activity) < std::numeric_limits<double>::epsilon() )
    throw std::logic_error( "activity_ratio_uncertainty: denominator activity is zero" );
  
  // Get covariance elements
  const double var_num = m_phys_units_cov[num_activity_param_index][num_activity_param_index];
  const double var_den = m_phys_units_cov[den_activity_param_index][den_activity_param_index];
  const double cov_num_den = m_phys_units_cov[num_activity_param_index][den_activity_param_index];

  // Check that variances are non-negative
  if( var_num < 0.0 || var_den < 0.0 )
    throw std::logic_error( "activity_ratio_uncertainty: negative variance in covariance matrix" );
  
  //const double ratio_variance = (var_num - 2.0 * num_activity * cov_num_den / den_activity
  //                              + num_activity * num_activity * var_den / (den_activity * den_activity))
  //                              / (den_activity * den_activity);

  //const double ratio_variance = (var_num
  //                                + ((num_activity*num_activity*var_den)/(den_activity*den_activity))
  //                                - ((2.0*num_activity*cov_num_den)/den_activity))
  //                              / (den_activity*den_activity);

  const double ratio_2 = std::pow( num_activity / den_activity, 2.0 );
  const double uncert_frac_2_num = (var_num/(num_activity*num_activity));
  const double uncert_frac_2_den = (var_den/(den_activity*den_activity));
  const double corr_term_2 = (2.0*cov_num_den)/(num_activity*den_activity);

  const double ratio_variance = ratio_2 * (uncert_frac_2_num + uncert_frac_2_den - corr_term_2);

  // Check for negative variance (numerical issues)
  if( ratio_variance < 0.0 )
  {
    // This could happen for a few reasons
    //  - Numerical precision issues - we are out of luck anyway
    //  - the Covariance matrix is calculated incorrectly; like if there is some spurious correlation - again out of luck.
    //  - Small denominator; the square or cube, or whatever of these valuse amplifies contributions of
    //    correlation too much.  In this case we can calc the uncertainty of the inverse ratio, and then convert
    //    this to what we want - but this is a bit sketch, and probably means the ratio in nearly meaningless anyway
    //    maybe.
    if( fabs(den_activity) < fabs(num_activity) )
    {
      const double inv_ratio_2 = std::pow( den_activity / num_activity, 2.0 );
      const double cov_den_num = m_phys_units_cov[den_activity_param_index][num_activity_param_index];
      const double inv_corr_term_2 = (2.0*cov_den_num)/(num_activity*den_activity);
      const double inv_ratio_variance = inv_ratio_2 * (uncert_frac_2_den + uncert_frac_2_num - inv_corr_term_2);
      if( inv_ratio_variance > 0.0 )
      {
        const double inv_ratio_uncert = sqrt( inv_ratio_variance );
        const double inv_ratio_val = den_activity / num_activity;
        const double orig_ratio_val = num_activity / den_activity;
        const double inv_inv_uncert = orig_ratio_val * (inv_ratio_uncert / inv_ratio_val);

        return inv_inv_uncert;
      }//if( inv_ratio_variance > 0.0 )
    }//if( den_activity < num_activity )

    return 0.0;
  }//if( ratio_variance < 0.0 )

  return std::sqrt( ratio_variance );
}//double activity_ratio_uncertainty(...)


double RelActAutoSolution::rel_activity( const SrcVariant &src, const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "rel_activity: invalid rel eff index" );

  const size_t nuc_index = nuclide_index( src, rel_eff_index );
  const vector<NuclideRelAct> &nuclides = m_rel_activities[rel_eff_index];
    
  // If Pu242, we will use the corrected mass ratio to Pu239 to get activity relative
  //  to Pu239, and return the multiple of this.
  //  The other Pu isotopes dont need a correction, I dont think
  const SandiaDecay::Nuclide * const nuc = nuclide(src);
  if( nuc && (nuc->atomicNumber == 94) && (nuc->massNumber == 242)
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

  return nuclides[nuc_index].rel_activity;
}//double rel_activity(...)
  

double RelActAutoSolution::nuclide_counts( const SrcVariant &src, const size_t rel_eff_index ) const
{
  assert( !RelActCalcAuto::is_null(src) );
  if( RelActCalcAuto::is_null(src) )
    throw std::logic_error( "RelActAutoSolution::nuclide_counts: invalid source" );
  
  double total_counts = 0.0;
  
  const double live_time = m_spectrum ? static_cast<double>(m_spectrum->live_time()) : 0.0;

  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "nuclide_counts: invalid rel eff index" ); 
  
  const vector<NuclideRelAct> &nuclides = m_rel_activities[rel_eff_index];

  for( const RoiRange &range : m_final_roi_ranges )
  {
    bool found_src = false;
    for( size_t nuc_index = 0; nuc_index < nuclides.size(); ++nuc_index )
    {
      const NuclideRelAct &nucinfo = nuclides[nuc_index];
      if( nucinfo.source != src )
        continue;

      found_src = true;

      const double rel_act = rel_activity( nucinfo.source, rel_eff_index ); //almost same as `nucinfo.rel_activity`, except with Pu corrections
        
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
                const size_t par_index = m_add_br_uncert_start_index + range_index;
                const double par_value = m_final_parameters[par_index];
                const double num_sigma_from_nominal = (par_value
                                                        - RelActCalcAutoImp::RelActAutoCostFcn::sm_peak_range_uncert_offset)
                                                          /RelActCalcAutoImp::RelActAutoCostFcn::sm_peak_range_uncert_par_scale;
                
                br_uncert_adj = max(0.0, 1.0 + num_sigma_from_nominal*m_options.additional_br_uncert);
                
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

    if( !found_src )
      throw runtime_error( "nuclide_counts: nuclide " + (is_null(src) ? string("null") : to_name(src)) + " not found in rel eff curve" );
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

  
pair<double,double> RelActAutoSolution::relative_efficiency_with_uncert( const double energy, const size_t rel_eff_index ) const
{
  typedef ceres::Jet<double,RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size> Jet;
  assert( m_rel_eff_forms.size() == m_rel_eff_coefficients.size() );
  assert( m_rel_eff_covariance.empty() || (m_rel_eff_covariance.size() == m_rel_eff_coefficients.size()) );
  
  assert( rel_eff_index < m_rel_eff_forms.size() );
  if( rel_eff_index >= m_rel_eff_forms.size() )
    throw std::logic_error( "relative_efficiency: invalid rel eff index" );
  
  if( m_covariance.empty() || (m_covariance.size() != m_final_parameters.size()) )
    throw runtime_error( "relative_efficiency_with_uncert: no valid covariance matrix." );
  
  const size_t num_par = m_final_parameters.size();
  double rel_eff = -1.0;
  vector<double> jacobian( num_par, 0.0 );
  for( size_t i = 0; i < num_par; i += RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size )
  {
    vector<Jet> x_local( begin(m_final_parameters), end(m_final_parameters) );
    for( size_t j = 0; (j < RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size) && (i+j < num_par); ++j )
      x_local[i+j].v[j] = 1.0;
    
    const Jet rel_eff_jet = m_cost_functor->relative_eff(energy, rel_eff_index, x_local);
    
    for( size_t j = 0; (j < RelActCalcAutoImp::RelActAutoCostFcn::sm_auto_diff_stride_size) && ((i+j) < num_par); ++j )
      jacobian[i+j] = rel_eff_jet.v[j];
    rel_eff = rel_eff_jet.a;
  }
  
  double uncertainty = 0.0;
  for( size_t i = 0; i < num_par; ++i )
  {
    assert( m_covariance[i].size() == num_par );
    
    if( m_covariance[i].size() != num_par )
      throw logic_error( "relative_efficiency_with_uncert: invalid Rel. Eff. covariance matrix row." );
    
    for( size_t j = 0; j < num_par; ++j )
      uncertainty += jacobian[i]*m_covariance[i][j]*jacobian[j];
  }
  
  uncertainty = std::sqrt( uncertainty );
  
#ifndef NDEBUG
  // We can double-check this rel eff calc, but more importantly, for non-physical models, we
  //  can check the uncertainty calculated with auto-differentiation (ie, above), against an
  //  analytical solution.
  const double manual_rel_eff = RelActAutoSolution::relative_efficiency( energy, rel_eff_index );
  const double releff_diff = fabs(manual_rel_eff - rel_eff);
  assert( (releff_diff < 1.0E-8) || (releff_diff < 1.0E-6*(std::max)(manual_rel_eff, rel_eff)) );
  
  if( m_rel_eff_forms[rel_eff_index] != RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    const RelActCalc::RelEffEqnForm eqn_form = m_rel_eff_forms[rel_eff_index];
    const vector<double> &coeffs = m_rel_eff_coefficients[rel_eff_index];
    assert( !m_rel_eff_covariance.empty() );
    assert( m_rel_eff_covariance.size() == m_rel_eff_coefficients.size() );
    if( rel_eff_index >= m_rel_eff_covariance.size() )
      throw runtime_error( "relative_efficiency_with_uncert: no valid Rel. Eff. covariance matrix." );
    const vector<vector<double>> &cov = m_rel_eff_covariance[rel_eff_index];
    assert( cov.empty() || (cov.size() == coeffs.size()) );
    if( cov.size() != coeffs.size() )
      throw logic_error( "relative_efficiency_with_uncert: invalid Rel. Eff. covariance matrix." );
    
    const double manual_uncert = RelActCalc::eval_eqn_uncertainty( energy, eqn_form, coeffs, cov );
    const double diff = fabs(manual_uncert - uncertainty);
    assert( (diff < 1.0E-8) || (diff < 1.0E-6*(std::max)(manual_uncert, uncertainty)) );
  }//if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
#endif //NDEBUG
  
  if( IsNan(rel_eff) || IsInf(rel_eff) )
    throw runtime_error( "relative_efficiency_with_uncert: invalid rel eff value." );
  
  if( IsNan(uncertainty) || IsInf(uncertainty) )
    throw runtime_error( "relative_efficiency_with_uncert: invalid uncertainty." );
  
  return {rel_eff, uncertainty};
}//pair<double,double> relative_efficiency_with_uncert( const double energy, const size_t rel_eff_index ) const

  
  
size_t RelActAutoSolution::nuclide_index( const SrcVariant &src, const size_t rel_eff_index ) const
{
  assert( m_rel_activities.size() == m_options.rel_eff_curves.size() );
  assert( !RelActCalcAuto::is_null(src) );
  if( RelActCalcAuto::is_null(src) )
    throw std::logic_error( "RelActAutoSolution::nuclide_index: invalid source" );
  
  assert( rel_eff_index < m_rel_activities.size() );
  if( rel_eff_index >= m_rel_activities.size() )
    throw std::logic_error( "nuclide_index: invalid rel eff index" );
  
  for( size_t i = 0; i < m_rel_activities[rel_eff_index].size(); ++i )
  {
    if( src == m_rel_activities[rel_eff_index][i].source )
    {
      return i;
    }
  }//for( size_t i = 0; i < m_rel_activities[rel_eff_index].size(); ++i )
  
  assert( 0 );
  throw runtime_error( "RelActAutoSolution: " + to_name(src) + " is an invalid nuclide." );
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
  
  const RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef input
                    = RelActCalcAutoImp::RelActAutoCostFcn::make_phys_eqn_input(
                            m_options.rel_eff_curves[rel_eff_index], m_drf, coeffs, 0 );
  
    
  return RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                            input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
}//string rel_eff_eqn_js_function() const
  
  
string RelActAutoSolution::rel_eff_eqn_js_uncert_fcn( const size_t rel_eff_index ) const
{ 
  if( m_final_roi_ranges.empty() )
    return "null";
  
  vector<RelActCalcAuto::RoiRange> rois = m_final_roi_ranges;
  
  std::sort( begin(rois), end(rois), []( const RelActCalcAuto::RoiRange &lhs, const RelActCalcAuto::RoiRange &rhs ){
    return lhs.lower_energy < rhs.lower_energy;
  } );

  vector<PeakDef> fit_peaks = m_fit_peaks;
  std::sort( begin(fit_peaks), end(fit_peaks), []( const PeakDef &lhs, const PeakDef &rhs ){
    return lhs.mean() < rhs.mean();
  } );
  
  double lower_energy = std::min( rois.front().lower_energy, fit_peaks.empty() ? 3000.0 : fit_peaks.front().mean() );
  double upper_energy = std::max( rois.back().upper_energy, fit_peaks.empty() ? 3000.0 : fit_peaks.back().mean() );
  
  if( lower_energy > 100 )
    lower_energy -=15;
  else if( lower_energy > 10 )
    lower_energy -= 5;
  else if( lower_energy > 1 )
    lower_energy -= 1;

  vector<double> energies;
  double current_energy = lower_energy;
  for( const PeakDef &peak : m_fit_peaks )
  {
    upper_energy = std::max( upper_energy, peak.mean() );

    double min_dx = 1.0;
    if( current_energy < 130 )
      min_dx = 1.0;
    else if( current_energy < 300 )
      min_dx = 5.0;
    else
      min_dx = 15.0;

    // We'll try to get in at least ~10 points between each peak
    if( peak.mean() > current_energy )
      min_dx = std::min( min_dx, 0.1*(peak.mean() - current_energy) );
    min_dx = std::max( min_dx, 1.0 ); //but less than a keV between points is just to small.
    
    for( ; current_energy < peak.mean(); current_energy += min_dx )
      energies.push_back( current_energy );
    
    if( !energies.empty() && (energies.back() < peak.mean()) )
      current_energy = peak.mean();
  }//for( const PeakDef &peak : m_fit_peaks )
  
  for( ; current_energy < upper_energy; current_energy += 15 )
    energies.push_back( current_energy );
  
  size_t num_points = 0;
  string fcn = "function(x){\n"
  "  const points = [";
  bool is_first_point = true;
  
  //cout << "RelEff: [";
  for( double x : energies )
  {
    try
    {
      const std::pair<double,double> rel_eff_uncert = relative_efficiency_with_uncert( x, rel_eff_index );
      //const double y = rel_eff_uncert.first;
      const double unc = rel_eff_uncert.second;
      //cout << "{" << rel_eff_uncert.first << ", " << unc << "}, ";

      //assert( (y >= 0.0) && !IsNan(y) && !IsInf(y) ); //Can happen when we are out of bounds of the physical model
      if( isnan(unc) || isinf(unc) )
        continue;
      
      fcn += is_first_point ? "" : ",";
      fcn += "[" + SpecUtils::printCompact(x, 4) + "," + SpecUtils::printCompact(unc, 4) + "]";
      
      num_points += 1;
      is_first_point = false;
    }catch( std::exception &e )
    {
      // This can happen when we are out of bounds of the physical model
    }
  }//for( double x : energies )
  
  //cout << "]" << endl;
  
  if( num_points < 2 )
    return "null";
  
  fcn += "];\n"
  "  if( x <= points[0][0] )\n"
  "    return points[0][1];\n"
  "  if( x >= points[points.length - 1][0] )\n"
  "    return points[points.length - 1][1];\n"
  "  for (let i = 0; i < points.length - 1; i++) {\n"
  "    const [x1, y1] = points[i];\n"
  "    const [x2, y2] = points[i + 1];\n"
  "    if( x >= x1 && x <= x2) {\n"
  "      const t = (x - x1) / (x2 - x1);\n"
  "      return y1 + t * (y2 - y1);\n"
  "    }\n"
  "  }\n"
  "console.assert(0,'Shouldnt get here in interpolating');\n"
  "return points[points.length - 1][1];"
  "}";
  
  return fcn;
}//string RelActAutoSolution::rel_eff_eqn_js_uncert_fcn(i)
  
  
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
  

  const RelActAutoSolution orig_sol = RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres(
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
      const double fwhm = eval_fwhm( energy, current_sol.m_fwhm_form, current_sol.m_fwhm_coefficients.data(), 
                                    current_sol.m_fwhm_coefficients.size(), current_sol.m_drf );
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
    
    
    vector<optional<RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef<double>>> phys_model_inputs( current_sol.m_rel_eff_coefficients.size() );
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
        
        optional<RelActCalcAutoImp::RelActAutoCostFcn::PhysModelRelEqnDef<double>> &phys_model_input = phys_model_inputs[rel_eff_index];
        if( (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel) && !phys_model_input.has_value() )
          phys_model_input = RelActCalcAutoImp::RelActAutoCostFcn::make_phys_eqn_input( rel_eff, current_sol.m_drf, rel_eff_coefs, 0 );
        
        for( const NuclideRelAct &rel_act : rel_acts )
        {
          assert( !is_null(rel_act.source) );
          if( !is_null(rel_act.source) )
            continue;

          const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(rel_act.source);
          
          if( (rel_eff.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable)
             && nuclide
             && (nuclide->atomicNumber == 94)
             && (nuclide->massNumber == 242) )
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
                                          current_sol.m_fwhm_coefficients.size(), 
                                          current_sol.m_drf );
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
            
            cout << "For " << energy << " keV " << to_name(rel_act.source)
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
      = RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres( updated_options, foreground, background, input_drf, all_peaks, cancel_calc );
      
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

  
string NuclideRelAct::name() const
{
  assert( !is_null(source) );
  if( is_null(source) )
    return "";
  return to_name(source);
}//NuclideRelAct::name()
  
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
  if( lhs.source != rhs.source )
    throw std::runtime_error( "Source in lhs and rhs are not the same" );
  
  if( fabs(lhs.age - rhs.age) > 1.0e-5 * std::max(lhs.age, rhs.age) )
    throw std::runtime_error( "Age in lhs and rhs are not the same" );
  
  if( lhs.fit_age != rhs.fit_age )
    throw std::runtime_error( "Fit age in lhs and rhs are not the same" );
  
  if( lhs.gammas_to_exclude.size() != rhs.gammas_to_exclude.size() )
    throw std::runtime_error( "Number of gammas to exclude in lhs and rhs are not the same" );
  
  if( lhs.fit_age_min.has_value() != rhs.fit_age_min.has_value() )
    throw std::runtime_error( "Fit age min in lhs and rhs are not the same" );
  if( lhs.fit_age_min.has_value() && rhs.fit_age_min.has_value() )
    if( fabs(lhs.fit_age_min.value() - rhs.fit_age_min.value()) > 1.0e-5 )
      throw std::runtime_error( "Fit age min in lhs and rhs are not the same" );

  if( lhs.fit_age_max.has_value() != rhs.fit_age_max.has_value() )
    throw std::runtime_error( "Fit age max in lhs and rhs are not the same" );
  if( lhs.fit_age_max.has_value() && rhs.fit_age_max.has_value() )
    if( fabs(lhs.fit_age_max.value() - rhs.fit_age_max.value()) > 1.0e-5 )
      throw std::runtime_error( "Fit age max in lhs and rhs are not the same" );

  if( lhs.min_rel_act.has_value() != rhs.min_rel_act.has_value() )
    throw std::runtime_error( "Min rel act in lhs and rhs are not the same" );
  if( lhs.min_rel_act.has_value() && rhs.min_rel_act.has_value() )
    if( fabs(lhs.min_rel_act.value() - rhs.min_rel_act.value()) > 1.0e-5 )
      throw std::runtime_error( "Min rel act in lhs and rhs are not the same" );

  if( lhs.max_rel_act.has_value() != rhs.max_rel_act.has_value() )
    throw std::runtime_error( "Max rel act in lhs and rhs are not the same" );
  if( lhs.max_rel_act.has_value() && rhs.max_rel_act.has_value() )
    if( fabs(lhs.max_rel_act.value() - rhs.max_rel_act.value()) > 1.0e-5 )
      throw std::runtime_error( "Max rel act in lhs and rhs are not the same" );

  if( lhs.starting_rel_act.has_value() != rhs.starting_rel_act.has_value() )
    throw std::runtime_error( "Starting rel act in lhs and rhs are not the same" );
  if( lhs.starting_rel_act.has_value() && rhs.starting_rel_act.has_value() )
    if( fabs(lhs.starting_rel_act.value() - rhs.starting_rel_act.value()) > 1.0e-5 )
      throw std::runtime_error( "Starting rel act in lhs and rhs are not the same" );
      

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

  if( lhs.act_ratio_constraints.size() != rhs.act_ratio_constraints.size() )
    throw std::runtime_error( "Number of nuclide constraints in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.act_ratio_constraints.size(); ++i )
    ActRatioConstraint::equalEnough( lhs.act_ratio_constraints[i], rhs.act_ratio_constraints[i] );

  if( lhs.mass_fraction_constraints.size() != rhs.mass_fraction_constraints.size() )
    throw std::runtime_error( "Number of mass fraction constraints in lhs and rhs are not the same" );
  
  for( size_t i = 0; i < lhs.mass_fraction_constraints.size(); ++i )
    MassFractionConstraint::equalEnough( lhs.mass_fraction_constraints[i], rhs.mass_fraction_constraints[i] );
}//RelEffCurveInput::equalEnough


void Options::equalEnough( const Options &lhs, const Options &rhs )
{
  if( lhs.fit_energy_cal != rhs.fit_energy_cal )
    throw std::runtime_error( "Fit energy calibration in lhs and rhs are not the same" );
  
  if( lhs.fwhm_form != rhs.fwhm_form )
    throw std::runtime_error( "FWHM form in lhs and rhs are not the same" );
  
  if( lhs.fwhm_estimation_method != rhs.fwhm_estimation_method )
    throw std::runtime_error( "FWHM estimation method in lhs and rhs are not the same" );

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



