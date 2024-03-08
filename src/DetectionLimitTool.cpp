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

#include <mutex>
#include <atomic>
#include <iostream>

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WTableCell>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>

#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Object>
#include <Wt/Json/Serializer>

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MinosError.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FumiliMinimizer.h"
#include "Minuit2/ScanMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMinimize.h"


#include "SpecUtils/SpecFile.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/D3SpectrumExport.h"


#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DetectionLimitTool.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;
using namespace Wt;

#define MU_CHARACTER "\xCE\xBC"

namespace
{

bool use_curie_units()
{
  InterSpec *interspec = InterSpec::instance();
  if( !interspec )
    return true;
  
  return !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", interspec );
}//bool use_curie_units()

  
string det_eff_geom_type_postfix( DetectorPeakResponse::EffGeometryType type )
{
  switch( type )
  {
    case DetectorPeakResponse::EffGeometryType::FarField:
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      return "";
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      return "/cm2";
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      return "/m2";
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      return "/g";
  }//switch( m_det_type )
  assert( 0 );
  return "";
}//string det_eff_geom_type_postfix( DetectorPeakResponse::EffGeometryType )
}//namespace

// We cant have MdaPeakRow in a private namespace since we forward declare it in the header.
class MdaPeakRow : public WContainerWidget
{
protected:
  DetectionLimitTool::MdaPeakRowInput m_input;
  
  const SandiaDecay::Nuclide *const m_nuclide;
  
  WText *m_title;
  WText *m_title_info;
  
  SwitchCheckbox *m_use_for_likelihood;
  
  /** The energy, in keV that the ROI starts at.
   When the user enters a value, it will be rounded to cover a whole channel; however, even if this wasnt done, I *think* things would
   work out okay (e.g., uncertainties should be calculated correctly, I think).
   */
  NativeFloatSpinBox *m_roi_start;
  
  /** Upper range of the ROI - again rounded when users change the value. */
  NativeFloatSpinBox *m_roi_end;
  
  /** These are the number of channels outside the ROI to use to compute the ROI line for Curie-style limits.  */
  NativeFloatSpinBox *m_num_side_channels;
  
  
  WComboBox *m_continuum;
  WText *m_poisonLimit;
  
  /** The upper-limits of activity.
   
   If computing for activity limits, then will be at the specified distance.
   If computing for distance limits, will be at 1-meter.
   */
  double m_simple_mda;
  double m_simple_excess_counts;
  
  /** The maximum detection distance - only filled out if computing for distance limits.  */
  double m_simple_max_det_dist;
  
  /** Wether of not to fix the continuum for the likelihood based limit to the bins on either side of the ROI, the full ROI, or let float
   and change dependent on the activity. */
  Wt::WComboBox *m_decon_cont_norm_method;
  
  Wt::Signal<> m_changed;
  
  
  void handleUseForLikelihoodChanged()
  {
    const bool use = m_use_for_likelihood->isChecked();
    m_decon_cont_norm_method->setEnabled( use );
    
    const bool fixed_at_edges = (m_decon_cont_norm_method->currentIndex() == static_cast<int>(DetectionLimitCalc::DeconContinuumNorm::FixedByEdges));
    m_continuum->setEnabled( use && !fixed_at_edges );
    if( fixed_at_edges )
    {
      assert( m_continuum->currentIndex() == 0 );
      assert( m_input.decon_continuum_type == PeakContinuum::OffsetType::Linear );
      m_continuum->setCurrentIndex( 0 );
      m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
    }
    
    m_input.use_for_likelihood = use;
    emitChanged();
  }

  
  void setSimplePoisonTxt()
  {
    try
    {
      const bool fixed_geom = m_input.drf->isFixedGeometry();
      const DetectorPeakResponse::EffGeometryType det_geom = m_input.drf->geometryType();
      const float &energy = m_input.energy;
      const double &distance = m_input.distance;
      const bool useCuries = use_curie_units();
      const double det_eff = fixed_geom ? m_input.drf->intrinsicEfficiency(energy)
                                        : m_input.drf->efficiency(energy, distance);
      const double counts_4pi = ((m_input.do_air_attenuation && !fixed_geom)
                                 ? m_input.counts_per_bq_into_4pi_with_air
                                 : m_input.counts_per_bq_into_4pi);
      const double gammas_per_bq = counts_4pi * det_eff;
      
      const DetectionLimitCalc::CurieMdaInput input = currieInput();
      
      const DetectionLimitCalc::CurieMdaResult result = DetectionLimitCalc::currie_mda_calc( input );
      //print_summary( cout, result, -1.0f );
      
      m_simple_excess_counts = result.source_counts;
      m_simple_mda = result.upper_limit / gammas_per_bq;
      m_simple_max_det_dist = -1;
      
      switch( m_input.limit_type )
      {
        case DetectionLimitTool::LimitType::Activity:
        {
          if( result.source_counts > result.decision_threshold )
          {
            // There is enough excess counts that we would reliably detect this activity, so we will
            //  give the activity range.
            const float lower_act = result.lower_limit / gammas_per_bq;
            const float upper_act = result.upper_limit / gammas_per_bq;
            const float nominal_act = result.source_counts / gammas_per_bq;
            
            const string lowerstr = PhysicalUnits::printToBestActivityUnits( lower_act, 2, useCuries )
                                    + det_eff_geom_type_postfix( det_geom );
            const string upperstr = PhysicalUnits::printToBestActivityUnits( upper_act, 2, useCuries )
                                    + det_eff_geom_type_postfix( det_geom );
            const string nomstr = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, useCuries )
                                  + det_eff_geom_type_postfix( det_geom );
            
            //cout << "At " << m_input.energy << "keV observed " << result.source_counts
            //     << " counts, at distance " << PhysicalUnits::printToBestLengthUnits(m_input.distance)
            //     << ", leading to nominal activity " << nomstr
            //     << endl;
            
            m_poisonLimit->setText( "<div>Observed " + nomstr + "</div>"
                                    "<div>Range [" + lowerstr + ", " + upperstr + "]</div>" );
            m_poisonLimit->setToolTip( "Detected activity, using just this Region Of Interests,"
                                      " and the observed excess of counts, as well as the"
                                      " statistical confidence interval." );
          }else if( result.upper_limit < 0 )
          {
            // This can happen when there are a lot fewer counts in the peak region than predicted
            //  from the sides - since this is non-sensical, we'll just say zero.
            const string unitstr = (useCuries ? "Ci" : "Bq") + det_eff_geom_type_postfix( det_geom );
            m_poisonLimit->setText( "<div>MDA: &lt; 0" + unitstr + "</div><div>(sig. fewer counts in ROI than predicted)</div>" );
            m_poisonLimit->setToolTip( "Significantly fewer counts were observed in the"
                                       " Region Of Interest, than predicted by the neighboring channels." );
          }else
          {
            // We will provide the upper bound on activity.
            const string mdastr = PhysicalUnits::printToBestActivityUnits( m_simple_mda, 2, useCuries )
                                  + det_eff_geom_type_postfix( det_geom );
            
            m_poisonLimit->setText( "MDA: " + mdastr );
            m_poisonLimit->setToolTip( "Minimum Detectable Activity, using just this Region Of"
                                      " Interests, assuming no signal is present.\n"
                                      "Basically uses the Currie method of calculating." );
          }
          
          break;
        }//case DetectionLimitTool::LimitType::Activity:
          
        case DetectionLimitTool::LimitType::Distance:
        {
          // We will iterate to find the distance in the below.
          //  If we arent taking into account attenuation in the air, its easy to solve for the
          //  distance.
          //  But, even if taking into account the air attenuation, the distance is still solvable
          //  in something like Mathematica, but the equation is like
          //   x^{-2}*exp(-0.000942*x) = A_{lim}  (where exp(-0.000942*x) is attenuation in air
          //    for 661 keV, and A_{lim} is for at one meter, and x is meters)
          //  so, for the moment its easier to just be consistent and use an iterative approach
          //  always.
          shared_ptr<const DetectorPeakResponse> drf = m_input.drf;
          if( !drf || !drf->isValid() )
            throw runtime_error( "DRF invalid" );
          
          // Round detection probability to nearest 1 decimal place.
          const double rnd_cl_percent = std::round(1000.0*input.detection_probability) / 10.0;
          
          const double activity = m_input.activity;
          if( activity <= 0.0 )
            throw runtime_error( "No activity specified" );
          
          // Make a convenience lambda that returns the efficiency taking into account both the
          //  geometric factor, and attenuation in the air.
          const double intrinsic_eff = m_input.drf->intrinsicEfficiency(m_input.energy);
          
          
          auto counts_at_distance = [this,intrinsic_eff,activity]( const double dist ) -> double {
            const double counts_4pi_no_air = m_input.counts_per_bq_into_4pi;
            const double geom_eff = m_input.drf->fractionalSolidAngle(m_input.drf->detectorDiameter(), dist);
            double air_eff = 1.0;
            if( m_input.do_air_attenuation )
            {
              const double coef = GammaInteractionCalc::transmission_coefficient_air( m_input.energy, dist );
              air_eff = exp( -1.0 * coef );
            }
            
            return activity * counts_4pi_no_air * geom_eff * intrinsic_eff * air_eff;
          };//counts_at_distance(...)
          
          if( result.source_counts > counts_at_distance(0.0) )
          {
            cout << "For " << input.gamma_energy << ", observe " << result.source_counts
                 << ", but expect " << counts_at_distance(0.0) << " at R=0" << endl;
            throw runtime_error( "More counts observed than would be on contact." );
          }
          
          if( result.source_counts > result.decision_threshold )
          {
            // There is enough excess counts that we would reliably detect this activity, at the
            //  measured distance, so we will calculate the range of distances the measurement
            //  may have happened at (e.x., when you take a measurement, and dont know how far
            //  away the source is; e.g., source behind wall).
            
            
            // We need to find the distance corresponding to the amount of counts given by both
            //  result.lower_limit and result.upper_limit
            // We will provide: "To the 95% CL, the measurement was taken between X and Y meters"
            
            if( result.lower_limit < 0.0 )
              throw runtime_error( "Lower limit of counts is less than zero, but counts observed is greater than L_d." );
            
            double max_distance = 1.0*PhysicalUnits::meter;
            
            while( max_distance < 100000.0*PhysicalUnits::meter ) //100 km is getting towards limits of numerical accuracy for DRF efficiency
            {
              const double n_expected_counts = counts_at_distance(max_distance);
              if( n_expected_counts < result.lower_limit )
                break;
              
              max_distance *= 2.0;
            }
            
            if( counts_at_distance(max_distance) >= result.upper_limit )
              throw runtime_error( "Maximum distance to large (>100km)" );
            
            auto lower_limit_dist = [&]( double dist ) -> double {
              return fabs( counts_at_distance(dist) - result.lower_limit );
            };
            
            auto upper_limit_dist = [&]( double dist ) -> double {
              return fabs( counts_at_distance(dist) - result.upper_limit );
            };
            
            auto nominal_dist = [&]( double dist ) -> double {
              return fabs( counts_at_distance(dist) - result.source_counts );
            };
            
            using boost::math::tools::brent_find_minima;
            
            /* brent_find_minima will use the following tolerances for the given bits:
             
             bits | Tolerance   | Remark
             _____|_____________|_________
             8    | 0.0078125   |
             10   | 0.00195312  |
             12   | 0.000488281 | Max for float (e.g., half std::numeric_limits<float>::digits)
             16   | 3.05176e-05 |
             20   | 1.90735e-06 |
             26   | 2.98023e-08 | Max for double (e.g., half std::numeric_limits<double>::digits)
             
             There might be a factor of 1/2 or something in brent_find_minima, I didnt follow its
             code super close.
             
             With 12 bits it looks like we are still almost always below 20 iterations.
             */
            const int bits = 12;
            
            boost::uintmax_t max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
            const pair<double, double> lower_crossing
                         = brent_find_minima( lower_limit_dist, 0.0, max_distance, bits, max_iter );
            
            max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
            const pair<double, double> upper_crossing
                         = brent_find_minima( upper_limit_dist, 0.0, max_distance, bits, max_iter );
            
            max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
            const pair<double, double> nominal_crossing
                             = brent_find_minima( nominal_dist, 0.0, max_distance, bits, max_iter );
            
            //cout << "Energy " << m_input.energy << "keV: nominal_dist={"
            //     << PhysicalUnits::printToBestLengthUnits(nominal_crossing.first)
            //     << ", " << nominal_crossing.second << "}, searched between 0 and "
            //     << PhysicalUnits::printToBestLengthUnits(max_distance)
            //     << " with " << max_iter << " iterations, and source activity "
            //     << PhysicalUnits::printToBestActivityUnits(activity, 2, useCuries)
            //     << "; source counts were " << result.source_counts
            //<< endl;
            
            const double lower_distance = lower_crossing.first;
            const double upper_distance = upper_crossing.first;
            const double nominal_distance = nominal_crossing.first;
            
            m_simple_max_det_dist = upper_distance;
            
            const string lowerstr = PhysicalUnits::printToBestLengthUnits( lower_distance, 2 );
            const string upperstr = PhysicalUnits::printToBestLengthUnits( upper_distance, 2 );
            const string nomstr = PhysicalUnits::printToBestLengthUnits( nominal_distance, 2 );
            
            // TODO: double check that this double-sided limit is correct, and we dont need to
            //  convert to a single sided or something
            
            char buffer[256];
            snprintf( buffer, sizeof(buffer),
                      "<div>Nominal distance %s</div><div>%.1f%% CL range [%s, %s]</div>",
                      nomstr.c_str(), rnd_cl_percent, upperstr.c_str(), lowerstr.c_str() );
            
            m_poisonLimit->setText( buffer );
          }else if( result.upper_limit < 0 )
          {
            // This can happen when there are a lot fewer counts in the peak region than predicted
            //  from the sides - since this is non-sensical ...
            m_poisonLimit->setText( "TODO: handle case where observed deficit of counts is really large" );
          }else
          {
            // We will provide a "At the 95% CL you would detect this source at X meters" answer
#ifndef _WIN32
  #warning "go back through this logic when I am less tired to make sure result.upper_limit is really what we want to use here"
#endif
            // TODO: go back through this logic when I am less tired to make sure result.upper_limit is really what we want to use here
            
            if( result.upper_limit < 0.0 )
              throw runtime_error( "Upper limit of counts is less than zero, but counts observed is greater than L_d." );
            
            double max_distance = 1.0*PhysicalUnits::meter;
            
            while( max_distance < 100000.0*PhysicalUnits::meter )
            {
              if( counts_at_distance(max_distance) < result.upper_limit )
                break;
              max_distance *= 2.0;
            }
            
            if( counts_at_distance(max_distance) >= result.upper_limit )
              throw runtime_error( "Maximum distance to large (>100km)" );
            
            auto upper_limit_dist = [&]( double dist ) -> double {
              return fabs( counts_at_distance(dist) - result.upper_limit );
            };
            
            using boost::math::tools::brent_find_minima;
            const int bits = 12; //Float has 24 bits of mantisa, so 12 would be its max precision
            boost::uintmax_t max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
            const pair<double, double> upper_crossing
                         = brent_find_minima( upper_limit_dist, 0.0, max_distance, bits, max_iter );
            
            const double upper_distance = upper_crossing.first;
            m_simple_max_det_dist = upper_distance;
            
            const string upperstr = PhysicalUnits::printToBestLengthUnits( upper_distance, 2 );
            
            char buffer[256];
            snprintf( buffer, sizeof(buffer),
                     "At %.1f%% CL you would detect source at %s",
                     rnd_cl_percent, upperstr.c_str() );
            
            m_poisonLimit->setText( buffer );
          }
          
          break;
        }//case DetectionLimitTool::LimitType::Distance:
      }//switch( m_input.limit_type )
      
      
      
    }catch( std::exception &e )
    {
      m_poisonLimit->setText( "Error calculating limit: " + string(e.what()) );
    }
  }//void setSimplePoisonTxt()
  
  
  void emitChanged()
  {
    setSimplePoisonTxt();
    m_changed.emit();
  }
  
  void handleContinuumNormMethodChange()
  {
    auto norm_type = static_cast<DetectionLimitCalc::DeconContinuumNorm>( m_decon_cont_norm_method->currentIndex() );
    switch( norm_type )
    {
      case DetectionLimitCalc::DeconContinuumNorm::Floating:
      case DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange:
        m_continuum->setEnabled( true );
        break;
        
      case DetectionLimitCalc::DeconContinuumNorm::FixedByEdges:
        m_continuum->setCurrentIndex( 0 );
        m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
        m_continuum->setEnabled( false );
        break;
        
      default:
        assert(0);
        norm_type = DetectionLimitCalc::DeconContinuumNorm::Floating;
        m_decon_cont_norm_method->setCurrentIndex( static_cast<int>(DetectionLimitCalc::DeconContinuumNorm::Floating ) );
        break;
    }//switch( m_decon_cont_norm_method->currentIndex() )
    
    m_input.decon_cont_norm_method = norm_type;
    
    emitChanged();
  }//void handleContinuumNormMethodChange()
  
  
  void roiChanged()
  {
    // We'll validate/check user quantities here, even though this should have already all been done
    setRoiStart( m_roi_start->value() );
    setRoiEnd( m_roi_end->value() );
    
    float numSideChannel = m_num_side_channels->value();
    if( numSideChannel < 1.0f )
    {
      m_num_side_channels->setValue(1.0f);
      numSideChannel = 1.0;
    }else if( std::round(numSideChannel) != numSideChannel )
    {
      numSideChannel = std::round(numSideChannel);
      m_num_side_channels->setValue(numSideChannel);
    }
    
    assert( numSideChannel > 0.0 );
    m_input.num_side_channels = static_cast<size_t>( std::round(numSideChannel) );
    
    switch( m_continuum->currentIndex() )
    {
      case 0:
        m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
        break;
      case 1:
        m_input.decon_continuum_type = PeakContinuum::OffsetType::Quadratic;
        break;
      default:
        assert( 0 );
        m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
        m_continuum->setCurrentIndex( 0 );
        break;
    }//switch( m_continuum->currentIndex() )
    
    auto norm_method = static_cast<DetectionLimitCalc::DeconContinuumNorm>( m_decon_cont_norm_method->currentIndex() );
    switch( norm_method )
    {
      case DetectionLimitCalc::DeconContinuumNorm::Floating:
      case DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange:
        break;
        
      case DetectionLimitCalc::DeconContinuumNorm::FixedByEdges:
        m_continuum->setCurrentIndex( 0 );
        m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
        break;
        
      default:
        m_input.decon_cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::Floating;
        m_continuum->setCurrentIndex( 0 );
        m_decon_cont_norm_method->setCurrentIndex( 0 );
        m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
        break;
    }//switch( m_decon_cont_norm_method->currentIndex() )
    
    m_input.decon_cont_norm_method = norm_method;
    
    emitChanged();
  }//void roiChanged()
  
  
  float roundToNearestChannelEdge( const float energy ) const
  {
    if( !m_input.measurement )
      return energy;
    
    auto cal = m_input.measurement->energy_calibration();
    
    if( !cal || !cal->valid() )
      return energy;
    
    const double channel = std::max( 0.0, cal->channel_for_energy(energy) ); //std::max isnt necessary, but JIC
    const double whole = std::floor(channel);
    const double frac = channel - whole;
    
    if( (frac >= 0.5) && ((whole+1) < cal->num_channels()) )
      return static_cast<float>( cal->energy_for_channel(whole + 1.0) );
    
    return static_cast<float>( cal->energy_for_channel(whole) );
  }//float roundToNearestChannelEdge( float energy )
  
  
public:
  MdaPeakRow( const DetectionLimitTool::MdaPeakRowInput &input,
             const SandiaDecay::Nuclide *nuclide,
             WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
  m_input( input ),
  m_nuclide( nuclide ),
  m_title( nullptr ),
  m_title_info( nullptr ),
  m_use_for_likelihood( nullptr ),
  m_roi_start( nullptr ),
  m_roi_end( nullptr ),
  m_num_side_channels( nullptr ),
  m_continuum( nullptr ),
  m_poisonLimit( nullptr ),
  m_simple_mda( -1.0 ),
  m_simple_excess_counts( -1.0 ),
  m_simple_max_det_dist( -1.0 ),
  m_decon_cont_norm_method( nullptr ),
  m_changed( this )
  {
    addStyleClass( "MdaPeakRow" );
    
    assert( input.measurement );
    assert( input.drf && input.drf->isValid() && input.drf->hasResolutionInfo() );
    assert( input.confidence_level > 0.5 && input.confidence_level < 1.0 );
    assert( input.roi_start < input.roi_end );
    assert( input.decon_continuum_type == PeakContinuum::OffsetType::Linear
           || input.decon_continuum_type == PeakContinuum::OffsetType::Quadratic );
    
    const bool fixed_geom = input.drf->isFixedGeometry();
    const DetectorPeakResponse::EffGeometryType det_geom = input.drf->geometryType();
    
    const bool do_air_atten = (input.do_air_attenuation && !fixed_geom);
    const double fwhm = input.drf->peakResolutionFWHM( input.energy );
    const double det_eff = fixed_geom ? input.drf->intrinsicEfficiency(input.energy)
                                      : input.drf->efficiency( input.energy, input.distance );
    const double counts_4pi = (do_air_atten ? input.counts_per_bq_into_4pi_with_air
                                            : input.counts_per_bq_into_4pi);
    
    char buffer[64] = { '\0' };
    
    WContainerWidget *leftColumn = new WContainerWidget( this );
    leftColumn->addStyleClass( "MdaRowSettings" );
    
    WContainerWidget *rightColumn = new WContainerWidget( this );
    rightColumn->addStyleClass( "MdaRowCurrieMda" );
    
    const string br_str = SpecUtils::printCompact(input.branch_ratio, 3);
    snprintf( buffer, sizeof(buffer), "%.2f keV, br=%s", input.energy, br_str.c_str() );
    m_title = new WText( buffer, leftColumn );
    m_title->addStyleClass( "MdaRowTitle GridFirstCol GridFirstRow GridSpanTwoCol" );
    
    
    const string fwhm_str = SpecUtils::printCompact(fwhm,3);
    const string cnts_per_bq_str = SpecUtils::printCompact(counts_4pi*det_eff, 3);
    const string act_postfix = det_eff_geom_type_postfix( det_geom );
    
    
    snprintf( buffer, sizeof(buffer), "FWHM=%s, %s cnts/bq%s",
             fwhm_str.c_str(), cnts_per_bq_str.c_str(), act_postfix.c_str() );
    m_title_info = new WText( buffer, rightColumn );
    m_title_info->addStyleClass( "MdaRowDetInfo" );
    
    
    m_use_for_likelihood = new SwitchCheckbox( "Use for multi-peak", leftColumn );
    m_use_for_likelihood->setChecked( input.use_for_likelihood );
    m_use_for_likelihood->addStyleClass( "GridFirstCol GridSecondRow GridSpanTwoCol UseForLiklihood" );
    
    
    WLabel *label = new WLabel( "ROI Lower:", leftColumn );
    label->addStyleClass( "GridFirstCol GridThirdRow" );
    m_roi_start = new NativeFloatSpinBox( leftColumn );
    m_roi_start->addStyleClass( "GridSecondCol GridThirdRow MdaRoiInput" );
    m_roi_start->setSpinnerHidden( true );
    m_roi_start->setFormatString( "%.2f" );
    label->setBuddy( m_roi_start );
    
    label = new WLabel( "ROI Upper:", leftColumn );
    label->addStyleClass( "GridFirstCol GridFourthRow" );
    m_roi_end = new NativeFloatSpinBox( leftColumn );
    m_roi_end->addStyleClass( "GridSecondCol GridFourthRow MdaRoiInput" );
    m_roi_end->setSpinnerHidden( true );
    m_roi_end->setFormatString( "%.2f" );
    label->setBuddy( m_roi_end );
    
    label = new WLabel( "Continuum", leftColumn );
    label->addStyleClass( "GridFirstCol GridFifthRow" );
    m_continuum = new WComboBox( leftColumn );
    m_continuum->addStyleClass( "GridSecondCol GridFifthRow" );
    m_continuum->addItem( "Linear" );
    m_continuum->addItem( "Quadratic" );
    m_continuum->setCurrentIndex( 0 );
    label->setBuddy( m_continuum );
    
    
    label = new WLabel( "Cont. Norm", leftColumn );
    label->addStyleClass( "GridFirstCol GridSixthRow" );
    m_decon_cont_norm_method = new WComboBox( leftColumn );
    m_decon_cont_norm_method->addStyleClass( "GridSecondCol GridSixthRow" );
    
    static_assert( static_cast<int>(DetectionLimitCalc::DeconContinuumNorm::Floating) == 0, "DeconContinuumNorm out of date" );
    static_assert( static_cast<int>(DetectionLimitCalc::DeconContinuumNorm::FixedByEdges) == 1, "DeconContinuumNorm out of date" );
    static_assert( static_cast<int>(DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange) == 2, "DeconContinuumNorm out of date" );
    m_decon_cont_norm_method->addItem( "Floating" );
    m_decon_cont_norm_method->addItem( "Fixed at edges" );
    m_decon_cont_norm_method->addItem( "Fixed full ROI" );
    m_decon_cont_norm_method->setCurrentIndex( static_cast<int>(m_input.decon_cont_norm_method) );
    
    const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
    const char *tooltip = "How the continuum normalization should be determined:"
    "<ul><li><b>Floating</b>: The polynomial continuum is fit for, at each given activity -"
         " the activity affects the continuum.</li>"
    "<li><b>Fixed at edges</b>: The channels on either side of the ROI are used to determine a"
         " linear continuum that is fixed, and not affected by the nuclides activity.</li>"
    "<li><b>Fixed full ROI</b>: The continuum is fit using the entire energy range of the ROI,"
         " assuming a Gaussian amplitude of zero; later, the Gaussian component of the peak will "
         " sit on top of this fixed continuum.<br/>"
         "This is effectively asserting that you know there is no signal peak present in the data."
         "  The continuum will not be affected by the nuclide activity value, and the Gaussian"
         " component of the peak will sit on top of this fixed continuum when evaluating the"
         " &chi;<sup>2</sup>.</li>"
    "</ul>";
      
    HelpSystem::attachToolTipOn( {label, m_decon_cont_norm_method},
                                tooltip, showToolTips, HelpSystem::ToolTipPosition::Right,
                                HelpSystem::ToolTipPrefOverride::RespectPreference );

    if( m_input.decon_cont_norm_method == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
    {
      assert( input.decon_continuum_type == PeakContinuum::OffsetType::Linear );
      m_input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
    }
    
    switch( m_input.decon_continuum_type )
    {
      case PeakContinuum::Linear:
        m_continuum->setCurrentIndex( 0 );
        break;
      case PeakContinuum::Quadratic:
        m_continuum->setCurrentIndex( 1 );
        break;
        
      case PeakContinuum::NoOffset:
      case PeakContinuum::Constant:
      case PeakContinuum::Cubic:
      case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep:
      case PeakContinuum::BiLinearStep:
      case PeakContinuum::External:
        assert( 0 );
        m_continuum->setCurrentIndex( 0 );
        m_input.decon_continuum_type = PeakContinuum::Linear;
        break;
    }//switch( input.offset_type )
    
    
    WContainerWidget *currieLimitContent = new WContainerWidget( rightColumn );
    currieLimitContent->addStyleClass( "MdaRowCurrieLimitContent" );
    
    WText *currie_label = new WText( "Single peak Limit", currieLimitContent );
    currie_label->addStyleClass( "MdaCurrieLimitTitle" );
    
    m_poisonLimit = new WText( "&nbsp;", currieLimitContent );
    m_poisonLimit->addStyleClass( "CurrieTxt");
    
    WContainerWidget *num_side_chan_row = new WContainerWidget( currieLimitContent );
    num_side_chan_row->addStyleClass( "CurrieNSideChanRow" );
    label = new WLabel( "Num Side Channels", num_side_chan_row );
    m_num_side_channels = new NativeFloatSpinBox( num_side_chan_row );
    m_num_side_channels->setFormatString( "%.0f" );
    m_num_side_channels->setSingleStep( 1.0 );
    m_num_side_channels->setRange( 1.0f, 50.0f );
    m_num_side_channels->setValue( input.num_side_channels );
    label->setBuddy( m_continuum );
  
    switch( m_input.limit_type )
    {
      case DetectionLimitTool::LimitType::Activity:
      {
        WPushButton *moreInfoButton = new WPushButton( currieLimitContent );
        moreInfoButton->setText( "further details..." );
        moreInfoButton->setStyleClass( "LinkBtn CurrieMoreInfoBtn" );
        moreInfoButton->clicked().connect( this, &MdaPeakRow::createMoreInfoWindow );
        
        break;
      }//case DetectionLimitTool::LimitType::Activity:
        
      case DetectionLimitTool::LimitType::Distance:
        // Not implemented
        break;
    }//switch( m_input.limit_type )
    
    
    m_use_for_likelihood->checked().connect( this, &MdaPeakRow::handleUseForLikelihoodChanged );
    m_use_for_likelihood->unChecked().connect( this, &MdaPeakRow::handleUseForLikelihoodChanged );
    m_roi_start->valueChanged().connect( this, &MdaPeakRow::roiChanged );
    m_roi_end->valueChanged().connect( this, &MdaPeakRow::roiChanged );
    m_num_side_channels->valueChanged().connect( this, &MdaPeakRow::roiChanged );
    m_continuum->activated().connect( this, &MdaPeakRow::roiChanged );
    m_continuum->changed().connect( this, &MdaPeakRow::roiChanged );
    m_decon_cont_norm_method->activated().connect( this, &MdaPeakRow::handleContinuumNormMethodChange );
    m_decon_cont_norm_method->changed().connect( this, &MdaPeakRow::handleContinuumNormMethodChange );
    
    
    setRoiStart( input.roi_start );
    setRoiEnd( input.roi_end );
    
    m_continuum->disable();
    m_decon_cont_norm_method->disable();
    
    setSimplePoisonTxt();
    
    //Right now we are just having DetectionLimitTool completely refresh on activity units change,
    //  but we could be a little more fine-grained about this.
    //InterSpec *viewer = InterSpec::instance();
    //InterSpecUser::addCallbackWhenChanged( viewer->m_user, "DisplayBecquerel",
    //                                      boost::bind(&MdaPeakRow::setSimplePoisonTxt, this) );
  }//MdaPeakRow constructor
  
  
  double simple_excess_counts() const
  {
    return m_simple_excess_counts;
  }
  
  double simple_mda() const
  {
    return m_simple_mda;
  }
  
  double simple_max_det_dist() const
  {
    return m_simple_max_det_dist;
  }
  
  Wt::Signal<> &changed()
  {
    return m_changed;
  }
  
  void setRoiStart( const float energy )
  {
    const float value = roundToNearestChannelEdge(energy);
    assert( roundToNearestChannelEdge(value) == value );
    const float oldval = m_roi_start->value();
    m_roi_start->setValue( value );
    m_input.roi_start = value;
    if( value != oldval )
      emitChanged();
  }//void setRoiStart( float energy )
  
  
  void setRoiEnd( const float energy )
  {
    const float value = roundToNearestChannelEdge(energy);
    assert( roundToNearestChannelEdge(value) == value );
    const float oldval = m_roi_end->value();
    m_roi_end->setValue( value );
    m_input.roi_end = value;
    if( value != oldval )
      emitChanged();
  }//void setRoiEnd( float energy )
  

  const DetectionLimitTool::MdaPeakRowInput &input() const
  {
    // Some sanity checks to make sure `m_input` has been kept in=sync with GUI widgets
    assert( fabs( m_input.roi_start - m_roi_start->value() ) < 0.1 );
    assert( fabs( m_input.roi_end - m_roi_end->value() ) < 0.1 );
    
    //roi_start/roi_end should already be rounded to nearest bin edge
    assert( fabs(m_input.roi_start - roundToNearestChannelEdge( m_roi_start->value() ) ) < 0.01 );
    assert( fabs(m_input.roi_end - roundToNearestChannelEdge( m_roi_end->value() ) ) < 0.01 );

    assert( static_cast<int>(m_input.decon_cont_norm_method) == m_decon_cont_norm_method->currentIndex() );
    assert( m_input.num_side_channels == m_num_side_channels->value() );
    assert( (m_continuum->currentIndex() == 0)
           || (m_continuum->currentIndex() == 1) );
    assert( ((m_continuum->currentIndex() == 0) && (m_input.decon_continuum_type == PeakContinuum::OffsetType::Linear))
           || ((m_continuum->currentIndex() == 1) && (m_input.decon_continuum_type == PeakContinuum::OffsetType::Quadratic)) );
    assert( m_input.use_for_likelihood == m_use_for_likelihood->isChecked() );
    
    return m_input;
  }
  
  DetectionLimitCalc::CurieMdaInput currieInput() const
  {
    auto m = m_input.measurement;
    if( !m || (m->num_gamma_channels() < 16) )
      throw runtime_error( "No measurement." );
    
    const float roi_lower_energy = m_roi_start->value();
    const float roi_upper_energy = m_roi_end->value();
    
    const size_t nsidebin = static_cast<size_t>( std::round( std::max( 1.0f, m_num_side_channels->value() ) ) );
    const size_t nchannels = m->num_gamma_channels();
    
    DetectionLimitCalc::CurieMdaInput input;
    input.spectrum = m;
    input.gamma_energy = m_input.energy;
    input.roi_lower_energy = m_roi_start->value();
    input.roi_upper_energy = m_roi_end->value();
    input.num_lower_side_channels = nsidebin;
    input.num_upper_side_channels = nsidebin;
    input.detection_probability = m_input.confidence_level;
    input.additional_uncertainty = 0.0f;  // TODO: can we get the DRFs contribution to form this?
    
    return input;
  }//DetectionLimitCalc::CurieMdaInput currieInput() const
  
  
  void createMoreInfoWindow()
  {
    char buffer[256];
    snprintf( buffer, sizeof(buffer), "%s %.2f keV Info",
             (m_nuclide ? m_nuclide->symbol.c_str() : "null"), m_input.energy );
    
    SimpleDialog *dialog = new SimpleDialog( buffer );
    dialog->addButton( "Close" );
    
    try
    {
      if( !m_input.measurement )
        throw runtime_error( "No measurement" );
        
      const bool useCuries = use_curie_units();
      const bool fixed_geom = m_input.drf->isFixedGeometry();
      const DetectorPeakResponse::EffGeometryType det_geom = m_input.drf->geometryType();
      const bool air_atten = (m_input.do_air_attenuation && !fixed_geom);
      const double &distance = m_input.distance;
      const float &energy = m_input.energy;
      const double intrinsic_eff = m_input.drf->intrinsicEfficiency( energy );
      const double geom_eff = m_input.drf->fractionalSolidAngle( m_input.drf->detectorDiameter(), distance );
      const double det_eff = fixed_geom ? intrinsic_eff : m_input.drf->efficiency(energy, distance);
      const double counts_4pi = (air_atten ? m_input.counts_per_bq_into_4pi_with_air
                                           : m_input.counts_per_bq_into_4pi);
      const double gammas_per_bq = counts_4pi * det_eff;
      
      snprintf( buffer, sizeof(buffer), "%.1f%%", 100.0*m_input.confidence_level );
      const string confidence_level = buffer;
      
      const DetectionLimitCalc::CurieMdaInput input = currieInput();
      const DetectionLimitCalc::CurieMdaResult result = DetectionLimitCalc::currie_mda_calc( input );
      

      //Gross Counts in Peak Foreground: 1389.813
      WTable *table = new WTable( dialog->contents() );
      table->addStyleClass( "MdaCurrieMoreInfoTable" );
      
      WTableCell *cell = nullptr;
      
      /** A helper lambda to add tool tips to both columns of the most recently created row of the table */
      auto addTooltipToRow = [table]( const string &tt ){
        WTableCell *cell = table->elementAt( table->rowCount() - 1, 2 );
        
        WImage *img = new WImage( cell );
        img->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
        img->setStyleClass("Wt-icon GridFourthRow GridThirdCol GridJustifyEnd");
        img->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
        
        HelpSystem::attachToolTipOn( img, tt, true, HelpSystem::ToolTipPosition::Right,
                                    HelpSystem::ToolTipPrefOverride::AlwaysShow );
      };//addTooltipToRow
      
      
      switch( m_input.limit_type )
      {
        case DetectionLimitTool::LimitType::Activity:
        {
          if( result.source_counts > result.decision_threshold )
          {
            // There is enough excess counts that we would reliably detect this activity, so we will
            //  give the activity range.
            const float lower_act = result.lower_limit / gammas_per_bq;
            const float upper_act = result.upper_limit / gammas_per_bq;
            const float nominal_act = result.source_counts / gammas_per_bq;
            
            const string lowerstr = PhysicalUnits::printToBestActivityUnits( lower_act, 2, useCuries )
                                    + det_eff_geom_type_postfix( det_geom );
            const string upperstr = PhysicalUnits::printToBestActivityUnits( upper_act, 2, useCuries )
                                    + det_eff_geom_type_postfix( det_geom );
            const string nomstr = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, useCuries )
                                  + det_eff_geom_type_postfix( det_geom );
            
            cell = table->elementAt( table->rowCount(), 0 );
            new WText( "Observed activity", cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( nomstr, cell );
            addTooltipToRow( "Greater than the &quot;critical level&quot;, L<sub>c</sub>,"
                            " counts were observed in the peak region." );
            
            cell = table->elementAt( table->rowCount(), 0 );
            new WText( "Activity range", cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( "[" + lowerstr + ", " + upperstr + "]", cell );
            addTooltipToRow( "The activity range estimate, to the " + confidence_level + " confidence level." );
          }else if( result.upper_limit < 0 )
          {
            // This can happen when there are a lot fewer counts in the peak region than predicted
            //  from the sides - since this is non-sensical, we'll just say zero.
            const string unitstr = useCuries ? "Ci" : "Bq";
            cell = table->elementAt( table->rowCount(), 1 );
            cell->setColumnSpan( 2 );
            new WText( "Activity &le; 0 " + unitstr, cell );
            addTooltipToRow( "Significantly fewer counts in peak region were observed,"
                            " than predicted by the neighboring regions." );
          }else
          {
            // We will provide the upper bound on activity.
            const string mdastr = PhysicalUnits::printToBestActivityUnits( m_simple_mda, 2, useCuries )
                                  + det_eff_geom_type_postfix( det_geom );
            
            cell = table->elementAt( table->rowCount(), 0 );
            new WText( "Activity upper bound", cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( mdastr , cell);
            
            addTooltipToRow( "The upper limit on how much activity could be present, to the "
                            + confidence_level + " confidence level." );
          }
          
          break;
        }//case DetectionLimitTool::LimitType::Activity:
          
        case DetectionLimitTool::LimitType::Distance:
        {
          WText *msg = new WText( "Distance not supported for additional info yet", dialog->contents() );
          msg->addStyleClass( "content" );
          msg->setInline( false );
          break;
        }
      }//switch( m_input.limit_type )
      
      // Add a blank row
      cell = table->elementAt( table->rowCount(), 0 );
      WText *txt = new WText( "&nbsp;", TextFormat::XHTMLText, cell );
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Lower region channels", cell );
      string val = "[" + std::to_string(result.first_lower_continuum_channel) + ", "
                    + std::to_string(result.last_lower_continuum_channel) + "]";
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      
      const double lower_lower_energy = m_input.measurement->gamma_channel_lower( result.first_lower_continuum_channel );
      const double lower_upper_energy = m_input.measurement->gamma_channel_lower( result.last_lower_continuum_channel );
      snprintf( buffer, sizeof(buffer), "The region above the peak in energy, that is being used"
               " to estimate the expected continuum counts in the peak region;"
               " corresponds to %.2f to %.2f keV", lower_lower_energy, lower_upper_energy );
      addTooltipToRow( buffer );
      
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Lower region counts", cell );
      val = SpecUtils::printCompact( result.lower_continuum_counts_sum, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The number of counts observed in the region below the peak region,"
                      " that is being used to estimate expected peak-region expected counts" );
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Upper region channels", cell );
      val = "[" + std::to_string(result.first_upper_continuum_channel) + ", "
                    + std::to_string(result.last_upper_continuum_channel) + "]";
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      const double upper_lower_energy = m_input.measurement->gamma_channel_lower( result.first_upper_continuum_channel );
      const double upper_upper_energy = m_input.measurement->gamma_channel_lower( result.last_lower_continuum_channel );
      snprintf( buffer, sizeof(buffer), "The region above the peak in energy, that is being used"
               " to estimate the expected continuum counts in the peak region;"
               " corresponds to %.2f to %.2f keV", upper_lower_energy, upper_upper_energy );
      addTooltipToRow( buffer );
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Upper region counts", cell );
      val = SpecUtils::printCompact( result.upper_continuum_counts_sum, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The number of counts observed in the region above the peak region,"
                      " that is being used to estimate expected peak-region expected counts" );
      
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Peak area channels", cell );
      val = "[" + std::to_string(result.last_lower_continuum_channel + 1) + ", "
                    + std::to_string(result.first_upper_continuum_channel - 1) + "]";
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      
      const double peak_lower_energy = m_input.measurement->gamma_channel_lower( result.last_lower_continuum_channel + 1 );
      const double peak_upper_energy = m_input.measurement->gamma_channel_lower( result.first_upper_continuum_channel - 1 );
      snprintf( buffer, sizeof(buffer), "The region the peak is being assumed to be within;"
               " corresponds to %.2f to %.2f keV", peak_lower_energy, peak_upper_energy );
      addTooltipToRow( buffer );
      
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Peak region counts", cell );
      val = SpecUtils::printCompact( result.peak_region_counts_sum, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The observed number of counts in the peak region" );
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Peak region null est.", cell );
      val = SpecUtils::printCompact( result.estimated_peak_continuum_counts, 5 )
            + " &plusmn; " + SpecUtils::printCompact( result.estimated_peak_continuum_uncert, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, TextFormat::XHTMLText, cell );
      addTooltipToRow( "An estimate of the expected number of counts, in the peak region, if it"
                      " is assumed that no signal is present.");
      
      // I believe this quantity corresponds to Currie's "critical level" ( L_c ),
      //   the net signal level (instrument response) above which an observed signal may be
      //   reliably recognized as "detected"
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Peak critical limit", cell );
      const double decision_threshold_act = result.decision_threshold / gammas_per_bq;
      val = SpecUtils::printCompact( result.decision_threshold, 4 )
            + " <span style=\"font-size: smaller;\">("
            + PhysicalUnits::printToBestActivityUnits( decision_threshold_act, 2, useCuries )
            + det_eff_geom_type_postfix( det_geom )
            + ")</span>";
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "Corresponds to Currie's &quot;critical level&quot;, L<sub>c</sub>,"
                      " that is the net signal level (instrument response) above which an"
                      " observed signal may be reliably recognized as &quot;detected&quot;." );
      
      
      // Note: I believe this quantity corresponds to Currie's "detection limit" (L_d) that
      //       is the “true” net signal level which may be a priori expected to lead to detection.
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Peak detection limit", cell );
      const double detection_limit_act = result.detection_limit / gammas_per_bq;
      val = SpecUtils::printCompact( result.detection_limit, 4 )
            + " <span style=\"font-size: smaller;\">("
            + PhysicalUnits::printToBestActivityUnits( detection_limit_act, 2, useCuries )
            + det_eff_geom_type_postfix( det_geom )
            + ")</span>";
      
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "Corresponds to Currie's &quot;detection limit&quot;, L<sub>d</sub>,"
                      " that is the &quot;true&quot; net signal level which may be, <i>a priori</i>."
                      " expected to lead to detection." );
      
      
      // Add a blank row
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "&nbsp;", TextFormat::XHTMLText, cell );
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Detector Intrinsic Eff.", cell );
      val = SpecUtils::printCompact( intrinsic_eff, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The efficiency for a gamma hitting the detector face,"
                      " to be detected in the full-energy peak." );
      
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Solid angle fraction", cell );
      val = SpecUtils::printCompact( geom_eff, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The fraction of the solid angle, the detector face takes up, at the specified distance." );
      
      const double shield_trans = m_input.counts_per_bq_into_4pi / m_input.branch_ratio / m_input.measurement->live_time();
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Shielding transmission", cell );
      val = SpecUtils::printCompact( shield_trans, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The fraction of gammas, at this energy, that will make it through the shielding without interacting." );
      
      const double air_trans = m_input.counts_per_bq_into_4pi_with_air / m_input.counts_per_bq_into_4pi;
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Air transmission", cell );
      val = SpecUtils::printCompact( air_trans, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The fraction of gammas, at this energy, that will make it through the air (assuming sea level) without interacting." );
      
      
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Nuclide branching ratio", cell );
      val = SpecUtils::printCompact( m_input.branch_ratio, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The number of gamma rays emitted at this energy, from the radioactive"
                      " source before any shielding, but accounting for nuclide age,"
                      " per decay of the parent nuclide." );
    }catch( std::exception &e )
    {
      WText *msg = new WText( "Error computing Currie limit information", dialog->contents() );
      msg->addStyleClass( "content" );
      msg->setInline( false );
    }//try / catch
  }//void createMoreInfoWindow()
  
};//class MdaPeakRow




DetectionLimitWindow::DetectionLimitWindow( InterSpec *viewer,
                                                     MaterialDB *materialDB,
                                                     WSuggestionPopup *materialSuggest )
: AuxWindow( "Detection Confidence Tool",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
   | AuxWindowProperties::SetCloseable
   | AuxWindowProperties::DisableCollapse) ),
  m_tool( nullptr )
{
  assert( viewer );
  rejectWhenEscapePressed( true );
  
  m_tool = new DetectionLimitTool( viewer, materialDB, materialSuggest );
  WContainerWidget *content = contents();
  content->addStyleClass( "DetectionLimitWindowContent" );
  content->addWidget( m_tool );
  
  AuxWindow::addHelpInFooter( footer(), "detection-confidence-tool" );
  
  //WContainerWidget *buttonDiv = footer();
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  const int width = ((screenW < 906) ? screenW : 0.75*screenW );
  const int height = ((screenH < 525) ? screenH : 0.95*screenH );
  resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//DetectionLimitWindow(...) constrctor


DetectionLimitWindow::~DetectionLimitWindow()
{
}




DetectionLimitTool::DetectionLimitTool( InterSpec *viewer,
                                                  MaterialDB *materialDB,
                                                  Wt::WSuggestionPopup *materialSuggest,
                                                  WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_needsUpdate( true ),
    m_chart( nullptr ),
    m_peakModel( nullptr ),
    m_nuclideEdit( nullptr ),
    m_ageEdit( nullptr ),
    m_currentNuclide( nullptr ),
    m_currentAge( -1.0 ),
    m_nuclideSuggest( nullptr ),
    m_detectorDisplay( nullptr ),
    m_distOrActivity( nullptr ),
    m_activityLabel( nullptr ),
    m_distanceLabel( nullptr ),
    m_distanceForActivityLimit( nullptr ),
    m_activityForDistanceLimit( nullptr ),
    m_materialDB( materialDB ),
    m_materialSuggest( materialSuggest ),
    m_shieldingSelect( nullptr ),
    m_minRelIntensity( nullptr ),
    m_attenuateForAir( nullptr ),
    m_displayActivityLabel( nullptr ),
    m_displayDistanceLabel( nullptr ),
    m_displayActivity( nullptr ),
    m_displayDistance( nullptr ),
    m_confidenceLevel( nullptr ),
    m_results( nullptr ),
    m_chi2Chart( nullptr ),
    m_bestChi2Act( nullptr ),
    m_upperLimit( nullptr ),
    m_errorMsg( nullptr ),
    m_fitFwhmBtn( nullptr )
{
  WApplication * const app = wApp;
  assert( app );
  
  app->useStyleSheet( "InterSpec_resources/DetectionLimitTool.css" );
  app->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  app->require( "InterSpec_resources/DetectionLimitTool.js" );
  
  addStyleClass( "DetectionLimitTool" );
  
  
  //new WLabel( "DetectionLimitTool", this );
  const WLength fieldWidth(4,WLength::FontEm);
  const WLength optionWidth(5.25,WLength::FontEm), buttonWidth(5.25,WLength::FontEm);
  
  // Container to hold and m_chart and m_results (which itself holds m_chi2Chart,  m_bestChi2Act,
  //  and m_upperLimit).
  WContainerWidget *upperCharts = new WContainerWidget( this );
  upperCharts->addStyleClass( "MdaCharts" );
  

  m_chart = new D3SpectrumDisplayDiv( upperCharts );
  
  /*
  //const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  //const int width = ((screenW < 906) ? screenW : 0.75*screenW ); // maybe cap width off at 1800 px?
  const int height = ((screenH < 525) ? screenH : 0.8*screenH );
  //m_chart->resize( WLength(width-20,WLength::Unit::Pixel), WLength(std::min(0.3*height,300.0),WLength::Unit::Pixel) );
  //m_chart->resize( WLength::Auto, WLength(300,WLength::Unit::Pixel) );
  cout << "height-->" << 1.0*height << endl;
  m_chart->setMinimumSize( WLength::Auto, WLength( std::max(200.0,0.4*height),WLength::Unit::Pixel) );
  */
  
  /** \TODO: Right now the user can drag the ROI around, and if the peak amplitude becaomes
             insignificant or under an area of 5, then the old peak is used, just with ROI range
             altered.  So what we'll do is just catch when the peak ROI drag is final, and then
             adjust the ROI, and re-due our fit - we should probably do something more intuitive
             so the user knows whats going on.
   */
  m_chart->existingRoiEdgeDragUpdate().connect( boost::bind( &DetectionLimitTool::roiDraggedCallback, this, _1, _2, _3, _4, _5, _6 ) );
  
  
  m_chart->setCompactAxis( true );
  m_chart->setXAxisTitle( "Energy (keV)" );
  m_chart->setYAxisTitle( "Counts/Channel" );
  m_chart->disableLegend();
  m_chart->showHistogramIntegralsInLegend( true );
  
  m_chart->applyColorTheme( m_interspec->getColorTheme() );
  m_interspec->colorThemeChanged().connect( m_chart, &D3SpectrumDisplayDiv::applyColorTheme );
  
  
  m_peakModel = new PeakModel( this );
  m_chart->setPeakModel( m_peakModel );

  // Create the user-input area under the spectrum and liklihood chart
  WContainerWidget *inputTable = new WContainerWidget( this );
  inputTable->addStyleClass( "Inputs" );
  
  
  // Create nuclide label and input
  WLabel *label = new WLabel( "Nuclide:", inputTable );
  label->addStyleClass( "GridFirstCol GridFirstRow" );
  
  
  m_nuclideEdit = new WLineEdit( "", inputTable );
  m_nuclideEdit->setMargin( 1 );
  m_nuclideEdit->addStyleClass( "GridSecondCol GridFirstRow" );
  m_nuclideEdit->setMinimumSize( fieldWidth, WLength::Auto );
  m_nuclideEdit->setAutoComplete( false );
  m_nuclideEdit->changed().connect( this, &DetectionLimitTool::handleUserNuclideChange );
  label->setBuddy( m_nuclideEdit );
  
  // Create nuclide suggestions
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( m_nuclideEdit );
  isoSuggestModel->excludeReactions( true );
  isoSuggestModel->excludeEscapes( true );
  isoSuggestModel->excludeXrays( true );
  m_nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
  m_nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
  m_nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  m_nuclideSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
  
  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_nuclideSuggest );
  
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  m_nuclideSuggest->forEdit( m_nuclideEdit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon
  
  
  // Create age input and time input validator
  label = new WLabel( "Age:", inputTable );
  label->addStyleClass( "GridFirstCol GridSecondRow" );
  
  m_ageEdit = new WLineEdit( "", inputTable );
  m_ageEdit->addStyleClass( "GridSecondCol GridSecondRow" );
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, m_ageEdit );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_ageEdit->setValidator(validator);
  m_ageEdit->setAutoComplete( false );
  label->setBuddy( m_ageEdit );
  
  m_ageEdit->changed().connect( this, &DetectionLimitTool::handleUserAgeChange );
  m_ageEdit->blurred().connect( this, &DetectionLimitTool::handleUserAgeChange );
  m_ageEdit->enterPressed().connect( this, &DetectionLimitTool::handleUserAgeChange );
  
  
  
  // Create distance input
  label = new WLabel( "Distance:", inputTable );
  label->addStyleClass( "GridFirstCol GridThirdRow" );
  m_distanceLabel = label;

  m_distanceForActivityLimit = new WLineEdit( "100 cm", inputTable );
  m_distanceForActivityLimit->addStyleClass( "GridSecondCol GridThirdRow" );
  label->setBuddy( m_distanceForActivityLimit );
  m_distanceForActivityLimit->changed().connect( this, &DetectionLimitTool::handleInputChange );
  m_distanceForActivityLimit->enterPressed().connect( this, &DetectionLimitTool::handleInputChange );
  
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, m_distanceForActivityLimit );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_distanceForActivityLimit->setValidator( distValidator );
  
  
  // We will but the activity label/input right next to the distance stuff, but since we default to
  //  calculating activity limit, we'll hide the activity stuff.
  label = new WLabel( "Activity:", inputTable );
  label->addStyleClass( "GridFirstCol GridThirdRow" );
  m_activityLabel = label;
  label->hide();
  
  m_activityForDistanceLimit = new WLineEdit( "0 uCi", inputTable );
  m_activityForDistanceLimit->addStyleClass( "GridSecondCol GridThirdRow" );
  label->setBuddy( m_activityForDistanceLimit );
  m_activityForDistanceLimit->changed().connect( this, &DetectionLimitTool::handleInputChange );
  m_activityForDistanceLimit->enterPressed().connect( this, &DetectionLimitTool::handleInputChange );
  
  WRegExpValidator *actvalidator = new WRegExpValidator( PhysicalUnits::sm_activityRegex, this );
  actvalidator->setFlags(Wt::MatchCaseInsensitive);
  m_activityForDistanceLimit->setValidator( actvalidator );
  m_activityForDistanceLimit->hide();
  
  
  
  SpectraFileModel *specFileModel = viewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( viewer, specFileModel, inputTable );
  m_detectorDisplay->addStyleClass( "GridThirdCol GridFirstRow" );
  viewer->detectorChanged().connect( boost::bind( &DetectionLimitTool::handleDrfSelected, this, boost::placeholders::_1 ) );
  viewer->detectorModified().connect( boost::bind( &DetectionLimitTool::handleDrfSelected, this, boost::placeholders::_1 ) );
  
  
  m_shieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest, inputTable );
  m_shieldingSelect->addStyleClass( "GridThirdCol GridSecondRow GridSpanTwoRows" );
  m_shieldingSelect->materialEdit()->setEmptyText( "<shielding material>" );
  m_shieldingSelect->materialChanged().connect( this, &DetectionLimitTool::handleInputChange );
  m_shieldingSelect->materialModified().connect( this, &DetectionLimitTool::handleInputChange );
  m_shieldingSelect->setMinimumSize( WLength(250), WLength::Auto );
  
  
  SwitchCheckbox *loglin = new SwitchCheckbox( "Log", "Lin", inputTable );
  loglin->addStyleClass( "MdaChartLogLin GridFourthCol GridFirstRow GridSpanTwoCol" );
  loglin->unChecked().connect( boost::bind( &D3SpectrumDisplayDiv::setYAxisLog, m_chart, true ) );
  loglin->checked().connect( boost::bind( &D3SpectrumDisplayDiv::setYAxisLog, m_chart, false ) );
  
  
  m_attenuateForAir = new WCheckBox( "Attenuate for air", inputTable );
  m_attenuateForAir->addStyleClass( "GridFourthCol GridSecondRow GridSpanTwoCol" );
  m_attenuateForAir->setChecked( true );
  m_attenuateForAir->checked().connect(this, &DetectionLimitTool::handleUserChangedUseAirAttenuate );
  m_attenuateForAir->unChecked().connect(this, &DetectionLimitTool::handleUserChangedUseAirAttenuate );
  
  
  m_distOrActivity = new SwitchCheckbox( "Activity Limit", "Distance Limit", inputTable );
  m_distOrActivity->addStyleClass( "GridFourthCol GridThirdRow GridSpanTwoCol" );

  m_distOrActivity->checked().connect( this, &DetectionLimitTool::handleUserChangedToComputeActOrDist );
  m_distOrActivity->unChecked().connect( this, &DetectionLimitTool::handleUserChangedToComputeActOrDist );
  
  
  
  label = new WLabel( "Peaks disp. act.:", inputTable );
  label->addStyleClass( "GridSixthCol GridFirstRow" );
  m_displayActivityLabel = label;
  
  m_displayActivity = new WLineEdit( inputTable );
  m_displayActivity->addStyleClass( "GridSeventhCol GridFirstRow" );
  label->setBuddy( m_displayActivity );
  
  m_displayActivity->setValidator( actvalidator );
  m_displayActivity->setTextSize( 10 );
  m_displayActivity->setText( "0 uCi" );
  m_displayActivity->enterPressed().connect( this, &DetectionLimitTool::updateShownPeaks );
  m_displayActivity->changed().connect( this, &DetectionLimitTool::updateShownPeaks );
  
  
  // Like with user input, we will put the put the display distance stuff right next to activity,
  //  and hide the display distance stuff
  label = new WLabel( "Peaks disp. dist.:", inputTable );
  label->addStyleClass( "GridSixthCol GridFirstRow" );
  m_displayDistanceLabel = label;
  
  m_displayDistance = new WLineEdit( inputTable );
  m_displayDistance->addStyleClass( "GridSeventhCol GridFirstRow" );
  label->setBuddy( m_displayDistance );
  
  m_displayDistance->setValidator( distValidator );
  m_displayDistance->setTextSize( 10 );
  m_displayDistance->setText( "1m" );
  m_displayDistance->enterPressed().connect( this, &DetectionLimitTool::updateShownPeaks );
  m_displayDistance->changed().connect( this, &DetectionLimitTool::updateShownPeaks );
  
  m_displayDistanceLabel->hide();
  m_displayDistance->hide();
  
  
  label = new WLabel( "Confidence Level:", inputTable );
  label->addStyleClass( "GridSixthCol GridSecondRow" );
  m_confidenceLevel = new WComboBox( inputTable );
  m_confidenceLevel->addStyleClass( "GridSeventhCol GridSecondRow" );
  
  for( auto cl = ConfidenceLevel(0); cl < NumConfidenceLevel; cl = ConfidenceLevel(cl+1) )
  {
    const char *txt = "";
    
    switch( cl )
    {
      case ConfidenceLevel::OneSigma:   txt = "68%";     break;
      case ConfidenceLevel::TwoSigma:   txt = "95%";     break;
      case ConfidenceLevel::ThreeSigma: txt = "99%";     break;
      case ConfidenceLevel::FourSigma:  txt = "4-sigma"; break;
      case ConfidenceLevel::FiveSigma:  txt = "5-sigma"; break;
      case ConfidenceLevel::NumConfidenceLevel:          break;
    }//switch( cl )
    
    m_confidenceLevel->addItem( txt );
  }//for( loop over confidence levels )
  
  m_confidenceLevel->setCurrentIndex( ConfidenceLevel::TwoSigma );
  
  m_confidenceLevel->activated().connect(this, &DetectionLimitTool::handleUserChangedConfidenceLevel );
  
  
  const auto primaryMeas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  auto spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( primaryMeas && spec )
  {
    //Make a copy of things so we dont mess the real stuff up - eventually we may want to
    //  Fully copy the SpecMeas 
    auto ourspec = make_shared<SpecUtils::Measurement>( *spec );
    ourspec->set_detector_name( "Aa1" );
    m_our_meas = make_shared<SpecMeas>();
    m_our_meas->setDetector( primaryMeas->detector() );
    m_our_meas->add_measurement( ourspec, true );
    m_peakModel->setPeakFromSpecMeas( m_our_meas, {ourspec->sample_number()} );
    m_chart->setData( ourspec, false );
  }//if( spec )

  m_results = new WContainerWidget( upperCharts );
  m_results->addStyleClass( "MdaResults" );
  m_results->hide();
  
  m_chi2Chart = new WContainerWidget( m_results );
  m_chi2Chart->addStyleClass( "MdaChi2Chart" );
  
  //////////////////////////////////////////////////////////////////////
  m_bestChi2Act = new WText( "&nbsp;", m_results );
  m_bestChi2Act->setInline( false );
  m_bestChi2Act->addStyleClass( "MdaResultTxt" );
    
  m_upperLimit = new WText( "&nbsp;", m_results );
  m_upperLimit->setInline( false );
  m_upperLimit->addStyleClass( "MdaResultTxt" );
  
  m_errorMsg = new WText("&nbsp;", this );
  m_errorMsg->addStyleClass( "MdaErrMsg" );
  m_errorMsg->hide();
  
  m_fitFwhmBtn = new WPushButton( "Fit FWHM...", this );
  m_fitFwhmBtn->addStyleClass( "MdaFitFwhm LightButton" );
  m_fitFwhmBtn->clicked().connect( this, &DetectionLimitTool::handleFitFwhmRequested );
  m_fitFwhmBtn->hide();
  
  WContainerWidget *peaksHolder = new WContainerWidget( this );
  peaksHolder->addStyleClass( "MdaPeaksArea" );
  
  WContainerWidget *titleBar = new WContainerWidget( peaksHolder );
  titleBar->addStyleClass( "MdaPeaksTitleBar" );
  
  WText *peaksTitle = new WText( "Gamma lines to use", titleBar );
  peaksTitle->addStyleClass( "MdaPeaksTitle" );
  
  WContainerWidget *filterDiv = new WContainerWidget( titleBar );
  filterDiv->addStyleClass( "MdaPeaksFilter" );
  label = new WLabel( "Min relative intensity:", filterDiv );
  label->addStyleClass( "GridFourthCol GridSecondRow" );
  m_minRelIntensity = new NativeFloatSpinBox( filterDiv );
  m_minRelIntensity->addStyleClass( "GridFifthCol GridSecondRow" );
  m_minRelIntensity->setMinimum( 0.0f );
  m_minRelIntensity->setMaximum( 0.999f );
  m_minRelIntensity->setSpinnerHidden( true );
  label->setBuddy( m_minRelIntensity );
  m_minRelIntensity->valueChanged().connect(this, &DetectionLimitTool::handleUserMinRelativeIntensityChange );
  
  
  
  m_peaks = new WContainerWidget( peaksHolder );
  m_peaks->addStyleClass( "MdaPeaks" );
  
  
  ReferencePhotopeakDisplay *reflines = viewer->referenceLinesWidget();
  if( reflines )
  {
    const ReferenceLineInfo &current = reflines->currentlyShowingNuclide();
    if( current.m_nuclide )
      m_nuclideEdit->setText( current.m_nuclide->symbol );
    
    m_ageEdit->setText( current.m_input.m_age ); // I think this should be non-empty
    
    const ShieldingSelect *shielding = reflines->shieldingSelect();
    if( shielding )
    {
      const bool isGeneric = shielding->isGenericMaterial();
      if( isGeneric && (shielding->arealDensity() > 0.0) )
      {
        const double an = shielding->atomicNumber();
        const double ad = shielding->arealDensity();
        m_shieldingSelect->setAtomicNumberAndArealDensity( an, ad );
      }
      
      if( !isGeneric && shielding->materialEdit() && shielding->thicknessEdit() )
      {
        const string material = shielding->materialEdit()->text().toUTF8();
        const string thick = shielding->thicknessEdit()->text().toUTF8();
        
        if( !material.empty() && !thick.empty() )
          m_shieldingSelect->setMaterialNameAndThickness( material, thick );
      }//if( generic ) / else
    }//if( shielding )
  }//if( reflines )
  
  initChi2Chart();
  
  // Incase we do have an initial nuclide set, treat as if user entered new info.
  handleUserNuclideChange();
  
  // Update the displayed activity units, when the user changes this preference.
  InterSpecUser::addCallbackWhenChanged( viewer->m_user, "DisplayBecquerel",
                                        boost::bind(&DetectionLimitTool::handleInputChange, this) );
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  // A quick sanity check to make sure looked-up definition of "Air" matches the statically compiled
  //  one
  const Material *air = (m_materialDB ? m_materialDB->material( "Air" ) : nullptr);
  assert( air );
  
  /*
   Expected coefs:{
     {25.0keV,0.0545}, {37.5keV,0.0301}, {56.3keV,0.0230}, {84.4keV,0.0201}, {126.5keV,0.0180},
     {189.8keV,0.0160}, {284.8keV,0.0139}, {427.2keV,0.0119}, {640.7keV,0.0101}, {961.1keV,0.00834},
     {1441.6keV,0.00680}, {2162.4keV,0.00549}, {3243.7keV,0.00442}, {4865.5keV,0.00359},
     {7298.2keV,0.00297}
   }
   */
  
  
  for( float energy = 25; energy < 10000; energy *= 1.5 )
  {
    const double air_coef = GammaInteractionCalc::transmission_coefficient_air( energy, 1.0*PhysicalUnits::m );
    const double mat_coef = GammaInteractionCalc::transmition_coefficient_material( air, energy, 1.0*PhysicalUnits::m );
    
    const double diff = fabs(air_coef - mat_coef);
    if( diff > 0.00001*(std::max)(air_coef,mat_coef) )
    {
      char msg[512];
      snprintf( msg, sizeof(msg),
               "Static 'Air' definition differed from database: static=%f, database=%f, diff=%f",
               air_coef, mat_coef, diff );
      log_developer_error( BOOST_CURRENT_FUNCTION, msg );
      assert( 0 );
    }
  }
#endif  //PERFORM_DEVELOPER_CHECKS
}//DetectionLimitTool constructor
  

  
DetectionLimitTool::~DetectionLimitTool()
{
}


void DetectionLimitTool::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  if( renderFull || m_needsUpdate )
    doCalc();
  
  m_needsUpdate = false;
  
  WContainerWidget::render( flags );
}//render( flags )


void DetectionLimitTool::initChi2Chart()
{
  m_chi2Chart->setJavaScriptMember( "chart", "new MdaChi2Chart(" + m_chi2Chart->jsRef() + ", {});");
  const string jsgraph = m_chi2Chart->jsRef() + ".chart";
  
  m_chi2Chart->setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + m_chi2Chart->id() + "') && "
             + m_chi2Chart->jsRef() + " && " + jsgraph + " )"
          + jsgraph + ".redraw();"
        "}"
      "});"
  );
  
  m_chi2Chart->callJavaScriptMember( "resizeObserver.observe", m_chi2Chart->jsRef() );
}



void DetectionLimitTool::roiDraggedCallback( double new_roi_lower_energy,
                                                  double new_roi_upper_energy,
                                                  double new_roi_lower_px,
                                                  double new_roi_upper_px,
                                                  double original_lower_energy,
                                                  bool is_final_range )
{
  if( !is_final_range )
    return;
  
  if( !m_peakModel )  //Shouldnt ever happen
    return;
  
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> origpeaks = m_peakModel->peaks();
  if( !origpeaks || origpeaks->empty() )
    return;  //shouldnt ever happen
  
  
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    
    const DetectionLimitTool::MdaPeakRowInput &input = rw->input();
    
    const float this_roi_start = input.roi_start;
    if( fabs(this_roi_start - original_lower_energy) < 0.1 )
    {
      rw->setRoiStart( new_roi_lower_energy );
      rw->setRoiEnd( new_roi_upper_energy );
      
      cout << "Updated ROI of peak at " << rw->input().energy << " keV to go from "
           << new_roi_lower_energy << " to " << new_roi_upper_energy << " keV" << endl;
      return;
    }
  }//for( auto w : m_peaks->children() )
  
  #if( PERFORM_DEVELOPER_CHECKS )
    char msg[512];
    snprintf( msg, sizeof(msg), "Failed to find a ROI that started at %f keV", new_roi_lower_energy );
    log_developer_error( BOOST_CURRENT_FUNCTION, msg );
  #endif

  cerr << "Failed to find a ROI that started at " << new_roi_lower_energy <<  " keV" << endl;
}//void roiDraggedCallback(...)


void DetectionLimitTool::handleUserAgeChange()
{
  handleNuclideChange( false );
}//void handleUserAgeChange();



void DetectionLimitTool::handleUserNuclideChange()
{
  handleNuclideChange( true );
}//void handleUserNuclideChange()


void DetectionLimitTool::handleNuclideChange( const bool update_to_default_age )
{
  double age = -1.0;
  const string isotxt = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotxt );

  if( !nuc || nuc->isStable() || nuc->decaysToChildren.empty() || nuc->decaysToStableChildren() )
  {
    age = 0.0;
    m_ageEdit->setText( "0y" );
    m_ageEdit->disable();
    
    if( !nuc || nuc->isStable() )
    {
      nuc = nullptr;
      m_ageEdit->setText( "" );
      m_nuclideEdit->setText( "" );
      
      if( nuc )
        passMessage( isotxt + " is stable", WarningWidget::WarningMsgHigh );
    }
  }else if( update_to_default_age )
  {
    m_ageEdit->enable();
    
    string agestr;
    age = PeakDef::defaultDecayTime( nuc, &agestr );
    m_ageEdit->setText( agestr );
  }else
  {
    // We will check that the input age is reasonable, and if not, modify it; this isnt 100%,
    //  but good enough for now.
    m_ageEdit->enable();
    const string agestr = m_ageEdit->text().toUTF8();
    
    try
    {
      age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
      if( age > 100.0*nuc->halfLife || age < 0.0 )
        throw std::runtime_error( "" );
    }catch(...)
    {
      string def_age_str;
      age = PeakDef::defaultDecayTime( nuc, &def_age_str );
      passMessage( "Changed age to a more reasonable value for " + nuc->symbol
            + " from '" + agestr + "' to '" + def_age_str + "'", WarningWidget::WarningMsgLow );
      
      m_ageEdit->setText( def_age_str );
    }//try / catch
  }//if( !nuc ... ) / else else
  
  const bool newNuclide = (nuc != m_currentNuclide);
  const bool newAge = (age != m_currentAge);
  if( !newNuclide && !newAge )
    return;
  
  m_currentNuclide = nuc;
  m_currentAge = age;
  
  if( newNuclide )
    calcAndSetDefaultMinRelativeIntensity();
  
  handleInputChange();
}//void handleNuclideChange()


void DetectionLimitTool::handleUserMinRelativeIntensityChange()
{
  if( m_minRelIntensity->valueText().empty() )
    calcAndSetDefaultMinRelativeIntensity();
  
  float value = m_minRelIntensity->value();
  if( value < 0 || value >= 1.0 )
  {
    calcAndSetDefaultMinRelativeIntensity();
    value = m_minRelIntensity->value();
  }
  
  handleInputChange();
}//void handleUserMinRelativeIntensityChange()



void DetectionLimitTool::handleUserChangedConfidenceLevel()
{
  handleInputChange();
}


void DetectionLimitTool::handleUserChangedUseAirAttenuate()
{
  handleInputChange();
}


void DetectionLimitTool::handleUserChangedToComputeActOrDist()
{
  const bool distanceLimit = m_distOrActivity->isChecked();
  
  m_activityLabel->setHidden( !distanceLimit );
  m_activityForDistanceLimit->setHidden( !distanceLimit );
  
  m_displayActivityLabel->setHidden( distanceLimit );
  m_displayActivity->setHidden( distanceLimit );
  
  m_distanceLabel->setHidden( distanceLimit );
  m_distanceForActivityLimit->setHidden( distanceLimit );
  
  m_displayDistanceLabel->setHidden( !distanceLimit );
  m_displayDistance->setHidden( !distanceLimit );
  
  const bool useCurie = use_curie_units();
  
  // A lamda to set the distance from one WLineEdit, to another, defaulting to "1m" if the
  //  input edit is invalid
  auto updateDistance = []( WLineEdit *fromEdit, WLineEdit *toEdit ){
    string diststr = fromEdit->text().toUTF8();
    
    try
    {
      if( PhysicalUnits::stringToDistance(diststr) <= PhysicalUnits::mm )
        throw std::exception();
    }catch( std::exception & )
    {
      diststr = "1 m";
    }
    
    toEdit->setText( diststr );
  };//updateDistance( from, to );
  
  
  // Similar to updateDistance, but for activity
  auto updateActivity = [useCurie]( WLineEdit *fromEdit, WLineEdit *toEdit ){
    string actstr = fromEdit->text().toUTF8();
    
    try
    {
      if( PhysicalUnits::stringToActivity(actstr) <= PhysicalUnits::bq )
        throw std::exception();
    }catch( std::exception & )
    {
      actstr = useCurie ? "1 uCi" : "37.0 kBq";
    }
    
    toEdit->setText( actstr );
  };//updateActivity( from, to );
  
  
  
  if( distanceLimit )
  {
    // Since we are switching from distance limit, lets start with the activity peaks are
    //  currently displayed for, if available, or 1 uCi otherwise.  Similarly for distance.
    updateActivity( m_displayActivity, m_activityForDistanceLimit );
    updateDistance( m_distanceForActivityLimit, m_displayDistance );
  }else //if( distanceLimit )
  {
    updateActivity( m_activityForDistanceLimit, m_displayActivity );
    updateDistance( m_displayDistance, m_distanceForActivityLimit );
  }// if( distanceLimit ) / else
  
  
  handleInputChange();
}//void handleUserChangedToComputeActOrDist()


void DetectionLimitTool::handleDrfSelected( std::shared_ptr<DetectorPeakResponse> new_drf )
{
  m_drf_cache = new_drf;
  
  // The DetectorDisplay::setDetector(...) function may not have been called yet because of
  //  the order of signal/slot connections - so we'll set that here to make sure we are up
  //  to date.
  m_detectorDisplay->setDetector( new_drf );
  
  handleInputChange();
}//void handleDrfChange()


void DetectionLimitTool::handleFitFwhmRequested()
{
  const bool use_auto_fit_peaks_too = true;
  m_interspec->fwhmFromForegroundWindow( use_auto_fit_peaks_too );
}//void handleFitFwhmRequested()


void DetectionLimitTool::calcAndSetDefaultMinRelativeIntensity()
{
  vector<GammaLineInfo> lines;
  
  try
  {
    lines = gammaLines();
  }catch( std::exception & )
  {
    // leave lines blank so value of 0.0 will be used
  }
  
  const size_t max_wanted_lines = 20;
  
  if( lines.size() <= max_wanted_lines )
  {
    m_minRelIntensity->setValue( 0.0f );
    return;
  }
  
  
  const shared_ptr<const DetectorPeakResponse> drf = detector();
  
  if( !drf || !drf->isValid() )
  {
    m_minRelIntensity->setValue( 0.0f );
    return;
  }
  
  double distance = -1.0;
  try
  {
    distance = currentDisplayDistance();
  }catch(...)
  {
  }
  
  const bool fixed_geom = drf->isFixedGeometry();
  const bool air_atten = (!fixed_geom && m_attenuateForAir->isChecked());
  
  if( (distance > 0.0) || fixed_geom )
  {
    for( GammaLineInfo &line : lines )
    {
      if( air_atten )
        line.gammas_into_4pi = line.gammas_4pi_after_air_attenuation;
      line.gammas_into_4pi *= (fixed_geom ? drf->intrinsicEfficiency(line.energy)
                                          : drf->efficiency( line.energy, distance ) );
    }//for( GammaLineInfo &line : lines )
  }//if( distance > 0.0 )
  
  // Sort 'lines' so the largest intensities come first
  std::sort( begin(lines), end(lines),
    []( const GammaLineInfo &lhs, const GammaLineInfo &rhs ) -> bool {
    return lhs.gammas_into_4pi > rhs.gammas_into_4pi;
  } );
  
  const double maxLineIntensity = lines.front().gammas_into_4pi;
  const double minWantedIntensity = lines[max_wanted_lines-1].gammas_into_4pi;
  
  // maxLineIntensity should always be valid, but jic
  if( IsInf(maxLineIntensity)
     || IsNan(maxLineIntensity)
     || (maxLineIntensity <= 0.0)
     || IsInf(minWantedIntensity)
     || IsNan(minWantedIntensity) )
  {
    m_minRelIntensity->setValue( 0.0f );
    return;
  }
  
  const double minRelIntensity = minWantedIntensity / maxLineIntensity;
   
  m_minRelIntensity->setValue( static_cast<float>(minRelIntensity) );
}//void calcAndSetDefaultMinRelativeIntensity()


void DetectionLimitTool::updateShownGammaLinesFromMinIntensity()
{
  handleInputChange();
}//void updateShownGammaLinesFromMinIntensity()


vector<DetectionLimitTool::GammaLineInfo> DetectionLimitTool::gammaLines() const
{
  if( !m_currentNuclide )
    throw runtime_error( "No nuclide specified" );
  
  double distance = 1.0 * PhysicalUnits::m;
  try
  {
    distance = currentDisplayDistance();
  }catch(...)
  {
  }
  
  //const shared_ptr<const DetectorPeakResponse> drf = m_our_meas->detector();
  //if( !drf || !drf->isValid() )
    //throw runtime_error( "DetectionLimitTool::gammaLines(): no DRF" );
  
  double air_distance = distance; //!< will potentially get corrected for shielding thickness
  
  float shielding_an = 0.0f, shielding_ad = 0.0f, shielding_thickness = 0.0f;
  shared_ptr<const Material> shielding_material;
  const bool generic_shielding = m_shieldingSelect->isGenericMaterial();
  if( generic_shielding )
  {
    shielding_an = m_shieldingSelect->atomicNumber();
    shielding_ad = m_shieldingSelect->arealDensity();
  }else
  {
    shielding_material = m_shieldingSelect->material();
    shielding_thickness = m_shieldingSelect->thickness();
    
    if( shielding_material )
      air_distance -= shielding_thickness;
  }//if( generic shielding ) / else
  
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( m_currentNuclide, 1.0E-3 * SandiaDecay::curie );
  
  const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( m_currentAge );
  const double parent_activity = mixture.activity( m_currentAge, m_currentNuclide );
  
  vector<SandiaDecay::EnergyRatePair> gammas = mixture.gammas( m_currentAge, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
  
  boost::function<double(float)> att_coef_fcn, air_atten_fcn;
    
  if( generic_shielding )
  {
    att_coef_fcn = boost::bind( &GammaInteractionCalc::transmition_coefficient_generic,
                                 shielding_an, shielding_ad, _1 );
  }else if( shielding_material && shielding_thickness > 0.0 )
  {
    att_coef_fcn = boost::bind( &GammaInteractionCalc::transmition_coefficient_material,
                                 shielding_material.get(), _1, shielding_thickness );
  }
    
  if( air_distance > 0.0 )
    air_atten_fcn = boost::bind( &GammaInteractionCalc::transmission_coefficient_air, _1, air_distance );
  
  
  vector<GammaLineInfo> lines;
  for( const auto &erp : gammas )
  {
    const double energy = erp.energy;
    const double decay_br = erp.numPerSecond / parent_activity;
    
    //br *= drf->efficiency( static_cast<float>(energy), static_cast<float>(distance) );
    const double shield_transmission = att_coef_fcn.empty() ? 1.0 : exp( -1.0*att_coef_fcn(energy) );
    const double air_transmission = air_atten_fcn.empty() ? 1.0 : exp( -1.0*air_atten_fcn(energy) );
      
    GammaLineInfo info;
    info.energy = energy;
    info.branching_ratio = decay_br;
    info.gammas_into_4pi = decay_br * shield_transmission;
    info.gammas_4pi_after_air_attenuation = decay_br * shield_transmission * air_transmission;
    info.air_transmission = air_transmission;
    info.shield_transmission = shield_transmission;
    
    lines.push_back( std::move(info) );
  }//for( const auto &erp : gammas )
  
  return lines;
}//vector<GammaLineInfo> gammaLines() const


void DetectionLimitTool::handleInputChange()
{
  // Grab currently displayed values incase its the same nuclide as before, and we are just updating
  //  things.
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    m_previousRoiValues[rw->input().energy] = rw->input();
  }
  
  // Clear out all the old results - we will replace all the widgets at the moment
  m_peaks->clear();
  m_results->hide();
  m_errorMsg->setText( "&nbsp;" );
  m_errorMsg->hide();
  
  auto spec = m_chart->data();
  if( !m_our_meas || !spec )
  {
    m_errorMsg->setText( "No Data Loaded" );
    m_errorMsg->show();
    return;
  }
  
  
  std::shared_ptr<const DetectorPeakResponse> drf = detector();
  
  if( !drf )
  {
    m_errorMsg->setText( "No DRF Loaded" );
    m_errorMsg->show();
    return;
  }
  
  if( !drf->hasResolutionInfo() || !drf->isValid() )
  {
    m_errorMsg->setText( "DRF does not have FWHM info - please fit for FWHM, or change DRF." );
    m_errorMsg->show();
    m_fitFwhmBtn->show();
    return;
  }
  m_fitFwhmBtn->hide();

  if( !m_currentNuclide )
  {
    m_errorMsg->setText( "No valid nuclide" );
    m_errorMsg->show();
    return;
  }
  
  const bool fixed_geom = drf->isFixedGeometry();
  const LimitType type = limitType();
  const bool do_air_atten = (!fixed_geom && m_attenuateForAir->isChecked());
  
  if( fixed_geom )
  {
    if( m_distOrActivity->isChecked() ) //checked means distance limit
    {
      m_distOrActivity->setUnChecked();
      m_distOrActivity->disable();
      handleUserChangedToComputeActOrDist();
      
      //`handleUserChangedToComputeActOrDist()` will call this function
      //  (`DetectionLimitTool::handleInputChange()`), but we wont get here, but we
      //  do have no need to continue on in our current function call.
      return;
    }//
  }else
  {
    m_distOrActivity->enable();
  }//if( m_distOrActivity->isChecked() ) / else
  
  double distance = 0.0, activity = 0.0;
  
  switch( type )
  {
    case LimitType::Distance:
    {
      // Currently using `MdaPeakRow::m_simple_mda` for distance @ 1 meter, is what we are doing
      //  so be careful if we change things.
      distance = 1.0*PhysicalUnits::meter;  //currentDisplayDistance();
      
      try
      {
        const string acttxt = m_activityForDistanceLimit->text().toUTF8();
        activity = PhysicalUnits::stringToActivity(acttxt);
        
        if( activity <= 0.0 )
          throw runtime_error( "Activity cant be zero or negative" );
      }catch( std::exception & )
      {
        m_errorMsg->setText( "Activity not valid" );
        m_errorMsg->show();
        return;
      }//try / catch
      
      break;
    }//case LimitType::Activity:
      
      
    case LimitType::Activity:
    {
      try
      {
        const string disttxt = m_distanceForActivityLimit->text().toUTF8();
        distance = PhysicalUnits::stringToDistance(disttxt);
        
        if( distance <= 0.0 )
          throw runtime_error( "Distance cant be zero or negative" );
      }catch( std::exception & )
      {
        m_errorMsg->setText( "Distance not valid" );
        m_errorMsg->show();
        return;
      }//try / catch
      
      break;
    }//case LimitType::Distance:
  }//switch( type )
  
  
  // We should have value min rel intensity value if we're here, but JIC, we'll check
  const float minDisplayIntensity = m_minRelIntensity->value();
  if( IsNan(minDisplayIntensity) || IsInf(minDisplayIntensity)
     || (minDisplayIntensity < 0.0f) || (minDisplayIntensity > 1.0f) )
  {
    m_errorMsg->setText( "Min relative intensity to show is invalid." );
    m_errorMsg->show();
    return;
  }

  
  setRefLinesAndGetLineInfo();
  
  vector<GammaLineInfo> lines; //{energy, br, gammas_into_4pi, gammas_4pi_after_air_attenuation, air_transmission, shield_transmission}

  try
  {
    lines = gammaLines();
  }catch(...)
  {
    m_errorMsg->setText( "Error Calculating effects of shielding" );
    m_errorMsg->show();
    return;
  }// try / catch
  
  double maxLineIntensity = 0.0;
  for( const auto &line : lines )
  {
    const double det_eff = fixed_geom ? drf->intrinsicEfficiency(line.energy)
                                      : drf->efficiency( line.energy, distance);
    
    // Note that line.gammas_into_4pi and line.gammas_4pi_after_air_attenuation already have
    //  shielding attenuation accounted for.
    const double intensity = det_eff * (do_air_atten ? line.gammas_4pi_after_air_attenuation : line.gammas_into_4pi);
    
    maxLineIntensity = std::max( maxLineIntensity, intensity );
  }//for( const auto &line : lines )
  
  if( maxLineIntensity <= 0.0 || IsInf(maxLineIntensity) || IsNan(maxLineIntensity) )
  {
    m_errorMsg->setText( "Error gamma yields - all intensities are zero." );
    m_errorMsg->show();
    return;
  }
  
  const float confLevel = currentConfidenceLevel();
  
  for( const auto &line : lines )
  {
    const double energy = line.energy;
    const double br = line.branching_ratio;
    const double det_eff = fixed_geom ? drf->intrinsicEfficiency(energy)
                                      : drf->efficiency( energy, distance);
    const double intensity = det_eff * (do_air_atten ? line.gammas_4pi_after_air_attenuation : line.gammas_into_4pi);
    
    const double relIntensity = intensity / maxLineIntensity;
    
    if( (relIntensity < minDisplayIntensity) || (relIntensity <= 0.0) )
    {
      continue;
    }
    
    const float fwhm = drf->peakResolutionFWHM(energy);
    
    MdaPeakRowInput input;
    input.use_for_likelihood = false;
    input.limit_type = type;
    input.do_air_attenuation = do_air_atten;
    input.measurement = spec;
    input.drf = drf;
    
    input.energy = energy;
    input.branch_ratio = br;
    input.counts_per_bq_into_4pi = line.gammas_into_4pi*spec->live_time();
    input.counts_per_bq_into_4pi_with_air = line.gammas_4pi_after_air_attenuation*spec->live_time();
    input.distance = distance;
    input.activity = activity;
    input.roi_start = energy - 1.25*fwhm; // recommended by ISO 11929:2010, could instead use 1.19
    input.roi_end = energy + 1.25*fwhm;
    input.num_side_channels = 4;
    input.confidence_level = confLevel;
    input.trans_through_air = line.air_transmission;
    input.trans_through_shielding = line.shield_transmission;
    input.decon_cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::Floating;
    input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
    
    const auto oldValPos = m_previousRoiValues.find( energy );
    if( oldValPos != end(m_previousRoiValues) )
    {
      // We will re-use the old user-selectable quantities from the MdaPeakRow widget, but we do not
      //  want to use re-use the options/values set outside of the MdaPeakRow widget.
      const MdaPeakRowInput &oldval = oldValPos->second;
      
      input.use_for_likelihood = oldval.use_for_likelihood;
      input.roi_start = oldval.roi_start;
      input.roi_end = oldval.roi_end;
      input.num_side_channels = oldval.num_side_channels;
      input.decon_cont_norm_method = oldval.decon_cont_norm_method;
      input.decon_continuum_type = oldval.decon_continuum_type;
      if( input.decon_cont_norm_method == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
      {
        assert( input.decon_continuum_type == PeakContinuum::OffsetType::Linear );
        input.decon_continuum_type = PeakContinuum::OffsetType::Linear;
      }
    }//if( we have seen this energy before )
    
    auto row = new MdaPeakRow( input, m_currentNuclide, m_peaks );
    
    row->changed().connect( this, &DetectionLimitTool::scheduleCalcUpdate );
  }//for( const auto &line : lines )
  
  m_results->show();
  
  scheduleCalcUpdate();
}//void handleInputChange()


float DetectionLimitTool::currentConfidenceLevel()
{
  const auto cl = ConfidenceLevel( m_confidenceLevel->currentIndex() );
  switch( cl )
  {
    case OneSigma:   return 0.682689492137086f;
    case TwoSigma:   return 0.954499736103642f;
    case ThreeSigma: return 0.997300203936740f;
    case FourSigma:  return 0.999936657516334f;
    case FiveSigma:  return 0.999999426696856f;
    case NumConfidenceLevel: break;
  }//switch( cl )
  
  assert( 0 );
  return 0.95f;
}//float DetectionLimitTool::currentConfidenceLevel()


double DetectionLimitTool::currentDisplayDistance() const
{
  string diststr;
  
  switch( limitType() )
  {
    case LimitType::Activity:
      diststr = m_distanceForActivityLimit->text().toUTF8();
      break;
   
    case LimitType::Distance:
      diststr = m_displayDistance->text().toUTF8();
      break;
  }//switch( limitType() )
  
  const double dist = PhysicalUnits::stringToDistance(diststr);
  if( dist <= 0.0 || IsNan(dist) || IsInf(dist) )
    throw runtime_error( "invalid distance" );
  
  return dist;
}//double currentDisplayDistance() const


double DetectionLimitTool::currentDisplayActivity() const
{
  string actstr;
  
  switch( limitType() )
  {
    case LimitType::Activity:
      actstr = m_displayActivity->text().toUTF8();
      break;
   
    case LimitType::Distance:
      actstr = m_activityForDistanceLimit->text().toUTF8();
      break;
  }//switch( limitType() )
  
  const double act = PhysicalUnits::stringToActivity(actstr);
  if( act <= 0.0 || IsNan(act) || IsInf(act) )
    throw runtime_error( "invalid activity" );
  
  return act;
}//double currentDisplayActivity() const


DetectionLimitTool::LimitType DetectionLimitTool::limitType() const
{
  return m_distOrActivity->isChecked() ? LimitType::Distance : LimitType::Activity;
}


void DetectionLimitTool::scheduleCalcUpdate()
{
  m_needsUpdate = true;
  scheduleRender();
}//void scheduleCalcUpdate()


void DetectionLimitTool::doCalc()
{
  try
  {
    const bool useCurie = use_curie_units();
  
    const bool is_dist_limit = (limitType() == DetectionLimitTool::LimitType::Distance);
  
    // If we are computing activity, then this next variable will be distance.
    // If we are computing distance, then it will be activity
    const double other_quantity = is_dist_limit ? currentDisplayActivity() : currentDisplayDistance();
    
    // We may be searching for either activity, or distance.  We'll shoehorn searching for either
    //  of these into this same function
    double min_search_quantity = std::numeric_limits<double>::infinity(), max_search_quantity = 0.0;
  
    int nused = 0;
    DetectorPeakResponse::EffGeometryType det_geom = DetectorPeakResponse::EffGeometryType::FarField;
    
    for( auto w : m_peaks->children() )
    {
      auto rw = dynamic_cast<MdaPeakRow *>( w );
      assert( rw );
      if( !rw->input().use_for_likelihood )
        continue;
    
      const MdaPeakRowInput &input = rw->input();
      const bool fixed_geom = input.drf->isFixedGeometry();
      const bool air_atten = (!fixed_geom && input.do_air_attenuation);
      det_geom = input.drf->geometryType(); //assumes all peaks have same DRF - which they should
      
      const double row_limit = std::max( 0.0, (is_dist_limit ? rw->simple_max_det_dist() : rw->simple_mda()) );
    
      ++nused;
    
      if( (row_limit > 0.0) && !IsInf(row_limit) && !IsNan(row_limit)  )
      {
        min_search_quantity = std::min( min_search_quantity, row_limit );
        max_search_quantity = std::max( max_search_quantity, row_limit );
        
        if( rw->simple_excess_counts() > 0.0 )
        {
          // `row_limit` is the upper limit, but we'll compute nominal value, to help increase range
          if( is_dist_limit )
          {
            // Ignore air attenuation for the moment; this is valid for finding the
            //  upper distance, and currently below we are setting the min distance
            //  to `min_allowed_quantity`, so it wont have any effect on that
            assert( !fixed_geom );
            const double def_dist = 1.0*PhysicalUnits::meter;
            const double det_eff_1m = fixed_geom ? input.drf->intrinsicEfficiency(input.energy)
            : input.drf->efficiency(input.energy, def_dist);
            const double counts_4pi = input.counts_per_bq_into_4pi;
            
            const double activity = other_quantity;
            
            const double nominal_dist = sqrt(counts_4pi * det_eff_1m * activity * def_dist * def_dist / rw->simple_excess_counts());
            
            min_search_quantity = std::min( min_search_quantity, nominal_dist );
            max_search_quantity = std::max( max_search_quantity, nominal_dist );
          }else
          {
            const double det_eff = fixed_geom ? input.drf->intrinsicEfficiency(input.energy)
            : input.drf->efficiency(input.energy, input.distance);
            const double counts_4pi = (air_atten ? input.counts_per_bq_into_4pi_with_air
                                       : input.counts_per_bq_into_4pi);
            const double gammas_per_bq = det_eff * counts_4pi;
            
            const double nominalActivity = rw->simple_excess_counts() / gammas_per_bq;
            min_search_quantity = std::min( min_search_quantity, nominalActivity );
            max_search_quantity = std::max( max_search_quantity, nominalActivity );
          }//if( is_dist_limit ) / else
        }
      }
    }//for( auto w : m_peaks->children() )
    
    if( !nused )
    {
      m_results->hide();
      if( m_errorMsg->isHidden() || m_errorMsg->text().empty() )
      {
        m_errorMsg->setText( "No peaks are selected" );
        m_errorMsg->show();
      }
      m_peakModel->setPeaks( vector<shared_ptr<const PeakDef>>{} );
      
      return;
    }//if( !nused )
    
    auto print_quantity = [is_dist_limit,det_geom,useCurie]( double quantity, int ndigits = 4 ) -> string {
      if( is_dist_limit )
        return PhysicalUnits::printToBestLengthUnits(quantity,ndigits);
      return PhysicalUnits::printToBestActivityUnits(quantity,ndigits,useCurie)
      + det_eff_geom_type_postfix(det_geom);
    };//print_quantity
    
    const bool shielding_has_thick = (!m_shieldingSelect->isGenericMaterial() && m_shieldingSelect->material());
    const double shield_thickness = (shielding_has_thick ? m_shieldingSelect->thickness() : 0.0);
    
    // A one week measurement, with a fully efficient detector, would be 0.6 counts
    //  for a uBq - a limit low enough its not likely to ever be encountered
    //const double min_allowed_quantity = 1.0E-6*PhysicalUnits::bq;
    const double min_allowed_quantity = is_dist_limit
                              ? (shielding_has_thick ? shield_thickness : 0.1*PhysicalUnits::mm)
                              : 0.0;
    
    if( min_search_quantity < min_allowed_quantity )
      min_search_quantity = min_allowed_quantity;
    
    if( IsInf(min_search_quantity) || (max_search_quantity == 0.0) )
    {
      min_search_quantity = min_allowed_quantity;
      //10 km or 1 kCi seem like reasonable large limits, but I'm sure at some point they wont be...
      max_search_quantity = is_dist_limit ? 10000.0*PhysicalUnits::meter : 1000.0*PhysicalUnits::curie;
    }else
    {
      // I *think* if we are here, there is a nominal activity
      //min_search_quantity *= 0.01;
      
      // It looks like multiplying `min_search_quantity` by something like 0.01 doesnt work super
      //  well since we are computing from very basic gross counts, but we want to find the limits
      //  using the peak-shape, which can vary wildly, say if we are right next to another peak,
      //  so for the time being, we'll just always start from `min_allowed_quantity`.
      min_search_quantity = min_allowed_quantity;
      max_search_quantity *= 100.0;
    }
    
    m_errorMsg->hide();
    
    const double yrange = 15;
    
    // TODO: we are scanning activity or distance, which is a single degree of freedom - but does it matter that
    //       we are marginalizing over (i.e., fitting for) the nuisance parameters of the peaks and
    //       stuff?  I dont *think* so.
    const boost::math::chi_squared chi_squared_dist( 1.0 );
    
    const float wantedCl = currentConfidenceLevel();
    
    // We want interval corresponding to 95%, where the quantile will give us CDF up to that
    //  point, so we actually want the quantile that covers 97.5% of cases.
    const float twoSidedCl = 0.5 + 0.5*wantedCl;
    
    const double cl_chi2_delta = boost::math::quantile( chi_squared_dist, twoSidedCl );
    
    const size_t nchi2 = 25;  //approx num chi2 to compute
    static_assert( nchi2 > 2, "We need at least two chi2" );
    
    vector<pair<double,double>> chi2s;
    double overallBestChi2 = 0.0, overallBestQuantity = 0.0;
    double upperLimit = 0.0, lowerLimit = 0.0;
    bool foundUpperCl = false, foundUpperDisplay = false, foundLowerCl = false, foundLowerDisplay = false;
    
    //
    double quantityRangeMin = 0.0, quantityRangeMax = 0.0;
    
    
    /// \TODO: currently all this stuff assumes a smooth continuously increasing Chi2 with increasing
    ///        activity, but this doesnt have to be the case, especially with quadratic continuums.
    
    std::atomic<size_t> num_iterations( 0 );
    
    //The boost::math::tools::bisect(...) function will make calls using the same value of activity,
    //  so we will cache those values to save some time.
    map<double,double> chi2cache;
    std::mutex chi2cache_mutex;
    
    const double mid_search_quantity = 0.5*(min_search_quantity + max_search_quantity);
    const double base_act = is_dist_limit ?  other_quantity : mid_search_quantity;
    const double base_dist = is_dist_limit ? mid_search_quantity : other_quantity;
    vector<PeakDef> base_peaks;
    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> base_input
                                  = getComputeForActivityInput( base_act, base_dist, base_peaks );
    if( !base_input )
      throw runtime_error( "missing input quantity." );
    
    
    auto compute_chi2 = [is_dist_limit,&base_input]( const double quantity, int *numDOF = nullptr ) -> double {
      assert( base_input );
      DetectionLimitCalc::DeconComputeInput input = *base_input;
      if( is_dist_limit )
        input.distance = quantity;
      else
        input.activity = quantity;
      const DetectionLimitCalc::DeconComputeResults results
                                            = DetectionLimitCalc::decon_compute_peaks( input );
      
      if( (results.num_degree_of_freedom == 0) && (results.chi2 == 0.0) )
        throw runtime_error( "No DOF" );
      
      if( numDOF )
        *numDOF = results.num_degree_of_freedom;
      
      return results.chi2;
    };//compute_chi2
    

    // This next lambda takes either distance or activity for its argument, depending which
    //  limit is being computed
    auto chi2ForQuantity = [this,&num_iterations,&chi2cache,&chi2cache_mutex,compute_chi2]( double const &quantity ) -> double {
      
      {//begin lock on chi2cache_mutex
        std::lock_guard<std::mutex> lock( chi2cache_mutex );
        const auto pos = chi2cache.find(quantity);
        if( pos != end(chi2cache) )
          return pos->second;
      }//end lock on chi2cache_mutex
      
      if( quantity < 0.0 )
        return std::numeric_limits<double>::max();
      
      const double chi2 = compute_chi2( quantity );
      
      ++num_iterations;
      
      {//begin lock on chi2cache_mutex
        std::lock_guard<std::mutex> lock( chi2cache_mutex );
        chi2cache.insert( std::pair<double,double>{quantity,chi2} );
      }//end lock on chi2cache_mutex
      
      return chi2;
    };//chi2ForQuantity
    
    
    /// `search_range` returns {best-chi2,best-quantity}
    auto search_range = [chi2ForQuantity]( double min_range, double max_range, boost::uintmax_t &max_iter ) -> pair<double, double> {
      // boost::brent_find_minima first evaluates the input at the range midpoint, then endpoints
      //  and so on - we could do a first pre-scan over the range to help make sure we dont miss
      //  a global minimum.
      
      
      //\TODO: if best activity is at min_search_quantity, it takes 50 iterations inside brent_find_minima
      //      to confirm; we could save this time by using just a little bit of intelligence...
      const int bits = 12; //Float has 24 bits of mantisa; should get us accurate to three significant figures
      
      return boost::math::tools::brent_find_minima( chi2ForQuantity, min_range, max_range, bits, max_iter );
    };//search_range lambda
    
    boost::uintmax_t max_iter = 100;  //this variable gets changed each use, so you need to reset it afterwards
    const pair<double, double> r = search_range( min_search_quantity, max_search_quantity, max_iter );
    
    overallBestChi2 = r.second;
    overallBestQuantity = r.first;
    
    cout << "Found min X2=" << overallBestChi2 << " with activity "
         << print_quantity(overallBestQuantity)
         << " and it took " << std::dec << num_iterations.load() << " iterations; searched from "
         << print_quantity(min_search_quantity)
         << " to " << print_quantity(max_search_quantity)
         << endl;
    
    //boost::math::tools::bracket_and_solve_root(...)
    auto chi2ForRangeLimit = [&chi2ForQuantity, overallBestChi2, yrange]( double const &quantity ) -> double {
      return chi2ForQuantity(quantity) - overallBestChi2 - yrange;
    };
    
    auto chi2ForCL = [&chi2ForQuantity, overallBestChi2, cl_chi2_delta]( double const &quantity ) -> double {
      return chi2ForQuantity(quantity) - overallBestChi2 - cl_chi2_delta;
    };
    
    //Tolerance is called with two values of quantity (activity or distance, depending which limit
    //  is being found); one with the chi2 below root, and one above
    auto tolerance = [chi2ForCL](double quantity_1, double quantity_2) -> bool{
      const double chi2_1 = chi2ForCL(quantity_1);
      const double chi2_2 = chi2ForCL(quantity_2);
      
      // \TODO: make sure tolerance is being used correctly - when printing info out for every call I'm not sure it is being used right... (but answers seem reasonable, so...)
      //cout << "Tolerance called with quantity_1=" << PhysicalUnits::printToBestActivityUnits(quantity_1,false)  //PhysicalUnits::printToBestLengthUnits(quantity_1)
      //     << ", quantity_2=" << PhysicalUnits::printToBestActivityUnits(quantity_2,4,false)
      //     << " ---> chi2_1=" << chi2_1 << ", chi2_2=" << chi2_2 << endl;
      
      return fabs(chi2_1 - chi2_2) < 0.025;
    };//tolerance(...)
    
    //cout << "chi2ForCL(min_search_quantity)=" << chi2ForCL(min_search_quantity) << endl;
    
    SpecUtilsAsync::ThreadPool pool;
    
    //Before trying to find lower-bounding activity, make sure the best value isnt the lowest
    //  possible value (i.e., zero in this case), and that if we go to the lowest possible value,
    //  that the chi2 will increase at least by cl_chi2_delta
    pool.post( [&lowerLimit,&quantityRangeMin,&foundLowerCl,&foundLowerDisplay,&num_iterations, //quantities we will modify
                 min_search_quantity,overallBestQuantity,overallBestChi2,yrange, //values we can capture by value
                 &tolerance,&chi2ForCL,&chi2ForQuantity,&print_quantity,&chi2ForRangeLimit //lamdas we will use
               ](){
      const double min_search_chi2 = chi2ForCL(min_search_quantity);
      if( (fabs(min_search_quantity - overallBestQuantity) > 0.001)
         && (min_search_chi2 > 0.0) )
      {
        pair<double,double> lower_val;
        
        boost::uintmax_t max_iter = 100;  //see note about needing to set before every use
        lower_val = boost::math::tools::bisect( chi2ForCL, min_search_quantity, overallBestQuantity, tolerance, max_iter );
        lowerLimit = 0.5*(lower_val.first + lower_val.second);
        foundLowerCl = true;
        cout << "lower_val CL activity="
             << print_quantity(lowerLimit)
             << " with Chi2(" << lowerLimit << ")=" << chi2ForQuantity(lowerLimit)
             << " (Best Chi2(" << overallBestQuantity << ")=" << overallBestChi2
             << "), num_iterations=" << std::dec << num_iterations.load() << " and search range from "
             << print_quantity(min_search_quantity)
             << " to "
             << print_quantity(overallBestQuantity)
             << endl;
        
        const double minActChi2 = chi2ForRangeLimit(min_search_quantity);
        if( minActChi2 < 0.0 )
        {
          quantityRangeMin = min_search_quantity;
          cout << "lower_val display activity being set to min_search_quantity (" 
               << min_search_quantity << "): minActChi2=" << minActChi2
               << ", with Chi2(" << quantityRangeMin << ")=" << chi2ForQuantity(quantityRangeMin) << endl;
        }else
        {
          try
          {
            boost::uintmax_t max_iter = 100;
            lower_val = boost::math::tools::bisect( chi2ForRangeLimit, min_search_quantity, lowerLimit, tolerance, max_iter );
            quantityRangeMin = 0.5*(lower_val.first + lower_val.second);
            foundLowerDisplay = true;
            cout << "lower_val display quantity=" << print_quantity(quantityRangeMin)
                 << " wih chi2=" << chi2ForQuantity(quantityRangeMin) << ", num_iterations=" << std::dec << num_iterations.load() << endl;
          }catch( std::exception &e )
          {
            const double delta_act = 0.1*(lowerLimit - quantityRangeMin);
            for( quantityRangeMin = lowerLimit; quantityRangeMin > 0; quantityRangeMin -= delta_act )
            {
              const double this_chi2 = chi2ForQuantity(quantityRangeMin);
              if( this_chi2 >= (overallBestChi2 + yrange) )
                break;
            }
            
            cout << "Couldnt find lower-limit of display range properly, so scanned down and found "
                 << print_quantity(quantityRangeMin)
                 << " where LowerLimit=" << print_quantity(lowerLimit)
                 << " and ActRangeMin=" << print_quantity(quantityRangeMin)
                 << " and BestActivity" << print_quantity(overallBestQuantity)
                 << endl;
          }//try / catch
        }//
      }else
      {
        lowerLimit = 0.0;
        //quantityRangeMin = overallBestQuantity;
        quantityRangeMin = min_search_quantity;
        cout << "lower_val activity/distance already at min" << endl;
      }//if( fabs(min_search_quantity - overallBestQuantity) > 0.001*PhysicalUnits::bq ) / else
    } );//pool.post( ... find lower limit ...)
    
    pool.post( [&upperLimit,&quantityRangeMax,&foundUpperCl,&foundUpperDisplay,&num_iterations, //quantities we will modify
                 max_search_quantity,overallBestQuantity,overallBestChi2,yrange,is_dist_limit,min_search_quantity, //values we can capture by value
                 &tolerance,&chi2ForCL,&chi2ForQuantity,&print_quantity,&chi2ForRangeLimit //lambdas we will use
               ](){
      const double max_search_chi2 = chi2ForCL(max_search_quantity);
      if( (fabs(max_search_quantity - overallBestQuantity) > 0.001)
         && (max_search_chi2 > 0.0)  )
      {
        pair<double,double> upper_val;
        boost::uintmax_t max_iter = 100;
        upper_val = boost::math::tools::bisect( chi2ForCL, overallBestQuantity, max_search_quantity, tolerance, max_iter );
        upperLimit = 0.5*(upper_val.first + upper_val.second);
        foundUpperCl = true;
        cout << "upper_val CL activity=" << print_quantity(upperLimit)
             << " wih chi2=" << chi2ForQuantity(upperLimit) << ", num_iterations=" << std::dec << num_iterations.load()
             << " and search range from " << print_quantity(overallBestQuantity)
             << " to "
             << print_quantity(max_search_quantity)
             << endl;
        
        const double maxSearchChi2 = chi2ForRangeLimit(max_search_quantity);
        if( maxSearchChi2 < 0.0 )
        {
          quantityRangeMax = max_search_quantity;
          cout << "upper_val display activity being set to max_search_quantity (" << max_search_quantity << "): maxSearchChi2=" << maxSearchChi2 << endl;
        }else
        {
          try
          {
            max_iter = 100;
            upper_val = boost::math::tools::bisect( chi2ForRangeLimit, upperLimit, max_search_quantity, tolerance, max_iter );
            quantityRangeMax = 0.5*(upper_val.first + upper_val.second);
            foundUpperDisplay = true;
            cout << "upper_val display quantity=" << print_quantity(quantityRangeMax)
                 << " wih chi2=" << chi2ForQuantity(quantityRangeMax) << ", num_iterations="
                 << std::dec << num_iterations.load() << endl;
          }catch( std::exception &e )
          {
            const double delta_act = std::max( 0.1*fabs(upperLimit - overallBestQuantity), 0.01*fabs(max_search_quantity - upperLimit) );
            for( quantityRangeMax = upperLimit; quantityRangeMax < max_search_quantity; quantityRangeMax -= delta_act )
            {
              const double this_chi2 = chi2ForQuantity(quantityRangeMax);
              if( this_chi2 >= (overallBestChi2 + yrange) )
                break;
            }
            
            cout << "Couldn't find upper-limit of display range properly, so scanned up and found "
                 << print_quantity(quantityRangeMax)
                 << " where UpperLimit Chi2(" << upperLimit << ")="
                 << print_quantity(upperLimit)
                 << " and ActRangeMax Chi2(" << quantityRangeMax << ")="
                 << print_quantity(quantityRangeMax)
                 << " and BestActivity Chi2(" << overallBestQuantity << ")="
                 << print_quantity(overallBestQuantity)
                 << endl;
          }//try / catch
        }
      }else
      {
        upperLimit = overallBestQuantity;
        quantityRangeMax = max_search_quantity;
        
        if( is_dist_limit )
        {
          // We might be at a huge distance, so lets find the distance at which we would start to
          //  kinda see something, ever so slightly
          try
          {
            auto chi2ForMinDelta = [&chi2ForQuantity, overallBestChi2, yrange]( double const &quantity ) -> double {
              return chi2ForQuantity(quantity) - overallBestChi2 - 0.01;
            };
            
            boost::uintmax_t max_iter = 100;
            const auto effective_upper_val = boost::math::tools::bisect( chi2ForMinDelta, min_search_quantity, overallBestQuantity, tolerance, max_iter );
            upperLimit = 0.5*(effective_upper_val.first + effective_upper_val.second);
            quantityRangeMax = upperLimit;
          }catch( std::exception &e )
          {
            assert( 0 );
          }
          //overallBestQuantity
        }//if( is_dist_limit )
        
        cout << "upper_val activity already at max" << endl;
      }
    } );//pool.post( ... find upper limit ...)
    
    pool.join();
    
    cout << "Found best chi2 and ranges with num_iterations=" << std::dec << num_iterations.load() << endl;

    assert( quantityRangeMin <= quantityRangeMax );
    if( quantityRangeMax < quantityRangeMin )
      std::swap( quantityRangeMin, quantityRangeMax );
    
    if( quantityRangeMax == quantityRangeMin )
    {
      assert( !foundLowerCl && !foundUpperCl );
      quantityRangeMin = 0.9*quantityRangeMin;
      quantityRangeMax = 1.1*quantityRangeMin;
    }
    
    const double initialRangeDelta = fabs(quantityRangeMax - quantityRangeMin);
    if( is_dist_limit && !foundUpperCl )
    {
      // This is when there are nearly zero or negative counts so the Chi2 will just stay flat
      //  at larger and larger distances; in this case we have set quantityRangeMax to be approx
      //  where you start getting a little effect, so now we'll add in a little area after
      //  this point so you can see the Chi2 curve is flattened out
      quantityRangeMax += 0.33 * initialRangeDelta;
    }
      
    if( foundLowerDisplay )
      quantityRangeMin = std::max( min_search_quantity, quantityRangeMin - 0.1*initialRangeDelta );
    
    if( foundUpperDisplay )
      quantityRangeMax = std::min( max_search_quantity, quantityRangeMax + 0.1*initialRangeDelta );
    
    const double quantity_delta = fabs(quantityRangeMax - quantityRangeMin) / nchi2;
    chi2s.resize( nchi2 );
    
    for( size_t i = 0; i < nchi2; ++i )
    {
      const double quantity = quantityRangeMin + quantity_delta*i;
      pool.post( [i, quantity, &chi2s, &chi2ForQuantity](){
        chi2s[i].first = quantity;
        chi2s[i].second = chi2ForQuantity(quantity);
      } );
    }
    pool.join();
    
    const double distance = is_dist_limit ? lowerLimit : other_quantity;
    const double activity = is_dist_limit ? other_quantity : upperLimit;
    
    int numDOF = 0;
    
    const string quantity_str = print_quantity(overallBestQuantity, 3);
    char buffer[128];
    
    string limit_str;
    if( foundLowerCl && foundUpperCl )
    {
      vector<PeakDef> peaks;
      double lowerQuantityChi2 = -999.9, upperQuantityChi2 = -999.9;
      if( is_dist_limit )
      {
        computeForActivity( other_quantity, lowerLimit, peaks, lowerQuantityChi2, numDOF );
        computeForActivity( other_quantity, upperLimit, peaks, upperQuantityChi2, numDOF );
      }else
      {
        computeForActivity( lowerLimit, other_quantity, peaks, lowerQuantityChi2, numDOF );
        computeForActivity( upperLimit, other_quantity, peaks, upperQuantityChi2, numDOF );
      }
      
      limit_str = print_quantity( overallBestQuantity, 3 );
      const string lower_limit_str = print_quantity( lowerLimit, 2 );
      const string upper_limit_str = print_quantity( upperLimit, 2 );
      
      // Chi2 at upper and lower limits *should* be the same, but since I dont totally trust
      //  everything yet, we'll allow showing a discrepancy so we can see something is up
      if( fabs(lowerQuantityChi2 - upperQuantityChi2) < 0.05*std::max(lowerQuantityChi2, upperQuantityChi2) )
        snprintf( buffer, sizeof(buffer), "%.1f", 0.5*(lowerQuantityChi2 + upperQuantityChi2) );
      else
        snprintf( buffer, sizeof(buffer), "%.1f and %.1f", lowerQuantityChi2, upperQuantityChi2 );
      
      const string chi2_str = buffer;
      
      //snprintf( buffer, sizeof(buffer), "%.1f%% coverage in [%s, %s], &chi;<sup>2</sup>=%s",
      //         0.1*std::round(1000.0*wantedCl), lower_limit_str.c_str(), upper_limit_str.c_str(),
      //         chi2_str.c_str() );
      
      snprintf( buffer, sizeof(buffer), "Between %s and %s at %.1f%% CL, &chi;<sup>2</sup>=%s",
               lower_limit_str.c_str(), upper_limit_str.c_str(),
               0.1*std::round(1000.0*wantedCl),
               chi2_str.c_str() );
    }else if( !foundLowerCl && !foundUpperCl )
    {
      limit_str = print_quantity( overallBestQuantity, 3 );
      snprintf( buffer, sizeof(buffer), "Error: failed upper and lower limits at %.1f%%",
               0.1*std::round(1000.0*wantedCl) );
    }else if( foundLowerCl )
    {
      vector<PeakDef> peaks;
      double lowerQuantityChi2 = -999.9;
      if( is_dist_limit )
      {
        computeForActivity( other_quantity, lowerLimit, peaks, lowerQuantityChi2, numDOF );
        
        limit_str = print_quantity( lowerLimit, 3 );
        const string print_limit_str = print_quantity( lowerLimit, 2 );
        
        //More stat-nerd-esk language, maybe, if its even correct, but lets print something
        //  easier to interpret, for the commoners, like myself.
        //snprintf( buffer, sizeof(buffer), "%.1f%% coverage at %s with &chi;<sup>2</sup>=%.1f",
        //         0.1*std::round(1000.0*wantedCl), print_limit_str.c_str(), lowerQuantityChi2 );
        
        snprintf( buffer, sizeof(buffer), "Can detect at %s at %.1f%% CL, &chi;<sup>2</sup>=%.1f",
                 print_limit_str.c_str(), 0.1*std::round(1000.0*wantedCl), lowerQuantityChi2 );
      }else
      {
        //computeForActivity( lowerLimit, other_quantity, peaks, lowerQuantityChi2, numDOF );
        snprintf( buffer, sizeof(buffer), "Error: Didn't find %.1f%% CL activity",
                 0.1*std::round(1000.0*wantedCl));
      }
    }else
    {
      assert( foundUpperCl );
      
      if( is_dist_limit )
      {
        snprintf( buffer, sizeof(buffer), "Error: Didn't find %.1f%% CL det. distance",
                 0.1*std::round(1000.0*wantedCl) );
      }else
      {
        vector<PeakDef> peaks;
        double upperQuantityChi2 = -999.9;
        computeForActivity( upperLimit, other_quantity, peaks, upperQuantityChi2, numDOF );
        limit_str = print_quantity(upperLimit,3);
        const string print_limit_str = print_quantity( upperLimit, 2 );
        //snprintf( buffer, sizeof(buffer), "%.1f%% coverage at %s with &chi;<sup>2</sup>=%.1f",
        //         0.1*std::round(1000.0*wantedCl), print_limit_str.c_str(), upperQuantityChi2 );
        snprintf( buffer, sizeof(buffer), "Can detect %s at %.1f%% CL, &chi;<sup>2</sup>=%.1f",
                 print_limit_str.c_str(), 0.1*std::round(1000.0*wantedCl), upperQuantityChi2 );
      }//if( is_dist_limit ) / else
    }
    
    // TODO: This `m_upperLimit` text could be moved to where the warnings is, ot the title area where the "Gamma lines to use" is
    //       This would give us more room, and maybe be more obvious to the user
    m_upperLimit->setText( WString::fromUTF8(buffer) );
    
    if( is_dist_limit )
    {
      m_displayDistance->setText( WString::fromUTF8(limit_str) );
      
      if( foundUpperCl )
      {
        snprintf( buffer, sizeof(buffer), "Best &chi;<sup>2</sup> of %.1f at %s, %i DOF",
                 overallBestChi2, quantity_str.c_str(), numDOF );
      }else
      {
        snprintf( buffer, sizeof(buffer), "&chi;<sup>2</sup> is %.1f at large distance, %i DOF",
                 overallBestChi2, numDOF );
      }
    }else
    {
      m_displayActivity->setText( WString::fromUTF8(limit_str) );
      
      snprintf( buffer, sizeof(buffer), "Best &chi;<sup>2</sup> of %.1f at %s, %i DOF",
               overallBestChi2, quantity_str.c_str(), numDOF );
    }//if( is_dist_limit ) / else
    
    m_bestChi2Act->setText( buffer );
    
    
    if( foundLowerCl && (lowerLimit > 0.0) )
    {
      //display lower limit to user...
    }

    
    const double avrg_quantity = 0.5*(chi2s.front().first + chi2s.back().first);
    const pair<string,double> &units = is_dist_limit
                    ? PhysicalUnits::bestLengthUnitHtml( avrg_quantity )
                    : PhysicalUnits::bestActivityUnitHtml( avrg_quantity, useCurie );

    double minchi2 = std::numeric_limits<double>::infinity();
    double maxchi2 = -std::numeric_limits<double>::infinity();
    
    Wt::Json::Object json;
    Wt::Json::Array chi2_xy;
    for( size_t i = 0; i < chi2s.size(); ++i )
    {
      minchi2 = std::min( minchi2, chi2s[i].second );
      maxchi2 = std::max( maxchi2, chi2s[i].second );
      
      Wt::Json::Object datapoint;
      datapoint["x"] = (chi2s[i].first / units.second);
      datapoint["y"] = chi2s[i].second;
      chi2_xy.push_back( std::move(datapoint) );
    }
    
    json["data"] = chi2_xy;
    
    
    double chi2range = maxchi2 - minchi2;
    if( IsInf(minchi2) || IsInf(maxchi2) ) //JIC
    {
      minchi2 = 0;
      maxchi2 = 100;
      chi2range = 0;
    }
    
    //json["xunits"] = WString::fromUTF8(units.first);
    
    string unit_str = units.first;
    SpecUtils::ireplace_all( unit_str, "&mu;", MU_CHARACTER );
    json["xtitle"] = WString::fromUTF8( (is_dist_limit ? "Distance (" : "Activity (") + unit_str + ")");
    json["ytitle"] = WString::fromUTF8( "\xCF\x87\xC2\xB2" ); //chi^2
    
    double xstart = (chi2s.front().first/units.second);
    double xend = (chi2s.front().first/units.second);
    //Blah blah blah, coarse labels to be round numbers
    //const double initialrange = xend - xstart;
    
    //json["xrange"] = Wt::Json::Array{ Wt::Json::Array( Wt::Json::Value(xstart), Wt::Json::Value(xend) ) };
    //json["yrange"] = Wt::Json::Array{ Wt::Json::Array( Wt::Json::Value(minchi2 - 0.1*chi2range), Wt::Json::Value(maxchi2 + 0.1*chi2range) ) };
      
    // TODO: draw line at upperLimit.
    
    std::shared_ptr<const ColorTheme> theme = m_interspec ? m_interspec->getColorTheme() : nullptr;
    if( theme )
    {
      auto csstxt = []( const Wt::WColor &c, const char * defcolor ) -> WString {
        return WString::fromUTF8( (c.isDefault() ? string(defcolor) : c.cssText()) );
      };
       
      json["lineColor"] = csstxt(theme->foregroundLine, "black");
      //json["axisColor"] = csstxt(theme->spectrumAxisLines, "black");
      json["chartBackgroundColor"] = csstxt(theme->spectrumChartBackground, "rgba(0,0,0,0)");
      //json["textColor"] = csstxt(theme->spectrumChartText, "black");
    }//if( theme )
    
    m_results->show();
    
    const string jsgraph = m_chi2Chart->jsRef() + ".chart";
    const string datajson = Wt::Json::serialize(json);
    m_chi2Chart->doJavaScript( jsgraph + ".setData(" + datajson + ");" );
    
    updateShownPeaks();
  }catch( std::exception &e )
  {
    m_bestChi2Act->setText( "" );
    m_upperLimit->setText( "" );
    m_results->hide();
    m_peakModel->setPeaks( vector<shared_ptr<const PeakDef>>{} );
    m_errorMsg->setText( "Error calculating Chi2: " + string(e.what()) );
    
    const string jsgraph = m_chi2Chart->jsRef() + ".chart";
    m_chi2Chart->doJavaScript( jsgraph + ".setData(null);" );
    
    m_errorMsg->show();
    return;
  }//try / catch
}//void doCalc()


void DetectionLimitTool::updateShownPeaks()
{
  cout << "Need to upgrade DetectionLimitTool::updateShownPeaks() for toggle between activity and distance" << endl;
#ifndef _WIN32
  #warning "Need to upgrade DetectionLimitTool::updateShownPeaks() for toggle between activity and distance"
#endif

  m_peakModel->removeAllPeaks();
  
  const bool distanceLimit = (limitType() == DetectionLimitTool::LimitType::Distance);
  
  string actStr;
  if( distanceLimit )
  {
    actStr = m_activityForDistanceLimit->text().toUTF8();
  }else
  {
    actStr = m_displayActivity->text().toUTF8();
  }
  
  double activity = 0;
  try
  {
    activity = PhysicalUnits::stringToActivity(actStr);
    
    if( activity <= 0.0 )
      throw runtime_error( "zero or negative activity." );
  }catch(...)
  {
    return;
  }//try / catch
  
  
  double distance = 0.0, air_distance = 0.0;
  
  try
  {
    distance = air_distance = currentDisplayDistance();
    
    if( !m_shieldingSelect->isGenericMaterial() && m_shieldingSelect->material() )
      air_distance -= m_shieldingSelect->thickness();
  }catch(...)
  {
    return;
  }
  
  
  int numDOF;
  double chi2;
  std::vector<PeakDef> peaks;
  computeForActivity( activity, distance, peaks, chi2, numDOF );
  m_peakModel->setPeaks( peaks );
    
  cout << "Done in updateShownPeaks()" << endl;
}//void DetectionLimitTool::updateShownPeaks()


shared_ptr<DetectionLimitCalc::DeconComputeInput> DetectionLimitTool::getComputeForActivityInput(
                                                                  const double activity,
                                                                  const double distance,
                                                                  std::vector<PeakDef> &peaks )
{
  peaks.clear();
  
  auto spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !spec )
  {
    cerr << "No displayed histogram!" << endl;
    return nullptr;
  }
  
  auto input = make_shared<DetectionLimitCalc::DeconComputeInput>();
  input->measurement = spec;
  input->drf = detector();
  
  input->include_air_attenuation = m_attenuateForAir->isChecked();
  input->distance = distance;
  input->activity = activity;
  if( !m_shieldingSelect->isGenericMaterial() && m_shieldingSelect->material() )
    input->shielding_thickness = m_shieldingSelect->thickness();
  
  for( auto w : m_peaks->children() )
  {
    auto rw = dynamic_cast<MdaPeakRow *>( w );
    assert( rw );
    
    if( !rw->input().use_for_likelihood )
      continue;
    
    const DetectionLimitTool::MdaPeakRowInput &row_input = rw->input();
    const float roi_start = row_input.roi_start;
    const float roi_end = row_input.roi_end;
    
    const size_t nsidebin = std::max( size_t(1), row_input.num_side_channels );
    
    shared_ptr<PeakContinuum> peak_continuum;
    
    const DetectionLimitCalc::DeconContinuumNorm cont_norm = row_input.decon_cont_norm_method;
    
    const PeakContinuum::OffsetType continuum_type = row_input.decon_continuum_type;
    switch( continuum_type )
    {
      case PeakContinuum::Linear:
      case PeakContinuum::Quadratic:
        break;
        
      case PeakContinuum::NoOffset:   case PeakContinuum::Constant:
      case PeakContinuum::Cubic:      case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep: case PeakContinuum::BiLinearStep:
      case PeakContinuum::External:
        throw logic_error( "DetectionLimitTool::computeForActivity: Unexpected continuum type" );
        break;
    }//switch( continuum_type )
    
    
    DetectionLimitCalc::DeconRoiInfo roi_info;
    
    roi_info.roi_start = roi_start;
    roi_info.roi_end = roi_end;
    roi_info.continuum_type = continuum_type;
    roi_info.cont_norm_method = cont_norm;
    roi_info.num_lower_side_channels = nsidebin;
    roi_info.num_upper_side_channels = nsidebin;
    
    // TODO: have all this peak-specific stuff go inside a loop
    DetectionLimitCalc::DeconRoiInfo::PeakInfo peak_info;
    
    peak_info.energy = row_input.energy;
    peak_info.fwhm = input->drf->peakResolutionFWHM(row_input.energy);
    peak_info.counts_per_bq_into_4pi = row_input.counts_per_bq_into_4pi;
    
    roi_info.peak_infos.push_back( peak_info );
    
    input->roi_info.push_back( roi_info );
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  if( input->roi_info.empty() )
  {
    cerr << "No peaks to do calc for!" << endl;
    return nullptr;
  }
  
  return input;
}//getComputeForActivityInput(...)
 

void DetectionLimitTool::computeForActivity( const double activity,
                                              const double distance,
                                                   std::vector<PeakDef> &peaks,
                                                   double &chi2, int &numDOF )
{
  chi2 = 0.0;
  numDOF = 0;
  
  auto input = getComputeForActivityInput(activity, distance, peaks );
  
  if( !input )
    return;
  
  const DetectionLimitCalc::DeconComputeResults results
                                        = DetectionLimitCalc::decon_compute_peaks( *input );
  
  peaks = results.fit_peaks;
  chi2 = results.chi2;
  numDOF = results.num_degree_of_freedom;
}//void DetectionLimitTool::computeForActivity(...)


std::shared_ptr<const DetectorPeakResponse> DetectionLimitTool::detector()
{
  std::shared_ptr<DetectorPeakResponse> drf = m_detectorDisplay->detector();
  if( drf && drf->isValid() && drf->hasResolutionInfo() )
    return drf;
  
  if( !m_our_meas )
    return drf;
  
  
  if( m_our_meas->detector()
     && m_our_meas->detector()->isValid()
     && m_our_meas->detector()->hasResolutionInfo() )
  {
    drf = m_our_meas->detector();
    m_detectorDisplay->setDetector( drf );
    return drf;
  }
  
  //m_detectorDisplay->setDetector( nullptr );
  return drf;
}//detector()


void DetectionLimitTool::setRefLinesAndGetLineInfo()
{
  if( !m_currentNuclide )
    return;
  
  auto spec = m_chart->data();
  if( !m_our_meas || !spec )
    return;
    
  std::shared_ptr<const DetectorPeakResponse> drf = detector();
  if( !drf || !drf->hasResolutionInfo() || !drf->isValid() )
    return;

  
  double distance = 0.0, air_distance = 0.0;
  try
  {
    air_distance = distance = currentDisplayDistance();
  }catch( std::exception & )
  {
    // For the sake of displaying *something*, we'll just use a meter if the user hasnt entered a
    //  valid distance
    air_distance = distance = 1.0*PhysicalUnits::meter;
  }//try / catch
  
  
  if( m_currentAge < 0.0 )
    return;
  
  const double brCutoff = 0.0;
  double shielding_an = 0.0, shielding_ad = 0.0, shielding_thickness = 0.0;
  std::shared_ptr<const Material> shielding_material;
  const bool generic_shielding = m_shieldingSelect->isGenericMaterial();
  if( generic_shielding )
  {
    shielding_an = m_shieldingSelect->atomicNumber();
    shielding_ad = m_shieldingSelect->arealDensity();
  }else
  {
    shielding_material = m_shieldingSelect->material();
    shielding_thickness = m_shieldingSelect->thickness();
    
    if( shielding_material )
      air_distance -= shielding_thickness;
  }//if( generic shielding ) / else
  
  auto theme = m_interspec->getColorTheme();
  
  
  
  RefLineInput ref_input;
  ref_input.m_input_txt = m_currentNuclide->symbol;
  ref_input.m_age = m_ageEdit->text().toUTF8();
    
  if( theme && theme->referenceLineColor.size() )
    ref_input.m_color = theme->referenceLineColor[0];
  
  ref_input.m_lower_br_cutt_off = 0.0;
  ref_input.m_promptLinesOnly = false;
    
  ref_input.m_showGammas = true;
  ref_input.m_showXrays = false;
  ref_input.m_showAlphas = false;
  ref_input.m_showBetas = false;
  ref_input.m_showCascades = false;
    
  ref_input.m_detector_name = drf->name();
  
  const std::function<float( float )> intrinsic_eff = drf->intrinsicEfficiencyFcn();
  
  if( m_attenuateForAir->isChecked() && (air_distance > 0.0) )
  {
    if( !m_materialDB )
      throw std::runtime_error( "DetectionLimitTool::setRefLinesAndGetLineInfo(): no material DB" );
    
    const auto air_atten_fcn = boost::bind( &GammaInteractionCalc::transmission_coefficient_air, _1, air_distance );
    
    ref_input.m_det_intrinsic_eff = [intrinsic_eff,air_atten_fcn]( float energy ) -> float {
      return intrinsic_eff(energy) * exp( -1.0 * air_atten_fcn( energy ) );
    };
  }else
  {
    ref_input.m_det_intrinsic_eff = intrinsic_eff;
  }

  
  if( generic_shielding )
  {
    ref_input.m_shielding_an = m_shieldingSelect->atomicNumberEdit()->text().toUTF8();
    ref_input.m_shielding_ad = m_shieldingSelect->arealDensityEdit()->text().toUTF8();
  }else
  {
    ref_input.m_shielding_name = m_shieldingSelect->materialEdit()->text().toUTF8();
    ref_input.m_shielding_thickness = m_shieldingSelect->thicknessEdit()->text().toUTF8();
  }
  
  
  std::shared_ptr<ReferenceLineInfo> reflines = ReferenceLineInfo::generateRefLineInfo( ref_input );
  if( reflines )
    m_chart->setReferncePhotoPeakLines( *reflines );
  else
    m_chart->clearAllReferncePhotoPeakLines();
}//std::vector<std::tuple<float,double,bool>> setRefLinesAndGetLineInfo();



void DetectionLimitTool::do_development()
{
  auto specfile = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !specfile )
  {
    new WLabel( "No Foreground File", this );
    return;
  }
  
  auto spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !spec )
  {
    new WLabel( "No foreground spectrum", this );
    return;
  }
  
  std::shared_ptr<DetectorPeakResponse> drf = specfile->detector();
  if( !drf )
  {
    new WLabel( "No DRF loaded", this );
    return;
  }
  
  if( !drf->hasResolutionInfo() || !drf->isValid() )
  {
    new WLabel( "DRF not valid or doesnt have FWHM info.", this );
    return;
  }
  
  
  vector<PeakDef> inputPeaks, fitPeaks;
  
  //vector<float> energies = { 80.9971, 276.4, 302.853, 356.017, 383.848 };
  vector<float> energies = { 356.017 };
  //vector<tuple<float,float>> peak_info = { {356.017f,0.0f} };
  
  //We provide <energy,gammas_per_uci> and optionally activity, and get out chi2, DOF, and activity if not provided
  //First, fit for preffered activity.
  //Then increase the activity until the Chi2 increases by quantile(chi_squared(n_pars),0.95)
  
  
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const float mean = energies[i];
    const float sigma = drf->peakResolutionSigma(mean);
    const float amplitude = 0.0f;
    PeakDef peak( mean, sigma, amplitude );
    peak.setFitFor( PeakDef::CoefficientType::Mean, false );
    peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
    peak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, true );
    shared_ptr<PeakContinuum> cont = peak.continuum();
    cont->setType( PeakContinuum::OffsetType::Linear );
    
    double lowerEnengy, upperEnergy;
    findROIEnergyLimits( lowerEnengy, upperEnergy, peak, spec );
    cont->setRange( lowerEnengy, upperEnergy );
    cont->calc_linear_continuum_eqn( spec, mean, lowerEnengy, upperEnergy, 4, 4 );
    cont->setPolynomialCoefFitFor( 0, true );
    cont->setPolynomialCoefFitFor( 1, true );
    
    inputPeaks.push_back( std::move(peak) );
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  ROOT::Minuit2::MnUserParameters inputFitPars;
  PeakFitChi2Fcn::addPeaksToFitter( inputFitPars, inputPeaks, spec, PeakFitChi2Fcn::kFitForPeakParameters );
  
      
  const int npeaks = static_cast<int>( inputPeaks.size() );
  PeakFitChi2Fcn chi2Fcn( npeaks, spec, nullptr );
  chi2Fcn.useReducedChi2( false );
      
  assert( inputFitPars.VariableParameters() != 0 );
          
  ROOT::Minuit2::MnUserParameterState inputParamState( inputFitPars );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2Fcn, inputParamState, strategy );
      
  unsigned int maxFcnCall = 5000;
  double tolerance = 2.5;
  tolerance = 0.5;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  const ROOT::Minuit2::MnUserParameters &fitParams = fitter.Parameters();
  //  minimum.IsValid()
  //      ROOT::Minuit2::MinimumState minState = minimum.State();
  //      ROOT::Minuit2::MinimumParameters minParams = minState.Parameters();
    
  //    cerr << endl << endl << "EDM=" << minimum.Edm() << endl;
  //    cerr << "MinValue=" <<  minimum.Fval() << endl << endl;
      
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
        
  if( !minimum.IsValid() )
  {
      //XXX - should we try to re-fit here? Or do something to handle the
      //      faliure in some reasonable way?
      cerr << endl << endl << "status is not valid"
            << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
            << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
            << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
            << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
            << "\n\tHasValidParameters: " << minimum.HasValidParameters()
            << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
            << endl;
      if( minimum.IsAboveMaxEdm() )
        cout << "\t\tEDM=" << minimum.Edm() << endl;
  }//if( !minimum.IsValid() )
      
  
  vector<double> fitpars = fitParams.Params();
  vector<double> fiterrors = fitParams.Errors();
  chi2Fcn.parametersToPeaks( fitPeaks, &fitpars[0], &fiterrors[0] );
      
  double initialChi2 = chi2Fcn.chi2( &fitpars[0] );
  
  //Lets try to keep whether or not to fit parameters should be the same for
  //  the output peaks as the input peaks.
  //Note that this doesnt account for peaks swapping with eachother in the fit
  assert( fitPeaks.size() == inputPeaks.size() );
  
  //for( size_t i = 0; i < near_peaks.size(); ++i )
  //  fitpeaks[i].inheritUserSelectedOptions( near_peaks[i], true );
  //for( size_t i = 0; i < fixedpeaks.size(); ++i )
  //  fitpeaks[i+near_peaks.size()].inheritUserSelectedOptions( fixedpeaks[i], true );
        
  set_chi2_dof( spec, fitPeaks, 0, fitPeaks.size() );
      
  
  auto label = new WLabel( "Initial Chi2=" + std::to_string(initialChi2), this );
  cout << "Intiial Chi2=" << initialChi2 << endl;
  label->setInline( false );
  
  using boost::math::chi_squared;
  using boost::math::quantile;
  using boost::math::complement;
  using boost::math::cdf;
  
  
  std::sort( fitPeaks.begin(), fitPeaks.end(), &PeakDef::lessThanByMean );
  for( auto &peak : fitPeaks )
  {
    const double amp = peak.amplitude();
    string msg = "Peak at " + std::to_string(peak.mean()) + " keV"
                 + " fit amplitude: " + std::to_string(amp)
                 + " and Chi2/dof=" + std::to_string(peak.chi2dof());
    auto label = new WLabel( msg, this );
    label->setInline( false );
    cout << msg << endl;
    
    std::vector<PeakDef *> peakptrs;
    for( auto &peak : fitPeaks )
      peakptrs.push_back( &peak );
      
    double chi2, dof;
    get_chi2_and_dof_for_roi( chi2, dof, spec, peakptrs );
    
    chi_squared dist( dof );
    const double prob = boost::math::cdf(dist,chi2); //Probability we would have seen a chi2 this large.
    double p_value = 1.0 - prob; //Probability we would have observed this good of a chi2, or better
    cout << "There are " << dof << " DOF" << endl;
    cout << "This chi2=" << chi2 << " and prob of observing a Chi2 at least this extreme " << p_value << endl;
    cout << "95% of the time, we would see a Chi2 value of " << quantile(dist, 0.95) << " or less" << endl;
    
    for( double alpha = 0.0; alpha < 1.0; alpha += 0.05 )
    {
      cout << "alpha=" << alpha << " --> quantile=" << quantile(dist, alpha) << endl;
    }
    
    //We are estimating 1 parameter (amplitude), so the amplitude of the peak that we are 95% confident
    //  it would be lower than is the amplitude that causes a chi2 change of 3.84
    //const double n_parameters_being_fit = 1;
    //chi_squared deltadist( n_parameters_being_fit );
    //cout << "quantile(0.95)=" << quantile(deltadist, 0.95) << " (should be 3.84)" << endl;
    
    //We wa
    
  }//for( const auto &peak : fitPeaks )
  
  
  
  
  
  
  
  /*
  vector<PeakDef> *all_peaks = &fitpeaks;
      std::unique_ptr< vector<PeakDef> > all_peaks_ptr;
      if( fixedpeaks.size() )
      {
        all_peaks = new vector<PeakDef>( fitpeaks );
        all_peaks_ptr.reset( all_peaks );
        all_peaks->insert( all_peaks->end(), fixedpeaks.begin(), fixedpeaks.end() );
        std::sort( all_peaks->begin(), all_peaks->end(), &PeakDef::lessThanByMean );
      }//if( fixedpeaks.size() )
      
      for( size_t i = 1; i <= fitpeaks.size(); ++i ) //Note weird convntion of index
      {
        const PeakDef *peak = &(fitpeaks[i-1]);
        const bool significant = chi2_significance_test(
                                                            *peak, stat_threshold, hypothesis_threshold,
                                                            *all_peaks, data );
        if( !significant )
        {
  #if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          DebugLog(cerr) << "\tPeak at mean=" << peak->mean()
                         << "is being discarded for not being significant"
                         << "\n";
  #endif
          fitpeaks.erase( fitpeaks.begin() + --i );
        }//if( !significant )
      }//for( size_t i = 1; i < fitpeaks.size(); ++i )
          
      bool removed_peak = false;
      for( size_t i = 1; i < fitpeaks.size(); ++i ) //Note weird convntion of index
      {
        PeakDef *this_peak = &(fitpeaks[i-1]);
        PeakDef *next_peak = &(fitpeaks[i+1-1]);
            
        const double min_sigma = min( this_peak->sigma(), next_peak->sigma() );
        const double mean_diff = next_peak->mean() - this_peak->mean();
            
        //In order to remove a gaussian, the peaks must both be within a sigma
        //  of eachother.  Note that this proccess doesnt care about the widths
        //  of the gaussians because we are assuming that the width of the gaussian
        //  should only be dependant on energy, so should only have one width of
        //  gaussian for a given energy
        if( (mean_diff/min_sigma) < 1.0 ) //XXX 1.0 chosen arbitrarily, and not checked
        {
  #if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          DebugLog(cerr) << "Removing duplicate peak at x=" << this_peak->mean() << " sigma="
              << this_peak->sigma() << " in favor of mean=" << next_peak->mean()
              << " sigma=" << next_peak->sigma() << "\n";
  #endif
              
          removed_peak = true;
              
          //Delete the peak with the worst chi2
          if( this_peak->chi2dof() > next_peak->chi2dof() )
            fitpeaks.erase( fitpeaks.begin() + i - 1 );
          else
            fitpeaks.erase( fitpeaks.begin() + i );
              
          i = i - 1; //incase we have multiple close peaks in a row
        }//if( (mean_diff/min_sigma) < 1.0 ) / else
      }//for( size_t i = 1; i < fitpeaks.size(); ++i )
          
      if( removed_peak )
        fitPeaks( fitpeaks, stat_threshold, hypothesis_threshold,
                  data, fitpeaks, fixedpeaks, false );
          
      if( datadefined_peaks.size() )
      {
        fitpeaks.insert( fitpeaks.end(),
                        datadefined_peaks.begin(), datadefined_peaks.end() );
        std::sort( fitpeaks.begin(), fitpeaks.end(), &PeakDef::lessThanByMean );
      }
          
      return;
    }catch( std::exception &e )
    {
      cerr << "fitPeaks(...)\n\tSerious programming logic error: caught"
           << " exception where I really shouldnt have.  what()=" << e.what()
           << endl;
    }catch(...)
    {
      cerr << "fitPeaks(...)\n\tSerious programming logic error: caught"
           << " unknown exception where I really shouldnt have." << endl;
    }//try/catch()
          
    //We will only reach here if there was no exception, so since never expect
    //  this to actually happen, just assign the results to be same as the input
    fitpeaks = all_near_peaks;
  }//vector<PeakDef> fitPeaks(...);
   */
}//void do_development()
