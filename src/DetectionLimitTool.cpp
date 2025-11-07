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
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UserPreferences.h"
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
#include "InterSpec/PhysicalUnitsLocalized.h"
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
  
  return !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec );
}//bool use_curie_units()

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
                                 : m_input.counts_per_bq_into_4pi__);
      const double gammas_per_bq = counts_4pi * det_eff;
      
      const DetectionLimitCalc::CurrieMdaInput input = currieInput();
      
      const DetectionLimitCalc::CurrieMdaResult result = DetectionLimitCalc::currie_mda_calc( input );
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
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            const string upperstr = PhysicalUnits::printToBestActivityUnits( upper_act, 2, useCuries )
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            const string nomstr = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, useCuries )
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            
            //cout << "At " << m_input.energy << "keV observed " << result.source_counts
            //     << " counts, at distance " << PhysicalUnits::printToBestLengthUnits(m_input.distance)
            //     << ", leading to nominal activity " << nomstr
            //     << endl;
            
            m_poisonLimit->setText( "<table><tr><td style=\"padding-right:5px\">Nominal:</td><td>" + nomstr + "</td></tr>"
                                   "<tr><td>Range:</td><td>[" + lowerstr + ", " + upperstr + "]</td></tr></table>" );
            m_poisonLimit->setToolTip( "Detected activity, using just this Region Of Interests,"
                                      " and the observed excess of counts, as well as the"
                                      " statistical confidence interval." );
          }else if( result.upper_limit < 0 )
          {
            // This can happen when there are a lot fewer counts in the peak region than predicted
            //  from the sides - since this is non-sensical, we'll just say zero.
            const string unitstr = (useCuries ? "Ci" : "Bq") + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            m_poisonLimit->setText( "<div>Limit: &lt; 0" + unitstr + "</div><div>(sig. fewer counts in ROI than predicted)</div>" );
            m_poisonLimit->setToolTip( "Significantly fewer counts were observed in the"
                                      " Region Of Interest, than predicted by the neighboring channels." );
          }else
          {
            // We will provide the upper bound on activity.
            const double simple_mda = result.upper_limit / gammas_per_bq;
            const string mdastr = PhysicalUnits::printToBestActivityUnits( simple_mda, 2, useCuries )
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            
            const double det_limit = result.detection_limit / gammas_per_bq;
            const string det_limit_str = PhysicalUnits::printToBestActivityUnits( det_limit, 2, useCuries )
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
           
            string txt = "<table><tr><td style=\"padding-right: 5px\">Upper Limit:</td><td>" + mdastr + "</td></tr>"
            "<tr><td style=\"padding-right: 5px\">Critical Limit:</td><td>" + det_limit_str + "</td></tr></table>";
            
            m_poisonLimit->setText( txt );
            m_poisonLimit->setToolTip( "The upper limit is the maximum activity present, for the given CL,"
                                      " given the observed number of counts in the ROI.\n"
                                      "The critical limit is the activity where the signal would reliably be detected.\n"
                                      "Both numbers calculated using the Currie method."
                                      );
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
            const double counts_4pi_no_air = m_input.counts_per_bq_into_4pi__;
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
            
            //char buffer[256];
            //snprintf( buffer, sizeof(buffer),
            //         "<div>Nominal distance %s</div><div>%.1f%% CL range [%s, %s]</div>",
            //         nomstr.c_str(), rnd_cl_percent, upperstr.c_str(), lowerstr.c_str() );
            char buffer[384];
            snprintf( buffer, sizeof(buffer),
                     "<table><tr><td style=\"padding-right:5px\">Nominal:</td><td>%s</td></tr><tr><td>Range:</td><td>[%s, %s]</td></tr></table>",
                     nomstr.c_str(), upperstr.c_str(), lowerstr.c_str() );
            
            m_poisonLimit->setText( buffer );
          }else if( result.upper_limit < 0 )
          {
            // This can happen when there are a lot fewer counts in the peak region than predicted
            //  from the sides - since this is non-sensical ...
            m_poisonLimit->setText( "Large deficit of counts in ROI" );
          }else
          {
            // We will provide a "Distance >=X meters at 95% CL" answer
            //
            // result.upper_limit, is the number of expected signal counts that we can be 95%
            //  certain the true signal is less than.
            // So we will find the distance where we would expect this many counts
            //  (result.upper_limit), from the source at the given activity, and this is the
            //  distance we can be 95% sure the detector is that distance, or farther from the
            //  source
            
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
            //snprintf( buffer, sizeof(buffer), "Distance ≥%s @%.1f%% CL", upperstr.c_str(), rnd_cl_percent );
            snprintf( buffer, sizeof(buffer), "Distance ≥%s", upperstr.c_str() );
            
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
                               : input.counts_per_bq_into_4pi__);
    
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
    const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
    
    
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
    
    const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
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
    
    WText *currie_label = new WText( "Single peak Currie Limit", currieLimitContent );
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
    //viewer->preferences()->addCallbackWhenChanged( "DisplayBecquerel",
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
  
  DetectionLimitCalc::CurrieMdaInput currieInput() const
  {
    auto m = m_input.measurement;
    if( !m || (m->num_gamma_channels() < 16) )
      throw runtime_error( "No measurement." );
    
    const float roi_lower_energy = m_roi_start->value();
    const float roi_upper_energy = m_roi_end->value();
    
    const size_t nsidebin = static_cast<size_t>( std::round( std::max( 1.0f, m_num_side_channels->value() ) ) );
    const size_t nchannels = m->num_gamma_channels();
    
    DetectionLimitCalc::CurrieMdaInput input;
    input.spectrum = m;
    input.gamma_energy = m_input.energy;
    input.roi_lower_energy = m_roi_start->value();
    input.roi_upper_energy = m_roi_end->value();
    input.num_lower_side_channels = nsidebin;
    input.num_upper_side_channels = nsidebin;
    input.detection_probability = m_input.confidence_level;
    input.additional_uncertainty = 0.0f;  // TODO: can we get the DRFs contribution to form this?
    
    return input;
  }//DetectionLimitCalc::CurrieMdaInput currieInput() const
  
  void createMoreInfoWindow()
  {
    try
    {
      if( !m_input.measurement )
        throw runtime_error( "No measurement" );
      
      if( !m_input.drf )
        throw runtime_error( "No detector efficiency function" );
      
      const DetectionLimitCalc::CurrieMdaInput input = currieInput();
      const DetectionLimitCalc::CurrieMdaResult result = DetectionLimitCalc::currie_mda_calc( input );
      
      
      DetectionLimitTool::createCurrieRoiMoreInfoWindow( m_nuclide,
                                 result,
                                 m_input.drf,
                                 m_input.limit_type,
                                 m_input.distance,
                                 m_input.do_air_attenuation,                           
                                 m_input.branch_ratio,
                                 m_input.shield_transmission );
    }catch( std::exception &e )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "Error computing %s, %.2f keV Info",
               (m_nuclide ? m_nuclide->symbol.c_str() : "null"), m_input.energy );
      
      SimpleDialog *dialog = new SimpleDialog( buffer,
                                              WString("Error computing Currie limit information: {1}").arg(e.what()) );
      dialog->addButton( "Close" );
    }//try / catch
  }//void createMoreInfoWindow()
  
};//class MdaPeakRow




DetectionLimitWindow::DetectionLimitWindow( InterSpec *viewer,
                                                     MaterialDB *materialDB,
                                                     WSuggestionPopup *materialSuggest )
: AuxWindow( "Detection Confidence Tool",
  (AuxWindowProperties::TabletNotFullScreen
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
}//DetectionLimitWindow(...) constructor


DetectionLimitWindow::~DetectionLimitWindow()
{
}


DetectionLimitTool *DetectionLimitWindow::tool()
{
  return m_tool;
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
  m_chart->disableLegend();
  m_chart->showHistogramIntegralsInLegend( true );

  //m_chart->xAxisSliderShown().connect(...)
  
  m_peakModel = new PeakModel( this );
  m_chart->setPeakModel( m_peakModel );

  // Create the user-input area under the spectrum and liklihood chart
  WContainerWidget *inputTable = new WContainerWidget( this );
  inputTable->addStyleClass( "Inputs" );
  
  
  // Create nuclide label and input
  WLabel *label = new WLabel( "Nuclide:&nbsp;&nbsp;", inputTable ); //The space so this will be the longest label, and not "Distance:"
  label->addStyleClass( "GridFirstCol GridFirstRow GridVertCenter" );
  
  
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
  label->addStyleClass( "GridFirstCol GridSecondRow GridVertCenter" );
  
  m_ageEdit = new WLineEdit( "", inputTable );
  m_ageEdit->addStyleClass( "GridSecondCol GridSecondRow" );
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), m_ageEdit );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_ageEdit->setValidator(validator);
  m_ageEdit->setAutoComplete( false );
  label->setBuddy( m_ageEdit );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  const char *tooltip =
  "<div>The age of the nuclide.</div>"
  "<br />"
  "<div>The age controls the amount of progeny in-growth; the activities of the parent"
  " (i.e., entered) nuclide are always reported for at the time of measurement."
  ".</div>";
  HelpSystem::attachToolTipOn( {label, m_ageEdit},
                              tooltip, showToolTips, HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::RespectPreference );
  
  m_ageEdit->changed().connect( this, &DetectionLimitTool::handleUserAgeChange );
  m_ageEdit->blurred().connect( this, &DetectionLimitTool::handleUserAgeChange );
  m_ageEdit->enterPressed().connect( this, &DetectionLimitTool::handleUserAgeChange );
  
  
  
  // Create distance input
  label = new WLabel( "Distance:", inputTable );
  label->addStyleClass( "GridFirstCol GridThirdRow GridVertCenter" );
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
  label->addStyleClass( "GridFirstCol GridThirdRow GridVertCenter" );
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
  m_chart->yAxisLogLinChanged().connect( boost::bind( &SwitchCheckbox::setUnChecked, loglin, boost::placeholders::_1 ) );
  
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
  label->addStyleClass( "GridSixthCol GridFirstRow GridVertCenter" );
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
  label->addStyleClass( "GridSixthCol GridFirstRow GridVertCenter" );
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
  label->addStyleClass( "GridSixthCol GridSecondRow GridVertCenter" );
  m_confidenceLevel = new WComboBox( inputTable );
  m_confidenceLevel->addStyleClass( "GridSeventhCol GridSecondRow GridVertCenter" );
  
  for( auto cl = ConfidenceLevel(0); cl < NumConfidenceLevel; cl = ConfidenceLevel(cl+1) )
  {
    const char *txt = "";
    
    switch( cl )
    {
      case ConfidenceLevel::NinetyFivePercent: txt = "95%";        break;
      case ConfidenceLevel::NinetyNinePercent: txt = "99%";        break;
      case ConfidenceLevel::OneSigma:          txt = "1σ (68.2%)"; break;
      case ConfidenceLevel::TwoSigma:          txt = "2σ (95.4%)"; break;
      case ConfidenceLevel::ThreeSigma:        txt = "3σ (99.7%)"; break;
      case ConfidenceLevel::FourSigma:         txt = "4σ";         break;
      case ConfidenceLevel::FiveSigma:         txt = "5σ";         break;
      case ConfidenceLevel::NumConfidenceLevel: break;
    }//switch( cl )
    
    m_confidenceLevel->addItem( txt );
  }//for( loop over confidence levels )
  
  m_confidenceLevel->setCurrentIndex( ConfidenceLevel::NinetyFivePercent );
  
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
    m_peakModel->setPeakFromSpecMeas( m_our_meas, {ourspec->sample_number()}, SpecUtils::SpectrumType::Foreground );
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
  viewer->preferences()->addCallbackWhenChanged( "DisplayBecquerel",
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
      log_developer_error( __func__, msg );
      assert( 0 );
    }
  }
#endif  //PERFORM_DEVELOPER_CHECKS
}//DetectionLimitTool constructor
  

  
DetectionLimitTool::~DetectionLimitTool()
{
}



void DetectionLimitTool::update_spectrum_for_currie_result( D3SpectrumDisplayDiv *chart,
                                       PeakModel *pmodel,
                                       const DetectionLimitCalc::CurrieMdaInput &input,
                                       const DetectionLimitCalc::CurrieMdaResult * const result,
                                       std::shared_ptr<const DetectorPeakResponse> drf,
                                       DetectionLimitTool::LimitType limitType,
                                       const double gammas_per_bq,
                                       const vector<DetectionLimitTool::CurrieResultPeak> &peaks )
{
  assert( chart );
  if( !chart )
    return;
  
  chart->removeAllDecorativeHighlightRegions();
  if( pmodel )
    pmodel->setPeaks( vector<shared_ptr<const PeakDef>>{} );
  
  shared_ptr<const SpecUtils::Measurement> spectrum = input.spectrum;
  if( !spectrum )
    spectrum = chart->data();
  
  assert( spectrum );
  if( !spectrum )
    return;
  
  InterSpec *viewer = InterSpec::instance();
  shared_ptr<const ColorTheme> theme = viewer ? viewer->getColorTheme() : nullptr;
  assert( theme );
  if( !theme )
    return;
  
  const bool useCuries = use_curie_units();
  
  const double confidence_level = input.detection_probability;
  
  // We should only get an exception if the ROI goes of the spectrum, or invalid energy cal or something
  try
  {
    double lower_continuum_counts_sum, peak_region_counts_sum, upper_continuum_counts_sum;
    double lower_lower_energy, lower_upper_energy, upper_lower_energy, upper_upper_energy;
    
    // If number of of side channels is zero, then we expect the following
    //  variables to all be zero
    assert( !result
           || (result->input.num_lower_side_channels != 0)
           || ((result->first_lower_continuum_channel == 0)
               && (result->last_lower_continuum_channel == 0)
               && (result->lower_continuum_counts_sum == 0)
               && (result->first_upper_continuum_channel == 0)
               && (result->last_upper_continuum_channel == 0)
               && (result->upper_continuum_counts_sum == 0))
           );
    
    const bool assertedNoSignal = (input.num_lower_side_channels == 0);
    
    // We will prefer to set the highlight regions from the results, but if we dont have results
    //  we'll do it from the input.
    if( result )
    {
      if( assertedNoSignal )
      {
        lower_upper_energy = spectrum->gamma_channel_upper( result->last_peak_region_channel );
        upper_upper_energy = lower_upper_energy;
        upper_lower_energy = spectrum->gamma_channel_lower( result->first_peak_region_channel );
        lower_lower_energy = upper_lower_energy;
        
        lower_continuum_counts_sum = 0;
        peak_region_counts_sum = result->peak_region_counts_sum;
        upper_continuum_counts_sum = 0;
      }else
      {
        lower_lower_energy = spectrum->gamma_channel_lower( result->first_lower_continuum_channel );
        lower_upper_energy = spectrum->gamma_channel_lower( result->first_upper_continuum_channel ); //need
        upper_lower_energy = spectrum->gamma_channel_upper( result->last_lower_continuum_channel );  //need
        upper_upper_energy = spectrum->gamma_channel_upper( result->last_upper_continuum_channel );
        
        lower_continuum_counts_sum = result->lower_continuum_counts_sum;
        peak_region_counts_sum = result->peak_region_counts_sum;
        upper_continuum_counts_sum = result->upper_continuum_counts_sum;
      }
    }else
    {
      const pair<size_t,size_t> channels 
            = DetectionLimitCalc::round_roi_to_channels( spectrum, input.roi_lower_energy, input.roi_upper_energy );
      
      const size_t first_peak_region_channel = channels.first;
      const size_t last_peak_region_channel = channels.second;
      
      if( assertedNoSignal )
      {
        lower_upper_energy = spectrum->gamma_channel_upper( last_peak_region_channel );
        upper_upper_energy = lower_upper_energy;
        upper_lower_energy = spectrum->gamma_channel_lower( first_peak_region_channel );
        lower_lower_energy = upper_lower_energy;
        
        lower_continuum_counts_sum = 0;
        peak_region_counts_sum = spectrum->gamma_channels_sum(first_peak_region_channel, last_peak_region_channel);
        upper_continuum_counts_sum = 0;
      }else
      {
        if( first_peak_region_channel < (input.num_lower_side_channels + 1) )
          throw std::runtime_error( "mda_counts_calc: lower peak region is outside spectrum energy range" );
        
        const size_t last_lower_continuum_channel = first_peak_region_channel - 1;
        const size_t first_lower_continuum_channel = last_lower_continuum_channel - input.num_lower_side_channels + 1;
        
        const size_t first_upper_continuum_channel = last_peak_region_channel + 1;
        const size_t last_upper_continuum_channel = first_upper_continuum_channel + input.num_upper_side_channels - 1;
        
        if( last_upper_continuum_channel >= spectrum->num_gamma_channels() )
          throw std::runtime_error( "mda_counts_calc: upper peak region is outside spectrum energy range" );
        
        lower_lower_energy = spectrum->gamma_channel_lower( first_lower_continuum_channel );
        lower_upper_energy = spectrum->gamma_channel_lower( first_upper_continuum_channel ); //need
        upper_lower_energy = spectrum->gamma_channel_upper( last_lower_continuum_channel ); //need
        upper_upper_energy = spectrum->gamma_channel_upper( last_upper_continuum_channel );
        
        lower_continuum_counts_sum = spectrum->gamma_channels_sum(first_lower_continuum_channel, last_lower_continuum_channel);
        peak_region_counts_sum =     spectrum->gamma_channels_sum(first_peak_region_channel, last_peak_region_channel);
        upper_continuum_counts_sum = spectrum->gamma_channels_sum(first_upper_continuum_channel, last_upper_continuum_channel);
      }//if( assertedNoSignal )
    }//if( result ) / else
    
    assert( !assertedNoSignal
           || ((lower_continuum_counts_sum == 0)
               && (upper_continuum_counts_sum == 0))
           );
    
    const double dx = upper_upper_energy - lower_lower_energy;
    chart->setXAxisRange( lower_lower_energy - 0.5*dx, upper_upper_energy + 0.5*dx );
    
    const int lower_ndec = 1 + static_cast<int>( std::ceil( fabs( std::log10( fabs(lower_continuum_counts_sum) ) ) ) );
    const int mid_ndec = 1 + static_cast<int>( std::ceil( fabs( std::log10( fabs(peak_region_counts_sum) ) ) ) );
    const int upper_ndec = 1 + static_cast<int>( std::ceil( fabs( std::log10( fabs(upper_continuum_counts_sum) ) ) ) );
    
    shared_ptr<const ColorTheme> theme = viewer->getColorTheme();
    assert( theme );
    if( !theme )
      throw runtime_error( "Invalid color theme" );
    
    if( !assertedNoSignal )
    {
      const string lower_txt = SpecUtils::printCompact( lower_continuum_counts_sum, lower_ndec );
      chart->addDecorativeHighlightRegion( lower_lower_energy, upper_lower_energy,
                                          theme->timeHistoryBackgroundHighlight,
                                          D3SpectrumDisplayDiv::HighlightRegionFill::BelowData,
                                          lower_txt );
    }//if( !assertedNoSignal )
    
    const string mid_txt = SpecUtils::printCompact( peak_region_counts_sum, mid_ndec );
    chart->addDecorativeHighlightRegion( upper_lower_energy, lower_upper_energy,
                                        theme->timeHistoryForegroundHighlight,
                                        D3SpectrumDisplayDiv::HighlightRegionFill::BelowData,
                                        mid_txt );
    
    if( !assertedNoSignal )
    {
      const string upper_txt = SpecUtils::printCompact( upper_continuum_counts_sum, upper_ndec );
      chart->addDecorativeHighlightRegion( lower_upper_energy, upper_upper_energy,
                                          theme->timeHistoryBackgroundHighlight,
                                          D3SpectrumDisplayDiv::HighlightRegionFill::BelowData,
                                          upper_txt );
    }//if( !assertedNoSignal )
    
    // We will only put the peak on the chart if there is a result
    if( result && pmodel )
    {
      const double peak_sigma = (drf && drf->hasResolutionInfo())
                                ? drf->peakResolutionSigma(input.gamma_energy)
                                : 0.25*(lower_upper_energy - upper_lower_energy);
      
      const DetectorPeakResponse::EffGeometryType det_geom = drf ? drf->geometryType() : DetectorPeakResponse::EffGeometryType::FarField;
      
      PeakDef generic_peak( input.gamma_energy, peak_sigma, 0.0 );
      generic_peak.continuum()->setRange( upper_lower_energy, lower_upper_energy );
      generic_peak.continuum()->setType(PeakContinuum::OffsetType::Linear);
      generic_peak.continuum()->setParameters( input.gamma_energy, {result->continuum_eqn[0], result->continuum_eqn[1]}, {} );
      
      double sum_counts_4pi = 0.0;
      vector<PeakDef> specific_peaks( peaks.size() );
      for( size_t i = 0; i < peaks.size(); ++i )
      {
        sum_counts_4pi += peaks[i].counts_4pi;
        specific_peaks[i].setMean( peaks[i].energy );
        specific_peaks[i].setSigma( peaks[i].fwhm/2.35482 );
        specific_peaks[i].setContinuum( generic_peak.continuum() );
      }
      
      char cl_str_buffer[64] = {'\0'};
      
      if( confidence_level < 0.999 )
        snprintf( cl_str_buffer, sizeof(cl_str_buffer), "%.1f%%", 100.0*confidence_level );
      else
        snprintf( cl_str_buffer, sizeof(cl_str_buffer), "1-%.2G", (1.0-confidence_level) );
     
      WString chart_title;
      const string cl_str = cl_str_buffer;
      
      switch( limitType )
      {
        case DetectionLimitTool::LimitType::Activity:
        {
          if( result->source_counts > result->decision_threshold )
          {
            assert( !assertedNoSignal );
            
            // There is enough excess counts that we would reliably detect this activity, so we will
            //  give the activity range.
            string lowerstr, upperstr, nomstr;
            
            if( gammas_per_bq > 0.0 )
            {
              const float lower_act = result->lower_limit / gammas_per_bq;
              const float upper_act = result->upper_limit / gammas_per_bq;
              const float nominal_act = result->source_counts / gammas_per_bq;
              
              lowerstr = PhysicalUnits::printToBestActivityUnits( lower_act, 2, useCuries )
              + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
              upperstr = PhysicalUnits::printToBestActivityUnits( upper_act, 2, useCuries )
              + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
              nomstr = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, useCuries )
              + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
              
              chart_title = "Estimated activity of " + nomstr + ".";
            }else
            {
              lowerstr = SpecUtils::printCompact(result->lower_limit, 4);
              upperstr = SpecUtils::printCompact(result->upper_limit, 4);
              nomstr = SpecUtils::printCompact(result->source_counts, 4);
              
              chart_title = "Excess counts of " + nomstr + ".";
            }
            
            //const string cl_str = SpecUtils::printCompact( 100.0*confidence_level, 3 );
            
            generic_peak.setPeakArea( result->source_counts );
            
            for( size_t i = 0; i < peaks.size(); ++i )
            {
              const double area = result->source_counts * peaks[i].counts_4pi / sum_counts_4pi;
              specific_peaks[i].setAmplitude( area );
            }
          }else if( result->upper_limit < 0 )
          {
            assert( !assertedNoSignal );
            
            // This can happen when there are a lot fewer counts in the peak region than predicted
            //  from the sides - since this is non-sensical, we'll just say zero.
            const string unitstr = useCuries ? "Ci" : "Bq";
            
            chart_title = "Activity < 0 " + unitstr;
            
            generic_peak.setPeakArea( 0.0 );
            for( size_t i = 0; i < peaks.size(); ++i )
              specific_peaks[i].setAmplitude( 0.0 );
          }else
          {
            if( assertedNoSignal )
            {
              // We will provide minimum counts reliably detectable
              string mdastr;
              if( gammas_per_bq > 0.0 )
              {
                const double simple_mda = result->detection_limit / gammas_per_bq;
                mdastr = PhysicalUnits::printToBestActivityUnits( simple_mda, 2, useCuries )
                      + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
              }else
              {
                mdastr = SpecUtils::printCompact( result->detection_limit, 4 ) + " counts";
              }
              
              chart_title = "Data + peak for " + mdastr;
              
              generic_peak.setPeakArea( result->detection_limit );
              for( size_t i = 0; i < peaks.size(); ++i )
              {
                const double area = result->detection_limit * peaks[i].counts_4pi / sum_counts_4pi;
                specific_peaks[i].setAmplitude( area );
              }
            }else
            {
              // We will provide the upper bound on activity.
              string mdastr;
              if( gammas_per_bq > 0.0 )
              {
                const double simple_mda = result->upper_limit / gammas_per_bq;
                mdastr = PhysicalUnits::printToBestActivityUnits( simple_mda, 2, useCuries )
                      + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
              }else
              {
                mdastr = SpecUtils::printCompact( result->upper_limit, 4 ) + " counts";
              }
              
              chart_title = "Peak for upper bound of " + mdastr + " @" + cl_str + " CL";
              
              generic_peak.setPeakArea( result->upper_limit );
              for( size_t i = 0; i < peaks.size(); ++i )
              {
                const double area = result->upper_limit * peaks[i].counts_4pi / sum_counts_4pi;
                specific_peaks[i].setAmplitude( area );
              }
            }//if( assertedNoSignal ) / else
          }//if( detected ) / else if( ....)
          
          break;
        }//case DetectionLimitTool::LimitType::Activity:
          
        case DetectionLimitTool::LimitType::Distance:
        {
          // Distance not supported yet
          break;
        }
      }//switch( limitType )
      
      chart->setChartTitle( chart_title );
      
      vector<PeakDef> peaks = (specific_peaks.size() > 0) ? std::move(specific_peaks) : vector<PeakDef>{generic_peak};
      
      if( assertedNoSignal )
      {
        for( PeakDef &p : peaks )
        {
          p.continuum()->setType( PeakContinuum::OffsetType::External );
          p.continuum()->setExternalContinuum( spectrum );
        }
      }//if( specIsBackground )
      
      set_chi2_dof( spectrum, peaks, 0, peaks.size() ); // Compute Chi2 for peaks
      
      pmodel->addPeaks( peaks );
    }//if( result )
  }catch( std::exception &e )
  {
    cerr << "update_spectrum_for_currie_result, Caught exception:" << e.what() << endl;
  }//try / catch
}//void update_spectrum_for_currie_result( D3SpectrumDisplayDiv *chart, PeakModel *pmodel )


SimpleDialog *DetectionLimitTool::createCurrieRoiMoreInfoWindow( const SandiaDecay::Nuclide *const nuclide,
                                const DetectionLimitCalc::CurrieMdaResult &result,
                                std::shared_ptr<const DetectorPeakResponse> drf,
                                DetectionLimitTool::LimitType limitType,
                                const double distance,
                                const bool do_air_attenuation,
                                const double branch_ratio,
                                double shield_transmission )
{
  const DetectionLimitCalc::CurrieMdaInput &input = result.input;
  const float &energy = input.gamma_energy;
  
  assert( (shield_transmission > 0.0) && (shield_transmission <= 1.0) );
  if( (shield_transmission <= 0.0) || (shield_transmission > 1.0) )
    shield_transmission = 1.0;
  
  char buffer[256];
  snprintf( buffer, sizeof(buffer), "%s%.2f keV Info",
           (nuclide ? (nuclide->symbol + " ").c_str() : ""), energy );
  
  SimpleDialog *dialog = new SimpleDialog( buffer );
  dialog->addButton( "Close" );
  
  if( drf && !drf->isValid() )
    drf.reset();
  
  const float live_time = result.input.spectrum ? result.input.spectrum->live_time() : 1.0f;
  
  
  try
  {
    wApp->useStyleSheet( "InterSpec_resources/DetectionLimitTool.css" ); // JIC
    
    if( !input.spectrum )
      throw runtime_error( "No measurement" );
      
    const double conf_level = input.detection_probability;
    
    const bool useCuries = use_curie_units();
    const bool fixed_geom = drf ? drf->isFixedGeometry() : false;
    const DetectorPeakResponse::EffGeometryType det_geom = drf ? drf->geometryType()
                                                : DetectorPeakResponse::EffGeometryType::FarField;
    const bool air_atten = (do_air_attenuation && !fixed_geom && (distance > 0.0));
    const double intrinsic_eff = drf ? drf->intrinsicEfficiency( energy ) : 1.0f;
    const double geom_eff = (drf && (distance >= 0.0)) ? drf->fractionalSolidAngle( drf->detectorDiameter(), distance ) : 1.0;
    const double det_eff = fixed_geom ? intrinsic_eff : (drf ? drf->efficiency(energy, distance) : 1.0);
     
    const double gammas_per_bq_into_4pi = branch_ratio * live_time * shield_transmission;
    
    double gammas_per_bq_into_4pi_with_air = gammas_per_bq_into_4pi;
    if( distance > 0.0 )
    {
      const double air_atten_coef = GammaInteractionCalc::transmission_coefficient_air( energy, distance );
      gammas_per_bq_into_4pi_with_air = gammas_per_bq_into_4pi * exp( -1.0*air_atten_coef );
    }
    
    const double counts_4pi = air_atten ? gammas_per_bq_into_4pi_with_air : gammas_per_bq_into_4pi;
    
    const double gammas_per_bq = counts_4pi * det_eff;
    
    if( conf_level < 0.999 )
      snprintf( buffer, sizeof(buffer), "%.1f%%", 100.0*conf_level );
    else
      snprintf( buffer, sizeof(buffer), "1-%.3g", (1.0-conf_level) );
    
    const string confidence_level = buffer;
    
    const bool assertedNoSignal = (result.input.num_lower_side_channels == 0);
    
    const double lower_upper_energy = input.spectrum->gamma_channel_lower( result.first_peak_region_channel );
    const double lower_lower_energy = assertedNoSignal ? lower_upper_energy : input.spectrum->gamma_channel_lower( result.first_lower_continuum_channel );
    const double upper_lower_energy = input.spectrum->gamma_channel_upper( result.last_peak_region_channel );
    const double upper_upper_energy = assertedNoSignal ? upper_lower_energy : input.spectrum->gamma_channel_upper( result.last_upper_continuum_channel );
    
    // Add chart
    InterSpec *viewer = InterSpec::instance();
    assert( viewer );
    D3SpectrumDisplayDiv *chart = new D3SpectrumDisplayDiv( dialog->contents() );
    chart->addStyleClass( "MdaCurrieChart" );
    chart->setXAxisTitle( "" );
    chart->setYAxisTitle( "", "" );
    shared_ptr<const SpecUtils::Measurement> hist = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    chart->setData( hist, true );
    chart->setYAxisLog( false );
    chart->disableLegend();
    const double dx = upper_upper_energy - lower_lower_energy;
    chart->setXAxisRange( lower_lower_energy - 0.5*dx, upper_upper_energy + 0.5*dx );
    chart->setShowPeakLabel( SpectrumChart::PeakLabels::kShowPeakUserLabel, true );
    
    //TODO: set no interaction, and implement drawing the various areas...
    
    PeakModel *pmodel = new PeakModel( chart );
    pmodel->setNoSpecMeasBacking();
    chart->setPeakModel( pmodel );
    pmodel->setForeground( hist );
    
    vector<CurrieResultPeak> dummy_peaks; //We'll pass empty vector, so the function will create a single peak for us
    update_spectrum_for_currie_result( chart, pmodel, input, &result, drf, limitType, gammas_per_bq, dummy_peaks );
    chart->setChartTitle( "" );
    
    
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
                                  HelpSystem::ToolTipPrefOverride::InstantAlways );
    };//addTooltipToRow
    
    
    switch( limitType )
    {
      case DetectionLimitTool::LimitType::Activity:
      {
        if( result.source_counts > result.decision_threshold )
        {
          assert( !assertedNoSignal );
          
          // There is enough excess counts that we would reliably detect this activity, so we will
          //  give the activity range.
          if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
          {
            WString obs_label = "Observed activity";
            const double nominal_act = result.source_counts / gammas_per_bq;
            const string nomstr = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, useCuries )
                                  + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );

            cell = table->elementAt( table->rowCount(), 0 );
            new WText( obs_label, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( nomstr, cell );
            addTooltipToRow( "Greater than the &quot;critical level&quot;, L<sub>c</sub>,"
                            " counts were observed in the peak region." );
          }

          {
            WString obs_label = "Observed counts";
            const string nomstr = SpecUtils::printCompact(result.source_counts, 4);

            cell = table->elementAt( table->rowCount(), 0 );
            new WText( obs_label, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( nomstr, cell );
            addTooltipToRow( "Greater than the &quot;critical level&quot;, L<sub>c</sub>,"
                            " counts were observed in the peak region." );
          }


          if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
          {
            WString range_label = "Activity range";
            const double lower_act = result.lower_limit / gammas_per_bq;
            const double upper_act = result.upper_limit / gammas_per_bq;
            const string lowerstr = PhysicalUnits::printToBestActivityUnits( lower_act, 2, useCuries )
                        + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            const string upperstr = PhysicalUnits::printToBestActivityUnits( upper_act, 2, useCuries )
                        + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
            cell = table->elementAt( table->rowCount(), 0 );
            new WText( range_label, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( "[" + lowerstr + ", " + upperstr + "]", cell );
            addTooltipToRow( "The signal activity range estimate, to the " + confidence_level + " confidence level." );
          }

          {
            WString range_label = "Counts range";
            const string lowerstr = SpecUtils::printCompact(result.lower_limit, 4);
            const string upperstr = SpecUtils::printCompact(result.upper_limit, 4);
            cell = table->elementAt( table->rowCount(), 0 );
            new WText( range_label, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( "[" + lowerstr + ", " + upperstr + "]", cell );
            addTooltipToRow( "The signal counts range estimate, to the " + confidence_level + " confidence level." );
          }//if( drf ) / else
        }else if( result.upper_limit < 0 )
        {
          assert( !assertedNoSignal );
          
          // This can happen when there are a lot fewer counts in the peak region than predicted
          //  from the sides - since this is non-sensical, we'll just say zero.
          const string unitstr = drf ? (useCuries ? "Ci" : "Bq") : " counts";
          cell = table->elementAt( table->rowCount(), 1 );
          cell->setColumnSpan( 2 );
          if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
          {
            new WText( "Activity &le; 0 " + unitstr, cell );
          }else
          {
            new WText( "Counts &le; 0", cell );
          }
          
          addTooltipToRow( "Significantly fewer counts in peak region were observed,"
                          " than predicted by the neighboring regions." );
        }else
        {
          // We will provide the upper bound on activity.
          if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
          {
            WString label_txt = "Activity upper bound";

            const double simple_mda = result.upper_limit / gammas_per_bq;
            const string mdastr = PhysicalUnits::printToBestActivityUnits( simple_mda, 2, useCuries )
                      + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );

            cell = table->elementAt( table->rowCount(), 0 );
            new WText( label_txt, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( mdastr , cell);

            addTooltipToRow( "The upper limit on how much activity could be present, to the "
                            + confidence_level + " confidence level." );
          }


          {
            WString label_txt = "Counts upper bound";
            const string mdastr = SpecUtils::printCompact( result.upper_limit, 4 ) + " counts";

            cell = table->elementAt( table->rowCount(), 0 );
            new WText( label_txt, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( mdastr , cell);

            addTooltipToRow( "The upper limit on how much activity could be present, to the "
                            + confidence_level + " confidence level." );
          }
          

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
    }//switch( limitType )
    
    // Add a blank row
    string val;
    cell = table->elementAt( table->rowCount(), 0 );
    WText *txt = new WText( "&nbsp;", TextFormat::XHTMLText, cell );
    
    if( !assertedNoSignal )
    {
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Lower region channels", cell );
      val = "[" + std::to_string(result.first_lower_continuum_channel) + ", "
            + std::to_string(result.last_lower_continuum_channel) + "]";
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      
      
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
    }else
    {
      cell = table->elementAt( table->rowCount(), 0 );
      cell->setColumnSpan( 2 );
      txt = new WText( "Using ROI area as background estimate", cell );
    }//if( !assertedNoSignal ) / else
    
    cell = table->elementAt( table->rowCount(), 0 );
    txt = new WText( assertedNoSignal ? "ROI area channels" : "Peak area channels", cell );
    val = "[" + std::to_string(result.first_peak_region_channel) + ", "
                  + std::to_string(result.last_peak_region_channel) + "]";
    cell = table->elementAt( table->rowCount() - 1, 1 );
    txt = new WText( val, cell );
    
    const double peak_lower_energy = input.spectrum->gamma_channel_lower( result.first_peak_region_channel );
    const double peak_upper_energy = input.spectrum->gamma_channel_upper( result.last_peak_region_channel );
    snprintf( buffer, sizeof(buffer), "The region the peak is being assumed to be within;"
             " corresponds to %.2f to %.2f keV", peak_lower_energy, peak_upper_energy );
    addTooltipToRow( buffer );
    
    
    cell = table->elementAt( table->rowCount(), 0 );
    txt = new WText( assertedNoSignal ? "ROI region counts" : "Peak region counts", cell );
    val = SpecUtils::printCompact( result.peak_region_counts_sum, 5 );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    txt = new WText( val, cell );
    addTooltipToRow( "The observed number of counts in the peak region" );
    
    cell = table->elementAt( table->rowCount(), 0 );
    txt = new WText( assertedNoSignal ? "ROI region null est.": "Peak region null est.", cell );
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
    txt = new WText( "Peak critical limit <span style=\"font-size: smaller;\">(L<sub>c</sub>)</span>", cell );
    if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
    {
      const double decision_threshold_act = result.decision_threshold / gammas_per_bq;
      val = SpecUtils::printCompact( result.decision_threshold, 4 )
            + " <span style=\"font-size: smaller;\">("
            + PhysicalUnits::printToBestActivityUnits( decision_threshold_act, 2, useCuries )
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom )
            + ")</span>";
    }else
    {
      val = SpecUtils::printCompact( result.decision_threshold, 4 ) + " counts";
    }
    
    cell = table->elementAt( table->rowCount() - 1, 1 );
    txt = new WText( val, cell );
    addTooltipToRow( "Corresponds to Currie's &quot;critical level&quot;, L<sub>c</sub>,"
                    " that is the net signal level (instrument response) above which an"
                    " observed signal may be reliably recognized as &quot;detected&quot;." );
    
    
    // Note: I believe this quantity corresponds to Currie's "detection limit" (L_d) that
    //       is the “true” net signal level which may be a priori expected to lead to detection.
    cell = table->elementAt( table->rowCount(), 0 );
    txt = new WText( "Peak detection limit <span style=\"font-size: smaller;\">(L<sub>d</sub>)</span>", cell );
    
    if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
    {
      const double detection_limit_act = result.detection_limit / gammas_per_bq;
      val = SpecUtils::printCompact( result.detection_limit, 4 )
            + " <span style=\"font-size: smaller;\">("
            + PhysicalUnits::printToBestActivityUnits( detection_limit_act, 2, useCuries )
            + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom )
            + ")</span>";
    }else
    {
      val = SpecUtils::printCompact( result.detection_limit, 4 ) + " counts";
    }
    
    cell = table->elementAt( table->rowCount() - 1, 1 );
    txt = new WText( val, cell );
    addTooltipToRow( "Corresponds to Currie's &quot;detection limit&quot;, L<sub>d</sub>,"
                    " that is the &quot;true&quot; net signal level which may be, <i>a priori</i>."
                    " expected to lead to detection." );
    
    
    switch( limitType )
    {
      case DetectionLimitTool::LimitType::Activity:
      {
        if( result.source_counts <= result.decision_threshold )
        {
          // There is NOT enough excess counts that we would reliably detect this activity, so we
          //  didnt give nominal activity above, so we'll do that here
          if( !assertedNoSignal )
          {
            WString obs_counts_label = "Observed counts";
            string nom_counts_str = SpecUtils::printCompact(result.source_counts, 4);
            if( result.source_counts < 0 )
              nom_counts_str = "&lt; 0 counts";

            nom_counts_str += " <span style=\"font-size: smaller;\">(below L<sub>c</sub>)</span>";
            cell = table->elementAt( table->rowCount(), 0 );
            new WText( obs_counts_label, cell );
            cell = table->elementAt( table->rowCount() - 1, 1 );
            new WText( WString::fromUTF8(nom_counts_str), cell );
            addTooltipToRow( "The observed signal counts is less than the &quot;critical level&quot;, L<sub>c</sub>,"
                            " so a detection can not be declared, but this is the excess over expected counts." );

            if( drf && (distance >= 0.0) && (gammas_per_bq > 0.0) )
            {
              WString obs_act_label = "Activity";
              const float nominal_act = result.source_counts / gammas_per_bq;

              string nom_act_str = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, useCuries )
                        + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
              if( nominal_act < 0 )
                nom_act_str = "&lt; " + PhysicalUnits::printToBestActivityUnits( 0.0, 2, useCuries )
                          + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );

              nom_act_str += " <span style=\"font-size: smaller;\">(below L<sub>c</sub>)</span>";
              cell = table->elementAt( table->rowCount(), 0 );
              new WText( obs_act_label, cell );
              cell = table->elementAt( table->rowCount() - 1, 1 );
              new WText( WString::fromUTF8(nom_act_str), cell );
              addTooltipToRow( "The observed signal activity is less than the &quot;critical level&quot;, L<sub>c</sub>,"
                              " so a detection can not be declared, but this is the excess over expected counts." );
            }//if( !nom_act_str.empty() )
          }//if( !assertedNoSignal )
        }//if( result.source_counts <= result.decision_threshold )
      }//case DetectionLimitTool::LimitType::Activity:
        
      case DetectionLimitTool::LimitType::Distance:
        break;
    }//switch( limitType )
    
    
    // Add a blank row
    cell = table->elementAt( table->rowCount(), 0 );
    txt = new WText( "&nbsp;", TextFormat::XHTMLText, cell );
    if( drf )
    {
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Detector Intrinsic Eff.", cell );
      val = SpecUtils::printCompact( intrinsic_eff, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The efficiency for a gamma hitting the detector face,"
                      " to be detected in the full-energy peak." );
      
      if( distance >= 0.0 )
      {
        cell = table->elementAt( table->rowCount(), 0 );
        txt = new WText( "Solid angle fraction", cell );
        val = SpecUtils::printCompact( geom_eff, 5 );
        cell = table->elementAt( table->rowCount() - 1, 1 );
        txt = new WText( val, cell );
        addTooltipToRow( "The fraction of the solid angle, the detector face takes up, at the specified distance." );
      }//if( distance >= 0.0 )
    }//if( drf )
    
    if( shield_transmission != 1.0 )
    {
      const double shield_trans = gammas_per_bq_into_4pi / branch_ratio / input.spectrum->live_time();
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Shielding transmission", cell );
      val = SpecUtils::printCompact( shield_transmission, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The fraction of gammas, at this energy, that will make it through the shielding without interacting." );
    }//if( branch_ratio != counts_per_bq_into_4pi )
    
    if( air_atten )
    {
      const double air_trans = gammas_per_bq_into_4pi_with_air / gammas_per_bq_into_4pi;
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Air transmission", cell );
      val = SpecUtils::printCompact( air_trans, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The fraction of gammas, at this energy, that will make it through the air (assuming sea level) without interacting." );
    }//if( air_atten )
    
    if( branch_ratio > 0.0 )
    {
      cell = table->elementAt( table->rowCount(), 0 );
      txt = new WText( "Nuclide branching ratio", cell );
      val = SpecUtils::printCompact( branch_ratio, 5 );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      txt = new WText( val, cell );
      addTooltipToRow( "The number of gamma rays emitted at this energy, from the radioactive"
                      " source before any shielding, but accounting for nuclide age,"
                      " per decay of the parent nuclide." );
    }//if( branch_ratio > 0.0 )
  }catch( std::exception &e )
  {
    cerr << "Error computing Currie limit information: " << e.what() << endl;
    WText *msg = new WText( "Error computing Currie limit information", dialog->contents() );
    msg->addStyleClass( "content" );
    msg->setInline( false );
  }//try / catch
  
  return dialog;
}//SimpleDialog *createCurrieRoiMoreInfoWindow()


Wt::Json::Object DetectionLimitTool::generateChartJson( const DetectionLimitCalc::DeconActivityOrDistanceLimitResult &result, const bool is_dist_limit )
{
  
  const bool useCurie = use_curie_units();
  const vector<pair<double,double>> &chi2s = result.chi2s;
  
  const double avrg_quantity = 0.5*(chi2s.front().first + chi2s.back().first);
  const pair<string,double> &units = is_dist_limit
  ? PhysicalUnits::bestLengthUnitHtml( avrg_quantity )
  : PhysicalUnits::bestActivityUnitHtml( avrg_quantity, useCurie );
  
  double minchi2 = std::numeric_limits<double>::infinity();
  for( size_t i = 0; i < chi2s.size(); ++i )
    minchi2 = std::min( minchi2, chi2s[i].second );
  
  double maxchi2 = -std::numeric_limits<double>::infinity();
  
  // The first/last Chi2 may be a huge (like when distance approaches 0), so we'll protect
  //  against this.
  // TODO: use `result.foundUpperCl` and `result.foundLowerCl` to see if we need this, and also, perhaps maybe we should be smarter in picking which points to compute the Chi2 for
  const double maxYRange = 15; //arbitrary
  
  
  Wt::Json::Object json;
  Wt::Json::Array chi2_xy;
  for( size_t i = 0; i < chi2s.size(); ++i )
  {
    if( chi2s[i].second > (minchi2 + maxYRange) )
      continue;
    
    maxchi2 = std::max( maxchi2, chi2s[i].second );
    
    Wt::Json::Object datapoint;
    datapoint["x"] = (chi2s[i].first / units.second);
    datapoint["y"] = chi2s[i].second;
    chi2_xy.push_back( std::move(datapoint) );
  }
  
  if( result.foundLowerCl && (result.lowerLimit > 0.0) )
  {
    //display lower limit to user...
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
  InterSpec *interspec = InterSpec::instance();
  std::shared_ptr<const ColorTheme> theme = interspec ? interspec->getColorTheme() : nullptr;
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
  
  return json;
};//generateChartJson lambda



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
                                            double new_roi_px,
                                            double original_lower_energy,
                                            const std::string &spectrum_type,
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
    log_developer_error( __func__, msg );
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
      age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
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
  
  const double confLevel = currentConfidenceLevel();
  
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
    input.shield_transmission = line.shield_transmission;
    input.counts_per_bq_into_4pi__ = line.gammas_into_4pi*spec->live_time();
    input.counts_per_bq_into_4pi_with_air = line.gammas_4pi_after_air_attenuation*spec->live_time();
    input.distance = distance;
    input.activity = activity;
    input.roi_start = energy - 1.25*fwhm; // recommended by ISO 11929:2010, could instead use 1.19
    input.roi_end = energy + 1.25*fwhm;
    input.num_side_channels = 4;
    input.confidence_level = confLevel;
    //input.trans_through_air = line.air_transmission;
    //input.trans_through_shielding = line.shield_transmission;
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


double DetectionLimitTool::currentConfidenceLevel()
{
  const auto cl = ConfidenceLevel( m_confidenceLevel->currentIndex() );
  switch( cl )
  {
    case NinetyFivePercent: return 0.95;
    case NinetyNinePercent: return 0.99;
    case OneSigma:          return 0.682689492137086;
    case TwoSigma:          return 0.954499736103642;
    case ThreeSigma:        return 0.997300203936740;
    case FourSigma:         return 0.999936657516334;
    case FiveSigma:         return 0.999999426696856;
    case NumConfidenceLevel: break;
  }//switch( cl )
  
  assert( 0 );
  return 0.95;
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


std::string DetectionLimitTool::encodeStateToUrl() const
{
  cerr << "`DetectionLimitTool::encodeStateToUrl()`: Not Implemented - returning empty string!!!" << endl;
  return "";
}


void DetectionLimitTool::handleAppUrl( std::string query_str )
{
  throw runtime_error( "DetectionLimitTool::handleAppUrl(...) not implemented!" );
}


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
            const double counts_4pi = input.counts_per_bq_into_4pi__;
            
            const double activity = other_quantity;
            
            const double nominal_dist = sqrt(counts_4pi * det_eff_1m * activity * def_dist * def_dist / rw->simple_excess_counts());
            
            min_search_quantity = std::min( min_search_quantity, nominal_dist );
            max_search_quantity = std::max( max_search_quantity, nominal_dist );
          }else
          {
            const double det_eff = fixed_geom ? input.drf->intrinsicEfficiency(input.energy)
            : input.drf->efficiency(input.energy, input.distance);
            const double counts_4pi = (air_atten ? input.counts_per_bq_into_4pi_with_air
                                       : input.counts_per_bq_into_4pi__);
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
    
    const double mid_search_quantity = 0.5*(min_search_quantity + max_search_quantity);
    const double base_act = is_dist_limit ?  other_quantity : mid_search_quantity;
    const double base_dist = is_dist_limit ? mid_search_quantity : other_quantity;
    vector<PeakDef> dummy_peaks{};
    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> base_input
                    = getComputeForActivityInput( base_act, base_dist, dummy_peaks );
    
    
    const float wantedCl = currentConfidenceLevel();
    DetectionLimitCalc::DeconActivityOrDistanceLimitResult result
                  = DetectionLimitCalc::get_activity_or_distance_limits( wantedCl, base_input,
                                                                  is_dist_limit, min_search_quantity,
                                                                    max_search_quantity, useCurie );
    
    // TODO: This `m_upperLimit` text could be moved to where the warnings is, ot the title area where the "Gamma lines to use" is
    //       This would give us more room, and maybe be more obvious to the user
    m_upperLimit->setText( WString::fromUTF8(result.limitText) );
    
    if( is_dist_limit )
      m_displayDistance->setText( WString::fromUTF8(result.quantityLimitStr) );
    else
      m_displayActivity->setText( WString::fromUTF8(result.quantityLimitStr) );
    
    m_bestChi2Act->setText( WString::fromUTF8(result.bestCh2Text) );
    
    m_results->show();
    
    const string jsgraph = m_chi2Chart->jsRef() + ".chart";
    const Wt::Json::Object chartJson = generateChartJson( result, is_dist_limit );
    const string datajson = Wt::Json::serialize(chartJson);
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
  
  /*
   vector<DetectionLimitTool::MdaPeakRowInput> input_rois;
   
   for( auto w : m_peaks->children() )
   {
     auto rw = dynamic_cast<MdaPeakRow *>( w );
     assert( rw );
     if( rw && rw->input().use_for_likelihood )
       input_rois.push_back( rw->input() );
   }//for( auto w : m_peaks->children() )
   */
  
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
    peak_info.counts_per_bq_into_4pi = row_input.counts_per_bq_into_4pi__;
    
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
  ref_input.m_showEscapes = false;
    
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
