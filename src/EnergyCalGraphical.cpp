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
#include <ctime>
#include <string>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WGroupBox>
#include <Wt/WCheckBox>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WPushButton>

#include "InterSpec/SpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/EnergyCalTool.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/EnergyCalGraphical.h"
#include "InterSpec/NativeFloatSpinBox.h"


using namespace std;
using namespace Wt;


EnergyCalGraphicalConfirm::EnergyCalGraphicalConfirm( double lowe, double highe,
                            EnergyCalTool *cal, time_t lastCal,
                            EnergyCalGraphicalConfirm::RecalTypes lastType,
                            float lastEnergy )
: AuxWindow( "Confirm Recalibration",
             (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
               | AuxWindowProperties::SetCloseable | AuxWindowProperties::DisableCollapse) ),
  m_calibrator( cal ), m_typeButtons( NULL), m_foregroundOnly( NULL ),
  m_startE( NULL ), m_finalE( NULL ),
  m_preserveLastCal( NULL ),
  m_lastType( lastType ),
  m_lastEnergy( lastEnergy )
{
  setResizable( false );
  rejectWhenEscapePressed();
  centerWindow();
  show();

  auto viewer = InterSpec::instance();
  assert( viewer );
  
  auto disp_meas = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const SpecUtils::EnergyCalibration> energycal = disp_meas
                                                             ? disp_meas->energy_calibration()
                                                             : nullptr;
  
  if( !disp_meas || !energycal || !energycal->valid()
      || (energycal->type() == SpecUtils::EnergyCalType::LowerChannelEdge) )
  {
    throw runtime_error( "You must have a dispayed spectrum to do a recalibration" );
    return;
  }//if( !spectrum->m_interspec->displayedHistogram(SpectrumType::Foreground) )
  
  WTable *table = new WTable( contents() );
  //    table->setHeaderCount( 1, Wt::Vertical );
  new WLabel( "Original Energy", table->elementAt(0,0) );
  m_startE = new NativeFloatSpinBox( table->elementAt(0,1) );
  new WLabel( "Modified Energy", table->elementAt(1,0) );
  m_finalE = new NativeFloatSpinBox( table->elementAt(1,1) );
  m_startE->setMaximum( 20000 );
  m_finalE->setMaximum( 20000 );
  setEnergies( lowe, highe );
  
  
  m_typeButtons = new WButtonGroup( contents() );
  WGroupBox *buttonBox = new WGroupBox( "Parameter to adjust", contents() );
  
  for( RecalTypes t = RecalTypes(0); t < NumRecalTypes; t = RecalTypes(t+1) )
  {
    const char *txt = "";
    switch( t )
    {
      case kOffset:       txt = "Offset";             break;
      case kLinear:       txt = "Linear";             break;
//      case kQuadratic:    txt = "Quadratic";          break;
      case kDeviation:    txt = "Add Deviation Pair"; break;
      case NumRecalTypes: txt = "";                   break;
    }//switch( t )
    
    WRadioButton *button = new WRadioButton( txt, buttonBox );
    button->setInline( false );
    m_typeButtons->addButton( button, t );
  }//for( loop over recal types )
  m_typeButtons->setSelectedButtonIndex( kLinear );
 
  time_t currenttime;
  time( &currenttime );
  const double secondsLastCal = difftime( currenttime, lastCal );
  const float energyrange = energycal->upper_energy() - energycal->lower_energy();
  const float prevdiff = std::max( fabs(m_lastEnergy - lowe), fabs(m_lastEnergy - highe) );
  
  if( secondsLastCal < 120.0
     && (m_lastType==kOffset || m_lastType==kLinear)
     && (m_lastEnergy > 0.0)
     && (prevdiff > 0.1*energyrange) ) //make sure the two points arent right next to eachother; we could probably increase this to like 0.2 without issue
  {
    char msg[128];
    snprintf(msg, sizeof(msg), "Preserve %.1f keV Cal.", m_lastEnergy );
    m_preserveLastCal = new WCheckBox( msg, contents() );
    m_preserveLastCal->setInline( false );
    m_preserveLastCal->setStyleClass( "PreserveLastCalCb" );
    
    const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", viewer );

    HelpSystem::attachToolTipOn( m_preserveLastCal,"This is only possible if a offset or"
                                " linear term adjustment was previously"
                                " made within the last 2 minutes.", showToolTips );
    m_preserveLastCal->setChecked();
    buttonBox->disable();
    m_preserveLastCal->checked().connect( buttonBox, &WWidget::disable );
    m_preserveLastCal->unChecked().connect( buttonBox, &WWidget::enable );
  }//if( preserve last cal possibly )
  
  const string applyToTxt = cal->applyToSummaryTxt();
  if( !applyToTxt.empty() )
  {
    WText *t = new WText( "Changes will be applied to<br />" + applyToTxt, XHTMLText, contents() );
    t->setInline( false );
    t->setAttributeValue( "style", "color: #737373; width: auto; text-align: center;" );
  }
  
  
  /*
  auto secondaryH  = viewer->displayedHistogram(SpecUtils::SpectrumType::SecondForeground);
   auto backgroundH = viewer->displayedHistogram(SpecUtils::SpectrumType::Background);
   
   if( secondaryH || backgroundH )
   {
     m_foregroundOnly = new WCheckBox( "Foreground Only", contents() );
     m_foregroundOnly->setInline( false );
     m_foregroundOnly->setAttributeValue( "style",
                 "margin-left:1.25em;margin-top:0.25em;margin-bottom=0.25em;" );
   }//if( secondaryH || backgroundH )
   
  
  const auto meas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const vector<string> detectors = viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  
  bool displayAll = true;
  if( meas )
  {
    for( const string &name : meas->detector_names() )
    {
      const auto pos = std::find( begin(detectors), end(detectors), name );
      displayAll = displayAll && (pos != end(detectors));
    }
  }//if( meas )
  
  if( !displayAll )
  {
    WText *t = new WText( "(Calibration only applied to<br />displayed detectors)",
                         XHTMLUnsafeText, contents() );
    t->setInline( false );
    t->setAttributeValue( "style", "color: #737373; width: auto; text-align: center;" );
  }//if( !displayAll )
*/
   
   
  AuxWindow::addHelpInFooter( footer(), "graphical-recal-dialog" );
  
  
  WPushButton *button = new WPushButton( "Cancel", footer() );
  button->clicked().connect( this, &AuxWindow::hide );
  
  button = new WPushButton( "Accept", footer()  );
  button->setIcon( "InterSpec_resources/images/accept.png" );
  button->clicked().connect( this, &EnergyCalGraphicalConfirm::apply );
  
  //Could place dialog near where the mouse is
  //  stringstream moveJs;
  //  moveJs << "var el=" << this->jsRef() << ";el.style.top='" << y << "px';"
  //         << "el.style.left='" << x << "px';";
  //  doJavaScript( moveJs.str();
}//EnergyCalGraphicalConfirm


EnergyCalGraphicalConfirm::~EnergyCalGraphicalConfirm()
{
}


void EnergyCalGraphicalConfirm::apply()
{
  const float startE = static_cast<float>(m_startE->value());
  const float finalE = static_cast<float>(m_finalE->value());
  const RecalTypes type = RecalTypes( m_typeButtons->selectedButtonIndex() );
  
  if( startE == finalE )
  {
    finished().emit(WDialog::Accepted);
    return;
  }//if( startE == finalE )
  
  
  //Need to implement recalibrating logic, and also logic to check how long
  //  ago user last calibrated this spectra...
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  shared_ptr<SpecMeas> foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const SpecUtils::Measurement> displ_foreground
                                  = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  if( !foreground || !displ_foreground )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  shared_ptr<const SpecUtils::EnergyCalibration> energycal = displ_foreground->energy_calibration();
  if( !energycal || !energycal->valid()
      || (energycal->type() == SpecUtils::EnergyCalType::LowerChannelEdge) )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  const float shift = finalE - startE;
  
  pair<float,float> added_dev_pair = {0.0f, 0.0f};
  const vector<pair<float,float>> orig_dev_pairs = energycal->deviation_pairs();
  const size_t nbin = foreground->num_gamma_channels();
    
  bool isOffsetOnly = false;
  
  //m_preserveLastCal will only be valid if m_lastType is kOffset or kLinear,
  // and m_lastEnergy > 0.0
  const bool preserveLast = (m_preserveLastCal && m_preserveLastCal->isChecked()
                             && (m_lastEnergy!=startE));
  
  vector<float> eqn = energycal->coefficients();
  assert( eqn.size() >= 2 );
  
  if( preserveLast )
  {
    const double lowchannel = energycal->channel_for_energy( m_lastEnergy );
    const double upchannel = energycal->channel_for_energy( startE );
    
    if( energycal->type() == SpecUtils::EnergyCalType::FullRangeFraction )
    {
      //From a quick empiracal test (so by no means definitive), the below
      //  behaves well for the case deviation pairs are defined.
      
      const double x1 = upchannel / nbin;
      const double x2 = lowchannel / nbin;
      const double E1 = eqn[0] + eqn[1]*x1;
      const double E2 = eqn[0] + eqn[1]*x2;
      eqn[1] = static_cast<float>( (E1 - E2 + shift) / (x1-x2) );
      eqn[0] = static_cast<float>( E2 - x2*eqn[1] );
    }else
    {
      assert( energycal->type() == SpecUtils::EnergyCalType::Polynomial
              || energycal->type() == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial );
      
      const float E1 = eqn[0] + eqn[1]*upchannel;
      const float E2 = eqn[0] + eqn[1]*lowchannel;
      eqn[1] = (E1 - E2 + shift) / (upchannel - lowchannel);
      eqn[0] = E2 - lowchannel*eqn[1];
    }
  }else
  {
    switch( type )
    {
      case kOffset:
        eqn[0] += shift;
        isOffsetOnly = true;
      break;
      
      case kLinear:
      {
        //From a quick empiracal test (so by no means definitive), the below
        //  behaves well for the case deviation pairs are defined.
        const double channel_num = energycal->channel_for_energy( startE );
        
        if( energycal->type() == SpecUtils::EnergyCalType::FullRangeFraction )
        {
          eqn[1] = 0.0;
          double otherterms = SpecUtils::fullrangefraction_energy( channel_num, eqn, nbin, {} );
          otherterms -= SpecUtils::correction_due_to_dev_pairs( finalE, orig_dev_pairs );
          
          eqn[1] = (finalE - otherterms) * nbin / channel_num ;
        }else
        {
          assert( energycal->type() == SpecUtils::EnergyCalType::Polynomial
                  || energycal->type() == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial );
          
          eqn[1] = 0.0;
          double otherterms = SpecUtils::polynomial_energy( channel_num, eqn, {} );
          otherterms -= SpecUtils::correction_due_to_dev_pairs( finalE, orig_dev_pairs );
          
          eqn[1] = (finalE - otherterms) / channel_num;
        }
        break;
      }//case kLinear:
            
//      case kQuadratic:
//      break;
      
      case kDeviation:
      {
        if( orig_dev_pairs.empty() )
        {
          added_dev_pair = { startE, shift };
        }else
        {
          //We are going to to find the bin number startE cooresponds to, then
          //  find the energy this would have cooresponded to if there were no
          //  deviation pairs, then we will add a deviation
          const double bin = energycal->channel_for_energy( startE );
          
          double originalE = 0.0;
          if( energycal->type() == SpecUtils::EnergyCalType::FullRangeFraction )
          {
            originalE = SpecUtils::fullrangefraction_energy( bin, energycal->coefficients(), nbin, {} );
          }else
          {
            assert( energycal->type() == SpecUtils::EnergyCalType::Polynomial
                    || energycal->type() == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial );
            originalE = SpecUtils::polynomial_energy( bin, energycal->coefficients(), {} );
          }
          
          const float offset = finalE - originalE;
          added_dev_pair = { static_cast<float>(finalE), static_cast<float>(offset) };
        }//if( orig_dev_pairs.empty() )
      
        break;
      }//case kDeviation:
      
      case NumRecalTypes:
      break;
    }//switch( type )
  }//if( preserveLast ) / else
  
  auto newcal = make_shared<SpecUtils::EnergyCalibration>();
  
  try
  {
    switch( energycal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        newcal->set_polynomial( nbin, eqn, orig_dev_pairs );
        break;
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
        newcal->set_full_range_fraction( nbin, eqn, orig_dev_pairs );
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        assert( 0 );
        break;
    }//switch( energycal->type() )
  }catch( std::exception & )
  {
    string msg;
    msg = "The coeffiecients {";
    for( size_t i = 0; i < eqn.size(); ++i )
      msg += (i ? ", " : "") + to_string( eqn[i] );
    msg += "} ";
    if( orig_dev_pairs.size() )
    {
      msg += " and deviation pairs {";
      for( size_t i = 0; i < orig_dev_pairs.size(); ++i )
        msg += (i ? ", [" : "[") + to_string(orig_dev_pairs[i].first)
              + "," + to_string(orig_dev_pairs[i].second) + "]";
      msg += "} ";
    }
    msg += "are invalid because they will cause higher numbered bins to have lower energy values."
           " Calibration not applied.";
    
    viewer->logMessage( WString::fromUTF8(msg), "", 3 );
    finished().emit(WDialog::Accepted);
    return;
  }//try / catch
  
  assert( newcal && newcal->valid() );
  
  try
  {
    
    if( preserveLast || (type == kOffset) || (type == kLinear) )
    {
      m_calibrator->applyCalChange( energycal, newcal, isOffsetOnly );
    }else if( type == kDeviation )
    {
      m_calibrator->addDeviationPair( added_dev_pair );
    }else
    {
      assert( 0 );
    }
   
    m_calibrator->setWasGraphicalRecal( type, finalE );
  }catch( std::exception &e )
  {
    viewer->logMessage( e.what(), "", 2 );
  }//try / catch
  
  finished().emit(WDialog::Accepted);
}//void apply()


void EnergyCalGraphicalConfirm::setEnergies( double xstart, double xfinish )
{
  m_startE->setValue( floor(100.0*xstart+0.5)/100.0 );
  m_finalE->setValue( floor(100.0*xfinish+0.5)/100.0 );
}//void setEnergies(...)

