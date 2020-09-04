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
#include <vector>

#include <Wt/WText>
#include <Wt/WGroupBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WButtonGroup>
#include <Wt/WApplication>
#include <Wt/WRadioButton>
#include <Wt/WContainerWidget>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/EnergyCalTool.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/EnergyCalMultiFile.h"
#include "InterSpec/EnergyCalAddActions.h"

// Forward declarations
class ConvertToPolyTool;

using namespace std;
using namespace Wt;


class ConvertToPolyTool : public WContainerWidget
{
  EnergyCalTool *m_cal;
  WButtonGroup *m_group;
  WText *m_fitCoefTxt;
  
  WPushButton *m_cancel;
  WPushButton *m_accept;
  
public:
  ConvertToPolyTool( EnergyCalTool *cal, AuxWindow *parent );
  ~ConvertToPolyTool();
  
  void fromLowerChannelEnergiesOptionChanged();
  void handleFinish( Wt::WDialog::DialogCode result );
};//class ConvertToPolyTool






EnergyCalAddActionsWindow::EnergyCalAddActionsWindow( const MoreActionsIndex actionType,
                             const std::vector<MeasToApplyCoefChangeTo> &measToChange,
                             EnergyCalTool *calibrator )
  : AuxWindow( "&nbsp", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                      | AuxWindowProperties::DisableCollapse) ),
    m_actionType( actionType ),
    m_measToChange( make_shared<vector<MeasToApplyCoefChangeTo>>(measToChange)),
    m_calibrator( calibrator )
{
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalAddActions.css");
  
  switch( m_actionType )
  {
    case MoreActionsIndex::Linearize:
      break;
      
    case MoreActionsIndex::Truncate:
      break;
      
    case MoreActionsIndex::CombineChannels:
      break;
      
    case MoreActionsIndex::ConvertToFrf:
      break;
      
    case MoreActionsIndex::ConvertToPoly:
      AuxWindow::setWindowTitle( "Convert to Polynomial Calibration" );
      new ConvertToPolyTool( m_calibrator, this );
      break;
      
    case MoreActionsIndex::MultipleFilesCal:
      AuxWindow::setWindowTitle( "Multi-File Calibration" );
      new EnergyCalMultiFile( m_calibrator, this );
      break;
      
    case MoreActionsIndex::NumMoreActionsIndex:
      break;
  }//switch( m_actionType )
  
  rejectWhenEscapePressed();
  AuxWindow::show();
  setWidth( 400.0 );
  AuxWindow::resizeToFitOnScreen();
  AuxWindow::centerWindow();
}//EnergyCalAddActionsWindow constructor
  
EnergyCalAddActionsWindow::~EnergyCalAddActionsWindow()
{
}
 





ConvertToPolyTool::ConvertToPolyTool( EnergyCalTool *cal, AuxWindow *parent )
  : WContainerWidget(),
    m_cal( cal ),
    m_group( nullptr ),
    m_fitCoefTxt( nullptr ),
    m_cancel( nullptr ),
    m_accept( nullptr )
{
  addStyleClass( "ConvertToPolyTool" );
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  if( parent )
    parent->stretcher()->addWidget( this, 0, 0  );
  
  set<size_t> coeforder, nchannels;
  set<SpecUtils::EnergyCalType> caltypes;
  vector<MeasToApplyCoefChangeTo> applicables = cal->measurementsToApplyCoeffChangeTo();
  for( const auto &toapplyto : applicables )
  {
    for( const auto &detname : toapplyto.detectors )
    {
      for( const int sample : toapplyto.sample_numbers )
      {
        auto m = toapplyto.meas->measurement( sample, detname );
        if( !m || !m->energy_calibration() || !m->energy_calibration()->valid()
           || (m->num_gamma_channels() < 5) )
          continue;
        
        nchannels.insert( m->num_gamma_channels() );
        caltypes.insert( m->energy_calibration()->type() );
        coeforder.insert( m->calibration_coeffs().size() );
      }//for( loop over samples )
    }//for( loop over detectors )
  }//for( const auto &m : applicables )
  
  if( (caltypes.size() != 1)
     || ((*begin(caltypes)) == SpecUtils::EnergyCalType::Polynomial)
     || ((*begin(caltypes)) == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial) )
  {
    if( parent )
    {
      auto close = parent->addCloseButtonToFooter();
      close->clicked().connect( boost::bind( &ConvertToPolyTool::handleFinish, this, WDialog::Rejected ) );
      parent->finished().connect( this, &ConvertToPolyTool::handleFinish );
    }
    
    const char *msgtxt = "";
    if( caltypes.empty() )
      msgtxt = "Conversion would not effect any spectra; try selecting more &quot;Apply Changes To&quot; criteria.";
    else if( caltypes.size() == 1 )
      msgtxt = "All spectra this would be applied to is already polynomial.";
    else
      msgtxt = "There is more than one calibration type of selected source spectra; try restricting &quot;Apply Changes To&quot; criteria.";
    
    WText *msg = new WText( msgtxt, this );
    msg->addStyleClass( "ConvertToNA" );
    msg->setInline( false );
    
    return;
  }//if( caltypes.size() != 1 || existing cal if polynomial )
  
  
  string applyToTxt = cal->applyToSummaryTxt();
  if( !applyToTxt.empty() )
  {
    applyToTxt = "Changes will be applied to: " + applyToTxt;
    WText *msg = new WText( applyToTxt, this );
    msg->addStyleClass( "ConvertToApplieTo" );
    msg->setInline( false );
  }//if( !applyToTxt.empty() )
  
  
  const SpecUtils::EnergyCalType origtype = *begin(caltypes);
  switch( origtype )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      break;
      
    case SpecUtils::EnergyCalType::FullRangeFraction:
    {
      assert( coeforder.size() );
      const size_t maxorder = *std::rbegin(coeforder);
      const char *msgtxt = "";
      if( maxorder > 4 )
        msgtxt = "Terms past the fourth term will be lost since there is no analougous terms.";
      else
        msgtxt = "No information will be lost if you continue.";
      WText *msg = new WText( msgtxt, this );
      msg->setInline( false );
      msg->addStyleClass( "ConvertToMsg" );
      break;
    }
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    {
      WText *msg = new WText( "How would you like the conversion to be performed?", this );
      msg->addStyleClass( "ConvertToMsg" );
      msg->setInline( false );
      
      m_group = new WButtonGroup( this );
      WGroupBox *box = new WGroupBox( this );
      
      WRadioButton *btn = new WRadioButton( "Linearize, rebinning counts to match new widths", box );
      btn->setToolTip( "The channels in the resulting spectra will all have the same width,"
                      " but counts in each channel will not be integer" );
      btn->setInline( false );
      m_group->addButton( btn, 0 );
      
      btn = new WRadioButton( "Linearize, leaving counts in channels the same", box );
      btn->setToolTip( "The channel energies will be reasigned so that each channel will have"
                      " same width, without taking into account the channels original width;"
                      " i.e., the counts in each channel will remain the same." );
      btn->setInline( false );
      m_group->addButton( btn, 1 );
      
      btn = new WRadioButton( "Fit channel energies to quadratic equation.", box );
      btn->setToolTip( "" );
      btn->setInline( false );
      m_group->addButton( btn, 2 );
      
      btn = new WRadioButton( "Fit channel energies to quartic equation.", box );
      btn->setToolTip( "" );
      btn->setInline( false );
      m_group->addButton( btn, 3 );
      
      m_group->setSelectedButtonIndex( 0 );
      m_group->checkedChanged().connect( this, &ConvertToPolyTool::fromLowerChannelEnergiesOptionChanged );
      
      break;
    }
  }//switch( origtype )
  
  
  WContainerWidget *buttonDiv = nullptr;
  if( parent )
  {
    buttonDiv = parent->footer();
  }else
  {
    buttonDiv = new WContainerWidget( this );
  }
  
  //AuxWindow::addHelpInFooter( buttonDiv, "convert-to-polynomial-dialog" );
  
  m_cancel = new WPushButton( "Cancel", buttonDiv );
  m_accept = new WPushButton( "Accept", buttonDiv );
  
  m_cancel->clicked().connect( boost::bind( &ConvertToPolyTool::handleFinish, this, WDialog::Rejected ) );
  m_accept->clicked().connect( boost::bind( &ConvertToPolyTool::handleFinish, this, WDialog::Accepted ) );
  
  if( parent )
  {
    parent->finished().connect( this, &ConvertToPolyTool::handleFinish );
    
    const int w = 600 < viewer->renderedWidth() ? 600 : viewer->renderedWidth();
    //const int h = static_cast<int>(0.8*viewer->renderedHeight());
    //parent->resizeWindow( w, h );
    parent->setWidth( w );
    
    parent->rejectWhenEscapePressed();
    
    parent->centerWindow();
  }//if( parent )
}//ConvertToPolyTool(...)
  

ConvertToPolyTool::~ConvertToPolyTool()
{
}
  
  
void ConvertToPolyTool::fromLowerChannelEnergiesOptionChanged()
{
  assert( m_group );
  
  if( !m_fitCoefTxt )
  {
    m_fitCoefTxt = new WText( "&nbsp;", this );
    m_fitCoefTxt->setInline( false );
    m_fitCoefTxt->addStyleClass( "ConvertFromLCERes" );
  }//if( !m_fitCoefTxt )
  
  const int index = m_group->selectedButtonIndex();
  
  try
  {
    m_accept->enable();
    
    switch( index )
    {
      case 0: //Linearize, rebinning counts to match
      case 1: //Linearize, leaving counts in channels the
        m_fitCoefTxt->setText( "" );
        return;
        
      case 2: //Fit channel energies to quadratic equation
      case 3: //Fit channel energies to quartic equation
        break;
        
      default:
        throw runtime_error( "Unexpected option selected" );
    }//switch( index )
    
    string msgsumm;
    const size_t ncoeffs = static_cast<size_t>(index+1);
    vector<MeasToApplyCoefChangeTo> applicables = m_cal->measurementsToApplyCoeffChangeTo();
    for( const auto &toapplyto : applicables )
    {
      if( !msgsumm.empty() )
        break;
      
      for( const auto &detname : toapplyto.detectors )
      {
        if( !msgsumm.empty() )
          break;
        
        for( const int sample : toapplyto.sample_numbers )
        {
          if( !msgsumm.empty() )
            break;
          
          auto m = toapplyto.meas->measurement( sample, detname );
          if( !m
             || !m->energy_calibration()
             || !m->energy_calibration()->valid()
             || (m->num_gamma_channels() < 5)
             || (m->energy_calibration()->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
            continue;
          
          assert( m->energy_calibration()->channel_energies() );
          
          const auto &lower_energies = *m->energy_calibration()->channel_energies();
          
          vector<float> coefs;
          const double avrgdiff = EnergyCal::fit_poly_from_lower_channel_energies( ncoeffs,
                                                                            lower_energies, coefs );
          
          char buffer[32];
          msgsumm = "Fit coefficients {";
          for( size_t i = 0; i < coefs.size(); ++i )
          {
            snprintf( buffer, sizeof(buffer), "%s%1.57g", (i ? ", " : ""), coefs[i] );
            msgsumm += buffer;
          }
          snprintf( buffer, sizeof(buffer), "%1.57g", avrgdiff );
          msgsumm += "} with an average error of " + string(buffer) + " keV in each channel";
        }//for( const int sample : toapplyto.sample_numbers )
      }//for( const auto &detname : toapplyto.detectors )
    }//for( const auto &toapplyto : applicables )
    
    if( msgsumm.empty() )
    {
      m_fitCoefTxt->setText( "No measurements would be effected?  Something may be a bit whack." );
    }else
    {
      m_fitCoefTxt->setText( WString::fromUTF8(msgsumm) );
    }// if( msgsumm.empty() ) / else
    
  }catch( std::exception &e )
  {
    m_fitCoefTxt->setText( "Error converting coefficients: " + string( e.what() ) );
    m_accept->disable();
  }//try / catch
}//void fromLowerChannelEnergiesOptionChanged()


void ConvertToPolyTool::handleFinish( Wt::WDialog::DialogCode result )
{
  
}//void handleFinish(...)

