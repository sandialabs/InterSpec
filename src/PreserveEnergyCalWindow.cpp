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

#include <vector>

#include <Wt/WText>
#include <Wt/WPushButton>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/PreserveEnergyCalWindow.h"

using namespace Wt;
using namespace std;

using SpecUtils::Measurement;
using SpecUtils::SpectrumType;


PreserveEnergyCalWindow::PreserveEnergyCalWindow(
                             std::shared_ptr<SpecMeas> newmeas,
                             const SpecUtils::SpectrumType newtype,
                             std::shared_ptr<SpecMeas> oldmeas,
                             const SpecUtils::SpectrumType oldtype,
                             EnergyCalTool *calibrator )
: AuxWindow( "Keep Calibration?",
             Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneModal) ),
  m_calibrator( calibrator ),
  m_newmeas( newmeas ),
  m_newtype( newtype ),
  m_oldtype( oldtype )
{
  if( !newmeas || !oldmeas || !calibrator )
    throw runtime_error( "PreserveEnergyCalWindow: invalid input" );
  
  std::shared_ptr<const Measurement> eqnoldmeas, eqnnewmeas;
  const vector< std::shared_ptr<const Measurement> > oldmeass = oldmeas->measurements();
  const vector< std::shared_ptr<const Measurement> > newmeass = newmeas->measurements();
  
  for( size_t i = 0; !eqnoldmeas && i < oldmeass.size(); ++i )
    if( oldmeass[i]->num_gamma_channels() )
      eqnoldmeas = oldmeass[i];
  for( size_t i = 0; !eqnnewmeas && i < newmeass.size(); ++i )
    if( newmeass[i]->num_gamma_channels() )
      eqnnewmeas = newmeass[i];
  
  if( !eqnoldmeas || !eqnnewmeas )
    throw runtime_error( "PreserveEnergyCalWindow: invalid input 2" );

  
  m_type = eqnoldmeas->energy_calibration_model();
  m_coeffs = eqnoldmeas->calibration_coeffs();
  m_devPairs = eqnoldmeas->deviation_pairs();
  
  switch( m_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      break;
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      throw runtime_error( "PreserveEnergyCalWindow: invalid type" );
    break;
  }//switch( m_devPairs )
  
  string msg = "<div>It looks like you have loaded a ";
  
  msg += descriptionText(newtype);
  msg += " spectrum";
  
  if( m_newtype == m_oldtype )
  {
    msg += " that is from the same detector as the previous spectrum";
  }else
  {
    msg += " that is from the same detector as the ";
    msg += descriptionText(m_oldtype);
    msg += " spectrum";
  }

  msg += ", but with a different calibration.</div>"
  "<table class=\"RecalibCoefTable\"><tr><th style=\"padding-right: 5px;\">Order</th><th style=\"padding-right: 5px;\">Previous (";
  
  switch( eqnnewmeas->energy_calibration_model() )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      msg += "Polynomial";
      break;
    case SpecUtils::EnergyCalType::FullRangeFraction: msg += "FRF"; break;
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      break;
  }//switch( newmeas->energy_calibration_model() )
  
  msg += ")</th><th>New (";
  switch( eqnoldmeas->energy_calibration_model() )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      msg += "Polynomial";
      break;
    case SpecUtils::EnergyCalType::FullRangeFraction: msg += "FRF"; break;
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      break;
  }//switch( oldmeas->energy_calibration_model() )
  msg += ")</th></tr>";
  

  const vector<float> &newcoefs = eqnnewmeas->calibration_coeffs();
  const size_t ncoefs = std::max( newcoefs.size(), m_coeffs.size() );
  
  for( size_t i = 0; i < ncoefs; ++i )
  {
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "<tr><td>%i</td>", static_cast<int>(i) );
    msg += buffer;
    
    if( i < m_coeffs.size() )
    {
      snprintf( buffer, sizeof(buffer), "<td>%.4g</td>", m_coeffs[i] );
      msg += buffer;
    }else
      msg += "<td></td>";
    
    if( i < newcoefs.size() )
    {
      snprintf( buffer, sizeof(buffer), "<td>%.4g</td>", newcoefs[i] );
      msg += buffer;
    }else
      msg += "<td></td>";
    msg += "</tr>";
  }//for( size_t i = 0; i < ncoefs; ++i )
  
  msg += "</table>";
  
  msg += "<div style=\"margin-top: 10px;\"><b>"
         "Would you like to use the previous calibration with the new file?"
         "</b></div>";
  
  new WText( msg, XHTMLUnsafeText, contents() );
  
  
  WContainerWidget *bottom = footer();
  //bottom->setStyleClass("modal-footer");
  WPushButton *cancel = new WPushButton( "No", bottom );
  cancel->clicked().connect( this, &PreserveEnergyCalWindow::emitReject );
  cancel->setWidth( WLength(47.5,WLength::Percentage) );
  cancel->setFloatSide(Wt::Left);
  
  WPushButton *use = new WPushButton( "Yes", bottom );
  use->clicked().connect( this, &PreserveEnergyCalWindow::doRecalibration );
  use->setWidth( WLength(47.5,WLength::Percentage) );
  use->setFloatSide( Wt::Right );
  
  rejectWhenEscapePressed();
  AuxWindow::disableCollapse();
  AuxWindow::show();
  AuxWindow::resize( 400.0, WLength::Auto );
  AuxWindow::centerWindow();
}//PreserveEnergyCalWindow


PreserveEnergyCalWindow::~PreserveEnergyCalWindow()
{
}


bool PreserveEnergyCalWindow::candidate( std::shared_ptr<SpecMeas> newmeas,
                                     std::shared_ptr<SpecMeas> oldmeas )
{
  if( !newmeas || !oldmeas || newmeas==oldmeas )
    return false;
  
  
  std::shared_ptr<const Measurement> eqnoldmeas, eqnnewmeas;
  const vector< std::shared_ptr<const Measurement> > oldmeass = oldmeas->measurements();
  const vector< std::shared_ptr<const Measurement> > newmeass = newmeas->measurements();
  
  for( size_t i = 0; !eqnoldmeas && i < oldmeass.size(); ++i )
    if( oldmeass[i]->num_gamma_channels() )
      eqnoldmeas = oldmeass[i];
  for( size_t i = 0; !eqnnewmeas && i < newmeass.size(); ++i )
    if( newmeass[i]->num_gamma_channels() )
      eqnnewmeas = newmeass[i];
  
  if( !eqnoldmeas || !eqnnewmeas )
    return false;
  
  if( newmeas->instrument_id() != oldmeas->instrument_id() )
    return false;
  
  if( newmeas->num_gamma_channels() != oldmeas->num_gamma_channels() )
    return false;
  
  
  
  const vector<float> &oldcoefs = eqnoldmeas->calibration_coeffs();
  const vector< pair<float,float> > &olddevpairs = eqnoldmeas->deviation_pairs();
  const SpecUtils::EnergyCalType oldtype = eqnoldmeas->energy_calibration_model();
  
  const vector<float> &newcoefs = eqnnewmeas->calibration_coeffs();
  const vector< pair<float,float> > &newdevpairs = eqnnewmeas->deviation_pairs();
  const SpecUtils::EnergyCalType newtype = eqnnewmeas->energy_calibration_model();

  if( oldcoefs.size() != newcoefs.size() )
    return true;
  if( olddevpairs.size() != newdevpairs.size() )
    return true;
  if( oldtype != newtype )
    return true;
  
  for( size_t i = 0; i < oldcoefs.size(); ++i )
  {
    const float a = oldcoefs[i];
    const float b = newcoefs[i];
    const float diff = fabs(a - b);
    
    //FLT_EPSILON the minimum positive number such that 1.0 + FLT_EPSILON != 1.0.
    //FLT_EPSILON==1.19209290E-07F
    //FLT_MIN==1.17549435E-38F
    if( (diff > (1.0E-5*std::max(fabs(a),fabs(b))))
        && diff > (std::pow(1.0E-3f,static_cast<float>(i+1))) )
      return true;
  }//for( size_t i = 0; i < oldcoefs.size(); ++i )
  
  for( size_t i = 0; i < olddevpairs.size(); ++i )
  {
    float a = olddevpairs[i].first;
    float b = newdevpairs[i].first;
    float diff = fabs(a - b);
    
    if( (diff > (1.0E-5*std::max(fabs(a),fabs(b))))
       && diff > (std::pow(1.0E-3f,static_cast<float>(i+1))) )
      return true;
    
    a = olddevpairs[i].second;
    b = newdevpairs[i].second;
    diff = fabs(a - b);
    
    if( (diff > (1.0E-5*std::max(fabs(a),fabs(b))))
       && diff > (std::pow(1.0E-3f,static_cast<float>(i+1))) )
      return true;
  }//for( size_t i = 0; i < olddevpairs.size(); ++i )
  
  return false;
}//candidate(...)


void PreserveEnergyCalWindow::doRecalibration()
{
#warning "Need to implement actually doing proper energy cal propogation"
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  std::shared_ptr<SpecMeas> meas = viewer->measurment(m_newtype);
  std::shared_ptr<const Measurement> displ_foreground
                                      = viewer->displayedHistogram(SpectrumType::Foreground);

  
  if( m_newmeas != meas )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  std::shared_ptr<const Measurement> eqnmeas;
  const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
  
  for( size_t i = 0; !eqnmeas && i < meass.size(); ++i )
    if( meass[i]->num_gamma_channels() )
      eqnmeas = meass[i];
  
  if( !eqnmeas )
  {
    finished().emit(WDialog::Accepted);
    return;
  }
  
  
  const vector<float> oldcoefs = eqnmeas->calibration_coeffs();
  const vector< pair<float,float> > olddevpairs = eqnmeas->deviation_pairs();
  const SpecUtils::EnergyCalType oldtype = eqnmeas->energy_calibration_model();
  
  const vector<string> detectors = viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  // Fix this for new calibration! //meas->recalibrate_by_eqn( m_coeffs, m_devPairs, m_type, detectors, false );

  assert( 0 );
  
#warning "Need to implement shifting peaks for calibration"
  //EnergyCalTool::shiftPeaksForEnergyCalibration( viewer->peakModel(),
  //                                              m_coeffs, m_devPairs, m_type,
  //                                              meas, m_newtype, oldcoefs, olddevpairs, oldtype );
  
  viewer->refreshDisplayedCharts();
  m_calibrator->refreshGuiFromFiles();
  
  finished().emit(WDialog::Accepted);
}//void doRecalibration()
