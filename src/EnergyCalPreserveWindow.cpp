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
#include <Wt/WMenu>
#include <Wt/WMenuItem>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/EnergyCal.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/EnergyCalPreserveWindow.h"

using namespace Wt;
using namespace std;

using SpecUtils::Measurement;
using SpecUtils::SpectrumType;
using SpecUtils::EnergyCalType;
using SpecUtils::EnergyCalibration;
using EnergyCal::translatePeaksForCalibrationChange;

namespace
{

// Returns map from detector name, to first energy calibration found for it.
//  Note that this doesnt account for energy calibration possibly being different for a detector for
//  different sample numbers
// TODO: should say currently displayed sample number, wo we can get energy calibration of current sample number
map<string,shared_ptr<const EnergyCalibration>> gamma_names_to_cals( shared_ptr<const SpecMeas> meas )
{
  map<string,shared_ptr<const EnergyCalibration>> gamma_cals;
  
  for( auto &name : meas->gamma_detector_names() )
  {
    shared_ptr<const EnergyCalibration> cal;
    for( auto sample : meas->sample_numbers() )
    {
      auto m = meas->measurement( sample, name );
      if( m && m->energy_calibration()
          && m->energy_calibration()->valid()
          && (m->num_gamma_channels() >= 5) )
      {
        cal = m->energy_calibration();
        break;
      }//if( we found energy calibration for this detector and sample )
    }//for( loop over sample numbers )
    
    if( cal )
      gamma_cals[name] = cal;
  }//for( loop over detector names )
  
  return gamma_cals;
}//gamma_names_to_cals(...)

}//namespace



EnergyCalPreserveWindow::EnergyCalPreserveWindow(
                             std::shared_ptr<SpecMeas> newmeas,
                             const SpecUtils::SpectrumType newtype,
                             std::shared_ptr<SpecMeas> oldmeas,
                             const SpecUtils::SpectrumType oldtype,
                             EnergyCalTool *calibrator )
: AuxWindow( "Keep Previous Calibration?",
             Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen) ),
  m_calibrator( calibrator ),
  m_newmeas( newmeas ),
  m_oldmeas( oldmeas ),
  m_newtype( newtype ),
  m_oldtype( oldtype )
{
  if( !newmeas || !oldmeas || !calibrator )
    throw runtime_error( "EnergyCalPreserveWindow: invalid input" );
  
  if( !candidate(newmeas,oldmeas) )
    throw runtime_error( "EnergyCalPreserveWindow: not candidate to transfer energy calibrations" );
 
  contents()->addStyleClass( "EnergyCalPreserve" );
  
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalPreserveWindow.css" );
  
  const auto new_cals = gamma_names_to_cals( newmeas );
  const auto old_cals = gamma_names_to_cals( oldmeas );
  
  string msg = "<div class=\"EnergyCalPreserveIntroTxt\">It looks like you have loaded a ";
  
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

  msg += ", but with a different calibration.</div>";
  
  auto txt = new WText( msg, XHTMLText, contents() );
  txt->setInline( false );
  
  auto create_table = []( const shared_ptr<const EnergyCalibration> &oldcal,
                         const shared_ptr<const EnergyCalibration> &newcal ) -> std::string {
    if( !oldcal || !newcal )
      return "";
    
    auto cal_type = []( const EnergyCalType type ) -> string {
      switch( type )
      {
        case EnergyCalType::Polynomial:
        case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          return "Polynomial";
          
        case EnergyCalType::FullRangeFraction:
          return "FullRangeFraction";
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
          return "Lower Channel Energies";
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
          break;
      }//switch( newmeas->energy_calibration_model() )
      
      return "";
    }; //cal_type(...) lambda
    
    
    string msg;
    msg += "<table class=\"RecalibCoefTable\">"
             "<tr>"
               "<th>Order</th>"
               "<th>Previous (" + cal_type(oldcal->type()) + ")</th>"
               "<th>New (" + cal_type(newcal->type()) + ")</th>"
             "</tr>";
  
    const auto &oldcoefs = oldcal->coefficients();
    const auto &newcoefs = newcal->coefficients();
    
    const size_t ncoefs = std::max( oldcoefs.size(), newcoefs.size() );
    
    for( size_t i = 0; i < ncoefs; ++i )
    {
      char buffer[64];
      snprintf( buffer, sizeof(buffer), "<tr><td>%i</td>", static_cast<int>(i) );
      msg += buffer;
      
      if( i < oldcoefs.size() )
      {
        snprintf( buffer, sizeof(buffer), "<td>%.4g</td>", oldcoefs[i] );
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
    
    const auto &olddev = oldcal->deviation_pairs();
    const auto &newdev = newcal->deviation_pairs();
    
    if( !olddev.empty() || !newdev.empty() )
    {
      msg += "<div>Old deviation pairs: [";
      for( size_t i = 0; i < olddev.size(); ++i )
        msg += (i?", {":"{") + to_string(olddev[i].first) + ", " + to_string(olddev[i].second) +"}";
      msg += "]</div>";
      
      msg += "<div>New deviation pairs: [";
      for( size_t i = 0; i < newdev.size(); ++i )
        msg += (i?", {":"{") + to_string(newdev[i].first) + ", " + to_string(newdev[i].second) +"}";
      msg += "]</div>";
    }//if( either have deviation pairs )
    
    return msg;
  };//create_table(...) lambda
  
  WContainerWidget *calcontent = new WContainerWidget( contents() );
  calcontent->addStyleClass( "EnergyCalPreserveCalContent" );
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( viewer->renderedHeight() > 250 )
    calcontent->setMaximumSize( WLength::Auto, viewer->renderedHeight() - 200 );
  else if( viewer->renderedHeight() > 10 )
    calcontent->setMaximumSize( WLength::Auto, 0.4*viewer->renderedHeight() );
  
  if( new_cals.size() == 1 )
  {
    const string &name = begin(new_cals)->first;
    const auto newcal = begin(new_cals)->second;
    const auto oldpos = old_cals.find(name);
    assert( oldpos != end(old_cals) );
    const auto oldcal = (oldpos == end(old_cals)) ? nullptr : oldpos->second;  //should always be valid
    
    msg = create_table( oldcal, newcal );
    auto txt = new WText( msg, XHTMLText, calcontent );
    txt->setInline( false );
  }else
  {
    WGridLayout *layout = new WGridLayout( calcontent );
    
    auto contentstack = new WStackedWidget();
    contentstack->addStyleClass( "EnergyCalPreserveStack" );
    WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
    contentstack->setTransitionAnimation( animation );
    
    auto menu = new WMenu( contentstack );
    menu->addStyleClass( "VerticalMenu" );
    //menu->itemSelected().connect( this, &... );
    
    layout->addWidget( menu, 0, 0 );
    layout->addWidget( contentstack, 0, 1 );
    layout->setColumnStretch( 1, 1 );
    
    for( const auto &nc : new_cals )
    {
      const string &name = nc.first;
      const auto newcal = nc.second;
      const auto old_pos = old_cals.find(name);
      assert( old_pos != end(old_cals) );
      const auto oldcal = (old_pos == end(old_cals)) ? nullptr : old_pos->second;  //should always be valid
      
      string tablestr = create_table( oldcal, newcal );
      auto tabltw = new WText( tablestr, XHTMLText );
      tabltw->setInline( false );
      auto item = menu->addItem( name, tabltw );
      
      //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
      item->clicked().connect( boost::bind(&WMenuItem::select, item) );
    }//for( const auto &nc : new_cals )
  }//if( new_cals.size == 1 ) / else
  
  
  
  msg = "<div class=\"EnergyCalPreserveBottomTxt\"><b>"
         "Would you like to use the previous calibration with the new file?"
         "</b></div>";
  
  txt = new WText( msg, XHTMLText, contents() );
  txt->setInline( false );
  
  
  WContainerWidget *bottom = footer();
  //bottom->setStyleClass("modal-footer");
  WPushButton *cancel = new WPushButton( "No", bottom );
  cancel->clicked().connect( this, &EnergyCalPreserveWindow::emitReject );
  cancel->setWidth( WLength(47.5,WLength::Percentage) );
  cancel->setFloatSide(Wt::Left);
  
  WPushButton *use = new WPushButton( "Yes", bottom );
  use->clicked().connect( this, &EnergyCalPreserveWindow::propogateCalibrations );
  use->setWidth( WLength(47.5,WLength::Percentage) );
  use->setFloatSide( Wt::Right );
  
  rejectWhenEscapePressed();
  AuxWindow::disableCollapse();
  AuxWindow::show();
  setWidth( 400.0 );
  AuxWindow::resizeToFitOnScreen();
  AuxWindow::centerWindow();
}//EnergyCalPreserveWindow


EnergyCalPreserveWindow::~EnergyCalPreserveWindow()
{
}


bool EnergyCalPreserveWindow::candidate( std::shared_ptr<SpecMeas> newmeas,
                                     std::shared_ptr<SpecMeas> oldmeas )
{
  if( !newmeas || !oldmeas || newmeas==oldmeas )
    return false;
  
  // If detectors have a different serial number, then they arent a potential candidate
  if( newmeas->instrument_id() != oldmeas->instrument_id() )
    return false;
  
  // Check if every detector in 'newmeas' with a calibration, also has a calibration in 'oldmeas'
  //  and the number of channels match
  const auto new_cals = gamma_names_to_cals( newmeas ); // TODO: should say currently displayed sample number, wo we can get energy calibration of current sample number
  const auto old_cals = gamma_names_to_cals( oldmeas );
  
  if( new_cals.empty() || old_cals.empty() || (new_cals.size() != old_cals.size()) )
    return false;
  
  bool allCalsMatch = true;
  
  for( const auto &nc : new_cals )
  {
    const auto old_pos = old_cals.find( nc.first );
    
    if( old_pos == end(old_cals) )
      return false;
    
    const auto &oldcal = old_pos->second;
    const auto &newcal = nc.second;
    
    if( oldcal->num_channels() != newcal->num_channels() )
      return false;
    
    allCalsMatch = (oldcal->type() == newcal->type());
    if( allCalsMatch )
    {
      const auto &oldcoefs = oldcal->coefficients();
      const auto &newcoefs = newcal->coefficients();
      
      allCalsMatch = (oldcoefs.size() == newcoefs.size());
      if( allCalsMatch )
      {
        for( size_t i = 0; allCalsMatch && (i < oldcoefs.size()); ++i )
        {
          const float a = oldcoefs[i], b = newcoefs[i];
          const float diff = fabs(a - b);

          // \TODO: make a better comparison method that is tuned to polynomial/FRF/LowerChannel
          allCalsMatch = ( (diff < (1.0E-5*std::max(fabs(a),fabs(b)))) || (diff < 1.0E-11) );
        }//for( size_t i = 0; i < oldcoefs.size(); ++i )
      }//if( allCalsMatch )
    }//if( allCalsMatch )
    
    if( !allCalsMatch )
      break;
    
    const auto &olddevpairs = oldcal->deviation_pairs();
    const auto &newdevpairs = newcal->deviation_pairs();
    allCalsMatch = (olddevpairs.size() == newdevpairs.size());
    if( allCalsMatch )
    {
      for( size_t i = 0; allCalsMatch && (i < olddevpairs.size()); ++i )
      {
        float a = olddevpairs[i].first;
        float b = newdevpairs[i].first;
        float diff = fabs(a - b);
        
        allCalsMatch = ( (diff < (1.0E-5*std::max(fabs(a),fabs(b)))) || (diff < 1.0E-4f) );
        
        if( allCalsMatch )
        {
          a = olddevpairs[i].second;
          b = newdevpairs[i].second;
          diff = fabs(a - b);
          
          allCalsMatch = ( (diff < (1.0E-5*std::max(fabs(a),fabs(b)))) || (diff < 1.0E-4f) );
        }//if( allCalsMatch )
      }//for( loop over deviation pairs )
    }//if( allCalsMatch )
  }//for( const auto &nc : new_cals )
  
  return !allCalsMatch;
}//candidate(...)


void EnergyCalPreserveWindow::propogateCalibrations()
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  const auto new_cals = gamma_names_to_cals( m_newmeas );
  const auto old_cals = gamma_names_to_cals( m_oldmeas );
  
  const auto dispSamples = viewer->displayedSamples(m_newtype);
  const auto dispMeas = viewer->measurment(m_newtype);
  
  map<set<int>,std::deque< std::shared_ptr<const PeakDef>>> new_peaks;
  
  const vector<string> gammaDetNames = m_newmeas->gamma_detector_names();
  const set<set<int>> samplesWithPeak = m_newmeas->sampleNumsWithPeaks();
  
  try
  {
    for( const set<int> &samples : samplesWithPeak )
    {
      auto oldPeaks = m_newmeas->peaks( samples );
      if( !oldPeaks || oldPeaks->empty() )
        continue;
      
      auto peak_cal = m_newmeas->suggested_sum_energy_calibration( samples, gammaDetNames );
      
      //We dont expect to ever actually not get a valid calibration, but JIC
      if( !peak_cal || !peak_cal->valid() )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Unexpectedly didnt get suggested energy cal for"
                            " samples/detectectso" );
#endif
        continue;
      }//if( no previous calibration )
      
      //Figure out what detector 'peak_cal' belongs too
      string detname;
      bool foundname = false;
      for( auto siter = begin(samples); !foundname && (siter != end(samples)); ++siter )
      {
        for( size_t i = 0; !foundname && (i < gammaDetNames.size()); ++i )
        {
          auto m = m_newmeas->measurement( *siter, gammaDetNames[i] );
          if( m && (m->energy_calibration() == peak_cal) )
          {
            foundname = true;
            detname = gammaDetNames[i];
          }//if( m && (m->energy_calibration() == peak_cal) )
        }//for( loop over names )
      }//for( loop over sample numbers )
      
      if( !foundname )  //probably shouldnt ever not find the calibration
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Unexpectedly couldnt find suggested energy cal" );
#endif
        continue;
      }//if( !foundname )
      
      auto old_cal_pos = old_cals.find( detname );
      if( old_cal_pos == end(old_cals) )  //shouldnt ever really happen
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Unexpectedly couldnt old calibration for detector" );
#endif
        continue;
      }//if( couldnt find old calibration )
      
      const auto newPeakCal = old_cal_pos->second;
      new_peaks[samples] = translatePeaksForCalibrationChange( *oldPeaks, peak_cal, newPeakCal );
    }//for( const set<int> &samples : samplesWithPeak )
  }catch( std::exception &e )
  {
    string msg = "Propagating old energy calibration to peaks in new file caused an error;"
                 " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#endif
            
    viewer->logMessage( msg, "", 2 );
    finished().emit(WDialog::Accepted);
    return;
  }//try / catch translate peaks
  
  
  
  try
  {
    // First loop over and make sure we can actually apply all the calibrations appropriately.
    for( const auto &m : m_newmeas->measurements() )
    {
      const auto new_file_cal = m ? m->energy_calibration() : nullptr;
      if( !new_file_cal || !new_file_cal->valid() || m->num_gamma_channels() < 5 )
        continue;
      
      const auto old_cal_pos = old_cals.find( m->detector_name() );
      if( old_cal_pos == end(old_cals) )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Unexcpected error finding calibration for a detector" );
#endif
        continue;
      }//if( old_cal_pos == end(old_cals) )
      
      const auto prev_file_cal = old_cal_pos->second;
      
      if( new_file_cal->num_channels() != prev_file_cal->num_channels() )
        throw runtime_error( "Detector '" + m->detector_name() + "' has a different number of"
                          " channels in old and new file, so not propogating energy calibration" );
    }//for( loop over measurements )
    
    
    for( const auto &m : m_newmeas->measurements() )
    {
      const auto new_file_cal = m ? m->energy_calibration() : nullptr;
      if( !new_file_cal || !new_file_cal->valid() || m->num_gamma_channels() < 5 )
        continue;
      
      const auto old_cal_pos = old_cals.find( m->detector_name() );
      if( old_cal_pos == end(old_cals) )
        continue;  //shouldnt happen
      
      const auto prev_file_cal = old_cal_pos->second;
      assert( new_file_cal->num_channels() == prev_file_cal->num_channels() );
      
      m_newmeas->set_energy_calibration( prev_file_cal, m );
    }//for( loop over measurements )
    
    
    //Now loop over and update peaks
    for( const auto &peakiter : new_peaks )
    {
      m_newmeas->setPeaks( peakiter.second, peakiter.first );
      
      //Check if we need to let the PeakModel know that dispalyed peaks have been updated.
      if( (m_newtype == SpectrumType::Foreground)
          && (dispMeas == m_newmeas)
          && (peakiter.first == dispSamples) )
      {
        viewer->peakModel()->setPeakFromSpecMeas( m_newmeas, dispSamples );
      }
    }//for( loop over new peaks mape )
  }catch( std::exception &e )
  {
    viewer->logMessage( e.what(), "", 2 );
    finished().emit(WDialog::Accepted);
    return;
  }//try / catch
  
  viewer->refreshDisplayedCharts();
  m_calibrator->refreshGuiFromFiles();
  
  finished().emit(WDialog::Accepted);
}//void propogateCalibrations()
