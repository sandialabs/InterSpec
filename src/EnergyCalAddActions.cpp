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
#include "InterSpec/PeakModel.h"
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
  AuxWindow *m_parent;
  WButtonGroup *m_group;
  WText *m_fitCoefTxt;
  
  WPushButton *m_cancel;
  WPushButton *m_accept;
  
public:
  ConvertToPolyTool( EnergyCalTool *cal, AuxWindow *parent );
  ~ConvertToPolyTool();
  
  void fromLowerChannelEnergiesOptionChanged();
  void convertLowerChannelEnegies( const size_t ncoeffs, const bool set_cals_to, std::string &msg );
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
    m_parent( parent ),
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
    }//case SpecUtils::EnergyCalType::FullRangeFraction:
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    {
      WText *msg = new WText( "How would you like the conversion to be performed?", this );
      msg->addStyleClass( "ConvertToMsg" );
      msg->setInline( false );
      
      m_group = new WButtonGroup( this );
      WGroupBox *box = new WGroupBox( this );
      
      WRadioButton *btn = new WRadioButton( "Linearize, rebinning counts to match new widths", box );
      btn->setToolTip( "The channels in the resulting spectra will all have the same width,"
                      " but counts in each channel will not be integer; this option preserves"
                      " spectral shape (counts per keV)." );
      btn->setInline( false );
      m_group->addButton( btn, 0 );
      
      btn = new WRadioButton( "Linearize, leaving counts in channels the same", box );
      btn->setToolTip( "The channel energies will be reasigned so that each channel will have"
                      " same width, without taking into account the channels original width;"
                      " i.e., the counts in each channel will remain the same, but spectral shape"
                      " (counts per keV) will be changed." );
      btn->setInline( false );
      m_group->addButton( btn, 1 );
      
      btn = new WRadioButton( "Fit channel energies to quadratic equation.", box );
      btn->setInline( false );
      m_group->addButton( btn, 2 );
      
      btn = new WRadioButton( "Fit channel energies to quartic equation.", box );
      btn->setInline( false );
      m_group->addButton( btn, 3 );
      
      m_group->setSelectedButtonIndex( 0 );
      m_group->checkedChanged().connect( this, &ConvertToPolyTool::fromLowerChannelEnergiesOptionChanged );
      
      m_group->setSelectedButtonIndex( 0 );
      WString msgstr = m_group->button(0) ? m_group->button(0)->toolTip() : WString();
      m_fitCoefTxt = new WText( msgstr, this );
      m_fitCoefTxt->setInline( false );
      m_fitCoefTxt->addStyleClass( "ConvertFromLCERes" );
      
      break;
    }//case SpecUtils::EnergyCalType::LowerChannelEdge:
  }//switch( origtype )
  
  
  WContainerWidget *buttonDiv = nullptr;
  if( parent )
    buttonDiv = parent->footer();
  else
    buttonDiv = new WContainerWidget( this );
  
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
  

void ConvertToPolyTool::convertLowerChannelEnegies( const size_t ncoeffs,
                                                    const bool set_cals, string &msgsumm )
{
  using SpecUtils::EnergyCalibration;
  
  msgsumm.clear();
  set<shared_ptr<SpecMeas>> shiftedPeaksFor;
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
  map<pair<shared_ptr<SpecMeas>,set<int>>, deque<shared_ptr<const PeakDef>>> updated_peaks;
  
  vector<MeasToApplyCoefChangeTo> applicables = m_cal->measurementsToApplyCoeffChangeTo();
  for( const auto &toapplyto : applicables )
  {
    for( const auto &detname : toapplyto.detectors )
    {
      for( const int sample : toapplyto.sample_numbers )
      {
        auto m = toapplyto.meas->measurement( sample, detname );
        auto cal = m ? m->energy_calibration() : nullptr;
        if( !cal
           || !cal->valid()
           || (cal->num_channels() < 5)
           || (cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
          continue;
        
        auto pos = updated_cals.find( cal );
        if( pos != end(updated_cals) )
          continue;
        
        const auto &lower_energies = *cal->channel_energies();
        const float chnlwidth = (cal->upper_energy() - cal->lower_energy()) / cal->num_channels();
        
        vector<float> coefs;
        const double avrgdiff = EnergyCal::fit_poly_from_lower_channel_energies( ncoeffs,
                                                                          lower_energies, coefs );
        
        auto newcal = make_shared<EnergyCalibration>();
        newcal->set_polynomial( cal->num_channels(), coefs, {} );
        updated_cals[cal] = newcal;
        
        if( msgsumm.empty() )
        {
          char buffer[32];
          msgsumm = "<p><div>Fit polynomial coefficients:</div><div>&nbsp;&nbsp;&nbsp;&nbsp;";
          for( size_t i = 0; i < coefs.size(); ++i )
          {
            snprintf( buffer, sizeof(buffer), "%s%1.5g", (i ? ", " : ""), coefs[i] );
            msgsumm += buffer;
          }
          msgsumm += "</div></p>";
          
          snprintf( buffer, sizeof(buffer), "%1.3g", avrgdiff );
          msgsumm += "<p>This has an average error of " + string(buffer) + " keV in each channel, which is ";
          snprintf( buffer, sizeof(buffer), "%1.3g", (100.0*avrgdiff/chnlwidth) );
          msgsumm += buffer;
          msgsumm += "% of the channels width.</p>";
          if( avrgdiff > 0.2*chnlwidth )
            msgsumm += "<p>The input energy calibration doesnt look amenable to convert to a"
            " polynomial; recomend to linearize with rebining.</p>";
        }//if( msgsumm.empty() )
      }//for( const int sample : toapplyto.sample_numbers )
    }//for( const auto &detname : toapplyto.detectors )
    
    
    if( !shiftedPeaksFor.count( toapplyto.meas) )
    {
      const auto gammaDetNames = toapplyto.meas->gamma_detector_names();
      const auto samplesWithPeak = toapplyto.meas->sampleNumsWithPeaks();
      for( const set<int> &samples : samplesWithPeak )
      {
        auto oldPeaks = toapplyto.meas->peaks( samples );
        if( !oldPeaks || oldPeaks->empty() )
            continue;
            
        auto peak_cal = toapplyto.meas->suggested_sum_energy_calibration( samples, gammaDetNames );
            
        //We dont expect to ever actually not get a valid calibration, but JIC
        if( !peak_cal || !peak_cal->valid() )
        {
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, "Unexpectedly didnt get suggested energy cal for"
                                         " samples/detectector" );
#endif
          continue;
        }//if( no previous calibration )
            
        auto pos = updated_cals.find(peak_cal);
        if( pos == end(updated_cals) )  //shouldnt ever really happen
        {
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, "Unexpectedly couldnt find cal for peaks" );
#endif
          continue;
        }//if( couldnt find old calibration )
          
        const auto newPeakCal = pos->second;
        try
        {
          updated_peaks[{toapplyto.meas,samples}]
                 = EnergyCal::translatePeaksForCalibrationChange( *oldPeaks, peak_cal, newPeakCal );
        }catch( std::exception &e )
        {
          throw runtime_error( "Error translating peaks for calibration change: "
                               + string(e.what()) );
        }//try / catch translate peaks
      }//for( const set<int> &samples : samplesWithPeak )
    }//if( !shiftedPeaksFor.count( toapplyto.meas) )
  }//for( const auto &toapplyto : applicables )
  
  if( !set_cals )
    return;
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  PeakModel *peakmodel = viewer->peakModel();
  assert( peakmodel );
  
  const auto foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const auto &foreSamples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  
  for( const auto &mp : updated_peaks )
  {
    const auto &spec = mp.first.first;
    const auto &samples = mp.first.second;
    const auto &newpeaks = mp.second;
    spec->setPeaks( newpeaks, samples );
    
    if( (spec == foreground) && (samples == foreSamples) )
      peakmodel->setPeakFromSpecMeas( spec, samples);
  }//for( const auto &mp : updated_peaks )
  
  
  for( const auto &toapplyto : applicables )
  {
    for( const auto &detname : toapplyto.detectors )
    {
      for( const int sample : toapplyto.sample_numbers )
      {
        auto m = toapplyto.meas->measurement( sample, detname );
        auto cal = m ? m->energy_calibration() : nullptr;
        if( !cal
           || !cal->valid()
           || (cal->num_channels() < 5)
           || (cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
          continue;
        
        auto pos = updated_cals.find( cal );
        if( pos != end(updated_cals) )
          toapplyto.meas->set_energy_calibration( pos->second, m );
      }//for( loop over samples )
    }//for( loop over detectors )
  }//for( loop over applicables )
}//void ConvertToPolyTool::convertLowerChannelEnegies(...)


void ConvertToPolyTool::fromLowerChannelEnergiesOptionChanged()
{
  assert( m_group );
  assert( m_fitCoefTxt );
  
  const int index = m_group->selectedButtonIndex();
  
  try
  {
    m_accept->enable();
    
    switch( index )
    {
      case 0: //Linearize, rebinning counts to match
      case 1: //Linearize, leaving counts in channels the
        if( m_group->checkedButton() )
          m_fitCoefTxt->setText( m_group->checkedButton()->toolTip() );
        else
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
    convertLowerChannelEnegies( ncoeffs, false, msgsumm );
    
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
  using namespace SpecUtils;
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  switch( result )
  {
    case WDialog::Rejected:
      cerr << "\nRejected ConvertToPolyTool" << endl;
    break;
      
    case WDialog::Accepted:
    {
      const int index = m_group ? m_group->selectedButtonIndex() : -1;
      
      if( index < 2 )
      {
        //Convert from FRF, or linearize spectra that is currently defined by lower energies
        map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
        vector<MeasToApplyCoefChangeTo> applicables = m_cal->measurementsToApplyCoeffChangeTo();
        for( const auto &toapplyto : applicables )
        {
          for( const auto &detname : toapplyto.detectors )
          {
            for( const int sample : toapplyto.sample_numbers )
            {
              auto m = toapplyto.meas->measurement( sample, detname );
              auto cal = m ? m->energy_calibration() : nullptr;
              if( !cal
                 || !cal->valid()
                 || (cal->num_channels() < 5)
                 || (cal->type() != EnergyCalType::LowerChannelEdge) )
                continue;
              
              const size_t nchannel = cal->num_channels();
              
              auto pos = updated_cals.find( cal );
              if( pos == end(updated_cals) )
              {
                vector<float> newcoefs;
                const auto &oldcoefs = cal->coefficients();
                switch( cal->type() )
                {
                  case EnergyCalType::Polynomial:
                  case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                    assert( 0 );
                    continue;
                    
                  case EnergyCalType::FullRangeFraction:
                    assert( index == -1 );
                    newcoefs = fullrangefraction_coef_to_polynomial( oldcoefs, nchannel );
                    break;
                    
                  case EnergyCalType::LowerChannelEdge:
                    assert( index == 0 || index == 1 );
                    newcoefs = { cal->lower_energy(),
                                 (cal->upper_energy() - cal->lower_energy())/nchannel };
                    break;
                    
                  case EnergyCalType::InvalidEquationType:
                    assert( 0 );
                    continue; //
                    break;
                }//switch( cal->type() )
                
                auto newcal = make_shared<EnergyCalibration>();
                newcal->set_polynomial( nchannel, newcoefs, cal->deviation_pairs() );
                pos = updated_cals.insert( {cal,newcal} ).first;
              }//if( we need to create a new calibration )
              
              if( index == -1 || index == 1 )
              {
                toapplyto.meas->set_energy_calibration( pos->second, m );
              }else if( index == 0 )
              {
                toapplyto.meas->rebin_measurement( pos->second, m );
              }else
              {
                assert( 0 );
              }
            }//for( loop over sample numbers )
          }//for( loop over detector names )
          
          if( index == 1 )
          {
            // \TODO: should create new peaks here since widths and stuff have changed
          }//if( index == 1 )
          
        }//for( const auto &toapplyto : applicables )
      }else
      {
        assert( index == 2 || index == 3 );
        string dummy;
        const size_t ncoeffs = static_cast<size_t>(index+1);
        convertLowerChannelEnegies( ncoeffs, true, dummy );
      }//if( convert from FRF or linearize ) / else ( fit polynomial from lower energies )
      
      cerr << "\nAccepted ConvertToPolyTool" << endl;
      m_cal->refreshGuiFromFiles();
      viewer->refreshDisplayedCharts();
      
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  if( m_parent )
    AuxWindow::deleteAuxWindow( m_parent );
}//void handleFinish(...)

