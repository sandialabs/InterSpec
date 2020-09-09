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
#include <tuple>
#include <vector>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WSpinBox>
#include <Wt/WGroupBox>
#include <Wt/WTableCell>
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
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/EnergyCalAddActions.h"

// Forward declarations
class ConvertCalTypeTool;

using namespace std;
using namespace Wt;


class ConvertCalTypeTool : public WContainerWidget
{
  /** The type of calibration to convert to */
  const SpecUtils::EnergyCalType m_targetType;
 
  shared_ptr<vector<MeasToApplyCoefChangeTo>> m_measToChange;
  
  /** The type of calibration we are converting from.
   If input SpecFiles has more than one type than this will be
   #SpecUtils::EnergyCalType::InvalidEquationType and the ConvertCalTypeTool constuctor will not
   allow conversion.
   */
  SpecUtils::EnergyCalType m_sourceType;
  
  EnergyCalTool * const m_cal;
  AuxWindow * const m_parent;
  
  /** Will be non-nullptr only if m_sourceType==SpecUtils::EnergyCalType::LowerChannelEdge: */
  WButtonGroup *m_group;
  WText *m_fitCoefTxt;
  
  WPushButton *m_cancel;
  WPushButton *m_accept;
  
public:
  ConvertCalTypeTool( const SpecUtils::EnergyCalType target,
                      shared_ptr<vector<MeasToApplyCoefChangeTo>> measToChange,
                      EnergyCalTool *cal, AuxWindow *parent );
  virtual ~ConvertCalTypeTool();
  
  void fromLowerChannelEnergiesOptionChanged();
  void convertLowerChannelEnegies( const size_t ncoeffs, const bool set_cals_to, std::string &msg );
  void handleFinish( Wt::WDialog::DialogCode result );
};//class ConvertCalTypeTool



class LinearizeCalTool : public WContainerWidget
{
  shared_ptr<vector<MeasToApplyCoefChangeTo>> m_measToChange;
  
  EnergyCalTool * const m_cal;
  AuxWindow * const m_parent;
  
  /** Number of orignal channels, lower energy, upper energy, number of channels. */
  map<size_t,tuple<NativeFloatSpinBox *,NativeFloatSpinBox *,WSpinBox *>> m_options;
  
  WPushButton *m_cancel;
  WPushButton *m_accept;
  
  bool checkValidInputs();
  void energyRangeUpdatedCallback();
  void numberOfChannelsUpdatedCallback();
public:
  LinearizeCalTool( shared_ptr<vector<MeasToApplyCoefChangeTo>> measToChange, EnergyCalTool *cal, AuxWindow *parent );
  virtual ~LinearizeCalTool();
  
  void handleFinish( Wt::WDialog::DialogCode result );
};//class LinearizeCalTool




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
      AuxWindow::setWindowTitle( "Linearize Spectrum" );
      new LinearizeCalTool( m_measToChange, m_calibrator, this );
      break;
      
    case MoreActionsIndex::Truncate:
      break;
      
    case MoreActionsIndex::CombineChannels:
      break;
      
    case MoreActionsIndex::ConvertToFrf:
      AuxWindow::setWindowTitle( "Convert to Full Range Fraction Calibration" );
      new ConvertCalTypeTool( SpecUtils::EnergyCalType::FullRangeFraction, m_measToChange, m_calibrator, this );
      break;
      
    case MoreActionsIndex::ConvertToPoly:
      AuxWindow::setWindowTitle( "Convert to Polynomial Calibration" );
      new ConvertCalTypeTool( SpecUtils::EnergyCalType::Polynomial, m_measToChange, m_calibrator, this );
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
 





ConvertCalTypeTool::ConvertCalTypeTool( const SpecUtils::EnergyCalType targetType,
                                       shared_ptr<vector<MeasToApplyCoefChangeTo>> measToChange,
                                      EnergyCalTool *cal, AuxWindow *parent )
  : WContainerWidget(),
    m_targetType( targetType ),
    m_measToChange( measToChange ),
    m_sourceType( SpecUtils::EnergyCalType::InvalidEquationType ),
    m_cal( cal ),
    m_parent( parent ),
    m_group( nullptr ),
    m_fitCoefTxt( nullptr ),
    m_cancel( nullptr ),
    m_accept( nullptr )
{
  using SpecUtils::EnergyCalType;
  
  assert( (m_targetType == SpecUtils::EnergyCalType::FullRangeFraction)
          || (m_targetType == SpecUtils::EnergyCalType::Polynomial)
          //LowerChannelEdge is untested as of yet, but should work in principle, but why would you do this?
          //|| (m_targetType == SpecUtils::EnergyCalType::LowerChannelEdge)
         );
  
  addStyleClass( "ConvertCalTypeTool" );
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  if( parent )
    parent->stretcher()->addWidget( this, 0, 0  );
  
  set<size_t> coeforder, nchannels;
  set<SpecUtils::EnergyCalType> caltypes;
  if( m_measToChange )
  {
    for( const auto &toapplyto : *m_measToChange )
    {
      for( const auto &detname : toapplyto.detectors )
      {
        for( const int sample : toapplyto.sample_numbers )
        {
          auto m = toapplyto.meas->measurement( sample, detname );
          if( !m || !m->energy_calibration() || !m->energy_calibration()->valid()
             || (m->num_gamma_channels() < 5) )
            continue;
          
          auto type = m->energy_calibration()->type();
          if( type == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
            type = SpecUtils::EnergyCalType::Polynomial;
          
          nchannels.insert( m->num_gamma_channels() );
          caltypes.insert( type );
          coeforder.insert( m->calibration_coeffs().size() );
        }//for( loop over samples )
      }//for( loop over detectors )
    }//for( const auto &m : applicables )
  }//if( m_measToChange )
  
  if( (caltypes.size() != 1) || ((*begin(caltypes)) == m_targetType) )
  {
    if( parent )
    {
      auto close = parent->addCloseButtonToFooter();
      close->clicked().connect( boost::bind( &ConvertCalTypeTool::handleFinish, this, WDialog::Rejected ) );
      parent->finished().connect( this, &ConvertCalTypeTool::handleFinish );
    }
    
    string msgtxt = "";
    if( caltypes.empty() )
      msgtxt = "Conversion would not effect any spectra; try selecting more"
               " &quot;Apply Changes To&quot; criteria.";
    else if( caltypes.size() == 1 )
      msgtxt = "All spectra this would be applied to is already "
               + string( m_targetType==EnergyCalType::Polynomial ? "polynomial." : "full range fraction.");
    else
      msgtxt = "There is more than one calibration type of selected source spectra;"
               " try restricting &quot;Apply Changes To&quot; criteria.";
    
    WText *msg = new WText( WString::fromUTF8(msgtxt), this );
    msg->addStyleClass( "ConvertToNA" );
    msg->setInline( false );
    
    return;
  }//if( caltypes.size() != 1 || cal is type wanted )
  
  
  string applyToTxt = cal->applyToSummaryTxt();
  if( !applyToTxt.empty() )
  {
    applyToTxt = "Changes will be applied to: " + applyToTxt;
    WText *msg = new WText( applyToTxt, this );
    msg->addStyleClass( "ConvertToApplieTo" );
    msg->setInline( false );
  }//if( !applyToTxt.empty() )
  
  
  m_sourceType = *begin(caltypes);
  switch( m_sourceType )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
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
      m_group->checkedChanged().connect( this, &ConvertCalTypeTool::fromLowerChannelEnergiesOptionChanged );
      
      m_group->setSelectedButtonIndex( 0 );
      WString msgstr = m_group->button(0) ? m_group->button(0)->toolTip() : WString();
      m_fitCoefTxt = new WText( msgstr, this );
      m_fitCoefTxt->setInline( false );
      m_fitCoefTxt->addStyleClass( "ConvertFromLCERes" );
      
      break;
    }//case SpecUtils::EnergyCalType::LowerChannelEdge:
      
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
    break;
  }//switch( m_sourceType )
  
  
  WContainerWidget *buttonDiv = nullptr;
  if( parent )
    buttonDiv = parent->footer();
  else
    buttonDiv = new WContainerWidget( this );
  
  //AuxWindow::addHelpInFooter( buttonDiv, "convert-to-polynomial-dialog" );
  
  m_cancel = new WPushButton( "Cancel", buttonDiv );
  m_accept = new WPushButton( "Accept", buttonDiv );
  
  m_cancel->clicked().connect( boost::bind( &ConvertCalTypeTool::handleFinish, this, WDialog::Rejected ) );
  m_accept->clicked().connect( boost::bind( &ConvertCalTypeTool::handleFinish, this, WDialog::Accepted ) );
  
  if( parent )
  {
    parent->finished().connect( this, &ConvertCalTypeTool::handleFinish );
    
    const int w = 600 < viewer->renderedWidth() ? 600 : viewer->renderedWidth();
    //const int h = static_cast<int>(0.8*viewer->renderedHeight());
    //parent->resizeWindow( w, h );
    parent->setWidth( w );
    
    parent->rejectWhenEscapePressed();
    
    parent->centerWindow();
  }//if( parent )
}//ConvertCalTypeTool(...)
  

ConvertCalTypeTool::~ConvertCalTypeTool()
{
}
  

void ConvertCalTypeTool::convertLowerChannelEnegies( const size_t ncoeffs,
                                                    const bool set_cals, string &msgsumm )
{
  using SpecUtils::EnergyCalibration;
  if( !m_measToChange )
  {
    assert( 0 );
    return;
  }
  
  msgsumm.clear();
  set<shared_ptr<SpecMeas>> shiftedPeaksFor;
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
  map<pair<shared_ptr<SpecMeas>,set<int>>, deque<shared_ptr<const PeakDef>>> updated_peaks;
  
  for( const auto &toapplyto : *m_measToChange )
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
        double avrgdiff = -999;
        
        auto newcal = make_shared<EnergyCalibration>();
        switch( m_targetType )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            avrgdiff = EnergyCal::fit_poly_from_channel_energies( ncoeffs,
                                                                        lower_energies, coefs );
            newcal->set_polynomial( cal->num_channels(), coefs, {} );
            break;
            
          case SpecUtils::EnergyCalType::FullRangeFraction:
            avrgdiff = EnergyCal::fit_full_range_fraction_from_channel_energies( ncoeffs,
                                                lower_energies, coefs );
            newcal->set_full_range_fraction( cal->num_channels(), coefs, {} );
            break;
            
          case SpecUtils::EnergyCalType::LowerChannelEdge:
          case SpecUtils::EnergyCalType::InvalidEquationType:
            assert( 0 );
            break;
        }//switch( m_targetType )
        
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
  
  
  for( const auto &toapplyto : *m_measToChange )
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
}//void ConvertCalTypeTool::convertLowerChannelEnegies(...)


void ConvertCalTypeTool::fromLowerChannelEnergiesOptionChanged()
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


void ConvertCalTypeTool::handleFinish( Wt::WDialog::DialogCode result )
{
  using namespace SpecUtils;
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  switch( result )
  {
    case WDialog::Rejected:
      cerr << "\nRejected ConvertCalTypeTool" << endl;
    break;
      
    case WDialog::Accepted:
    {
      const int index = m_group ? m_group->selectedButtonIndex() : -1;
      
      if( !m_measToChange )
      {
        assert( 0 );
        break;
      }
      
      if( index < 2 )
      {
        //Convert from FRF, or linearize spectra that is currently defined by lower energies
        map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
        for( const auto &toapplyto : *m_measToChange )
        {
          for( const auto &detname : toapplyto.detectors )
          {
            for( const int sample : toapplyto.sample_numbers )
            {
              auto m = toapplyto.meas->measurement( sample, detname );
              auto cal = m ? m->energy_calibration() : nullptr;
              if( !cal || !cal->valid() || (cal->num_channels() < 5) )
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
                  {
                    switch( m_targetType )
                    {
                      case SpecUtils::EnergyCalType::Polynomial:
                      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                        continue;
                        
                      case SpecUtils::EnergyCalType::FullRangeFraction:
                        newcoefs = polynomial_coef_to_fullrangefraction( oldcoefs, nchannel );
                        break;
                        
                      case SpecUtils::EnergyCalType::LowerChannelEdge:
                        newcoefs = *cal->channel_energies();
                        break;
                        
                      case SpecUtils::EnergyCalType::InvalidEquationType:
                        assert( 0 );
                        continue;
                    }//switch( m_targetType )
                    break;
                  }// cal->type() == EnergyCalType::Polynomial
                    
                  case EnergyCalType::FullRangeFraction:
                  {
                    switch( m_targetType )
                    {
                      case SpecUtils::EnergyCalType::Polynomial:
                      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                        newcoefs = fullrangefraction_coef_to_polynomial( oldcoefs, nchannel );
                        break;
                        
                      case SpecUtils::EnergyCalType::FullRangeFraction:
                        continue;
                        
                      case SpecUtils::EnergyCalType::LowerChannelEdge:
                        newcoefs = *cal->channel_energies();
                        break;
                        
                      case SpecUtils::EnergyCalType::InvalidEquationType:
                        assert(0);
                        continue;
                    }//switch( m_targetType )
                    assert( index == -1 );
                    
                    break;
                  }// cal->type() == EnergyCalType::FullRangeFraction
                    
                  case EnergyCalType::LowerChannelEdge:
                  {
                    assert( index == 0 || index == 1 );
                    switch( m_targetType )
                    {
                      case SpecUtils::EnergyCalType::Polynomial:
                      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                        newcoefs = { cal->lower_energy(),
                                     (cal->upper_energy() - cal->lower_energy())/nchannel };
                        break;
                        
                      case SpecUtils::EnergyCalType::FullRangeFraction:
                        newcoefs = { cal->lower_energy(),
                                     (cal->upper_energy() - cal->lower_energy()) };
                        break;
                        
                      case SpecUtils::EnergyCalType::LowerChannelEdge:
                        newcoefs = *cal->channel_energies();
                        break;
                        
                      case SpecUtils::EnergyCalType::InvalidEquationType:
                        assert(0);
                        continue;
                    }//switch( m_targetType )
                    break;
                  }// cal->type() == EnergyCalType::LowerChannelEdge:
                    
                  case EnergyCalType::InvalidEquationType:
                    assert( 0 );
                    continue; //
                    break;
                }//switch( cal->type() )
                
                auto newcal = make_shared<EnergyCalibration>();
                switch( m_targetType )
                {
                  case SpecUtils::EnergyCalType::Polynomial:
                  case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                    newcal->set_polynomial( nchannel, newcoefs, cal->deviation_pairs() );
                    break;
                  
                  case SpecUtils::EnergyCalType::FullRangeFraction:
                    newcal->set_full_range_fraction( nchannel, newcoefs, cal->deviation_pairs() );
                    break;
                    
                  case SpecUtils::EnergyCalType::LowerChannelEdge:
                    newcal->set_lower_channel_energy( nchannel, std::move(newcoefs) );
                    break;
                    
                  case SpecUtils::EnergyCalType::InvalidEquationType:
                    assert( 0 );
                    break;
                }//switch( m_targetType )
                    
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
      }//if( convert from FRF or linearize ) / else ( fit polynomial/FRF from lower energies )
      
      cerr << "\nAccepted ConvertCalTypeTool" << endl;
      m_cal->refreshGuiFromFiles();
      viewer->refreshDisplayedCharts();
      
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  if( m_parent )
    AuxWindow::deleteAuxWindow( m_parent );
}//void handleFinish(...)




LinearizeCalTool::LinearizeCalTool( shared_ptr<vector<MeasToApplyCoefChangeTo>> measToChange,
                                    EnergyCalTool *cal, AuxWindow *parent )
  : WContainerWidget(),
    m_measToChange( measToChange ),
    m_cal( cal ),
    m_parent( parent ),
    m_options{},
    m_cancel( nullptr ),
    m_accept( nullptr )
{
  using SpecUtils::EnergyCalType;
        
  addStyleClass( "LinearizeCalTool" );
    
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
    
  if( parent )
    parent->stretcher()->addWidget( this, 0, 0  );
    
  //Create a mapping that for each unique number of channels, gives the min/max energy of any
  //  measurements with that energy range
  map<size_t,std::pair<float,float>> nchannelToRange;
  if( m_measToChange )
  {
    for( const auto &toapplyto : *m_measToChange )
    {
      for( const auto &detname : toapplyto.detectors )
      {
        for( const int sample : toapplyto.sample_numbers )
        {
          auto m = toapplyto.meas->measurement( sample, detname );
          auto cal = m ? m->energy_calibration() : nullptr;
          const size_t nchannel = cal ? cal->num_channels() : size_t(0);
          if( !cal || !cal->valid() || (nchannel < 5) )
            continue;
          
          const float lower = cal->lower_energy();
          const float upper = cal->upper_energy();
          
          const auto pos = nchannelToRange.find(nchannel);
          if( pos == end(nchannelToRange) )
          {
            nchannelToRange.insert( {nchannel,{lower,upper}} );
          }else
          {
            pos->second.first = std::min( pos->second.first, lower );
            pos->second.second = std::max( pos->second.second, upper );
          }
        }//for( loop over samples )
      }//for( loop over detectors )
    }//for( const auto &m : applicables )
  }//if( m_measToChange )
  
  if( nchannelToRange.empty() )
  {
    if( parent )
    {
      auto close = parent->addCloseButtonToFooter();
      close->clicked().connect( boost::bind( &LinearizeCalTool::handleFinish, this, WDialog::Rejected ) );
      parent->finished().connect( this, &LinearizeCalTool::handleFinish );
    }
    
    const string msgtxt = "Conversion would not effect any spectra; try selecting more"
    " &quot;Apply Changes To&quot; criteria.";
    
    WText *msg = new WText( WString::fromUTF8(msgtxt), this );
    msg->addStyleClass( "ConvertToNA" );
    msg->setInline( false );
    
    return;
  }//if( nchannelToRange.empty() )
  
  
  //Create widgets
  for( const auto &chanToRange : nchannelToRange )
  {
    WTable *table = new WTable( this );
    table->addStyleClass( "LinearizeInputTbl" );
    
    const size_t nchannel = chanToRange.first;
    const float lower_energy = chanToRange.second.first;
    const float upper_energy = chanToRange.second.second;
    
    const int multiple = (nchannelToRange.size() > 1);
    if( multiple )
    {
      auto cell = table->elementAt(0, 0);
      cell->setColumnSpan( 2 );
      string msg = "For " + to_string(nchannel) + " channel spectra:";
      WText *txt = new WText( WString::fromUTF8(msg), cell );
      txt->addStyleClass( "NChannelSpecTxt" );
    }//if( multiple )
    
    const int input_nchar = 8;
    const int input_width_px = 75;
    
    auto cell = table->elementAt( 0 + multiple, 0 );
    auto label = new WLabel( "Lower Energy", cell );
    label->addStyleClass( "LinearizeLabel" );
    
    cell = table->elementAt( 0 + multiple, 1 );
    NativeFloatSpinBox *lower = new NativeFloatSpinBox( cell );
    lower->setValue( lower_energy );
    lower->setRange( -3000.0f, upper_energy );
    lower->setTextSize( input_nchar ); //doesnt seem to have an effect, but is showing up in the DOM
    lower->setWidth( input_width_px );
    lower->valueChanged().connect( this, &LinearizeCalTool::energyRangeUpdatedCallback );
    label->setBuddy( lower );
    
    cell = table->elementAt( 1 + multiple, 0 );
    label = new WLabel( "Upper Energy", cell );
    label->addStyleClass( "LinearizeLabel" );
    
    cell = table->elementAt( 1 + multiple, 1 );
    NativeFloatSpinBox *upper = new NativeFloatSpinBox( cell );
    upper->setValue( upper_energy );
    upper->setRange( lower_energy, 32000000.0f );
    upper->setTextSize( input_nchar );  //doesnt seem to have an effect, but is showing up in the DOM
    upper->setWidth( input_width_px );
    upper->valueChanged().connect( this, &LinearizeCalTool::energyRangeUpdatedCallback );
    label->setBuddy( upper );
    
    
    cell = table->elementAt( 2 + multiple, 0);
    label = new WLabel( "Num Channels", cell );
    label->addStyleClass( "LinearizeLabel" );
    
    cell = table->elementAt( 2 + multiple, 1 );
    WSpinBox *nchannelSB = new WSpinBox( cell );
    nchannelSB->setMinimum( 5 );
    nchannelSB->setMaximum( 65536 );
    nchannelSB->setValue( static_cast<int>(nchannel) );
    //nchannelSB->setNativeControl( true );
    //nchannelSB->removeStyleClass( "Wt-spinbox" ); //
    nchannelSB->setTextSize( input_nchar );
    nchannelSB->setWidth( input_width_px - 16 );
    nchannelSB->valueChanged().connect( this, &LinearizeCalTool::numberOfChannelsUpdatedCallback );
    label->setBuddy( nchannelSB );
    
    m_options[nchannel] = { lower, upper, nchannelSB };
  }//for( const auto &ranges : nchannelToRange )
  
  string applyToTxt = cal->applyToSummaryTxt();
  if( !applyToTxt.empty() )
  {
    applyToTxt = "Changes will be applied to: " + applyToTxt;
    WText *msg = new WText( applyToTxt, this );
    msg->addStyleClass( "ConvertToApplieTo" );
    msg->setInline( false );
  }//if( !applyToTxt.empty() )
  
  WContainerWidget *buttonDiv = nullptr;
  if( parent )
    buttonDiv = parent->footer();
  else
    buttonDiv = new WContainerWidget( this );
  
  //AuxWindow::addHelpInFooter( buttonDiv, "linearize-dialog" );
  
  m_cancel = new WPushButton( "Cancel", buttonDiv );
  m_accept = new WPushButton( "Accept", buttonDiv );
  
  m_cancel->clicked().connect( boost::bind( &LinearizeCalTool::handleFinish, this, WDialog::Rejected ) );
  m_accept->clicked().connect( boost::bind( &LinearizeCalTool::handleFinish, this, WDialog::Accepted ) );
  
  if( parent )
  {
    parent->finished().connect( this, &LinearizeCalTool::handleFinish );
    
    const int w = 600 < viewer->renderedWidth() ? 600 : viewer->renderedWidth();
    //const int h = static_cast<int>(0.8*viewer->renderedHeight());
    //parent->resizeWindow( w, h );
    parent->setWidth( w );
    
    parent->rejectWhenEscapePressed();
    
    parent->centerWindow();
  }//if( parent )
}//ConvertCalTypeTool(...)
  
  
LinearizeCalTool::~LinearizeCalTool()
{
    
}


bool LinearizeCalTool::checkValidInputs()
{
  // \TODO: indicate to user what is wrong
  bool allValid = true;
  
  for( auto &m : m_options )
  {
    const size_t nOrigChann = m.first;
    const float lower = get<0>(m.second)->value();
    const float upper = get<1>(m.second)->value();
    const size_t nCurrentChan = get<2>(m.second)->value();
    
    if( nCurrentChan < 5 || lower >= upper )
    {
      allValid = false;
      
      auto viewer = InterSpec::instance();
      if( viewer )
      {
        string msg = "The input for spectra with " + to_string(nOrigChann)
                     + " channels are invalid.  You must have at least 6 channels and upper energy"
                       " must be greater than lower energy";
        viewer->logMessage( WString::fromUTF8(msg), "", 2 );
      }
    }//if( nCurrentChan < 5 || lower >= upper )
  }//for( auto &m : m_options )
  
  return allValid;
}//void checkValidInputs()


void LinearizeCalTool::energyRangeUpdatedCallback()
{
  m_accept->setEnabled( checkValidInputs() );
}//void energyRangeUpdatedCallback()
  
  
void LinearizeCalTool::numberOfChannelsUpdatedCallback()
{
  m_accept->setEnabled( checkValidInputs() );
}//void numberOfChannelsUpdatedCallback()
  
  
void LinearizeCalTool::handleFinish( Wt::WDialog::DialogCode result )
{
  using namespace SpecUtils;
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  switch( result )
  {
    case WDialog::Rejected:
      cerr << "\nRejected ConvertCalTypeTool" << endl;
      break;
      
    case WDialog::Accepted:
    {
      cout << "\nAccepted LinearizeCalTool" << endl;
      
      if( !m_measToChange )
      {
        assert( 0 );
        break;
      }
      
      if( !checkValidInputs() || !m_measToChange )
        break;  //shouldnt ever happen
      
      map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
      for( const auto &toapplyto : *m_measToChange )
      {
        for( const auto &detname : toapplyto.detectors )
        {
          for( const int sample : toapplyto.sample_numbers )
          {
            auto m = toapplyto.meas->measurement( sample, detname );
            auto cal = m ? m->energy_calibration() : nullptr;
            const size_t nchannel = cal ? cal->num_channels() : size_t(0);
            if( !cal || !cal->valid() || (nchannel < 5) )
              continue;
            
            auto pos = updated_cals.find(cal);
            if( pos == end(updated_cals) )
            {
              if( !m_options.count(nchannel) )
              {
#if( PERFORM_DEVELOPER_CHECKS )
                log_developer_error( __func__, "Unexpected failure to find linearization info for"
                                               " a nchhannel" );
#endif
                continue;
              }//if( !m_options.count(nchannel) )
              
              const size_t nchan = static_cast<size_t>( get<2>(m_options[nchannel])->value() );
              const float lower = get<0>(m_options[nchannel])->value();
              const float upper = get<1>(m_options[nchannel])->value();
              
              auto newcal = make_shared<EnergyCalibration>();
              
              try
              {
                newcal->set_polynomial( nchan, {lower,(upper-lower)/nchan}, {} );
              }catch( std::exception &e ) //Shouldnt ever actually happen, but JIC
              {
                if( viewer )
                  viewer->logMessage( "Error creating new cal: " + string(e.what()), "", 3 );
                continue;
              }//try / catch to create new calibration
              
              pos = updated_cals.insert( {cal,newcal} ).first;
            }//if( we need to create a new cal )
            
            try
            {
              toapplyto.meas->rebin_measurement( pos->second, m );
            }catch( std::exception &e )  //Shouldnt ever actually happen, but JIC
            {
              if( viewer )
                viewer->logMessage( "Error rebinning whiel linearizing: " + string(e.what()), "", 3 );
            }//try / catch to create new calibration
          }//for( const int sample : toapplyto.sample_numbers )
        }//for( const auto &detname : toapplyto.detectors )
      }//for( const auto &toapplyto : *m_measToChange )
      
      m_cal->refreshGuiFromFiles();
      viewer->refreshDisplayedCharts();
      
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  if( m_parent )
    AuxWindow::deleteAuxWindow( m_parent );
}//void handleFinish( WDialog::DialogCode result )


