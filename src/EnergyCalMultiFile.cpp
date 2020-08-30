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
#include <deque>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WGroupBox>
#include <Wt/WTextArea>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>


#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ReactionGamma.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/EnergyCalMultiFile.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"

using namespace Wt;
using namespace std;

using SpecUtils::Measurement;
using SpecUtils::SpectrumType;

namespace
{
  const size_t ns_min_num_coef = 4;
}


/*

EnergyCalMultiFile::EnergyCalMultiFile( EnergyCalTool *cal, AuxWindow *parent )
: WContainerWidget(),
  m_calibrator( cal ),
  m_model( nullptr ),
  m_fitFor(),
  m_coefvals(),
  m_use( nullptr ),
  m_cancel( nullptr ),
  m_fit( nullptr ),
  m_fitSumary( nullptr ),
  m_calVal(),
  m_calUncert()
{
  InterSpec *viewer = InterSpec::instance();
 
  if( parent )
    parent->stretcher()->addWidget( this, 0, 0  );
  
  m_model = new EnergyCalMultiFileModel( cal, this );

  RowStretchTreeView *tree = new RowStretchTreeView();
  tree->setSortingEnabled(false);
  tree->setModel( m_model );
  tree->setColumn1Fixed( false );
  tree->setHeaderHeight( 2 );
  tree->setColumnWidth( 0, 200 );
  tree->setColumnWidth( 1, 150 );
  tree->setColumnWidth( 2, 150 );
  
  WGroupBox *fitFor = new WGroupBox( "Coefficents to fit for" );
  WGridLayout *fitForLayout = new WGridLayout( fitFor );
  
  
  for( size_t i = 0; i < ns_min_num_coef; ++i )
  {
    m_calVal.push_back( -1.0 );
    m_calUncert.push_back( -1.0 );
    
    WLabel *label = 0;
    switch( i )
    {
      case 0: label = new WLabel( "Offset" ); break;
      case 1: label = new WLabel( "Linear" ); break;
      case 2: label = new WLabel( "Quadratic" ); break;
      case 3: label = new WLabel( "Quartic" ); break;
      default: label = new WLabel( std::to_string(i+1) + "th" ); break;
    }//switch( i )
    
    auto coefval = new WLineEdit();
    auto fitcb = new WCheckBox( "Fit" );
    
    m_coefvals.push_back( coefval );
    m_fitFor.push_back( fitcb );
    coefval->disable();
    fitcb->setChecked( (i < 2 ) );
    
    fitForLayout->addWidget( label,   i, 0 );
    fitForLayout->addWidget( coefval, i, 1 );
    fitForLayout->addWidget( fitcb,   i, 2 );
  }//for( int i = 0; i < sm_numCoefs; ++i )
  
  assert( m_calVal.size() == m_calUncert.size() );
  assert( m_calVal.size() == m_fitFor.size() );
  assert( m_calVal.size() == m_coefvals.size() );
                      
  fitForLayout->setColumnStretch( 1, 1 );
  
  m_fitSumary = new WTextArea();
  m_fitSumary->setHeight( 75 );
  m_fitSumary->setMaximumSize( WLength::Auto, 75 );
  
  WContainerWidget *instructions = new WContainerWidget();
  WText *line = new WText( "Select peaks to use from each file then click &quot;Fit&quot;.", instructions );
  line->setInline( false );
  line = new WText( "If satisfied, click &quot;Use&quot; to set calibration for involved files.", instructions );
  line->setInline( false );
  line = new WText( "Calibration will be applied to all files with at least one selected peak.", instructions );
  line->setInline( false );

  
  WGridLayout *layout = new WGridLayout( this );
  layout->setContentsMargins( 0, 0, 0, 0 );
  
  layout->addWidget( instructions, 0, 0 );
  layout->addWidget( tree, 1, 0 );
  layout->setRowStretch( 1, 1 );
  layout->addWidget( fitFor, 2, 0 );
  layout->addWidget( m_fitSumary, 3, 0 );
  
  WContainerWidget *buttonDiv = nullptr;
  
  if( parent )
  {
    buttonDiv = parent->footer();
  }else
  {
    buttonDiv = new WContainerWidget();
    layout->addWidget( buttonDiv, 4, 0 );
  }
  
  AuxWindow::addHelpInFooter( buttonDiv, "multi-file-calibration-dialog" );
  
  m_cancel = new WPushButton( "Cancel", buttonDiv );
  m_fit    = new WPushButton( "Fit", buttonDiv );
  m_use    = new WPushButton( "Use", buttonDiv );
  
  m_use->disable();
  m_cancel->clicked().connect( boost::bind( &EnergyCalMultiFile::handleFinish, this, WDialog::Rejected ) );
  m_use->clicked().connect( boost::bind( &EnergyCalMultiFile::handleFinish, this, WDialog::Accepted ) );
  m_fit->clicked().connect( this, &EnergyCalMultiFile::doFit );
  
  m_fitSumary->disable();
  m_fitSumary->hide();
  
  updateCoefDisplay();
  
  if( parent )
  {
    parent->finished().connect( this, &EnergyCalMultiFile::handleFinish );
  
    const int w = 600 < viewer->renderedHeight() ? 600 : viewer->renderedHeight();
    const int h = static_cast<int>(0.8*viewer->renderedHeight());
    parent->resizeWindow( w, h );
  
    parent->rejectWhenEscapePressed();

    parent->centerWindow();
  }//if( parent )
}//EnergyCalMultiFile constructor


EnergyCalMultiFile::~EnergyCalMultiFile()
{
  
}//~EnergyCalMultiFile()


void EnergyCalMultiFile::doFit()
{
  //XXX - The logic of this function is very similar to
  //      Recalibrator::recalibrateByPeaks() - so a refactoriztion should be
  //      done
  auto interspec = InterSpec::instance();
  assert( interspec );
  
  vector< std::shared_ptr<const PeakDef> > peakstouse;
  for( size_t i = 0; i < m_model->m_peaks.size(); ++i )
  {
     const vector< pair<bool,std::shared_ptr<const PeakDef> > > &peaks
                                                          = m_model->m_peaks[i];
     for( size_t j = 0; j < peaks.size(); ++j )
     {
       if( peaks[j].first && peaks[j].second )
         peakstouse.push_back( peaks[j].second );
     }//for( size_t j = 0; j < m_peaks.size(); ++j )
  }//for( size_t i = 0; i < m_peaks.size(); ++i )
  
  try
  {
    const size_t npeaks = peakstouse.size();
    const int num_coeff_fit = m_fitFor[0]->isChecked()
                              + m_fitFor[1]->isChecked()
                              + m_fitFor[2]->isChecked();
    
    if( num_coeff_fit < 1 )
    {
      const char *msg = "You must select at least one coefficient to fit for";
      interspec->logMessage( msg, "", 3 );
      return;
    }//if( num_coeff_fit < 1 )
    
    if( num_coeff_fit > static_cast<int>(npeaks) )
    {
      const char *msg = "You must select at least as many peaks as coeficents to fit for";
      interspec->logMessage( msg, "", 3 );
      return;
    }//if( num_coeff_fit < 1 )
    
    
    std::shared_ptr<const SpecMeas> meas = m_calibrator->m_interspec->measurment(SpectrumType::Foreground);
    
    if( !meas )
    {
      const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit";
      interspec->logMessage( msg, "", 3 );
      return;
    }
    
    
    std::shared_ptr<const Measurement> eqnmeas;
    const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
    
    for( size_t i = 0; !eqnmeas && i < meass.size(); ++i )
      if( meass[i]->num_gamma_channels() )
        eqnmeas = meass[i];
    if( !eqnmeas )
    {
      const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit (unexpected error)";
      interspec->logMessage( msg, "", 3 );
      return;
    }//if( !eqnmeas )
    
    const std::shared_ptr<const std::vector<float>> &binning = eqnmeas->channel_energies();
    const size_t nchannel = eqnmeas->num_gamma_channels();
    
    if( !binning || nchannel < 16 )
    {
      const char *msg = "The spectrum isnt high enough resolution to fit for a calibration";
      interspec->logMessage( msg, "", 3 );
      return;
    }
    
    const SpecUtils::EnergyCalType calibration_type = eqnmeas->energy_calibration_model();
    
    vector<float> calib_coefs = eqnmeas->calibration_coeffs();
    
    switch( calibration_type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
      {
        const char *msg = "Invalid starting calibration type: unknown or lower bin edge energy not"
                          " allowed";
        interspec->logMessage( msg, "", 3 );
        return;
        break;
      }//case LowerChannelEdge or InvalidEquationType
    }//switch( calibration_type )
    
    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contian the error of the fit mean for that peak
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( size_t peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakDef &peak = *peakstouse[peakn];
      
      const double wantedEnergy = peak.gammaParticleEnergy();
      
      EnergyCal::RecalPeakInfo peakInfo;
      peakInfo.peakMean = peak.mean();
      peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
      if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
        peakInfo.peakMeanUncert = 0.5;
        
      peakInfo.photopeakEnergy = wantedEnergy;
      if( m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::FullRangeFraction )
      {
        peakInfo.peakMeanBinNumber = SpecUtils::find_fullrangefraction_channel( peak.mean(),
                                                                  calib_coefs, nchannel,
                                                                  eqnmeas->deviation_pairs(),
                                                                  0.001f );
      }else
      {
        const vector<float> fwfcoef = SpecUtils::polynomial_coef_to_fullrangefraction( calib_coefs, nchannel );
        peakInfo.peakMeanBinNumber = SpecUtils::find_fullrangefraction_channel( peak.mean(),
                                                                  fwfcoef, nchannel,
                                                                  eqnmeas->deviation_pairs(),
                                                                  0.001f );
      }//if( FullRangeFraction ) / else
        
      if( IsNan(peakInfo.peakMeanBinNumber)
          || IsInf(peakInfo.peakMeanBinNumber) )
        throw runtime_error( "Invalid result fromm "
                             "find_fullrangefraction_channel(...)" );
      peakInfos.push_back( peakInfo );
    }//for( int col = 0; col < numModelCol; ++col )
    
    
    bool fit_coefs = false;
    vector<double> parValues, parErrors;
    
    try
    {
      vector<float> lls_fit_coefs = calib_coefs, lls_fit_coefs_uncert;
      
      //Convert eqation type to polynomial incase there are any fixed paramters.
      if( m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::FullRangeFraction )
        lls_fit_coefs = SpecUtils::fullrangefraction_coef_to_polynomial(calib_coefs, nchannel );
      lls_fit_coefs.resize( sm_numCoefs, 0.0f );
      
      vector<bool> fitfor( lls_fit_coefs.size(), false );
      
      for( size_t i = 0; i < calib_coefs.size() && i < sm_numCoefs; ++i )
        fitfor[i] = m_fitFor[i]->isChecked();
      
      const double chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor,
                                              eqnmeas->num_gamma_channels(),
                                              eqnmeas->deviation_pairs(),
                                              lls_fit_coefs, lls_fit_coefs_uncert );
      
      stringstream msg;
      msg << "\nfit_energy_cal_poly gave chi2=" << chi2 << " with coefs={";
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
        msg << lls_fit_coefs[i] << "+-" << lls_fit_coefs_uncert[i] << ", ";
      msg << "}\n";
      cout << msg.str() << endl;
      
      parValues.resize( lls_fit_coefs.size() );
      parErrors.resize( lls_fit_coefs.size() );
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
      {
        parValues[i] = lls_fit_coefs[i];
        parErrors[i] = lls_fit_coefs_uncert[i];
      }
      
      switch( m_calibrator->m_coeffEquationType )
      {
        case SpecUtils::EnergyCalType::Polynomial:
          break;
          
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          m_calibrator->m_coeffEquationType = SpecUtils::EnergyCalType::Polynomial;
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          lls_fit_coefs = SpecUtils::polynomial_coef_to_fullrangefraction(lls_fit_coefs, nchannel);
          lls_fit_coefs_uncert = SpecUtils::polynomial_coef_to_fullrangefraction(lls_fit_coefs_uncert, nchannel);
          break;
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
        case SpecUtils::EnergyCalType::LowerChannelEdge:
          assert( 0 );  //shouldnt ever get here
          break;
      }//switch( m_coeffEquationType )
      
      fit_coefs = true;
    }catch( std::exception &e )
    {
      cerr << "fit_energy_cal_poly threw: " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[512] = { '\0' };
      snprintf( buffer, sizeof(buffer)-1, "fit_energy_cal_poly threw: %s", e.what() );
      log_developer_error( __func__, buffer );
#endif
      fit_coefs = false;
    }//try / catch fit for coefficents using least linear squares
    
    
    if( !fit_coefs )
    {
      //This Minuit based fitting methodolgy is depreciated I think; the LLS code
      //  should work better, and seems to be releiabel, but leaving this code
      //  in for a while as a backup
      const auto eqntype = m_calibrator->m_coeffEquationType;
      const auto &devpairs = eqnmeas->deviation_pairs();
      if( calib_coefs.size() < sm_numCoefs )
        calib_coefs.resize( sm_numCoefs, 0.0 );
      
      vector<bool> fitfor( calib_coefs.size(), false );
      for( size_t i = 0; i < sm_numCoefs; ++i )
        fitfor[i] = m_fitFor[i]->isChecked();
        
      
      string warning_msg;
      vector<float> coefs, coefs_uncert;
      EnergyCal::fit_energy_cal_iterative( peakInfos, nchannel, eqntype, fitfor, calib_coefs,
                                 devpairs, coefs, coefs_uncert, warning_msg );
      
      if( warning_msg.size() )
        passMessage( warning_msg, "", WarningWidget::WarningMsgHigh );
      
      assert( coefs.size() == coefs_uncert.size() );
      parValues.resize( coefs.size(), 0.0 );
      parErrors.resize( coefs_uncert.size(), 0.0 );
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        parValues[i] = coefs[i];
        parErrors[i] = coefs_uncert[i];
      }
    }//if( !fit_coefs )
    
    for( size_t i = 0; i < parValues.size(); ++i )
      if( IsInf(parValues[i]) || IsNan(parValues[i]) )
        throw runtime_error( "Invalid calibration parameter from fit :(" );
    
    assert( parValues.size() >= sm_numCoefs );
    
    m_eqnType = calibration_type;
    for( int i = 0; i < sm_numCoefs; ++i )
    {
      m_calVal[i] = parValues[i];
      m_calUncert[i] = parErrors[i];
    }//for( size_t i = 0; i < parValues.size(); ++i )
    
    //Try to loop over peaks to give chi2 values and such
    vector<float> float_coef;
    for( const double d : parValues )
      float_coef.push_back( static_cast<float>(d) );
    
    if( m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::Polynomial
       || m_calibrator->m_coeffEquationType == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
      float_coef = SpecUtils::polynomial_coef_to_fullrangefraction( float_coef, nchannel );
    
    stringstream msg;
    for( const EnergyCal::RecalPeakInfo &info : peakInfos )
    {
      const double predictedMean = SpecUtils::fullrangefraction_energy( info.peakMeanBinNumber, float_coef, nchannel, eqnmeas->deviation_pairs() );
      double uncert = ((info.peakMeanUncert<=0.0) ? 1.0 : info.peakMeanUncert );
      double chi2 = pow(predictedMean - info.photopeakEnergy, 2.0 ) / (uncert*uncert);
      
      msg << "-Peak originally at " << info.peakMean << " +- "
          << info.peakMeanUncert << " keV for photopeak at "
          << info.photopeakEnergy << " keV ended up at " << predictedMean
          << " keV and contributed " << chi2 << " towards the chi2.\n";
    }//for( const &RecalPeakInfo info : peakInfo )
    
    m_fitSumary->setText( msg.str() );
    m_fitSumary->show();
    updateCoefDisplay();
    m_use->enable();
  }catch( std::exception &e )
  {
    m_use->disable();
    string exceptionmsg = e.what();
    string msg = "Failed calibration by fitting peak means.";
    
    if( SpecUtils::starts_with( exceptionmsg, "ErrorMsg" ) )
      msg = exceptionmsg.substr(8);
    
    cerr << "EnergyCalMultiFile::doFit(): \n\tCaught: " << exceptionmsg << endl;
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void doFit()


void EnergyCalMultiFile::updateCoefDisplay()
{
  assert( m_calVal.size() == m_calUncert.size() );
  assert( m_calVal.size() == m_fitFor.size() );
  assert( m_calVal.size() == m_coefvals.size() );
  
  for( size_t i = 0; i < m_calVal.size(); ++i )
  {
    char msg[32];
    snprintf( msg, sizeof(msg), "%.4g", m_calVal[i] );
    
    WString val(msg);
    if( m_calUncert[i] > 0.0 )
    {
#ifndef WT_NO_STD_WSTRING
      val += L" \x00B1 ";  //plusminus
#else
      val += " +- ";
#endif
      snprintf( msg, sizeof(msg), "%.4g", m_calUncert[i] );
      val += msg;
    }//if( uncertainty is available )
    
    m_coefvals[i]->setText( val );
  }//for( int i = 0; i < sm_numCoefs; ++i )
}//void updateCoefDisplay()


void EnergyCalMultiFile::handleFinish( WDialog::DialogCode result )
{
  switch( result )
  {
    case WDialog::Rejected:
      cerr << "\nRejected EnergyCalMultiFile" << endl;
    break;
      
    case WDialog::Accepted:
    {
      InterSpec *viewer = m_calibrator->m_interspec;
      
      // @TODO: the equation could totally be invalid
      vector<float> eqn( m_calVal, m_calVal + sm_numCoefs );
      while( !eqn.empty() && eqn.back()==0.0 )
        eqn.resize( eqn.size()-1 );
      
      static_assert( sm_numCoefs >= 3, "" );
      cerr << "\n\nm_calVal={" << m_calVal[0] << ", " << m_calVal[1] << ", " << m_calVal[2] << "}" << endl;
      
      vector<float> oldcalibpars;
      std::vector<std::pair<float,float>> devpairs;
      SpecUtils::EnergyCalType oldEqnType;
      
      std::shared_ptr<SpecMeas> fore = viewer->measurment(SpectrumType::Foreground);
      std::shared_ptr<SpecMeas> back = viewer->measurment(SpectrumType::Background);
      std::shared_ptr<SpecMeas> second = viewer->measurment(SpectrumType::SecondForeground);
      std::shared_ptr<const Measurement> displ_foreground
                                     = viewer->displayedHistogram(SpectrumType::Foreground);

      
      std::shared_ptr<const Measurement> eqnmeass[3];
      for( int i = 0; i < 3; ++i )
      {
        const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
        std::shared_ptr<SpecMeas> meas = viewer->measurment(type);
        if( meas )
        {
          const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
          for( size_t j = 0; !eqnmeass[i] && j < meass.size(); ++j )
            if( meass[j]->num_gamma_channels() )
              eqnmeass[i] = meass[j];
        }
      }//for( int i = 0; i < 3; ++i )
      
      if( fore && eqnmeass[ForeIndex] )
      {
        oldcalibpars = eqnmeass[ForeIndex]->calibration_coeffs();
        oldEqnType = eqnmeass[ForeIndex]->energy_calibration_model();
        devpairs = eqnmeass[ForeIndex]->deviation_pairs();
      }else if( second && eqnmeass[SecondIndex] )
      {
        oldcalibpars = eqnmeass[SecondIndex]->calibration_coeffs();
        oldEqnType = eqnmeass[SecondIndex]->energy_calibration_model();
        devpairs = eqnmeass[SecondIndex]->deviation_pairs();
      }else if( back && eqnmeass[BackIndex] )
      {
        oldcalibpars = eqnmeass[BackIndex]->calibration_coeffs();
        oldEqnType = eqnmeass[BackIndex]->energy_calibration_model();
        devpairs = eqnmeass[BackIndex]->deviation_pairs();
      }else
      {
        break;
      }
      
      vector<string> displayed_detectors = viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
      
      
      if( !m_calibrator->m_applyTo[ForeIndex]->isChecked() )
        fore.reset();
      if( !m_calibrator->m_applyTo[BackIndex]->isChecked() )
        back.reset();
      if( !m_calibrator->m_applyTo[SecondIndex]->isChecked() )
        second.reset();
      
      //XXX - applying the calibration parameters does not shift peaks!
      for( size_t i = 0; i < m_model->m_peaks.size(); ++i )
      {
        const vector< pair<bool,std::shared_ptr<const PeakDef> > > &peaks
                                                          = m_model->m_peaks[i];
        bool used = false;
        for( size_t j = 0; j < peaks.size(); ++j )
          used |= peaks[j].first;
        
        if( !used )
          continue;
        
        SpectraFileModel *fileModel = viewer->fileManager()->model();
        std::shared_ptr<SpectraFileHeader> header
                                = fileModel->fileHeader( static_cast<int>(i) );
        if( !header )
          continue;
        
        std::shared_ptr<SpecMeas> meas = header->parseFile();
        if( !meas )
          continue;
        
        // Fix this for new calibration! //meas->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
        cerr << "\n\nRecalled " << header->displayName().toUTF8()
             << " to m_calVal={" << m_calVal[0] << ", " << m_calVal[1]
             << ", " << m_calVal[2] << "}" << endl;
        
        if( meas == fore )
        {
          fore.reset();
          Recalibrator::shiftPeaksForEnergyCalibration(
                                            m_calibrator->m_peakModel,
                                            eqn, devpairs, m_eqnType,
                                            meas, SpectrumType::Foreground,
                                            oldcalibpars, devpairs, oldEqnType );
        }else
        {
          assert( 0 );
          //meas->shiftPeaksForRecalibration( oldcalibpars, devpairs, oldEqnType,
          //                                 eqn, devpairs, m_eqnType );
        }//if( meas == fore ) / else
        
        if( meas == back )
          back.reset();
        if( meas == second )
          second.reset();
      }//for( size_t i = 0; i < m_peaks.size(); ++i )
      
      // Fix this for new calibration! //if( fore )
      // Fix this for new calibration! //  fore->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      // Fix this for new calibration! //if( back )
      // Fix this for new calibration! //  back->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      // Fix this for new calibration! //if( second )
      // Fix this for new calibration! //  second->recalibrate_by_eqn( eqn, devpairs, m_eqnType, displayed_detectors, false );
      
      viewer->refreshDisplayedCharts();
      m_calibrator->refreshRecalibrator();
      
      cerr << "\nAccepted EnergyCalMultiFile" << endl;
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  delete this;
}//void handleFinish(...)



EnergyCalMultiFileModel::EnergyCalMultiFileModel( EnergyCalTool *calibrator, Wt::WObject *parent )
: WAbstractItemModel( parent ),
  m_calibrator( calibrator ),
  m_fileModel( nullptr )
{
  refreshData();
  
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  
  SpecMeasManager *measManager = interspec->fileManager();
  assert( SpecMeasManager );
  
  m_fileModel = measManager->model();
  assert( m_fileModel );
  
  m_fileModel->rowsInserted().connect( boost::bind(&EnergyCalMultiFileModel::refreshData, this) );
  m_fileModel->rowsRemoved().connect( boost::bind(&EnergyCalMultiFileModel::refreshData, this) );
  
  //Should add in a listener here to the PeakModel to see if peaks are
  //  added/removed; this might neccesitate tracking
}//EnergyCalMultiFileModel constructor


EnergyCalMultiFileModel::~EnergyCalMultiFileModel()
{
}//~EnergyCalMultiFileModel()


WModelIndex EnergyCalMultiFileModel::index( int row, int column, const WModelIndex &parent ) const
{
  if( parent.isValid() )
  {
    if( column < 0 || column >= 3 )
      return WModelIndex();
    
    if( parent.row() < 0 || parent.row() >= static_cast<int>(m_peaks.size()) )
      return WModelIndex();
    
    return createIndex( row, column, ::uint64_t(parent.row()) );
  }//if( parent.isValid() )
  
  if( row < 0 || row >= static_cast<int>(m_peaks.size()) )
    return WModelIndex();
  
  if( column != 0 )
    return WModelIndex();
  return createIndex( row, column, numeric_limits<unsigned int>::max() );
}//index(...)


WModelIndex EnergyCalMultiFileModel::parent( const WModelIndex &index ) const
{
  if( !index.isValid() )
    return WModelIndex();
  
  const ::uint64_t filenum = index.internalId();
  if( filenum == numeric_limits<unsigned int>::max() )
    return WModelIndex();
  
  if( filenum < m_peaks.size() )
    return createIndex( filenum, 0, numeric_limits<unsigned int>::max() );
  
  return WModelIndex();
}//parent(...)


int EnergyCalMultiFileModel::rowCount( const WModelIndex &parent ) const
{
  if( !parent.isValid() )
    return static_cast<int>( m_peaks.size() );
  
  if( parent.parent().isValid() )
    return 0;
  
  const int filenum = parent.row();
  if( filenum >= 0 && filenum < static_cast<int>(m_peaks.size()) )
    return static_cast<int>( m_peaks[filenum].size() );
  
  return 0;
}//rowCount(...)


int EnergyCalMultiFileModel::columnCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 3;
  return 1;
}//columnCount(...)

  
boost::any EnergyCalMultiFileModel::data( const Wt::WModelIndex &index, int role ) const
{
  if( (role != Wt::DisplayRole) && (role != Wt::CheckStateRole) )
    return boost::any();
  
  if( !index.isValid() )
    return boost::any();
  
  if( !index.parent().isValid() && index.column() != 0 )
    return boost::any();
  
  const int row = index.row();
  const int col = index.column();
  
  if( index.parent().isValid() )
  {
    const int filenum = index.parent().row();
    if( filenum < 0 || filenum >= static_cast<int>(m_peaks.size()) )
      return boost::any();
    if( row < 0 || row >= static_cast<int>(m_peaks[filenum].size()) )
      return boost::any();
    
    const bool checked = m_peaks[filenum][row].first;
    std::shared_ptr<const PeakDef> peak = m_peaks[filenum][index.row()].second;
    
    if( !peak )
      return boost::any();
    
    if( !peak->parentNuclide() && !peak->xrayElement() && !peak->reaction() )
      return boost::any();
    
    if( role == Wt::CheckStateRole )
    {
      if( col == 0 )
        return boost::any( checked );
      return boost::any();
    }
    
    if( col == 0 )
    {
      if( peak->parentNuclide() )
        return boost::any( WString::fromUTF8(peak->parentNuclide()->symbol) );
      if( peak->xrayElement() )
        return boost::any( WString::fromUTF8(peak->xrayElement()->symbol) );
      if( peak->reaction() )
        return boost::any( WString::fromUTF8(peak->reaction()->name()) );
      return boost::any( WString("--") );
    }//if( col == 0 )
    
    char msg[128];
    msg[0] = '\0';
    
    if( col == 1 )
      snprintf(msg, sizeof(msg), "%.2f keV Detected", peak->mean() );
    
    if( col == 2 )
    {
      try
      {
        const float wantedEnergy = peak->gammaParticleEnergy();
        snprintf(msg, sizeof(msg), "%.2f keV Expected", wantedEnergy );
      }catch(...)
      {
        snprintf(msg, sizeof(msg), "%s", "--" );
      }
    }//if( col == 2 )
    
    if( !strlen(msg) )
      return boost::any();
    return boost::any( WString(msg) );
  }//if( index.parent().isValid() )
  
  if( role != Wt::DisplayRole )
    return boost::any();
  
  if( row < 0 || row >= static_cast<int>(m_peaks.size()) )
    return boost::any();
  
  std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( row );
  if( header )
    return boost::any( header->displayName() );
  return boost::any();
}//data(...)


bool EnergyCalMultiFileModel::setData( const WModelIndex &index, const boost::any &value, int role )
{
  if( role != Wt::CheckStateRole
      || !index.isValid() || !index.parent().isValid() )
    return false;
  
  const int filenum = index.parent().row();
  const int row = index.row();
  const int col = index.column();
  
  if( col != 0 )
    return false;
  if( filenum >= static_cast<int>(m_peaks.size()) )
    return false;
  if( row >= static_cast<int>(m_peaks[filenum].size()) )
    return false;
  
  try
  {
    m_peaks[filenum][row].first = boost::any_cast<bool>( value );
  }catch(...)
  {
    cerr << "Got bad boost::any_cast<bool>( value )" << endl;
    return false;
  }
  
  return true;
}//setData(...)



WFlags<ItemFlag> EnergyCalMultiFileModel::flags( const WModelIndex &index ) const
{
  if( index.isValid() && index.parent().isValid() && index.column()==0 )
    return WFlags<ItemFlag>(ItemIsUserCheckable);
  return WFlags<ItemFlag>(ItemIsXHTMLText);
}//flags(...)


boost::any EnergyCalMultiFileModel::headerData( int section, Wt::Orientation orientation,
                                                int role ) const
{
  if( orientation != Wt::Horizontal )
    return boost::any();
  if( role != DisplayRole )
    return boost::any();
  
  if( section == 0 )
    boost::any( WString("Nuclide") );
  if( section == 1 )
    boost::any( WString("Obs. Energy") );
  if( section == 2 )
    boost::any( WString("Photopeak Energy") );
  return boost::any();
}//headerData(...);


void EnergyCalMultiFileModel::refreshData()
{
  if(!m_peaks.empty())
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_peaks.size())-1 );
    m_peaks.clear();
    endRemoveRows();
  }//if( m_peaks.size() )
  
  vector<std::vector< std::pair<bool,std::shared_ptr<const PeakDef> > > > newdata;
  
  const int nfile = m_fileModel->rowCount();
  for( int filenum = 0; filenum < nfile; ++filenum )
  {
    vector< pair<bool,std::shared_ptr<const PeakDef> > > peaks;
    std::shared_ptr<SpectraFileHeader> header
                                           = m_fileModel->fileHeader( filenum );
    if( !header )
    {
      newdata.push_back( peaks );
      continue;
    }//if( !header )
    
    std::shared_ptr<SpecMeas> spec = header->parseFile();
    
    if( !spec )
    {
      newdata.push_back( peaks );
      continue;
    }//if( !header )
    
    typedef set<int> IntSet;
    typedef std::shared_ptr<const PeakDef> PeakPtr;
    const set< set<int> > peaksamplenums = spec->sampleNumsWithPeaks();
    for( const IntSet &samplnums : peaksamplenums )
    {
      std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > measpeaks;
      measpeaks = spec->peaks( samplnums );
      if( !measpeaks )
        continue;
      
      for( const PeakPtr &peak : *measpeaks )
      {
        const bool use = peak->useForCalibration();
        pair<bool,std::shared_ptr<const PeakDef> > peakpair( use, peak );
        peaks.push_back( peakpair );
      }//for( const PeakPtr &peak : peaks )
      
      newdata.push_back( peaks );
    }//for( const IntSet &samplnums : peaksamplenums )
  }//for( int filenum = 0; filenum < nfile; ++filenum )

  if(!newdata.empty())
  {
    beginInsertRows( WModelIndex(), 0, static_cast<int>(newdata.size()) -1 );
    m_peaks.swap( newdata );
    endInsertRows();
  }//if( newdata.size() )
}//refreshData()


*/
