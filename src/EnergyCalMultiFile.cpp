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
#include <map>
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
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
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
using SpecUtils::EnergyCalibration;

namespace
{
  const size_t ns_min_num_coef = 4;
}

EnergyCalMultiFile::EnergyCalMultiFile( EnergyCalTool *cal, AuxWindow *parent )
: WContainerWidget(),
  m_calibrator( cal ),
  m_parent( parent ),
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
    WLabel *label = 0;
    switch( i )
    {
      case 0: label = new WLabel( "Offset" ); break;
      case 1: label = new WLabel( "Linear" ); break;
      case 2: label = new WLabel( "Quadratic" ); break;
      case 3: label = new WLabel( "Cubic" ); break;
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
  
    const int w = 600 < viewer->renderedWidth() ? 600 : viewer->renderedWidth();
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
  auto interspec = InterSpec::instance();
  assert( interspec );
  
  m_calVal.clear();
  m_calUncert.clear();
  m_devPairs.clear();
  
  shared_ptr<const SpecMeas> meas = interspec->measurment(SpectrumType::Foreground);
  shared_ptr<const Measurement> dispmeas = interspec->displayedHistogram(SpectrumType::Foreground);
  shared_ptr<const EnergyCalibration> disp_cal = dispmeas ? dispmeas->energy_calibration() : nullptr;
  
  if( !disp_cal || !disp_cal->valid() || disp_cal->num_channels() < 16 )
  {
    const char *msg = "You need to be displaying a foreground spectrum to do a calibration fit";
    interspec->logMessage( msg, "", 3 );
    return;
  }
  
  
  try
  {
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( size_t filenum = 0; filenum < m_model->m_peaks.size(); ++filenum )
    {
      const vector<EnergyCalMultiFileModel::PeakInfo_t> &peaks = m_model->m_peaks[filenum];
      for( size_t j = 0; j < peaks.size(); ++j )
      {
        if( get<0>(peaks[j]) && get<1>(peaks[j]) && get<2>(peaks[j]) )
        {
          const PeakDef &peak = *get<2>(peaks[j]);
          const shared_ptr<const EnergyCalibration> &peakcal = get<0>(peaks[j]);
          const double wantedEnergy = peak.gammaParticleEnergy();
          
          EnergyCal::RecalPeakInfo peakInfo;
          peakInfo.peakMean = peak.mean();
          peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
          if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
            peakInfo.peakMeanUncert = 0.5;
          peakInfo.photopeakEnergy = wantedEnergy;
          peakInfo.peakMeanBinNumber = peakcal->channel_for_energy( peak.mean() );
          
          if( IsNan(peakInfo.peakMeanBinNumber) || IsInf(peakInfo.peakMeanBinNumber) )
            throw runtime_error( "Invalid result from EnergyCalibration::channel_for_energy(...)" );
          peakInfos.push_back( peakInfo );
        }//if( energy cal and peak ptrs are valid, and we should use this peak for fitting )
      }//for( loop over peaks for a file )
    }//for( int col = 0; col < numModelCol; ++col )
    
    
    const size_t npeaks = peakInfos.size();
    const size_t ncoeffs = m_fitFor.size();
    const size_t nchannel = disp_cal->num_channels();
    
    int num_coeff_fit = 0;
    vector<bool> fitfor( ncoeffs, false );
    for( size_t i = 0; i < m_fitFor.size(); ++i )
    {
      fitfor[i] = m_fitFor[i]->isChecked();
      num_coeff_fit += m_fitFor[i]->isChecked();
    }
    
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
    
    
    bool fit_coefs = false;
    try
    {
      vector<float> lls_fit_coefs( ncoeffs, 0.0f ), lls_fit_coefs_uncert( ncoeffs, 0.0f );
      
      const double chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor,
                                               disp_cal->num_channels(),
                                               disp_cal->deviation_pairs(),
                                              lls_fit_coefs, lls_fit_coefs_uncert );
      
      stringstream msg;
      msg << "\nfit_energy_cal_poly gave chi2=" << chi2 << " with coefs={";
      for( size_t i = 0; i < lls_fit_coefs.size(); ++i )
        msg << lls_fit_coefs[i] << "+-" << lls_fit_coefs_uncert[i] << ", ";
      msg << "}\n";
      cout << msg.str() << endl;
      
      fit_coefs = true;
      m_calVal = lls_fit_coefs;
      m_calUncert = lls_fit_coefs_uncert;
      m_devPairs = disp_cal->deviation_pairs();
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
      const auto &devpairs = disp_cal->deviation_pairs();
      vector<float> starting_coefs( ncoeffs, 0.0 );
      starting_coefs[1] = disp_cal->upper_energy() / disp_cal->num_channels();
      
      string warning_msg;
      vector<float> coefs, coefs_uncert;
      EnergyCal::fit_energy_cal_iterative( peakInfos, nchannel,
                              SpecUtils::EnergyCalType::Polynomial, fitfor, starting_coefs,
                              devpairs, coefs, coefs_uncert, warning_msg );
      
      if( warning_msg.size() )
        interspec->logMessage( warning_msg, "", 3 );
      
      assert( coefs.size() == ncoeffs );
      assert( coefs.size() == coefs_uncert.size() );
      
      for( size_t i = 0; i < coefs.size(); ++i )
        if( IsInf(coefs[i]) || IsNan(coefs[i]) )
          throw runtime_error( "Invalid calibration parameter from fit :(" );
      
      fit_coefs = true;
      m_calVal = coefs;
      m_calUncert = coefs_uncert;
      m_devPairs = disp_cal->deviation_pairs();
    }//if( !fit_coefs )
    
    //Try to loop over peaks to give chi2 values and such
    //  \TODO: Put this information in the table (e.g., modify EnergyCalMultiFileModel to hold it)
    stringstream msg;
    for( const EnergyCal::RecalPeakInfo &info : peakInfos )
    {
      const double predictedMean = SpecUtils::polynomial_energy( info.peakMeanBinNumber, m_calVal, m_devPairs );
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
    interspec->logMessage(msg, "", 3);
  }//try / catch
}//void doFit()


void EnergyCalMultiFile::updateCoefDisplay()
{
  if( m_calVal.empty() )
  {
    for( size_t i = 0; i < m_coefvals.size(); ++i )
      m_coefvals[i]->setText( "--" );
    return;
  }
  
  assert( m_calVal.size() == m_calUncert.size() );
  assert( m_calVal.size() == m_fitFor.size() );
  assert( m_calVal.size() == m_coefvals.size() );
  
  for( size_t i = 0; i < m_calVal.size(); ++i )
  {
    char msg[64] = { '\0' };
    if( m_calUncert[i] > 0.0 )
      snprintf( msg, sizeof(msg), "%.4g \xC2\xB1 %.4g", m_calVal[i], m_calUncert[i] );
    else
      snprintf( msg, sizeof(msg), "%.4g", m_calVal[i] );
    
    m_coefvals[i]->setText( WString::fromUTF8(msg) );
  }//for( int i = 0; i < sm_numCoefs; ++i )
}//void updateCoefDisplay()


void EnergyCalMultiFile::applyCurrentFit()
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  if( m_calVal.size() < 2 )
  {
    viewer->logMessage( "Currently fit calibration is invalid", "", 3 );
    return;
  }
  
  vector<pair<string,string>> error_msgs;
  const auto foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const auto dispsamples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  
  for( size_t filenum = 0; filenum < m_model->m_peaks.size(); ++filenum )
  {
    shared_ptr<SpectraFileHeader> header;
    const vector<EnergyCalMultiFileModel::PeakInfo_t> &peaks = m_model->m_peaks[filenum];
    for( size_t j = 0; !header && (j < peaks.size()); ++j )
    {
      if( get<0>(peaks[j]) && get<1>(peaks[j]) && get<2>(peaks[j]) && get<3>(peaks[j]) )
        header = get<3>(peaks[j]);
    }
    
    if( !header )
      continue;
    
    // We will apply the calibrations on a file-by-file level; this has the side-effect that if we
    //  do run into an issue, that one file wont pick up the change, but all other files will...
    //  not sure if this is the right way or not.
    try
    {
      shared_ptr<SpecMeas> spec = header->parseFile();
      if( !spec ) //shouldnt happen
        continue;
      
      set<shared_ptr<const EnergyCalibration>> newcals;
      map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
      map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
      
      // We will first update the energy calibration for Measurements associated with peaks, since
      //  this is kinda un-ambigous for how to do it for multi-deector systems.
      const auto &detnames = spec->gamma_detector_names();
      const auto peaksets = spec->sampleNumsWithPeaks();
      for( const set<int> &samples : peaksets )
      {
        auto peaks = spec->peaks( samples );
        if( !peaks || peaks->empty() )
          continue;
        
        auto dispcal = spec->suggested_sum_energy_calibration( samples, detnames );
        if( !dispcal || !dispcal->valid() )
          continue;
        
        shared_ptr<const EnergyCalibration> newcal;
        auto calpos = updated_cals.find( dispcal );
        if( calpos != end(updated_cals) )
        {
          newcal = calpos->second;
        }else
        {
          auto cal = make_shared<EnergyCalibration>();
          cal->set_polynomial( dispcal->num_channels(), m_calVal, m_devPairs );
          newcal = cal;
          calpos = updated_cals.insert( {dispcal, newcal} ).first;
        }
        
        
        for( const int sample : samples )
        {
          for( const auto m : spec->sample_measurements(sample) )
          {
            auto cal = m ? m->energy_calibration() : nullptr;
            if( !cal || !cal->valid() || (cal->num_channels() < 5) )
              continue;
            
            if( cal == dispcal )
            {
              updated_cals[cal] = newcal;
              newcals.insert( newcal );
            }else
            {
              auto pos = updated_cals.find( dispcal );
              if( pos == end(updated_cals) )
              {
                auto thisnewcal = EnergyCal::propogate_energy_cal_change( dispcal, newcal, cal );
                updated_cals[cal] = thisnewcal;
                newcals.insert( thisnewcal );
              }
            }//
          }//for( const auto m : spec->sample_measurements(sample) )
        }//for( const int sample : samples )
        
        updated_peaks[peaks] = EnergyCal::translatePeaksForCalibrationChange( *peaks, dispcal, newcal );
      }//for( const set<int> &samples : peaksets )
      
      
      for( shared_ptr<const SpecUtils::Measurement> m : spec->measurements() )
      {
        shared_ptr<const EnergyCalibration> oldcal = m ? m->energy_calibration() : nullptr;
        if( !oldcal || !oldcal->valid() || (oldcal->num_channels() < 5) )
          continue;
        
        auto pos = updated_cals.find( oldcal );
        if( pos != end(updated_cals) )
          continue;

        // Uhg, for multi-detector systems, I dont really know what to do, so heres something...
        //  \TODO: figure out a better way to apply the correct calibration; maybe match up detect
        //         names or somethign?
        auto dispcal = spec->suggested_sum_energy_calibration( {m->sample_number()}, detnames );
        if( !dispcal )
          dispcal = oldcal;
        
        auto newdispcal = make_shared<EnergyCalibration>();
        newdispcal->set_polynomial( oldcal->num_channels(), m_calVal, m_devPairs);
        
        auto newcal = EnergyCal::propogate_energy_cal_change( dispcal, newdispcal, oldcal );
        updated_cals[oldcal] = newcal;
        newcals.insert( newcal );
      }//for( auto m : spec->measurements() )
      
      //We have made all the new energy calibrations, and translated all the peaks; lets set them
      for( shared_ptr<const SpecUtils::Measurement> m : spec->measurements() )
      {
        shared_ptr<const EnergyCalibration> oldcal = m ? m->energy_calibration() : nullptr;
        auto pos = updated_cals.find( oldcal );
        if( pos != end(updated_cals) )
          spec->set_energy_calibration( pos->second, m );
      }//for( auto m : spec->measurements() )
      
      
      for( const set<int> &samples : peaksets )
      {
        auto peaks = spec->peaks( samples );
        if( !peaks || peaks->empty() )
          continue;
        auto pos = updated_peaks.find( peaks );
        if( pos != end(updated_peaks) )
          spec->setPeaks( pos->second, samples );
      }//for( const set<int> &samples : peaksets )
      
      assert( foreground );
      //Trigger an update of peak views if we are on the foreground SpecFile
      if( (spec == foreground) && foreground->peaks(dispsamples) )
        viewer->peakModel()->setPeakFromSpecMeas(foreground,dispsamples);
    }catch( std::exception &e )
    {
      const string filename = header ? header->displayName().toUTF8() : std::string("unknown");
      error_msgs.emplace_back( SpecUtils::filename(filename), e.what() );
    }//try / catch to set the calibrations
  }//for( loop over files )
  
  viewer->refreshDisplayedCharts();
  m_calibrator->refreshGuiFromFiles();
  
  string errmsg;
  for( size_t i = 0; (i < error_msgs.size()) && (i < 5); ++i )
  {
    const string &filename = error_msgs[i].first;
    const string &msg = error_msgs[i].second;
    errmsg += "<p>" + error_msgs[i].first + ":"
              " <span style=\"font-style: italic; font-family: monospace; color: red;\">"
              + error_msgs[i].second
              + "</span></p>";
  }//for( loop over err messages )
  
  if( !errmsg.empty() )
  {
    errmsg = "There was an error an error changing calibration of "
             + std::to_string(error_msgs.size()) + " file:\n" + errmsg;
    viewer->logMessage( WString::fromUTF8(errmsg), "", 3 );
  }else
  {
    viewer->logMessage( "Have updated calibration for all files with any peaks selected.", "", 1 );
  }
}//void applyCurrentFit()


void EnergyCalMultiFile::handleFinish( WDialog::DialogCode result )
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  switch( result )
  {
    case WDialog::Rejected:
      cerr << "\nRejected EnergyCalMultiFile" << endl;
    break;
      
    case WDialog::Accepted:
    {
      applyCurrentFit(); 
      cerr << "\nAccepted EnergyCalMultiFile" << endl;
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  if( m_parent )
    AuxWindow::deleteAuxWindow( m_parent );
}//void handleFinish(...)



EnergyCalMultiFileModel::EnergyCalMultiFileModel( EnergyCalTool *calibrator, Wt::WObject *parent )
: WAbstractItemModel( parent ),
  m_calibrator( calibrator ),
  m_fileModel( nullptr )
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  
  SpecMeasManager *measManager = interspec->fileManager();
  assert( measManager );
  
  m_fileModel = measManager->model();
  assert( m_fileModel );
  
  refreshData();
  
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
    
    const PeakInfo_t &info = m_peaks[filenum][row];
    
    const bool checked = get<1>(info);
    shared_ptr<const PeakDef> peak = get<2>(info);
    
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
    get<1>(m_peaks[filenum][row]) = boost::any_cast<bool>( value );
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
  
  vector<vector<PeakInfo_t>> newdata;
  
  const int nfile = m_fileModel->rowCount();
  for( int filenum = 0; filenum < nfile; ++filenum )
  {
    vector<PeakInfo_t> peaks;
    shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( filenum );
    
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
    
    const auto gammaDetNames = spec->gamma_detector_names();
    
    typedef std::shared_ptr<const PeakDef> PeakPtr;
    const set<set<int>> peaksamplenums = spec->sampleNumsWithPeaks();
    for( const set<int> &samplnums : peaksamplenums )
    {
      shared_ptr<const EnergyCalibration> cal;
      try
      {
        cal = spec->suggested_sum_energy_calibration( samplnums, gammaDetNames );
      }catch( std::exception & )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Unexpected failure of suggested_sum_energy_calibration" );
#endif
        continue;
      }//try / catch
      
      auto measpeaks = spec->peaks( samplnums );
      if( !measpeaks )
        continue;
      
      for( const PeakPtr &peak : *measpeaks )
        peaks.emplace_back( cal, peak->useForCalibration(), peak, header );
      
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
