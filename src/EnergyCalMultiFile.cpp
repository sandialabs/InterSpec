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

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/SpectraFileModel.h"
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


/** For the EnergyCalMultiFileModel class, we will create a unique WModelIndex internal ID ptr by multiplying each level if information
 by some multiples, and using the bottom three bits to indicate what type of index this is.
 
 Maybe not super-solid, but whatever.
 */
namespace
{
const uint64_t ns_file_level_multiple   = 0x0001000000000000;
const uint64_t ns_sample_level_multiple = 0x0000000010000000;
const uint64_t ns_peak_level_multiple   = 0x0000000000000010;

const uint64_t ns_is_file_bit           = 0x0000000000000004;
const uint64_t ns_is_sample_bit         = 0x0000000000000002;
const uint64_t ns_is_peak_bit           = 0x0000000000000001;

enum class ModelLevel
{
  File, SampleSet, Peak, Invalid
};

uint64_t to_internal_id( ModelLevel level, int filenum, int samplesnum, int peaknum )
{
  uint64_t file_ptr_part = filenum;
  file_ptr_part *= ns_file_level_multiple;
  
  uint64_t sample_ptr_part = samplesnum;
  sample_ptr_part *= ns_sample_level_multiple;
  
  uint64_t peak_ptr_part = peaknum;
  peak_ptr_part *= ns_peak_level_multiple;
  
  uint64_t internalptr = file_ptr_part + sample_ptr_part + peak_ptr_part;
  
  switch( level )
  {
    case ModelLevel::File:       internalptr |= ns_is_file_bit;   break;
    case ModelLevel::SampleSet:  internalptr |= ns_is_sample_bit; break;
    case ModelLevel::Peak:       internalptr |= ns_is_peak_bit;   break;
    case ModelLevel::Invalid:                                     break;
  }//switch( level )
  
  return internalptr;
}//to_internal_id(...)


void from_internal_id( const uint64_t internal_id, ModelLevel &level, int &filenum, int &samplesnum, int &peaknum )
{
  filenum = samplesnum = peaknum = 0;
  
  const bool is_filelevel   = (internal_id & ns_is_file_bit);
  const bool is_samplelevel = (internal_id & ns_is_sample_bit);
  const bool is_peaklevel   = (internal_id & ns_is_peak_bit);
  
  const int nlevels = is_filelevel + is_samplelevel + is_peaklevel;
  assert( nlevels <= 1 );
  
  if( nlevels > 1 )
  {
    level = ModelLevel::Invalid;
    return;
  }
  
  if( is_filelevel )
    level = ModelLevel::File;
  else if( is_samplelevel )
    level = ModelLevel::SampleSet;
  else if( is_peaklevel )
    level = ModelLevel::Peak;
  else
    level = ModelLevel::Invalid;
  
  uint64_t file_part   = internal_id / ns_file_level_multiple;
  uint64_t sample_part = ((internal_id % ns_file_level_multiple)   / ns_sample_level_multiple);
  uint64_t peak_part   = ((internal_id % ns_sample_level_multiple) / ns_peak_level_multiple);
  
  switch( level )
  {
    case ModelLevel::Invalid:
      assert( file_part == 0 );
      // fall-through intentional
      
    case ModelLevel::File:
      assert( sample_part == 0 );
      // fall-through intentional
      
    case ModelLevel::SampleSet:
      assert( peak_part == 0 );
      break;
      
    case ModelLevel::Peak:
      break;
  }//switch( level )
  
  filenum = static_cast<int>(file_part);
  samplesnum = static_cast<int>(sample_part);
  peaknum = static_cast<int>(peak_part);
}//void from_internal_id(...)


}//namespace


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
  tree->setHeaderHeight( 20 );
  tree->setColumnWidth( 0, 200 );
  tree->setColumnWidth( 1, 150 );
  tree->setColumnWidth( 2, 150 );
  tree->expandToDepth( 2 );
  
  WGroupBox *fitFor = new WGroupBox( "Coefficients to fit for" );
  WGridLayout *fitForLayout = new WGridLayout( fitFor );
  
  
  for( int i = 0; i < static_cast<int>(ns_min_num_coef); ++i )
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
    
    coefval->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    coefval->setAttributeValue( "autocorrect", "off" );
    coefval->setAttributeValue( "spellcheck", "off" );
#endif
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
    interspec->logMessage( msg, 3 );
    return;
  }
  
  
  try
  {
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( const auto &samplesinfos : m_model->m_data )
    {
      for( const EnergyCalMultiFileModel::SamplesPeakInfo_t &samplesinfo : samplesinfos )
      {
        const shared_ptr<const SpecUtils::EnergyCalibration> &peakcal = get<2>(samplesinfo);
        const vector<EnergyCalMultiFileModel::UsePeakInfo_t> &peakinfos = get<3>(samplesinfo);
        
        if( !peakcal || !peakcal->valid() )
          continue;
        
        for( const EnergyCalMultiFileModel::UsePeakInfo_t &info : peakinfos )
        {
          const bool use = get<0>(info);
          const shared_ptr<const PeakDef> &peakptr = get<1>(info);
          
          if( use && peakptr )
          {
            const PeakDef &peak = *peakptr;
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
    }//for( const EnergyCalMultiFileModel::SamplesPeakInfo_t &samplesinfo : samplesinfos )
    
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
      interspec->logMessage( msg, 3 );
      return;
    }//if( num_coeff_fit < 1 )
    
    if( num_coeff_fit > static_cast<int>(npeaks) )
    {
      const char *msg = "You must select at least as many peaks as coefficients to fit for";
      interspec->logMessage( msg, 3 );
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
        interspec->logMessage( warning_msg, 3 );
      
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
    interspec->logMessage(msg, 3);
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
    viewer->logMessage( "Currently fit calibration is invalid", 3 );
    return;
  }
  
  vector<pair<string,string>> error_msgs;
  const auto foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const auto dispsamples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  
  for( size_t filenum = 0; filenum < m_model->m_data.size(); ++filenum )
  {
    const vector<EnergyCalMultiFileModel::SamplesPeakInfo_t> &samplesinfos = m_model->m_data[filenum];
    
    if( samplesinfos.empty() )
      continue;
    
    const shared_ptr<SpectraFileHeader> header = get<0>(samplesinfos[0]);
    if( !header )
      continue;
    
    shared_ptr<SpecMeas> spec;
    try
    {
      spec = header->parseFile();
    }catch( std::exception &e )
    {
      const string filename = header ? header->displayName().toUTF8() : std::string("unknown");
      error_msgs.emplace_back( SpecUtils::filename(filename), e.what() );
    }
    
    if( !spec ) //shouldnt happen, but JIC
      continue;
    
    //  TODO: decide if we should apply changes only for samples that have a peak participating in the fit...
    int num_samples_used = 0;
    
    for( const EnergyCalMultiFileModel::SamplesPeakInfo_t &samplesinfo : samplesinfos )
    {
      for( const EnergyCalMultiFileModel::UsePeakInfo_t &peakInfo : get<3>(samplesinfo) )
      {
        const bool use = get<0>(peakInfo);
        const shared_ptr<const PeakDef> &peak = get<1>(peakInfo);
        
        if( use && peak && peak->hasSourceGammaAssigned() )
          num_samples_used += 1;
      }
    }//for( const EnergyCalMultiFileModel::SamplesPeakInfo_t &samplesinfo : samplesinfos )
    
    if( !num_samples_used )
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
      map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks; //old peaks, to new peaks
      map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_auto_search_peaks;

      // We will first update the energy calibration for Measurements associated with peaks, since
      //  this is kinda un-ambiguous for how to do it for multi-detector systems.
      const auto &detnames = spec->gamma_detector_names();
      const auto peaksets = spec->sampleNumsWithPeaks();
      const set<set<int>> samplesWithAutoSearchPeak = spec->sampleNumsWithAutomatedSearchPeaks();

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


      for( const set<int> &samples : samplesWithAutoSearchPeak )
      {
        shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = spec->automatedSearchPeaks(samples);
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

        updated_auto_search_peaks[peaks] = EnergyCal::translatePeaksForCalibrationChange( *peaks, dispcal, newcal );
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


      for( const set<int> &samples : samplesWithAutoSearchPeak )
      {
        shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = spec->automatedSearchPeaks( samples );
        if( !peaks || peaks->empty() )
          continue;
        const auto pos = updated_auto_search_peaks.find( peaks );
        if( pos != end(updated_auto_search_peaks) )
        {
          auto peaks = make_shared<std::deque<std::shared_ptr<const PeakDef>>>(pos->second);
          spec->setAutomatedSearchPeaks( samples, peaks );
        }
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
    viewer->logMessage( WString::fromUTF8(errmsg), 3 );
  }else
  {
    viewer->logMessage( "Have updated calibration for all files with any peaks selected.", 1 );
  }
}//void applyCurrentFit()


void EnergyCalMultiFile::handleFinish( WDialog::DialogCode result )
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  switch( result )
  {
    case WDialog::Rejected:
    {
      cerr << "\nRejected EnergyCalMultiFile" << endl;
      
      UndoRedoManager *undoManager = viewer->undoRedoManager();
      if( m_parent && undoManager )
      {
        auto undo = [](){
          InterSpec *viewer = InterSpec::instance();
          EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
          if( tool )
            tool->moreActionBtnClicked( MoreActionsIndex::MultipleFilesCal );
        };
        
        auto redo = [](){
          InterSpec *viewer = InterSpec::instance();
          EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
          if( tool )
            tool->cancelMoreActionWindow();
        };
        
        undoManager->addUndoRedoStep( undo, redo, "Cancel multi-file energy cal" );
      }//if( m_parent && undoManager )
      
      break;
    }//case WDialog::Rejected:
      
    case WDialog::Accepted:
    {
      applyCurrentFit(); 
      cerr << "\nAccepted EnergyCalMultiFile" << endl;
      
      UndoRedoManager *undoManager = viewer->undoRedoManager();
      if( undoManager )
        undoManager->clearUndoRedu();
      
      break;
    }//case WDialog::Accepted:
  }//switch( result )
  
  if( m_parent )
    m_parent->hide();
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
  //  added/removed; this might necessitate tracking
}//EnergyCalMultiFileModel constructor


EnergyCalMultiFileModel::~EnergyCalMultiFileModel()
{
}//~EnergyCalMultiFileModel()





WModelIndex EnergyCalMultiFileModel::index( int row, int column, const WModelIndex &parent ) const
{
  // Make sure we dont create anything below the peak info level
  
  ModelLevel parent_level, index_level = ModelLevel::Invalid;
  int filenum, samplesnum, peaknum;
  from_internal_id( parent.internalId(), parent_level, filenum, samplesnum, peaknum );
  
  switch( parent_level )
  {
    case ModelLevel::Invalid:
      filenum = row;
      index_level = ModelLevel::File;
      if( column != 0 )
        return WModelIndex();
      break;
      
    case ModelLevel::File:
      samplesnum = row;
      index_level = ModelLevel::SampleSet;
      if( column != 0 )
        return WModelIndex();
      break;
      
    case ModelLevel::SampleSet:
      peaknum = row;
      index_level = ModelLevel::Peak;
      if( (column < 0) || (column >= 3) )
        return WModelIndex();
      break;
      
    case ModelLevel::Peak:
      return WModelIndex();
  }//switch( parent_level )
  
  const uint64_t internal_index = to_internal_id( index_level, filenum, samplesnum, peaknum );
  
  // We could check valid filenum, samplesnum, and peaknum, but we'll just rely on the rest of
  //  the code actually working correctly...
  
  return createIndex( row, column, internal_index );
}//index(...)


WModelIndex EnergyCalMultiFileModel::parent( const WModelIndex &index ) const
{
  if( !index.isValid() )
    return WModelIndex();
  
  const uint64_t internal_id = index.internalId();
  
  ModelLevel parent_level = ModelLevel::Invalid, index_level;
  int filenum, samplesnum, peaknum;
  from_internal_id( internal_id, index_level, filenum, samplesnum, peaknum );
  
  switch( index_level )
  {
    case ModelLevel::Invalid:
    case ModelLevel::File:
      parent_level = ModelLevel::Invalid;
      return WModelIndex();
      
    case ModelLevel::SampleSet:
      samplesnum = peaknum = 0;
      parent_level = ModelLevel::File;
      break;
      
    case ModelLevel::Peak:
      peaknum = 0;
      parent_level = ModelLevel::SampleSet;
      break;
  }//switch( index_level )
  
  const uint64_t parent_index = to_internal_id( parent_level, filenum, samplesnum, peaknum );
  
  switch( parent_level )
  {
    case ModelLevel::File:
      return createIndex( filenum, 0, parent_index );
      
    case ModelLevel::SampleSet:
      return createIndex( samplesnum, 0, parent_index );
      
    case ModelLevel::Peak:
    case ModelLevel::Invalid:
      assert( 0 );
      return WModelIndex();
  }//switch( parent_level )
  
  // We could validate filenum and samplesnum are valid...
  
  return WModelIndex();
}//parent(...)


int EnergyCalMultiFileModel::rowCount( const WModelIndex &parent ) const
{
  if( !parent.isValid() ) // Number of files.
    return static_cast<int>( m_data.size() );
  
  const uint64_t parent_id = parent.internalId();
  
  ModelLevel parent_level;
  int filenum, samplesnum, peaknum;
  from_internal_id( parent_id, parent_level, filenum, samplesnum, peaknum );
  
  const int row = parent.row();
  const int col = parent.column();
  
  switch( parent_level )
  {
    case ModelLevel::Invalid:
      assert( 0 );
      return static_cast<int>( m_data.size() );
      
    case ModelLevel::File:
      assert( row == filenum );
      if( (col != 0) || (row < 0) || (row >= static_cast<int>(m_data.size())) )
        return 0;
      return static_cast<int>( m_data[row].size() );
      
    case ModelLevel::SampleSet:
      if( (col != 0) || (filenum < 0) || (filenum >= static_cast<int>(m_data.size()))
         || (samplesnum < 0) || (samplesnum >= m_data[filenum].size() ) )
        return 0;
      
      return static_cast<int>( get<3>(m_data[filenum][samplesnum]).size() );
      
    case ModelLevel::Peak:
      return 0;
  }//switch( parent_level )
  
  return 0;
}//rowCount(...)


int EnergyCalMultiFileModel::columnCount( const WModelIndex &parent ) const
{
  ModelLevel level;
  int filenum, samplesnum, peaknum;
  from_internal_id( parent.internalId(), level, filenum, samplesnum, peaknum );
  
  switch( level )
  {
    case ModelLevel::Invalid:
    case ModelLevel::File:
      return 1;
      
    case ModelLevel::SampleSet:
      return 3;
      
    case ModelLevel::Peak:
      return 0;
  }//switch( level )
  
  return 0;
}//columnCount(...)

  
boost::any EnergyCalMultiFileModel::data( const Wt::WModelIndex &index, int role ) const
{
  if( (role != Wt::DisplayRole) && (role != Wt::CheckStateRole) )
    return boost::any();
  
  if( !index.isValid() )
    return boost::any();
  
  
  ModelLevel level;
  int filenum, samplesnum, peaknum;
  from_internal_id( index.internalId(), level, filenum, samplesnum, peaknum );
  
  const int row = index.row();
  const int col = index.column();
  
  if( (filenum < 0) || (filenum >= static_cast<int>(m_data.size())) )
    return boost::any();
  
  const vector<SamplesPeakInfo_t> &fileinfos = m_data[filenum];
  
  
  switch( level )
  {
    case ModelLevel::Invalid:
      return boost::any();

    case ModelLevel::File:
    {
      // Return filename here
      if( (col != 0) || (role != Wt::DisplayRole) )
        return boost::any();
        
      assert( filenum == row );
      
      if( fileinfos.empty() )
        return boost::any();
      
      shared_ptr<SpectraFileHeader> header = get<0>( fileinfos[0] );
      if( !header )
        return boost::any();
      
      Wt::WString value = header->displayName();
      if( value.empty() )
        value = WString::fromUTF8( "File " + std::to_string(filenum) );
      
      return boost::any( value );
    }//case ModelLevel::File:
      
      
    case ModelLevel::SampleSet:
    {
      // Return sample numbers + title here
      assert( row == samplesnum );
      if( (col != 0) || (role != Wt::DisplayRole) )
        return boost::any();
      
      if( (samplesnum < 0) || (samplesnum >= static_cast<int>(fileinfos.size())) )
        return boost::any();
      
      const SamplesPeakInfo_t &setinfo = fileinfos[samplesnum];
      
      const shared_ptr<SpectraFileHeader> &header = get<0>(setinfo);
      const set<int> &samplenums = get<1>(setinfo);
      
      string title;
      
      if( header && (samplenums.size() == 1) )
      {
        shared_ptr<SpecMeas> meas = header ? header->measurementIfInMemory() : nullptr;
        if( meas )
        {
          const int sample = *begin(samplenums);
          for( shared_ptr<const SpecUtils::Measurement> &m : meas->sample_measurements(sample) )
          {
            if( m && !m->title().empty() )
            {
              title = m->title();
              break;
            }
          }//for( loop over measurements of this sample number )
          
          if( title.empty() && (meas->sample_numbers().size() > 1) )
            title = "Sample " + std::to_string(sample);
        }//if( we have the measurement )
        
        if( title.empty() )
          title = "Peaks in File";
      }else
      {
        title = "Samples " + SpecUtils::sequencesToBriefString(samplenums);
      }
      
      //if( title.size() > 80 )
      //  title = title.substr(0,77) + "...";
      
      return boost::any( WString::fromUTF8(title) );
    }//case ModelLevel::SampleSet:
      
      
    case ModelLevel::Peak:
    {
      assert( row == peaknum );
      if( col < 0 || col > 2 )
        return boost::any();
      
      if( (samplesnum < 0) || (samplesnum >= static_cast<int>(fileinfos.size())) )
        return boost::any();
      
      const vector<UsePeakInfo_t> &peakinfos = get<3>(fileinfos[samplesnum]);
      if( (peaknum < 0) || (peaknum >= static_cast<int>(peakinfos.size())) )
        return boost::any();
      
      const bool checked = get<0>(peakinfos[peaknum]);
      shared_ptr<const PeakDef> peak = get<1>(peakinfos[peaknum]);
      
      if( role == Wt::CheckStateRole )
      {
        if( col == 0 )
          return boost::any( checked );
        return boost::any();
      }
      
      if( !peak )
        return boost::any();
      
      if( !peak->parentNuclide() && !peak->xrayElement() && !peak->reaction() )
        return boost::any();
      
      if( role == Wt::CheckStateRole )
      {
        if( index.column() == 0 )
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
        
        return boost::any( WString::fromUTF8("--") );
      }//if( col == 0 )
      
      char msg[128] = { '\0' };
      
      if( col == 1 )
      {
        snprintf(msg, sizeof(msg), "%.2f keV Detected", peak->mean() );
      }else if( col == 2 )
      {
        try
        {
          const float wantedEnergy = peak->gammaParticleEnergy();
          snprintf(msg, sizeof(msg), "%.2f keV Expected", wantedEnergy );
        }catch(...)
        {
          snprintf(msg, sizeof(msg), "%s", "--" );
        }
      }else
      {
        return boost::any();
      } //if( col == 1 ) / else ...
      
      return boost::any( WString::fromUTF8(msg) );
      
      break;
    }//case ModelLevel::Peak:
  }//switch( level )
  
    
  return boost::any();
}//data(...)


bool EnergyCalMultiFileModel::setData( const WModelIndex &index, const boost::any &value, int role )
{
  if( (role != Wt::CheckStateRole) || (index.column() != 0) )
    return false;
  
  ModelLevel level;
  int filenum, samplesnum, peaknum;
  from_internal_id( index.internalId(), level, filenum, samplesnum, peaknum );
  
  if( level != ModelLevel::Peak )
    return false;
  
  if( (filenum < 0) || filenum >= static_cast<int>(m_data.size()) )
    return false;
  
  vector<SamplesPeakInfo_t> &filedata = m_data[filenum];

  if( (samplesnum < 0) || (samplesnum >= static_cast<int>(filedata.size())) )
    return false;
  
  vector<UsePeakInfo_t> &peakinfos = get<3>(filedata[samplesnum]);
  
  if( (peaknum < 0) || (peaknum >= static_cast<int>(peakinfos.size())) )
    return false;
  
  try
  {
    get<0>(peakinfos[peaknum]) = boost::any_cast<bool>( value );
  }catch(...)
  {
    cerr << "Got bad boost::any_cast<bool>( value )" << endl;
    return false;
  }
  
  return true;
}//setData(...)



WFlags<ItemFlag> EnergyCalMultiFileModel::flags( const WModelIndex &index ) const
{
  ModelLevel level;
  int filenum, samplesnum, peaknum;
  from_internal_id( index.internalId(), level, filenum, samplesnum, peaknum );
  
  if( (level == ModelLevel::Peak) && (index.column() == 0) )
    return WFlags<ItemFlag>(ItemIsUserCheckable);
  
  //return WFlags<ItemFlag>(ItemIsXHTMLText);
  return WFlags<ItemFlag>();
}//flags(...)


boost::any EnergyCalMultiFileModel::headerData( int section, Wt::Orientation orientation,
                                                int role ) const
{
//  if( role == DisplayRole )
//  {
//    switch( section )
//    {
//      case 0: return boost::any( WString("Nuclide") );
//      case 1: return boost::any( WString("Obs. Energy") );
//      case 2: return boost::any( WString("Photopeak Energy") );
//    }
//  }//if( role == DisplayRole )
  
  return boost::any();
}//headerData(...);


void EnergyCalMultiFileModel::refreshData()
{
  if( !m_data.empty() )
  {
    beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_data.size()) - 1 );
    m_data.clear();
    endRemoveRows();
  }//if( m_peaks.size() )
  
  vector<vector<SamplesPeakInfo_t>> newdata;
  
  const int nfile = m_fileModel->rowCount();
  for( int filenum = 0; filenum < nfile; ++filenum )
  {
    vector<SamplesPeakInfo_t> peaks;
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
      SamplesPeakInfo_t samplesinfo;
      
      get<0>(samplesinfo) = header;
      get<1>(samplesinfo) = samplnums;
      shared_ptr<const SpecUtils::EnergyCalibration> &cal = get<2>(samplesinfo);
      vector<UsePeakInfo_t> &samplespeaks = get<3>(samplesinfo);
      
      try
      {
        cal = spec->suggested_sum_energy_calibration( samplnums, gammaDetNames );
      }catch( std::exception & )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Unexpected failure of suggested_sum_energy_calibration" );
#endif
      }//try / catch
      
      auto measpeaks = spec->peaks( samplnums );
      if( cal && measpeaks )
      {
        for( const PeakPtr &peak : *measpeaks )
        {
          if( peak && peak->hasSourceGammaAssigned() )
            samplespeaks.emplace_back( peak->useForEnergyCalibration(), peak );
        }
      }
      
      if( !samplespeaks.empty() )
        peaks.push_back( samplesinfo );
    }//for( const IntSet &samplnums : peaksamplenums )
    
    newdata.push_back( peaks );
  }//for( int filenum = 0; filenum < nfile; ++filenum )

  m_data.swap( newdata );
  reset();
  
  //if( !newdata.empty() )
  //{
  //  beginInsertRows( WModelIndex(), 0, static_cast<int>(newdata.size()) -1 );
  //  m_data.swap( newdata );
  //  endInsertRows();
  //}//if( newdata.size() )
}//refreshData()
